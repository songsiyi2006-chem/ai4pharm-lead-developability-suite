#!/usr/bin/env python3
"""Task 8: executed synthetic ADC engineering, with explicit mass accounting.

Units: conjugation time h; cell species molecules/cell; tissue amounts molecules;
length um; concentration nM. The HIC retention and PD parameters are hypotheses.
No clinical efficacy or product-specific parameters are estimated by this script.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import nnls
from scipy.special import voigt_profile
from scipy.sparse import lil_matrix, eye
from scipy.sparse.linalg import splu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NA = 6.02214076e23
NM_TO_MOLECULES_UM3 = NA * 1e-24
STATES = [(m, n) for m in (0, 2, 4, 6, 8) for n in range(m + 1)]
INDEX = {state: i for i, state in enumerate(STATES)}
DEFAULT = {
    "seed": 82026, "reduction_rate_h": 3.0, "reduction_time_h": 1.0,
    "conjugation_rate_per_equivalent_h": 0.5, "steric_penalty": 0.18,
    "conjugation_time_h": 8.0, "payload_equivalents": [2.0, 4.0, 6.0, 10.0],
    "downstream_payload_equivalents": 4.0, "ssa_antibodies": 1000, "ssa_replicates": 24,
    "hic_gradient_min": 35.0, "hic_initial_salt_m": 1.5,
    "hic_retention_base_min": 2.0, "hic_retention_scale_min": 2.0,
    "hic_hydrophobic_exponent": 1.15, "hic_sigma_min": 0.35, "hic_gamma_min": 0.12,
    "receptors_per_cell": 1e6, "adc_bath_nm": 10.0, "bath_volume_l_per_cell": 1e-9,
    "cell_volume_um3": 2000.0, "kon_per_nm_h": 0.1, "koff_h": 0.1,
    "internalization_h": 0.05, "receptor_recycling_h": 0.15, "endosome_sorting_h": 0.15,
    "cleavage_vmax_adc_per_cell_h": 50000.0, "cleavage_km_adc_per_cell": 50000.0,
    "payload_metabolism_h": 0.03, "high_permeability_h": 0.2, "low_permeability_h": 0.002,
    "tissue_radius_um": 200.0, "core_radius_um": 20.0, "extracellular_fraction": 0.3,
    "diffusion_um2_h": 300.0, "extracellular_clearance_h": 0.04,
    "pd_ec50_nm": 50.0, "pd_max_kill_h": 0.05, "horizon_h": 72.0,
    "radial_cells": 80, "time_step_h": 0.1, "output_interval_h": 0.5,
}


def validate_config(c):
    if set(c) != set(DEFAULT):
        raise ValueError(f"Config keys must match defaults; differences={set(c)^set(DEFAULT)}")
    for key, val in c.items():
        if key == "payload_equivalents":
            if not isinstance(val, list) or len(val) < 2 or len(set(val)) != len(val):
                raise ValueError("payload_equivalents requires at least two distinct values")
            if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0 for x in val):
                raise ValueError("payload equivalents must be finite positive")
            continue
        if isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val) or val < 0:
            raise ValueError(f"{key} must be a finite nonnegative number")
    for k in ("seed", "ssa_antibodies", "ssa_replicates", "radial_cells"):
        if int(c[k]) != c[k]:
            raise ValueError(f"{k} must be integer")
    for k in ("conjugation_time_h", "hic_gradient_min", "hic_initial_salt_m", "hic_sigma_min", "hic_gamma_min", "receptors_per_cell", "bath_volume_l_per_cell", "cell_volume_um3", "cleavage_km_adc_per_cell", "tissue_radius_um", "core_radius_um", "pd_ec50_nm", "horizon_h", "time_step_h", "output_interval_h"):
        if c[k] <= 0:
            raise ValueError(f"{k} must be positive")
    if not 0 < c["extracellular_fraction"] < 1 or not 0 < c["core_radius_um"] < c["tissue_radius_um"]:
        raise ValueError("Invalid tissue geometry")
    if c["radial_cells"] < 10 or c["ssa_antibodies"] < 100 or c["ssa_replicates"] < 2:
        raise ValueError("Insufficient spatial grid or stochastic replicates")
    if c["downstream_payload_equivalents"] not in c["payload_equivalents"]:
        raise ValueError("downstream_payload_equivalents must be a simulated conjugation condition")
    for ratio in (c["output_interval_h"] / c["time_step_h"], c["horizon_h"] / c["output_interval_h"]):
        if not math.isclose(ratio, round(ratio), abs_tol=1e-9):
            raise ValueError("horizon/output interval/time step must divide exactly")
    rtmax = c["hic_retention_base_min"] + c["hic_retention_scale_min"] * 8**c["hic_hydrophobic_exponent"]
    if c["hic_retention_scale_min"] <= 0 or rtmax >= c["hic_gradient_min"]:
        raise ValueError("HIC peak centers must be distinct and inside gradient")


def reduction_prob(c):
    return -np.expm1(-c["reduction_rate_h"] * c["reduction_time_h"])


def initial_conjugation(c, equivalents):
    p = reduction_prob(c)
    y = np.zeros(len(STATES) + 1)
    for b in range(5):
        y[INDEX[(2*b, 0)]] = math.comb(4, b) * p**b * (1-p)**(4-b)
    y[-1] = equivalents
    return y


def conjugation_rhs(t, y, c):
    del t
    d = np.zeros_like(y)
    linker = max(y[-1], 0.0)
    for i, (m, n) in enumerate(STATES):
        if n < m:
            v = c["conjugation_rate_per_equivalent_h"] * (m-n) * np.exp(-c["steric_penalty"]*n) * linker * max(y[i], 0.0)
            d[i] -= v
            d[INDEX[(m, n+1)]] += v
            d[-1] -= v
    return d


def conjugate(c, equivalents):
    t = np.linspace(0, c["conjugation_time_h"], 161)
    sol = solve_ivp(conjugation_rhs, (t[0], t[-1]), initial_conjugation(c, equivalents), args=(c,), t_eval=t, method="DOP853", rtol=2e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    f = np.zeros((len(t), 9))
    for i, (_, n) in enumerate(STATES):
        f[:, n] += sol.y[i]
    mean = f @ np.arange(9)
    return {"time": t, "fractions": f, "mean": mean, "free_linker": sol.y[-1],
            "mass_error": float(np.max(np.abs(mean + sol.y[-1] - equivalents))),
            "probability_error": float(np.max(np.abs(f.sum(axis=1)-1))), "min_fraction": float(f.min())}


def gillespie(c, equivalents):
    rng = np.random.default_rng(int(c["seed"]))
    nr, nabs = int(c["ssa_replicates"]), int(c["ssa_antibodies"])
    output = np.zeros((nr, 9))
    coefficients = np.array([c["conjugation_rate_per_equivalent_h"]*(m-n)*np.exp(-c["steric_penalty"]*n) for m,n in STATES])
    destinations = np.array([INDEX[(m, min(n+1,m))] for m,n in STATES])
    for rep in range(nr):
        counts = np.zeros(len(STATES), dtype=int)
        bonds = rng.binomial(4, reduction_prob(c), nabs)
        for b in range(5):
            counts[INDEX[(2*b,0)]] = np.count_nonzero(bonds == b)
        linkers = round(equivalents*nabs)
        t = 0.0
        while linkers:
            prop = coefficients * counts * linkers/nabs
            total = prop.sum()
            if total <= 0:
                break
            t += rng.exponential(1/total)
            if t > c["conjugation_time_h"]:
                break
            i = min(int(np.searchsorted(np.cumsum(prop), rng.random()*total)), len(prop)-1)
            counts[i] -= 1
            counts[destinations[i]] += 1
            linkers -= 1
        for i, (_, n) in enumerate(STATES):
            output[rep,n] += counts[i]/nabs
    return output


def hic(c, fractions):
    time_grid = np.linspace(0, c["hic_gradient_min"], 7001)
    centers = c["hic_retention_base_min"] + c["hic_retention_scale_min"]*np.arange(9)**c["hic_hydrophobic_exponent"]
    basis = voigt_profile(time_grid[:,None]-centers, c["hic_sigma_min"], c["hic_gamma_min"])
    # Each species has equal detector response per antibody, normalized within run.
    basis /= np.trapezoid(basis, time_grid, axis=0)
    components = basis*fractions
    signal = components.sum(axis=1)
    fit, residual = nnls(basis, signal)
    fwhm = []
    for j in range(9):
        half = basis[:,j].max()/2
        selected = np.flatnonzero(basis[:,j] >= half)
        fwhm.append(float(time_grid[selected[-1]] - time_grid[selected[0]]))
    rs = 1.18*np.diff(centers)/(np.array(fwhm[:-1])+np.array(fwhm[1:]))
    left, right = (centers[3]+centers[4])/2, (centers[4]+centers[5])/2
    mask = (time_grid>=left)&(time_grid<=right)
    collected = np.trapezoid(components[mask], time_grid[mask], axis=0)
    purity = collected[4]/collected.sum() if collected.sum()>0 else None
    recovery = collected[4]/fractions[4] if fractions[4]>0 else None
    return {"time":time_grid,"centers":centers,"components":components,"signal":signal,
            "fitted_fractions":fit,"fit_residual":float(residual),"fwhm":fwhm,"rs":rs,
            "dar4_window":[float(left),float(right)],"dar4_purity":None if purity is None else float(purity),
            "dar4_recovery":None if recovery is None else float(recovery)}


CELL_NAMES = ["external_adc", "free_surface_receptor", "bound_adc", "endosome_adc", "lysosome_adc", "recycling_receptor", "degraded_adc", "cytoplasm_payload", "exported_payload", "metabolized_payload"]


def cell_rhs(t, y, c, mean_dar, permeability):
    del t
    x, r, b, e, l, ri, z, p, q, met = np.maximum(y,0)
    concentration = x/(NA*1e-9*c["bath_volume_l_per_cell"])
    bind = c["kon_per_nm_h"]*concentration*r
    unbind = c["koff_h"]*b
    internalize = c["internalization_h"]*b
    recycle = c["receptor_recycling_h"]*ri
    sort = c["endosome_sorting_h"]*e
    cleave = c["cleavage_vmax_adc_per_cell_h"]*l/(c["cleavage_km_adc_per_cell"]+l)
    export = permeability*p
    metabolism = c["payload_metabolism_h"]*p
    return np.array([-bind+unbind,-bind+unbind+recycle,bind-unbind-internalize,internalize-sort,sort-cleave,internalize-recycle,cleave,mean_dar*cleave-export-metabolism,export,metabolism])


def cell_model(c, mean_dar, permeability):
    y0 = np.zeros(10)
    y0[0] = c["adc_bath_nm"]*NA*1e-9*c["bath_volume_l_per_cell"]
    y0[1] = c["receptors_per_cell"]
    sol = solve_ivp(cell_rhs,(0,c["horizon_h"]),y0,args=(c,mean_dar,permeability),method="Radau",rtol=2e-9,atol=1e-5,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    t = np.linspace(0,c["horizon_h"],int(round(c["horizon_h"]/.1))+1)
    y=sol.sol(t)
    adc = y[[0,2,3,4,6]].sum(axis=0)
    receptor=y[[1,2,5]].sum(axis=0)
    payload=mean_dar*y[[0,2,3,4]].sum(axis=0)+y[[7,8,9]].sum(axis=0)
    return {"solution":sol,"time":t,"values":y,"initial_adc":float(y0[0]),
            "adc_relative_error":float(np.max(np.abs(adc-y0[0]))/max(1,y0[0])),
            "receptor_relative_error":float(np.max(np.abs(receptor-y0[1]))/max(1,y0[1])),
            "payload_relative_error":float(np.max(np.abs(payload-mean_dar*y0[0]))/max(1,mean_dar*y0[0])),
            "min_species":float(y.min())}


def radial_operator(c, permeability, n=None, absorbing=True):
    n=int(n or c["radial_cells"])
    edges=np.linspace(0,c["tissue_radius_um"],n+1)
    centers=(edges[:-1]+edges[1:])/2
    vol=4*np.pi/3*np.diff(edges**3)
    corevol=4*np.pi/3*np.diff(np.minimum(edges,c["core_radius_um"])**3)
    ve=c["extracellular_fraction"]*vol
    vi=(1-c["extracellular_fraction"])*(vol-corevol)
    matrix=lil_matrix((2*n,2*n))
    dr=edges[1]-edges[0]
    for i in range(n-1):
        g=c["diffusion_um2_h"]*c["extracellular_fraction"]*4*np.pi*edges[i+1]**2/dr
        matrix[i,i]-=g/ve[i]; matrix[i,i+1]+=g/ve[i+1]
        matrix[i+1,i]+=g/ve[i]; matrix[i+1,i+1]-=g/ve[i+1]
    boundary=0.0
    if absorbing:
        boundary=c["diffusion_um2_h"]*c["extracellular_fraction"]*4*np.pi*edges[-1]**2/(dr/2)/ve[-1]
        matrix[n-1,n-1]-=boundary
    for i in range(n):
        uptake=permeability*vi[i]/ve[i]
        matrix[i,i]-=uptake+c["extracellular_clearance_h"]
        matrix[i,n+i]+=permeability
        matrix[n+i,i]+=uptake
        matrix[n+i,n+i]-=permeability+c["payload_metabolism_h"]
    return {"matrix":matrix.tocsc(),"r":centers,"ve":ve,"vi":vi,"core_weights":corevol/corevol.sum(),
            "core_cells":float((1-c["extracellular_fraction"])*corevol.sum()/c["cell_volume_um3"]),"boundary_rate":boundary}


def tissue_model(c, cell, permeability, n=None, dt=None):
    grid=radial_operator(c,permeability,n)
    n=len(grid["r"]); dt=float(dt or c["time_step_h"])
    steps=int(round(c["horizon_h"]/dt)); stride=int(round(c["output_interval_h"]/dt))
    if not np.isclose(steps*dt,c["horizon_h"]) or stride<1 or not np.isclose(stride*dt,c["output_interval_h"]):
        raise ValueError("Refinement steps must divide the output interval and horizon")
    solve=splu(eye(2*n,format="csc")-dt*grid["matrix"]).solve
    amounts=np.zeros(2*n); hazard=np.zeros(n); sink=0.0; supplied=0.0; max_error=0.0
    times=[0.0]; snapshots=[amounts.copy()]; hazards=[hazard.copy()]; ledgers=[[0.0,0.0,0.0,0.0]]
    q=cell["solution"].sol(np.arange(steps+1)*dt)[8]*grid["core_cells"]
    for step in range(1,steps+1):
        added=max(0.0,float(q[step]-q[step-1]))
        rhs=amounts.copy(); rhs[:n]+=added*grid["core_weights"]
        amounts=solve(rhs)
        supplied+=added
        lost=dt*(c["extracellular_clearance_h"]*amounts[:n].sum()+c["payload_metabolism_h"]*amounts[n:].sum()+grid["boundary_rate"]*amounts[n-1])
        sink+=lost
        ci=np.divide(amounts[n:],grid["vi"]*NM_TO_MOLECULES_UM3,out=np.zeros(n),where=grid["vi"]>0)
        hazard+=dt*c["pd_max_kill_h"]*ci/(c["pd_ec50_nm"]+ci)
        err=abs(amounts.sum()+sink-supplied)/max(1.0,supplied)
        max_error=max(max_error,err)
        if step%stride==0:
            times.append(step*dt); snapshots.append(amounts.copy()); hazards.append(hazard.copy())
            ledgers.append([supplied,float(amounts.sum()),sink,err])
    a=np.array(snapshots); hazard=np.array(hazards)
    ce=a[:,:n]/grid["ve"]/NM_TO_MOLECULES_UM3
    ci=np.divide(a[:,n:],grid["vi"]*NM_TO_MOLECULES_UM3,out=np.zeros_like(a[:,n:]),where=grid["vi"]>0)
    survival=np.exp(-hazard)
    mean_survival=(survival[-1]*grid["vi"]).sum()/grid["vi"].sum()
    t=np.array(times)
    source_c=cell["solution"].sol(t)[7]/c["cell_volume_um3"]/NM_TO_MOLECULES_UM3
    source_survival=np.exp(-np.trapezoid(c["pd_max_kill_h"]*source_c/(c["pd_ec50_nm"]+source_c),t))
    above=grid["r"][(ce[-1]>c["pd_ec50_nm"]) & (grid["vi"]>0)]
    return {"time":t,"grid":grid,"ce":ce,"ci":ci,"survival":survival,"ledger":np.array(ledgers),"horizon_h":c["horizon_h"],
            "bystander_survival_at_horizon":float(mean_survival),"source_survival_at_horizon":float(source_survival),
            "threshold_radius_um":float(above.max()) if len(above) else None,"mass_relative_error":max_error,
            "minimum_amount":float(a.min()),"final_amounts":a[-1]}


def json_write(path, value):
    def scalar(x):
        if isinstance(x,np.generic): return x.item()
        raise TypeError(f"Unsupported JSON type: {type(x)}")
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False,default=scalar)+"\n",encoding="utf-8")


def csv_write(path, header, rows):
    with Path(path).open("w",newline="",encoding="utf-8") as stream:
        writer=csv.writer(stream); writer.writerow(header); writer.writerows(rows)


def figures(out, c, conjugations, ss, h, cell_results, tissues):
    folder=out/"figures_task8"; folder.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
    base=conjugations[c["downstream_payload_equivalents"]]
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    for n in range(9): ax[0].plot(base["time"],base["fractions"][:,n],label=f"DAR {n}")
    ax[0].set(xlabel="Conjugation time (h)",ylabel="Antibody fraction",title="A. Finite linker pool; partial reduction")
    ax[0].legend(ncol=3,fontsize=8)
    width=.8/len(conjugations)
    for j,(eq,r) in enumerate(conjugations.items()):
        ax[1].bar(np.arange(9)-.4+width*(j+.5),r["fractions"][-1],width,label=f"{eq:g} eq; mean {r['mean'][-1]:.2f}")
    ax[1].set(xlabel="Drug-to-antibody ratio",ylabel="Final fraction",title="B. All integer DAR states retained",xticks=np.arange(9))
    ax[1].legend(fontsize=8)
    fig.suptitle("Task 8A | Synthetic reduction/conjugation master equations")
    fig.savefig(folder/"fig1_adc_conjugation_dar_distribution.png"); plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4.8),layout="constrained",gridspec_kw={"width_ratios":[2,1]})
    for n in range(9):
        ax[0].plot(h["time"],h["components"][:,n],label=f"DAR {n}",ls="-" if n%2==0 else "--",lw=1.8 if n%2==0 else 1)
    ax[0].plot(h["time"],h["signal"],color="black",lw=2,label="Total")
    ax[0].axhline(0,color="gray",lw=.8)
    ax[0].axvspan(*h["dar4_window"],alpha=.12,color="green",label="DAR4 collection window")
    ax[0].set(xlabel="Retention time (min)",ylabel="Normalized detector response",title="A. Empirical salt-gradient Voigt twin")
    ax[0].legend(ncol=3,fontsize=7)
    ax[1].plot(h["time"],c["hic_initial_salt_m"]*(1-h["time"]/c["hic_gradient_min"]),color="#507997")
    ax[1].set(xlabel="Retention time (min)",ylabel="Salt (M)",title="B. Assumed linear gradient")
    purity="NA" if h["dar4_purity"] is None else f"{100*h['dar4_purity']:.1f}%"
    ax[1].text(.05,.45,f"DAR4 window purity: {purity}\nAdjacent Gaussian-equivalent Rs:\n3/4 = {h['rs'][3]:.2f}; 4/5 = {h['rs'][4]:.2f}\n\nKnown-kernel noiseless inversion\nNot an experimental HIC validation",transform=ax[1].transAxes,fontsize=9,bbox={"facecolor":"white","alpha":.9,"edgecolor":"none"})
    fig.suptitle("Task 8B | All DAR 0–8 peaks, including odd species")
    fig.savefig(folder/"fig2_analytical_hic_chromatogram_twin.png"); plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for label,result in cell_results.items():
        t=result["time"]; y=result["values"]; mask=t<=48
        axes[0,0].plot(t[mask],y[2,mask]/1e6,label=label)
        axes[0,1].plot(t[mask],y[3,mask]/1e6,label=f"{label}: endosome")
        axes[0,1].plot(t[mask],y[4,mask]/1e6,ls="--",label=f"{label}: lysosome")
        flux=c["cleavage_vmax_adc_per_cell_h"]*y[4]/(c["cleavage_km_adc_per_cell"]+y[4])*base["mean"][-1]
        axes[1,0].plot(t[mask],flux[mask]/1e6,label=label)
        axes[1,1].plot(t[mask],y[7,mask]/c["cell_volume_um3"]/NM_TO_MOLECULES_UM3,label=label)
    titles=[("Surface-bound ADC","ADC (million/cell)"),("Endosomal and lysosomal ADC","ADC (million/cell)"),("Lysosomal payload liberation","Payload (million/cell/h)"),("Cytoplasmic free payload","Concentration (nM)")]
    for ax,(title,ylab) in zip(axes.flat,titles): ax.set(title=title,xlabel="Time (h)",ylabel=ylab); ax.legend(fontsize=8)
    fig.suptitle(f"Task 8C | Finite bath; incoming mean DAR = {base['mean'][-1]:.3f}; hypothetical permeability")
    fig.savefig(folder/"fig3_intracellular_lysosomal_release_ode.png"); plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    maximum=max(ti["ce"].max() for ti in tissues.values())
    for j,(label,ti) in enumerate(tissues.items()):
        ax=axes[0,j]
        mesh=ax.pcolormesh(ti["grid"]["r"],ti["time"],ti["ce"],shading="auto",cmap="magma",vmin=0,vmax=maximum)
        if ti["ce"].max()>c["pd_ec50_nm"]:
            contour=ax.contour(ti["grid"]["r"],ti["time"],ti["ce"],levels=[c["pd_ec50_nm"]],colors="cyan",linewidths=1)
            ax.clabel(contour,fmt={c["pd_ec50_nm"]:"assumed EC50"},fontsize=8,manual=[(c["core_radius_um"],c["horizon_h"]*.55)])
        ax.axvline(c["core_radius_um"],color="white",ls="--",lw=.8)
        ax.set(title=f"{label}: extracellular payload",xlabel="Radius (um)",ylabel="Time (h)")
        fig.colorbar(mesh,ax=ax,label="nM")
        mask=ti["grid"]["vi"]>0
        axes[1,0].plot(ti["grid"]["r"][mask],ti["survival"][-1,mask],label=f"{label}: Ag-negative")
        axes[1,1].plot(ti["grid"]["r"],ti["ce"][-1],label=f"{label}: extracellular")
        axes[1,1].plot(ti["grid"]["r"][mask],ti["ci"][-1,mask],ls="--",label=f"{label}: intracellular")
    axes[1,0].set(title=f"Hypothetical response at {c['horizon_h']:g} h",xlabel="Radius (um)",ylabel="Surviving fraction",ylim=(0,1.02))
    axes[1,1].set(title="Final spatial concentration",xlabel="Radius (um)",ylabel="Concentration (nM)")
    axes[1,0].legend(fontsize=8); axes[1,1].legend(fontsize=8)
    fig.suptitle("Task 8D | Conservative spherical transport; threshold contour is not a validated killing radius")
    fig.savefig(folder/"fig4_bystander_killing_spatiotemporal_contour.png"); plt.close(fig)


def numeric_checks(c):
    checks=[]
    def record(name,error,tolerance):
        checks.append({"name":name,"error":float(error),"tolerance":tolerance,"passed":bool(error<=tolerance)})
    conj=conjugate(c,c["downstream_payload_equivalents"])
    record("finite_linker_stoichiometry",conj["mass_error"],1e-8)
    record("antibody_probability",conj["probability_error"],1e-9)
    record("dar_fraction_positivity",max(0,-conj["min_fraction"]),1e-10)
    zero=conjugate(dict(c,conjugation_rate_per_equivalent_h=0),4)
    record("zero_conjugation_limit",np.max(np.abs(zero["fractions"][:,0]-1)),1e-12)
    unreduced=conjugate(dict(c,reduction_rate_h=0),4)
    record("unreduced_no_conjugation",unreduced["mean"].max(),1e-12)
    # Independent-site analytic solution with constant linker (not the finite bath).
    rate=.3; duration=2.0; p=1-np.exp(-rate*duration)
    y0=np.zeros(9);y0[0]=1
    def rhs(t,y):
        v=rate*np.arange(8,0,-1)*y[:-1]
        return np.r_[-v[0],v[:-1]-v[1:],v[-1]]
    result=solve_ivp(rhs,(0,duration),y0,rtol=1e-11,atol=1e-13).y[:,-1]
    exact=np.array([math.comb(8,n)*p**n*(1-p)**(8-n) for n in range(9)])
    record("constant_reagent_binomial_limit",np.max(abs(result-exact)),1e-9)
    cell=cell_model(c,float(conj["mean"][-1]),c["high_permeability_h"])
    for name in ("adc_relative_error","receptor_relative_error","payload_relative_error"):
        record(name,cell[name],1e-8)
    record("cell_species_positivity",max(0,-cell["min_species"]),1e-5)
    blocked=cell_model(dict(c,cleavage_vmax_adc_per_cell_h=0),4,c["high_permeability_h"])
    record("zero_cleavage_zero_payload",np.max(abs(blocked["values"][7:10])),1e-9)
    cc=dict(c,extracellular_clearance_h=0,payload_metabolism_h=0)
    grid=radial_operator(cc,c["high_permeability_h"],20,absorbing=False)
    record("closed_finite_volume_column_balance",np.max(abs(np.asarray(grid["matrix"].sum(axis=0)))),1e-10)
    # Equal concentrations at equilibrium have no net exchange or diffusion.
    uniform=np.r_[grid["ve"],grid["vi"]]
    record("uniform_concentration_equilibrium",np.max(abs(grid["matrix"]@uniform))/max(uniform),1e-10)
    return checks


def run(c,out):
    validate_config(c)
    out=Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {out}")
    out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter(); data=out/"data";data.mkdir()
    json_write(out/"config.json",c)
    conjs={eq:conjugate(c,eq) for eq in c["payload_equivalents"]}
    base=conjs[c["downstream_payload_equivalents"]]; md=float(base["mean"][-1])
    ss=gillespie(c,c["downstream_payload_equivalents"])
    h=hic(c,base["fractions"][-1])
    cells={label:cell_model(c,md,c[key]) for label,key in [("High permeability","high_permeability_h"),("Low permeability","low_permeability_h")]}
    tissues={label:tissue_model(c,cells[label],c[key]) for label,key in [("High permeability","high_permeability_h"),("Low permeability","low_permeability_h")]}
    csv_write(data/"dar_timecourses.csv",["payload_equivalents","time_h",*[f"dar_{n}_fraction" for n in range(9)],"mean_dar","free_linker_equivalents"],([eq,t,*r["fractions"][i],r["mean"][i],r["free_linker"][i]] for eq,r in conjs.items() for i,t in enumerate(r["time"])))
    csv_write(data/"gillespie_replicates.csv",["replicate",*[f"dar_{n}_fraction" for n in range(9)],"mean_dar"],([i,*f,float(f@np.arange(9))] for i,f in enumerate(ss)))
    csv_write(data/"hic_chromatogram.csv",["time_min","salt_m",*[f"dar_{n}_signal" for n in range(9)],"total_signal"],([t,c["hic_initial_salt_m"]*(1-t/c["hic_gradient_min"]),*h["components"][i],h["signal"][i]] for i,t in enumerate(h["time"])))
    csv_write(data/"cell_trafficking.csv",["permeability_scenario","time_h",*[n+"_molecules_per_cell" for n in CELL_NAMES]],([label,t,*r["values"][:,i]] for label,r in cells.items() for i,t in enumerate(r["time"])))
    csv_write(data/"radial_timecourses.csv",["scenario","time_h","radius_um","extracellular_nm","neighbor_intracellular_nm","hypothetical_neighbor_survival","neighbor_cell_volume_um3"],([label,t,r,ti["ce"][i,j],ti["ci"][i,j],ti["survival"][i,j],ti["grid"]["vi"][j]] for label,ti in tissues.items() for i,t in enumerate(ti["time"]) for j,r in enumerate(ti["grid"]["r"])))
    csv_write(data/"tissue_mass_ledger.csv",["scenario","time_h","supplied_molecules","resident_molecules","lost_molecules","relative_error"],([label,t,*ti["ledger"][i]] for label,ti in tissues.items() for i,t in enumerate(ti["time"])))
    checks=numeric_checks(c)
    for label,ti in tissues.items():
        checks.append({"name":label+"_tissue_mass_balance","error":ti["mass_relative_error"],"tolerance":1e-8,"passed":ti["mass_relative_error"]<=1e-8})
        checks.append({"name":label+"_tissue_positivity","error":max(0,-ti["minimum_amount"]),"tolerance":1e-8,"passed":ti["minimum_amount"]>=-1e-8})
    high=tissues["High permeability"]
    refined_t=tissue_model(c,cells["High permeability"],c["high_permeability_h"],dt=c["time_step_h"]/2)
    refined_r=tissue_model(c,cells["High permeability"],c["high_permeability_h"],n=2*c["radial_cells"],dt=c["time_step_h"]/2)
    refinement={"time_step_halving_survival_absolute_change":abs(high["bystander_survival_at_horizon"]-refined_t["bystander_survival_at_horizon"]),"spatial_doubling_at_fine_time_survival_absolute_change":abs(refined_t["bystander_survival_at_horizon"]-refined_r["bystander_survival_at_horizon"]),"fine_time_step_h":c["time_step_h"]/2,"fine_radial_cells":2*c["radial_cells"]}
    for key in list(refinement)[:2]: checks.append({"name":key,"error":refinement[key],"tolerance":0.01,"passed":refinement[key]<=0.01})
    fine_profile=np.interp(refined_t["grid"]["r"],refined_r["grid"]["r"],refined_r["ce"][-1])
    profile_error=float(np.max(abs(refined_t["ce"][-1]-fine_profile))/max(1e-20,np.max(refined_r["ce"][-1])))
    refinement["extracellular_profile_max_error_normalized_by_peak"]=profile_error
    checks.append({"name":"spatial_concentration_profile_refinement","error":profile_error,"tolerance":0.02,"passed":profile_error<=0.02})
    sensitivity=[]
    for key in ("diffusion_um2_h","extracellular_clearance_h"):
        for factor in (.5,2.0):
            cc=dict(c,**{key:c[key]*factor})
            ti=tissue_model(cc,cells["High permeability"],c["high_permeability_h"])
            sensitivity.append({"parameter":key,"factor":factor,"bystander_survival":ti["bystander_survival_at_horizon"],"threshold_radius_um":ti["threshold_radius_um"]})
    json_write(out/"sensitivity.json",sensitivity)
    ssmeans=ss@np.arange(9); stderr=float(ssmeans.std(ddof=1)/np.sqrt(len(ssmeans)))
    summary={"evidence_status":"executed synthetic calculation; no clinical validation", "mean_dar_by_equivalents":{str(eq):float(r["mean"][-1]) for eq,r in conjs.items()},"downstream_mean_dar":md,"reduced_bond_probability":float(reduction_prob(c)),"downstream_final_dar_fractions":base["fractions"][-1].tolist(),"gillespie":{"antibodies_per_replicate":c["ssa_antibodies"],"replicates":c["ssa_replicates"],"mean_dar":float(ssmeans.mean()),"standard_error":stderr,"absolute_ode_difference":abs(float(ssmeans.mean())-md),"fraction_max_absolute_difference":float(np.max(abs(ss.mean(axis=0)-base["fractions"][-1])))},"hic":{"dar4_window_min":h["dar4_window"],"dar4_purity":h["dar4_purity"],"dar4_recovery":h["dar4_recovery"],"adjacent_gaussian_equivalent_resolution":h["rs"].tolist(),"noiseless_known_basis_max_fraction_error":float(max(abs(h["fitted_fractions"]-base["fractions"][-1])))},"cell_conservation":{label:{k:r[k] for k in ("adc_relative_error","receptor_relative_error","payload_relative_error","min_species")} for label,r in cells.items()},"tissue":{label:{k:ti[k] for k in ("bystander_survival_at_horizon","source_survival_at_horizon","threshold_radius_um","mass_relative_error","minimum_amount")} for label,ti in tissues.items()},"refinement":refinement,"computation_counts":{"conjugation_conditions":len(conjs),"gillespie_replicates":len(ss),"nominal_cell_conditions":2,"nominal_tissue_conditions":2,"refinement_tissue_conditions":2,"one_at_a_time_tissue_sensitivities":4,"radial_data_rows":sum(len(ti["time"])*len(ti["grid"]["r"]) for ti in tissues.values())},"verification_passed":all(x["passed"] for x in checks)}
    summary["horizon_h"]=c["horizon_h"]
    json_write(out/"summary.json",summary);json_write(out/"verification.json",checks)
    figures(out,c,conjs,ss,h,cells,tissues)
    json_write(out/"run_metadata.json",{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__,"elapsed_seconds":time.perf_counter()-start,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    files={p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob("*")) if p.is_file() and p.name!="manifest.json"}
    json_write(out/"manifest.json",{"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"files_sha256":files})
    print(json.dumps(summary,indent=2))
    if not summary["verification_passed"]: raise RuntimeError("Numerical checks failed; outputs retained for diagnosis")


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path(__file__).resolve().parent/"results")
    parser.add_argument("--config",type=Path)
    parser.add_argument("--write-example",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args(argv)
    if args.write_example:
        if args.write_example.exists(): raise FileExistsError(args.write_example)
        json_write(args.write_example,DEFAULT);return 0
    c=copy.deepcopy(DEFAULT)
    if args.config:
        supplied=json.loads(args.config.read_text(encoding="utf-8"));c.update(supplied)
    validate_config(c)
    if args.self_test:
        checks=numeric_checks(c);print(json.dumps(checks,indent=2));return 0 if all(r["passed"] for r in checks) else 1
    run(c,args.out);return 0


if __name__=="__main__":
    raise SystemExit(main())
