"""Task 6: reproducible hypothetical ASD thermodynamics and kinetics.

Offline CPU computation; no fitted compound, clinical prediction, or wet-lab data.
Run from repo root: python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --out NEW_DIR
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import least_squares
from scipy.special import expit, log_expit, logit

R = 8.31446261815324  # J mol^-1 K^-1
KB = 1.380649e-23  # J K^-1
NA = 6.02214076e23
ROOT = Path(__file__).resolve().parent
SOURCES = [
    {"id": "SHINETSU", "url": "https://www.setylose.com/fileadmin/download_pfmd/49.pdf", "type": "manufacturer primary brochure", "locator": "Shin-Etsu AQOAT page, Tg table", "supports": "HPMCAS Tg 122 C, DSC second heating 10 C/min; grades have different opening pH", "accessed": "2026-09-20"},
    {"id": "BASF_VA64", "url": "https://download.basf.com/p1/EN_StaticDocuments_6145/en/Technical_Information_Technical_Information_English.pdf", "type": "manufacturer primary technical information", "locator": "page 13, melt extrusion", "supports": "Kollidon VA 64 Tg 101 C; ASD carrier; hygroscopicity context", "accessed": "2026-09-20"},
    {"id": "BASF_SOLUPLUS", "url": "https://pharmaceutical.basf.com/global/en/pharma-solutions/products/soluplus", "type": "manufacturer primary product page", "supports": "Carrier identity and solubilizer/matrix roles; not numerical Tg or Hansen parameters", "accessed": "2026-09-20"},
    {"id": "EVONIK", "url": "https://www.evonik.com/en/products/hc/pr_52000884.html", "type": "manufacturer primary product page", "supports": "EUDRAGIT L100 identity and release above pH 6.0; not numerical Tg", "accessed": "2026-09-20"},
    {"id": "FH_EXPERIMENT", "url": "https://pubmed.ncbi.nlm.nih.gov/31430958/", "type": "original research", "supports": "Measured AAPS boundaries and composition/temperature dependent chi; warns against universal constant-Hansen chi", "accessed": "2026-09-20"},
    {"id": "FH_THERMAL", "url": "https://pubmed.ncbi.nlm.nih.gov/21416468/", "type": "original research", "supports": "Thermal phase diagram modelling and experimental validation in a specific drug/copovidone system", "accessed": "2026-09-20"},
    {"id": "FREE_DRUG", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5972073/", "doi": "10.1016/j.jconrel.2018.04.014", "type": "original research", "supports": "Polymer mixing can reduce amorphous chemical potential and maximum free concentration; free drug drives passive transport", "accessed": "2026-09-20"},
]

DEFAULT = {
    "evidence_class": "hypothetical_parameterized_mechanistic_scenarios_not_experiment",
    "seed": 20260920,
    "drug": {"name": "Hypothetical neutral poorly soluble API", "molar_mass_g_mol": 390.0, "density_g_cm3": 1.30, "molar_volume_cm3_mol": 300.0, "hansen_MPa_half": [20.0, 9.0, 10.0], "Tg_K": 318.15, "N_segments": 1.0, "crystalline_solubility_mg_ml": 0.01, "amorphous_solubility_ratio": 10.0, "water_sorption_a": 0.003},
    "polymers": [
        {"name": "HPMC-AS", "hansen_MPa_half": [18.5, 11.0, 12.0], "density_g_cm3": 1.20, "N_segments": 100.0, "Tg_K": 395.15, "Tg_evidence": "manufacturer typical 122 C; SHINETSU, not a batch measurement", "water_sorption_a": 0.016},
        {"name": "PVP-VA", "hansen_MPa_half": [18.0, 12.0, 7.0], "density_g_cm3": 1.18, "N_segments": 80.0, "Tg_K": 374.15, "Tg_evidence": "manufacturer typical 101 C; BASF_VA64, not a batch measurement", "water_sorption_a": 0.035},
        {"name": "Soluplus", "hansen_MPa_half": [18.0, 8.0, 8.0], "density_g_cm3": 1.10, "N_segments": 120.0, "Tg_K": 343.15, "Tg_evidence": "assumed 70 C for scenario, not verified from cited product page", "water_sorption_a": 0.023},
        {"name": "Eudragit L100", "hansen_MPa_half": [17.5, 12.0, 8.0], "density_g_cm3": 1.25, "N_segments": 150.0, "Tg_K": 468.15, "Tg_evidence": "assumed 195 C for scenario, not verified from cited product page", "water_sorption_a": 0.012},
    ],
    "parameter_provenance": "All numeric drug, Hansen, density, segment, sorption, VFT and kinetic values are analyst assumptions, except the two explicitly attributed typical polymer Tgs. Carrier names do not make these calibrated carrier predictions.",
    "storage": {"conditions": [[298.15, 0.60], [313.15, 0.75]], "loadings": [0.1, 0.2, 0.3, 0.5], "water_Tg_K": 136.0, "water_GT_K": 5.0, "sorption_model": "q_dry=a*RH/(1-0.5*RH), a linear dry-mass blend of component coefficients", "vft_Tg_minus_T0_K": 50.0, "vft_log10_tau0_s": -14.0, "vft_log10_tau_at_Tg_s": 2.0, "induction_to_relaxation_ratio": 1e6, "VFT_status": "hypothetical anchored curve; no measured relaxation or crystallization calibration"},
    "kinetics": {"dose_mg": 100.0, "volume_ml": 900.0, "pH": 6.8, "temperature_K": 310.15, "medium": "assumed pH 6.8 intestinal-style medium; no specified bile salt composition or measured micellar partitioning; NOT a validated FaSSIF experiment", "hours": 6.0, "samples": 721, "diffusivity_cm2_s": 5e-6, "boundary_layer_cm": 0.005, "area_cm2": {"crystalline": 120.0, "amorphous": 800.0, "ASD": 650.0}, "absorption_h_inv": 0.30, "gamma0_J_m2": 0.004, "gamma_increment_J_m2": 0.004, "polymer_half_effect_mg_ml": 0.10, "nucleation_h_inv": 5.0, "growth_h_inv": 0.80, "growth_inhibition_ml_mg": 8.0, "site_loss_h_inv": 0.30, "seed_floor": 0.002, "induction_threshold_fraction_dose": 0.01, "supersaturation_threshold": 1.05, "solver": "DOP853", "rtol": 2e-8, "atol": 1e-10, "max_step_h": 0.02},
}


def validate_config(cfg):
    """Validate the fixed Task-6 experiment design and editable physical parameters."""
    def shape(value, example, where="config"):
        if isinstance(example, dict):
            if not isinstance(value, dict) or set(value) != set(example):
                raise ValueError(f"{where}: keys must match the exported example")
            for key in example:
                shape(value[key], example[key], where+"."+key)
        elif isinstance(example, list):
            if not isinstance(value, list) or len(value) != len(example):
                raise ValueError(f"{where}: list length must be {len(example)}")
            for i,(v,e) in enumerate(zip(value,example)):
                shape(v,e,f"{where}[{i}]")
        elif isinstance(example, (int,float)):
            if isinstance(value, bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError(f"{where}: finite numeric value required")
        elif not isinstance(value,str):
            raise ValueError(f"{where}: text value required")
    shape(cfg,DEFAULT)
    drug,storage,kinetics=cfg["drug"],cfg["storage"],cfg["kinetics"]
    if drug["N_segments"] != 1:
        raise ValueError("Nd=1 is required: the reference site is one drug molecular volume")
    if not math.isclose(drug["density_g_cm3"]*drug["molar_volume_cm3_mol"],drug["molar_mass_g_mol"],rel_tol=1e-6):
        raise ValueError("drug molar volume must equal molecular weight/density")
    positive_drug=("density_g_cm3","molar_volume_cm3_mol","molar_mass_g_mol","Tg_K","crystalline_solubility_mg_ml")
    if any(drug[x]<=0 for x in positive_drug) or drug["amorphous_solubility_ratio"]<1 or drug["water_sorption_a"]<0:
        raise ValueError("invalid drug physical parameters")
    for p,expected in zip(cfg["polymers"],DEFAULT["polymers"]):
        if p["name"] != expected["name"]:
            raise ValueError("the four carrier names/order define this Task-6 experiment")
        if min(p["density_g_cm3"],p["Tg_K"],p["N_segments"])<=0 or p["water_sorption_a"]<0:
            raise ValueError("invalid polymer physical parameters")
    if any(x<0 for component in [drug]+cfg["polymers"] for x in component["hansen_MPa_half"]):
        raise ValueError("Hansen components must be nonnegative")
    if storage["loadings"] != [.1,.2,.3,.5]:
        raise ValueError("Task-6 loading comparison is fixed at 10,20,30,50 wt%")
    if any(t<=0 or not 0<=rh<1 for t,rh in storage["conditions"]):
        raise ValueError("storage conditions require Kelvin >0 and RH in [0,1)")
    if min(storage[x] for x in ("water_Tg_K","water_GT_K","vft_Tg_minus_T0_K","induction_to_relaxation_ratio"))<=0:
        raise ValueError("invalid storage model parameters")
    if storage["vft_log10_tau_at_Tg_s"]<=storage["vft_log10_tau0_s"]:
        raise ValueError("VFT Tg anchor must exceed the prefactor")
    positive=("dose_mg","volume_ml","temperature_K","hours","diffusivity_cm2_s","boundary_layer_cm","polymer_half_effect_mg_ml","rtol","atol","max_step_h")
    nonnegative=("absorption_h_inv","gamma0_J_m2","gamma_increment_J_m2","nucleation_h_inv","growth_h_inv","growth_inhibition_ml_mg","site_loss_h_inv","seed_floor")
    if any(kinetics[x]<=0 for x in positive) or any(kinetics[x]<0 for x in nonnegative):
        raise ValueError("kinetic rates/physical scales outside permitted domain")
    if any(v<=0 for v in kinetics["area_cm2"].values()) or not 0<kinetics["induction_threshold_fraction_dose"]<1 or kinetics["supersaturation_threshold"]<=1:
        raise ValueError("invalid area, induction fraction or supersaturation threshold")
    if kinetics["induction_threshold_fraction_dose"] != .01 or kinetics["supersaturation_threshold"] != 1.05:
        raise ValueError("the defined Task-6 endpoints use 1% dose induction and S>1.05")
    if kinetics["absorption_h_inv"]==0:
        raise ValueError("positive uptake is required for the ER_abs comparison; closed-vessel runs disable uptake separately")
    if kinetics["solver"] != "DOP853" or not isinstance(kinetics["samples"],int) or kinetics["samples"]<3:
        raise ValueError("this experiment uses DOP853 and at least 3 output samples")
    if kinetics["max_step_h"]>kinetics["hours"]:
        raise ValueError("max_step_h must not exceed the observation horizon")
    return cfg


def mass_to_volume_fraction(w, rho_drug, rho_poly):
    w = np.asarray(w)
    if np.any((w < 0) | (w > 1)) or min(rho_drug, rho_poly) <= 0:
        raise ValueError("mass fraction must be in [0,1] and densities positive")
    return (w / rho_drug) / (w / rho_drug + (1-w) / rho_poly)


def volume_to_mass_fraction(phi, rho_drug, rho_poly):
    phi = np.asarray(phi)
    return phi * rho_drug / (phi*rho_drug + (1-phi)*rho_poly)


def chi_hansen(drug, polymer, temperature):
    """MPa * cm3 = J, so the resulting interaction parameter is dimensionless."""
    if temperature <= 0:
        raise ValueError("temperature must be positive Kelvin")
    d = np.asarray(drug["hansen_MPa_half"]) - polymer["hansen_MPa_half"]
    return float(drug["molar_volume_cm3_mol"] * np.dot(d*d, [1.0, .25, .25]) / (R*temperature))


def free_energy(phi, chi, n_poly, n_drug=1.0):
    phi = np.asarray(phi, dtype=float)
    if np.any((phi < 0) | (phi > 1)):
        raise ValueError("volume fraction outside [0,1]")
    a = np.zeros_like(phi)
    b = np.zeros_like(phi)
    np.log(phi, out=a, where=phi > 0)
    np.log(1-phi, out=b, where=phi < 1)
    return phi*a/n_drug + (1-phi)*b/n_poly + chi*phi*(1-phi)


def chemical_potentials(z, chi, n_poly, n_drug=1.0):
    # Logit coordinates retain log(1-phi) when phi numerically rounds to one.
    p, q = expit(z), expit(-np.asarray(z))
    c = 1/n_drug-1/n_poly
    return np.array([log_expit(z)/n_drug+c*q+chi*q*q,
                     log_expit(-np.asarray(z))/n_poly-c*p+chi*p*p])


def phase_boundaries(chi, n_poly, n_drug=1.0):
    critical = 0.5*(1/math.sqrt(n_drug)+1/math.sqrt(n_poly))**2
    if chi <= critical + 1e-9:
        return {"two_phase": False, "chi_critical": critical, "binodal": None, "spinodal": None, "common_tangent_residual": None}
    # 1/(Nd*p)+1/(Np*(1-p))-2*chi = 0.
    roots = np.sort(np.roots([2*chi, 1/n_poly-1/n_drug-2*chi, 1/n_drug]))
    lo, hi = map(float, roots)
    zl, zh = logit(lo), logit(hi)
    def residual(z):
        return chemical_potentials(z[0], chi, n_poly, n_drug)-chemical_potentials(z[1], chi, n_poly, n_drug)
    sol = least_squares(residual, [zl-1.0, zh+2.0], bounds=([-500, zh+1e-9], [zl-1e-9, 500]), xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=1500)
    err = float(max(abs(residual(sol.x))))
    if not sol.success or err > 2e-8:
        raise RuntimeError(f"common-tangent solve failed: chi={chi}, residual={err}")
    bino = expit(sol.x)
    if not (bino[0] < lo < hi < bino[1]):
        raise RuntimeError("binodal does not enclose spinodal")
    return {"two_phase": True, "chi_critical": critical, "binodal": bino.tolist(), "binodal_logit": sol.x.tolist(), "spinodal": [lo, hi], "common_tangent_residual": err}


def gordon_taylor(w, t1, t2, k):
    if min(t1, t2, k) <= 0:
        raise ValueError("Tg must be Kelvin and K positive")
    return (w*t1+k*(1-w)*t2)/(w+k*(1-w))


def storage_result(drug, poly, loading, temp, rh, storage):
    if not 0 <= rh < 1:
        raise ValueError("RH must be a fraction, never water mass fraction")
    k = drug["density_g_cm3"]*drug["Tg_K"]/(poly["density_g_cm3"]*poly["Tg_K"])
    dry = float(gordon_taylor(loading, drug["Tg_K"], poly["Tg_K"], k))
    a = loading*drug["water_sorption_a"]+(1-loading)*poly["water_sorption_a"]
    q = a*rh/(1-.5*rh)
    water = q/(1+q)
    wet = float(gordon_taylor(1-water, dry, storage["water_Tg_K"], storage["water_GT_K"]))
    t0 = wet-storage["vft_Tg_minus_T0_K"]
    # Never evaluate beyond the pole or report a fabricated finite shelf life.
    valid = temp > t0
    logtau = None
    logind = None
    if valid:
        logtau = storage["vft_log10_tau0_s"] + (storage["vft_log10_tau_at_Tg_s"]-storage["vft_log10_tau0_s"])*(wet-t0)/(temp-t0)
        logind = logtau+math.log10(storage["induction_to_relaxation_ratio"])-math.log10(3600)
    return {"polymer": poly["name"], "drug_loading": loading, "storage_K": temp, "RH_fraction": rh, "assumed_water_mass_fraction": water, "dry_Tg_K": dry, "wet_Tg_K": wet, "Tg_minus_storage_K": wet-temp, "heuristic_below_30K_margin": wet-temp < 30, "VFT_T0_K": t0, "VFT_domain_valid": valid, "log10_tau_s_hypothetical": logtau, "log10_induction_h_hypothetical": logind, "shelf_life_status": "not_calibrated_not_predictable"}


def cnt_barrier(supersaturation, gamma, molecular_volume_m3, temp):
    if gamma < 0 or molecular_volume_m3 <= 0 or temp <= 0:
        raise ValueError("invalid CNT parameter")
    if supersaturation <= 1:
        return math.inf
    return 16*math.pi*gamma**3*molecular_volume_m3**2/(3*(KB*temp)**3*math.log(supersaturation)**2)


def simulation(cfg, formulation, loading=.2, absorb=False, gamma_scale=1.0, nucleation_scale=1.0, absorption_scale=1.0, tolerance_scale=1.0):
    if formulation not in ("crystalline", "amorphous", "ASD"):
        raise ValueError("unknown formulation")
    if not 0 < loading <= 1:
        raise ValueError("loading must be in (0,1]")
    d, k = cfg["drug"], cfg["kinetics"]
    dose, volume = k["dose_mg"], k["volume_ml"]
    cs = d["crystalline_solubility_mg_ml"]
    cap = cs if formulation == "crystalline" else cs*d["amorphous_solubility_ratio"]
    polymer_c = dose*(1-loading)/loading/volume if formulation == "ASD" else 0.0
    if formulation == "ASD":
        p = cfg["polymers"][0]
        phi = float(mass_to_volume_fraction(loading, d["density_g_cm3"], p["density_g_cm3"]))
        activity = math.exp(float(chemical_potentials(logit(phi), chi_hansen(d, p, k["temperature_K"]), p["N_segments"])[0])) if loading < 1 else 1.0
        cap *= activity  # Assumes congruent release of a homogeneous dry blend.
    else:
        activity = 1.0
    transport = k["diffusivity_cm2_s"]*k["area_cm2"][formulation]/k["boundary_layer_cm"]*3600  # mL/h
    crystal_transport = k["diffusivity_cm2_s"]*k["area_cm2"]["crystalline"]/k["boundary_layer_cm"]*3600
    gamma = (k["gamma0_J_m2"]+k["gamma_increment_J_m2"]*polymer_c/(polymer_c+k["polymer_half_effect_mg_ml"]))*gamma_scale
    growth = k["growth_h_inv"]/(1+k["growth_inhibition_ml_mg"]*polymer_c)
    ka = k["absorption_h_inv"]*absorption_scale if absorb else 0.0
    molecular_volume = d["molar_volume_cm3_mol"]*1e-6/NA
    def rhs(_t, y):
        undiss, dissolved, precip, absorbed, sites = np.maximum(y, 0)
        concentration = dissolved/volume
        saturation = concentration/cs
        area_fraction = (undiss/dose)**(2/3)
        dissolve = transport*area_fraction*max(cap-concentration, 0)
        redissolve = crystal_transport*(precip/dose)**(2/3)*max(cs-concentration, 0)
        barrier = cnt_barrier(saturation, gamma, molecular_volume, k["temperature_K"])
        # Sites are a phenomenological growth-site activation variable, not particle counts.
        activate = k["nucleation_h_inv"]*nucleation_scale*max(saturation-1, 0)*math.exp(-min(barrier, 745))
        site_loss = k["site_loss_h_inv"]*max(1-saturation, 0)*sites
        crystallize = growth*(sites+k["seed_floor"]*(precip/dose)**(2/3))*max(concentration-cs, 0)*volume
        uptake = ka*dissolved
        return [-dissolve, dissolve-crystallize+redissolve-uptake, crystallize-redissolve, uptake, activate-site_loss]
    time = np.linspace(0, k["hours"], k["samples"])
    sol = solve_ivp(rhs, (0,k["hours"]), [dose,0,0,0,0], t_eval=time, method=k["solver"], rtol=k["rtol"]*tolerance_scale, atol=k["atol"]*tolerance_scale, max_step=k["max_step_h"]*math.sqrt(tolerance_scale))
    if not sol.success:
        raise RuntimeError(sol.message)
    mass_error = float(np.max(abs(sol.y[:4].sum(axis=0)-dose)))
    min_mass = float(np.min(sol.y[:4]))
    if mass_error > 1e-6 or min_mass < -1e-6:
        raise RuntimeError(f"mass/positivity failure: {mass_error}, {min_mass}")
    c = sol.y[1]/volume
    auc = float(trapezoid(c,time))
    excess = float(trapezoid(np.maximum(c-cs,0),time))
    threshold = dose*k["induction_threshold_fraction_dose"]
    hits = np.flatnonzero(sol.y[2] >= threshold)
    induction = None
    if len(hits):
        i = hits[0]
        induction = float(np.interp(threshold, sol.y[2, i-1:i+1], time[i-1:i+1]))
    # Linear crossings on the saved grid, not a count of points times dt.
    duration = 0.0
    threshold_c = cs*k["supersaturation_threshold"]
    for t1,t2,c1,c2 in zip(time[:-1],time[1:],c[:-1],c[1:]):
        if c1 >= threshold_c and c2 >= threshold_c:
            duration += t2-t1
        elif (c1-threshold_c)*(c2-threshold_c) < 0:
            fraction = (threshold_c-c1)/(c2-c1)
            duration += (t2-t1)*((1-fraction) if c2 > c1 else fraction)
    summary = {"formulation": formulation, "drug_loading": loading if formulation == "ASD" else 1.0, "absorption_enabled": absorb, "gamma_scale": gamma_scale, "nucleation_scale": nucleation_scale, "absorption_scale": absorption_scale, "drug_activity_in_dry_ASD": activity, "source_solubility_mg_ml": cap, "polymer_concentration_mg_ml_assumed": polymer_c, "gamma_J_m2_assumed": gamma, "AUC_mg_h_ml": auc, "excess_AUC_mg_h_ml": excess, "Cmax_mg_ml": float(c.max()), "absorbed_mg": float(sol.y[3,-1]), "fraction_absorbed_proxy": float(sol.y[3,-1]/dose), "induction_1pct_precip_h": induction, "induction_right_censored": induction is None, "supersaturation_over_1p05_duration_h": duration, "max_mass_error_mg": mass_error, "minimum_compartment_mg": min_mass, "solver_nfev": sol.nfev}
    return {"time": time, "states": sol.y, "concentration": c, "summary": summary}


def write_csv(path, rows):
    if not rows:
        raise ValueError("empty table")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")


def make_figures(cfg, out, phase_rows, storage_rows, curves, loading_runs):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white"})
    folder = out/"figures_task6"
    folder.mkdir()
    colors = ["#166b8f", "#de7044", "#39866d", "#815998"]
    fig, axs = plt.subplots(2,4,figsize=(16,8),layout="constrained")
    drug = cfg["drug"]
    for i,p in enumerate(cfg["polymers"]):
        w = np.linspace(0,1,501)
        phi = mass_to_volume_fraction(w,drug["density_g_cm3"],p["density_g_cm3"])
        ch = chi_hansen(drug,p,298.15)
        axs[0,i].plot(w,free_energy(phi,ch,p["N_segments"]),color=colors[i])
        axs[0,i].set(title=f"{p['name']} | chi(25 C)={ch:.3f}",xlabel="Drug mass fraction",ylabel="Mixing free energy / RT per site")
        rows = [r for r in phase_rows if r["polymer"]==p["name"] and r["two_phase"]]
        if rows:
            t = np.array([r["temperature_K"]-273.15 for r in rows])
            b1,b2,s1,s2 = [np.array([r[k] for r in rows]) for k in ("binodal_low_w","binodal_high_w","spinodal_low_w","spinodal_high_w")]
            ax=axs[1,i]
            ax.fill_betweenx(t,b1,b2,color="#f1ddaf",alpha=.7,label="Two-phase; metastable margins")
            ax.fill_betweenx(t,s1,s2,color="#df9380",alpha=.75,label="Unstable homogeneous state")
            ax.plot(b1,t,"-",color=colors[i],label="Binodal (common tangent)")
            ax.plot(b2,t,"-",color=colors[i])
            ax.plot(s1,t,"--",color="#5b3d40",label="Spinodal")
            ax.plot(s2,t,"--",color="#5b3d40")
        else:
            axs[1,i].text(.5,.5,"One amorphous phase\nthroughout this T range",ha="center",va="center",transform=axs[1,i].transAxes)
        axs[1,i].set(xlim=(0,1),ylim=(0,150),xlabel="Drug mass fraction",ylabel="Temperature (C)")
    axs[1,1].legend(fontsize=7,loc="lower left")
    fig.suptitle("Hypothetical dry binary ASD: equilibrium miscibility, not crystallization stability",fontsize=15)
    fig.savefig(folder/"fig1_flory_huggins_miscibility_phase_diagram.png",dpi=300)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for ax,(temp,rh) in zip(axs,cfg["storage"]["conditions"]):
        for color,p in zip(colors,cfg["polymers"]):
            w=np.linspace(0,1,101)
            vals=[storage_result(drug,p,float(x),temp,rh,cfg["storage"]) for x in w]
            ax.plot(100*w,[v["wet_Tg_K"]-273.15 for v in vals],color=color,label=p["name"])
            ax.plot(100*w,[v["dry_Tg_K"]-273.15 for v in vals],color=color,ls=":",alpha=.5)
        ax.axhspan(-30,temp+30-273.15,color="#f5d7d0",alpha=.6,label="Tg - storage <30 K heuristic")
        ax.axhline(temp-273.15,color="#555555",lw=1)
        ax.set(xlim=(0,100),ylim=(0,205),xlabel="Drug loading (wt%)",ylabel="Glass transition (C)",title=f"{temp-273.15:.0f} C / {rh*100:.0f}% RH | assumed sorption")
    axs[0].legend(fontsize=8)
    fig.suptitle("Moisture plasticization scenarios | solid: humid; dotted: dry",fontsize=14)
    fig.savefig(folder/"fig2_gordon_taylor_tg_depression.png",dpi=300)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for color,form in zip(colors,("crystalline","amorphous","ASD")):
        for i,absorb in enumerate((False,True)):
            sim=curves[(form,absorb)]
            axs[i].plot(sim["time"],sim["concentration"]*1000,color=color,label="20% HPMC-AS" if form=="ASD" else form.capitalize())
    for ax in axs:
        ax.axhline(drug["crystalline_solubility_mg_ml"]*1000,color="black",lw=1,ls="--",label="Crystalline solubility")
        ax.set(xlim=(0,cfg["kinetics"]["hours"]),xlabel="Time (h)",ylabel="Free molecular concentration (ug/mL)")
    axs[0].set_title(f"Closed {cfg['kinetics']['volume_ml']:g} mL dissolution vessel")
    axs[1].set_title("Same vessel + assumed absorption sink")
    axs[0].legend(fontsize=8)
    fig.suptitle("Hypothetical spring-and-parachute kinetics | no validated FaSSIF or clinical PK",fontsize=13)
    fig.savefig(folder/"fig3_spring_and_parachute_dissolution.png",dpi=300)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    x=np.arange(len(loading_runs))
    labels=[f"{s['summary']['drug_loading']*100:.0f}%\nP:D {(1-s['summary']['drug_loading'])/s['summary']['drug_loading']:.2g}:1" for s in loading_runs]
    times=[s["summary"]["induction_1pct_precip_h"] or cfg["kinetics"]["hours"] for s in loading_runs]
    bars=axs[0].bar(x,times,color=colors)
    for b,s in zip(bars,loading_runs):
        if s["summary"]["induction_right_censored"]:
            b.set_hatch("//")
            horizon=cfg["kinetics"]["hours"]
            axs[0].text(b.get_x()+b.get_width()/2,horizon*1.01,f">{horizon:g} h",ha="center",fontsize=9)
    base=curves[("crystalline",False)]["summary"]["AUC_mg_h_ml"]
    ratios=[s["summary"]["AUC_mg_h_ml"]/base for s in loading_runs]
    axs[1].bar(x,ratios,color=colors)
    axs[1].axhline(1,color="black",lw=1,ls="--")
    for ax in axs:
        ax.set_xticks(x,labels)
        ax.set_xlabel("HPMC-AS drug loading and polymer:drug mass ratio")
    axs[0].set(ylabel=f"First precipitated mass >={cfg['kinetics']['induction_threshold_fraction_dose']*100:g}% dose (h)",ylim=(0,cfg["kinetics"]["hours"]*1.12),title="Operational solution induction; hatched = censored")
    axs[1].set(ylabel=f"AUC / crystalline AUC (0-{cfg['kinetics']['hours']:g} h)",title="Closed-vessel molecular concentration AUC")
    fig.suptitle("Hypothetical loading comparison: precipitation inhibition vs reduced drug activity",fontsize=13)
    fig.savefig(folder/"fig4_polymeric_precipitation_inhibition_efficiency.png",dpi=300)
    plt.close(fig)


def run(out, cfg=None):
    cfg=validate_config(json.loads(json.dumps(DEFAULT if cfg is None else cfg)))
    if out.exists():
        raise FileExistsError("Use a new output directory; existing results are never overwritten")
    out.mkdir(parents=True)
    inputs=out/"inputs"
    results=out/"results"
    inputs.mkdir(); results.mkdir()
    write_json(inputs/"parameters.json",cfg)
    write_json(inputs/"sources.json",SOURCES)
    phase_rows=[]
    max_tangent=0.0
    for p in cfg["polymers"]:
        for temp in np.linspace(273.15,423.15,81):
            ch=chi_hansen(cfg["drug"],p,float(temp))
            phases=phase_boundaries(ch,p["N_segments"])
            bd=[None,None]; sp=[None,None]
            if phases["two_phase"]:
                bd=volume_to_mass_fraction(np.array(phases["binodal"]),cfg["drug"]["density_g_cm3"],p["density_g_cm3"]).tolist()
                sp=volume_to_mass_fraction(np.array(phases["spinodal"]),cfg["drug"]["density_g_cm3"],p["density_g_cm3"]).tolist()
                max_tangent=max(max_tangent,phases["common_tangent_residual"])
            phase_rows.append({"polymer":p["name"],"temperature_K":float(temp),"chi":ch,"chi_critical":phases["chi_critical"],"two_phase":phases["two_phase"],"binodal_low_w":bd[0],"binodal_high_w":bd[1],"spinodal_low_w":sp[0],"spinodal_high_w":sp[1],"common_tangent_residual":phases["common_tangent_residual"]})
    write_csv(results/"phase_boundaries.csv",phase_rows)
    storage_rows=[storage_result(cfg["drug"],p,w,t,rh,cfg["storage"]) for p,w,(t,rh) in itertools.product(cfg["polymers"],cfg["storage"]["loadings"],cfg["storage"]["conditions"])]
    write_csv(results/"storage_scenarios.csv",storage_rows)
    curves={(f,a):simulation(cfg,f,absorb=a) for f,a in itertools.product(("crystalline","amorphous","ASD"),(False,True))}
    rows=[]
    for (form,absorb),sim in curves.items():
        for i,t in enumerate(sim["time"]):
            rows.append({"formulation":form,"absorption_enabled":absorb,"time_h":float(t),"undissolved_mg":float(sim["states"][0,i]),"dissolved_mg":float(sim["states"][1,i]),"precipitated_mg":float(sim["states"][2,i]),"absorbed_mg":float(sim["states"][3,i]),"growth_site_activation":float(sim["states"][4,i]),"free_concentration_mg_ml":float(sim["concentration"][i])})
    write_csv(results/"concentration_mass_timeseries.csv",rows)
    summaries=[sim["summary"] for sim in curves.values()]
    write_csv(results/"formulation_summary.csv",summaries)
    loading_runs=[curves[("ASD",False)] if w==.2 else simulation(cfg,"ASD",loading=w) for w in cfg["storage"]["loadings"]]
    write_csv(results/"loading_summary.csv",[x["summary"] for x in loading_runs])
    sensitivity=[]
    for gamma,nuc,uptake in itertools.product((.8,1.,1.2),(.1,1.,10.),(.5,1.,2.)):
        test=simulation(cfg,"ASD",absorb=True,gamma_scale=gamma,nucleation_scale=nuc,absorption_scale=uptake)
        baseline=simulation(cfg,"crystalline",absorb=True,absorption_scale=uptake)
        s=test["summary"].copy()
        s["ER_abs_vs_matched_crystal"]=s["absorbed_mg"]/baseline["summary"]["absorbed_mg"]
        sensitivity.append(s)
    write_csv(results/"sensitivity_27_scenarios.csv",sensitivity)
    checks=[]
    for key,sim in curves.items():
        refined=simulation(cfg,key[0],absorb=key[1],tolerance_scale=.1)
        delta=float(np.max(abs(refined["concentration"]-sim["concentration"])))
        checks.append({"formulation":key[0],"absorption_enabled":key[1],"max_concentration_difference_mg_ml":delta,"AUC_relative_difference":abs(refined["summary"]["AUC_mg_h_ml"]-sim["summary"]["AUC_mg_h_ml"])/sim["summary"]["AUC_mg_h_ml"]})
    write_csv(results/"solver_refinement.csv",checks)
    make_figures(cfg,out,phase_rows,storage_rows,curves,loading_runs)
    crystal_auc=curves[("crystalline",False)]["summary"]["AUC_mg_h_ml"]
    crystal_abs=curves[("crystalline",True)]["summary"]["absorbed_mg"]
    summary={"evidence_class":cfg["evidence_class"],"phase_grid_rows":len(phase_rows),"two_phase_grid_rows":sum(r["two_phase"] for r in phase_rows),"storage_rows":len(storage_rows),"VFT_outside_domain_rows":sum(not r["VFT_domain_valid"] for r in storage_rows),"primary_ODE_runs":6,"additional_loading_runs":3,"sensitivity_paired_ODE_runs":54,"solver_refinement_ODE_runs":6,"total_ODE_runs":69,"random_sampling":False,"seed_recorded_but_unused":cfg["seed"],"max_common_tangent_residual":max_tangent,"max_mass_error_mg":max(s["max_mass_error_mg"] for s in summaries+sensitivity+[r["summary"] for r in loading_runs]),"max_solver_refinement_concentration_difference_mg_ml":max(c["max_concentration_difference_mg_ml"] for c in checks),"ASD20_closed_vessel":curves[("ASD",False)]["summary"],"ASD20_AUC_ratio":curves[("ASD",False)]["summary"]["AUC_mg_h_ml"]/crystal_auc,"ASD20_ER_abs_proxy":curves[("ASD",True)]["summary"]["absorbed_mg"]/crystal_abs,"sensitivity_ER_abs_range":[min(s["ER_abs_vs_matched_crystal"] for s in sensitivity),max(s["ER_abs_vs_matched_crystal"] for s in sensitivity)],"clinical_exposure_status":"not_identified; no permeability/PK/first-pass data","optimization_status":"no validated optimum; 20% HPMC-AS is a prespecified candidate"}
    write_json(results/"summary.json",summary)
    manifest={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"deterministic":True,"files":{p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob("*")) if p.is_file()}}
    write_json(out/"manifest.json",manifest)
    print(json.dumps(summary,indent=2))


def self_test():
    suite=unittest.defaultTestLoader.discover(str(ROOT.parents[1]/"tests"),pattern="test_task6.py")
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if result.testsRun == 0:
        raise RuntimeError("No Task 6 tests discovered")
    return result.wasSuccessful()


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,help="new output directory")
    parser.add_argument("--config",type=Path,help="complete JSON configuration matching --write-example")
    parser.add_argument("--write-example",type=Path,help="write editable default configuration to a new file and exit")
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.write_example is not None:
        if args.write_example.exists():
            parser.error("example configuration path already exists")
        args.write_example.parent.mkdir(parents=True,exist_ok=True)
        write_json(args.write_example,DEFAULT)
        sys.exit(0)
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if args.out is None:
        parser.error("--out NEW_DIR is required")
    configuration=json.loads(args.config.read_text(encoding="utf-8")) if args.config else None
    run(args.out.resolve(),configuration)
