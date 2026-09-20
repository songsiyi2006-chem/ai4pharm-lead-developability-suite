#!/usr/bin/env python3
"""Task 5: evidence-labelled CYP inhibition scenarios coupled to the actual Task 4.

All numeric compound/PK inputs are synthetic demonstration assumptions. Drug names
are qualitative anchors, not validated drug models. Offline; Python >=3.10,
NumPy/SciPy/Matplotlib. --out requires a new directory; --self-test writes nothing.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
import platform
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

HERE = Path(__file__).resolve().parent
TASK4_PATH = HERE.parent / 'task4_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py'
spec = importlib.util.spec_from_file_location('task4_for_task5', TASK4_PATH)
t4 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = t4
spec.loader.exec_module(t4)
ISO = ('CYP3A4', 'CYP2D6', 'CYP2C9')
KDEG = np.array([0.019, 0.025, 0.020])
KDEG_GUT = np.array([0.030, 0.040, 0.030])
LAST_DOSE = 324.0  # 28 BID doses at 0, 12, ... 324 h, fourteen treatment days
END = LAST_DOSE + 240.0
SOURCE_M12 = 'https://www.fda.gov/media/161199/download'
SOURCE_TABLE = 'https://www.fda.gov/drugs/drug-interactions-labeling/drug-development-and-drug-interactions-table-substrates-inhibitors-and-inducers'


def panel():
    # Each numeric entry is a constructed scenario, NEVER an extracted Ki/PK value.
    rows = [
        ('Ketoconazole', 531.43, 200., .04, 12., [.03, 8., 8.], [0., 0., 0.], [1., 1., 1.], 'CYP3A reversible-inhibitor archetype'),
        ('Clarithromycin', 747.95, 500., .15, 25., [8., 80., 80.], [.45, 0., 0.], [2., 1., 1.], 'CYP3A time-dependent-inhibitor archetype'),
        ('Ritonavir', 720.94, 100., .02, 15., [.025, 1., 20.], [.8, .02, 0.], [.15, 3., 1.], 'CYP3A mixed-inhibition archetype'),
        ('Fluconazole', 306.27, 100., .20, 10., [2., 30., .5], [0., 0., 0.], [1., 1., 1.], 'moderate CYP3A/CYP2C9 qualitative anchor'),
        ('Quinidine', 324.42, 200., .15, 25., [15., .02, 30.], [0., 0., 0.], [1., 1., 1.], 'CYP2D6 reversible-inhibitor archetype'),
        ('Amoxicillin', 365.40, 250., .20, 15., [1e12, 1e12, 1e12], [0., 0., 0.], [1., 1., 1.], 'assumed near-zero CYP-inhibition control; no clinical exclusion claim'),
    ]
    return [dict(name=n, mw_g_mol=mw, dose_mg=d, fu_p=fu,
                 clint_mic_ul_min_mg=cl, ki_uM=ki, kinact_h=kin, KI_uM=KI,
                 role=role, evidence='synthetic_scenario', ka_h=1., fa=.9,
                 fg=1., gfr_l_h=7.5, qgut_l_h=18., fu_gut=1.)
            for n,mw,d,fu,cl,ki,kin,KI,role in rows]


def probes():
    return [dict(name='Midazolam', dose_mg=2., mw_g_mol=325.77, fu_p=.05,
                 clint_mic_ul_min_mg=150., ka_h=1.5, fa=.95, fg=.5,
                 hepatic_intrinsic_fractions=[.9,0.,0.], gut_fractions=[1.,0.,0.]),
            dict(name='Metoprolol', dose_mg=50., mw_g_mol=267.37, fu_p=.20,
                 clint_mic_ul_min_mg=100., ka_h=1.2, fa=.95, fg=1.,
                 hepatic_intrinsic_fractions=[0.,.9,0.], gut_fractions=[0.,0.,0.])]


def validate_compound(d):
    if set(d)!=set(panel()[0]):
        raise ValueError('compound keys must match --write-example schema exactly')
    if d['evidence']!='synthetic_scenario':
        raise ValueError('This demonstration schema retains synthetic_scenario evidence; measured-input validation is not implemented')
    for key in ('mw_g_mol','dose_mg','fu_p','clint_mic_ul_min_mg','ka_h','qgut_l_h','fu_gut'):
        if not math.isfinite(d[key]) or d[key] <= 0:
            raise ValueError(f'{key} must be finite and positive')
    if not .01 <= d['fu_p'] <= .2 or not 0 < d['fa'] <= 1 or not 0 < d['fg'] <= 1 or not 0 < d['fu_gut'] <= 1:
        raise ValueError('invalid binding/absorption fraction (Task 4 fu_p range is 0.01..0.20)')
    for key in ('ki_uM','kinact_h','KI_uM'):
        a = np.asarray(d[key], float)
        if a.shape != (3,) or not np.all(np.isfinite(a)) or np.any(a < 0) or (key != 'kinact_h' and np.any(a == 0)):
            raise ValueError(f'{key}: three nonnegative values, with positive Ki and KI, required')


def pbpk_model(d):
    # Reuse actual Task 4 Config, PBPK partitioning, circulation and derivative.
    c = t4.Config(compound_name=d['name']+'_SYNTHETIC',
                  mw_g_mol=d['mw_g_mol'], fu_p=d['fu_p'],
                  clint_mic_ul_min_mg=d['clint_mic_ul_min_mg'],
                  ka_h=d['ka_h'], fa=d['fa'], fg=d['fg'],
                  gfr_l_h=d.get('gfr_l_h',7.5),
                  kp_overrides={n: 1. for n in t4.TISSUES})
    return t4.PBPK(c)


def mg_l_to_uM(value, mw):
    return np.asarray(value)*1000./mw


def kobs(i, kinact, KI):
    return np.asarray(kinact)*np.asarray(i)/(np.asarray(KI)+np.asarray(i))


def activity(e, i, ki):
    return np.asarray(e)/(1.+np.asarray(i)/np.asarray(ki))


def official_static_aucr(ahep, agut, fm=.9, fg=.5):
    """M12 inhibition-only static limit, negligible extrahepatic clearance.

    A*B here is active-enzyme fraction/(1+I/Ki); induction factor C=1.
    fm is hepatic clearance fraction, not an independently estimated systemic fm.
    """
    if not 0 <= fm <= 1 or not 0 < fg <= 1 or not 0 <= ahep <= 1 or not 0 <= agut <= 1:
        raise ValueError('invalid static fraction')
    den = (ahep*fm+1-fm)*(agut*(1-fg)+fg)
    return math.inf if den == 0 else 1./den


def severity(r):
    if not math.isfinite(r) or r <= 0:
        raise ValueError('finite positive AUCR required')
    return 'strong_range' if r >= 5 else 'moderate_range' if r >= 2 else 'weak_range' if r >= 1.25 else 'below_weak_threshold'


def basic_ratios(cmax_uM, ki_uM, kinact_h, KI_uM, kdeg_h):
    """M12 basic hepatic screens, distinct from both static and dynamic AUCR."""
    return 1+cmax_uM/ki_uM, 1+float(kobs(5*cmax_uM,kinact_h,KI_uM))/kdeg_h


class Exposure:
    """Piecewise dense 17-state solution: Task 4 eleven + six normalized CYP pools."""
    def __init__(self, d, kdeg_factor=1., kinact_factor=1., rtol=2e-7, atol=1e-9):
        validate_compound(d)
        self.d, self.model = d, pbpk_model(d)
        self.kdeg = KDEG*kdeg_factor
        self.kgut = KDEG_GUT*kdeg_factor
        self.kinact = np.asarray(d['kinact_h'])*kinact_factor
        self.segments=[]
        state=np.r_[np.zeros(11),np.ones(6)]
        edges=list(np.arange(0,LAST_DOSE+1,12.))+[END]
        mass_error=0.
        min_amount=0.
        for j,(start,end) in enumerate(zip(edges[:-1],edges[1:])):
            state[t4.GUT] += d['dose_mg']
            expected=(j+1)*d['dose_mg']
            sol=solve_ivp(self.rhs,(start,end),state,method='BDF',rtol=rtol,atol=atol,max_step=2.,dense_output=True)
            if not sol.success or not np.all(np.isfinite(sol.y)):
                raise RuntimeError(sol.message)
            mass_error=max(mass_error,float(np.max(np.abs(sol.y[:10].sum(axis=0)-expected))))
            min_amount=min(min_amount,float(np.min(sol.y[:7])))
            if np.min(sol.y[11:]) < -1e-7 or np.max(sol.y[11:]) > 1+1e-7:
                raise RuntimeError('enzyme state outside [0,1]')
            state=sol.y[:,-1].copy()
            self.segments.append(sol)
        if mass_error > 2e-6*28*d['dose_mg'] or min_amount < -1e-7:
            raise RuntimeError('perpetrator mass balance/positivity failed')
        self.mass_error=mass_error
        self.min_amount=min_amount
        self.ends=np.asarray([s.t[-1] for s in self.segments])

    def concentrations(self,y):
        m=self.model;d=self.d
        ih=mg_l_to_uM(d['fu_p']*np.maximum(y[t4.LIVER],0)/(m.c.volumes_l['liver']*m.kp['liver']),d['mw_g_mol'])
        ip=mg_l_to_uM(d['fu_p']*np.maximum(y[t4.BLOOD],0)*m.cp_factor,d['mw_g_mol'])
        # Explicit gut inlet proxy; Task 4 has no perfused enterocyte compartment.
        ig=ip+mg_l_to_uM(d['fu_gut']*d['fa']*d['ka_h']*np.maximum(y[t4.GUT],0)/d['qgut_l_h'],d['mw_g_mol'])
        return ih,ig,ip

    def rhs(self,t,y):
        ih,ig,_=self.concentrations(y)
        return np.r_[self.model.derivative(y[:11]),
                     self.kdeg*(1-y[11:14])-kobs(ih,self.kinact,self.d['KI_uM'])*y[11:14],
                     self.kgut*(1-y[14:17])-kobs(ig,self.kinact,self.d['KI_uM'])*y[14:17]]

    def at(self,t):
        scalar=np.ndim(t)==0
        times=np.atleast_1d(t).astype(float)
        if np.any(times < 0) or np.any(times > END):
            raise ValueError('time outside exposure horizon')
        # Dose discontinuities are right-continuous; enzymes and blood are continuous.
        idx=np.searchsorted(self.ends,times,side='right').clip(max=len(self.segments)-1)
        y=np.empty((17,len(times)))
        for i in np.unique(idx):
            mask=idx==i;y[:,mask]=self.segments[i].sol(times[mask])
        return y[:,0] if scalar else y

    def factors(self,t):
        y=self.at(t);ih,ig,_=self.concentrations(y)
        return activity(y[11:14],ih,self.d['ki_uM']),activity(y[14:17],ig,self.d['ki_uM'])


def victim_rhs(model,p,y,ah,ag):
    """Task 4 derivative, changing only hepatic metabolism and gut first pass."""
    dy=model.derivative(y)
    w=np.asarray(p['hepatic_intrinsic_fractions'])
    scale=1-w.sum()+np.dot(w,ah)
    cu=model.c.fu_p*y[t4.LIVER]/(model.c.volumes_l['liver']*model.kp['liver'])
    delta=model.clint*cu*(scale-1)
    dy[t4.LIVER]-=delta;dy[t4.HEP]+=delta
    wg=np.asarray(p['gut_fractions'])
    gscale=1-wg.sum()+np.dot(wg,ag)
    fg=1/(1+(1/p['fg']-1)*gscale)
    change=model.c.fa*(fg-p['fg'])*model.c.ka_h*y[t4.GUT]
    dy[t4.LIVER]+=change;dy[t4.FEC]-=change
    return dy


def generalized_static(model,p,ah,ag):
    """Exact integrated oral AUC in the constant-activity linear Task 4 limit."""
    w=np.asarray(p['hepatic_intrinsic_fractions']);wg=np.asarray(p['gut_fractions'])
    h=1-w.sum()+np.dot(w,ah);g=1-wg.sum()+np.dot(wg,ag)
    fg=1/(1+(1/p['fg']-1)*g)
    q=model.c.flows_l_h['liver']
    x=model.c.fu_p/model.c.rb*model.clint*h
    fh=q/(q+x);clh=q*x/(q+x)*model.c.rb
    return p['dose_mg']*p['fa']*fg*fh/(clh+model.clrenal)


def simulate_victim(p,exposure=None,start=312.,rtol=2e-7,atol=1e-10):
    model=pbpk_model(p)
    y=np.zeros(11);y[t4.GUT]=p['dose_mg']
    end=END
    edges=sorted(set([start]+[float(x) for x in np.arange(0,LAST_DOSE+1,12.) if start<x<end]+[end]))
    segments=[];mass_error=0.
    for a,b in zip(edges[:-1],edges[1:]):
        def rhs(t,s):
            ah,ag=(np.ones(3),np.ones(3)) if exposure is None else exposure.factors(t)
            return victim_rhs(model,p,s,ah,ag)
        sol=solve_ivp(rhs,(a,b),y,method='BDF',rtol=rtol,atol=atol,max_step=2.,dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        mass_error=max(mass_error,float(np.max(np.abs(sol.y[:10].sum(axis=0)-p['dose_mg']))))
        if np.min(sol.y[:7]) < -1e-7: raise RuntimeError('negative probe amount')
        y=sol.y[:,-1].copy();segments.append(sol)
    if mass_error>2e-6*p['dose_mg']: raise RuntimeError('probe mass conservation failed')
    # The label is finite-horizon AUC; residual fraction is reported, no unjustified tail extrapolation.
    time=np.linspace(start,min(start+96,END),769)
    ends=np.array([s.t[-1] for s in segments]);idx=np.searchsorted(ends,time,side='right').clip(max=len(segments)-1)
    cp=np.array([model.cp(segments[k].sol(t)) for t,k in zip(time,idx)])
    return dict(auc_mg_h_l=float(y[t4.AUC]),cmax_mg_l=float(cp.max()),
                residual_fraction=float(y[:7].sum()/p['dose_mg']),mass_error_mg=mass_error,
                times_h=time-start,cp_mg_l=cp,model=model,horizon_h=END-start)


def recovery(exposure,compartment,index,target=.9):
    state_index=(11 if compartment=='hepatic' else 14)+index
    grid=np.linspace(LAST_DOSE,END,2401)
    values=exposure.at(grid)[state_index]
    below=np.flatnonzero(values<target)
    if len(below)==0: return 0.
    k=int(below[-1])
    if k==len(grid)-1: return None
    # Last bracketed crossing; all later 0.1 h samples must remain above target.
    crossing=brentq(lambda t:exposure.at(t)[state_index]-target,grid[k],grid[k+1],xtol=1e-9)
    return float((crossing-LAST_DOSE)/24.)


def save_csv(path,rows):
    if not rows: raise ValueError('empty output table')
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def save_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def parameter_ledger(out,ps):
    rows=[]
    for d in ps:
        for key,value in d.items():
            if isinstance(value,(int,float,list)):
                for j,v in enumerate(value if isinstance(value,list) else [value]):
                    rows.append(dict(entity=d['name'],parameter=key,isoform=ISO[j] if isinstance(value,list) else '',value=v,
                                     evidence='synthetic_scenario_assumption',source='analyst_defined_not_extracted',
                                     limitation='Drug name does not validate this numeric value'))
    for p in probes():
        for key,v in p.items():
            if isinstance(v,(int,float,list)):
                rows.append(dict(entity=p['name'],parameter=key,isoform='',value=json.dumps(v),evidence='user_specified' if key=='dose_mg' else 'synthetic_scenario_assumption',source='user_task_prompt' if key=='dose_mg' else 'analyst_defined_not_extracted',limitation='No clinical PK calibration'))
    for j,n in enumerate(ISO):
        for key,v in [('hepatic_kdeg_h',KDEG[j]),('gut_kdeg_h',KDEG_GUT[j])]:
            rows.append(dict(entity='enzyme',parameter=key,isoform=n,value=v,evidence='user_example_assumption' if j==0 and key=='hepatic_kdeg_h' else 'synthetic_scenario_assumption',source='user_task_prompt' if j==0 and key=='hepatic_kdeg_h' else 'analyst_defined_not_extracted',limitation='Not individually measured turnover'))
    save_csv(out/'parameter_provenance.csv',rows)


def literature_reference(out):
    """Separate published MODEL parameter response, never plugged into clinical PK."""
    source='https://pmc.ncbi.nlm.nih.gov/articles/PMC2812061/'
    parameters=[dict(compound='Clarithromycin',parameter='KI',value=5.3,unit='uM',site='liver_and_gut',table='2'),
                dict(compound='Clarithromycin',parameter='kinact',value=.4,unit='h^-1',site='liver',table='2'),
                dict(compound='Clarithromycin',parameter='kinact',value=4.,unit='h^-1',site='gut',table='2')]
    for row in parameters:
        row.update(source_url=source,doi='10.1124/dmd.109.028746',evidence='published_in_vivo_model_estimate_not_measured_in_vitro_kinetics',limitation='Point estimates only; binding convention not reinterpreted as unbound; not transplanted into Task5 PBPK')
    save_csv(out/'literature_parameter_reference.csv',parameters)
    rows=[]
    for site,kin in [('liver',.4),('gut',4.)]:
        for concentration in (.1,1.,10.):
            rate=float(kobs(concentration,kin,5.3))
            for time in np.linspace(0,4,41):
                rows.append(dict(source_doi='10.1124/dmd.109.028746',site=site,assumed_constant_I_uM=concentration,time_h=float(time),computed_kobs_h=rate,computed_survival_no_resynthesis=float(np.exp(-rate*time)),evidence='calculated_transfer_function_of_published_model_parameters_not_assay_data'))
    save_csv(out/'literature_parameter_response.csv',rows)
    save_json(out/'sources.json',{'access_date':'2026-09-20','sources':[
        {'id':'FDA_M12_2024','url':SOURCE_M12,'supports':'basic screens and inhibition-only mechanistic static model; not numeric compound inputs'},
        {'id':'FDA_CYP_TABLE','url':SOURCE_TABLE,'supports':'qualitative inhibitor/substrate examples; no numeric Ki or PK values inferred'},
        {'id':'Quinney2010','url':source,'doi':'10.1124/dmd.109.028746','supports':'Table 2 published clarithromycin model estimates; separate transfer-function calculation'},
        {'id':'Rock2014','url':'https://pubmed.ncbi.nlm.nih.gov/25274602/','doi':'10.1124/mol.114.094862','supports':'experimental ritonavir-derived CYP3A4 apoprotein Lys257 adduct'},
        {'id':'Lin2013','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3781371/','supports':'ritonavir heme alteration in a reconstituted enzyme system'},
        {'id':'Bergamottin2012','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3336797/','supports':'furan bioactivation and CYP3A4 covalent residue modification'}]})


def figures(out,exposures,curves,ps):
    figdir=out/'figures_task5';figdir.mkdir()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
    inc_rows=[];fits=[]
    for ax,name in zip(axes,['Clarithromycin','Ritonavir']):
        d=next(x for x in ps if x['name']==name)
        concentrations=np.array([.02,.1,.5,2.,10.]);tt=np.linspace(0,4,101)
        rates=[]
        for c in concentrations:
            rate=float(kobs(c,d['kinact_h'][0],d['KI_uM'][0]));yy=np.exp(-rate*tt)
            rates.append(-np.polyfit(tt,np.log(yy),1)[0])
            ax.plot(tt,yy,label=f'{c:g} µM')
            inc_rows += [dict(compound=name,concentration_uM=c,time_h=t,relative_activity=a,evidence='noise_free_synthetic_assay') for t,a in zip(tt,yy)]
        popt,_=curve_fit(lambda x,ki,KI:ki*x/(KI+x),concentrations,rates,p0=[.3,1.],bounds=(0,np.inf),xtol=1e-12,ftol=1e-12,gtol=1e-12)
        fits.append(dict(compound=name,true_kinact_h=d['kinact_h'][0],fit_kinact_h=float(popt[0]),true_KI_uM=d['KI_uM'][0],fit_KI_uM=float(popt[1]),evidence='noise_free_parameter_recovery_not_measurement'))
        ax.set(title=name+' scenario',xlabel='Preincubation (h)',ylabel='Remaining activity');ax.legend(fontsize=8)
    fig.suptitle('Constructed MBI assay: no enzyme synthesis and no experimental observations')
    fig.savefig(figdir/'fig1_cyp_inactivation_curves.png');plt.close(fig)
    save_csv(out/'synthetic_preincubation.csv',inc_rows);save_csv(out/'synthetic_kinetic_recovery.csv',fits)

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for ax,p in zip(axes,probes()):
        for name in ('Control','Clarithromycin','Ritonavir','Quinidine'):
            c=curves[(name,p['name'],312.)]
            mask=c['times_h']<=48
            ax.plot(c['times_h'][mask],c['cp_mg_l'][mask]*1000,label=name)
        ax.set(title=p['name']+f" {p['dose_mg']:g} mg oral",xlabel='Hours after day-14 probe',ylabel='Total plasma (ng/mL)');ax.legend(fontsize=8)
    fig.suptitle('Synthetic Task 4-coupled PK; last inhibitor dose 12 h after probe')
    fig.savefig(figdir/'fig2_dynamic_ddi_midazolam_pk.png');plt.close(fig)

    # Fixed I=1uM, KI=1uM and kinact derived from y. This is a conditional
    # decision diagram, not a universal safety boundary or clinical FDA matrix.
    x=np.geomspace(1e-4,100,220);y=np.geomspace(1e-4,10,220);xx,yy=np.meshgrid(x,y)
    r1=1+xx;r2=1+(yy*5/(1+5))/.019
    ah=(.019/(.019+yy/2))/(1+xx)
    aucr=1/(.9*ah+.1) # Fg=1, negligible extrahepatic clearance
    cat=np.where((r1<1.02)&(r2<1.25),0,np.where(aucr>=5,2,1))
    fig,ax=plt.subplots(figsize=(9,6))
    im=ax.pcolormesh(x,y,cat,cmap=ListedColormap(['#73b69d','#efc85b','#d96868']),norm=BoundaryNorm([-.5,.5,1.5,2.5],3),shading='auto')
    ax.set(xscale='log',yscale='log',xlabel='Cmax,u / Ki (dimensionless)',ylabel='kinact / KI (h⁻¹ µM⁻¹)',title='Conditional synthetic screen: Cmax,u = KI = 1 µM, fm = 0.9, Fg = 1')
    ax.axvline(.02,color='black',ls='--',lw=.8)
    cb=fig.colorbar(im,ax=ax,ticks=[0,1,2]);cb.ax.set_yticklabels(['Below basic cutoffs*','Further evaluation','Model AUCR ≥ 5'])
    fig.subplots_adjust(left=.11,right=.77,bottom=.19,top=.88)
    fig.text(.10,.025,'*Constructed hepatic scenario only; no clinical safety or contraindication claim.',fontsize=9)
    fig.savefig(figdir/'fig3_fda_ddi_risk_heatmap.png');plt.close(fig)
    save_csv(out/'conditional_risk_grid.csv',[dict(imax_u_over_ki=float(a),kinact_over_KI_h_uM=float(b),r1=float(c),r2=float(d),illustrative_aucr=float(e),category=int(f)) for a,b,c,d,e,f in zip(xx.ravel(),yy.ravel(),r1.ravel(),r2.ravel(),aucr.ravel(),cat.ravel())])

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    tt=np.linspace(LAST_DOSE,END,1001)
    for name in ('Clarithromycin','Ritonavir'):
        e=exposures[name];v=e.at(tt)
        axes[0].plot((tt-LAST_DOSE)/24,v[11],label=name+' hepatic')
        axes[0].plot((tt-LAST_DOSE)/24,v[14],ls='--',label=name+' gut')
        ih,ig,ip=e.concentrations(v)
        axes[1].semilogy((tt-LAST_DOSE)/24,np.maximum(ih,1e-9),label=name+' liver unbound')
    axes[0].axhline(.9,color='grey',ls=':');axes[0].set(ylabel='Active CYP3A4 / baseline',ylim=(0,1.05),xlabel='Days since LAST dose');axes[0].legend(fontsize=8)
    axes[1].set(ylabel='Residual inhibitor (µM)',xlabel='Days since LAST dose');axes[1].legend(fontsize=8)
    fig.suptitle('Recovery retains inhibitor PK and ongoing inactivation after the final dose')
    fig.savefig(figdir/'fig4_cyp3a4_resynthesis_timeline.png');plt.close(fig)


def execute(out,ps=None):
    out.mkdir(parents=True,exist_ok=False)
    ps=panel() if ps is None else ps
    if len(ps)!=6 or {d['name'] for d in ps}!={d['name'] for d in panel()}:
        raise ValueError('panel requires the six named scenario compounds, each exactly once')
    for d in ps: validate_compound(d)
    parameter_ledger(out,ps)
    literature_reference(out)
    save_json(out/'parameters_used.json',{'compound_panel':ps,'probe_scenarios':probes(),'hepatic_kdeg_h':KDEG.tolist(),'gut_kdeg_h':KDEG_GUT.tolist(),'dose_times_h':np.arange(0,LAST_DOSE+1,12.).tolist(),'last_dose_h':LAST_DOSE,'end_h':END,'numeric_evidence':'synthetic_except_user_probe_doses','all_tissue_kp_overrides':1.,'rtol':2e-7,'atol':1e-9,'full_task4_configs':{d['name']:asdict(pbpk_model(d).c) for d in ps+probes()}})
    exposures={};curves={};rows=[];wash=[];screen=[];trajectory=[]
    controls={}
    for p in probes():
        for start in (0.,312.):
            controls[(p['name'],start)]=simulate_victim(p,start=start)
            curves[('Control',p['name'],start)]=controls[(p['name'],start)]
    for d in ps:
        print('Computing '+d['name'],flush=True)
        ex=Exposure(d);exposures[d['name']]=ex
        tt=np.arange(0,END+.25,.25);v=ex.at(tt);ih,ig,ip=ex.concentrations(v)
        for k,t in enumerate(tt):
            trajectory.append(dict(compound=d['name'],time_h=t,plasma_unbound_uM=ip[k],liver_unbound_uM=ih[k],gut_proxy_uM=ig[k],**{f'hep_{n}':v[11+j,k] for j,n in enumerate(ISO)},**{f'gut_{n}':v[14+j,k] for j,n in enumerate(ISO)}))
        cycle=(tt>=312)&(tt<324);cmax=float(ip[cycle].max())
        mean_ih=float(np.mean(ih[cycle]));mean_ig=float(np.mean(ig[cycle]))
        ah=activity(KDEG/(KDEG+kobs(mean_ih,d['kinact_h'],d['KI_uM'])),mean_ih,d['ki_uM'])
        ag=activity(KDEG_GUT/(KDEG_GUT+kobs(mean_ig,d['kinact_h'],d['KI_uM'])),mean_ig,d['ki_uM'])
        for j,n in enumerate(ISO):
            r1,r2=basic_ratios(cmax,d['ki_uM'][j],d['kinact_h'][j],d['KI_uM'][j],KDEG[j])
            gut_ratio=d['dose_mg']*1000/d['mw_g_mol']/.250/d['ki_uM'][j]
            screen.append(dict(compound=d['name'],isoform=n,day14_cmax_uM=cmax,r1=r1,r2_m12_5x=r2,gut_dose250_over_ki=gut_ratio if j==0 else '',hepatic_further_evaluation=r1>=1.02 or r2>=1.25,gut_further_evaluation=gut_ratio>=10 if j==0 else 'not_assessed',evidence='synthetic_screen_not_clinical_risk'))
            for comp in ('hepatic','gut'):
                value=recovery(ex,comp,j)
                wash.append(dict(compound=d['name'],compartment=comp,isoform=n,recovery_90_days_after_last_dose=value,censored_at_10days=value is None,kdeg_half_life_h=math.log(2)/(KDEG[j] if comp=='hepatic' else KDEG_GUT[j]),definition='last_upward_crossing_with_all_later_0.1h_samples_above_90pct'))
        for p in probes():
            model=pbpk_model(p)
            static=generalized_static(model,p,ah,ag)/generalized_static(model,p,np.ones(3),np.ones(3))
            j=0 if p['name']=='Midazolam' else 1
            official=official_static_aucr(ah[j],ag[j] if j==0 else 1.,.9,p['fg'])
            for start in (0.,312.):
                c=simulate_victim(p,ex,start);curves[(d['name'],p['name'],start)]=c
                control=controls[(p['name'],start)];ratio=c['auc_mg_h_l']/control['auc_mg_h_l']
                rows.append(dict(compound=d['name'],probe=p['name'],probe_time_h=start,auc_horizon_h=c['horizon_h'],auc_control_mg_h_l=control['auc_mg_h_l'],auc_inhibited_mg_h_l=c['auc_mg_h_l'],dynamic_aucr=ratio,illustrative_threshold_band=severity(ratio),generalized_constant_mean_static_aucr=static,m12_static_negligible_extrahepatic_aucr=official,probe_residual_fraction=c['residual_fraction'],mass_error_mg=c['mass_error_mg']))
    print('Computing numerical sensitivity and turnover scenarios',flush=True)
    # Independent tighter solve checks a genuinely coupled TDI trajectory and AUC.
    clari=next(d for d in ps if d['name']=='Clarithromycin')
    fine=Exposure(clari,rtol=2e-9,atol=1e-11)
    fp=simulate_victim(probes()[0],fine,rtol=2e-9,atol=1e-12)
    base=curves[('Clarithromycin','Midazolam',312.)]
    convergence=abs(fp['auc_mg_h_l']/base['auc_mg_h_l']-1)
    sensitivity=[]
    for factor in (.5,1.,2.):
        ex=exposures['Clarithromycin'] if factor==1 else Exposure(clari,kdeg_factor=factor)
        c=simulate_victim(probes()[0],ex)
        sensitivity.append(dict(kdeg_factor=factor,midazolam_day14_aucr=c['auc_mg_h_l']/controls[('Midazolam',312.)]['auc_mg_h_l'],hepatic_cyp3a4_recovery_days=recovery(ex,'hepatic',0)))
    save_csv(out/'dynamic_static_comparison.csv',rows);save_csv(out/'washout_recovery.csv',wash)
    save_csv(out/'basic_regulatory_screen.csv',screen);save_csv(out/'enzyme_and_inhibitor_trajectories.csv',trajectory)
    save_csv(out/'turnover_sensitivity.csv',sensitivity)
    pkrows=[]
    for (name,probe,start),c in curves.items():
        pkrows += [dict(compound=name,probe=probe,probe_time_h=start,time_after_probe_h=float(t),plasma_mg_l=float(x)) for t,x in zip(c['times_h'],c['cp_mg_l'])]
    save_csv(out/'probe_pk_curves.csv',pkrows)
    figures(out,exposures,curves,ps)
    summary={'evidence':'executed_synthetic_mechanistic_scenarios_not_clinical_predictions','compounds':6,'isoforms':list(ISO),'dose_count':28,'main_probe_comparisons':len(rows),'washout_entries':len(wash),'main_perpetrator_mass_error_mg':max(e.mass_error for e in exposures.values()),'max_probe_mass_error_mg':max(r['mass_error_mg'] for r in rows),'max_probe_residual_fraction':max(r['probe_residual_fraction'] for r in rows),'tighter_solver_clarithromycin_midazolam_auc_relative_change':convergence,'task4_derivative_imported':True,'task4_sha256':hashlib.sha256(TASK4_PATH.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},'limits':['No measured numeric compound input or fitted clinical PK','No self-inhibition, induction, inhibitory metabolites or transporter kinetics','Gut concentration is a flow-based proxy, not a new validated gut-wall model','Two probes are independent hypothetical administrations, not a validated cocktail','Day 14 is not automatically claimed to be steady state','AUC is finite horizon, with residual drug fractions reported','Recovery is enzyme abundance, not combined reversible-inhibition clearance recovery','Basic screens are not AUCR and no below-cutoff result is a clinical safety claim']}
    if convergence>1e-4: raise RuntimeError('solver convergence gate failed')
    save_json(out/'results_summary.json',summary)
    files={p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
    save_json(out/'manifest.json',{'algorithm':'sha256','files':files})
    print(json.dumps(summary,indent=2),flush=True)


def self_test():
    import unittest
    suite=unittest.defaultTestLoader.discover(str(HERE.parent/'tests'),pattern='test_task5.py')
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--config',type=Path,help='JSON compound panel, matching --write-example schema')
    parser.add_argument('--write-example',type=Path,help='Write editable synthetic panel to a new file')
    args=parser.parse_args()
    if args.self_test: return self_test()
    if args.write_example is not None:
        with args.write_example.open('x',encoding='utf-8') as f:
            json.dump({'compound_panel':panel()},f,indent=2,ensure_ascii=False,allow_nan=False)
            f.write('\n')
        return 0
    if args.out is None: parser.error('--out NEW_DIRECTORY is required')
    ps=None
    if args.config is not None:
        obj=json.loads(args.config.read_text(encoding='utf-8'))
        if set(obj)!={'compound_panel'}: parser.error('config must contain only compound_panel')
        ps=obj['compound_panel']
    execute(args.out.resolve(),ps);return 0


if __name__=='__main__':
    raise SystemExit(main())
