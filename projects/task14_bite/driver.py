"""Finite-receptor BiTE equilibrium, contact-volume conversion and cell lysis."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE='fig14_bite_synapse_crosslinking_curve.png'
NA=6.02214076e23
CONFIG={'KD_CD3_nM':100.,'KD_TAA_nM':1.,'cooperativity':1.,'accessible_CD3_per_effector':5000.,'accessible_TAA_per_target':10000.,'effective_contact_volume_L_per_initial_target':1e-12,'initial_target_cells':1e5,'effector_target_ratio':1.,'lysis_max_per_h':.12,'synapse_half_effect_complexes_per_cell':500.,'hill_exponent':2.,'horizon_h':48.,'evidence':'All biological parameters are illustrative assumptions; contact volume is a coarse-graining convention, not measured cleft geometry.'}

def equilibrium(total,et,tt,ke=100.,kt=1.,alpha=1.):
    """Concentrations in nM. Species returned E,T,B,EB,TB,EBT."""
    if min(total,et,tt)<0 or min(ke,kt)<=0 or alpha<0:raise ValueError('Nonnegative totals and positive KDs are required')
    if total==0:return np.array([et,tt,0,0,0,0],dtype=float)
    def at(b):
        def values(e):
            t=tt/(1+b/kt+alpha*e*b/(ke*kt));x=alpha*e*t*b/(ke*kt)
            return t,x
        e=brentq(lambda e:e*(1+b/ke)+values(e)[1]-et,0,et,xtol=1e-13) if et else 0.
        t,x=values(e)
        return np.array([e,t,b,e*b/ke,t*b/kt,x])
    b=brentq(lambda b:at(b)[[2,3,4,5]].sum()-total,0,total,xtol=1e-13)
    return at(b)

def totals(fraction=1.,ratio=1.,expression=1.):
    factor=1e9/(NA*CONFIG['effective_contact_volume_L_per_initial_target'])
    return CONFIG['accessible_CD3_per_effector']*ratio*factor,CONFIG['accessible_TAA_per_target']*expression*fraction*factor

def synapses(dose,fraction=1.,ratio=1.,expression=1.,ke=100.,kt=1.,alpha=1.):
    e,t=totals(fraction,ratio,expression);x=equilibrium(dose,e,t,ke,kt,alpha)[5]
    return float(x*1e-9*NA*CONFIG['effective_contact_volume_L_per_initial_target']/max(fraction,1e-15))

def simulate(dose=10.,ratio=1.,expression=1.,ke=100.,kt=1.,alpha=1.,rtol=2e-7):
    c=CONFIG
    def rhs(t,y):
        s=synapses(dose,max(y[0],0.),ratio,expression,ke,kt,alpha)
        frac=(s/c['synapse_half_effect_complexes_per_cell'])**c['hill_exponent']
        flux=c['lysis_max_per_h']*frac/(1+frac)*max(y[0],0.)
        return [-flux,flux]
    sol=solve_ivp(rhs,(0,48),[1.,0.],t_eval=np.linspace(0,48,193),rtol=rtol,atol=rtol*1e-3)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def csv_write(path,header,rows):
    with path.open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)

def run(out:Path)->dict:
    out=Path(out)
    if out.exists():raise FileExistsError('Use a new output directory')
    (out/'figures').mkdir(parents=True)
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2)+'\n',encoding='utf-8')
    rows=[];grid=[];maxbal=0.;doses=np.logspace(-4,5,109)
    e,t=totals()
    for ke in [10.,100.,1000.]:
        for kt in [.1,1.,10.]:
            for dose in doses:
                eq=equilibrium(dose,e,t,ke,kt);s=synapses(dose,ke=ke,kt=kt)
                error=max(abs(eq[[0,3,5]].sum()-e),abs(eq[[1,4,5]].sum()-t),abs(eq[[2,3,4,5]].sum()-dose));maxbal=max(maxbal,error)
                rows.append([dose,ke,kt,*eq,s,error])
    for ratio in [.1,1.,5.]:
        for expression in [.1,1.,3.]:
            for dose in np.r_[0.,np.logspace(-3,4,29)]:
                sol=simulate(dose,ratio,expression)
                grid.append([dose,ratio,expression,sol.y[0,-1],sol.y[1,-1],synapses(dose,ratio=ratio,expression=expression),float(np.max(np.abs(sol.y.sum(axis=0)-1)))])
    csv_write(out/'binding_grid.csv',['total_BiTE_nM','KD_CD3_nM','KD_TAA_nM','free_CD3_nM','free_TAA_nM','free_BiTE_nM','CD3_BiTE_nM','TAA_BiTE_nM','trimer_nM','trimer_complexes_per_target','max_mass_balance_error_nM'],rows)
    csv_write(out/'lysis_grid.csv',['total_BiTE_nM','effector_target_ratio','TAA_expression_multiplier','surviving_fraction_48h','lysed_fraction_48h','initial_trimer_complexes_per_target','cell_ledger_error'],grid)
    trajectories=[]
    for dose in [0.,.01,1.,10.,1000.,1e4]:
        s=simulate(dose)
        for i,time in enumerate(s.t):trajectories.append([dose,time,*s.y[:,i],synapses(dose,s.y[0,i])])
    csv_write(out/'trajectories.csv',['total_BiTE_nM','time_h','surviving_fraction','lysed_fraction','trimer_complexes_per_surviving_target'],trajectories)
    base=simulate();fine=simulate(rtol=2e-10);zero=simulate(0.);noeff=simulate(ratio=0.);noalpha=simulate(alpha=0.)
    checks={'binding_grid_max_balance_error_nM':maxbal,'cell_ledger_max_error':max(r[-1] for r in grid),'refinement_max_fraction_difference':float(np.max(np.abs(base.y-fine.y))),'zero_BiTE_lysis':float(zero.y[1,-1]),'zero_effector_lysis':float(noeff.y[1,-1]),'zero_cooperativity_lysis':float(noalpha.y[1,-1]),'minimum_cell_fraction':float(base.y.min())}
    checks['passed']=bool(maxbal<1e-6 and checks['cell_ledger_max_error']<1e-8 and checks['refinement_max_fraction_difference']<1e-5 and max(checks['zero_BiTE_lysis'],checks['zero_effector_lysis'],checks['zero_cooperativity_lysis'])==0 and checks['minimum_cell_fraction']>=0)
    nominal=[r for r in rows if r[1]==100. and r[2]==1.];peak=max(nominal,key=lambda r:r[9]);high=nominal[-1]
    summary={'task':14,'evidence':'executed closed receptor equilibrium with clamped BiTE exposure and phenomenological lysis','binding_grid_points':len(rows),'lysis_ode_runs':len(grid),'nominal_peak_BiTE_nM_in_grid':peak[0],'peak_trimer_complexes_per_target':peak[9],'high_dose_100uM_fraction_of_peak':high[9]/peak[9],'hook_detected_in_nominal_grid':bool(high[9]<.5*peak[9]),'10nM_surviving_fraction_48h':float(base.y[0,-1]),'limitation':'The hook follows this reversible bispecific equilibrium under fixed totals; not universal across TCE architectures or clinical exposure windows. No cytokine toxicity, serial killing kinetics, exhaustion or calibrated clinical response.'}
    for name,obj in [('summary.json',summary),('verification.json',checks)]: (out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False});fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ke in [10.,100.,1000.]:
        r=np.array([r for r in rows if r[1]==ke and r[2]==1.]);ax[0,0].semilogx(r[:,0],r[:,9],label=f'CD3 KD {ke:g} nM')
    ax[0,0].set(xlabel='Total BiTE (nM)',ylabel='Trimer complexes / target cell',title='A  Finite binding pools and asymmetric affinity');ax[0,0].legend(fontsize=8)
    for ratio in [.1,1.,5.]:
        r=np.array([r for r in grid if r[1]==ratio and r[2]==1. and r[0]>0]);ax[0,1].semilogx(r[:,0],1-r[:,3],label=f'E:T = {ratio:g}')
    ax[0,1].set(xlabel='Clamped total BiTE (nM)',ylabel='Lysed fraction at 48 h',title='B  Dose × effector availability');ax[0,1].legend(fontsize=8)
    for dose in [.01,1.,10.,1000.,1e4]:
        s=simulate(dose);ax[1,0].plot(s.t,s.y[0],label=f'{dose:g} nM')
    ax[1,0].set(xlabel='Time (h)',ylabel='Target surviving fraction',title='C  Nonlinear lysis with target depletion');ax[1,0].legend(fontsize=8)
    r=np.array(nominal);ax[1,1].semilogx(r[:,0],r[:,5],label='Free BiTE');ax[1,1].semilogx(r[:,0],r[:,6],label='CD3–BiTE');ax[1,1].semilogx(r[:,0],r[:,7],label='TAA–BiTE');ax[1,1].semilogx(r[:,0],r[:,8],label='Trimer');ax[1,1].set(yscale='log',xlabel='Total BiTE (nM)',ylabel='Species (nM)',title='D  Binary occupation competes with crosslinking');ax[1,1].legend(fontsize=8)
    fig.suptitle('Task 14 | BiTE crosslinking and lysis: calculated hook, not a universal clinical law',fontweight='bold');fig.savefig(out/'figures'/FIGURE,dpi=300);plt.close(fig)
    if not checks['passed']:raise AssertionError(checks)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);print(json.dumps(run(p.parse_args().out),indent=2))
