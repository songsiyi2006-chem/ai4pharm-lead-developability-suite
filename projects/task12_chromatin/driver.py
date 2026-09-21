"""Competitive free/chromatin target degradation with explicit mass action."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE='fig12_epigenetic_chromatin_depletion.png'
CONFIG={'target_nM':100.,'e3_nM':20.,'degrader_nM':100.,'free_target_initial_nM':30.,'chromatin_target_initial_nM':70.,'KD_E_nM':50.,'KD_T_nM':100.,'kon_per_nM_h':.01,'alpha_free':10.,'chromatin_penalty_kcal_per_mol':2.,'RT_kcal_per_mol':.61633,'chromatin_capture_per_h':.14,'chromatin_release_per_h':.06,'target_turnover_per_h':.025,'ternary_degradation_per_h':.5,'horizon_h':48.,'evidence':'All numerical biological parameters are hypothetical; no fitted BRD4/EZH2/HDAC dataset.'}
NAMES=['F','C','P','E','PF','PC','EP','EFP','ECP','target_degraded','target_synthesized']

def equilibrium(p_total, f_total=30., c_total=70., penalty=2.):
    """Exact binding-only mass balances for two competing target pools (nM)."""
    ke,kt,et=50.,100.,20.
    af=10.; ac=af*np.exp(-penalty/CONFIG['RT_kcal_per_mol'])
    def at(p):
        def values(e):
            f=f_total/(1+p/kt+af*e*p/(ke*kt)); c=c_total/(1+p/kt+ac*e*p/(ke*kt))
            return f,c,af*e*p*f/(ke*kt),ac*e*p*c/(ke*kt)
        def residual(e):
            f,c,tf,tc=values(e)
            return e*(1+p/ke)+tf+tc-et
        e=brentq(residual,0,et,xtol=1e-12)
        f,c,tf,tc=values(e)
        return np.array([f,c,p,e,p*f/kt,p*c/kt,p*e/ke,tf,tc])
    if p_total==0: return at(0.)
    p=brentq(lambda p:at(p)[[2,4,5,6,7,8]].sum()-p_total,0,p_total,xtol=1e-12)
    return at(p)

def rhs(t,y,penalty=2.,kdeg=.5):
    c=CONFIG; d=np.zeros(11)
    def flux(rate,reactants,products):
        for i in reactants: d[i]-=rate
        for i in products: d[i]+=rate
    k=c['kon_per_nM_h']; ke=c['KD_E_nM']; kt=c['KD_T_nM']; af=c['alpha_free']; ac=af*np.exp(-penalty/c['RT_kcal_per_mol'])
    for target,binary,ternary,alpha in [(0,4,7,af),(1,5,8,ac)]:
        flux(k*y[2]*y[target]-k*kt*y[binary],[2,target],[binary])
        flux(k*y[6]*y[target]-k*kt/alpha*y[ternary],[6,target],[ternary])
        flux(k*y[binary]*y[3]-k*ke/alpha*y[ternary],[binary,3],[ternary])
        flux(kdeg*y[ternary],[ternary],[2,3,9])
    flux(k*y[2]*y[3]-k*ke*y[6],[2,3],[6])
    flux(c['chromatin_capture_per_h']*y[0]-c['chromatin_release_per_h']*y[1],[0],[1])
    # Every target-containing species turns over and releases surviving E3/degrader.
    for i,products in [(0,[9]),(1,[9]),(4,[2,9]),(5,[2,9]),(7,[2,3,9]),(8,[2,3,9])]:
        flux(c['target_turnover_per_h']*y[i],[i],products)
    synthesis=c['target_turnover_per_h']*c['target_nM']
    d[0]+=synthesis; d[10]+=synthesis
    return d

def simulate(dose=100.,penalty=2.,kdeg=.5,rtol=1e-8):
    y0=np.array([30.,70.,dose,20.,0,0,0,0,0,0,0])
    # All constitutive synthesis enters the free pool. The initial 30:70 split
    # relaxes even without drug; treated pools require matched untreated controls.
    sol=solve_ivp(lambda t,y:rhs(t,y,penalty,kdeg),(0,48),y0,t_eval=np.linspace(0,48,193),method='BDF',rtol=rtol,atol=rtol*1e-3)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

def pools(y): return y[[0,4,7]].sum(axis=0),y[[1,5,8]].sum(axis=0)
def balances(y,dose):
    return [float(np.max(np.abs(y[[2,4,5,6,7,8]].sum(axis=0)-dose))),float(np.max(np.abs(y[[3,6,7,8]].sum(axis=0)-20))),float(np.max(np.abs(y[[0,1,4,5,7,8]].sum(axis=0)+y[9]-y[10]-100)))]
def csv_write(path,head,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(rows)

def run(out:Path)->dict:
    out=Path(out)
    if out.exists():raise FileExistsError('Use a new output directory')
    (out/'figures').mkdir(parents=True)
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2)+'\n',encoding='utf-8')
    base=simulate(); zero=simulate(0); disabled=simulate(kdeg=0); fine=simulate(rtol=1e-10)
    trajectories=[]
    for label,s in [('baseline',base),('no_degrader',zero),('no_induced_degradation',disabled),('no_penalty',simulate(penalty=0))]:
        for i,t in enumerate(s.t):trajectories.append([label,t,*s.y[:,i]])
    csv_write(out/'trajectories.csv',['scenario','time_h',*[n+'_nM' for n in NAMES]],trajectories)
    rows=[]; eqrows=[]; griderr=0.
    for penalty in [0.,1.,2.,4.]:
        for dose in np.r_[0.,np.logspace(-2,4,37)]:
            s=simulate(dose,penalty); f,c=pools(s.y); griderr=max(griderr,*balances(s.y,dose))
            rows.append([dose,penalty,f[-1],c[-1],(f+c)[-1],float(s.y[9,-1]),float(s.y[10,-1])])
            eq=equilibrium(dose,penalty=penalty);eqrows.append([dose,penalty,*eq])
    csv_write(out/'dose_penalty_grid.csv',['degrader_total_nM','penalty_kcal_mol','free_pool_48h_nM','chromatin_pool_48h_nM','total_48h_nM','degraded_nM','synthesized_nM'],rows)
    csv_write(out/'binding_equilibria.csv',['degrader_total_nM','penalty_kcal_mol',*[n+'_nM' for n in NAMES[:9]]],eqrows)
    b=balances(base.y,100.); f,c=pools(base.y);fz,cz=pools(zero.y)
    checks={'degrader_balance_error_nM':b[0],'e3_balance_error_nM':b[1],'target_ledger_error_nM':b[2],'grid_max_balance_error_nM':griderr,'refinement_max_abs_difference_nM':float(np.max(np.abs(base.y-fine.y))),'no_degrader_total_deviation_nM':float(np.max(np.abs(fz+cz-100))),'no_induced_degradation_total_deviation_nM':float(np.max(np.abs(sum(pools(disabled.y))-100))),'minimum_species_nM':float(min(base.y.min(),zero.y.min()))}
    checks['passed']=bool(griderr<1e-5 and checks['refinement_max_abs_difference_nM']<1e-4 and checks['no_degrader_total_deviation_nM']<1e-5 and checks['no_induced_degradation_total_deviation_nM']<1e-5 and checks['minimum_species_nM']>-1e-7)
    summary={'task':12,'evidence':'executed competitive mass-action hypothesis; no measured chromatin occlusion energies','grid_runs':len(rows),'free_pool_48h_nM':float(f[-1]),'chromatin_pool_48h_nM':float(c[-1]),'total_remaining_fraction_48h':float((f[-1]+c[-1])/100),'no_drug_free_pool_48h_nM':float(fz[-1]),'no_drug_chromatin_pool_48h_nM':float(cz[-1]),'alpha_chromatin':float(10*np.exp(-2/CONFIG['RT_kcal_per_mol'])),'limitation':'No nucleosome structure or measured binding/degradation rates; energy penalty is an exposed assumption. Compare treated pools to matched untreated trajectories.'}
    for name,obj in [('summary.json',summary),('verification.json',checks)]: (out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    ax[0,0].plot(base.t,f,label='Free pool');ax[0,0].plot(base.t,c,label='Chromatin pool');ax[0,0].plot(zero.t,fz,'--',label='Free: untreated');ax[0,0].plot(zero.t,cz,'--',label='Chromatin: untreated');ax[0,0].set(xlabel='Time (h)',ylabel='Target (nM)',title='A  Explicit synthesis and two pools');ax[0,0].legend(fontsize=8)
    for penalty in [0.,1.,2.,4.]:
        r=np.array([r for r in rows if r[1]==penalty and r[0]>0]);ax[0,1].semilogx(r[:,0],r[:,4],label=f'Penalty {penalty:g} kcal/mol')
    ax[0,1].set(xlabel='Total degrader (nM)',ylabel='Total target at 48 h (nM)',title='B  Dose × chromatin-access penalty');ax[0,1].legend(fontsize=8)
    ax[1,0].plot(base.t,base.y[7],label='E–P–free target');ax[1,0].plot(base.t,base.y[8],label='E–P–chromatin target');ax[1,0].set(xlabel='Time (h)',ylabel='Ternary complex (nM)',title='C  Competition for finite E3 and degrader');ax[1,0].legend()
    ax[1,1].plot(base.t,f+c,label='Remaining');ax[1,1].plot(base.t,base.y[9],label='Cumulative degraded');ax[1,1].plot(base.t,100+base.y[10],label='Initial + synthesized');ax[1,1].set(xlabel='Time (h)',ylabel='Target-equivalent (nM)',title='D  Open-system target ledger');ax[1,1].legend(fontsize=8)
    fig.suptitle('Task 12 | Chromatin degradation: hypothesis parameters, exact stoichiometry',fontweight='bold');fig.savefig(out/'figures'/FIGURE,dpi=300);plt.close(fig)
    if not checks['passed']:raise AssertionError(checks)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);print(json.dumps(run(p.parse_args().out),indent=2))
