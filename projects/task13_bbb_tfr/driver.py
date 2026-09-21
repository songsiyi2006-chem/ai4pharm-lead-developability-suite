"""Finite antibody/receptor BBB trafficking and conditional affinity tradeoffs."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE='fig13_bbb_tfr_transcytosis_profile.png'
CONFIG={'plasma_volume_L':3.,'brain_extravascular_volume_L':.3,'endosome_volume_L':.001,'brain_vascular_volume_L':.015,'initial_plasma_nM':100.,'receptor_total_nmol':.03,'apical_KD_nM':100.,'kon_per_nM_h':.02,'pH_dependence_log10_KD_per_pH':.7,'endosomal_pH_limit':5.8,'acidification_per_h':1.,'internalization_per_h':.3,'unbound_receptor_internalization_per_h':.02,'receptor_recycling_per_h':.3,'bound_cargo_lysosome_per_h':.08,'free_cargo_lysosome_per_h':.03,'free_cargo_to_brain_per_h':.2,'free_cargo_to_plasma_per_h':.08,'plasma_clearance_per_h':.03,'brain_clearance_per_h':.02,'horizon_h':72.,'evidence':'All numerical system parameters are hypothetical; no human or mouse PK calibration.'}
NAMES=['plasma_antibody','surface_receptor','surface_complex','endosomal_receptor','endosomal_complex','endosomal_free_antibody','brain_free_antibody','systemic_cleared','lysosomal_degraded','brain_cleared']

def simulate(kd=100.,dose=100.,ph_slope=.7,receptor=.03,rtol=2e-8,internalization=.3):
    c=CONFIG
    def rhs(t,y):
        ap,rs,cs,re,ce,ae,ab,_,_,_=y
        ph=c['endosomal_pH_limit']+(7.4-c['endosomal_pH_limit'])*np.exp(-c['acidification_per_h']*t)
        kendo=kd*10**(ph_slope*(7.4-ph))
        bind=c['kon_per_nM_h']*(ap/c['plasma_volume_L']*rs-kd*cs)
        eb=c['kon_per_nM_h']*(ae/c['endosome_volume_L']*re-kendo*ce)
        ci=internalization*cs;ri=c['unbound_receptor_internalization_per_h']*rs;rr=c['receptor_recycling_per_h']*re
        ld=c['bound_cargo_lysosome_per_h']*ce;lf=c['free_cargo_lysosome_per_h']*ae
        tb=c['free_cargo_to_brain_per_h']*ae;tp=c['free_cargo_to_plasma_per_h']*ae
        cp=c['plasma_clearance_per_h']*ap;cb=c['brain_clearance_per_h']*ab
        return [-bind+tp-cp,-bind-ri+rr,bind-ci,ri-rr-eb+ld,ci+eb-ld,-eb-lf-tb-tp,tb-cb,cp,ld+lf,cb]
    # Unbound receptor begins at its recycling/internalization steady state.
    rs=receptor*.3/(.3+.02);re=receptor-rs
    y0=[dose*c['plasma_volume_L'],rs,0,re,0,0,0,0,0,0]
    sol=solve_ivp(rhs,(0,72),y0,t_eval=np.linspace(0,72,289),method='BDF',rtol=rtol,atol=rtol*1e-4)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def exposure(s):
    cp=s.y[0]/CONFIG['plasma_volume_L'];cb=s.y[6]/CONFIG['brain_extravascular_volume_L']
    ap=float(np.trapezoid(cp,s.t));ab=float(np.trapezoid(cb,s.t))
    vascular=cp*CONFIG['brain_vascular_volume_L']/CONFIG['brain_extravascular_volume_L']
    return ap,ab,ab/ap,float(np.trapezoid(cb+vascular,s.t)/ap)

def balances(s,dose=100.,receptor=.03):
    ab=s.y[[0,2,4,5,6,7,8,9]].sum(axis=0)
    rec=s.y[[1,2,3,4]].sum(axis=0)
    return float(np.max(np.abs(ab-dose*3))),float(np.max(np.abs(rec-receptor)))

def csv_write(path,head,rows):
    with path.open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(head);w.writerows(rows)

def run(out:Path)->dict:
    out=Path(out)
    if out.exists():raise FileExistsError('Use a new output directory')
    (out/'figures').mkdir(parents=True)
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2)+'\n',encoding='utf-8')
    rows=[];traces=[];maxerr=0.;maxrec=0.
    kds=np.logspace(-1,4,41)
    for slope,dose,rec,label in [(s,100.,.03,'pH_sensitivity') for s in [0.,.7,1.4]]+[(.7,d,.03,'dose_sensitivity') for d in [1.,10.,1000.]]+[(.7,100.,r,'receptor_sensitivity') for r in [.01,.1]]:
        for kd in kds:
            s=simulate(kd,dose,slope,rec);a,b,kp,app=exposure(s);err,rerr=balances(s,dose,rec);maxerr=max(err,maxerr);maxrec=max(rerr,maxrec)
            rows.append([label,kd,dose,slope,rec,a,b,kp,app,float(s.y[8,-1])])
    for kd in [.1,100.,1e4]:
        s=simulate(kd)
        for i,t in enumerate(s.t):traces.append([kd,t,*s.y[:,i]])
    csv_write(out/'affinity_sensitivity.csv',['scenario','KD_apical_nM','initial_plasma_nM','pH_slope','receptor_nmol','plasma_AUC_nM_h','brain_AUC_nM_h','Kp_extravascular','Kp_with_vascular_contamination','lysosomal_loss_nmol'],rows)
    csv_write(out/'trajectories.csv',['KD_apical_nM','time_h',*[n+'_nmol' for n in NAMES]],traces)
    base=simulate();fine=simulate(rtol=2e-10);zero=simulate(receptor=0.);noentry=simulate(internalization=0.)
    ab,rec=balances(base);ap,abr,kp,app=exposure(base)
    checks={'antibody_ledger_max_error_nmol':ab,'receptor_balance_max_error_nmol':rec,'grid_antibody_max_error_nmol':maxerr,'grid_receptor_max_error_nmol':maxrec,'refinement_relative_state_error':float(np.max(np.abs(base.y-fine.y))/300),'zero_receptor_brain_nmol':float(np.max(np.abs(zero.y[6]))),'no_internalization_brain_nmol':float(np.max(np.abs(noentry.y[6]))),'minimum_state_nmol':float(base.y.min()),'vascular_Kp_increment_error':float(abs(app-kp-.05))}
    checks['passed']=bool(maxerr<1e-5 and maxrec<1e-8 and checks['refinement_relative_state_error']<1e-6 and checks['zero_receptor_brain_nmol']<1e-12 and checks['no_internalization_brain_nmol']<1e-12 and checks['minimum_state_nmol']>-1e-8 and checks['vascular_Kp_increment_error']<1e-12)
    nominal=[r for r in rows if r[0]=='pH_sensitivity' and r[3]==.7];best=max(nominal,key=lambda r:r[7])
    summary={'task':13,'evidence':'executed finite-pool, uncalibrated antibody trafficking','grid_runs':len(rows),'nominal_Kp_extravascular_0_72h':kp,'nominal_Kp_with_vascular_contamination':app,'nominal_brain_AUC_nM_h':abr,'nominal_best_KD_nM_in_grid':best[1],'best_at_grid_boundary':bool(best[1] in [kds[0],kds[-1]]),'best_within_requested_50_500_nM':bool(50<=best[1]<=500),'limitation':'No species-specific receptor abundance, endogenous transferrin competition, avidity or brain target binding. All brain antibody is free; no Kp,uu inference for a target-binding antibody.'}
    for name,obj in [('summary.json',summary),('verification.json',checks)]: (out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False});fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for slope in [0.,.7,1.4]:
        r=np.array([r[1:] for r in rows if r[0]=='pH_sensitivity' and r[3]==slope]);ax[0,0].semilogx(r[:,0],r[:,6],label=f'pH slope {slope:g}')
    ax[0,0].set(xlabel='Apical KD (nM)',ylabel='Extravascular brain / plasma AUC',title='A  Affinity optimum is a calculated outcome');ax[0,0].legend(fontsize=8)
    for kd in [.1,100.,1e4]:
        s=simulate(kd);ax[0,1].plot(s.t,s.y[6]/.3,label=f'KD {kd:g} nM');ax[1,0].plot(s.t,s.y[4],label=f'KD {kd:g} nM')
    ax[0,1].set(xlabel='Time (h)',ylabel='Free brain antibody (nM)',title='B  Finite-dose brain pharmacokinetics');ax[0,1].legend(fontsize=8)
    ax[1,0].set(xlabel='Time (h)',ylabel='Endosomal receptor complex (nmol)',title='C  Endosomal retention and cargo loss');ax[1,0].legend(fontsize=8)
    ax[1,1].bar(['Extravascular','With vascular carryover'],[kp,app],color=['#0e7490','#eab308']);ax[1,1].set(ylabel='Concentration AUC ratio',title='D  Measurement definition changes apparent Kp')
    fig.suptitle('Task 13 | BBB–TfR shuttle: conserved amounts, conditional transport tradeoffs',fontweight='bold');fig.savefig(out/'figures'/FIGURE,dpi=300);plt.close(fig)
    if not checks['passed']:raise AssertionError(checks)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);print(json.dumps(run(p.parse_args().out),indent=2))
