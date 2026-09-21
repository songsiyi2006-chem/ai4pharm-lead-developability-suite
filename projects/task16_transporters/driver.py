"""Actual Task 4 PBPK extension with finite intracellular transporter amounts."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE = 'fig16_transporter_oatp_pgp_kinetics.png'
TASK4_PATH = Path(__file__).resolve().parents[1] / 'task04_pbpk' / 'run_task4_pbpk_pharmacokinetics_dose_prediction.py'
spec = importlib.util.spec_from_file_location('task16_existing_task4', TASK4_PATH)
task4 = importlib.util.module_from_spec(spec); sys.modules[spec.name] = task4; spec.loader.exec_module(task4)
HC, RC, EC = 11, 12, 13
CONFIG = {'dose_mg': 100., 'horizon_h': 96., 'OATP1B1_Vmax_mg_h': 12., 'OATP1B3_Vmax_mg_h': 4.,
          'OATP_Km_mg_l': .1, 'hepatic_passive_l_h': 2., 'hepatocyte_volume_l': 1.3,
          'intracellular_fu': .2, 'hepatic_metabolic_CL_l_h': 5., 'biliary_Vmax_mg_h': 3.,
          'biliary_Km_mg_l': .1, 'renal_cell_volume_l': .15, 'renal_uptake_Vmax_mg_h': 3.,
          'renal_uptake_Km_mg_l': .2, 'renal_passive_l_h': .5,
          'renal_Pgp_BCRP_Vmax_mg_h': [.8, .6], 'renal_efflux_Km_mg_l': .15,
          'enterocyte_volume_l': .25, 'enterocyte_basolateral_h': .8,
          'gut_Pgp_BCRP_Vmax_mg_h': [2., 1.], 'gut_efflux_Km_mg_l': .2,
          'gut_fecal_transit_h': .1, 'inhibitor_I_over_Ki': [0., 1., 10., 100.],
          'provenance': 'transporter and cell parameters assumed; inherited Task4 physiology is scenario data',
          'clinical_recommendations': None}

def write_json(p, v): p.write_text(json.dumps(v, indent=2, allow_nan=False) + '\n', encoding='utf-8')
def write_csv(p, h, r):
    with p.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(h); w.writerows(r)

def michaelis(c, vmax, km, inhibition=0.):
    if vmax < 0 or km <= 0 or inhibition < 0: raise ValueError('invalid transport parameters')
    c = max(float(c), 0.)
    return vmax * c / (km * (1. + inhibition) + c)

class TransportPBPK:
    def __init__(self, inhibition=0.):
        if inhibition < 0: raise ValueError('inhibition cannot be negative')
        self.inhibition = inhibition
        self.base = task4.PBPK(task4.Config(fa=1., fg=1.))

    def derivative(self, t, y):
        b = self.base; c = b.c
        d = np.zeros(14); d[:11] = b.derivative(y[:11])
        # Task4 circulation and glomerular filtration remain in force.
        # Remove its complete perfusion-limited hepatic metabolic flux first.
        cu_liver = c.fu_p * y[task4.LIVER] / (c.volumes_l['liver'] * b.kp['liver'])
        old_hepatic = b.clint * cu_liver
        d[task4.LIVER] += old_hepatic; d[task4.HEP] -= old_hepatic
        cu_kidney = c.fu_p * y[task4.KIDNEY] / (c.volumes_l['kidney'] * b.kp['kidney'])
        cu_h = .2 * y[HC] / 1.3; cu_r = .2 * y[RC] / .15; cu_e = .2 * y[EC] / .25
        influx = (michaelis(cu_liver, 12., .1, self.inhibition) +
                  michaelis(cu_liver, 4., .1, self.inhibition) + 2. * (cu_liver - cu_h))
        metabolism = 5. * cu_h; bile = michaelis(cu_h, 3., .1)
        d[task4.LIVER] -= influx; d[HC] += influx - metabolism - bile
        d[task4.HEP] += metabolism; d[task4.FEC] += bile
        renal_in = michaelis(cu_kidney, 3., .2) + .5 * (cu_kidney - cu_r)
        renal_out = michaelis(cu_r, .8, .15) + michaelis(cu_r, .6, .15)
        d[task4.KIDNEY] -= renal_in; d[RC] += renal_in - renal_out; d[task4.REN] += renal_out
        # Replace the existing direct gut-depot -> liver shortcut with gut cells.
        uptake = c.ka_h * y[task4.GUT]
        d[task4.LIVER] -= uptake
        gut_out = michaelis(cu_e, 2., .2) + michaelis(cu_e, 1., .2)
        basolateral = .8 * y[EC]; fecal = .1 * y[task4.GUT]
        d[EC] += uptake - gut_out - basolateral
        d[task4.GUT] += gut_out - fecal; d[task4.FEC] += fecal
        d[task4.LIVER] += basolateral
        return d

    def simulate(self, dose=100., route='IV', rtol=2e-8, samples=961):
        if dose <= 0 or route not in {'IV', 'oral'}: raise ValueError('positive dose and IV/oral route required')
        y0 = np.zeros(14); y0[task4.BLOOD if route == 'IV' else task4.GUT] = dose
        t = np.linspace(0., 96., samples)
        sol = solve_ivp(self.derivative, (0., 96.), y0, method='Radau', t_eval=t,
                        rtol=rtol, atol=rtol * .01)
        if not sol.success: raise RuntimeError(sol.message)
        return sol.t, sol.y

def mass_ledger(y): return y[:10].sum(axis=0) + y[11:].sum(axis=0)

def run(out: Path) -> dict:
    out = Path(out)
    if out.exists(): raise FileExistsError(f'Output must be new: {out}')
    out.mkdir(parents=True); (out / 'figures').mkdir()
    write_json(out / 'config.json', CONFIG)
    rows = []; stats = []; runs = {}; err = 0.; minimum = 0.
    for route in ['IV', 'oral']:
        for inhib in CONFIG['inhibitor_I_over_Ki']:
            t,y = TransportPBPK(inhib).simulate(route=route); runs[(route,inhib)] = (t,y)
            err = max(err, float(np.max(np.abs(mass_ledger(y) - 100.))))
            minimum = min(minimum, float(y.min()))
            stats.append({'route': route, 'I_over_Ki': inhib, 'AUC_0_96_mg_h_l': float(y[10,-1]),
                          'Cmax_mg_l': float(y[0].max()/5), 'remaining_amount_mg': float(y[:7,-1].sum()+y[11:,-1].sum()),
                          'metabolism_mg': float(y[7,-1]), 'urinary_mg': float(y[8,-1]), 'fecal_biliary_mg': float(y[9,-1])})
            rows.extend([[route,inhib,ti,*y[:,j]] for j,ti in enumerate(t)])
    header = ['route','I_over_Ki','time_h',*task4.NAMES,'metabolic_loss_mg','urinary_loss_mg','fecal_biliary_loss_mg',
              'plasma_AUC_mg_h_l','hepatocyte_mg','renal_cell_mg','enterocyte_mg']
    write_csv(out / 'pbpk_transporter_trajectories.csv', header, rows)
    sweep=[]
    for dose in [10.,100.,300.]:
        baseline = None
        for inhib in [0.,.1,.3,1.,3.,10.,30.,100.]:
            t,y = TransportPBPK(inhib).simulate(dose=dose, samples=49)
            auc=float(y[10,-1]); baseline=auc if baseline is None else baseline
            sweep.append([dose,inhib,auc,auc/baseline,float(y[0].max()/5.)])
    write_csv(out / 'dose_inhibition_sensitivity.csv',['dose_mg','I_over_Ki','AUC_0_96','AUCR','Cmax_mg_l'],sweep)
    _,ref = TransportPBPK(10.).simulate(rtol=2e-10)
    refine=float(np.max(np.abs(ref-runs[('IV',10.)][1])) / max(1.,np.max(np.abs(ref))))
    check={'mass_balance_error_mg':err,'minimum_state':minimum,'refinement_scaled_error':refine,
           'task4_derivative_reused':True,'old_hepatic_flux_removed_before_extension':True,
           'passed':bool(err<1e-5 and minimum>-1e-7 and refine<2e-5)}
    if not check['passed']: raise RuntimeError(str(check))
    for route in ['IV','oral']:
        base=next(x['AUC_0_96_mg_h_l'] for x in stats if x['route']==route and x['I_over_Ki']==0.)
        for x in stats:
            if x['route']==route: x['AUCR_0_96']=x['AUC_0_96_mg_h_l']/base
    summary={'task':16,'evidence':'synthetic extended PBPK with explicit mass conservation','scenarios':stats,
             'clinical_recommendations':None,'reused_module':TASK4_PATH.name,
             'limitations':['No measured transporter kinetics or tissue calibration.',
                            'Inhibitor is a constant imposed unbound I/Ki ratio, not perpetrator PK.',
                            'Inherited tissue volumes are effective exchange volumes; added cell volumes are an unvalidated split.',
                            'AUC is truncated at 96 hours; no patient dose prediction.']}
    write_json(out/'summary.json',summary); write_json(out/'verification.json',check)
    fig,ax=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    conc=np.geomspace(.001,10.,120)
    for inhib in [0.,1.,10.,100.]: ax[0,0].semilogx(conc,[michaelis(x,16.,.1,inhib) for x in conc],label=f'I/Ki={inhib:g}')
    ax[0,0].set(xlabel='Unbound donor concentration (mg/L)',ylabel='OATP1B1 + 1B3 flux (mg/h)',title='A  Competitive uptake inhibition');ax[0,0].legend()
    for inhib in CONFIG['inhibitor_I_over_Ki']:
        t,y=runs[('IV',inhib)]; ax[0,1].semilogy(t,y[0]/5.,label=f'I/Ki={inhib:g}')
    ax[0,1].set(xlabel='Time (h)',ylabel='Plasma concentration (mg/L)',title='B  Same 100 mg synthetic IV input');ax[0,1].legend()
    for dose in [10.,100.,300.]:
        a=np.asarray([r for r in sweep if r[0]==dose]);ax[1,0].semilogx(a[1:,1],a[1:,3],'o-',label=f'{dose:g} mg')
    ax[1,0].set(xlabel='Imposed I/Ki',ylabel='AUC0–96 ratio',title='C  Saturation and exposure interaction');ax[1,0].legend()
    t,y=runs[('oral',0.)]
    for idx,label in [(HC,'hepatocyte'),(RC,'renal cell'),(EC,'enterocyte')]:ax[1,1].plot(t,y[idx],label=label)
    ax[1,1].set(xlabel='Time (h)',ylabel='Intracellular amount (mg)',title='D  Finite cellular donor pools, oral');ax[1,1].legend()
    for a in ax.flat:a.grid(alpha=.2)
    fig.suptitle('Task 16 | Extended Task 4 circulation — assumed transporter kinetics',fontsize=14)
    fig.savefig(out/'figures'/FIGURE,dpi=300);plt.close(fig)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path(__file__).parent/'outputs')
    print(json.dumps(run(p.parse_args().out),indent=2))
