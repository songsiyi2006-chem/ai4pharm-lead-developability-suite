"""CPIC phenotype translation plus explicitly hypothetical parent/metabolite PK."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIGURE='fig18_pgx_cyp2d6_phenotype_kinetics.png'
CYP2D6_ALLELES={'*1':1.,'*2':1.,'*4':0.,'*5':0.,'*10':.25,'*17':.5,'*41':.25}
CYP2C19_ALLELES={'*1':'normal','*2':'none','*3':'none','*17':'increased','*9':'decreased','*12':'unknown'}
COHORTS={'CYP2D6':[('*4','*4'),('*1','*4'),('*1','*1'),('*1x2','*1')],
         'CYP2C19':[('*2','*2'),('*1','*2'),('*1','*1'),('*1','*17'),('*17','*17')]}
FACTORS={'PM':0.,'IM':.35,'NM':1.,'RM':1.4,'UM':1.8}
CONFIG={'CYP2D6_allele_subset_activity':CYP2D6_ALLELES,'CYP2C19_allele_subset_function':CYP2C19_ALLELES,
        'hypothetical_substrate_phenotype_clearance_factors':FACTORS,
        'parent_volume_l':50.,'metabolite_volume_l':30.,'non_CYP_clearance_l_h':2.,
        'NM_CYP_clearance_l_h':8.,'metabolite_clearance_l_h':4.,'metabolite_molar_yield':.8,
        'oral_bioavailability':.85,'absorption_h':1.,'dose_nmol':1000.,'dose_interval_h':24.,'number_doses':7,
        'phenoconversion_inhibitor_I_over_Ki':[0.,1.,100.],
        'classification_source':'CPIC CYP2D6 2026 update; CYP2C19 2022 clopidogrel guideline',
        'PK_provenance':'all PK parameters and factor mapping assumed; not CPIC dosing recommendations',
        'clinical_recommendations':None}

def write_json(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def write_csv(p,h,r):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(h);w.writerows(r)

def cyp2d6_phenotype(score):
    if score is None:return 'indeterminate'
    if not np.isfinite(score) or score<0:raise ValueError('activity score must be finite and nonnegative')
    if score==0:return 'PM'
    if score<1.25:return 'IM'
    if score<=2.25:return 'NM'
    return 'UM'

def cyp2d6_diplotype(alleles):
    if len(alleles)!=2:raise ValueError('two phased allele entries required')
    score=0.
    for entry in alleles:
        parts=entry.split('x');base=parts[0]
        if len(parts)>2:raise ValueError('unsupported copy notation')
        copies=int(parts[1]) if len(parts)==2 else 1
        if copies<1:raise ValueError('gene deletion must be represented as *5')
        if base not in CYP2D6_ALLELES:return {'activity_score':None,'phenotype':'indeterminate'}
        score+=copies*CYP2D6_ALLELES[base]
    return {'activity_score':score,'phenotype':cyp2d6_phenotype(score)}

def cyp2c19_diplotype(alleles):
    if len(alleles)!=2:raise ValueError('two allele entries required')
    kinds=[CYP2C19_ALLELES.get(a,'unknown') for a in alleles]
    if 'unknown' in kinds:return 'indeterminate'
    if kinds.count('none')==2:return 'PM'
    if 'decreased' in kinds:return 'likely_PM' if 'none' in kinds else 'likely_IM'
    if 'none' in kinds:return 'IM' # *17 does not compensate for a no-function allele.
    if kinds.count('increased')==2:return 'UM'
    if 'increased' in kinds:return 'RM'
    return 'NM'

def pk_matrix(factor=1.,inhibitor=0.,enzyme_fraction=.8):
    if factor<0 or inhibitor<0 or not 0<=enzyme_fraction<=1:raise ValueError('invalid PK parameters')
    cl_cyp=10.*enzyme_fraction*factor/(1.+inhibitor);cl_other=10.*(1-enzyme_fraction)
    k=np.zeros((6,6));ka=1.;f=.85;vp=50.;vm=30.;clm=4.;yield_m=.8
    k[0,0]=-ka;k[1,0]=f*ka;k[3,0]=(1-f)*ka
    k[1,1]=-(cl_cyp+cl_other)/vp;k[2,1]=yield_m*cl_cyp/vp
    k[2,2]=-clm/vm;k[3,1]=(cl_other+(1-yield_m)*cl_cyp)/vp;k[3,2]=clm/vm
    k[4,1]=1/vp;k[5,2]=1/vm
    return k

def simulate(factor=1.,inhibitor=0.,enzyme_fraction=.8,step=.25):
    nstep=round(24./step)
    if nstep<1 or not np.isclose(nstep*step,24.):raise ValueError('step must divide 24 h')
    k=pk_matrix(factor,inhibitor,enzyme_fraction);p=expm(k*step)
    y=np.zeros(6);times=[];states=[];dose_count=[]
    for day in range(7):
        y[0]+=1000.
        # At a dosing boundary, retain the post-dose state, not duplicate time rows.
        for j in range(nstep):
            times.append(day*24+j*step);states.append(y.copy());dose_count.append(day+1)
            y=p@y
    times.append(168.);states.append(y.copy());dose_count.append(7)
    return np.array(times),np.array(states).T,np.array(dose_count)

def run(out:Path)->dict:
    out=Path(out)
    if out.exists():raise FileExistsError(f'Output must be new: {out}')
    out.mkdir(parents=True);(out/'figures').mkdir();write_json(out/'config.json',CONFIG)
    classifications=[];stats=[];rows=[];runs={};err=0.;minimum=0.
    for gene,diplotypes in COHORTS.items():
        for alleles in diplotypes:
            result=cyp2d6_diplotype(alleles) if gene=='CYP2D6' else {'activity_score':None,'phenotype':cyp2c19_diplotype(alleles)}
            pheno=result['phenotype'];label='/'.join(alleles);factor=FACTORS[pheno]
            classifications.append([gene,label,result['activity_score'],pheno,factor,'assumed PK factor, not AS-to-clearance calibration'])
            for inhibitor in [0.,1.,100.]:
                t,y,n=simulate(factor,inhibitor);runs[(gene,pheno,inhibitor)]=(t,y)
                err=max(err,float(np.max(np.abs(y[:4].sum(axis=0)-1000*n))))
                minimum=min(minimum,float(y.min()))
                record={'gene':gene,'diplotype':label,'genetic_phenotype':pheno,'I_over_Ki':inhibitor,
                        'simulated_CYP_factor':factor/(1+inhibitor),'parent_AUC_0_168_nM_h':float(y[4,-1]),
                        'active_metabolite_AUC_0_168_nM_h':float(y[5,-1]),'parent_Cmax_nM':float((y[1]/50).max()),
                        'clinical_recommendations':None}
                stats.append(record)
                rows.extend([[gene,label,pheno,inhibitor,tj,*y[:,j],n[j]] for j,tj in enumerate(t)])
    write_csv(out/'genotype_translation.csv',['gene','diplotype','activity_score','phenotype','assumed_PK_factor','factor_evidence'],classifications)
    write_csv(out/'parent_metabolite_trajectories.csv',['gene','diplotype','genetic_phenotype','I_over_Ki','time_h','depot_nmol','parent_nmol','metabolite_nmol','eliminated_nmol_equivalent','parent_AUC_nM_h','metabolite_AUC_nM_h','doses_received'],rows)
    factors=[]
    for gene in COHORTS:
        ref=next(s for s in stats if s['gene']==gene and s['genetic_phenotype']=='NM' and s['I_over_Ki']==0.)
        for s in stats:
            if s['gene']!=gene:continue
            p=ref['parent_AUC_0_168_nM_h']/s['parent_AUC_0_168_nM_h']
            m=ref['active_metabolite_AUC_0_168_nM_h']/s['active_metabolite_AUC_0_168_nM_h'] if s['active_metabolite_AUC_0_168_nM_h']>1e-14 else None
            factors.append([gene,s['genetic_phenotype'],s['I_over_Ki'],p,m,'mathematical exposure scaling only; no dose recommendation'])
    write_csv(out/'exposure_normalization_factors.csv',['gene','phenotype','I_over_Ki','parent_normalization','active_metabolite_normalization','interpretation'],factors)
    sweep=[]
    for fm in [.3,.6,.8,.95]:
        _,ref,_=simulate(1.,0.,fm)
        for factor in np.linspace(0,2.,21):
            _,y,_=simulate(float(factor),0.,fm)
            sweep.append([fm,factor,y[4,-1]/ref[4,-1],y[5,-1]/ref[5,-1]])
    write_csv(out/'enzyme_fraction_sensitivity.csv',['NM_enzyme_clearance_fraction','activity_factor','parent_AUCR','metabolite_AUCR'],sweep)
    t,y,_=simulate(1.,1.,step=.125);_,coarse,_=simulate(1.,1.)
    refine=float(np.max(np.abs(y[:,::2]-coarse))/max(1.,np.max(np.abs(y))))
    verification={'molar_equivalent_ledger_error_nmol':err,'minimum_state':minimum,'matrix_step_refinement_scaled_error':refine,
                  'PM_active_metabolite_AUC':float(runs[('CYP2D6','PM',0.)][1][5,-1]),
                  'passed':bool(err<1e-5 and minimum>-1e-8 and refine<1e-10)}
    if not verification['passed']:raise RuntimeError(str(verification))
    summary={'task':18,'evidence':'source-backed phenotype labels; hypothetical substrate PK',
             'classification_date':'2026-09-21','scenarios':stats,'clinical_recommendations':None,
             'limitations':['No real patient genotype or drug-specific dose model.',
                            'Activity score classifies CYP2D6; clearance factors are independent assumptions.',
                            'CYP2C19 uses allele-function combinations, including a separate rapid phenotype.',
                            'Inhibitor changes model clearance without rewriting inherited genotype.',
                            'Active parent and bioactivated metabolite can imply opposite exposure changes.']}
    write_json(out/'summary.json',summary);write_json(out/'verification.json',verification)
    fig,ax=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    colors={'PM':'#b91c1c','IM':'#d97706','NM':'#059669','UM':'#2563eb'}
    for pheno in ['PM','IM','NM','UM']:
        t,y=runs[('CYP2D6',pheno,0.)]
        ax[0,0].plot(t,y[1]/50.,label=pheno,color=colors[pheno]);ax[0,1].plot(t,y[2]/30.,label=pheno,color=colors[pheno])
    ax[0,0].set(xlabel='Time (h)',ylabel='Parent concentration (nM)',title='A  Same hypothetical oral input');ax[0,0].legend()
    ax[0,1].set(xlabel='Time (h)',ylabel='Active metabolite concentration (nM)',title='B  Bioactivation changes the interpretation');ax[0,1].legend()
    for inhibitor in [0.,1.,100.]:
        t,y=runs[('CYP2D6','NM',inhibitor)];ax[1,0].plot(t,y[1]/50.,label=f'NM, I/Ki={inhibitor:g}')
    ax[1,0].set(xlabel='Time (h)',ylabel='Parent concentration (nM)',title='C  Phenoconversion scenario');ax[1,0].legend()
    arr=[next(r for r in factors if r[0]=='CYP2D6' and r[1]==ph and r[2]==0.) for ph in ['PM','IM','NM','UM']]
    x=np.arange(4);ax[1,1].bar(x-.18,[r[3] for r in arr],.36,label='parent AUC equalization')
    ax[1,1].bar(x+.18,[r[4] if r[4] is not None else np.nan for r in arr],.36,label='metabolite AUC equalization')
    ax[1,1].text(0,.27,'PM metabolite:\nno finite factor',fontsize=8,ha='center')
    ax[1,1].set(xticks=x,xticklabels=['PM','IM','NM','UM'],ylabel='Input scaling factor',title='D  Model exposure normalization, not dosing');ax[1,1].legend(fontsize=8)
    for a in ax.flat:a.grid(alpha=.2)
    fig.suptitle('Task 18 | Current CPIC categories + uncalibrated PK scenarios',fontsize=14)
    fig.savefig(out/'figures'/FIGURE,dpi=300);plt.close(fig)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path(__file__).parent/'outputs')
    print(json.dumps(run(p.parse_args().out),indent=2))
