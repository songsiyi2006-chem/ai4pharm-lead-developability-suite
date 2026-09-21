#!/usr/bin/env python3
"""Single entry point for Tasks 11–20; archive, CPU recomputation and ZIP delivery.

Independent project folders keep methods and reports locatable. --bundle writes
one self-contained source/data/report archive; installed scientific runtimes are
external dependencies. No environment installation, network or Git side effects.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, shutil, signal, subprocess, sys, time, zipfile
from pathlib import Path
for _key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
os.environ['MPLBACKEND']='Agg'
ROOT=Path(__file__).resolve().parent
PROJECTS=('task12_chromatin','task13_bbb_tfr','task14_bite','task15_receptor_signaling',
          'task16_transporters','task18_pgx','task19_metadynamics','task20_denovo')
PUBLISHED_TASKS=(12,13,14,15,16,18,19,20)
TITLES=('染色质靶蛋白降解','BBB–TfR 跨胞转运','BiTE 三元平衡','受体校验与细胞因子',
        '主动转运与 PBPK','药物基因组学','全原子三元界面采样','三维口袋分子生成')
FIGURES=('fig12_epigenetic_chromatin_depletion.png','fig13_bbb_tfr_transcytosis_profile.png',
         'fig14_bite_synapse_crosslinking_curve.png','fig15_receptor_proofreading_cytokine_balance.png',
         'fig16_transporter_oatp_pgp_kinetics.png','fig18_pgx_cyp2d6_phenotype_kinetics.png',
         'fig19_ternary_metadynamics_pmf_landscape.png','fig20_denovo_pareto_lead_optimization.png')
BOUNDARIES=('模型区分自由态和染色质结合态，解耦罚能尚无实验校准。',
 '质量守恒的转运情景；最佳亲和力由参数决定，尚无临床外推。',
 '有限总量质量作用平衡；hook 与杀伤饱和分别计算。',
 'IL6R 阻断降低信号，游离 IL6 可升高；没有休克预测。',
 '实际复用 Task 4 循环；载体参数与 DDI 比值为情景结果。',
 '撞击器质量分布与肺区域沉积分开；吸收和清除需实测校准。',
 '采用当前基因分型标准；暴露匹配只是模型输出，无临床处方。',
 '真实全原子短程 pilot；自由能未收敛，协同性 alpha 保留未知。',
 '实际图生成、三维评分与 Vina 对照；nM 亲和力和合成路线未证实。')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def prepare_output(path):
    out=Path(path).resolve()
    if out==ROOT or ROOT.is_relative_to(out):raise ValueError('Output overlaps repository or ancestor')
    for folder in ('projects','docs','tests','tools','.git','figures_advanced','data_advanced','omnibus','figures_omnibus','data_omnibus'):
        if out.is_relative_to(ROOT/folder):raise ValueError('Use a new work/ directory, outside scientific archives')
    if out.exists():raise FileExistsError('Output exists; select a new directory')
    out.mkdir(parents=True);return out
def task_command(task,target,python,atomistic_python=None,vina=None,preparation_python=None):
    if task not in range(11,21):raise ValueError('Task must be 11..20')
    executable=atomistic_python if task==19 and atomistic_python else python
    command=[str(executable),str(ROOT/'projects'/PROJECTS[task-11]/'driver.py'),'--out',str(target)]
    if task==20 and vina:
        command+=['--vina',str(Path(vina).resolve())]
        if preparation_python:command+=['--preparation-python',str(Path(preparation_python).resolve())]
    return command
def execute(command,log,timeout):
    opts={'creationflags':subprocess.CREATE_NEW_PROCESS_GROUP} if os.name=='nt' else {'start_new_session':True}
    with log.open('w',encoding='utf-8') as handle, subprocess.Popen(command,cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT,**opts) as proc:
        try:return proc.wait(timeout=timeout)
        except (subprocess.TimeoutExpired,KeyboardInterrupt):
            if os.name=='nt':subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=False)
            else:
                try:os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError:pass
            if proc.poll() is None:proc.kill()
            proc.wait();raise
def project_files(project):
    base=ROOT/'projects'/project
    return sorted(p for p in base.rglob('*') if p.is_file() and
                  not any(x in ('work','__pycache__','.local_run_history','.pytest_cache') for x in p.relative_to(base).parts))
def collect(out,tasks,mode,calculations,execution):
    sources={};records=[]
    for task in tasks:
        index=PUBLISHED_TASKS.index(task);project=PROJECTS[index];base=calculations[task];filename=FIGURES[index]
        summary=json.loads((base/'summary.json').read_text(encoding='utf-8-sig'))
        if task==19 and not summary.get('atomistic',{}).get('status','').startswith('completed_'):
            raise ValueError('Task 19 atomistic pilot was not completed; inspect its execution log')
        destination=out/'figures_advanced'/filename;destination.parent.mkdir(exist_ok=True)
        shutil.copy2(base/'figures'/filename,destination)
        dump(out/'data_advanced'/'summaries'/f'task{task}.json',summary)
        archived_files=project_files(project)
        for p in archived_files:sources[p.relative_to(ROOT).as_posix()]=sha(p)
        records.append({'task':task,'project':project,'title':TITLES[index],'file':filename,
                        'sha256':sha(destination),'summary':summary,'evidence_boundary':BOUNDARIES[index],
                        'calculation_output':str(base.resolve()),'output_sha256':{p.relative_to(base).as_posix():sha(p) for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts}})
    for p in (ROOT/Path(__file__).name,ROOT/'projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py'):
        sources[p.relative_to(ROOT).as_posix()]=sha(p)
    for lang in ('ZH','EN'):
        p=ROOT/f'AI4PHARM_ADVANCED_TREATISE_{lang}.md'
        if not p.exists():raise FileNotFoundError(f'Master treatise absent: {p.name}')
        shutil.copy2(p,out/p.name)
    dump(out/'data_advanced/figure_provenance.json',records)
    dump(out/'data_advanced/execution.json',execution)
    dump(out/'data_advanced/run_summary.json',{'tasks':tasks,'mode':mode,'python':platform.python_version(),
         'numerical_threads':1,'scientific_status':'Executed calculations, with task-specific calibration and acceptance limits; no clinical recommendations'})
    dump(out/'manifest_advanced.json',{'schema':'ai4pharm-advanced-1','tasks':tasks,'mode':mode,
         'source_sha256':sources,'files_sha256':{p.relative_to(out).as_posix():sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest_advanced.json'}})
    return records
def bundle(path):
    path=Path(path).resolve()
    if path.exists():raise FileExistsError('Bundle already exists')
    if path.suffix!='.zip':raise ValueError('Bundle must end in .zip')
    files=set()
    for project in PROJECTS:files.update(project_files(project))
    files.update(project_files('task04_pbpk'))
    files.update(ROOT.glob('AI4PHARM_ADVANCED_TREATISE_*.md'))
    for folder in ('figures_advanced','data_advanced','docs/advanced'):
        files.update(p for p in (ROOT/folder).rglob('*') if p.is_file())
    files.update([Path(__file__).resolve(),ROOT/'requirements.txt',ROOT/'LICENSE',ROOT/'manifest_advanced.json'])
    path.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files):z.write(p,p.relative_to(ROOT).as_posix())
        z.writestr('BUNDLE_README.txt','AI4Pharm Tasks 11–20 source, inputs, executed outputs and bilingual reports.\nExtract and run: python run_ai4pharm_advanced_modalities_suite.py --out work/new_archive\nSee docs/advanced/README.md for runtimes, optional all-atom dependencies and Vina.\n')
    return {'bundle':str(path),'bytes':path.stat().st_size,'sha256':sha(path)}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path);p.add_argument('--mode',choices=('archived','recompute'),default='archived')
    p.add_argument('--tasks',nargs='+',type=int,default=PUBLISHED_TASKS)
    p.add_argument('--task-timeout',type=float,default=3600);p.add_argument('--atomistic-python',type=Path)
    p.add_argument('--vina',type=Path);p.add_argument('--preparation-python',type=Path);p.add_argument('--bundle',type=Path)
    a=p.parse_args()
    if a.bundle:
        print(json.dumps(bundle(a.bundle),indent=2));return
    if not a.out:p.error('--out is required unless --bundle is selected')
    tasks=sorted(set(a.tasks))
    if not tasks or any(t not in PUBLISHED_TASKS for t in tasks):p.error('Only publishable tasks 12,13,14,15,16,18,19,20 are enabled')
    if a.task_timeout<=0:p.error('--task-timeout must be positive')
    out=prepare_output(a.out);execution=[];calculations={}
    for task in tasks:
        project=PROJECTS[PUBLISHED_TASKS.index(task)]
        if a.mode=='archived':calculations[task]=ROOT/'projects'/project/'outputs';continue
        target=out/'calculations'/project/'outputs';target.parent.mkdir(parents=True)
        command=task_command(task,target,sys.executable,a.atomistic_python,a.vina,a.preparation_python)
        log=out/'logs'/f'task{task}.txt';log.parent.mkdir(exist_ok=True);start=time.perf_counter()
        print(f'Task {task}: running on CPU; log={log}',flush=True)
        status='failed';code=None
        try:code=execute(command,log,a.task_timeout);status='passed' if code==0 else 'failed'
        except subprocess.TimeoutExpired:status='timeout'
        execution.append({'task':task,'status':status,'returncode':code,'seconds':time.perf_counter()-start,'command':command,'log':log.relative_to(out).as_posix()})
        dump(out/'data_advanced/execution.json',execution)
        if status!='passed':raise RuntimeError(f'Task {task} {status}; partial outputs and log retained')
        calculations[task]=target
    collect(out,tasks,a.mode,calculations,execution)
    print(json.dumps({'status':'delivered','mode':a.mode,'tasks':tasks,'out':str(out)},indent=2))
if __name__=='__main__':main()
