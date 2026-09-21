"""Verify advanced delivery bytes, figure metadata and executed evidence fields."""
from __future__ import annotations
import argparse,json,sys,hashlib
from pathlib import Path
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from run_ai4pharm_advanced_modalities_suite import PROJECTS,FIGURES,PUBLISHED_TASKS
def inspect(bundle=ROOT,source=ROOT):
    errors=[];count=0
    manifest=json.loads((bundle/'manifest_advanced.json').read_text(encoding='utf-8'))
    if manifest.get('schema')!='ai4pharm-advanced-1':errors.append('Unexpected schema')
    tasks=manifest.get('tasks',[])
    for field,base in [('files_sha256',bundle),('source_sha256',source)]:
        mapping=manifest.get(field,{})
        if not mapping:errors.append('Missing manifest mapping: '+field)
        for name,value in mapping.items():
            path=(base/name).resolve();count+=1
            if not path.is_relative_to(base.resolve()) or not path.is_file():errors.append('Missing or escaping path: '+name)
            elif hashlib.sha256(path.read_bytes()).hexdigest()!=value:errors.append('Checksum mismatch: '+name)
    records=json.loads((bundle/'data_advanced/figure_provenance.json').read_text(encoding='utf-8'))
    if sorted(x['task'] for x in records)!=tasks:errors.append('Task coverage differs')
    images=[]
    for t in tasks:
        filename=FIGURES[PUBLISHED_TASKS.index(t)]
        with Image.open(bundle/'figures_advanced'/filename) as im:
            dpi=im.info.get('dpi',())
            if len(dpi)!=2 or any(abs(v-300)>.02 for v in dpi):errors.append('DPI mismatch: '+filename)
            if min(im.size)<1500:errors.append('Image too small: '+filename)
            images.append({'task':t,'pixels':list(im.size),'dpi':list(dpi)});im.verify()
        for lang in ('ZH','EN'):
            name=f'AI4PHARM_ADVANCED_TREATISE_{lang}.md'
            if 'figures_advanced/'+filename not in (bundle/name).read_text(encoding='utf-8'):errors.append('Report missing figure: '+name+' '+filename)
    if 19 in tasks:
        s=json.loads((bundle/'data_advanced/summaries/task19.json').read_text())['atomistic']
        if not s['status'].startswith('completed_'):errors.append('Task 19 atomistic pilot incomplete')
        if s['pmf_converged'] or s['cooperativity_alpha'] is not None:errors.append('Task 19 unsupported thermodynamic claim')
    if 20 in tasks:
        s=json.loads((bundle/'data_advanced/summaries/task20.json').read_text())
        if s['population']!=100 or s['generations']!=20 or s['history_rows']!=2000:errors.append('Task 20 full dimensions not executed')
        if s['affinity_nM'] is not None:errors.append('Unsupported calibrated affinity')
        if s['random_search_control']['unique_graph_budget']!=s['unique_graphs_evaluated']:errors.append('Random-control budget differs')
    return {'status':'passed' if not errors else 'failed','tasks':tasks,'hashes_checked':count,'figures':images,'errors':errors,
            'scope':'Artifact integrity, executed dimensions and scientific acceptance fields; not experimental validation'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path);p.add_argument('--bundle',type=Path,default=ROOT);a=p.parse_args()
    result=inspect(a.bundle);text=json.dumps(result,indent=2,ensure_ascii=False)
    if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(text+'\n',encoding='utf-8')
    print(text);raise SystemExit(0 if result['status']=='passed' else 1)
