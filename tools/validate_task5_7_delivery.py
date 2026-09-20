"""Read-only integration checks for Tasks 5-7; no scientific outputs overwritten."""
from pathlib import Path
import argparse,ast,hashlib,json,re,subprocess
from urllib.parse import unquote,urlsplit
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
BASE='ced2ff75baa0810b1a61db54a9386b1af81b3664'
TASKS=('task5_cyp_ddi','task6_asd','task7_qsp')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inspect(out,reproductions=None):
    out.mkdir(parents=True,exist_ok=True)
    changed=subprocess.check_output(['git','diff',BASE,'--name-only'],cwd=ROOT,text=True).splitlines()
    allowed={'.gitattributes','README.md','manifest_task2.json'}
    included=set(subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines())
    legacy_changes=[p for p in changed if p not in allowed and not p.startswith(('docs/','tools/validate_task5_7_delivery.py','tests/test_task5','tests/test_task6','tests/test_task7',*TASKS))]
    assert not legacy_changes,legacy_changes
    basefiles=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    preserved=0
    for name in basefiles:
        if name in allowed:continue
        # Git autocrlf may materialize working bytes; compare staged/index-free Git diff first.
        assert name not in changed,name
        preserved+=1
    py=[];figs=[];links=[];broken=[]
    docs=[ROOT/'README.md']+list((ROOT/'docs').rglob('*.md'))
    for folder in TASKS:
        for path in (ROOT/folder).rglob('*.py'):
            if path.relative_to(ROOT).as_posix() not in included:continue
            ast.parse(path.read_text(encoding='utf-8-sig'));py.append(str(path.relative_to(ROOT)))
        docs += [p for p in (ROOT/folder).rglob('*.md') if p.relative_to(ROOT).as_posix() in included]
        images=[p for p in (ROOT/folder).rglob('fig[1-4]_*.png') if p.relative_to(ROOT).as_posix() in included]
        assert len(images)==4,(folder,len(images))
        for path in images:
            with Image.open(path) as im:
                assert all(abs(d-300)<.02 for d in im.info.get('dpi',[])),path
                assert len(im.info.get('dpi',[]))==2,path
                im.verify()
            figs.append({'path':path.relative_to(ROOT).as_posix(),'sha256':sha(path),'dpi':[300,300]})
    for path in docs:
        # Original attachment may contain commands, placeholders and incomplete fences.
        if path.name.startswith('TASK') and path.parent.name=='task5_7_prompts':continue
        content=path.read_text(encoding='utf8')
        content=re.sub(r'```.*?```','',content,flags=re.S)
        for url in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',content):
            raw=url.strip().split(' "')[0].strip('<>')
            if not raw or raw.startswith(('#','http:','https:','mailto:')):continue
            target=path.parent/unquote(urlsplit(raw).path)
            if not target.exists():broken.append({'document':path.relative_to(ROOT).as_posix(),'target':raw})
            links.append(raw)
    assert not broken,broken
    manifest_checks=[]; reproduction_checks=[]
    manifests=(('task5_cyp_ddi/results','files'),('task6_asd','files'),('task7_qsp','files_sha256'))
    drivers=('run_task5_cyp_ddi_mechanism_based_inhibition.py','run_task6_asd_formulation_supersaturation_kinetics.py','run_task7_qsp_tumor_immune_pkpd_synergy.py')
    for i,((folder,key),driver) in enumerate(zip(manifests,drivers)):
        base=ROOT/folder
        manifest=json.loads((base/'manifest.json').read_text(encoding='utf8'))
        for name,expected in manifest[key].items():
            actual=sha(base/name)
            assert actual==expected,(folder,name,'frozen manifest mismatch')
            manifest_checks.append({'path':(base/name).relative_to(ROOT).as_posix(),'sha256':actual})
            if reproductions and name!='run_log.json':
                other=reproductions[i]/name
                assert other.exists() and sha(other)==actual,(folder,name,'independent rerun mismatch')
                reproduction_checks.append({'task':TASKS[i],'artifact':name,'byte_identical':True})
        source_info=manifest
        if i==0:
            source_info=json.loads((base/'results_summary.json').read_text(encoding='utf8'))
            assert source_info['task4_sha256']==sha(ROOT/'task4_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py')
        assert source_info['script_sha256']==sha(ROOT/TASKS[i]/driver),(folder,'source hash mismatch')
    result={'base_commit':BASE,'preexisting_files_preserved':preserved,'permitted_existing_metadata_changes':sorted(set(changed)&allowed),
            'new_python_sources_parsed':py,'figure_metadata':figs,'local_markdown_links_checked':len(links),'broken_links':broken,
            'manifest_artifacts_verified':manifest_checks,'independent_reproduction_comparisons':reproduction_checks,
            'excluded_from_reproduction_equality':['run_log.json: wall-clock timing differs','manifest.json: derived timing-log hash differs'],
            'scope':'Structural, provenance and optional byte-for-byte reproduction checks; numerical checks and visual review recorded separately'}
    (out/'integration.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reproductions',nargs=3,type=Path,metavar=('TASK5','TASK6','TASK7'))
    a=p.parse_args();inspect(a.out,a.reproductions)
