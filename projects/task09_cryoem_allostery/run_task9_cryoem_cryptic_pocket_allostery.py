"""Public-structure anchored SYNTHETIC ensemble, MSM, pocket geometry and PRS.

The 10 ns clock, state populations and interpolated coordinates are assumptions.
No particle images, cryoDRGN training, MD, binding validation or clinical claims.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import unittest
import numpy as np
from scipy import linalg, ndimage, spatial
from scipy.sparse.csgraph import dijkstra
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
R = 0.00198720425864083  # kcal mol-1 K-1
DEFAULT = dict(seed=20260920, n_conformers=500, temperature_K=298.15,
               lag_ns=10.0, assumed_state_free_energies_kcal_mol=[0,.35,.65,1.05,1.6],
               jump_probability=.35, reversible_edge_pseudocount=.5,
               bootstrap_replicates=200, grid_spacing_A=.75, roi_halfwidth_A=9.,
               pocket_probe_A=1., sasa_probe_A=1.4, sasa_sphere_points=96,
               density_voxel_A=1., density_psf_fwhm_A=2.5,
               synthetic_loop_displacement_A=.25, anm_cutoff_A=13.,
               prs_force_counts=[64,256,1024], active_site_residue=11)
RADII = {'C':1.70,'N':1.55,'O':1.52,'S':1.80}
Z = {'C':6,'N':7,'O':8,'S':16}

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n',encoding='utf-8')

def table(path, names, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f); writer.writerow(names); writer.writerows(rows)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_config(c):
    if set(c)!=set(DEFAULT): raise ValueError('Configuration keys must match the exported example.')
    for k in ['seed','n_conformers','bootstrap_replicates','sasa_sphere_points','active_site_residue']:
        if not isinstance(c[k],int) or isinstance(c[k],bool): raise ValueError(f'{k} must be integer')
    if not 100<=c['n_conformers']<=2000: raise ValueError('n_conformers must be 100..2000')
    if not 10<=c['bootstrap_replicates']<=2000: raise ValueError('bootstrap_replicates must be 10..2000')
    for k in ['temperature_K','lag_ns','grid_spacing_A','roi_halfwidth_A','pocket_probe_A',
              'sasa_probe_A','density_voxel_A','density_psf_fwhm_A','anm_cutoff_A']:
        if not np.isfinite(c[k]) or c[k]<=0: raise ValueError(f'{k} must be finite and positive')
    if not .4<=c['grid_spacing_A']<=1.5 or not 6<=c['roi_halfwidth_A']<=15: raise ValueError('Grid exceeds CPU budget')
    if not .1<c['jump_probability']<.9: raise ValueError('jump_probability must be .1.. .9')
    if not 24<=c['sasa_sphere_points']<=1024: raise ValueError('sasa_sphere_points must be 24..1024')
    if c['seed']<0: raise ValueError('seed must be nonnegative')
    if not .5<=c['density_voxel_A']<=2 or not 1<=c['density_psf_fwhm_A']<=6: raise ValueError('Density grid/PSF exceeds supported budget')
    if not .05<=c['reversible_edge_pseudocount']<=5: raise ValueError('positive regularization required')
    if not 0<=c['synthetic_loop_displacement_A']<=1: raise ValueError('loop displacement must be 0..1 A')
    if len(c['assumed_state_free_energies_kcal_mol'])!=5 or not np.isfinite(c['assumed_state_free_energies_kcal_mol']).all(): raise ValueError('Five finite state energies required')
    if len(c['prs_force_counts'])!=3 or any(not isinstance(n,int) or not 16<=n<=8192 for n in c['prs_force_counts']): raise ValueError('Three bounded integer force counts required')
    if c['prs_force_counts']!=sorted(set(c['prs_force_counts'])): raise ValueError('Force counts must be strictly increasing')

def read_pdb(path, hetero=False):
    atoms={}
    for line in path.read_text().splitlines():
        if line[:6].strip()!=('HETATM' if hetero else 'ATOM') or line[21]!='A': continue
        if line[16] not in ' A': continue  # coherent A alternate, not a per-atom occupancy mixture
        element=line[76:78].strip() or line[12:16].strip()[0]
        if element not in RADII: continue
        key=(int(line[22:26]),line[12:16].strip())
        atoms[key]=(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]),line[17:20].strip(),element)
    return atoms

def align_anchors():
    a=read_pdb(HERE/'inputs/4W51.pdb'); b=read_pdb(HERE/'inputs/4W59.pdb')
    keys=sorted(set(a)&set(b)); x=np.array([a[k][0] for k in keys]); y=np.array([b[k][0] for k in keys])
    ca=np.array([i for i,k in enumerate(keys) if k[1]=='CA'])
    # Fit the core; the rearranging 100-120 region is excluded from alignment.
    fit=[i for i in ca if not 100<=keys[i][0]<=120]
    cx=x[fit].mean(0); cy=y[fit].mean(0); u,_,vt=linalg.svd((y[fit]-cy).T@(x[fit]-cx))
    correction=np.diag([1.,1.,np.linalg.det(u@vt)]); rot=u@correction@vt
    y=(y-cy)@rot+cx
    lig=read_pdb(HERE/'inputs/4W59.pdb',True)
    ligand=np.array([v[0] for v in lig.values() if v[1]=='3GZ'])
    if not len(ligand): raise ValueError('Expected 4W59 3GZ ligand missing')
    center=((ligand-cy)@rot+cx).mean(0)
    return keys,x,y,ca,[a[k][1] for k in keys],np.array([a[k][2] for k in keys]),center

def generator(c):
    g=np.asarray(c['assumed_state_free_energies_kcal_mol']); pi=np.exp(-(g-g.min())/(R*c['temperature_K'])); pi/=pi.sum()
    t=np.zeros((5,5))
    for i in range(5):
        for j in [i-1,i+1]:
            if 0<=j<5: t[i,j]=.5*c['jump_probability']*min(1,pi[j]/pi[i])
        t[i,i]=1-t[i].sum()
    return t,pi

def simulate_states(t,n,rng,pi=None):
    states=np.zeros(n,dtype=int); states[0]=rng.choice(len(t),p=pi) if pi is not None else 0
    for i in range(1,n): states[i]=rng.choice(len(t),p=t[states[i-1]])
    return states

def estimate_msm(states, prior=.5, lag=1):
    counts=np.zeros((5,5)); np.add.at(counts,(states[:-lag],states[lag:]),1)
    support=np.eye(5)+np.eye(5,k=1)+np.eye(5,k=-1)
    # Symmetric flux/count estimator with explicit regularization, not reversible MLE.
    flux=(counts+counts.T)/2+prior*support
    mass=flux.sum(1); t=flux/mass[:,None]; pi=mass/mass.sum()
    return t,pi,counts

def mfpt(t,target=4,lag_ns=10.):
    keep=np.array([i for i in range(len(t)) if i!=target]); m=np.zeros(len(t))
    m[keep]=linalg.solve(np.eye(len(keep))-t[np.ix_(keep,keep)],np.ones(len(keep))*lag_ns)
    return m

def nonlinear_embedding(coords,ca):
    flat=coords[:,ca].reshape(len(coords),-1)
    dist=spatial.distance.squareform(spatial.distance.pdist(flat,'sqeuclidean'))/len(ca)
    bandwidth=np.median(dist[dist>0]); kernel=np.exp(-dist/(2*bandwidth))
    kernel-=kernel.mean(0)[None,:]+kernel.mean(1)[:,None]-kernel.mean()
    values,vectors=linalg.eigh(kernel,subset_by_index=[len(coords)-2,len(coords)-1])
    z=vectors[:,::-1]*np.sqrt(np.maximum(values[::-1],0))
    for j in range(2):
        if z[np.argmax(abs(z[:,j])),j]<0: z[:,j]*=-1
    return z,float(bandwidth)

def sphere_points(n):
    i=np.arange(n); z=1-2*(i+.5)/n; theta=i*np.pi*(3-np.sqrt(5))
    return np.column_stack([np.sqrt(1-z*z)*np.cos(theta),np.sqrt(1-z*z)*np.sin(theta),z])

def surface_clearance(tree,points,radii):
    k=min(8,len(radii)); d,ix=tree.query(points,k=k)
    if k==1: return d-radii[ix]
    gap=(d-radii[ix]).min(1)
    # Nearest center need not have the nearest vdW surface. Resolve every ambiguous point exactly.
    ambiguous=d[:,-1]-radii.max()<gap
    if ambiguous.any(): gap[ambiguous]=(spatial.distance.cdist(points[ambiguous],tree.data)-radii).min(1)
    return gap

def lining_sasa(x,radii,selected,probe,n):
    if not len(selected): return 0.
    directions=sphere_points(n); area=0.
    neighbors=spatial.cKDTree(x); rr=radii+probe
    for i in selected:
        pts=x[i]+rr[i]*directions
        js=[j for j in neighbors.query_ball_point(x[i],rr[i]+rr.max()) if j!=i]
        exposed=np.ones(n,dtype=bool) if not js else np.all(spatial.distance.cdist(pts,x[js])>=rr[js][None,:]-1e-10,axis=1)
        area+=4*np.pi*rr[i]**2*exposed.mean()
    return float(area)

def pocket_metrics(x,elements,resnames,center,c,spacing=None,with_sasa=True):
    h=spacing or c['grid_spacing_A']; half=c['roi_halfwidth_A']; axis=np.arange(-half,half+h*.1,h)
    shape=(len(axis),)*3; points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3)+center
    radii=np.array([RADII[e] for e in elements]); tree=spatial.cKDTree(x)
    clearance=surface_clearance(tree,points,radii).reshape(shape)
    occupied=clearance<c['pocket_probe_A']; enclosure=np.zeros(shape,dtype=np.uint8)
    for ax in range(3):
        enclosure+=(np.cumsum(occupied,axis=ax)>0)
        enclosure+=(np.flip(np.cumsum(np.flip(occupied,axis=ax),axis=ax),axis=ax)>0)
    candidate=(~occupied)&(enclosure>=5)
    labels,n=ndimage.label(candidate)
    if n:
        # Ligand-seeded component: select the eligible voxel nearest the bound-ligand centroid.
        idx=np.flatnonzero(candidate); nearest=idx[np.argmin(np.linalg.norm(points[idx]-center,axis=1))]
        component=labels==labels.ravel()[nearest]
    else: component=np.zeros(shape,dtype=bool)
    pocket=points[component.ravel()]; volume=float(len(pocket)*h**3)
    if len(pocket):
        ptree=spatial.cKDTree(pocket); d,_=ptree.query(x); selected=np.flatnonzero(d<radii+2.5)
        aromatic=float(np.mean(np.isin(np.array(resnames)[selected],['PHE','TYR','TRP','HIS'])))
        polarity=float(np.mean(np.isin(elements[selected],['N','O','S'])))
        hydrophobic=float(np.mean(elements[selected]=='C'))
        depth=float(ndimage.distance_transform_edt(component).max()*h)
        enc=float(enclosure[component].mean()/6)
        boundary=bool(component[0].any() or component[-1].any() or component[:,0].any() or component[:,-1].any() or component[:,:,0].any() or component[:,:,-1].any())
        raw=.49*np.log(volume)+.78*np.log(max(aromatic,1e-6))-.22*polarity-1.2
        score=float(1/(1+np.exp(-raw)))
        area=lining_sasa(x,radii,selected,c['sasa_probe_A'],c['sasa_sphere_points']) if with_sasa else None
    else:
        selected=np.array([],dtype=int); aromatic=polarity=hydrophobic=depth=enc=score=0.; raw=None; area=0.; boundary=False
    classification='indeterminate_roi_truncated' if boundary else ('meets_sampled_gate' if volume>500 and score>.7 else 'does_not_meet_sampled_gate')
    return dict(volume_A3=volume,full_component_volume_A3=None if boundary else volume,threshold_classification=classification,
                lining_sasa_A2=area,aromatic_residue_atom_fraction=aromatic,
                polar_atom_fraction=polarity,enclosure_fraction=enc,hydrophobic_enclosure_depth_A=depth*hydrophobic,
                raw_unbounded_score=raw,bounded_uncalibrated_score=score,roi_boundary_contact=boundary),selected

def anm_hessian(x,cutoff):
    n=len(x); h=np.zeros((3*n,3*n)); distance=spatial.distance.squareform(spatial.distance.pdist(x))
    for i,j in zip(*np.where(np.triu((distance<cutoff)&(distance>0),1))):
        d=x[j]-x[i]; block=np.outer(d,d)/(d@d)
        si=slice(3*i,3*i+3); sj=slice(3*j,3*j+3)
        h[si,si]+=block; h[sj,sj]+=block; h[si,sj]-=block; h[sj,si]-=block
    return h,distance

def prs_calculation(x,cutoff,force_counts,seed):
    h,distance=anm_hessian(x,cutoff); values,vectors=linalg.eigh(h)
    positive=values>max(values[-1]*1e-9,1e-12); zero=int((~positive).sum())
    if zero!=6: raise ValueError(f'ANM has {zero} null modes; six rigid modes required')
    pinv=(vectors[:,positive]/values[positive])@vectors[:,positive].T
    n=len(x); blocks=pinv.reshape(n,3,n,3).transpose(0,2,1,3)
    exact=(blocks**2).sum((2,3))/3
    rng=np.random.default_rng(seed); f=rng.normal(size=(n,max(force_counts),3)); f/=np.linalg.norm(f,axis=2)[:,:,None]
    convergence=[]; sampled=None
    for count in force_counts:
        cov=np.einsum('jfa,jfb->jab',f[:,:count],f[:,:count])/count
        sampled=np.einsum('ijka,jab,ijkb->ij',blocks,cov,blocks,optimize=True)
        convergence.append(dict(forces_per_residue=count,relative_frobenius_error=float(np.linalg.norm(sampled-exact)/np.linalg.norm(exact))))
    normalized=exact/np.diag(exact)[None,:]
    coupling=exact/np.sqrt(np.diag(exact)[:,None]*np.diag(exact)[None,:])
    return dict(hessian=h,pseudoinverse=pinv,eigenvalues=values,exact=exact,sampled=sampled,
                normalized=normalized,coupling=coupling,distance=distance,zero_modes=zero,
                convergence=convergence,pseudoinverse_residual=float(np.linalg.norm(h@pinv@h-h)/np.linalg.norm(h)))

def coupling_path(prs,source_indices,target,cutoff):
    d=prs['distance']; coupling=np.clip(prs['coupling'],1e-12,1)
    graph=np.where((d<cutoff)&(d>0),-np.log(coupling)+.01*d,0)
    distances,predecessor=dijkstra(graph,directed=False,indices=target,return_predecessors=True)
    candidates=[int(i) for i in source_indices if i!=target and np.isfinite(distances[i])]
    if not candidates: raise ValueError('No cavity-to-active-site path')
    source=min(candidates,key=lambda i:distances[i]); path=[source]
    while path[-1]!=target:
        path.append(int(predecessor[path[-1]]))
        if path[-1]<0 or len(path)>len(d): raise ValueError('Broken shortest path')
    return path,float(sum(d[i,j] for i,j in zip(path[:-1],path[1:]))),float(d[source,target])

def density_maps(coords,elements,states,c):
    h=c['density_voxel_A']; sigma=c['density_psf_fwhm_A']/(2*np.sqrt(2*np.log(2)))/h
    lower=np.floor(coords.min((0,1))-5); upper=np.ceil(coords.max((0,1))+5)
    shape=np.ceil((upper-lower)/h).astype(int)+1; maps=np.zeros((5,*shape),dtype=np.float64)
    weights=np.array([Z[e] for e in elements]); individual_peak=[]
    for xyz,s in zip(coords,states):
        indices=np.rint((xyz-lower)/h).astype(int); volume=np.zeros(shape)
        np.add.at(volume,tuple(indices.T),weights/h**3)
        volume=ndimage.gaussian_filter(volume,sigma=sigma,mode='constant')
        maps[s]+=volume; individual_peak.append(float(volume.max()))
    counts=np.bincount(states,minlength=5)
    if (counts==0).any(): raise ValueError('Density state has no sampled conformer; select a better sampled simulation')
    average=maps.sum(0)/len(coords); maps/=counts[:,None,None,None]
    return maps.astype(np.float32),average.astype(np.float32),lower,h,float(np.mean(individual_peak))

def save_figures(out,z,states,q,pockets,maps,average,lower,h,center,msm,prs,ca_xyz,resids,path,c):
    figdir=out/'figures_task9'; figdir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    fig=plt.figure(figsize=(12,8)); gs=fig.add_gridspec(2,3)
    ax=fig.add_subplot(gs[:,0]); scat=ax.scatter(z[:,0],z[:,1],c=states,cmap='viridis',s=15,alpha=.7)
    ax.set(xlabel='Kernel PCA coordinate 1',ylabel='Kernel PCA coordinate 2',title=f'{len(states)} synthetic coordinate snapshots')
    fig.colorbar(scat,ax=ax,label='Assigned geometric macrostate',ticks=range(5),shrink=.6)
    cut=int(np.clip(round((center[2]-lower[2])/h),0,average.shape[2]-1))
    extent=[lower[0],lower[0]+h*average.shape[0],lower[1],lower[1]+h*average.shape[1]]
    panels=[(maps[0],'Closed-bin mean'),(maps[4],'Open-bin mean'),(average,'Population-weighted ensemble'),(abs(maps[4]-maps[0]),'Absolute open-minus-closed difference')]
    shared_max=max(float(grid[:,:,cut].max()) for grid,_ in panels[:3])
    for k,(grid,title) in enumerate(panels):
        ax=fig.add_subplot(gs[k//2,k%2+1]); im=ax.imshow(grid[:,:,cut].T,origin='lower',extent=extent,cmap='magma',aspect='equal',vmin=0,vmax=shared_max if k<3 else None)
        fig.colorbar(im,ax=ax,shrink=.65,label='Z-weighted proxy / Å³' if k<3 else '|Δ proxy| / Å³')
        ax.scatter(*center[:2],s=24,facecolors='none',edgecolors='cyan'); ax.set(title=title,xlabel='x (Å)',ylabel='y (Å)')
    fig.suptitle(f"Structure-anchored synthetic density: {c['density_psf_fwhm_A']} Å Gaussian PSF FWHM; {c['density_voxel_A']} Å voxels\nNo cryo-EM particles, CTF, reconstruction or cryoDRGN training",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.94)); fig.savefig(figdir/'fig1_cryoem_latent_conformational_manifold.png',dpi=300,bbox_inches='tight',pad_inches=.15); plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(14,4.8)); g=np.array(msm['estimated_free_energies_kcal_mol']); lo,hi=np.array(msm['bootstrap_free_energy_95_interval']).T
    axs[0].plot(range(5),c['assumed_state_free_energies_kcal_mol'],'--',color='gray',label='Assumed generator')
    axs[0].errorbar(range(5),g,yerr=[np.maximum(g-lo,0),np.maximum(hi-g,0)],marker='o',capsize=4,label='Regularized count estimate')
    axs[0].set(xlabel='Macrostate: 0 closed → 4 open',ylabel='Population free energy (kcal/mol)',title='State free energies are not TS barriers'); axs[0].legend(fontsize=8)
    im=axs[1].imshow(msm['estimated_transition_matrix'],vmin=0,vmax=1,cmap='Blues'); fig.colorbar(im,ax=axs[1],label=f"Transition probability / assumed {c['lag_ns']} ns")
    axs[1].set(xlabel='Destination state',ylabel='Source state',title='Reversible five-state model')
    for i in range(5):
        for j in range(5): axs[1].text(j,i,f"{msm['estimated_transition_matrix'][i][j]:.2f}",ha='center',va='center',fontsize=8)
    axs[2].bar(range(5),msm['observed_counts'],color=plt.cm.viridis(np.linspace(0,1,5))); axs[2].set(xlabel='Macrostate',ylabel='Sampled snapshots',title='Finite sampling and assumed clock')
    axs[2].text(.04,.97,f"ΔG(open−closed) = {g[-1]:.3f} kcal/mol\nMFPT(0→4) = {msm['mfpt_closed_to_open_ns']:.1f} ns\nk = {msm['opening_rate_per_ns']:.5f} ns⁻¹\nConditional on the synthetic generator",transform=axs[2].transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='gray',alpha=.97,pad=5))
    fig.suptitle('Regularization and parametric-bootstrap uncertainty are explicit; no target energy was imposed on the estimate')
    fig.tight_layout(); fig.savefig(figdir/'fig2_msm_free_energy_pathway.png',dpi=300,bbox_inches='tight',pad_inches=.15); plt.close(fig)

    vol=np.array([p['volume_A3'] for p in pockets]); score=np.array([p['bounded_uncalibrated_score'] for p in pockets]); sasa=np.array([p['lining_sasa_A2'] for p in pockets]); order=np.argsort(q)
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    clipped=np.array([p['roi_boundary_contact'] for p in pockets])
    axs[0,0].plot(q[order],vol[order],lw=.8); axs[0,0].scatter(q[clipped],vol[clipped],s=9,color='orange',label='ROI truncated: full volume unknown'); axs[0,0].axhline(500,ls='--',color='gray',label='Requested 500 Å³ gate'); axs[0,0].set(xlabel='Interpolation coordinate q',ylabel='ROI cavity volume (Å³)',title='Geometric volume with censoring'); axs[0,0].legend(fontsize=7)
    axs[0,1].plot(q[order],score[order],color='#d17c24'); axs[0,1].axhline(.7,ls='--',color='gray'); axs[0,1].set(xlabel='Interpolation coordinate q',ylabel='Uncalibrated logistic score',ylim=(0,1),title='A bounded heuristic is not binding evidence')
    axs[1,0].plot(np.arange(len(q))*c['lag_ns']/1000,vol,lw=.8); axs[1,0].set(xlabel='Assumed synthetic time (µs)',ylabel='Cavity volume (Å³)',title='Sample ordering is supplied by the generator')
    axs[1,1].scatter(vol,sasa,c=q,cmap='viridis',s=10); axs[1,1].set(xlabel='Cavity volume (Å³)',ylabel='Lining probe-accessible area (Å²)',title='Ligand removed; waters and ions excluded')
    fig.tight_layout(); fig.savefig(figdir/'fig3_dynamic_pocket_volume_druggability.png',dpi=300,bbox_inches='tight',pad_inches=.15); plt.close(fig)

    fig=plt.figure(figsize=(13,5.5)); ax=fig.add_subplot(121)
    im=ax.imshow(np.log10(np.maximum(prs['normalized'],1e-5)),origin='lower',cmap='magma',vmin=-3,vmax=0); fig.colorbar(im,ax=ax,label='log10 response / source self-response')
    ax.set(xlabel='Perturbed residue index',ylabel='Responding residue index',title='ANM pseudoinverse: isotropic PRS')
    ax=fig.add_subplot(122,projection='3d'); ax.plot(*ca_xyz.T,color='silver',lw=.8,alpha=.8)
    xyz=ca_xyz[path]; ax.plot(*xyz.T,color='#d55e00',marker='o',lw=3,markersize=4)
    for idx in path: ax.text(*ca_xyz[idx],str(resids[idx]),fontsize=8)
    ax.scatter(*ca_xyz[path[0]],color='navy',s=70,label='Cavity lining source'); ax.scatter(*ca_xyz[path[-1]],color='green',s=70,label=f'Residue {c["active_site_residue"]} target')
    ax.set(xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)',title='Contact-constrained coupling path\nMechanical hypothesis; no energy-flow validation'); ax.legend(fontsize=8,loc='lower left')
    fig.tight_layout(); fig.savefig(figdir/'fig4_prs_allosteric_network_matrix.png',dpi=300,bbox_inches='tight',pad_inches=.2); plt.close(fig)

def run(out,c):
    started=time.perf_counter(); validate_config(c); out.mkdir(parents=True,exist_ok=True); results=out/'results'; results.mkdir(exist_ok=True)
    for name in ['4W51.pdb','4W59.pdb']:
        dst=out/'inputs'/name; dst.parent.mkdir(exist_ok=True)
        if dst.resolve()!=(HERE/'inputs'/name).resolve(): shutil.copyfile(HERE/'inputs'/name,dst)
    source_provenance=HERE/'inputs/provenance.json'
    if source_provenance.exists() and out.resolve()!=HERE.resolve(): shutil.copyfile(source_provenance,out/'inputs/provenance.json')
    dump(out/'inputs/config.json',c)
    keys,x,y,ca,resnames,elements,center=align_anchors(); rng=np.random.default_rng(c['seed']); t_true,pi_true=generator(c)
    states=simulate_states(t_true,c['n_conformers'],rng,pi_true); q=(states+rng.uniform(.001,.999,len(states)))/5
    loop=np.array([max(0,1-abs(k[0]-110)/12) for k in keys]); directions=rng.normal(size=(len(q),3)); directions/=np.linalg.norm(directions,axis=1)[:,None]
    amplitude=c['synthetic_loop_displacement_A']*np.sin(np.pi*q)
    coords=x[None]+q[:,None,None]*(y-x)[None]+amplitude[:,None,None]*loop[None,:,None]*directions[:,None,:]
    z,bandwidth=nonlinear_embedding(coords,ca)
    np.savez_compressed(results/'synthetic_conformers.npz',coordinates_A=coords.astype(np.float32),interpolation_q=q,states=states,ca_indices=ca,residue_ids=np.array([k[0] for k in keys]),atom_names=np.array([k[1] for k in keys]),elements=elements,latent_kernel_pca=z)
    table(results/'conformer_latent_trajectory.csv',['frame','assumed_time_ns','state','interpolation_q','z1','z2'],[(i,i*c['lag_ns'],int(states[i]),q[i],*z[i]) for i in range(len(q))])
    t,pi,counts=estimate_msm(states,c['reversible_edge_pseudocount']); free=-R*c['temperature_K']*np.log(pi/pi[0]); opening=mfpt(t,lag_ns=c['lag_ns'])[0]
    bootstrap=[]; unseen=0
    for _ in range(c['bootstrap_replicates']):
        s=simulate_states(t,len(states),rng,pi); tb,pib,_=estimate_msm(s,c['reversible_edge_pseudocount'])
        unseen+=int(len(np.unique(s))<5); gb=-R*c['temperature_K']*np.log(pib/pib[0]); bootstrap.append([*gb,mfpt(tb,lag_ns=c['lag_ns'])[0]])
    bootstrap=np.array(bootstrap); t2,_,_=estimate_msm(states,c['reversible_edge_pseudocount'],lag=2)
    msm=dict(estimator='symmetrized transition counts plus explicit diagonal/nearest-neighbor pseudocount; not reversible MLE',clock_status='assumed synthetic; no MD or experimental chronology',
             generator_transition_matrix=t_true.tolist(),generator_stationary=pi_true.tolist(),observed_counts=np.bincount(states,minlength=5).tolist(),transition_counts=counts.tolist(),
             estimated_transition_matrix=t.tolist(),estimated_stationary=pi.tolist(),estimated_free_energies_kcal_mol=free.tolist(),
             bootstrap_free_energy_95_interval=np.quantile(bootstrap[:,:5],[.025,.975],axis=0).T.tolist(),
             mfpt_closed_to_open_ns=float(opening),opening_rate_per_ns=float(1/opening),bootstrap_mfpt_95_interval_ns=np.quantile(bootstrap[:,5],[.025,.975]).tolist(),
             bootstrap_replicates_missing_at_least_one_state=unseen,bootstrap_replicates=c['bootstrap_replicates'],
             detailed_balance_max_error=float(np.max(abs(pi[:,None]*t-pi[None,:]*t.T))),stationarity_max_error=float(np.max(abs(pi@t-pi))),
             ck_lag2_frobenius_disagreement=float(np.linalg.norm(t@t-t2)),opening_penalty_meets_requested_2_to_5=bool(2<=free[-1]<=5))
    dump(results/'msm.json',msm); table(results/'msm_bootstrap.csv',['G0','G1','G2','G3','G4','mfpt_ns'],bootstrap)
    print(f'MSM and {len(coords)}-coordinate ensemble complete',flush=True)
    pockets=[]
    for i,xyz in enumerate(coords):
        p,_=pocket_metrics(xyz,elements,resnames,center,c); pockets.append(p)
        if (i+1)%100==0: print(f'Pocket/SASA {i+1}/{len(coords)}',flush=True)
    names=list(pockets[0]); table(results/'pocket_geometry.csv',['frame','q',*names],[[i,q[i],*[p[n] for n in names]] for i,p in enumerate(pockets)])
    grid_check=[]
    for label,xyz in [('closed_anchor',x),('open_anchor',y)]:
        for spacing in [.5,.75,1.]:
            p,_=pocket_metrics(xyz,elements,resnames,center,c,spacing,False); grid_check.append(dict(anchor=label,spacing_A=spacing,**p))
    dump(results/'grid_convergence.json',grid_check)
    roi_check=[]
    for label,xyz in [('closed_anchor',x),('mid_interpolation',(x+y)/2),('open_anchor',y)]:
        for halfwidth in [9.,12.,15.]:
            expanded=dict(c); expanded['roi_halfwidth_A']=halfwidth
            p,_=pocket_metrics(xyz,elements,resnames,center,expanded,with_sasa=False)
            roi_check.append(dict(anchor=label,roi_halfwidth_A=halfwidth,**p))
    dump(results/'roi_sensitivity.json',roi_check)
    print('Pocket grid convergence complete',flush=True)
    maps,average,lower,h,peak=density_maps(coords,elements,states,c)
    np.savez_compressed(results/'synthetic_density_maps.npz',state_mean_density=maps,ensemble_mean_density=average,origin_A=lower,voxel_A=h,psf_fwhm_A=c['density_psf_fwhm_A'],state_counts=np.bincount(states,minlength=5))
    prs=prs_calculation(x[ca],c['anm_cutoff_A'],c['prs_force_counts'],c['seed']+1)
    _,lining=pocket_metrics(y,elements,resnames,center,c,with_sasa=False); lining_resids=set(keys[i][0] for i in lining)
    resids=[keys[i][0] for i in ca]; source_indices=[i for i,rid in enumerate(resids) if rid in lining_resids]; target=resids.index(c['active_site_residue'])
    path,path_length,straight=coupling_path(prs,source_indices,target,c['anm_cutoff_A'])
    np.savez_compressed(results/'anm_prs.npz',hessian=prs['hessian'],pseudoinverse=prs['pseudoinverse'],eigenvalues=prs['eigenvalues'],prs_isotropic=prs['exact'],prs_random_force=prs['sampled'],prs_normalized=prs['normalized'],coupling=prs['coupling'],residue_ids=resids,ca_coordinates_A=x[ca])
    table(results/'allosteric_path.csv',['path_order','residue_id','x_A','y_A','z_A'],[[i,resids[j],*x[ca][j]] for i,j in enumerate(path)])
    neighbor_dist=np.linalg.norm(np.diff(coords[:,ca],axis=1),axis=2)
    # Exclude sequence neighbors; this is an inexpensive coarse geometric screen, not energy minimization.
    close=[]
    for xyz in coords[:,ca]:
        d=spatial.distance.squareform(spatial.distance.pdist(xyz)); d[np.abs(np.subtract.outer(np.arange(len(ca)),np.arange(len(ca))))<=1]=np.inf; close.append(float(d.min()))
    gate=[i for i,p in enumerate(pockets) if p['threshold_classification']=='meets_sampled_gate']
    first=min(gate,key=lambda i:q[i]) if gate else None
    summary=dict(evidence='public X-ray endpoints plus synthetic interpolation and assumed Markov time; not experimental cryo-EM/MD/drug discovery',
                 conformers=len(coords),heavy_atoms=len(keys),residues=len(ca),kernel_pca_bandwidth_A2=bandwidth,
                 anchor_ca_rmsd_A=float(np.sqrt(np.mean(np.sum((x[ca]-y[ca])**2,axis=1)))),
                 volume_range_A3=[min(p['volume_A3'] for p in pockets),max(p['volume_A3'] for p in pockets)],
                 uncalibrated_score_range=[min(p['bounded_uncalibrated_score'] for p in pockets),max(p['bounded_uncalibrated_score'] for p in pockets)],
                 high_volume_high_score_frames=len(gate),first_sampled_gate_q=None if first is None else float(q[first]),
                 threshold_indeterminate_frames=sum(p['threshold_classification']=='indeterminate_roi_truncated' for p in pockets),
                 exact_threshold_status='not established; finite sampled geometries and unvalidated score',
                 roi_boundary_contact_frames=sum(p['roi_boundary_contact'] for p in pockets),
                 ca_adjacent_distance_range_A=[float(neighbor_dist.min()),float(neighbor_dist.max())],
                 minimum_nonadjacent_ca_distance_A=min(close),
                 anm_zero_modes=prs['zero_modes'],anm_pseudoinverse_residual=prs['pseudoinverse_residual'],prs_force_convergence=prs['convergence'],
                 path_residue_ids=[resids[i] for i in path],path_length_A=path_length,path_end_to_end_A=straight,path_end_to_end_above_30_A=straight>30,
                 mean_individual_density_peak=peak,ensemble_density_peak=float(average.max()),
                 density_units='atomic-number weighted Gaussian proxy / A^3; not calibrated electrostatic potential',
                 density_population_weighting='actual synthetic state counts / total frames',msm=msm)
    dump(results/'summary.json',summary)
    save_figures(out,z,states,q,pockets,maps,average,lower,h,center,msm,prs,x[ca],resids,path,c)
    write_reports(out,c,summary)
    verification=dict(all_finite_coordinates=bool(np.isfinite(coords).all()),probabilities_normalized=bool(np.allclose(t.sum(1),1)),detailed_balance=msm['detailed_balance_max_error']<1e-12,
                      six_rigid_modes=prs['zero_modes']==6,pseudoinverse_identity=prs['pseudoinverse_residual']<1e-8,
                      raw_formula_not_claimed_bounded=True,synthetic_clock_explicit=True,random_force_relative_error_below_0_12=prs['convergence'][-1]['relative_frobenius_error']<.12)
    dump(out/'verification.json',verification)
    dump(out/'run_log.json',dict(runtime_seconds=time.perf_counter()-started,python=sys.version,numpy=np.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__))
    refresh_manifest(out)
    print(json.dumps(summary,indent=2),flush=True)
    if not all(verification.values()): raise RuntimeError('A numerical verification gate failed')
    return summary

def refresh_manifest(out):
    files={str(p.relative_to(out)).replace('\\','/'):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    dump(out/'manifest.json',dict(schema_version=1,script_sha256=sha(Path(__file__)),files_sha256=files))

def write_reports(out,c,s):
    m=s['msm']; g=m['estimated_free_energies_kcal_mol'][-1]; lo,hi=m['bootstrap_free_energy_95_interval'][-1]
    common=f"""| Quantity / 指标 | Executed result / 已执行结果 |
|---|---:|
| Conformers / 合成构象 | {s['conformers']} |
| Common heavy atoms / 匹配重原子 | {s['heavy_atoms']} |
| Cα residues / 残基 | {s['residues']} |
| Core-aligned Cα endpoint RMSD / 端点偏差 | {s['anchor_ca_rmsd_A']:.6f} Å |
| State counts 0–4 / 状态计数 | {m['observed_counts']} |
| Estimated open-minus-closed population ΔG | {g:.6f} kcal/mol |
| Parametric-bootstrap 95% interval | [{lo:.6f}, {hi:.6f}] kcal/mol |
| Assumed-clock closed-to-open MFPT | {m['mfpt_closed_to_open_ns']:.6f} ns |
| Conditional opening rate 1/MFPT | {m['opening_rate_per_ns']:.9f} ns⁻¹ |
| MFPT bootstrap 95% interval | {m['bootstrap_mfpt_95_interval_ns']} ns |
| Bootstrap draws missing a state / 缺态重抽样 | {m['bootstrap_replicates_missing_at_least_one_state']}/{m['bootstrap_replicates']} |
| Cavity volume range / 空腔体积 | {s['volume_range_A3']} Å³ |
| Uncalibrated bounded score range / 未校准评分 | {s['uncalibrated_score_range']} |
| Joint volume >500 Å³ and score >0.7 / 同时达阈值 | {s['high_volume_high_score_frames']} frames |
| ROI boundary-contact frames / 触边帧 | {s['roi_boundary_contact_frames']} |
| Indeterminate joint-gate frames / 联合阈值未知帧 | {s['threshold_indeterminate_frames']} |
| Adjacent Cα distance range / 相邻 Cα 距离 | {s['ca_adjacent_distance_range_A']} Å |
| Minimum nonadjacent Cα distance / 非相邻最近距 | {s['minimum_nonadjacent_ca_distance_A']:.6f} Å |
| ANM null modes / 零模 | {s['anm_zero_modes']} |
| Relative H H⁺ H − H norm / 伪逆残差 | {s['anm_pseudoinverse_residual']:.3g} |
| PRS force direction convergence / 相对误差 | {[(v['forces_per_residue'],round(v['relative_frobenius_error'],6)) for v in s['prs_force_convergence']]} |
| Contact path residue IDs / 耦合路径 | {' → '.join(map(str,s['path_residue_ids']))} |
| Path contour length / 路径长度 | {s['path_length_A']:.6f} Å |
| Endpoint separation / 端点直线距离 | {s['path_end_to_end_A']:.6f} Å |
"""
    links="""
## Primary sources / 原始来源

- [4W51 apo coordinates](https://www.rcsb.org/structure/4W51) and [4W59 n-hexylbenzene-bound coordinates](https://www.rcsb.org/structure/4W59); [Merski et al., 2015](https://doi.org/10.1073/pnas.1500806112). These are X-ray structures of T4 lysozyme L99A.
- [Zhong et al., cryoDRGN, 2021](https://doi.org/10.1038/s41592-020-01049-4): methodological context only; its neural model is not executed here.
- [Trendelkamp-Schroer et al., reversible MSM estimation, 2015](https://arxiv.org/abs/1507.05990): reference for reversibility and uncertainty; this code uses a simpler explicitly stated estimator.
- [Atilgan et al., ANM, 2001](https://doi.org/10.1016/S0006-3495(01)76033-X) and [Atilgan & Atilgan, PRS, 2009](https://doi.org/10.1371/journal.pcbi.1000544).
- [Kuroki et al., T4 lysozyme catalytic-site study](https://pmc.ncbi.nlm.nih.gov/articles/PMC17713/): supports the Glu11 active-site reference; does not validate the computed pathway.
- [wwPDB/RCSB data usage policy](https://www.rcsb.org/pages/usage-policy): deposited PDB coordinate data are CC0; retain structure attribution.

## Inspectable outputs / 可检查产物

- [Coordinate archive](results/synthetic_conformers.npz), [ordered latent trajectory](results/conformer_latent_trajectory.csv), [density voxel maps](results/synthetic_density_maps.npz).
- [MSM estimate](results/msm.json), [bootstrap replicates](results/msm_bootstrap.csv), [pocket geometry](results/pocket_geometry.csv), [grid convergence](results/grid_convergence.json), [ROI sensitivity](results/roi_sensitivity.json).
- [ANM/PRS matrices](results/anm_prs.npz), [path coordinates](results/allosteric_path.csv), [summary](results/summary.json), [numerical verification](verification.json), [byte manifest](manifest.json).
- [Figure 1](figures_task9/fig1_cryoem_latent_conformational_manifold.png), [Figure 2](figures_task9/fig2_msm_free_energy_pathway.png), [Figure 3](figures_task9/fig3_dynamic_pocket_volume_druggability.png), [Figure 4](figures_task9/fig4_prs_allosteric_network_matrix.png).
"""
    en=f"""# Task 9 — Public-structure anchored synthetic heterogeneity and mechanical allostery

This executable benchmark combines two public X-ray endpoints with synthetic coordinate interpolation, Gaussian scattering proxies, a five-state Markov generator, geometric cavity measurements and an ANM/PRS calculation. It does **not** establish a new druggable pocket, a cancer-target mechanism, a physical kinetic rate, or a validated allosteric route. The model protein is bacteriophage T4 lysozyme L99A, not a human kinase or GPCR.

## Executed results

{common}

No sampled, fully contained component reaches the requested joint cavity gate; {s['threshold_indeterminate_frames']} frames have unknown full volume and indeterminate gate classification because they touch the ROI boundary. This does not exclude a larger or druggable pocket. The contact path does not span 30 Å end to end in the archived default run. Only {m['observed_counts'][-1]} open frames were sampled; the generator's declared opening penalty is {c['assumed_state_free_energies_kcal_mol'][-1]-c['assumed_state_free_energies_kcal_mol'][0]} kcal/mol. Sampling noise and explicit pseudocount regularization explain the difference from the estimate. No seed, state energy or score calibration was tuned to force the requested outcome. The estimated ΔG is a **state population difference**, not an activation barrier or a transition-state free energy.

## 9A. Structures, synthetic coordinates and density

Chain A heavy atoms common to 4W51 and 4W59 are matched by residue number and atom name. Alternate A is selected coherently; alternate B, hydrogens, waters, ions and all ligands are excluded from the protein calculation. The bound ligand 3GZ supplies only a cavity-center coordinate and is then removed. Core Cα atoms outside residues 100–120 define a proper Kabsch rotation. The 4W59 alternate-A endpoint is ligand-conditioned; it does not establish a spontaneously populated apo-open structure.

For each generated state s, q=(s+u)/5 with u uniform on (0.001,0.999). Coordinates are (1−q)x_closed+q x_bound plus a tapered random loop displacement, largest around residue 110. The loop amplitude is {c['synthetic_loop_displacement_A']} Å times sin(πq). This produces continuous geometric variability and an explicit ordering, but it is neither an atomistic dynamical integrator nor an energy-minimized transition pathway. Adjacent Cα compression to the reported minimum reveals an interpolation artifact; the coarse contact check does not test covalent stereochemistry or all heavy-atom clashes.

Kernel PCA uses the centered RBF kernel of aligned Cα Cartesian distances. Its two latent axes are geometry-derived nonlinear coordinates, not learned cryoDRGN embeddings. Atomic-number weighted coordinate deposits are Gaussian-filtered on a {c['density_voxel_A']} Å grid, with **{c['density_psf_fwhm_A']} Å FWHM** point spread. Pixel spacing and PSF width are different quantities. There are no electron micrographs, orientation inference, CTF, detector noise, FSC resolution estimates or Coulomb-potential calibration. Maps average the actual generated snapshots, with state-count weighting for the global mean. Density smoothing/attenuation is an illustrative consequence of coordinate averaging; no experimentally observed loop disappearance is claimed.

## 9B. Reversible synthetic MSM and uncertainty

The declared generator uses nearest-neighbor Metropolis transitions with proposed step probability {c['jump_probability']}/2 per direction and Boltzmann weights from the configured state energies at {c['temperature_K']} K. Reflecting endpoints are implemented by rejected proposals. The initial state is drawn from the generator's stationary distribution. One saved sample is assigned **{c['lag_ns']} ns by assumption**. These time labels cannot be recovered from unordered cryo-EM snapshots.

From {s['conformers']-1} adjacent-frame transitions, the symmetric flux estimator is F=(C+Cᵀ)/2+αS, α={c['reversible_edge_pseudocount']}, where S includes diagonal and nearest-neighbor entries. T_ij=F_ij/Σ_jF_ij and π_i=Σ_jF_ij/Σ_ijF_ij exactly satisfy detailed balance. This is a regularized moment estimator, not a reversible maximum-likelihood fit. The PMF is −RT ln(π_i/π_0). MFPT solves (I−T_nonabsorbing)m=τ1 with state 4 absorbing; k_open=1/m_0 is a conditional first-passage rate, not an elementary transition rate.

{c['bootstrap_replicates']} parametric replicate trajectories use the fitted T and are refitted with the same prior. Percentile intervals measure finite-sampling variation under this synthetic model; they exclude force-field, structural and clock uncertainty. Missing states in bootstrap draws remain visible and are regularized, never silently discarded. The 2τ Chapman–Kolmogorov Frobenius discrepancy is {m['ck_lag2_frobenius_disagreement']:.6f}; this is a finite-data diagnostic, not a successful physical Markovianity test.

## 9C. Geometry, area and the score boundary

The ligand-centered cubic ROI has half-width {c['roi_halfwidth_A']} Å. A voxel is accessible when an assumed {c['pocket_probe_A']} Å probe clears all vdW spheres (C/N/O/S radii 1.70/1.55/1.52/1.80 Å). At least five of six axis directions must encounter a protein obstacle within the ROI. The connected eligible component nearest the ligand centroid defines the seeded cavity. Volume is occupied-voxel count times spacing cubed. Thus this is a local sphere-filling/occlusion algorithm, not a global pocket detector, alpha-shape method or validated commercial score.

Lining atoms lie within their vdW radius +2.5 Å of cavity voxels. Fibonacci-sphere exposure with {c['sasa_sphere_points']} points and a {c['sasa_probe_A']} Å probe supplies lining probe-accessible surface area. Enclosed surfaces are included; access from bulk solvent is not proven. Aromaticity is the fraction of lining atoms belonging to Phe/Tyr/Trp/His residues, **not** a measured aromatic surface fraction. Polarity is the N/O/S atom fraction. Enclosure is the mean occluded-direction fraction; hydrophobic depth is maximum component distance-transform radius multiplied by the lining-carbon fraction, an explicitly geometric proxy.

The supplied formula is dimensionally interpreted as raw=0.49 ln[V/(1 Å³)]+0.78 ln(aromaticity)−0.22 polarity−1.2, with a 10⁻⁶ aromaticity floor. It is unbounded. We retain raw and additionally export sigmoid(raw), a bounded **uncalibrated heuristic** with no binding-probability interpretation. Empty cavities receive zero score and null raw. No exact transition threshold can be inferred from these samples. Boundary-contact volumes are ROI-limited descriptions of the detected component and cannot establish a closed pocket volume; full volume is null and classification indeterminate. The 0.5/0.75/1.0 Å endpoint audit measures voxel sensitivity. A separate 9/12/15 Å half-width audit evaluates closed/mid/open representative structures. Since changing ROI can also change directional enclosure and connectivity, a numerically larger component is not automatically the same physical cavity. The per-frame censoring remains, and fine-grid agreement is not global-cavity convergence.

## 9D. ANM, PRS and a contact-constrained hypothesis

The closed-endpoint Cα network uses unit springs and a {c['anm_cutoff_A']} Å contact cutoff. Each pair adds uuᵀ diagonal blocks and −uuᵀ off-diagonal blocks to the 3N×3N Hessian. Six rigid translation/rotation modes are removed with a relative eigenvalue tolerance; all positive internal modes contribute to H⁺. The singular Hessian is never inverted directly. Absolute displacements have arbitrary spring/force units.

For a unit isotropic force at residue j, M_ij=||H⁺_ij||²_F/3. Monte Carlo force sets of {c['prs_force_counts']} independent directions per source are evaluated via their exact empirical covariance, algebraically identical to averaging squared individual displacement responses. Source self-response normalizes each column. Symmetric coupling uses M_ij/sqrt(M_ii M_jj). Contact edges cost −ln(coupling)+0.01 distance/Å, an explicitly chosen routing metric. Dijkstra selects the cheapest route from any detected lining residue to reference active-site residue {c['active_site_residue']}; this endpoint selection can favor a nearby lining residue. The resulting path is not energy transport, causality, mutational validation or proof of long-range allostery.

## Industry relevance and next acceptance gates

This reusable calculation can expose pipeline failure modes before spending resources on structure-enabled hit discovery: conformational support, cavity truncation, score calibration and pathway sensitivity are separate questions. A real program needs target-specific cryo-EM particles or validated MD, state population/kinetic measurements, side-chain and solvent refinement, unbiased pocket detection, ligandability benchmarks, multiple-structure ANM sensitivity and perturbation experiments. Only then should predicted pockets prioritize fragment screening or medicinal chemistry. The archived negative gates are retained as outcomes rather than converted into discovery claims.
{links}"""
    zh=f"""# 任务九：公开结构锚定的合成构象异质性与机械变构计算

本项目已执行两端公开 X 射线结构配准、{s['conformers']} 个合成坐标构象、非线性嵌入、三维高斯散射代理、五态可逆马尔可夫模型、网格空腔/SASA 及 ANM/PRS。它**没有证明新的可成药口袋、癌症靶点机制、真实动力学速率或已验证变构路径**。模型蛋白是噬菌体 T4 溶菌酶 L99A，不是人源激酶或 GPCR。

## 已执行结果

{common}

在完整包含于ROI的已采样连通区域中，没有观察到“体积大于500 Å³且评分大于0.7”的联合条件；另有{s['threshold_indeterminate_frames']}帧触边，其完整体积和联合阈值分类为未知，**不能排除更大或可成药口袋的存在**。默认路径端点未跨越30 Å。开态仅有{m['observed_counts'][-1]}个样本，生成器预先设定的开态代价为{c['assumed_state_free_energies_kcal_mol'][-1]-c['assumed_state_free_energies_kcal_mol'][0]} kcal/mol。与估计的差别体现有限采样及显式伪计数正则化，不能宣称发现或验证目标能量区间。未为满足阈值调整随机种子、状态能量或评分参数。这里的ΔG是**状态占有率差异**，并非过渡态活化自由能或真实势垒。

## 9A：结构与连续几何构象

使用4W51和4W59的A链共有重原子，以残基编号和原子名逐一匹配。统一取alternate A构象；不混拼不同占有率原子。排除B构象、氢、水、离子和配体。4W59中的3GZ配体仅提供口袋中心，随后从几何计算删除。使用100–120残基以外的Cα做Kabsch配准，保证旋转矩阵行列式为+1。4W59端点属于配体条件下的结构，不能等同于已证实的无配体自发开态。

生成器先给出状态s，再取q=(s+u)/5，其中u服从(0.001,0.999)均匀分布。坐标为两端线性插值，叠加以110残基为中心逐渐衰减的随机环区位移；最大幅度为{c['synthetic_loop_displacement_A']} Å×sin(πq)。这使构象几何连续，但不能保证路径动力学可达。相邻Cα最短距离已显示插值压缩，且没有做力场能量最小化、共价立体化学或全重原子冲突审查；因此这些帧不是MD轨迹。

非线性嵌入使用配准Cα坐标距离构建RBF核，再做核PCA，输出两个潜变量。没有训练cryoDRGN，也没有真实单颗粒图像。各原子按原子序数加权投到{c['density_voxel_A']} Å体素，再用FWHM={c['density_psf_fwhm_A']} Å的高斯点扩散函数平滑。**FWHM是设置的PSF宽度，不是体素间距，也不是FSC验证的实验分辨率**。未模拟CTF、取向推断、探测噪声或定量库仑势。各状态平均和总体平均均依据实际生成帧数加权，散射单位为原子序数加权代理/Å³。图示只能说明坐标平均造成的平滑和局部衰减，不能声称真实环区密度消失。

## 9B：五态可逆模型与有限采样

在{c['temperature_K']} K下用配置中的五个假设自由能构造Boltzmann权重，以每方向{c['jump_probability']}/2的提议概率构造最近邻Metropolis转移，边界拒绝移出状态空间的步长。初态按平衡分布抽取。每帧**人为赋予{c['lag_ns']} ns**；未从结构或冷冻电镜推断时间。

从{s['conformers']-1}对相邻帧计数得到C，定义对称流量F=(C+Cᵀ)/2+αS，α={c['reversible_edge_pseudocount']}，S包含对角线和相邻态。归一化行得到T，行和归一化得到π，从代数上满足详细平衡。这是带正则的矩估计，不是假称的可逆最大似然估计。PMF为−RT ln(π_i/π_0)。将状态4吸收，解(I−T_nonabsorbing)m=τ1得到首次到达时间；1/m_0为假设时钟下的条件开启速率。

实际完成{c['bootstrap_replicates']}次参数自助法：用拟合T重抽样等长轨迹，并使用同一正则重新拟合。95%百分位区间仅包含该合成模型的采样变化，未包括结构、力场、时间标尺等不确定性。缺态重抽样次数单独保留，不静默删除。2τ的Chapman–Kolmogorov Frobenius差为{m['ck_lag2_frobenius_disagreement']:.6f}，只是有限样本诊断，不能据此宣布真实蛋白满足马尔可夫性。

## 9C：从几何计算空腔、表面积和评分

围绕配体中心建立半宽{c['roi_halfwidth_A']} Å的局部立方网格。以C/N/O/S=1.70/1.55/1.52/1.80 Å为假设vdW半径，检测{c['pocket_probe_A']} Å探针是否与蛋白相交。候选体素要求六个轴向至少五向在ROI内遇到蛋白；选择距离配体中心最近的候选连通区域。体积为体素数量×间距³。该算法是局部球填充/遮挡方法，不是全蛋白无偏口袋搜索、alpha-shape或商用可成药评分。

距空腔体素小于自身vdW半径+2.5 Å的原子定义衬里。使用{c['sasa_sphere_points']}点Fibonacci球与{c['sasa_probe_A']} Å探针计算衬里可接近面积，包含内部封闭表面，不证明与体相溶剂连通。芳香性定义为属于Phe/Tyr/Trp/His残基的衬里原子比例，非精确芳香表面积；极性为N/O/S比例。遮挡比例是六方向平均遮挡，疏水围合深度是连通区域最大距离变换半径×衬里碳比例，均明确为几何代理。

原提示词公式并不在[0,1]内。这里保留raw=0.49 ln[V/(1 Å³)]+0.78 ln(芳香比例)−0.22极性−1.2；对零芳香比例使用10⁻⁶下限，再额外计算sigmoid(raw)。后者仅为**未校准有界描述符**，不等于结合概率或成药概率。空腔为空时raw记录null、评分记0。有限帧数不足以给出精确连续阈值，未观察到达到的门槛不能擅自声明达到。

触及ROI边界的帧不能作为完整闭合口袋体积；完整体积置null、阈值分类标为indeterminate，仅保留ROI内计算值。已在两端完成0.5/0.75/1.0 Å网格比较，另在闭/中间/开三个代表结构上执行9/12/15 Å半宽ROI敏感性。ROI扩大可能同时改变方向遮挡与连通性，因此扩大后的数值并非自动代表同一个物理空腔。逐帧截断标志仍保留；体素结果接近不代表全空腔体积收敛。

## 9D：ANM伪逆与扰动响应

在闭端点Cα上，用{c['anm_cutoff_A']} Å接触截断和单位弹簧构造3N×3N Hessian。接触对在对角块加uuᵀ、交叉块减uuᵀ。按相对本征值阈值剔除六个平移/旋转零模，使用其余全部内部模构造H⁺，没有直接对奇异H求逆。弹簧常数和外力幅度未经物理标定，位移仅有任意响应单位。

对源残基j施加单位各向同性力，精确平均平方响应为M_ij=||H⁺_ij||²_F/3。每残基实际抽取{c['prs_force_counts']}个随机方向，通过经验力协方差计算，代数上等价于逐方向计算位移后求平方均值。闭式结果与随机方向结果比较，最后的相对误差见上表。每列按源自身响应归一化；图中列是施力残基，行是响应残基。

耦合定义为M_ij/sqrt(M_ii M_jj)。只保留接触边，代价为−ln(耦合)+0.01距离/Å，用Dijkstra从全部被检测衬里残基到参考活性位点{c['active_site_residue']}选择最低代价路径。该定义可能偏向近邻衬里端点；结果不能冒充跨越30 Å的信号传递，也不表示真实能量流、因果调控或突变验证。必须在报告保留实际端距和路径长度。

## 业界使用场景与验收边界

本项目适合作为结构驱动先导发现前的计算流程试验：分别审查状态支持度、ROI截断、评分校准和机械路径可靠性。要推进真实靶点，应取得原始冷冻电镜颗粒或经验证MD、人口与动力学测量、侧链和水环境优化、无偏口袋检测、片段结合基准，以及多结构/接触参数敏感性和突变实验。只有这些证据逐步补齐后，才能把结果用于片段筛选或药化优先级。本次未满足的门槛已作为负结果保留。
{links}"""
    (out/'CRYOEM_CRYPTIC_POCKET_REPORT_EN.md').write_text(en,encoding='utf-8')
    (out/'CRYOEM_CRYPTIC_POCKET_REPORT_ZH.md').write_text(zh,encoding='utf-8')

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--out',type=Path); parser.add_argument('--config',type=Path); parser.add_argument('--write-example',type=Path); parser.add_argument('--self-test',action='store_true'); a=parser.parse_args()
    if a.write_example:
        if a.write_example.exists(): parser.error('--write-example destination already exists')
        dump(a.write_example,DEFAULT); return 0
    if a.self_test:
        suite=unittest.defaultTestLoader.discover(str(HERE.parents[1]/'tests'),pattern='test_task9.py')
        if not suite.countTestCases(): raise RuntimeError('No Task9 tests discovered')
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    config=dict(DEFAULT)
    if a.config: config.update(json.loads(a.config.read_text(encoding='utf-8')))
    if a.out is None: parser.error('--out NEW_DIRECTORY is required; archived artifacts are never overwritten by the CLI')
    if a.out.exists(): parser.error('--out destination already exists; use a new directory')
    run(a.out.resolve(),config); return 0

if __name__=='__main__': raise SystemExit(main())
