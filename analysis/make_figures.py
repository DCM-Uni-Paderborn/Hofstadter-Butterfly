"""Create publication figures without altering the supplied DFT image data."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
 'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.6})

def xyz(path):
    rows=path.read_text().splitlines()[2:]
    return np.array([r.split()[0] for r in rows]), np.array([[float(v) for v in r.split()[1:4]] for r in rows])

records=[]
files=sorted((ROOT/'data/structures_d2p15').glob('*/atoms.xyz'),key=lambda p:float(p.parent.name.split('angle_')[1]))
for f in files:
    el,p=xyz(f)
    records.append({'theta_deg':float(f.parent.name.split('angle_')[1]),'total':len(p),
     'bottom_C':int(np.sum(p[:,2]==0)),'top_C':int(np.sum((p[:,2]>0)&(el=='C'))),
     'top_H':int(np.sum(el=='H'))})
with (ROOT/'data/geometry_counts.csv').open('w') as out:
    writer=csv.DictWriter(out,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)

fig,ax=plt.subplots(1,2,figsize=(3.4,2.0),gridspec_kw={'width_ratios':[1.65,1]},layout='constrained')
el,p=xyz(files[25]);low=p[:,2]==0;top=(p[:,2]>0)&(el=='C');hyd=el=='H'
ax[0].scatter(p[low,0]/10,p[low,1]/10,s=.45,c='#777777',label='lower C',rasterized=True)
ax[0].scatter(p[top,0]/10,p[top,1]/10,s=.5,c='#0072B2',label='upper C',rasterized=True)
ax[0].scatter(p[hyd,0]/10,p[hyd,1]/10,s=2.5,c='#D55E00',label='edge H',rasterized=True)
A=np.array([61.4878,0])/10;B=np.array([30.7439,53.25])/10
box=np.array([[0,0],A,A+B,B,[0,0]])
ax[0].plot(box[:,0],box[:,1],lw=.6,c='black');ax[0].set_aspect('equal')
ax[0].set_xlabel('$x$ (nm)');ax[0].set_ylabel('$y$ (nm)');ax[0].set_title(r'(a) DFT, $45^\circ$',loc='left')
R=4.5;v1=np.array([np.sqrt(3),0]);v2=np.array([np.sqrt(3)/2,1.5]);basis=np.array([[0,0],[0,1]])
P=np.array([i*v1+j*v2+b for i in range(-10,11) for j in range(-10,11) for b in basis if np.linalg.norm(i*v1+j*v2+b)<=R])
rot=np.array([[np.cos(np.pi/4),-np.sin(np.pi/4)],[np.sin(np.pi/4),np.cos(np.pi/4)]]);Q=P@rot.T
ax[1].scatter(P[:,0],P[:,1],s=3,c='#777777');ax[1].scatter(Q[:,0],Q[:,1],s=3,c='#0072B2')
ax[1].add_patch(plt.Circle((0,0),R,fill=False,color='black',lw=.6));ax[1].set_aspect('equal');ax[1].axis('off')
ax[1].set_title('(b) Model',loc='left');ax[1].text(.5,-.04,'Coincident disks\n(schematic)',ha='center',va='top',transform=ax[1].transAxes)
fig.savefig(FIG/'geometry.pdf');plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(6.8,2.2),layout='constrained')
th=[r['theta_deg'] for r in records]
axes[0].plot(th,[r['top_C'] for r in records],'.-',ms=3,lw=.7,color='#0072B2');axes[0].set_ylabel('Upper-layer C atoms')
axes[1].plot(th,[r['top_H'] for r in records],'.-',ms=3,lw=.7,color='#D55E00');axes[1].set_ylabel('Terminating H atoms')
for i,a in enumerate(axes):a.set_xlabel(r'$\theta$ (deg)');a.set_title(f'({chr(97+i)})',loc='left');a.grid(alpha=.15)
fig.savefig(FIG/'geometry_counts.pdf');plt.close(fig)

from assemble_dft_panels import assemble
assemble()

def tb_dos(file):
    data=np.load(file);positive=data['theta'];es=data['eigenvalues']
    # Layer exchange plus rigid rotation proves equality at +theta and -theta.
    angles=np.r_[-positive[:0:-1],positive]*180/np.pi
    es=np.concatenate([es[:0:-1],es],axis=0)
    lo,hi=es.min(),es.max();grid=np.linspace(lo,hi,2001);sigma=(hi-lo)/1000
    rho=np.empty((len(grid),len(angles)))
    for i,ev in enumerate(es):
        rho[:,i]=np.exp(-.5*((grid[:,None]-ev[None,:])/sigma)**2).mean(axis=1)/(np.sqrt(2*np.pi)*sigma)
    return data,angles,grid,rho,sigma

tbfiles=[ROOT/f'data/tb_reproduction/R25_d{d}_xi0.03.npz' for d in ['0.99','0.98','0.97']]
if all(p.exists() for p in tbfiles):
    fig,axes=plt.subplots(1,3,figsize=(6.8,2.45),layout='constrained')
    report=[]
    for i,(ax,file) in enumerate(zip(axes,tbfiles)):
        data,angles,grid,rho,sigma=tb_dos(file)
        im=ax.pcolormesh(angles,grid,rho,shading='auto',cmap='viridis',norm=Normalize(0,.5),rasterized=True)
        ax.set_facecolor(plt.get_cmap('viridis')(0))
        ax.set_xlabel(r'$\theta$ (deg)');ax.set_xticks([-60,-30,0,30,60]);ax.set_ylim(-6,6)
        if i==0:ax.set_ylabel(r'$E/t_{\rm nn}$')
        ax.set_title(f'({chr(97+i)}) $d/a={float(data["d"]):.2f}$',loc='left')
        report.append({'d':float(data['d']),'dimension':len(data['eigenvalues'][0]),'sigma':float(sigma),'min':float(grid.min()),'max':float(grid.max()),'max_dos':float(rho.max())})
    fig.colorbar(im,ax=axes,label=r'$t_{\rm nn}\rho(E)$',shrink=.86,pad=.015,extend='max')
    fig.savefig(FIG/'tb_main.pdf',dpi=400);plt.close(fig)
    (ROOT/'data/tb_rendering.json').write_text(json.dumps(report,indent=2))

rangefiles=[ROOT/f'data/tb_reproduction/R25_d0.97_xi{xi}.npz' for xi in ['0.03','0.1','0.3']]
if all(p.exists() for p in rangefiles):
    fig,axes=plt.subplots(1,3,figsize=(6.8,2.8),layout='constrained')
    for i,(ax,file) in enumerate(zip(axes,rangefiles)):
        data,angles,grid,rho,sigma=tb_dos(file)
        im=ax.pcolormesh(angles,grid,rho,shading='auto',cmap='viridis',norm=Normalize(0,.5),rasterized=True)
        ax.set_facecolor(plt.get_cmap('viridis')(0))
        ax.set_xlabel(r'$\theta$ (deg)');ax.set_xticks([-60,-30,0,30,60]);ax.set_ylim(-6,6)
        if i==0:ax.set_ylabel(r'$E/t_{\rm nn}$')
        ax.set_title(rf'({chr(97+i)}) $\xi/a={float(data["xi"]):g}$',loc='left')
    fig.colorbar(im,ax=axes,label=r'$t_{\rm nn}\rho(E)$',shrink=.86,pad=.015,extend='max')
    fig.savefig(FIG/'tb_range.pdf',dpi=400);plt.close(fig)
print('Figures built; DFT panels are unchanged source images.')
