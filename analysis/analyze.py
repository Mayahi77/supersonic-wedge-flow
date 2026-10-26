"""Extract cell data, shock comparisons and synchronized plots from ASCII OpenFOAM.

Uses the actual z=0 boundary face topology and its owner mapping; no assumed
Gmsh/OpenFOAM cell ordering. Run artifacts are kept outside the source case.
"""
import argparse
import gzip
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.path import Path as PolygonPath
from PIL import Image


def foam_list(path):
    text = path.read_text()
    match = re.search(r'\n\s*(\d+)\s*\n\(', text)
    if not match:
        raise ValueError(str(path))
    return int(match[1]), text[match.end():].split('\n)', 1)[0]


def field(path, n, vector=False):
    text = path.read_text()
    match = re.search(r'internalField\s+uniform\s+([^;]+);', text)
    if match:
        value = np.fromstring(match[1].strip('() '), sep=' ')
        return np.tile(value, (n, 1)) if vector else np.full(n, value[0])
    match = re.search(r'internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\(', text)
    assert int(match[1]) == n
    values = np.fromstring(text[match.end():].split('\n)', 1)[0].replace('(', '').replace(')', ''), sep=' ')
    return values.reshape(n, 3) if vector else values


def mesh(case):
    base = case / 'constant/polyMesh'
    _, s = foam_list(base / 'points')
    points = np.fromstring(s.replace('(', '').replace(')', ''), sep=' ').reshape(-1, 3)
    _, s = foam_list(base / 'faces')
    faces = [np.fromstring(m, sep=' ', dtype=int) for m in re.findall(r'\d+\(([^)]*)\)', s)]
    _, s = foam_list(base / 'owner')
    owners = np.fromstring(s, sep=' ', dtype=int)
    n = owners.max() + 1
    polygons = [None] * n
    for face, owner in zip(faces, owners):
        if np.all(np.abs(points[face, 2]) < 1e-10):
            polygons[owner] = points[face, :2]
    assert all(p is not None and len(p) == 4 for p in polygons)
    polygons = np.array(polygons)
    centers = polygons.mean(axis=1)
    areas = np.abs(np.sum(polygons[:, :, 0] * np.roll(polygons[:, :, 1], -1, axis=1)
                          - polygons[:, :, 1] * np.roll(polygons[:, :, 0], -1, axis=1), axis=1)) / 2
    return polygons, centers, areas


def theory():
    gamma, M, theta = 1.4, 2., np.deg2rad(10)
    def f(b):
        return 2 / np.tan(b) * (M*M*np.sin(b)**2 - 1) / (M*M*(gamma+np.cos(2*b))+2) - np.tan(theta)
    a, b = np.arcsin(1/M), np.deg2rad(60)
    for _ in range(70):
        c = (a+b)/2
        if f(c) > 0: b = c
        else: a = c
    beta = (a+b)/2
    mn = M*np.sin(beta)
    pr = 1+2*gamma/(gamma+1)*(mn*mn-1)
    rr = (gamma+1)*mn*mn/((gamma-1)*mn*mn+2)
    m2n = np.sqrt((1+(gamma-1)/2*mn*mn)/(gamma*mn*mn-(gamma-1)/2))
    return dict(beta_deg=float(np.rad2deg(beta)), pressure_ratio=float(pr),
                p_Pa=float(1e5*pr), T_K=float(300*pr/rr), M=float(m2n/np.sin(beta-theta)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--animate', action='store_true')
    args = parser.parse_args()
    case, out = args.case.resolve(), args.output.resolve()
    for name in ['results', 'media', 'media/frames']: (out/name).mkdir(parents=True, exist_ok=True)
    polygons, centers, areas = mesh(case)
    n = len(centers)
    times = sorted((float(p.name), p) for p in case.iterdir() if p.is_dir() and re.fullmatch(r'[\d.eE+-]+', p.name) and (p/'p').exists())
    # Recover the effective gas constant from saved thermodynamic state, avoiding
    # small differences among OpenFOAM releases' physical constants.
    last = times[-1][1]
    R_cells = field(last/'p',n)/(field(last/'rho',n)*field(last/'T',n))
    R = float(np.median(R_cells))
    assert np.max(abs(R_cells/R-1)) < 1e-7
    cp = 1004.703
    gamma = cp/(cp-R)
    th = theory()
    beta = np.deg2rad(th['beta_deg'])
    # Isolated near-corner post-shock region, away from the smeared shock/wall/outlet.
    downstream = (centers[:,0]>.7)&(centers[:,0]<.95)&(centers[:,1]>.09)&(centers[:,1]<.14)
    upstream = (centers[:,0]>.15)&(centers[:,0]<.4)&(centers[:,1]>.3)&(centers[:,1]<.7)
    assert downstream.sum()>5 and upstream.sum()>5
    # Cell intersections of requested line x=.8, z=.005; cell values, not interpolated.
    intersections=[]
    for i, poly in enumerate(polygons):
        ys=[]
        for a,b in zip(poly,np.roll(poly,-1,axis=0)):
            if min(a[0],b[0]) <= .8 < max(a[0],b[0]):
                ys.append(a[1]+(.8-a[0])/(b[0]-a[0])*(b[1]-a[1]))
        if len(ys)==2:
            lo,hi=sorted(ys)
            lo,hi=max(lo,.15),min(hi,.8)
            if hi>lo: intersections.append((i,lo,hi))
    intersections.sort(key=lambda v:v[1])
    idx=np.array([r[0] for r in intersections]); y=np.array([(r[1]+r[2])/2 for r in intersections])
    cache={}; history=[]; profiles=[]
    for t,directory in times:
        p=field(directory/'p',n); T=field(directory/'T',n); U=field(directory/'U',n,True)
        M=np.linalg.norm(U,axis=1)/np.sqrt(gamma*R*T)
        assert all(np.isfinite(v).all() for v in [p,T,U,M]) and T.min()>0 and p.min()>0
        mean=lambda v,mask:float(np.average(v[mask],weights=areas[mask]))
        history.append([t,mean(p,upstream),mean(p,downstream),mean(T,downstream),mean(M,downstream)])
        profiles.append(p[idx]); cache[t]=(p,T,M)
    history=np.array(history); profiles=np.array(profiles)
    np.savetxt(out/'results/history.csv',history,delimiter=',',header='time_s,upstream_p_Pa,downstream_p_Pa,downstream_T_K,downstream_M',comments='')
    np.savetxt(out/'results/pressure-profiles.csv',np.column_stack([history[:,0],profiles]),delimiter=',',header='time_s,'+','.join('y_%.8g_m'%v for v in y),comments='')
    final_t=max(cache); p,T,M=cache[final_t]
    np.savez_compressed(out/'results/selected-fields.npz',polygons=polygons,centers=centers,areas=areas,
        p_001=cache[.01][0],T_001=cache[.01][1],M_001=cache[.01][2],p_002=p,T_002=T,M_002=M)
    metrics={'theory_nominal_M2':th,'gas':{'R_J_kgK':R,'gamma':gamma,'inlet_M':694.45/np.sqrt(gamma*R*300)},
             'downstream_roi':{'x_m':[.7,.95],'y_m':[.09,.14],'cells':int(downstream.sum())},
             'upstream_roi':{'x_m':[.15,.4],'y_m':[.3,.7],'cells':int(upstream.sum())},
             'final_time_s':final_t,'time_directories':len(times),'final_roi':dict(zip(['time_s','upstream_p_Pa','downstream_p_Pa','downstream_T_K','downstream_M'],history[-1].tolist()))}
    metrics['saved_snapshot_change_001_to_002']={}
    for label,k in [('p',0),('T',1),('M',2)]:
        old,new=cache[.01][k],cache[.02][k]
        metrics['saved_snapshot_change_001_to_002'][label]={'area_weighted_relative_L1':float(np.sum(areas*abs(new-old))/np.sum(areas*abs(old))), 'max_absolute':float(abs(new-old).max())}
    verification=out/'.work/verification'
    verification_history=[]
    if verification.exists():
        for t in np.arange(1,21)*.001:
            directory=verification/f'{t:.3f}'.rstrip('0')
            pv=field(directory/'p',n);tv=field(directory/'T',n);uv=field(directory/'U',n,True)
            mv=np.linalg.norm(uv,axis=1)/np.sqrt(gamma*R*tv)
            verification_history.append([t,mean(pv,upstream),mean(pv,downstream),mean(tv,downstream),mean(mv,downstream)])
        verification_history=np.array(verification_history)
        np.savetxt(out/'results/verification-history.csv',verification_history,delimiter=',',header='time_s,upstream_p_Pa,downstream_p_Pa,downstream_T_K,downstream_M',comments='')
        metrics['verification_stationarity']={}
        for t0,t1 in [('0.01','0.02'),('0.019','0.02')]:
            change={}
            for label in ['p','T']:
                old=field(verification/t0/label,n);new=field(verification/t1/label,n)
                change[label]={'area_weighted_relative_L1':float(np.sum(areas*abs(new-old))/np.sum(areas*abs(old))),'max_absolute':float(abs(new-old).max())}
            metrics['verification_stationarity'][t0+'_to_'+t1]=change
    # Shock width between 10%/90% pressure-rise thresholds, using descending profile.
    pf=profiles[-1]; thresholds=1e5+(th['p_Pa']-1e5)*np.array([.9,.5,.1])
    crossings=[]
    for threshold in thresholds:
        hits=np.where((pf[:-1]>=threshold)&(pf[1:]<threshold))[0]
        crossings.append(float(np.interp(threshold,[pf[hits[0]+1],pf[hits[0]]],[y[hits[0]+1],y[hits[0]]])) if len(hits) else None)
    metrics['profile']={'x_m':.8,'z_m':.005,'y_range_m':[.15,.8], 'sampling':'piecewise-constant cell intersections; thresholds linearly interpolated between segment midpoints',
                        'theoretical_shock_y_m':float(.3*np.tan(beta)), 'y_at_90_50_10_percent':crossings,
                        'width_10_90_m':crossings[2]-crossings[0] if None not in crossings else None}
    log=(out/'logs/original-startup-rhoCentralFoam.log')
    if log.exists() or log.with_suffix('.log.gz').exists():
        log_text=log.read_text() if log.exists() else gzip.open(log.with_suffix('.log.gz'),'rt').read()
        co=re.findall(r'Mean and max Courant Numbers = ([\d.eE+-]+) ([\d.eE+-]+)',log_text)
        metrics['observed_max_Co']=max(float(v[1]) for v in co)
    (out/'results/metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for label,values,unit,clim in [('pressure',p/1000,'p [kPa]',(100,180)),('temperature',T,'T [K]',(300,360)),('mach',M,'Mach number',(1.5,2.05))]:
        fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
        coll=PolyCollection(polygons,array=values,cmap='viridis',clim=clim,edgecolors='none',antialiased=False);ax.add_collection(coll)
        ax.plot([.5,1.5],[0,np.tan(beta)],'w--',lw=1,label='Nominal oblique-shock theory')
        ax.plot([.8,.8],[.15,.8],'r-',lw=1,label='Pressure sampling line')
        ax.set(xlim=(0,1.5),ylim=(0,1),xlabel='x [m]',ylabel='y [m]',title=f'OpenFOAM {label}, t = {final_t:.3f} s');ax.set_aspect('equal')
        ax.legend(loc='upper left',fontsize=8);fig.colorbar(coll,ax=ax,label=unit,shrink=.8)
        fig.savefig(out/f'media/{label}.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
    for t in [.001,.002,.005,.01,.02]:
        i=np.argmin(abs(history[:,0]-t));ax.plot(y,profiles[i]/1000,'o-',ms=2,label=f'{history[i,0]:g} s')
    ax.axhline(th['p_Pa']/1000,color='.3',ls='--',label='Theory downstream')
    ax.axhline(100,color='.5',ls=':');ax.axvline(.3*np.tan(beta),color='.3',ls=':',label='Theory shock crossing')
    ax.set(xlabel='y [m] at x = 0.8 m, z = 0.005 m',ylabel='p [kPa]',title='Shock formation and numerical smearing');ax.legend(fontsize=8,ncol=2)
    fig.savefig(out/'media/pressure-profile.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(8,7),sharex=True,constrained_layout=True)
    for ax,col,ref,label in zip(axes,[2,3,4],[th['p_Pa'],th['T_K'],th['M']],['p [Pa]','T [K]','Mach']):
        ax.plot(history[:-1,0],history[:-1,col],label='Saved startup ROI mean');ax.plot(history[-1,0],history[-1,col],'o',label='Saved 0.02 s snapshot')
        if len(verification_history): ax.plot(verification_history[:,0],verification_history[:,col],':',label='Independent verification run')
        ax.axhline(ref,color='.3',ls='--',label='Nominal theory');ax.set_ylabel(label);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);axes[0].set_title('Downstream region: 0.7 < x < 0.95 m, 0.09 < y < 0.14 m');axes[-1].set_xlabel('Physical time [s]')
    fig.savefig(out/'media/stationarity.png',dpi=170);plt.close(fig)
    if args.animate:
        ts=[t for t,_ in times if t<=.01]
        fig=plt.figure(figsize=(12,7))
        ax=fig.add_axes([.12,.48,.76,.46]);prof=fig.add_axes([.08,.09,.39,.28]);hist=fig.add_axes([.58,.09,.38,.28])
        coll=PolyCollection(polygons,array=cache[0.][0]/1000,cmap='viridis',clim=(100,180),edgecolors='none',antialiased=False);ax.add_collection(coll)
        ax.set(xlim=(0,1.5),ylim=(0,1),xlabel='x [m]',ylabel='y [m]');ax.set_aspect('equal');ax.plot([.8,.8],[.15,.8],color='red',lw=1)
        fig.colorbar(coll,cax=fig.add_axes([.88,.50,.018,.4]),label='p [kPa]')
        line,=prof.plot(y,profiles[0]/1000,color='tab:blue');prof.axhline(th['p_Pa']/1000,color='.4',ls='--');prof.axvline(.3*np.tan(beta),color='.5',ls=':')
        prof.set(xlim=(.15,.8),ylim=(95,185),xlabel='y [m] at x = 0.8 m',ylabel='p [kPa]',title='Synchronized pressure profile')
        hist.plot(history[:len(ts),0],history[:len(ts),2]/1000,color='.6');hist.axhline(th['p_Pa']/1000,color='.4',ls='--');cursor=hist.axvline(0,color='tab:red')
        dot,=hist.plot([0],[100],'ro');hist.set(xlim=(0,.01),ylim=(95,185),xlabel='Physical time [s]',ylabel='Downstream ROI p [kPa]',title='ROI history; dashed line = nominal theory')
        def update(i):
            t=ts[i];coll.set_array(cache[t][0]/1000);line.set_ydata(profiles[i]/1000);cursor.set_xdata([t,t]);dot.set_data([t],[history[i,2]/1000]);ax.set_title(f'OpenFOAM startup | t = {t:.5f} s | cell data (Matplotlib)')
            return coll,line,cursor,dot
        for i in range(len(ts)):
            update(i);fig.savefig(out/f'media/frames/frame-{i:04d}.png',dpi=100)
        plt.close(fig)
        # One global palette prevents frame-dependent palette flicker and keeps
        # the GIF small. Fixed axes prevent text/camera movement during playback.
        palette_source=Image.new('RGB',(600,350*3),'white')
        for j,i in enumerate([0,20,len(ts)-1]):
            palette_source.paste(Image.open(out/f'media/frames/frame-{i:04d}.png').convert('RGB').resize((600,350)),(0,350*j))
        palette=palette_source.quantize(colors=96)
        gif_frames=[]
        for i in range(0,len(ts),4):
            frame=Image.open(out/f'media/frames/frame-{i:04d}.png').convert('RGB').resize((600,350),Image.LANCZOS)
            gif_frames.append(frame.quantize(palette=palette,dither=Image.NONE))
        gif_frames[0].save(out/'media/startup-plots-preview.gif',save_all=True,append_images=gif_frames[1:],duration=200,loop=0,optimize=False,disposal=1)
    print(json.dumps(metrics,indent=2))


if __name__=='__main__': main()
