import numpy as np, sys, rasterio, pandas as pd
from rasterio.features import geometry_mask
from shapely.geometry import Polygon, mapping
from scipy import ndimage as ndi
sys.path.insert(0,"/home/claude/Doutorado/paper3-lai-drone/analise")
from indices import compute_indices
tag,ortho,dsm=sys.argv[1:4]; FR=0.75   # retângulo útil = 75% do passo da grade
cx,cy,th,pu,pv,nu,nv=np.load(f"{tag}_grid.npy"); nu,nv=int(nu),int(nv)
u=np.array([np.cos(th),np.sin(th)]); v=np.array([-np.sin(th),np.cos(th)])
def rect(c,fu,fv):
    hu,hv=fu*pu/2,fv*pv/2
    return Polygon([c+su*hu*u+sv*hv*v for su,sv in [(-1,-1),(1,-1),(1,1),(-1,1)]])
so,sd=rasterio.open(ortho),rasterio.open(dsm)
assert so.transform==sd.transform
rows=[]; allexg=[]
cache={}
for i in range(nu):
    for j in range(nv):
        c=np.array([cx,cy])+(i-(nu-1)/2)*pu*u+(j-(nv-1)/2)*pv*v
        P=rect(c,FR,FR); Pz=rect(c,1.0,1.0)   # zona com bordadura (inclui metade das ruas) p/ solo
        win=rasterio.windows.from_bounds(*Pz.bounds,transform=so.transform).round_offsets().round_lengths()
        tr=so.window_transform(win)
        a=so.read(window=win).astype(float); z=sd.read(1,window=win).astype(float)
        shp=a.shape[1:]
        m_in=~geometry_mask([mapping(P)],shp,tr); m_z=~geometry_mask([mapping(Pz)],shp,tr)
        ok=(a[3]>0)&(z>-9000)
        R,G,B=a[0],a[1],a[2]
        cache[(i,j)]=(R,G,B,z,m_in&ok,m_z&ok)
        S_=R+G+B; exgn=np.where(S_>0,(2*G-R-B)/np.maximum(S_,1),np.nan)
        allexg.append(exgn[m_in&ok])
# limiar de Otsu do ExG normalizado, global
vals=np.concatenate(allexg); vals=vals[np.isfinite(vals)]
h,e=np.histogram(vals,256); cc=(e[:-1]+e[1:])/2; w0=np.cumsum(h); w1=w0[-1]-w0
m0=np.cumsum(h*cc)/np.maximum(w0,1); m1=(np.sum(h*cc)-np.cumsum(h*cc))/np.maximum(w1,1); thr=cc[np.argmax(w0*w1*(m0-m1)**2)]
for (i,j),(R,G,B,z,mi,mz) in cache.items():
    r_,g_,b_=R[mi],G[mi],B[mi]
    idx=compute_indices(r_,g_,b_); veg=idx["ExG"]>thr
    row=dict(gi=i,gj=j,n_pix=mi.sum(),CC=veg.mean(),ExG_DN=np.mean(2*g_-r_-b_))
    for k,vv in idx.items(): row[k]=np.nanmean(vv)
    mR,mG,mB=r_.mean(),g_.mean(),b_.mean()
    for k,vv in compute_indices(np.array([mR]),np.array([mG]),np.array([mB])).items(): row[k+"_bm"]=float(vv[0])
    zz=z[mi]; zground=np.percentile(z[mz],2)
    # solo local: pixels de solo (não-vegetação) na zona ampliada
    Sz=R+G+B; exz=np.where(Sz>0,(2*G-R-B)/np.maximum(Sz,1),np.nan)
    soil=mz&(exz<thr)
    zsoil=np.median(z[soil]) if soil.sum()>100 else np.nan
    row.update(z_solo_p2=zground,z_solo_med=zsoil)
    for q in (50,90,95,99): row[f"H_p{q}"]=np.percentile(zz,q)-zsoil
    vz=zz[veg]
    row["H_veg_media"]=(vz.mean()-zsoil) if vz.size else np.nan
    hh=np.clip(zz-zsoil,0,None); row["Vol_m3_m2"]=hh.mean()     # volume do dossel por m² (m)
    row["H_dp"]=zz.std()
    rows.append(row)
df=pd.DataFrame(rows); df.to_csv(f"{tag}_parcelas.csv",index=False)
print(tag,"Otsu",thr); print(df[["gi","gj","n_pix","CC","ExG_DN","ExG","H_p90","H_p99","H_veg_media","Vol_m3_m2"]].round(3).to_string(index=False))
