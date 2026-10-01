import numpy as np, sys, itertools
from scipy import ndimage as ndi
from scipy.optimize import minimize
tag,x0,y0=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]); cx0,cy0=float(sys.argv[4]),float(sys.argv[5]); res=0.05
d=np.load(f"{tag}_dens.npy")
sm=ndi.uniform_filter(d,size=int(2.4/res))
def samp(x,y):
    c=(x-x0)/res-0.5; r=(y0-y)/res-0.5
    return ndi.map_coordinates(sm,[r,c],order=1,mode="constant",cval=0)
def centers(p,nu,nv):
    cx,cy,th,pu,pv=p
    u=np.array([np.cos(th),np.sin(th)]); v=np.array([-np.sin(th),np.cos(th)])
    I,J=np.meshgrid(np.arange(nu)-(nu-1)/2,np.arange(nv)-(nv-1)/2,indexing="ij")
    return cx+I*pu*u[0]+J*pv*v[0], cy+I*pu*u[1]+J*pv*v[1], u, v
def score(p,nu,nv):
    X,Y,u,v=centers(p,nu,nv)
    inside=samp(X.ravel(),Y.ravel())
    al=[]
    for s in (-0.5,0.5):
        al.append(samp((X+s*p[3]*u[0]).ravel(),(Y+s*p[3]*u[1]).ravel()))
        al.append(samp((X+s*p[4]*v[0]).ravel(),(Y+s*p[4]*v[1]).ravel()))
    return -(inside.mean()-np.concatenate(al).mean())
best=None
for nu,nv in [(6,4),(4,6)]:
    for th0 in np.deg2rad(np.arange(0,180,5)):
        for pu0,pv0 in [(4.0,5.5),(5.5,4.0),(4.5,4.5),(5,5)]:
            for dx,dy in itertools.product([-1,0,1],[-1,0,1]):
                r=minimize(score,[cx0+dx,cy0+dy,th0,pu0,pv0],args=(nu,nv),method="Nelder-Mead",options=dict(xatol=1e-3,fatol=1e-5,maxiter=3000))
                if best is None or r.fun<best[0]: best=(r.fun,r.x,nu,nv)
f,p,nu,nv=best
print(tag,"score",-f,"params cx cy th(deg) pu pv",p[0],p[1],np.rad2deg(p[2])%360,p[3],p[4],"nu nv",nu,nv)
np.save(f"{tag}_grid.npy",np.r_[p,nu,nv])
X,Y,u,v=centers(p,nu,nv)
inside=samp(X.ravel(),Y.ravel()); print("densidade nas 24 parcelas min/mediana",inside.min().round(2),np.median(inside).round(2))
