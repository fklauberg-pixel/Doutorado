import numpy as np, sys, rasterio, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rasterio.enums import Resampling
from matplotlib.patches import Polygon
tag,f=sys.argv[1],sys.argv[2]; fu=float(sys.argv[3]) if len(sys.argv)>3 else 0.75
cx,cy,th,pu,pv,nu,nv=np.load(f"{tag}_grid.npy"); nu,nv=int(nu),int(nv)
u=np.array([np.cos(th),np.sin(th)]); v=np.array([-np.sin(th),np.cos(th)])
with rasterio.open(f) as s:
    sc=0.03/s.res[0]; H,W=int(s.height/sc),int(s.width/sc)
    a=s.read(out_shape=(4,H,W),resampling=Resampling.average); b=s.bounds
img=np.moveaxis(a[:3],0,-1).copy(); img[a[3]==0]=255
fig,ax=plt.subplots(figsize=(12,13)); ax.imshow(img,extent=[b.left,b.right,b.bottom,b.top])
k=0
for i in range(nu):
    for j in range(nv):
        c=np.array([cx,cy])+(i-(nu-1)/2)*pu*u+(j-(nv-1)/2)*pv*v
        hu,hv=fu*pu/2,fu*pv/2
        P=[c+su*hu*u+sv*hv*v for su,sv in [(-1,-1),(1,-1),(1,1),(-1,1)]]
        ax.add_patch(Polygon(P,fill=False,ec="red",lw=1.2)); ax.text(*c,f"{i},{j}",color="yellow",fontsize=8,ha="center")
ax.set_xlim(cx-15,cx+15); ax.set_ylim(cy-15,cy+15)
fig.savefig(f"{tag}_overlay.png",dpi=70,bbox_inches="tight")
