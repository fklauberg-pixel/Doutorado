import rasterio, numpy as np, sys
from rasterio.enums import Resampling
from scipy import ndimage as ndi
f, tag = sys.argv[1], sys.argv[2]
res = 0.05
with rasterio.open(f) as s:
    sc = res / s.res[0]
    H, W = int(s.height / sc), int(s.width / sc)
    a = s.read(out_shape=(4, H, W), resampling=Resampling.average).astype(float)
    tr = s.transform * s.transform.scale(s.width / W, s.height / H)
R, G, B, A = a
S = R + G + B + 1e-9
exg = 2 * G / S - R / S - B / S
veg = (exg > 0.08) & (A > 0)
# vegetação densa em blocos: média móvel 1 m
dens = ndi.uniform_filter(veg.astype(float), size=int(1.0 / res))
np.save(f"{tag}_dens.npy", dens)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 8)); ax.imshow(dens, cmap="Greens"); fig.savefig(f"{tag}_dens.png", dpi=90); plt.close(fig)
print(tag, H, W, tr)
