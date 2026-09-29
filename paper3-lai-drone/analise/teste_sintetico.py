"""Gera um ortomosaico e parcelas SINTÉTICOS e roda o pipeline inteiro.

Serve só para verificar que o código funciona antes de chegarem os dados reais.
Nenhum número produzido aqui vai para o manuscrito.

Uso:  python analise/teste_sintetico.py <pasta_de_saida>
"""
import subprocess
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.transform import from_origin
from shapely.geometry import box

BASE = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/teste_lai").resolve()
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(1)

gsd, x0, y0 = 0.02, 170000.0, 9710000.0     # 2 cm para o teste ficar leve
W, H = int(6 * 4.5 / gsd), int(4 * 5.5 / gsd)
img = np.zeros((3, H, W), dtype="uint8")
soil = np.array([150, 120, 95])
leaf = np.array([60, 120, 45])
polys, lai_rows = [], []
pid = 0
for row in range(4):
    for col in range(6):
        pid += 1
        lai = rng.uniform(1.0, 4.5)
        cover = 1 - np.exp(-0.6 * lai)                 # Beer-Lambert simplificado
        c0, r0 = int((col * 4.5 + 0.25) / gsd), int((row * 5.5 + 0.75) / gsd)
        c1, r1 = c0 + int(4 / gsd), r0 + int(4 / gsd)
        veg = rng.random((r1 - r0, c1 - c0)) < cover
        for b in range(3):
            base = np.where(veg, leaf[b] * (1.1 - 0.05 * lai), soil[b])
            img[b, r0:r1, c0:c1] = np.clip(base + rng.normal(0, 8, veg.shape), 0, 255)
        polys.append(box(x0 + c0 * gsd, y0 - r1 * gsd, x0 + c1 * gsd, y0 - r0 * gsd))
        lai_rows.append({"Parcela": pid, "Bloco": row + 1,
                         "TLA_cm2_planta": lai / 5.0 * 1e4 * rng.normal(1, 0.08)})

with rasterio.open(out / "ortomosaico_rgb.tif", "w", driver="GTiff", width=W, height=H, count=3,
                   dtype="uint8", crs="EPSG:31980", transform=from_origin(x0, y0, gsd, gsd)) as dst:
    dst.write(img)
gpd.GeoDataFrame({"Parcela": range(1, 25)}, geometry=polys, crs="EPSG:31980").to_file(out / "parcelas.gpkg")
pd.DataFrame(lai_rows).to_csv(out / "lai_ci202.csv", index=False)

cfg = yaml.safe_load(open(BASE / "config.yaml"))
cfg.update(ortomosaico=str(out / "ortomosaico_rgb.tif"), parcelas=str(out / "parcelas.gpkg"),
           lai_campo=str(out / "lai_ci202.csv"), densidade_plantas_m2=5.0,
           pasta_resultados=str(out / "resultados"))
yaml.safe_dump(cfg, open(out / "config_teste.yaml", "w"))

py = sys.executable
subprocess.run([py, BASE / "analise/01_extrair_indices.py", out / "config_teste.yaml"], check=True)
subprocess.run([py, BASE / "analise/02_modelagem.py", out / "config_teste.yaml", "--rapido",
                ], check=True)
