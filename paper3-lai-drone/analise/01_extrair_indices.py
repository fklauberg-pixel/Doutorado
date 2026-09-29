"""Etapa 1: extrai os índices RGB por parcela a partir do ortomosaico.

Uso:  python analise/01_extrair_indices.py [config.yaml]
Saída: resultados/indices_por_parcela.csv

Para cada parcela: aplica bordadura interna, calcula os índices pixel a pixel e
tira a média (a) de todos os pixels da parcela e (b) só dos pixels de vegetação,
separados do solo por limiar de Otsu (1979) sobre o ExG. Também registra a
cobertura do dossel (CC, fração de pixels de vegetação).
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.mask import mask

sys.path.insert(0, str(Path(__file__).parent))
from indices import compute_indices  # noqa: E402

BASE = Path(__file__).resolve().parents[1]


def otsu(values, bins=256):
    values = values[np.isfinite(values)]
    hist, edges = np.histogram(values, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers) / np.maximum(w0, 1)
    m1 = (np.sum(hist * centers) - np.cumsum(hist * centers)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return centers[np.argmax(between)]


def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    ortho = BASE / cfg["ortomosaico"]
    plots = gpd.read_file(BASE / cfg["parcelas"])
    idcol = cfg["campo_id_parcela"]

    with rasterio.open(ortho) as src:
        if plots.crs != src.crs:
            plots = plots.to_crs(src.crs)
        if plots.crs.is_geographic and cfg["buffer_interno_m"]:
            raise SystemExit("Reprojete as parcelas para um CRS métrico (ex.: EPSG:31980) "
                             "para aplicar a bordadura em metros.")
        geoms = plots.geometry.buffer(-cfg["buffer_interno_m"])

        # Limiar de solo/vegetação calculado uma única vez sobre todas as parcelas
        exg_all, pix = [], {}
        for pid, geom in zip(plots[idcol], geoms):
            arr, _ = mask(src, [geom], crop=True, filled=False, indexes=[1, 2, 3])
            valid = ~np.any(np.ma.getmaskarray(arr), axis=0)
            R, G, B = (arr[i].data[valid].astype("float64") for i in range(3))
            pix[pid] = (R, G, B)
            exg_all.append(compute_indices(R, G, B)["ExG"])
    thr = otsu(np.concatenate(exg_all)) if cfg["mascara_solo"] == "otsu_exg" else -np.inf

    rows = []
    for pid, (R, G, B) in pix.items():
        idx = compute_indices(R, G, B)
        veg = idx["ExG"] > thr
        row = {idcol: pid, "n_pixels": R.size, "CC": veg.mean()}
        for name, v in idx.items():
            row[f"{name}_todo"] = np.nanmean(v)
            row[f"{name}_veg"] = np.nanmean(v[veg]) if veg.any() else np.nan
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(idcol)
    out.attrs["limiar_otsu_exg"] = thr
    dest = BASE / cfg.get("pasta_resultados", "resultados") / "indices_por_parcela.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False)
    print(f"Limiar Otsu (ExG normalizado): {thr:.4f}")
    print(f"{len(out)} parcelas -> {dest}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else BASE / "config.yaml")
