"""Análise preliminar com o que já existe: ExG por parcela (QGIS, DN 0-255)
contra a área foliar por planta do CI-202 (Livro2.xlsx / EXG.csv).

Só um índice, então só cabe regressão linear simples. Serve para dimensionar o
problema (força da relação e amplitude do LAI) antes de extrair os 13 índices.

Uso:  python analise/00_preliminar_exg.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut

BASE = Path(__file__).resolve().parents[1]
d = pd.read_csv(BASE / "dados" / "preliminar_exg_tla.csv")
x, y, n = d[["ExG_DN_QGIS"]], d["TLA_cm2_planta"].values, len(d)

r, p = stats.pearsonr(x.iloc[:, 0], y)
z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
lr = LinearRegression().fit(x, y)


def cv(splitter, groups=None):
    pred = np.empty(n)
    for tr, te in splitter.split(x, y, groups):
        pred[te] = LinearRegression().fit(x.iloc[tr], y[tr]).predict(x.iloc[te])
    rmse = np.sqrt(np.mean((y - pred) ** 2))
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2), rmse, 100 * rmse / y.mean()


lines = [
    f"n = {n} parcelas",
    f"Área foliar (cm² planta⁻¹): média {y.mean():.0f}, DP {y.std(ddof=1):.0f}, "
    f"CV {100 * y.std(ddof=1) / y.mean():.1f}%, amplitude {y.min():.0f}-{y.max():.0f}",
    f"ExG (DN): média {d.ExG_DN_QGIS.mean():.1f}, amplitude {d.ExG_DN_QGIS.min():.1f}-{d.ExG_DN_QGIS.max():.1f}",
    f"r de Pearson = {r:.3f} (IC95% {np.tanh(z - 1.96 * se):.3f} a {np.tanh(z + 1.96 * se):.3f}), p = {p:.4f}",
    f"Ajuste: TLA = {lr.intercept_:.1f} + {lr.coef_[0]:.2f} * ExG; R² = {lr.score(x, y):.3f}",
    "LOOCV:               R² = {:.3f}; RMSE = {:.0f} cm²; rRMSE = {:.1f}%".format(*cv(LeaveOneOut())),
    "Deixa um bloco fora: R² = {:.3f}; RMSE = {:.0f} cm²; rRMSE = {:.1f}%".format(
        *cv(LeaveOneGroupOut(), d["Bloco"].values)),
]
rho, prho = stats.spearmanr(d.ExG_DN_QGIS, y)
lines.append(f"Spearman rho = {rho:.3f}, p = {prho:.4f}")
out = BASE / "resultados" / "preliminar"
out.mkdir(parents=True, exist_ok=True)
(out / "preliminar_exg_tla.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
