"""Etapa 2: compara regressão convencional e Random Forest na estimativa do LAI.

Uso:
  python analise/02_modelagem.py [config.yaml] [--sufixo todo|veg] [--rapido]

Entradas: resultados/indices_por_parcela.csv (etapa 1) e o CSV do CI-202.
Saídas em resultados/<sufixo>/: tabelas CSV e figuras PNG.

Modelos (todos ajustados SÓ com os dados de treino de cada partição, para que a
seleção de variáveis não vaze informação da validação):
  Nulo     média do treino (linha de base)
  ULR      regressão linear simples com o índice de maior |r| no treino
  MLR      regressão linear múltipla com índices |r| > limiar (Du et al., 2022),
           podados por colinearidade (|r| entre preditores > 0,95) e limitados
           a 3 preditores, pois n = 24
  RF_sel   Random Forest com os mesmos índices do MLR
  RF_todos Random Forest com todos os índices
Hiperparâmetros do RF (max_features, min_samples_leaf) escolhidos pelo erro
out-of-bag dentro do treino.

Esquemas de validação:
  kfold    k-fold repetido (principal)
  loocv    leave-one-out
  lobo     leave-one-block-out (deixa um bloco inteiro de fora)
  split    N partições aleatórias 70/30, como em Du et al. (2022)
"""
import argparse
import sys
import warnings
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy import stats  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.ensemble import RandomForestRegressor  # noqa: E402
from sklearn.inspection import permutation_importance  # noqa: E402
from sklearn.linear_model import LinearRegression  # noqa: E402
from sklearn.model_selection import (LeaveOneGroupOut, LeaveOneOut,  # noqa: E402
                                     RepeatedKFold, ShuffleSplit)

warnings.filterwarnings("ignore", category=UserWarning)
BASE = Path(__file__).resolve().parents[1]
MODELOS = ["Nulo", "ULR", "MLR", "RF_sel", "RF_todos"]


# ----------------------------------------------------------------------------
# Seleção de variáveis e modelos
# ----------------------------------------------------------------------------
def corr_with_y(X, y):
    return pd.Series({c: np.corrcoef(X[c], y)[0, 1] for c in X.columns}).fillna(0)


def select_features(X, y, thr, max_k=3, collin=0.95):
    r = corr_with_y(X, y).abs().sort_values(ascending=False)
    cand = list(r[r > thr].index) or list(r.index[:1])
    keep = []
    for c in cand:
        if all(abs(np.corrcoef(X[c], X[k])[0, 1]) <= collin for k in keep):
            keep.append(c)
        if len(keep) == max_k:
            break
    return keep


def fit_rf(X, y, cfg_rf, seed, fast):
    grid = list(product(cfg_rf["max_features"], cfg_rf["min_samples_leaf"]))
    if fast:
        grid = grid[::3]
    best, best_err = None, np.inf
    for mf, msl in grid:
        rf = RandomForestRegressor(n_estimators=150, max_features=mf, min_samples_leaf=msl,
                                   oob_score=True, bootstrap=True, random_state=seed, n_jobs=1)
        rf.fit(X, y)
        err = np.mean((rf.oob_prediction_ - y) ** 2)
        if err < best_err:
            best, best_err = (mf, msl), err
    rf = RandomForestRegressor(n_estimators=cfg_rf["n_estimators"], max_features=best[0],
                               min_samples_leaf=best[1], random_state=seed, n_jobs=1)
    return rf.fit(X, y), best


def fit_predict_all(Xtr, ytr, Xte, cfg, seed, fast):
    thr = cfg["limiar_r_selecao"]
    pred, info = {}, {}
    pred["Nulo"] = np.full(len(Xte), ytr.mean())

    best1 = corr_with_y(Xtr, ytr).abs().idxmax()
    pred["ULR"] = LinearRegression().fit(Xtr[[best1]], ytr).predict(Xte[[best1]])
    info["ULR"] = [best1]

    sel = select_features(Xtr, ytr, thr)
    pred["MLR"] = LinearRegression().fit(Xtr[sel], ytr).predict(Xte[sel])
    info["MLR"] = sel

    rf, hp = fit_rf(Xtr[sel].values, ytr, cfg["random_forest"], seed, fast)
    pred["RF_sel"] = rf.predict(Xte[sel].values)
    info["RF_sel"] = sel + [f"hp={hp}"]

    rf, hp = fit_rf(Xtr.values, ytr, cfg["random_forest"], seed, fast)
    pred["RF_todos"] = rf.predict(Xte.values)
    info["RF_todos"] = [f"hp={hp}"]
    return pred, info


# ----------------------------------------------------------------------------
# Métricas
# ----------------------------------------------------------------------------
def metrics(y, p):
    y, p = np.asarray(y), np.asarray(p)
    sse = np.sum((y - p) ** 2)
    rmse = np.sqrt(sse / len(y))
    return {
        "R2": 1 - sse / np.sum((y - y.mean()) ** 2),
        "RMSE": rmse,
        "rRMSE_pct": 100 * rmse / y.mean(),
        "MAE": np.mean(np.abs(y - p)),
        "vies": np.mean(p - y),
    }


def run_pooled(X, y, splitter, groups, cfg, fast, n_rep):
    """Para kfold/loocv/lobo: junta as predições de todas as dobras de uma
    repetição e calcula as métricas sobre o conjunto completo."""
    n = len(y)
    splits = list(splitter.split(X, y, groups))
    per_rep = len(splits) // n_rep
    rows, oof_first, chosen = [], None, []
    fold_mse = {m: [] for m in MODELOS}
    seed0 = cfg["validacao"]["semente"]
    fitted = Parallel(n_jobs=-1)(
        delayed(fit_predict_all)(X.iloc[tr], y[tr], X.iloc[te], cfg, seed0 + i, fast)
        for i, (tr, te) in enumerate(splits))
    for rep in range(n_rep):
        oof = {m: np.full(n, np.nan) for m in MODELOS}
        for k, (tr, te) in enumerate(splits[rep * per_rep:(rep + 1) * per_rep]):
            pred, info = fitted[rep * per_rep + k]
            for m in MODELOS:
                oof[m][te] = pred[m]
                fold_mse[m].append(np.mean((y[te] - pred[m]) ** 2))
            chosen.append(info)
        for m in MODELOS:
            rows.append({"repeticao": rep, "modelo": m, **metrics(y, oof[m])})
        if rep == 0:
            oof_first = oof
    return pd.DataFrame(rows), oof_first, chosen, fold_mse


def run_split(X, y, cfg, fast):
    v = cfg["validacao"]
    ss = ShuffleSplit(n_splits=v["splits_70_30"], test_size=0.3, random_state=v["semente"])
    splits = list(ss.split(X))
    fitted = Parallel(n_jobs=-1)(
        delayed(fit_predict_all)(X.iloc[tr], y[tr], X.iloc[te], cfg, v["semente"] + i, fast)
        for i, (tr, te) in enumerate(splits))
    rows = []
    for i, ((tr, te), (pred, _)) in enumerate(zip(splits, fitted)):
        for m in MODELOS:
            rows.append({"repeticao": i, "modelo": m, **metrics(y[te], pred[m])})
    return pd.DataFrame(rows)


def corrected_t(diffs, n_train, n_test):
    """Teste t corrigido para reamostragem (Nadeau & Bengio, 2003)."""
    d = np.asarray(diffs)
    k = len(d)
    var = np.var(d, ddof=1)
    t = d.mean() / np.sqrt((1 / k + n_test / n_train) * var)
    return t, 2 * stats.t.sf(abs(t), df=k - 1)


def summarise(df):
    g = df.groupby("modelo")[["R2", "RMSE", "rRMSE_pct", "MAE", "vies"]]
    s = g.mean().add_suffix("_media").join(g.std().add_suffix("_dp"))
    return s.loc[[m for m in MODELOS if m in s.index]].round(4)


# ----------------------------------------------------------------------------
# Figuras
# ----------------------------------------------------------------------------
INK, MUTED, MARK = "#1f2328", "#8c959f", "#2f6db5"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(alpha=0.25, lw=0.6)


def fig_corr(r_table, dest):
    r_table = r_table.sort_values("r")
    fig, ax = plt.subplots(figsize=(5, 0.28 * len(r_table) + 1))
    ax.barh(r_table.index, r_table["r"], color=MARK, height=0.6)
    ax.errorbar(r_table["r"], r_table.index,
                xerr=[r_table["r"] - r_table["IC95_inf"], r_table["IC95_sup"] - r_table["r"]],
                fmt="none", ecolor=INK, lw=0.8, capsize=2)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("r de Pearson com o LAI (IC 95%)")
    style(ax)
    fig.tight_layout()
    fig.savefig(dest, dpi=300)
    plt.close(fig)


def fig_obs_pred(y, oof, label_y, dest):
    ms = ["ULR", "MLR", "RF_sel", "RF_todos"]
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.3), sharex=True, sharey=True)
    lo, hi = min(y.min(), min(oof[m].min() for m in ms)), max(y.max(), max(oof[m].max() for m in ms))
    pad = 0.05 * (hi - lo)
    for ax, m in zip(axes, ms):
        met = metrics(y, oof[m])
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=MUTED, lw=1, ls="--")
        ax.scatter(y, oof[m], s=36, color=MARK, edgecolor="white", linewidth=1)
        ax.set_title(m, fontsize=10, color=INK)
        ax.text(0.04, 0.96, f"R² = {met['R2']:.2f}\nRMSE = {met['RMSE']:.3g}",
                transform=ax.transAxes, va="top", fontsize=9, color=INK)
        ax.set_xlabel(f"{label_y} observado (CI-202)")
        style(ax)
    axes[0].set_ylabel(f"{label_y} estimado (validação cruzada)")
    fig.tight_layout()
    fig.savefig(dest, dpi=300)
    plt.close(fig)


def fig_rmse_box(df, dest):
    ms = [m for m in MODELOS if m in df["modelo"].unique()]
    data = [df.loc[df.modelo == m, "RMSE"] for m in ms]
    fig, ax = plt.subplots(figsize=(5.5, 3.3))
    ax.boxplot(data, tick_labels=ms, widths=0.5, patch_artist=True,
               boxprops=dict(facecolor="#dbe6f4", edgecolor=MARK),
               medianprops=dict(color=INK), whiskerprops=dict(color=MUTED),
               capprops=dict(color=MUTED), flierprops=dict(markeredgecolor=MUTED, markersize=3))
    ax.set_ylabel("RMSE de validação cruzada")
    style(ax)
    fig.tight_layout()
    fig.savefig(dest, dpi=300)
    plt.close(fig)


def fig_importance(imp, dest):
    imp = imp.sort_values("media")
    fig, ax = plt.subplots(figsize=(5, 0.28 * len(imp) + 1))
    ax.barh(imp.index, imp["media"], xerr=imp["dp"], color=MARK, height=0.6,
            error_kw=dict(ecolor=INK, lw=0.8, capsize=2))
    ax.set_xlabel("Importância por permutação (queda no R²)")
    style(ax)
    fig.tight_layout()
    fig.savefig(dest, dpi=300)
    plt.close(fig)


# ----------------------------------------------------------------------------
def load_data(cfg, sufixo, indices_csv=None):
    idcol = cfg["campo_id_parcela"]
    ind = pd.read_csv(indices_csv or BASE / cfg.get("pasta_resultados", "resultados") / "indices_por_parcela.csv")
    lai = pd.read_csv(BASE / cfg["lai_campo"])
    df = lai.merge(ind, on=idcol, how="inner", validate="one_to_one")
    if len(df) != len(lai):
        print(f"AVISO: {len(lai) - len(df)} parcelas do CI-202 sem índices correspondentes")
    feats = [c for c in df.columns if c.endswith(f"_{sufixo}")]
    if "CC" in df.columns:
        feats.append("CC")
    X = df[feats].rename(columns=lambda c: c.replace(f"_{sufixo}", ""))
    X = X.loc[:, X.std() > 0]
    area = df[cfg["coluna_area_foliar"]].astype(float).values
    dens = cfg.get("densidade_plantas_m2")
    if dens:
        y, label = area / 1e4 * dens, "LAI (m² m⁻²)"
    else:
        y, label = area, "Área foliar (cm² planta⁻¹)"
        print("AVISO: densidade de plantas não informada; usando área foliar por planta "
              "(r e R² não mudam; RMSE fica em cm² planta⁻¹).")
    groups = df["Bloco"].values if "Bloco" in df.columns else None
    return df, X, y, label, groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", nargs="?", default=str(BASE / "config.yaml"))
    ap.add_argument("--sufixo", default="todo", choices=["todo", "veg"])
    ap.add_argument("--rapido", action="store_true", help="menos repetições, para testar")
    ap.add_argument("--indices", default=None, help="CSV de índices alternativo")
    ap.add_argument("--saida", default=None)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    v = cfg["validacao"]
    if a.rapido:
        v["repeticoes"], v["splits_70_30"] = 5, 20
    out = Path(a.saida) if a.saida else BASE / cfg.get("pasta_resultados", "resultados") / a.sufixo
    out.mkdir(parents=True, exist_ok=True)

    df, X, y, label, groups = load_data(cfg, a.sufixo, a.indices)
    n = len(y)
    print(f"n = {n} parcelas, {X.shape[1]} preditores: {', '.join(X.columns)}")

    # 1. Correlações (descritivas, com todos os dados)
    rows = []
    for c in X.columns:
        r, p = stats.pearsonr(X[c], y)
        z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
        rows.append({"indice": c, "r": r, "p": p,
                     "IC95_inf": np.tanh(z - 1.96 * se), "IC95_sup": np.tanh(z + 1.96 * se)})
    rt = pd.DataFrame(rows).set_index("indice").sort_values("r", key=abs, ascending=False)
    rt.round(4).to_csv(out / "tab_correlacoes.csv")
    fig_corr(rt, out / "fig_correlacoes.png")

    # 2. Regressão linear simples índice a índice (LOOCV), tabela descritiva
    rows = []
    for c in X.columns:
        p = np.empty(n)
        for tr, te in LeaveOneOut().split(X):
            p[te] = LinearRegression().fit(X.iloc[tr][[c]], y[tr]).predict(X.iloc[te][[c]])
        lr = LinearRegression().fit(X[[c]], y)
        rows.append({"indice": c, "a": lr.intercept_, "b": lr.coef_[0],
                     "R2_ajuste": lr.score(X[[c]], y), **{f"{k}_LOOCV": v_ for k, v_ in metrics(y, p).items()}})
    pd.DataFrame(rows).set_index("indice").sort_values("RMSE_LOOCV").round(4).to_csv(out / "tab_ULR_por_indice.csv")

    # 3. Validação cruzada dos cinco modelos
    rkf = RepeatedKFold(n_splits=v["kfold"], n_repeats=v["repeticoes"], random_state=v["semente"])
    res_kf, oof_kf, chosen_kf, fold_mse = run_pooled(X, y, rkf, None, cfg, a.rapido, v["repeticoes"])
    res_lo, oof_lo, _, _ = run_pooled(X, y, LeaveOneOut(), None, cfg, a.rapido, 1)
    tabs = {"kfold_repetido": summarise(res_kf), "LOOCV": summarise(res_lo)}
    if groups is not None and len(np.unique(groups)) > 2:
        res_bl, _, _, _ = run_pooled(X, y, LeaveOneGroupOut(), groups, cfg, a.rapido, 1)
        tabs["deixa_um_bloco_fora"] = summarise(res_bl)
    res_sp = run_split(X, y, cfg, a.rapido)
    tabs["split_70_30"] = summarise(res_sp)
    pd.concat(tabs, names=["esquema"]).to_csv(out / "tab_desempenho_modelos.csv")
    res_kf.to_csv(out / "kfold_metricas_por_repeticao.csv", index=False)

    # 4. Comparação pareada RF x melhor linear (mesmas dobras)
    n_te = n // v["kfold"]
    comp = []
    for rf in ("RF_sel", "RF_todos"):
        for lin in ("ULR", "MLR"):
            d = np.array(fold_mse[lin]) - np.array(fold_mse[rf])
            t, p = corrected_t(d, n - n_te, n_te)
            r_rf = res_kf[res_kf.modelo == rf]["RMSE"].values
            r_li = res_kf[res_kf.modelo == lin]["RMSE"].values
            comp.append({"comparacao": f"{rf} vs {lin}", "dif_RMSE_media": (r_li - r_rf).mean(),
                         "pct_repeticoes_RF_melhor": 100 * np.mean(r_rf < r_li),
                         "t_corrigido": t, "p_corrigido": p})
    pd.DataFrame(comp).round(4).to_csv(out / "tab_comparacao_RF_vs_linear.csv", index=False)

    # 5. Frequência com que cada índice foi selecionado nas dobras
    sel = pd.Series([f for c in chosen_kf for f in c["MLR"]]).value_counts() / len(chosen_kf) * 100
    sel.round(1).rename("pct_dobras").to_csv(out / "tab_frequencia_selecao_MLR.csv")
    ulr = pd.Series([c["ULR"][0] for c in chosen_kf]).value_counts() / len(chosen_kf) * 100
    ulr.round(1).rename("pct_dobras").to_csv(out / "tab_frequencia_selecao_ULR.csv")

    # 6. Modelos finais com todos os dados
    final_sel = select_features(X, y, cfg["limiar_r_selecao"])
    mlr = LinearRegression().fit(X[final_sel], y)
    with open(out / "modelos_finais.txt", "w") as f:
        f.write(f"Resposta: {label}; n = {n}\n")
        f.write(f"MLR: y = {mlr.intercept_:.4f} " +
                " ".join(f"{b:+.4f}*{c}" for b, c in zip(mlr.coef_, final_sel)) +
                f"   R2_ajuste = {mlr.score(X[final_sel], y):.3f}\n")
        if len(final_sel) > 1:
            Z = X[final_sel].values
            for j, c in enumerate(final_sel):
                others = np.delete(Z, j, axis=1)
                r2 = LinearRegression().fit(others, Z[:, j]).score(others, Z[:, j])
                f.write(f"  VIF {c}: {1 / max(1 - r2, 1e-9):.1f}\n")
        rf, hp = fit_rf(X.values, y, cfg["random_forest"], v["semente"], a.rapido)
        f.write(f"RF_todos: max_features, min_samples_leaf = {hp}\n")
    imp = permutation_importance(rf, X.values, y, n_repeats=50, random_state=v["semente"], scoring="r2")
    imp = pd.DataFrame({"media": imp.importances_mean, "dp": imp.importances_std}, index=X.columns)
    imp.round(4).to_csv(out / "tab_importancia_RF.csv")
    fig_importance(imp, out / "fig_importancia_RF.png")

    fig_obs_pred(y, oof_lo, label.split(" (")[0], out / "fig_observado_vs_estimado_LOOCV.png")
    fig_rmse_box(res_kf, out / "fig_RMSE_kfold.png")
    print(pd.concat(tabs, names=["esquema"])[["R2_media", "RMSE_media", "rRMSE_pct_media"]])
    print(f"Resultados em {out}")


if __name__ == "__main__":
    main()
