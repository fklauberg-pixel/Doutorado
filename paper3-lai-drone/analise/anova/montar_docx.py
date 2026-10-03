"""Monta o .docx com tabelas e figuras (padrão Biochar / Field Crops Research).
Lê as saídas de analise/anova/anova_tukey_regressao.R (dados do Livro2)."""
import pandas as pd, sys
from docx import Document
from docx.shared import Pt, Mm

R = "resultados/anova_tukey"; out = sys.argv[1]
A = pd.read_csv(f"{R}/Table_ANOVA.csv"); M = pd.read_csv(f"{R}/Table_means_Tukey.csv"); G = pd.read_csv(f"{R}/Table_regression.csv")
VAR = {"TLA (cm2 plant-1)": "TLA (cm² plant⁻¹)", "LAI (m2 m-2)": "LAI (m² m⁻²)", "SPAD index": "SPAD index", "ExG (DN)": "ExG (DN)"}
ORD = list(VAR)
TRAT = ("BC0 and BC12, 0 and 12 Mg ha⁻¹ of açaí biochar; L0, L75 and L100, 0, 75 and 100% of the recommended "
        "dolomitic lime rate.")

doc = Document(); st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(10)
def par(t, b=False, i=False, size=None):
    p = doc.add_paragraph(); r = p.add_run(t); r.bold = b; r.italic = i
    if size: r.font.size = Pt(size)
def tabela(df):
    t = doc.add_table(rows=1, cols=len(df.columns)); t.style = "Table Grid"
    for j, c in enumerate(df.columns):
        t.rows[0].cells[j].text = str(c); t.rows[0].cells[j].paragraphs[0].runs[0].bold = True
    for _, row in df.iterrows():
        cs = t.add_row().cells
        for j, c in enumerate(df.columns): cs[j].text = "" if pd.isna(row[c]) else str(row[c])
    for row in t.rows:
        for c in row.cells:
            for p in c.paragraphs:
                for r in p.runs: r.font.size = Pt(8)
def fp(p): return "< 0.001" if p < 0.001 else f"{p:.3f}"

doc.add_heading("Maize leaf area, LAI, SPAD and ExG under biochar and lime: ANOVA, Tukey and regression", 1)
par("Data: Livro2.xlsx (24 plots; TLA from the CI-202, SPAD, and ExG extracted in QGIS), the dataset confirmed as valid "
    "on 3 Oct 2026. LAI = TLA × 10 plants m⁻² / 10 000.", i=True, size=8)

doc.add_heading("Statistical analysis (Methods text)", 2)
par("Data were analysed as a randomized complete block design in a 2 × 3 factorial arrangement "
    "(biochar: 0 and 12 Mg ha⁻¹; lime: 0, 75 and 100% of the recommended rate) with four blocks. "
    "Normality of residuals and homogeneity of variances were checked with the Shapiro–Wilk and Levene tests, respectively. "
    "When the biochar × lime interaction was significant (p < 0.05), biochar rates were compared within each lime rate and "
    "lime rates within each biochar rate; otherwise, main effects were compared. Means were separated by Tukey's HSD test "
    "(p < 0.05). Linear relationships among leaf area, SPAD and ExG were evaluated by Pearson's correlation and ordinary "
    "least-squares regression (n = 24). Analyses were performed in R (R Core Team, 2025) with the packages emmeans, "
    "multcomp and car, and figures with ggplot2. [VERIFICAR: R version and citations]")

doc.add_heading("Table 1. Analysis of variance", 2)
w = A[A.Source != "Residuals"].copy(); w["Fp"] = w.apply(lambda r: f"{r.F:.2f} ({fp(r.p)})", axis=1)
T1 = w.pivot(index="Source", columns="Variable", values="Fp").reindex(["Block", "Biochar", "Lime", "Biochar:Lime"])[ORD]
T1.index = ["Block", "Biochar (BC)", "Lime (L)", "BC × L"]
res = A[A.Source == "Residuals"].set_index("Variable")
for lab, col, f in [("CV (%)", "CV_pct", "{:.1f}"), ("Shapiro–Wilk p", "Shapiro_p", "{:.3f}"), ("Levene p", "Levene_p", "{:.3f}")]:
    T1.loc[lab] = [f.format(res.loc[c, col]) for c in ORD]
T1.columns = [VAR[c] for c in ORD]; T1 = T1.reset_index().rename(columns={"index": "Source of variation"})
T1.insert(1, "df", ["3", "1", "2", "2", "", "", ""])
tabela(T1)
par("Values are F (p). Residual df = 15. LAI is TLA multiplied by a constant, so both have the same F, p and CV.", size=8)

doc.add_heading("Table 2. Means ± standard deviation and Tukey groups", 2)
M["v"] = M.apply(lambda r: r.Mean_SD if r.Tukey == "ns" else f"{r.Mean_SD} {r.Tukey}", axis=1)
T2 = M.pivot(index=["Effect", "Level"], columns="Variable", values="v")
T2 = T2.reindex([("Biochar", "BC0"), ("Biochar", "BC12"), ("Lime", "L0"), ("Lime", "L75"), ("Lime", "L100")])[ORD]
T2.columns = [VAR[c] for c in ORD]; T2 = T2.reset_index().rename(columns={"Effect": "Factor", "Level": "Rate"})
tabela(T2)
par("No biochar × lime interaction was significant (Table 1), so main effects are shown. Within each factor, means followed "
    "by the same letter do not differ (Tukey, p < 0.05); no letters, F test not significant. n = 12 (biochar) and n = 8 (lime). "
    + TRAT + " 95% confidence half-widths are in Table_means_Tukey.csv.", size=8)

doc.add_heading("Table 3. Linear regressions (n = 24)", 2)
G2 = G.copy(); G2["p"] = G2.p.apply(fp)
G2["Response"] = G2.Response.map({"TLA": "TLA", "LAI": "LAI", "SPAD": "SPAD"})
G2 = G2[["Response", "Predictor", "Equation", "r", "r_CI95", "R2", "R2_adj", "F", "p", "RMSE"]]
G2.columns = ["Response", "Predictor", "Equation", "r", "95% CI of r", "R²", "Adjusted R²", "F", "p", "RMSE"]
tabela(G2)
par("RMSE in the units of the response (TLA, cm² plant⁻¹; LAI, m² m⁻²; SPAD, index units). In TLA (or LAI) ~ ExG + SPAD, "
    "neither predictor was significant (partial p: ExG = 0.093; SPAD = 0.773; VIF = 1.23).", size=8)

doc.add_heading("Figures", 2)
doc.add_picture(f"{R}/Fig_bars_TLA_SPAD_ExG.png", width=Mm(140))
par("Fig. 1 Main effects of biochar and lime on (a, b) leaf area per plant (TLA), (c, d) SPAD index and (e, f) excess green "
    "index (ExG) of maize. Bars are means and error bars are ± standard deviation (n = 12 for biochar and n = 8 for lime). "
    "Means followed by the same letter do not differ by Tukey's test (p < 0.05); no letters, F test not significant. "
    "F and p of each factor and of the biochar × lime interaction are shown above each panel. " + TRAT, size=8)
doc.add_picture(f"{R}/Fig_regression.png", width=Mm(170))
par("Fig. 2 Relationships between (a) leaf area per plant (TLA) and ExG, (b) TLA and SPAD index and (c) SPAD index and ExG "
    "in maize plots (n = 24). Lines are ordinary least-squares fits with 95% confidence bands. Symbol shape, biochar rate; "
    "fill, lime rate. " + TRAT, size=8)
doc.add_picture(f"{R}/Fig_bars_LAI_SPAD_ExG.png", width=Mm(140))
par("Fig. S1 Same as Fig. 1 with leaf area index (LAI = TLA × 10 plants m⁻²) in panels a and b.", size=8)

doc.add_heading("Format notes", 2)
par("Figures are exported at 600 dpi (PNG) and as vector PDF, sans-serif font, double-column width (≤ 180 mm), "
    "in line with Elsevier (Field Crops Research) and Springer (Biochar) artwork guidelines [VERIFICAR: current author "
    "guidelines of each journal]. Units follow SI. SPAD and ExG are dimensionless (ExG computed on 0–255 digital numbers).", size=8)
doc.save(out)
