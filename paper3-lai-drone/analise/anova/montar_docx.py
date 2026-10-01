"""Monta o .docx com tabelas e figuras (padrão Biochar / Field Crops Research)."""
import pandas as pd, sys
from docx import Document
from docx.shared import Pt, Mm
from docx.enum.text import WD_ALIGN_PARAGRAPH
R="resultados/anova_tukey"; out=sys.argv[1]
A=pd.read_csv(f"{R}/Table_ANOVA.csv"); M=pd.read_csv(f"{R}/Table_means_Tukey.csv"); G=pd.read_csv(f"{R}/Table_regression.csv")
doc=Document(); st=doc.styles["Normal"]; st.font.name="Arial"; st.font.size=Pt(10)
def par(t,b=False,i=False,size=None):
    p=doc.add_paragraph(); r=p.add_run(t); r.bold=b; r.italic=i
    if size: r.font.size=Pt(size)
    return p
def tabela(df):
    t=doc.add_table(rows=1,cols=len(df.columns)); t.style="Table Grid"
    for j,c in enumerate(df.columns): t.rows[0].cells[j].text=str(c); t.rows[0].cells[j].paragraphs[0].runs[0].bold=True
    for _,row in df.iterrows():
        cs=t.add_row().cells
        for j,c in enumerate(df.columns): cs[j].text="" if pd.isna(row[c]) else str(row[c])
    for row in t.rows:
        for c in row.cells:
            for p in c.paragraphs:
                for r in p.runs: r.font.size=Pt(8)
def fp(p): return "< 0.001" if p<0.001 else f"{p:.3f}"
doc.add_heading("Maize leaf area, SPAD and ExG under biochar and lime: ANOVA, Tukey and regression",1)
par("Data: table sent by Filipe on 1 Oct 2026 (24 plots). LAI = TLA × 10 plants m⁻² / 10 000. "
    "Note: this table differs from Livro2.xlsx/EXG.csv in plots 6, 19, 20 and 22 (TLA) and plot 19 (ExG); "
    "in plots 19 and 20 the reported LAI does not match TLA × 10 (3.24 vs 3.20; 2.30 vs 2.33). "
    "[VERIFICAR: confirm these values against the original CI-202 and QGIS records before submission.]", i=True, size=8)

doc.add_heading("Statistical analysis (Methods text)",2)
par("Data were analysed as a randomized complete block design in a 2 × 3 factorial arrangement "
    "(biochar: 0 and 12 Mg ha⁻¹; lime: 0, 75 and 100% of the recommended rate) with four blocks. "
    "Normality of residuals and homogeneity of variances were checked with the Shapiro–Wilk and Levene tests, respectively. "
    "When the biochar × lime interaction was significant (p < 0.05), biochar rates were compared within each lime rate and "
    "lime rates within each biochar rate; otherwise, main effects were compared. Means were separated by Tukey's HSD test "
    "(p < 0.05). Linear relationships among LAI, SPAD and ExG were evaluated by Pearson's correlation and ordinary least-squares "
    "regression (n = 24). Analyses were performed in R (R Core Team, 2025) with the packages emmeans, multcomp and car, "
    "and figures with ggplot2. [VERIFICAR: R version and citations]")

doc.add_heading("Table 1. Analysis of variance",2)
w=A[A.Source!="Residuals"].copy(); w["Fp"]=w.apply(lambda r: f"{r.F:.2f} ({fp(r.p)})",axis=1)
T1=w.pivot(index="Source",columns="Variable",values="Fp").reindex(["Block","Biochar","Lime","Biochar:Lime"])
T1.index=["Block (B)","Biochar (BC)","Lime (L)","BC × L"]; T1=T1[[c for c in ["LAI (m2 m-2)","SPAD index","ExG (DN)"]]]
res=A[A.Source=="Residuals"].set_index("Variable")
T1.loc["CV (%)"]=[f"{res.loc[c,'CV_pct']:.1f}" for c in T1.columns]
T1.loc["Shapiro–Wilk p"]=[f"{res.loc[c,'Shapiro_p']:.3f}" for c in T1.columns]
T1.loc["Levene p"]=[f"{res.loc[c,'Levene_p']:.3f}" for c in T1.columns]
T1.columns=["LAI (m² m⁻²)","SPAD index","ExG (DN)"]; T1=T1.reset_index().rename(columns={"index":"Source of variation"})
T1.insert(1,"df",["3","1","2","2","","",""])
tabela(T1)
par("Values are F (p). Residual df = 15. TLA gives the same F and p as LAI (LAI is TLA × constant).",size=8)

doc.add_heading("Table 2. Means ± standard deviation and Tukey groups",2)
li=M[M.Variable.str.startswith(("LAI","TLA"))].copy()
li["v"]=li.Mean_SD+" "+li.Tukey
T2=li.pivot(index="Level",columns="Variable",values="v").reindex(["BC0 L0","BC0 L75","BC0 L100","BC12 L0","BC12 L75","BC12 L100"])
T2.columns=["LAI (m² m⁻²)","TLA (cm² plant⁻¹)"]; tabela(T2.reset_index().rename(columns={"Level":"Treatment"}))
par("Biochar × lime interaction significant (p = 0.001). Uppercase letters compare biochar rates within each lime rate; "
    "lowercase letters compare lime rates within each biochar rate (Tukey, p < 0.05). n = 4. BC, biochar (Mg ha⁻¹); L, lime (% of recommended rate).",size=8)
se=M[~M.Variable.str.startswith(("LAI","TLA"))].copy(); se["v"]=se.Mean_SD+" "+se.Tukey
T3=se.pivot(index=["Effect","Level"],columns="Variable",values="v").reindex([("Biochar","BC0"),("Biochar","BC12"),("Lime","L0"),("Lime","L75"),("Lime","L100")])
T3.columns=["ExG (DN)","SPAD index"]; T3=T3[["SPAD index","ExG (DN)"]].reset_index()
tabela(T3)
par("No biochar × lime interaction (SPAD p = 0.841; ExG p = 0.110): main effects shown. Means followed by the same letter "
    "within a factor do not differ (Tukey, p < 0.05). n = 12 (biochar) and n = 8 (lime). 95% confidence half-widths are in Table_means_Tukey.csv.",size=8)

doc.add_heading("Table 3. Linear regressions (n = 24)",2)
G2=G[G.Response!="TLA"].copy()
G2["p"]=G2.p.apply(fp); G2=G2[["Response","Predictor","Equation","r","r_CI95","R2","R2_adj","F","p","RMSE"]]
G2.columns=["Response","Predictor","Equation","r","95% CI of r","R²","Adjusted R²","F","p","RMSE"]
tabela(G2)
par("In LAI ~ ExG + SPAD, only ExG was significant (partial p: ExG < 0.001; SPAD = 0.311; VIF = 1.23). RMSE in the units of the response.",size=8)

doc.add_heading("Figures",2)
doc.add_picture(f"{R}/Fig_bars_LAI_SPAD_ExG.png",width=Mm(140))
par("Fig. 1 (a) Leaf area index (LAI) of maize as affected by the interaction of biochar and lime rates; uppercase letters compare "
    "biochar rates within each lime rate and lowercase letters compare lime rates within each biochar rate. (b–e) Main effects of "
    "biochar and lime on SPAD index (b, c) and excess green index, ExG (d, e). Bars are means and error bars are ± standard deviation. "
    "Means followed by the same letter do not differ by Tukey's test (p < 0.05).",size=8)
doc.add_picture(f"{R}/Fig_regression.png",width=Mm(170))
par("Fig. 2 Relationships between (a) LAI and ExG, (b) LAI and SPAD index and (c) SPAD index and ExG in maize plots (n = 24). "
    "Lines are ordinary least-squares fits with 95% confidence bands. Symbols: biochar rate (shape) and lime rate (fill).",size=8)

doc.add_heading("Format notes",2)
par("Figures are exported at 600 dpi (PNG) and as vector PDF, Arial/sans-serif, single-column width ≤ 90 mm and double-column ≤ 180 mm, "
    "which fits Elsevier (Field Crops Research) and Springer (Biochar) artwork guidelines [VERIFICAR: current author guidelines of each journal]. "
    "Units follow SI: Mg ha⁻¹, m² m⁻², cm² plant⁻¹. SPAD and ExG are dimensionless (ExG computed on 0–255 digital numbers).",size=8)
doc.save(out)
