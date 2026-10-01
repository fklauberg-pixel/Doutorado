"""Análise da tabela por parcela (TLA, LAI, SPAD, ExG) enviada pelo Filipe em 01/10/2026,
com GY, AGB e H do Livro2. Versão A = tabela do Filipe (parcelas 20 e 22 com TLA 2333,30 e 2161,96);
versão B = Livro2 (2933,30 e 2861,96)."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, statsmodels.api as sm
from scipy import stats
L=pd.read_excel("/mnt/project-files/Livro2.xlsx")
base=pd.DataFrame(dict(Parcela=L.Parcela,BC=L.BC,Lm=L.L,BL=L.BL.astype(str),SPAD=L.SPAD,ExG=L.ExG,
    GY=L.GY,AGBse=L.AGBse,AGBce=L.AGBce,H=L.H,Dc=L.Dc,NF=L.n_folhas,TLA_B=L["TLA cm2 plant-1"]))
base["TLA_A"]=base.TLA_B; base.loc[base.Parcela==20,"TLA_A"]=2333.30; base.loc[base.Parcela==22,"TLA_A"]=2161.96
base["Lnum"]=base.Lm.str[1:].astype(int)
out=[]
def p(*a): s=" ".join(str(x) for x in a); print(s); out.append(s)
for ver in ["A","B"]:
    d=base.copy(); d["LAI"]=d[f"TLA_{ver}"]*10/1e4; d["CCC"]=d.LAI*d.SPAD   # conteúdo de clorofila do dossel (proxy)
    p(f"\n===== Versão {ver} =====")
    p(f"LAI média {d.LAI.mean():.2f}, CV {100*d.LAI.std()/d.LAI.mean():.1f}%")
    for v in ["LAI","SPAD","ExG","CCC","GY","AGBce","H"]:
        a=sm.stats.anova_lm(smf.ols(f"{v}~C(BL)+C(BC)*C(Lm)",d).fit())
        m=d.groupby(["BC","Lm"])[v].mean().unstack()[["L0","L75","L100"]]
        p(f"{v:6s} F_BC={a.F.iloc[1]:.1f} p={a['PR(>F)'].iloc[1]:.4f} | F_L={a.F.iloc[2]:.1f} p={a['PR(>F)'].iloc[2]:.4f} | F_int={a.F.iloc[3]:.2f} p={a['PR(>F)'].iloc[3]:.3f} | CVexp={100*np.sqrt(a.mean_sq.iloc[-1])/d[v].mean():.1f}%")
        p("       médias BC0:",m.loc['BC0'].round(2).tolist()," BC12:",m.loc['BC12'].round(2).tolist())
    p("Correlações de Pearson (n=24):")
    vs=["LAI","SPAD","ExG","CCC","GY","AGBce","H"]
    CR=d[vs].corr().round(2); p(CR.to_string())
    # ExG explicado por SPAD e LAI (padronizado)
    z=d[vs].apply(stats.zscore)
    f=smf.ols("ExG~SPAD+LAI",z).fit(); p(f"ExG ~ SPAD + LAI (padronizado): beta_SPAD={f.params.SPAD:.2f} (p={f.pvalues.SPAD:.3f}), beta_LAI={f.params.LAI:.2f} (p={f.pvalues.LAI:.3f}), R2={f.rsquared:.2f}")
    f2=smf.ols("ExG~C(BL)+SPAD+LAI",d).fit(); p(f"   com bloco: p_SPAD={f2.pvalues.SPAD:.3f}, p_LAI={f2.pvalues.LAI:.3f}, R2={f2.rsquared:.2f}")
    for v in ["LAI","SPAD","CCC"]:
        f=smf.ols(f"ExG~{v}",d).fit(); p(f"ExG ~ {v}: R2={f.rsquared:.3f}, AIC={f.aic:.1f}")
    # médias de tratamento
    mt=d.groupby(["BC","Lm"])[["LAI","SPAD","ExG","CCC","GY"]].mean()
    for a_,b_ in [("ExG","LAI"),("ExG","SPAD"),("ExG","CCC"),("GY","LAI"),("GY","SPAD"),("GY","ExG"),("GY","CCC")]:
        r,pv=stats.pearsonr(mt[a_],mt[b_]); p(f"médias de tratamento r({a_},{b_})={r:.2f} p={pv:.3f}")
    # rendimento: o que explica GY
    for form in ["GY~LAI","GY~SPAD","GY~ExG","GY~CCC","GY~LAI+SPAD"]:
        f=smf.ols(form,d).fit(); p(f"{form}: R2={f.rsquared:.2f}, AIC={f.aic:.1f}, p={f.f_pvalue:.4f}")
    # contraste linear de calcário (dose)
    for v in ["SPAD","ExG","LAI"]:
        f=smf.ols(f"{v}~C(BL)+C(BC)+Lnum",d).fit(); p(f"{v}: inclinação por % de dose de calcário = {f.params.Lnum:.4f} (p={f.pvalues.Lnum:.4f})")
open("resultados/tabela_q1.txt","w").write("\n".join(out))
# correlação dentro de tratamento (resíduos de bloco + tratamento): o drone enxerga variação entre parcelas além do tratamento?
for ver in ["A","B"]:
    d=base.copy(); d["LAI"]=d[f"TLA_{ver}"]*10/1e4; d["CCC"]=d.LAI*d.SPAD
    R=pd.DataFrame({v:smf.ols(f"{v}~C(BL)+C(BC)*C(Lm)",d).fit().resid for v in ["LAI","SPAD","ExG","CCC","H","GY"]})
    p(f"\nVersão {ver}, correlação dos resíduos (dentro de tratamento):"); p(R.corr().round(2).to_string())
open("resultados/tabela_q1.txt","w").write("\n".join(out))
