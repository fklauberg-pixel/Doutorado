"""Simulação: com n = 24 parcelas, regressão linear x Random Forest (LOOCV).
Verdade: LAI = f(altura, cobertura) + ruído, com R² verdadeiro fixado.
Cenário A: relação linear. Cenário B: relação com saturação (exponencial)."""
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut, cross_val_predict
def r2(y,p): return 1-np.sum((y-p)**2)/np.sum((y-y.mean())**2)
def um(cen,R2true,rep):
    rng=np.random.default_rng(rep*31+int(R2true*10)+(cen=="linear")*1000)
    n=24; H=rng.normal(0.9,0.2,n); CC=np.clip(rng.normal(0.5,0.15,n),0.05,0.95)
    sig=0.6*H+1.5*CC if cen=="linear" else 3*(1-np.exp(-1.2*(H+CC)))
    y=sig+rng.normal(0,np.sqrt(np.var(sig)*(1-R2true)/R2true),n); X=np.c_[H,CC]
    pl=cross_val_predict(LinearRegression(),X,y,cv=LeaveOneOut())
    pr=cross_val_predict(RandomForestRegressor(100,min_samples_leaf=2,random_state=rep,n_jobs=1),X,y,cv=LeaveOneOut())
    return dict(cenario=cen,R2_verdadeiro=R2true,R2_linear=r2(y,pl),R2_RF=r2(y,pr))
jobs=[(c,r,k) for c in ["linear","saturacao"] for r in [0.3,0.5,0.7] for k in range(40)]
d=pd.DataFrame(Parallel(n_jobs=4)(delayed(um)(*j) for j in jobs))
s=d.groupby(["cenario","R2_verdadeiro"]).agg(R2_linear=("R2_linear","mean"),R2_RF=("R2_RF","mean"))
s["pct_RF_melhor"]=d.assign(w=d.R2_RF>d.R2_linear).groupby(["cenario","R2_verdadeiro"]).w.mean()*100
print(s.round(3)); s.round(3).to_csv("simulacao_linear_vs_RF.csv")
for cvp in [0.10,0.15,0.20]:
    se=cvp/np.sqrt(5); tot=0.115; rel=max(0,(tot**2-se**2)/tot**2)
    print(f"CV entre plantas {cvp:.0%}: confiabilidade do LAI de campo = {rel:.2f}")
