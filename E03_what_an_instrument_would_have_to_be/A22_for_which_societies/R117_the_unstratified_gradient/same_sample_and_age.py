"""E03·A22·R117(同轮混淆控制 + G4 第二协变量)

**最强混淆(先写下来的那个):** 两组题来自**不同样本、不同年份** ——
性四题单独 n=15,000 / 21 个年份,性别角色单独 n=29,591 / 23 个年份,交集损失 22% 与 61%。
⇒ **限定到七题全答的同一批人**重跑。
**G4:** `#671` 两个协变量都测过(年龄与教育),只报一个就是选择性报告 ⇒ 两个都跑。
结果(n=11,576):**教育上性四题 +0.1349 ✅ / 性别角色 −0.0322 落地板 —— 4.19× 反号;
年龄上性四题 +0.0318 落地板 / 性别角色 +0.0906 ✅ —— 0.35× 同号。**
**一个「作答一致性特质」不可能产生这个交叉。**
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ","age"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
j=df.dropna(subset=SEX+FEM+["educ","age"])
print(f"七题全答 + educ + age:**n = {len(j):,}** · {int(j.year.min())}–{int(j.year.max())}")
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def cn(rx,ry,s):
    a=np.sort(rx);b=np.sort(ry); return wc(a,b if s>0 else b[::-1],np.ones_like(a))
def weak(R,W):
    v=[]
    for i,k in combinations(range(R.shape[1]),2):
        r=wc(R[:,i],R[:,k],W)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        c=cn(R[:,i],R[:,k],1 if r>0 else -1)
        if np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
    return min(v) if v else np.nan
RK={nm:np.column_stack([pd.Series(j[c]).rank().to_numpy(float) for c in it])
    for nm,it in [("性四题",SEX),("性别角色三题",FEM)]}
res={}
for cov,grid,bw in [("educ",np.arange(8,19.01,1.0),2.5),("age",np.arange(22,79.01,5.0),7.0)]:
    x=j[cov].to_numpy(float); rng=np.random.default_rng(20260806); res[cov]={}
    print(f"\n=== {cov} · 带宽 {bw} ===")
    for nm,R in RK.items():
        f=lambda e:np.array([weak(R,np.exp(-0.5*((e-g)/bw)**2))
                             if np.exp(-0.5*((e-g)/bw)**2).sum()>=200 else np.nan for g in grid])
        d=lambda y:(lambda m:float(np.mean(y[m][-3:])-np.mean(y[m][:3])))(np.isfinite(y))
        obs=d(f(x)); nul=np.array([abs(d(f(rng.permutation(x)))) for _ in range(300)])
        q=float(np.quantile(nul,0.95)); p=float((nul>=abs(obs)).mean()); res[cov][nm]=(obs,q,p)
        print(f"  {nm:8s} Δ = **{obs:+.4f}**  零 95% 分位 {q:.4f}  p = **{p:.4f}**  {'✅' if p<0.05 else '⛔ 落在地板上'}")
    m,pl=res[cov]["性四题"][0],res[cov]["性别角色三题"][0]
    print(f"  ⇒ 幅度比 **{abs(m)/abs(pl):.2f}×** · {'**反号**' if m*pl<0 else '同号'}")
json.dump({"n":int(len(j)),"res":{k:{a:list(b) for a,b in v.items()} for k,v in res.items()},
           "unchallenged":True},open(OUT/"same_sample_and_age.json","w"),indent=1,ensure_ascii=False)
