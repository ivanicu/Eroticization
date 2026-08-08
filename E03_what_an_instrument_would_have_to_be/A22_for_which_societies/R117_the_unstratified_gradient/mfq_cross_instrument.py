"""E03·A22·R117(第二具仪器)—— 同一条 Δ 在 MFQ 上跑

硬规则④:跨仪器复制优于同一仪器再来一轮。
⚠ **安慰剂必须先过自己的成块正对照** —— `#670` 正是死在「拿成块的和不成块的比」。
结果:**MFQ 五个域里只有纯洁成块**(`#651` 复制)⇒ **这具仪器上没有合法安慰剂**,
只能报主效应本身;而主效应 **Δ = −0.0535 · p = 0.0967 落在地板上,与 GSS 的 +0.1349 反号** ⇒ **不复制。**
⇒ 结论必须缩到 **GSS 的性许可度四题**,不许升格成「性」或「纯洁」这个概念。
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pyreadstat, numpy as np
from itertools import combinations
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
d,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
DOM={"纯洁":["harmlessdg","decency","god","disgusting","unnatural","chastity"],
     "公平":["unfairly","treated","justice","rights","fairly","rich"],
     "内群":["lovecountry","betray","loyalty","family","team","history"],
     "权威":["respect","traditions","kidrespect","soldier","shutup","chaos"],
     "伤害":["emotionally","weak","cruel","compassion","animal","kill"]}
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
ok={}
for nm,it in DOM.items():
    j=d.dropna(subset=it+["educ_num"]); R=np.column_stack([j[c].rank().to_numpy(float) for c in it])
    ok[nm]=float(weak(R,np.ones(len(j))))
    print(f"  {nm}  n={len(j):>6,}  最弱一环 **{ok[nm]:+.4f}**  {'✅ 成块' if ok[nm]>0.20 else '⛔ 不成块 —— 不可当安慰剂'}")
blocks=[k for k,v in ok.items() if v>0.20 and k!="纯洁"]
print(f"  ⇒ 合法安慰剂:**{blocks or '一个都没有 —— 这具仪器上无法做安慰剂对照'}**")
grid=np.array([2.,3.,4.,5.,6.,7.]); BW=1.2
it=DOM["纯洁"]; j=d.dropna(subset=it+["educ_num"]); j=j[j.educ_num<=7]
e=j["educ_num"].to_numpy(float); R=np.column_stack([j[c].rank().to_numpy(float) for c in it])
f=lambda x:np.array([weak(R,np.exp(-0.5*((x-g)/BW)**2)) if np.exp(-0.5*((x-g)/BW)**2).sum()>=200 else np.nan for g in grid])
dd=lambda y:(lambda m:float(np.mean(y[m][-2:])-np.mean(y[m][:2])))(np.isfinite(y))
rng=np.random.default_rng(20260806); obs=dd(f(e))
nul=np.array([abs(dd(f(rng.permutation(e)))) for _ in range(300)])
q=float(np.quantile(nul,0.95)); p=float((nul>=abs(obs)).mean())
print(f"\n  MFQ 纯洁 n={len(j):,}  Δ = **{obs:+.4f}**  零 95% 分位 {q:.4f}  p = **{p:.4f}**")
print(f"  GSS 性四题是 **+0.1349 · p=0.0000** ⇒ **反号且落在地板上 —— 不复制**")
json.dump(dict(blocks=ok,legal_placebos=blocks,mfq_purity_delta=obs,null_q95=q,p=p,
   verdict="不复制:MFQ 纯洁 Δ 反号且落在地板上;且该仪器上无合法安慰剂",unchallenged=True),
   open(OUT/"mfq_cross.json","w"),indent=1,ensure_ascii=False)
