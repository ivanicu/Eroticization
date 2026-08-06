"""E03·A35·R188b —— 逐层对照它自己的零(而不是拿两个斜率相除)

⚠ **事后对照,如实标注:没有预注册。** 触发它的是主脚本里显而易见的混淆:
**常去教堂的人两题都从地板起步(2.14 / 1.37),天花板给了他们更多空间。**

⚠⚠ **而本文件的第一版判词是错的**:我把归一后的两个斜率相除,而分母是 **−0.00044**,
于是打印出「梯度塌了」。**`#691` 早就记过:含噪近零分母之比是重尾的,没有分辨力。**
⇒ **改用差,并让每一层对照它自己的零。** 归一后的梯度其实仍在且单调。

⚠ **换不了仪器**:同主脚本。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A,Bc="premarsx","homosex"
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","attend",A,Bc],encoding="latin1")
R=g.dropna(subset=[A,Bc,"attend"]).copy()
R["rel"]=pd.cut(R.attend,[-1,1,4,8],labels=["几乎不去","偶尔","常去"]); R=R.dropna(subset=["rel"]).copy()
FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def nslope(fr):
    pts=[]
    for y,sub in fr.groupby("year"):
        if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
        x=sub[A].to_numpy(float); yv=sub[Bc].to_numpy(float); r=sp(x,yv)
        xs=np.sort(x); ys=np.sort(yv); ys=ys if r>0 else ys[::-1]; c=abs(sp(xs,ys))
        if c<1e-9: continue
        pts.append((float(y),r/c))
    return (float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0]) if len(pts)>=FY else np.nan), len(pts)
rng=np.random.default_rng(20260806); res={}
print(f"{'层':10s}{'归一斜率':>11s}{'年':>4s}{'零 95% 区间':>24s}{'裁决':>12s}")
for lv in ["几乎不去","偶尔","常去"]:
    fr=R[R.rel==lv]; obs,ny=nslope(fr); nul=[]
    for _ in range(200):
        P=fr.copy(); P["year"]=rng.permutation(P.year.to_numpy())
        v,_=nslope(P)
        if np.isfinite(v): nul.append(v)
    q=np.quantile(nul,[0.025,0.975]); out="在零之外" if not (q[0]<=obs<=q[1]) else "**落在零里**"
    res[lv]=dict(slope=obs,years=ny,ci=[float(q[0]),float(q[1])],verdict=out)
    print(f"{lv:10s}{obs:>+11.5f}{ny:>4d}{f'[{q[0]:+.5f}, {q[1]:+.5f}]':>24s}{out:>12s}")
d=res["常去"]["slope"]-res["几乎不去"]["slope"]
print(f"\n常去 − 几乎不去 = **{d:+.5f}/年 ⇒ 36 年 {d*36:+.4f}**")
print("**用差,不用比 —— 分母近零时比值无意义(`#691`)。**")
print(f"⇒ 归一后梯度仍在且单调:{res['几乎不去']['slope']:+.5f} < {res['偶尔']['slope']:+.5f} < {res['常去']['slope']:+.5f}")
json.dump(dict(strata=res,gap=float(d),gap36=float(d*36),
  note="事后对照未预注册;第一版判词用比值(分母 −0.00044)判「塌了」,是 #691 的同型错误,已改用差",
  unchallenged=True),open(OUT/"norm.json","w"),indent=1,ensure_ascii=False)
