"""#751 第三臂 —— 两具仪器不同意,分歧的假设是「一般道德」指的是哪一半

GSS 的 G = 逃税 · 骗福利(内容上靠近**伤害-公平**)      -> 偏掉后性保留 97.2%
MFQ 的 G = 伤害 · 公平 · 内群体 · 权威                  -> 偏掉后只保留 59.5%
⚠ 两个 G 不是同一个构念。MFQ 的里面装着**权威与内群体**,GSS 的没有。
⚠⚠ 而 `#750` 在**社会**层量到的亚军品格,恰好是**服从 +0.3030** —— 十个里唯一逼近性的那个。

**预注册的分歧假设**:MFQ 的 G 里若只留伤害+公平(GSS 那一半),保留率应跳回 GSS 那一侧;
若只留内群体+权威,保留率应更低。⇒ 「性是不是自成一条线」**取决于线的另一侧站的是谁**。

预注册判词:
  W-A 分歧由「权威/内群体」造成:偏掉 HARM+FAIRNESS 的保留率 ≥ 85%,且偏掉 INGROUP+AUTHORITY 的 ≤ 75%
  W-B 分歧不由此造成:两个半边的保留率相差 < 10 个百分点
  W-C 之间 ⇒ 判不了
⚠ 阈值写在跑之前。每个规格报它自己的零。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp
RNG=np.random.default_rng(2194)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
SEX="chastity"; P=["disgusting","decency","god","harmlessdg","unnatural"]
F=["HARM_AVG","FAIRNESS_AVG","INGROUP_AVG","AUTHORITY_AVG"]
D=d[[SEX]+P+F].apply(pd.to_numeric,errors="coerce").dropna()
def rho(a,b): return sp(pd.Series(a),pd.Series(b))
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    return float(np.corrcoef(resid(r(a),np.column_stack([r(C[:,j]) for j in range(C.shape[1])])),
                             resid(r(b),np.column_stack([r(C[:,j]) for j in range(C.shape[1])])))[0,1])
mb=float(np.median([rho(D[SEX],D[p]) for p in P]))
print(f"偏前中位数 {mb:+.4f} · n={len(D)}\n")
print(f"⚠ 收敛那一半,逐个基础(这是分歧假设的直接证据):")
for f in F: print(f"  ρ(chastity, {f:14s}) = {rho(D[SEX],D[f]):+.4f}")
SPECS={"全部四个(= MFQ 臂)":F,"只伤害+公平(≈ GSS 那一半)":F[:2],"只内群体+权威":F[2:],
       "只权威":["AUTHORITY_AVG"],"只内群体":["INGROUP_AVG"],"只伤害":["HARM_AVG"],"只公平":["FAIRNESS_AVG"]}
res={}
print(f"\n{'规格':28s} {'偏后中位':>9s} {'保留':>7s} {'该规格自己的零 95% 上界':>22s}")
for name,cols in SPECS.items():
    a=float(np.median([prho(D[SEX].to_numpy(),D[p].to_numpy(),D[cols].to_numpy()) for p in P]))
    nul=[float(np.median([prho(D[SEX].to_numpy(),D[p].to_numpy(),
          D[cols].to_numpy()[RNG.permutation(len(D))]) for p in P])) for _ in range(400)]
    q=float(np.quantile(nul,.975))
    res[name]=dict(after=a,retention=a/mb,null_q975=q,cols=cols)
    print(f"  {name:26s} {a:+9.4f} {a/mb*100:6.1f}% {q:+22.4f}")
rHF=res["只伤害+公平(≈ GSS 那一半)"]["retention"]; rIA=res["只内群体+权威"]["retention"]
print("\n"+"="*64)
if rHF>=0.85 and rIA<=0.75:
    v=(f"**W-A:分歧由「权威/内群体」造成 —— 偏掉伤害+公平后性保留 {rHF*100:.1f}%(≈ GSS 的 97.2%),"
       f"偏掉内群体+权威后只剩 {rIA*100:.1f}% ⇒ 「性是不是自成一条线」取决于线的另一侧站的是谁**")
elif abs(rHF-rIA)<0.10: v=f"**W-B:两个半边差 {abs(rHF-rIA)*100:.1f} 个百分点 ⇒ 分歧不由权威/内群体造成,另找假设**"
else: v=f"**W-C:伤害+公平 {rHF*100:.1f}% · 内群体+权威 {rIA*100:.1f}% —— 落在两条预注册之间,判不了**"
print(v)
json.dump(dict(n=len(D),med_before=mb,specs=res,verdict=v,
               conv={f:rho(D[SEX],D[f]) for f in F}),open(OUT/"which_half.json","w"),ensure_ascii=False,indent=1)
