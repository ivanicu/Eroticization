"""#751 第二臂 · 同一个问题,换一具人层仪器:MFQ(Graham/Haidt/Nosek 2009 Study 3)

⚠ **为什么必须有这一臂,以及闸是怎么逼出来的。**
GSS 那一臂给出「人层上性也自成一条线」,但它的控制量只有**两条题、且只覆盖「对国家不诚实」**。
`readme_gate` 的 `single_instrument`(`#658`)因此阻断了那一轮 —— **而它指出的正是我自己写进 NEXT ① 的弱点。**
正确的修法不是往文件里塞第二个仪器名(那正是该规则代理账里写明的失败模式),
是**把同一个问题放到第二具人层仪器上**。

⚠⚠ **而这一臂有真实的机会推翻第一臂,这是选它的理由:**
MFQ 的纯洁基础与权威/内群体是文献里出了名的高相关(binding foundations),
**所以 MFQ 很可能给出「性只是一般道德的一个侧面」** —— 与 GSS 那一臂相反。
按 realstat §2.5:**两具仪器不同号 ⇒ 框架本身就是发现,不许平均,去找它们分歧的那个假设。**

G1 估计量:
   S = `chastity`(MFQ 里唯一明确关于性的道德条目)
   G = 四个**非性**道德基础的均值(HARM · FAIRNESS · INGROUP · AUTHORITY)—— 比 GSS 的两条税务题厚得多
   P = 纯洁基础里**除 chastity 外**的其余条目(⚠ 必须剔除 chastity 本身,否则是部分-整体污染)
   ① 收敛:ρ(chastity, G)
   ② 判别:偏掉 G 之后,chastity 与 P 之间还剩多少
"""
import pandas as pd, numpy as np, json, pathlib, sys, itertools
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp
RNG=np.random.default_rng(1194)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)

PUR=[c for c in ["disgusting","decency","god","harmlessdg","chastity","unnatural"] if c in d.columns]
NONSEX=["HARM_AVG","FAIRNESS_AVG","INGROUP_AVG","AUTHORITY_AVG"]
print("=== 硬规则①:变量名不是测量 ===")
for c in PUR+NONSEX+["PURITY_AVG"]:
    if c in d.columns:
        v=pd.to_numeric(d[c],errors="coerce")
        print(f"  {c:14s} n={int(v.notna().sum()):6d} 取值 [{v.min():.2f}, {v.max():.2f}] 均值 {v.mean():.3f}")
    else: print(f"  {c:14s} 不存在")
SEXITEM="chastity"; P=[c for c in PUR if c!=SEXITEM]
print(f"\n  性条目 = {SEXITEM} · 其余纯洁条目 P = {P}")
print(f"  ⚠ chastity 已从 P 中剔除 —— 否则是部分-整体污染(`#749` 那一族)")

D=d[[SEXITEM]+P+NONSEX].apply(pd.to_numeric,errors="coerce").dropna()
print(f"\n完整个案 n={len(D)}(全样本 {len(d)})")
z=lambda s:(s-s.mean())/s.std(ddof=1)
D=D.copy(); D["G"]=z(D[NONSEX]).mean(axis=1)

def rho(a,b): return sp(pd.Series(a),pd.Series(b))
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
def alpha(df):
    k=df.shape[1]; return k/(k-1)*(1-df.var(ddof=1).sum()/df.sum(axis=1).var(ddof=1))
relG=alpha(z(D[NONSEX]))
print(f"\n=== 信度:G(四个非性基础)Cronbach α = {relG:.4f} ===")
print(f"  ⚠ 对比 GSS 那一臂的两条题 0.7031 —— **这一臂的控制量更厚,校正不足更少,对我的结论更不利**")

rSG=rho(D[SEXITEM],D.G)
print(f"\n=== ① 收敛:ρ({SEXITEM}, 四个非性道德基础) = {rSG:+.4f} ===")
print(f"  ⚠ GSS 那一臂的对应数是 +0.1317。**若这里显著更高,两具仪器不同意,框架就是发现。**")

before={p:rho(D[SEXITEM],D[p]) for p in P}
after ={p:prho(D[SEXITEM].to_numpy(),D[p].to_numpy(),D[NONSEX].to_numpy()) for p in P}
print(f"\n=== ② 判别:偏掉四个非性基础之后,{SEXITEM} 与其余纯洁条目还剩多少 ===")
print(f"  {'与':14s} {'偏前':>8s} {'偏后':>8s} {'保留':>7s}")
for p in P: print(f"  {p:14s} {before[p]:+8.4f} {after[p]:+8.4f} {after[p]/before[p]*100:6.1f}%")
mb=float(np.median(list(before.values()))); ma=float(np.median(list(after.values())))
print(f"  {'中位数':14s} {mb:+8.4f} {ma:+8.4f} {ma/mb*100:6.1f}%")

NP=2000
nul=[float(np.median([prho(D[SEXITEM].to_numpy(),D[p].to_numpy(),
       D[NONSEX].to_numpy()[RNG.permutation(len(D))]) for p in P])) for _ in range(NP)]
nq=(float(np.quantile(nul,.025)),float(np.quantile(nul,.975)))
print(f"\n=== 零:打乱「谁的非性道德配谁的纯洁态度」{NP} 次,重算偏相关 ===")
print(f"  零的 95% 区间 {nq[0]:+.4f} … {nq[1]:+.4f}")
gz=RNG.normal(0,1e-9,(len(D),1))
pc=float(np.median([prho(D[SEXITEM].to_numpy(),D[p].to_numpy(),gz) for p in P]))
print(f"=== 正控:控制量为常数时须退回偏前 {mb:+.4f},实测 {pc:+.4f} ⇒ {'通过' if abs(pc-mb)<0.01 else '不通过'} ===")

print("\n"+"="*64)
gss_conv=0.1317; gss_ret=0.972
if abs(pc-mb)>=0.01: v="**UNVERIFIED:正控没过**"
elif rSG>=2*gss_conv and ma/mb<0.75:
    v=(f"**W2 / 两具仪器不同意:MFQ 上 ρ(性, 一般道德) = {rSG:+.4f}(GSS 是 +0.1317),"
       f"偏掉后只保留 {ma/mb*100:.1f}%(GSS 是 97.2%)⇒ 框架就是发现,不许平均**")
elif ma/mb>=0.75:
    v=(f"**W1 / 两具仪器同向:MFQ 上偏掉四个非性基础后仍保留 {ma/mb*100:.1f}%(GSS 97.2%),"
       f"收敛 {rSG:+.4f} vs GSS +0.1317 ⇒ 人这个单位上性自成一条线,在两具仪器上都成立**")
else: v=f"**W3:收敛 {rSG:+.4f}、保留 {ma/mb*100:.1f}% —— 落在两条预注册之间,这个设计判不了**"
print(v)
json.dump(dict(n=len(D),sex_item=SEXITEM,purity_rest=P,rel_G=relG,rho_SG=rSG,
               before=before,after=after,med_before=mb,med_after=ma,null_ci=nq,
               pos_control=pc,gss_conv=gss_conv,gss_retention=gss_ret,verdict=v),
          open(OUT/"mfq_arm.json","w"),ensure_ascii=False,indent=1)
