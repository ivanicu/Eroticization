"""#753 第二臂 —— 人层那个「族」,用切社会层的同一把刀切一次

第一臂在**社会**层量到:没有「服从那一族」,只有服从这一个(责任只值 +0.0046)。
`#752` 在**人**层说性属于**权威-内群体这一族**,证据是两个 ρ 都高(+0.5001 / +0.4194)。
⚠ **但「两个都高」可能是一件事,不是两件** —— 若权威与内群体彼此高度相关。
本臂做**与社会层完全同一个操作**:把族里的成员互相偏掉,看谁还站得住。

预注册(两支我都不预先偏好,写清楚它们各自意味着什么):
  W-1 **族是假的,两层其实一致**:偏掉权威后 ρ(chastity, 内群体) 掉到 ≤ 一半,
      且偏掉内群体后 ρ(chastity, 权威) 保留 ≥ 75% ⇒ 人层也是「权威这一个」,
      ⇒ **两个单位说的是同一件事,而 `#752` 的族读法两层都不成立。**
  W-2 **族是真的,两层确实不同**:两个偏相关都保留 ≥ 75%
      ⇒ 人层有两个独立成分,社会层只有一个 —— 单位本身改变了结构。
  W-3 之间 ⇒ 判不了
⚠ MFQ n=7267,偏相关的零上界约 0.02 ⇒ **这不会是 n 的问题,会是构念的问题。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp
RNG=np.random.default_rng(1195)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
SEX="chastity"; F=["HARM_AVG","FAIRNESS_AVG","INGROUP_AVG","AUTHORITY_AVG"]
D=d[[SEX]+F].apply(pd.to_numeric,errors="coerce").dropna()
def rho(a,b): return sp(pd.Series(a),pd.Series(b))
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
print(f"n={len(D)}")
print(f"\n=== 先问最便宜的那个问题:权威与内群体彼此相关多少? ===")
rAI=rho(D.AUTHORITY_AVG,D.INGROUP_AVG)
print(f"  ρ(权威, 内群体) = {rAI:+.4f}   ⚠ 若很高,「两个都高」就是一件事")
for a,b in [("HARM_AVG","FAIRNESS_AVG"),("AUTHORITY_AVG","HARM_AVG"),("INGROUP_AVG","FAIRNESS_AVG")]:
    print(f"  ρ({a:13s}, {b:13s}) = {rho(D[a],D[b]):+.4f}")

bA=rho(D[SEX],D.AUTHORITY_AVG); bI=rho(D[SEX],D.INGROUP_AVG)
pA=prho(D[SEX].to_numpy(),D.AUTHORITY_AVG.to_numpy(),D.INGROUP_AVG.to_numpy())
pI=prho(D[SEX].to_numpy(),D.INGROUP_AVG.to_numpy(),D.AUTHORITY_AVG.to_numpy())
print(f"\n=== 把族里的两个成员互相偏掉(与社会层同一把刀)===")
print(f"  {'':34s} {'偏前':>8s} {'偏后':>8s} {'保留':>7s}")
print(f"  ρ(chastity, 权威   | 偏掉内群体) {bA:+8.4f} {pA:+8.4f} {pA/bA*100:6.1f}%")
print(f"  ρ(chastity, 内群体 | 偏掉权威  ) {bI:+8.4f} {pI:+8.4f} {pI/bI*100:6.1f}%")
NP=2000
nA=[abs(prho(D[SEX].to_numpy(),RNG.permutation(D.AUTHORITY_AVG.to_numpy()),D.INGROUP_AVG.to_numpy())) for _ in range(NP)]
nI=[abs(prho(D[SEX].to_numpy(),RNG.permutation(D.INGROUP_AVG.to_numpy()),D.AUTHORITY_AVG.to_numpy())) for _ in range(NP)]
qA,qI=float(np.quantile(nA,.95)),float(np.quantile(nI,.95))
print(f"  各自的**零的 95% 分位**:权威 {qA:.4f} · 内群体 {qI:.4f}")
gz=RNG.normal(0,1e-9,(len(D),1))
pc=prho(D[SEX].to_numpy(),D.AUTHORITY_AVG.to_numpy(),gz)
print(f"  正控:控制量为常数须退回 {bA:+.4f},实测 {pc:+.4f} ⇒ {'通过' if abs(pc-bA)<0.01 else '不通过'}")
print("\n"+"="*66)
kA,kI=pA/bA,pI/bI
if abs(pc-bA)>=0.01: v="**UNVERIFIED:正控没过**"
elif kI<=0.5 and kA>=0.75:
    v=(f"**W-1:偏掉权威后内群体只剩 {kI*100:.1f}%,而偏掉内群体后权威保留 {kA*100:.1f}% "
       f"⇒ 人层也是「权威这一个」,不是一个族 —— **两个单位说的是同一件事,`#752` 的族读法两层都不成立**")
elif kA>=0.75 and kI>=0.75:
    v=(f"**W-2:两个都保留(权威 {kA*100:.1f}% · 内群体 {kI*100:.1f}%)⇒ 人层真有两个独立成分,"
       f"而社会层只有一个 —— **单位本身改变了结构**")
else: v=f"**W-3:权威保留 {kA*100:.1f}% · 内群体保留 {kI*100:.1f}% —— 落在两条预注册之间,判不了**"
print(v)
json.dump(dict(n=len(D),rho_auth_ingroup=rAI,biv_auth=bA,biv_ingroup=bI,
               partial_auth=pA,partial_ingroup=pI,null95_auth=qA,null95_ingroup=qI,
               pos_control=pc,verdict=v),open(OUT/"family_is_it.json","w"),ensure_ascii=False,indent=1)
