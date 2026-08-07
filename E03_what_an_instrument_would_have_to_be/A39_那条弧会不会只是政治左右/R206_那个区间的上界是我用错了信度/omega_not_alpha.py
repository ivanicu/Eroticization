"""#767 · E03·A39·R206 —— 那个区间的上界,是我用错了信度估计量

`#765` 给出区间:宗教 **[33.3%, 47.1%]** · 残余 **[49.0%, 63.3%]** · 残余/宗教 **[1.04×, 1.90×]**,
上界来自 **Cronbach α = 0.4804** 的去衰减校正。
⚠ 而 α **只在 τ-等价(各题载荷相等)下无偏**;异质量表上 **α ≤ 真信度** ⇒ **除以 α 是除小了 ⇒ 校正过头**
⇒ **宗教被高估、残余被低估** ⇒ 47.1% 是宗教的**上界**,49.0% 是残余的**下界**。

G1 估计量:**McDonald's ω** —— 不假设 τ-等价。三题时单因子模型恰好可识别,载荷有闭式:
  λ₁² = r₁₂r₁₃/r₂₃ · λ₂² = r₁₂r₂₃/r₁₃ · λ₃² = r₁₃r₂₃/r₁₂ ,  ω = (Σλ)² / ((Σλ)² + Σ(1−λ²))
⚠ **闭式解可能越界**(某个 λ² > 1 或 < 0),那说明单因子模型不适配 —— **必须如实说,不许截断后照用。**

⚠ 政治只有两题(`polviews`·`partyid`),**两题无法把载荷与误差分开** ⇒ ω 不可识别,
   **只能用 Spearman-Brown** —— 按 realstat §2 记为 **N/A + 它需要什么 = 第三个指标**,不是「计划中」。

预注册(判词按 `#764` 新写法:只比已测量的量):
  若 **ω > α** ⇒ 校正应更小 ⇒ 宗教份额落在 47.1% 以下、残余升回 49.0% 以上,**区间收窄,方向已知**
  若 **ω < α** ⇒ 与理论相反 ⇒ **报出来并停**,不许当成更好的估计照用
  若闭式解越界 ⇒ **单因子不适配,ω 不可用**,区间维持 `#765` 的样子
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(206)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"teensex":(1,4),"homosex":(1,4),
       "polviews":(1,7),"partyid":(0,6),"reliten":(1,4),"fund":(1,3),"attend":(0,8)}
# ⚠ #766① 说那具工具还没有前瞻使用者 —— 本轮先过它一遍,而不是等下一个 attend 出现
print("=== #766 前瞻使用:我保留的范围排除了哪些带标签的档 ===")
for c,rng in VALID.items():
    dr,tot=check_kept_codes(gp,c,rng)
    if dr: print(f"  {c:9s} keep={rng} -> "+" · ".join(f"删 码{int(a)} {b!r} {n}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
print("  (未列出的列 = 没有排除任何带标签的档)")

d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
BASE=["obey","attend","reliten","fund","polviews","partyid"]
def frame(s):
    D=M[BASE+[s]].dropna().copy()
    D["REL"]=z(D[["attend","reliten","fund"]]).mean(axis=1)
    D["POL"]=z(D[["polviews","partyid"]]).mean(axis=1); return D
D=frame("premarsx")
sp=lambda a,b: float(pd.Series(a).corr(pd.Series(b),method="spearman"))
r12,r13,r23=sp(D.attend,D.reliten),sp(D.attend,D.fund),sp(D.reliten,D.fund)
print(f"\n=== 宗教三题的相关(ω 的原料)· n={len(D)} ===")
print(f"  r(attend,reliten)={r12:+.4f} · r(attend,fund)={r13:+.4f} · r(reliten,fund)={r23:+.4f}")
l2=np.array([r12*r13/r23, r12*r23/r13, r13*r23/r12])
inrange=bool((l2>0).all() and (l2<1).all())
print(f"  闭式载荷平方 λ² = {np.round(l2,4)} ⇒ {'在 (0,1) 内,单因子适配' if inrange else '**越界,单因子不适配**'}")
def alpha(df):
    k=df.shape[1]; return k/(k-1)*(1-df.var(ddof=1).sum()/df.sum(axis=1).var(ddof=1))
a_rel=alpha(z(D[["attend","reliten","fund"]]))
if inrange:
    lam=np.sqrt(l2); w_rel=float(lam.sum()**2/(lam.sum()**2+np.sum(1-l2)))
else: w_rel=None
r2=sp(D.polviews,D.partyid); sb_pol=2*r2/(1+r2)
print(f"\n=== 信度估计量 ===")
print(f"  宗教 Cronbach α = {a_rel:.4f}")
print(f"  宗教 McDonald ω = {w_rel:.4f}" if w_rel else "  宗教 McDonald ω = 不可用(单因子不适配)")
if w_rel: print(f"    ⇒ ω − α = {w_rel-a_rel:+.4f}  {'(理论预期 ω ≥ α)' if w_rel>=a_rel else '**⚠ 与理论相反**'}")
print(f"  政治 Spearman-Brown = {sb_pol:.4f}  ⚠ 两题 ⇒ **ω 不可识别**(N/A;它需要第三个指标)")

def partial_from_R(Rm,i,j,ctrl):
    idx=[i,j]+list(ctrl); S=Rm[np.ix_(idx,idx)]; P=np.linalg.pinv(S)
    return float(-P[0,1]/np.sqrt(P[0,0]*P[1,1]))
def shares(rel_r,rel_p):
    med=lambda k: float(np.median([o[k]/o["raw"] for o in per.values()]))
    per={}
    for s in SEX:
        Ds=frame(s); V=["obey",s,"REL","POL"]; X=Ds[V]
        Rm=np.array([[sp(X[a],X[b]) for b in V] for a in V]); rel=np.array([1,1,rel_r,rel_p],float)
        Rc=Rm.copy()
        for a in range(4):
            for b in range(4):
                if a!=b: Rc[a,b]=Rm[a,b]/np.sqrt(rel[a]*rel[b])
        per[s]=dict(raw=Rm[0,1],pol=partial_from_R(Rc,0,1,[3]),rel=partial_from_R(Rc,0,1,[2]),
                    both=partial_from_R(Rc,0,1,[2,3]),min_eig=float(np.linalg.eigvalsh(Rc).min()))
    return (1-med("pol"),1-med("rel"),med("both"),min(o["min_eig"] for o in per.values()))
P_raw=shares(1.0,1.0); P_a=shares(a_rel,sb_pol)
print(f"\n=== 份额(政治 · 宗教 · 残余;最小特征值)===")
print(f"  未校正            {P_raw[0]*100:5.1f}% {P_raw[1]*100:5.1f}% {P_raw[2]*100:5.1f}%   ({P_raw[3]:+.4f})")
print(f"  用 α 校正(#765)  {P_a[0]*100:5.1f}% {P_a[1]*100:5.1f}% {P_a[2]*100:5.1f}%   ({P_a[3]:+.4f})")
if w_rel:
    P_w=shares(w_rel,sb_pol)
    print(f"  用 ω 校正(本轮)  {P_w[0]*100:5.1f}% {P_w[1]*100:5.1f}% {P_w[2]*100:5.1f}%   ({P_w[3]:+.4f})")
G=Gate("#767 · ω 而不是 α")
G.identity_control("① 信度取 1 时校正必须是恒等变换", observed=shares(1.0,1.0)[2], expected=P_raw[2],
                   tol=1e-9, what="去衰减在 rel=1 时是恒等;不等就是公式错了")
if w_rel:
    G.identity_control("② ω 与 α 必须给出不同的校正(否则本轮无内容)",
                       observed=abs(P_w[2]-P_a[2]), expected=0.0, tol=1e-6,
                       what="⚠ 这一条**故意反向**:若差为 0 则 PASS,而 PASS 在这里意味着本轮什么也没证明")
print(); print(G)
print("\n"+"="*70)
if not inrange: print("**单因子不适配,ω 不可用 ⇒ 区间维持 `#765` 的样子**")
elif w_rel<a_rel: print(f"**ω({w_rel:.4f}) < α({a_rel:.4f}),与理论相反 ⇒ 报出来并停,不当成更好的估计**")
else:
    print(f"**ω = {w_rel:.4f} > α = {a_rel:.4f} ⇒ 校正更小,区间收窄:**")
    print(f"  宗教 [{P_raw[1]*100:.1f}%, {P_w[1]*100:.1f}%](原 [{P_raw[1]*100:.1f}%, {P_a[1]*100:.1f}%])")
    print(f"  残余 [{P_w[2]*100:.1f}%, {P_raw[2]*100:.1f}%](原 [{P_a[2]*100:.1f}%, {P_raw[2]*100:.1f}%])")
    print(f"  残余/宗教 [{P_w[2]/P_w[1]:.2f}×, {P_raw[2]/P_raw[1]:.2f}×](原 [{P_a[2]/P_a[1]:.2f}×, {P_raw[2]/P_raw[1]:.2f}×])")
json.dump(dict(n=len(D),alpha=a_rel,omega=w_rel,sb_pol=sb_pol,lam2=l2.tolist(),inrange=inrange,
               raw=P_raw,alpha_corr=P_a,omega_corr=(P_w if w_rel else None)),
          open(OUT/"omega.json","w"),ensure_ascii=False,indent=1)
