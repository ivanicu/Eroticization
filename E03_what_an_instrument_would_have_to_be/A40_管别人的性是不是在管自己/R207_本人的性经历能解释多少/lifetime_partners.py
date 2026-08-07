"""#768 第二臂 —— 主臂用的是「去年」,而那个民间理论说的是「一生」

第一臂:`partners`(去年伴侣数)· `sexfreq`(去年频率)· `evstray` 合起来只解释 **1.9%**。
⚠ **但这是真缺口**:对一个已婚五十岁的人,「去年伴侣数」恒等于 1,与他的**历史**无关。
   而「管别人的性 = 管自己做过的事」这个说法讲的是**一生**,不是去年。
⇒ 换成一生伴侣数 `numwomen + nummen`。

⚠ **码位陷阱,跑前摊开**:两列都带**哨兵码 989–997**('several' · 'many/lots' · 'n.a.' · refused),
   而 **99 分位就是 990** —— 直接当数用会把「refused」读成 997 个伴侣。⇒ 有效码取 **0–988**。

预注册(判词按 `#764` 新写法):
  一生伴侣数解释的份额,与**政治 [6.8%, 11.6%]** 和**宗教 [33.3%, 46.7%]** 并排比。
  若 ≥ 政治下界 6.8% ⇒ **主臂用错了仪器**,「去年」不是这个理论要的量,第一臂的否定要放宽;
  若 < 6.8% ⇒ 否定更强,而且是在这个理论**自己要求的**变量上。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(1207)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"homosex":(1,4),"teensex":(1,4),
       "numwomen":(0,988),"nummen":(0,988),"age":(18,89),"sex":(1,2)}
print("=== #766 前瞻使用 ===")
for c,rng in VALID.items():
    dr,tot=check_kept_codes(gp,c,rng)
    if dr: print(f"  {c:9s} keep={rng} -> "+" · ".join(f"删 码{int(a)} {b!r} {n}({sh*100:.2f}%)" for a,b,n,sh in dr[:3]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["lifetime"]=M.numwomen+M.nummen
sub0=M[["obey","lifetime","age"]].dropna()
print(f"\n一生伴侣数 n={len(sub0)} · 分位 "
      f"{[float(np.nanquantile(sub0.lifetime,q)) for q in (.25,.5,.75,.9,.99)]} · 最大 {float(sub0.lifetime.max())}")
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
SPECS={"只一生伴侣数":["lifetime"],"只年龄":["age"],"一生伴侣数 + 年龄":["lifetime","age"]}
print(f"\n=== G3 网格 · 保留率 ===")
print(f"  {'控制量':18s}"+"".join(f"{s[:9]:>11s}" for s in SEX)+f"{'中位':>9s}")
res={}; ns={}
for nm,cols in SPECS.items():
    row={}
    for s in SEX:
        sub=M[["obey","lifetime","age",s]].dropna()
        raw=prho(sub.obey.to_numpy(),sub[s].to_numpy())
        v=prho(sub.obey.to_numpy(),sub[s].to_numpy(),sub[cols].to_numpy())
        fl=3*1.65/np.sqrt(len(sub)); row[s]=dict(val=v,raw=raw,n=len(sub),keep=(v/raw if abs(raw)>=fl else None)); ns[s]=len(sub)
    res[nm]=row; ks=[row[s]["keep"] for s in SEX if row[s]["keep"] is not None]
    res[nm]["_med"]=float(np.median(ks))
    print(f"  {nm:18s}"+"".join(f"{row[s]['keep']*100:10.1f}%" for s in SEX)+f"{res[nm]['_med']*100:8.1f}%")
print(f"  {'n':18s}"+"".join(f"{ns[s]:11d}" for s in SEX))
sub=M[["obey","lifetime","age","premarsx"]].dropna()
raw=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy())
pc=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),RNG.normal(0,1e-9,(len(sub),1)))
nul=[prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),
          sub[["lifetime"]].to_numpy()[RNG.permutation(len(sub))]) for _ in range(300)]
G=Gate("#768 第二臂 · 一生伴侣数")
G.identity_control("① 常数控制须回到偏前", observed=pc, expected=raw, tol=0.005, what="仪器活着吗")
G.identity_control("② 打乱一生伴侣数须回到偏前", observed=float(np.median(nul)), expected=raw, tol=0.005,
                   what="打乱只毁掉「谁的性史配谁的态度」")
G.offset_control("③ 偏掉真的一生伴侣数后须显著低于偏前",
                 effect=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[["lifetime"]].to_numpy()),
                 offset=raw, spread=float(np.std(nul)),
                 null_kind="同一批人、同一对题,唯一差别是控制量是真的一生伴侣数还是打乱后的")
print(); print(G)
life=1-res["只一生伴侣数"]["_med"]; age=1-res["只年龄"]["_med"]; both=1-res["一生伴侣数 + 年龄"]["_med"]
print("\n"+"="*70)
if not all(r[2] for r in G.rows): print("**UNVERIFIED:闸没全过**")
else:
    print(f"  政治 [6.8%, 11.6%] · 宗教 [33.3%, 46.7%]")
    print(f"  **一生伴侣数 {life*100:.1f}%** · 只年龄 {age*100:.1f}% · 两者一起 {both*100:.1f}%")
    print(f"  ⇒ {'**≥ 政治下界 ⇒ 主臂用错了仪器**' if life>=0.068 else '**< 政治下界 ⇒ 否定更强,而且是在这个理论自己要求的变量上**'}")
json.dump(dict(n=ns,med={k:res[k]["_med"] for k in SPECS},life=life,age=age,both=both,
               gate_ok=all(r[2] for r in G.rows)),open(OUT/"lifetime.json","w"),ensure_ascii=False,indent=1)
