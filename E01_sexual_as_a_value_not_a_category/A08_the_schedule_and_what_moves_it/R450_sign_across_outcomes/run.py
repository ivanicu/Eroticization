import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A155 R450 -- `S` 在所有结局上朝同一个方向吗

`#405a` 的表里掉出一件事:**在「能不能改」上 `z(S)` 是负的(−0.093),
而在其它三个结局上都是正的。**
而**页面上关于 `S` 的每一句话都来自「羞耻」那一个结局** —— **没有一句检查过它在别处是否同号。**

⚠ 这**不是** `#396b` 那一类空分离器:**`S` 在多结局上的符号一致性从没被估过。**

ESTIMAND        项目自己的 **29 个结局**上逐一拟合 `结局 ~ z(S) + c3⁻ + 类别数`(同一控制集);
                主量 = **`S` 系数的符号一致性**(取多数符号的比例)。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**;offset 零非退化。
                【非零支】一致性**越过** offset 零的上侧 -> `S` 确实朝一个方向;
                【零支】落在零里 -> **`S` 不是一个有统一方向的量**,而页面上每一句关于它的话
                        都只在**它被发现的那个结局**上成立。
⚠ 零的种类     `offset_control`:**符号一致性的零绝不是零** ——
                29 个结局里,即使 `S` 与它们全都无关,多数符号也会占到一半以上。
                零 = **`lib.nulls.perm_in` 打乱 `S`** 后重算一致性的分布(保住结局之间的相关)。
⚠ 多重性       29 个结局 -> **报分布**,不报单格。
IMPOSSIBLE      ① 29 个结局彼此相关 -> 一致性会被**抬高**,**偏向「朝一个方向」** -> 方向上**不保守**;
                ② 符号一致不代表**幅度**一致;③ 结局的**编码方向**由问卷决定,不由我 ——
                   所以「一致性低」也可能只是**问卷的编码方向不统一**,这一条无法在本轮排除。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R401=(ROOT/'E01_sexual_as_a_value_not_a_category/A131_breadth_leakage/R401_which_coordinate/run.py').read_text()
exec(_R401.split('"""',2)[2].split('inv=pd.read_csv')[0])
# ⚠ R401 的 splice 来自 R333,提供的是 `coords` 与 R333 自己算好的 `S`,**没有** `fit_apply`。
# (第一版照抄了别轮的写法 —— **splice 提供什么,要看它自己,不要凭记忆。**)
ALLR=np.flatnonzero(ok); CO=coords(ALLR); C3=-CO[4]
# `ncat`(起始类别数)也不在这条 splice 里 —— 就地建,并按 `#392e` 先看它自己。
BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
_inv=pd.read_csv('data/derived/inventory.csv')
_ons=[c for c in _inv[_inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
ncat=np.column_stack([np.isfinite(d[c].map(BINo).values.astype(float)) for c in _ons]).sum(1).astype(float)
print(f"⚠ `#392e`:`ncat`(起始类别数)取值 0–{int(ncat.max())} · 众数 **{int(pd.Series(ncat).mode().iloc[0])}** · "
      f"与 `S` 相关 **{np.corrcoef(ncat[ok&np.isfinite(S)],S[ok&np.isfinite(S)])[0,1]:+.4f}**")
print(f"项目自己的结局集:**{len(OUT)}** 个")
M0=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(ncat)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def coefS(y,Sv=None):
    Sv=S if Sv is None else Sv
    g=M0&np.isfinite(y)&np.isfinite(Sv)
    k=int(g.sum())
    if k<500: return np.nan
    X=np.column_stack([np.ones(k),z(Sv,g),z(C3,g),z(ncat,g)]); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); return float(b[1])
rows=[]
for nm,y in OUT:
    b=coefS(np.asarray(y,dtype=float))
    if np.isfinite(b): rows.append(dict(v_out=str(nm)[:44],v_b=b))
T=pd.DataFrame(rows); check_columns(T,'R450')
T.to_csv(pathlib.Path(__file__).parent/'results'/'S_signs.csv',index=False)
K=len(T); npos=int((T.v_b>0).sum())
CONS=max(npos,K-npos)/K
print(f"可算的结局 **{K}** 个 · `S` 系数为正 **{npos}** · 为负 **{K-npos}**")
print(f"**符号一致性(多数符号的比例)= {CONS:.4f}**")
print(f"   最大 **{T.v_b.max():+.4f}**({T.loc[T.v_b.idxmax(),'v_out'][:36]})· "
      f"最小 **{T.v_b.min():+.4f}**({T.loc[T.v_b.idxmin(),'v_out'][:36]})")
NP_=400
nul=[]
for s_ in range(NP_):
    Sp=perm_in(S,M0,9300+s_)
    bs=[coefS(np.asarray(y,dtype=float),Sp) for nm,y in OUT]
    bs=[x for x in bs if np.isfinite(x)]
    if not bs: continue
    p=sum(1 for x in bs if x>0); nul.append(max(p,len(bs)-p)/len(bs))
nul=np.array(nul); HI=float(np.percentile(nul,95))
print(f"\n⚠ offset 零(**`lib.nulls.perm_in` 打乱 `S`**,保住结局之间的相关;"
      f"**29 个结局里多数符号本来就占一半以上 -> 这个零绝不是零**):")
print(f"   **{nul.mean():.4f} ± {nul.std():.4f}** · 95 分位(上侧)**{HI:.4f}**")
print(f"   实测 **{CONS:.4f}** -> **{(CONS-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈 -> `S` 确实朝一个方向**' if CONS>HI else '**落在零里 -> `S` 不是一个有统一方向的量**'}")
negs=[]
for s_ in range(200):
    Sp=perm_in(S,M0,99000+s_)
    bs=[coefS(np.asarray(y,dtype=float),Sp) for nm,y in OUT]
    bs=[x for x in bs if np.isfinite(x)]
    p=sum(1 for x in bs if x>0); negs.append(max(p,len(bs)-p)/len(bs))
negs=np.array(negs); rate=float((negs>HI).mean())
print(f"\n负对照(**越阈率**,打乱 `S` 200 次):**{100*rate:.1f}%**(合格 1–12%)")
MDE=None
print(f"\nguard 26 = **MDE 扫描**(种一个真实的统一方向),每级 20 次:")
for gg in (0.02,0.04,0.06,0.10):
    hit=0
    for s_ in range(20):
        rg=np.random.default_rng(120+int(gg*100)*101+s_)
        bs=[]
        for nm,y in OUT:
            yy=np.asarray(y,dtype=float).copy()
            g=M0&np.isfinite(yy)
            yy[g]=yy[g]+gg*np.nanstd(yy[g])*z(S,g)          # 给每个结局都加同向的 S
            bs.append(coefS(yy))
        bs=[x for x in bs if np.isfinite(x)]
        p=sum(1 for x in bs if x>0)
        if max(p,len(bs)-p)/len(bs)>HI: hit+=1
    print(f"   统一方向强度 **{gg:.2f}** -> 检出 **{hit}/20 = {hit*5:>3d}%**")
    if MDE is None and hit>=16: MDE=gg
MDE_=MDE if MDE else 0.15
NONNULL=CONS>HI
g=Gate('`S` 在所有结局上朝同一个方向吗')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 有意义的统一方向强度 0.05',MDE_,0.05,True,
    what='MDE 扫描 80% 检出',branch='non_null' if NONNULL else 'null')
g.asserted('★【两支】offset 零非退化(多数符号本来就占一半以上)',nul.std()>0,
           f"{nul.mean():.4f} ± {nul.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】一致性越过 offset 零的上侧 -> `S` 确实朝一个方向',NONNULL,
               f"{CONS:.4f} vs 上侧 {HI:.4f} · 为正 {npos}/{K}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **IMPOSSIBLE ①**:29 个结局彼此相关 -> 一致性被**抬高**,**偏向「朝一个方向」**"
      f" -> 方向上**不保守**。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
