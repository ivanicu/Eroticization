import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A158 R453 -- `S` 是不是两个方向的混合

`#406a`:`S` 在 29 个结局里 8 个反号。`#408a`:那 8 个的共同点**不是**卡片的题型分类。
**⇒ 换个问法:也许问题不在「哪些结局」,而在 **`S` 本身不是一个方向**。**
`S` 是**块内所选项的平均稀有度**,而 `#401`/`#402` 已显示「稀有」与「常规也管用」**可以并存**。

两个活着的世界:
**A `S` 是两个方向的混合** -> 拆成**越轨块**与**普通块**两半后,**两半各自的符号一致性都上升**;
**B `S` 就是一个含混的量** -> 拆开后不升 -> **页面上每一句关于 `S` 的话都必须带结局限定**(`#406d` 已如此)。

⚠ 块的分法来自 `#327` 的**载荷自助**(`TRG` / `ORD`)—— **数据给的,不是我挑的**。

ESTIMAND        `S_TRG`(越轨块内平均稀有度)· `S_ORD`(普通块内);
                各自在 29 个结局上拟合 `结局 ~ S_half + c3⁻ + 类别数`;
                主量 = **两半各自的符号一致性**,对着合起来的 `S` 的 **0.7241**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;
                        **guard 26 显式传 `main_quantity='discrete_count'` 与 `sweep_detection`**(`#407c`)。
                【非零支】**两半都**越过 offset 零 -> 世界 A;
                【零支】未越 -> 世界 B。
⚠ 零的种类     `offset_control`:**拆开后一致性上升的零绝不是零** ——
                **任意**把块分成两组,两半各自都可能比合起来更一致(每半的题更少、更同质)。
                零 = **随机把块分成同样大小的两组**后重算两半一致性的分布。
IMPOSSIBLE      ① 两半的 n 比 `S` 小 -> 系数更不稳,**而随机零用同样的块数,已吸收这一点**;
                ② `TRG`/`ORD` 是在**同一份数据**上定的(`#327`)-> 这是**不同仪器**,不是独立样本;
                ③ 「两个方向」不等于「两个心理构念」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R401=(ROOT/'E01_sexual_as_a_value_not_a_category/A131_breadth_leakage/R401_which_coordinate/run.py').read_text()
exec(_R401.split('"""',2)[2].split('inv=pd.read_csv')[0])
ALLR=np.flatnonzero(ok); CO=coords(ALLR); C3=-CO[4]
BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
_inv=pd.read_csv('data/derived/inventory.csv')
_ons=[c for c in _inv[_inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
ncat=np.column_stack([np.isfinite(d[c].map(BINo).values.astype(float)) for c in _ons]).sum(1).astype(float)
print(f"⚠ 块的分法来自 `#327` 的载荷自助:**越轨块 {len(TRG)}** · **普通块 {len(ORD)}** · 共 {len(MB)} 块")
def S_over(blocks,minc=3):
    """在指定的块集合上算「所选项的平均稀有度」——与 `S` 同一构造,只是块集合变了。"""
    cv=np.zeros(NN); ps=np.zeros(NN)
    for bi in blocks:
        M_,ppl=MB[bi]
        rr=-np.log(np.clip(M_.mean(0),1e-4,1.)); nb=M_.sum(1)
        v=np.where(nb>0,(M_@rr)/np.maximum(nb,1),np.nan)
        g=np.isfinite(v); cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(cv>=minc,ps/np.maximum(cv,1),np.nan)
S_T=S_over(TRG); S_O=S_over(ORD)
print(f"⚠ `#392e`:`S_TRG` 有限 **{int(np.isfinite(S_T).sum()):,}** · `S_ORD` 有限 **{int(np.isfinite(S_O).sum()):,}** · "
      f"corr **{np.corrcoef(S_T[np.isfinite(S_T)&np.isfinite(S_O)],S_O[np.isfinite(S_T)&np.isfinite(S_O)])[0,1]:+.4f}**")
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def consistency(Sv):
    bs=[]
    for nm,y in OUT:
        y=np.asarray(y,dtype=float)
        g=ok&np.isfinite(Sv)&np.isfinite(C3)&np.isfinite(ncat)&np.isfinite(y)
        k=int(g.sum())
        if k<500: continue
        X=np.column_stack([np.ones(k),z(Sv,g),z(C3,g),z(ncat,g)])
        b,*_=np.linalg.lstsq(X,z(y,g),rcond=None); bs.append(float(b[1]))
    if not bs: return np.nan,0
    p=sum(1 for x in bs if x>0); return max(p,len(bs)-p)/len(bs),len(bs)
cT,kT=consistency(S_T); cO,kO=consistency(S_O)
S_ALL=0.7241                                   # `#406a`,从 results/ 读过的那个数
print(f"\n两半各自的符号一致性:")
print(f"   `S_TRG`(越轨块)**{cT:.4f}**({kT} 个结局)")
print(f"   `S_ORD`(普通块)**{cO:.4f}**({kO} 个结局)")
print(f"   合起来的 `S`(`#406a`)**{S_ALL:.4f}**")
NP_=400; rg=np.random.default_rng(13); nul=[]
allb=list(range(len(MB)))
for _ in range(NP_):
    p=rg.permutation(allb); a=sorted(p[:len(TRG)]); b=sorted(p[len(TRG):len(TRG)+len(ORD)])
    ca,_=consistency(S_over(a)); cb,_=consistency(S_over(b))
    if np.isfinite(ca) and np.isfinite(cb): nul.append(min(ca,cb))
nul=np.array(nul); HI=float(np.percentile(nul,95))
OBS=min(cT,cO)
print(f"\n⚠ offset 零(**随机把块分成同样大小的两组** {NP_} 次;"
      f"**任意分法都可能让两半各自更一致 -> 这个零绝不是零**):")
print(f"   两半中**较小**的一致性:**{nul.mean():.4f} ± {nul.std():.4f}** · 95 分位 **{HI:.4f}**")
print(f"   实测较小者 **{OBS:.4f}** -> **{(OBS-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈 -> 世界 A(`S` 是两个方向的混合)**' if OBS>HI else '**落在零里 -> 世界 B**'}")
negs=[]
rg2=np.random.default_rng(91)
for _ in range(200):
    p=rg2.permutation(allb); a=sorted(p[:len(TRG)]); b=sorted(p[len(TRG):len(TRG)+len(ORD)])
    ca,_=consistency(S_over(a)); cb,_=consistency(S_over(b))
    if np.isfinite(ca) and np.isfinite(cb): negs.append(min(ca,cb))
negs=np.array(negs); rate=float((negs>HI).mean())
print(f"\n负对照(**越阈率**,随机分块 {len(negs)} 次):**{100*rate:.1f}%**(合格 1–12%)")
print(f"\nguard 26 = **MDE 扫描**(⚠ 离散主量,**显式声明**),每级 20 次:")
det=[]
for gg in (0.2,0.4,0.6,0.8):
    hit=0
    for s_ in range(20):
        rg3=np.random.default_rng(30+int(gg*10)*107+s_)
        # 种一个「两半各自更一致」的世界:把两半的系数按比例推成同号
        bsT=[];bsO=[]
        for Sv,acc in ((S_T,bsT),(S_O,bsO)):
            for nm,y in OUT:
                y=np.asarray(y,dtype=float)
                g=ok&np.isfinite(Sv)&np.isfinite(C3)&np.isfinite(ncat)&np.isfinite(y)
                if int(g.sum())<500: continue
                X=np.column_stack([np.ones(int(g.sum())),z(Sv,g),z(C3,g),z(ncat,g)])
                b,*_=np.linalg.lstsq(X,z(y,g),rcond=None); acc.append(float(b[1]))
        def push(bs):
            bs=np.array(bs); f=rg3.random(len(bs))<gg
            bs[f]=np.abs(bs[f]); p=int((bs>0).sum()); return max(p,len(bs)-p)/len(bs)
        if min(push(bsT),push(bsO))>HI: hit+=1
    det.append(hit/20); print(f"   两半各按 **{gg:.0%}** 推成同号 -> 检出 **{hit}/20 = {hit*5:>3d}%**")
MDE=next((g for g,dv in zip((0.2,0.4,0.6,0.8),det) if dv>=0.8),1.0)
print(f"   **MDE = {MDE:.0%}** · 有意义 **50%**")
pd.DataFrame([dict(v_half='S_TRG',v_cons=cT,v_k=kT),dict(v_half='S_ORD',v_cons=cO,v_k=kO),
              dict(v_half='S_all',v_cons=S_ALL,v_k=29),
              dict(v_half='_null_hi',v_cons=HI,v_k=len(nul))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'split_S.csv',index=False)
NONNULL=OBS>HI
g=Gate('`S` 是不是两个方向的混合')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26(显式声明离散主量 + 扫描)',MDE,0.50,True,what='MDE 扫描',
    branch='non_null' if NONNULL else 'null',main_quantity='discrete_count',sweep_detection=det)
g.asserted('★【两支】offset 零非退化(任意分法都可能让两半更一致)',nul.std()>0,
           f"{nul.mean():.4f} ± {nul.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】两半中**较小**的一致性越过 offset 零 -> 世界 A',NONNULL,
               f"{OBS:.4f} vs 上侧 {HI:.4f} · TRG {cT:.4f} · ORD {cO:.4f} · 合起来 {S_ALL:.4f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
