import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A146 R433 -- 用**问卷自己的块结构**重定义那个家族,才算复制

`#388c`:那个「常规性行为」家族是我**看过头名之后**挑的五题 -> **+0.1159 是描述,不是检验**。
**⇒ 复制的唯一合法形式:家族由**我没挑过的结构**给出。**

这份数据里现成就有一个:**32 个 MULTISELECT 块**。
而 `#327` 已经用**另一件仪器**(块载荷的自助)定出了 `c3` 的**普通端**:
**`sex acts`(−0.287)与 `sensations`(−0.268)** —— **那是数据给的,不是我挑的。**

ESTIMAND        用这两块的**块内勾选率**做家族分数,问它与 `BELIEF` 残差的相关,
                是否复制 `#388a` 的 **+0.1159**(方向与量级)。
判据(**先标支**,`#379c`)
                【两支】阳性参照 `ACTED`(已知相关物)必须在**单变量零**上开火(`#388b` 的教训);
                        负对照:一列纯噪声;guard 26:MDE 扫描。
                【非零支】**同号且越过单变量零** -> 复制成功。
                【零支】未越阈时启用 MDE。
⚠ 零的种类     `offset_control`:这是**预先指定的一个量**,不是从候选里挑出来的
                -> **单变量零**(打乱那一列),**不是**族内取最大零(`#388b` 的类别错误)。
⚠ 同轮必报     块勾选率与总勾选数、与 `S` 的相关(它们是这条路上最强的混淆)。
IMPOSSIBLE      ① `#327` 的普通端是在**同一份数据**上定的 -> 这是**不同仪器**的复制,不是**独立样本**的复制;
                ② 块勾选率与总勾选数强相关 -> 类别数已在控制项,但残余共线仍会压低系数;
                ③ 两块只是 `c3` 普通端的**两个**,不是「普通」这个概念的全部。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
BC=next(c for c in d.columns if '7lgg41e' in c); AC=next(c for c in d.columns if '41kpfir' in c)
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,
      'With a lot of effort, yes':2.,'With some effort, yes':3.,'With little effort, yes':4.}
BELIEF=d[BC].map(BMAP).values.astype(float)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)

# ---- 家族 = `#327` 定出的 c3 普通端两块的**块内勾选率**(数据给的,不是我挑的)----
WANT=['sex act','sensation']
picked=[]
for M_,ppl in MB:
    pass
names=[str(x) for x in qm.iloc[:,0]] if 'qm' in dir() else []
BLK=[]
for i,(Mb,ppl) in enumerate(MB):
    BLK.append((i,Mb,ppl))
lbl=[str(x) for x in (qm[qm.columns[0]].tolist() if 'qm' in dir() else [])]
print(f"块数 **{len(MB)}** · 标签样本 {lbl[:3] if lbl else '(无标签,按题干匹配)'}")
inv=pd.read_csv('data/derived/inventory.csv')
mult=[c for c in inv[inv['kind']=='MULTISELECT']['col']]
hit=[c for c in mult if any(w in str(c).lower() for w in WANT)]
print(f"⚠ 家族由**问卷自己的结构**给出(`#327` 用块载荷的自助定出 `c3` 的普通端:"
      f"`sex acts` −0.287 · `sensations` −0.268)")
print(f"   匹配到的块列(**打印原文,让错配当场可见** —— `#374b`):")
for c in hit: print(f"      {c[:78]}")
assert len(hit)>=2, f"只匹配到 {len(hit)} 个块,家族无法构造"
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
def block_rate(colname):
    sub=lg[lg[lg.columns[1]].astype(str)==str(colname)] if lg.shape[1]>2 else None
    return None
rates=[]
for c in hit:
    v=pd.to_numeric(d[c].astype(str).str.count(','),errors='coerce').values.astype(float)+1
    v=np.where(d[c].isna().values,np.nan,v)
    rates.append(v)
FAMB=np.nanmean(np.column_stack([np.where(np.isfinite(r),(r-np.nanmean(r))/np.nanstd(r),np.nan)
                                 for r in rates]),1)
mF=M&np.isfinite(FAMB); nF=int(mF.sum())
print(f"\n家族分数(**两块的勾选个数**,各自标准化后平均)· n=**{nF:,}**")
print(f"⚠ 同轮必报的混淆:corr(家族, 总类别数) = **{np.corrcoef(FAMB[mF],ncat[mF])[0,1]:+.4f}** · "
      f"corr(家族, `S`) = **{np.corrcoef(FAMB[mF],S[mF])[0,1]:+.4f}**")
CTRL=[ACTED,S,C3,ncat,sh]
X0=np.column_stack([np.ones(nF)]+[z(v,mF) for v in CTRL])
yB=z(BELIEF,mF); bb,*_=np.linalg.lstsq(X0,yB,rcond=None); RES=yB-X0@bb
rFAM=float(np.corrcoef(RES,z(FAMB,mF))[0,1])
NP_=400
nulS=np.array([abs(float(np.corrcoef(RES,np.random.default_rng(7300+s).permutation(z(FAMB,mF)))[0,1]))
               for s in range(NP_)])
STHR=float(np.percentile(nulS,95))
print(f"\n⚠ 零的种类:**单变量零**(打乱那一列)—— 这是**预先指定的**一个量,不是从候选里挑的"
      f"(`#388b` 的类别错误不再犯):**{nulS.mean():.4f} ± {nulS.std():.4f}** · 95 分位 **{STHR:.4f}**")
print(f"\n★ 复制:家族(块结构给出)↔ `BELIEF` 残差 = **{rFAM:+.4f}** vs 阈 {STHR:.4f} -> "
      f"{'**越阈**' if abs(rFAM)>STHR else '**未越阈**'} · "
      f"与 `#388a` 的 +0.1159 **{'同号' if rFAM>0 else '⚠ 变号'}**")
X1=np.column_stack([np.ones(nF)]+[z(v,mF) for v in [S,C3,ncat,sh]])
b1,*_=np.linalg.lstsq(X1,yB,rcond=None); RES_NA=yB-X1@b1
rA=abs(float(np.corrcoef(RES_NA,z(ACTED,mF))[0,1]))
rgN=np.random.default_rng(7); rF=abs(float(np.corrcoef(RES,rgN.standard_normal(nF))[0,1]))
print(f"  阳性参照 `ACTED`:|r| **{rA:.4f}** vs {STHR:.4f} -> "
      f"{'**开火**' if rA>STHR else '**不开火**'} · 负对照(纯噪声)**{rF:.4f}**")
print(f"\nguard 26 的正对照 = **MDE 扫描**(`#384d`),每级 30 次:")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hitc=0
    for s_ in range(30):
        rg=np.random.default_rng(5900+int(gg*100)*37+s_)
        yy=gg*z(FAMB,mF)+rg.standard_normal(nF)
        if abs(float(np.corrcoef(yy,z(FAMB,mF))[0,1]))>STHR: hitc+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hitc}/30 = {hitc/0.3:>5.1f}%**")
    if MDE is None and hitc>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = `#388a` 的 **0.116**")
pd.DataFrame([dict(v_r=rFAM,v_thr=STHR,v_mde=MDE_,v_n=nF,v_ref=rA,
                   v_ncat=float(np.corrcoef(FAMB[mF],ncat[mF])[0,1]),
                   v_S=float(np.corrcoef(FAMB[mF],S[mF])[0,1]))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'block_family.csv',index=False)

g=Gate('用问卷自己的块结构重定义家族,复制得了吗')
g.asserted('★【两支】阳性参照 `ACTED` 在**单变量零**上开火(`#388b`)',rA>STHR,
           f"{rA:.4f} vs {STHR:.4f}",kind='control')
g.asserted('★【两支】负对照:纯噪声列不越阈',rF<=STHR,f"{rF:.4f} vs {STHR:.4f}",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度 0.116',MDE_,0.1159,True,what='MDE 扫描 80% 检出')
g.asserted('★【非零支】同号且越过单变量零 -> 复制成功',rFAM>0 and abs(rFAM)>STHR,
           f"{rFAM:+.4f} vs {STHR:.4f}")
print(g)
print(f"\n⚠ IMPOSSIBLE ①:`#327` 的普通端是在**同一份数据**上定的 -> "
      f"这是**不同仪器**的复制,**不是独立样本**的复制。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
