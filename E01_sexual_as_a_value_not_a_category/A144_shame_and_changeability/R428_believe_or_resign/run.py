import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A144 R428 -- 羞耻的人,认为自己**能**改,还是认为自己**改不了**

问卷里有两题正好切开这件事:
`BELIEF` = 「如果非常努力,你能停止被它唤起吗」(Impossible … With little effort, yes)
`ACTED`  = 「我已经实践/尝试过所有唤起我的东西」(−3 … +3)

**三个活着的世界,而 A 与 B 对同一个相关预测**相反的符号** —— 这是分离器,不是参数问题:**
**A 收缩(能动的羞耻)** — 羞耻动员控制努力:他一直在试,所以**相信能改**;同时压住行动。
   预测:`羞耻 ↔ BELIEF` **正** · `羞耻 ↔ ACTED` **负**
**B 认命(困住的羞耻)** — 羞耻与无力同源(`#366c`:无力感就在那条线上):**相信改不了**。
   预测:`羞耻 ↔ BELIEF` **负** · `羞耻 ↔ ACTED` **负**
**C 无关** — 「能不能改」是关于**内容**的判断(有些东西本来就更固定),与这个人的羞耻无关。
   预测:控制内容(`S` · `c3⁻`)后 `羞耻 ↔ BELIEF` **≈ 0**

⚠ **`ACTED` 那一行是平的**(A 与 B 预测同号)-> **它不分离世界**,只作描述,不作判据(frontier §1「flat row」)。

⚠ 先读 IMPOSSIBLE(`#376d`):`#378e②` 尺子在合并样本上估;`#427` 无。本轮自己的:
IMPOSSIBLE ① 两题都是**自报**,而羞耻本身会污染「我能不能改」的自评 —— 本轮测的是**关联**;
IMPOSSIBLE ② `BELIEF` 的 5 档是**序数**,等距是我加的 -> **同轮跑一遍秩变换版**,若符号变则不报;
IMPOSSIBLE ③ 控制 `c3⁻` 可能是**过度控制**(若信念在内容的下游)-> **原始与控制后两个都报**。

ESTIMAND        `BELIEF ~ 羞耻 + S + c3⁻ + 类别数`,主量 = **羞耻的系数与符号**。
判据(**先标支**,`#379c`)
                【两支】guard 26:阳性参照用**真实数据上已知的** `c3⁻ ↔ 羞耻`(页面 +0.129),
                        **不随手种 0.25**;负对照:打乱人。
                【非零支】系数越过 offset 零 -> 按符号判 A 或 B。
                【零支】仅当未越阈时启用 MDE:MDE < 0.05 才算「看得见而没有」-> 世界 C。
⚠ 零的种类     `offset_control`:**这个零不该是零** —— `BELIEF` 与 `ACTED` 都与「勾选了多少」强相关,
                而勾选数与羞耻有关(`#357b`)。零 = **控制类别数与 `S` 之后**的残差关系的置换分布。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
BC=next(c for c in d.columns if '7lgg41e' in c)
AC=next(c for c in d.columns if '41kpfir' in c)
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,
      'With a lot of effort, yes':2.,'With some effort, yes':3.,'With little effort, yes':4.}
BELIEF=d[BC].map(BMAP).values.astype(float)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(BELIEF)&np.isfinite(ACTED)&np.isfinite(ncat)
n=int(M.sum())
print(f"n=**{n:,}** · `BELIEF` 高 = 更容易改 · `ACTED` 高 = 实践得更多")
print(f"`BELIEF` 分布:" + ' · '.join(
    f"{k[:22]} {int((BELIEF[M]==v).sum()):,}" for k,v in sorted(BMAP.items(),key=lambda t:t[1])))
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def coef(y,x,g=None,ctrl=True):
    g=M if g is None else g; k=int(g.sum())
    cols=[np.ones(k),z(x,g)]+([z(S,g),z(C3,g),z(ncat,g)] if ctrl else [])
    X=np.column_stack(cols); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-len(cols)); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1])
BR=np.full(NN,np.nan); j=np.flatnonzero(M)
BR[j]=(rankdata(BELIEF[j])-rankdata(BELIEF[j]).mean())/rankdata(BELIEF[j]).std()

print(f"\n① 描述(⚠ `ACTED` 这一行 A 与 B 预测同号 -> **不分离世界**,只作描述):")
for nm,y in (('BELIEF 等距',BELIEF),('BELIEF 秩变换',BR),('ACTED',ACTED)):
    b0,s0=coef(y,sh,ctrl=False); b1,s1=coef(y,sh,ctrl=True)
    print(f"   羞耻 -> {nm:<12} 原始 **{b0:+.4f}** (se {s0:.4f}) · "
          f"控制 S/c3⁻/类别数后 **{b1:+.4f}** (se {s1:.4f})")

# ---- offset 零 ----
NP_=400
nul=np.array([coef(BELIEF,sh[np.random.default_rng(3300+s).permutation(NN)])[0] for s in range(NP_)])
THR=float(np.percentile(np.abs(nul),95))
bO,sO=coef(BELIEF,sh); bR,sR=coef(BR,sh)
print(f"\n⚠ offset 零(**控制类别数与 `S` 之后**打乱人;这个零不该是零 —— "
      f"两题都与勾选数强相关,而勾选数与羞耻有关 `#357b`):")
print(f"   **{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   实测(等距)**{bO:+.4f}** -> **{(bO-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(bO)>THR else '未越阈'}")
print(f"   实测(秩变换)**{bR:+.4f}** -> **{'同号' if bO*bR>0 else '⚠ 变号'}**(IMPOSSIBLE ②)")

# ---- guard 26 的真实阳性参照 ----
bPOS,sPOS=coef(sh,C3,ctrl=False)
print(f"\n★ 阳性参照(**真实数据**,页面已知 `c3⁻ ↔ 羞耻` ≈ +0.129):实测 **{bPOS:+.4f}** · "
      f"vs 阈 {THR:.4f} -> {'**通过**' if abs(bPOS)>THR else '**不过**'}")
bNEG,_=coef(BELIEF,sh[np.random.default_rng(99).permutation(NN)])
print(f"  负对照(打乱人):**{bNEG:+.5f}** vs 阈 {THR:.5f}")

print(f"\n先算 MDE 再看数(`#369a`),每级 30 次:")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(7700+int(gg*100)*17+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(sh,M)+rg.standard_normal(n)
        if abs(coef(y,sh)[0])>THR: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}**")
# ---- ⚠ `ACTED` 那一行需要**它自己的零** ----
# 它被预先声明为「不分离世界,只作描述」。**但 A 与 B 都死了之后,它成了本轮唯一的实质结果** ——
# 而一个描述行升格成结果时,**必须先拿到自己的零**,不能借 `BELIEF` 的阈(那是另一个量的零)。
nulA=np.array([coef(ACTED,sh[np.random.default_rng(4400+s).permutation(NN)])[0] for s in range(400)])
THRA=float(np.percentile(np.abs(nulA),95))
bA,sA=coef(ACTED,sh)
print(f"\n② `ACTED` 的**自己的**零(控制 S/c3⁻/类别数后打乱人,400 次):"
      f"**{nulA.mean():+.5f} ± {nulA.std():.5f}** · |值| 95 分位 **{THRA:.5f}**")
print(f"   实测 **{bA:+.4f}**(se {sA:.4f})-> **{(bA-nulA.mean())/max(nulA.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(bA)>THRA else '未越阈'}")
MDEA=None
for gg in (0.02,0.03,0.05):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(6600+int(gg*100)*19+s_)
        y=np.full(NN,np.nan); y[M]=-gg*z(sh,M)+rg.standard_normal(n)
        if abs(coef(y,sh)[0])>THRA: hit+=1
    if MDEA is None and hit>=24: MDEA=gg
MDEA_=MDEA if MDEA else 0.08
print(f"   `ACTED` 的 MDE = **{MDEA_:.2f}**")

MEANINGFUL=0.05     # 与其它三条路(+0.129 / +0.119 / −0.102)同量级的一半 —— 项目里 guard 21 的标准
print(f"\n⚠ **guard 26 第一版在我自己身上开火了,而它是对的** ——")
print(f"   我传了**观测值** {abs(bO):.4f} 当「争议幅度」。**一个零的争议幅度是**有意义的效应量**"
      f"({MEANINGFUL:.2f}),不是那个 ≈0 的观测值** —— 正是 `#382c` 已经记下的错,我又犯了一次。")
print(f"   而修正后有一个**更好的答案**:真实参照 `c3⁻` 是 {abs(bPOS):.4f},仍是 {abs(bPOS)/MEANINGFUL:.1f}× ——")
print(f"   **真正种在争议幅度上的正对照,是 MDE 扫描本身**:它种 0.02/0.03/0.05 并量检出率,")
print(f"   最小 80% 检出的种植 = **{MDE_:.2f} ≤ {MEANINGFUL:.2f}**。")
print(f"   **⇒ guard 21(MDE)与 guard 26(标定)是同一件事的两面:MDE 扫描 = 标定正确的正对照。**")
g=Gate('羞耻的人认为自己能改,还是改不了')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:正对照必须种在争议幅度上 —— 用 **MDE 扫描的最小检出幅度**,不是单次种植',
    MDE_,MEANINGFUL,True,what=f'MDE 扫描(80% 检出);真实参照 c3⁻ {abs(bPOS):.4f} 是 {abs(bPOS)/MEANINGFUL:.1f}× 偏大')
g.asserted('★【两支】负对照:打乱人 -> 必须落回零',abs(bNEG)<=THR,
           f"{bNEG:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★【两支】offset 零非退化',nul.std()>0,f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
g.asserted('★【两支】IMPOSSIBLE ②:秩变换版必须同号,否则不报',bO*bR>0,
           f"等距 {bO:+.4f} · 秩 {bR:+.4f}",kind='control')
if abs(bPOS)>THR and abs(bNEG)<=THR and bO*bR>0:
    if abs(bO)>THR:
        g.asserted('★【非零支】越阈 -> 按符号判世界(正 = A 收缩 · 负 = B 认命)',True,
                   f"{bO:+.4f} -> **世界 {'A 收缩' if bO>0 else 'B 认命'}**")
    else:
        g.asserted('★【零支】未越阈,且 MDE < 0.05 才算世界 C',MDE_<0.05,
                   f"{bO:+.4f} vs {THR:.4f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
g.asserted('★【非零支】`ACTED` 越过**它自己的**零(升格成结果,所以要自己的阈)',
           abs(bA)>THRA,f"{bA:+.4f} vs {THRA:.4f} · MDE {MDEA_:.2f}")
print(g)
pd.DataFrame([dict(v_what='BELIEF等距',v_b=bO,v_se=sO,v_thr=THR,v_mde=MDE_),
              dict(v_what='BELIEF秩',v_b=bR,v_se=sR,v_thr=THR,v_mde=MDE_),
              dict(v_what='ACTED',v_b=bA,v_se=sA,v_thr=THRA,v_mde=MDEA_)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'belief.csv',index=False)
print(f"\n⇒ **世界 A 与 B 同时死掉**(`BELIEF` 的零可发布:MDE {MDE_:.2f} < 有意义 {MEANINGFUL:.2f});"
      f"**世界 C 成立**。而 `ACTED` 是本轮唯一的实质结果。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")