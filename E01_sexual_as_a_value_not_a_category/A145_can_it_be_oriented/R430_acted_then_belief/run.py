import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A145 R430 -- 「做过了所以不那么羞耻」这条路,数据支持吗

`#385` 拿掉了「起始时间造出来的」这个对手。**剩下的对手是 `#384f④` 的反向:
「做过了所以不那么羞耻」。**

这份数据能对它说一句话的地方只有一个,而且**它不需要羞耻这个变量本身**
(所以不会被羞耻的自报污染 —— `#384f①` 的那条 IMPOSSIBLE 在这里不适用):

**若「做过 -> 羞耻降低」成立**,那条路的机制是**试过了、发现它还在** ——
所以做得多的人应当更倾向于认为**「这东西改不掉」**。预测:`ACTED → BELIEF` **负**。
**若「羞耻 -> 少做」成立**,`ACTED ↔ BELIEF` **没有这个方向上的理由**。预测:**≈ 0 或正**。

⚠ 而 `#384a` 已经证明**羞耻本身与 `BELIEF` 无关**(+0.0135,零可发布)——
**所以这里的零是有内容的:它不是「什么都测不到」,而是「那条机制该留下的痕迹没有留下」。**

ESTIMAND        `BELIEF ~ ACTED + S + c3⁻ + 类别数 + 羞耻`,主量 = **`ACTED` 的系数与符号**。
判据(**先标支**,`#379c`)
                【两支】guard 26 用 **MDE 扫描**(`#384d`)· 负对照(**层内/掩码内**打乱,`#385c`)·
                        offset 零非退化 · 秩变换版同号(沿用 `#384` 的 IMPOSSIBLE ②)。
                【非零支】**负号且越阈** -> 支持「做过 -> 羞耻降低」;**正号且越阈** -> 反对它。
                【零支】未越阈时启用 MDE:MDE < 0.05 才算「该留的痕迹没留下」。
⚠ 零的种类     `offset_control`:**这个零不该是零** —— 两题都与「勾选了多少」强相关(`#357b`),
                零 = **控制 `S`/`c3⁻`/类别数/羞耻之后**、在掩码内打乱 `ACTED` 的分布。
IMPOSSIBLE      ① 两题仍是自报,只是**羞耻不再是自变量**,所以羞耻的自报偏差不直接进入主量;
                ② `BELIEF` 的等距是我加的 -> 秩变换版必须同号;
                ③ **这仍不是因果**:`ACTED` 与 `BELIEF` 同时测量,谁在前无法由本设计定。
                   本轮能说的只有「**那条机制预测的痕迹在不在**」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
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
BR=np.full(NN,np.nan); j=np.flatnonzero(M)
r_=rankdata(BELIEF[j]); BR[j]=(r_-r_.mean())/r_.std()
def coef(y,x,g=None):
    g=M if g is None else g; k=int(g.sum())
    X=np.column_stack([np.ones(k),z(x,g),z(S,g),z(C3,g),z(ncat,g),z(sh,g)])
    yy=z(y,g); b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-6); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1])
def perm_in(v,g,seed):
    o=v.copy(); jj=np.flatnonzero(g&np.isfinite(v))
    o[jj]=v[jj][np.random.default_rng(seed).permutation(len(jj))]; return o
print(f"n=**{n:,}** · `BELIEF` 高 = 更容易改 · `ACTED` 高 = 实践得更多")
print(f"⚠ 主量里**没有**羞耻作为自变量(它只是控制项)-> `#384f①` 的自报污染不直接进入。\n")
bO,sO=coef(BELIEF,ACTED); bRk,sRk=coef(BR,ACTED)
NP_=400
nul=np.array([coef(BELIEF,perm_in(ACTED,M,3900+s))[0] for s in range(NP_)])
THR=float(np.percentile(np.abs(nul),95))
print(f"⚠ offset 零(**控制 S/c3⁻/类别数/羞耻之后**,在掩码内打乱 `ACTED`;"
      f"这个零不该是零 —— 两题都与勾选数强相关 `#357b`):")
print(f"   **{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"\n实测 `ACTED → BELIEF`(等距)**{bO:+.4f}** (se {sO:.4f}) -> "
      f"**{(bO-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(bO)>THR else '**未越阈**'} · 符号 **{'负' if bO<0 else '正'}**")
print(f"   秩变换版 **{bRk:+.4f}** -> **{'同号' if bO*bRk>0 else '⚠ 变号'}**(IMPOSSIBLE ②)")
print(f"\nguard 26 的正对照 = **MDE 扫描**(`#384d`),每级 30 次:")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(4700+int(gg*100)*29+s_)
        y=np.full(NN,np.nan); y[M]=-gg*z(ACTED,M)+rg.standard_normal(n)
        if abs(coef(y,ACTED)[0])>THR: hit+=1
    print(f"   种植 **{-gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
MEANINGFUL=0.05
print(f"   **MDE = {MDE_:.2f}** · 争议幅度(有意义)= **{MEANINGFUL:.2f}**")
bNEG,_=coef(BELIEF,perm_in(ACTED,M,77))
pd.DataFrame([dict(v_what='BELIEF等距',v_b=bO,v_se=sO,v_thr=THR,v_mde=MDE_),
              dict(v_what='BELIEF秩',v_b=bRk,v_se=sRk,v_thr=THR,v_mde=MDE_)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'acted_belief.csv',index=False)

g=Gate('「做过了所以不那么羞耻」这条路,数据支持吗')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:正对照用 **MDE 扫描**,不是单次种植',MDE_,MEANINGFUL,True,
    what='MDE 扫描 80% 检出 vs 有意义效应量')
g.asserted('★【两支】负对照:掩码内打乱 `ACTED` -> 必须落回零',abs(bNEG)<=THR,
           f"{bNEG:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★【两支】offset 零非退化',nul.std()>0,f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
g.asserted('★【两支】IMPOSSIBLE ②:秩变换版必须同号',bO*bRk>0,
           f"等距 {bO:+.4f} · 秩 {bRk:+.4f}",kind='control')
if MDE_<=MEANINGFUL and abs(bNEG)<=THR and bO*bRk>0:
    if abs(bO)>THR:
        g.asserted('★【非零支】越阈 -> 按符号判(负 = 支持「做过->羞耻降低」· 正 = 反对它)',True,
                   f"{bO:+.4f} -> **{'支持' if bO<0 else '**反对**'}**")
    else:
        g.asserted('★【零支】未越阈,且 MDE < 有意义 -> 那条机制该留的痕迹没留下',MDE_<MEANINGFUL,
                   f"{bO:+.4f} vs {THR:.4f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **仍不是因果**:`ACTED` 与 `BELIEF` 同时测量,谁在前本设计定不了。"
      f"本轮只说「**那条机制预测的痕迹在不在**」。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
