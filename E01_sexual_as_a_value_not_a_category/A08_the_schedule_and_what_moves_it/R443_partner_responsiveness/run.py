import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A151 R443 -- 「被伴侣带动」:两道题互为复制,还是同一道题的两次改写

`#398c` 撞见两道题,都与「治疗性」的残差正相关,且都不是已知的四条路:
`P1` = 「如果伴侣被某样东西唤起,我也能被它唤起」**+0.1080**
`P2` = 「我能被伴侣做我喜欢的事唤起,即使我知道他自己并不被唤起」**+0.0761**

⚠ **先量它们彼此的相关**(`#392e`,跑在建模之前):
**若 r 高(>0.7)-> 它们是同一道题的改写,「互为复制」不算复制,本轮到此为止;**
**若 r 中等 -> 各自存活才有信息(`#391` 的留一逻辑,换到两题上)。**

ESTIMAND        ① `corr(P1, P2)`;
                ② `治疗性 ~ P1 + 羞耻 + S + c3⁻ + 类别数 + ACTED`(**P2 完全不进模型**),再反过来;
                ③ 两者同时进模型时的偏系数。
判据(**先标支**,`#379c`)
                【两支】`corr(P1,P2) ≤ 0.7`(否则互为复制无效,不判)·
                        阳性参照(**羞耻自己**,在**未去掉它**的残差上,单变量零)·
                        负对照用**越阈率**(`#395b`)· guard 26 用 **MDE 扫描**。
                【非零支】**两道各自都同号且越过自己的单变量零** -> 一个构念的两个指标,复制成立;
                          **只有一道** -> 那不是构念,是那一道题。
⚠ 零的种类     `offset_control`:两者都与 `ACTED`/类别数相关 -> **这个零不该是零**;
                零 = 控制之后、**`lib.nulls.perm_in`** 在掩码内打乱结局的分布。
⚠ 读数         **从 `results/` 读,不从终端读**(`#398d`)。
IMPOSSIBLE      ① 两道题都自报且措辞相近 -> 同源方差会抬高它们的相关,**这会让判据偏向「不算复制」**,
                   即**保守**方向;② 扫描发现的题在**同一份数据**上复制 -> 这是**内部**一致性,不是外部复制。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
TC=next(c for c in d.columns if 'vmq8jqw' in c)
THER=pd.to_numeric(d[TC],errors='coerce').values.astype(float)
AC=next(c for c in d.columns if '41kpfir' in c)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
P1C=next(c for c in d.columns if 'If my partner is aroused by something' in str(c))
P2C=next(c for c in d.columns if 'even if I know they' in str(c))
P1=pd.to_numeric(d[P1C],errors='coerce').values.astype(float)
P2=pd.to_numeric(d[P2C],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(THER)&np.isfinite(ACTED)&np.isfinite(P1)&np.isfinite(P2)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
print("⚠ **`#392e`:两道题进模型前先各自看清楚**")
for nm,v,cc in (('P1',P1,P1C),('P2',P2,P2C)):
    u=np.unique(v[np.isfinite(v)]); ga=M&np.isfinite(anc)
    print(f"   {nm} `{str(cc)[:56]}`")
    print(f"      取值 {u.tolist()} · 众数 **{float(pd.Series(v[np.isfinite(v)]).mode().iloc[0]):g}** · "
          f"与锚 `Totalsexacts` 相关 **{np.corrcoef(v[ga],anc[ga])[0,1]:+.4f}**")
R12=float(np.corrcoef(P1[M],P2[M])[0,1])
print(f"\n① **`corr(P1, P2)` = {R12:+.4f}** -> "
      f"{'⚠ **>0.7:它们是同一道题的改写,互为复制不算复制**' if abs(R12)>0.7 else '**中等 —— 各自存活才有信息**'}")
CTRL=[sh,S,C3,ncat,ACTED]
def fit(y,xs):
    X=np.column_stack([np.ones(n)]+[z(v,M) for v in xs]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-len(xs)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:],se[1:]
b1,s1=fit(THER,[P1]+CTRL); b2,s2_=fit(THER,[P2]+CTRL); bJ,sJ=fit(THER,[P1,P2]+CTRL)
NP_=400
n1=np.array([fit(perm_in(THER,M,4100+s),[P1]+CTRL)[0][0] for s in range(NP_)])
n2=np.array([fit(perm_in(THER,M,4600+s),[P2]+CTRL)[0][0] for s in range(NP_)])
T1=float(np.percentile(np.abs(n1),95)); T2=float(np.percentile(np.abs(n2),95))
print(f"\n② 各自单独(**另一道完全不进模型**),控制 羞耻·`S`·`c3⁻`·类别数·`ACTED`:")
print(f"   P1 **{b1[0]:+.4f}** (se {s1[0]:.4f}) vs 单变量零 **{T1:.4f}** · {'**越阈**' if abs(b1[0])>T1 else '未越阈'}")
print(f"   P2 **{b2[0]:+.4f}** (se {s2_[0]:.4f}) vs 单变量零 **{T2:.4f}** · {'**越阈**' if abs(b2[0])>T2 else '未越阈'}")
print(f"\n③ 一起放:P1 **{bJ[0]:+.4f}** (se {sJ[0]:.4f}) · P2 **{bJ[1]:+.4f}** (se {sJ[1]:.4f})")
negs=np.array([fit(perm_in(THER,M,70000+s),[P1]+CTRL)[0][0] for s in range(200)])
rate=float((np.abs(negs)>T1).mean())
X1=np.column_stack([np.ones(n)]+[z(v,M) for v in [S,C3,ncat,ACTED]]); yT=z(THER,M)
bb,*_=np.linalg.lstsq(X1,yT,rcond=None); RES_NS=yT-X1@bb
nulS=np.array([abs(float(np.corrcoef(RES_NS,perm_in(sh,M,7700+s)[M])[0,1])) for s in range(400)])
STHR=float(np.percentile(nulS,95)); rS=abs(float(np.corrcoef(RES_NS,z(sh,M))[0,1]))
print(f"\n阳性参照(**羞耻自己**,未去掉它的残差):|r| **{rS:.4f}** vs {STHR:.4f} -> "
      f"{'**开火**' if rS>STHR else '不开火'} · 负对照**越阈率** **{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**,每级 30 次(种在 P1 上):")
MDE=None
for gg in (0.02,0.03,0.05,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(1200+int(gg*100)*67+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(P1,M)+rg.standard_normal(n)
        if abs(fit(y,[P1]+CTRL)[0][0])>T1: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = **{max(abs(b1[0]),abs(b2[0])):.4f}**")
pd.DataFrame([dict(v_item='P1',v_alone=b1[0],v_joint=bJ[0],v_thr=T1),
              dict(v_item='P2',v_alone=b2[0],v_joint=bJ[1],v_thr=T2),
              dict(v_item='corr_P1P2',v_alone=R12,v_joint=np.nan,v_thr=np.nan)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'partner.csv',index=False)
OK1=abs(b1[0])>T1; OK2=abs(b2[0])>T2; SAME=np.sign(b1[0])==np.sign(b2[0])
g=Gate('两道「被伴侣带动」的题,互为复制还是同一道题的两次改写')
g.asserted('★【两支】`corr(P1,P2)` ≤ 0.7(否则互为复制无效)',abs(R12)<=0.7,
           f"{R12:+.4f}",kind='control')
g.asserted('★【两支】阳性参照:羞耻自己必须开火',rS>STHR,f"{rS:.4f} vs {STHR:.4f}",kind='control')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,max(abs(float(b1[0])),abs(float(b2[0]))),True,
    what='MDE 扫描 80% 检出')
if abs(R12)<=0.7 and rS>STHR and 0.01<=rate<=0.12:
    g.asserted('★【非零支】两道各自同号且各自越阈 -> 一个构念的两个指标',SAME and OK1 and OK2,
               f"P1 {b1[0]:+.4f}/{T1:.4f} · P2 {b2[0]:+.4f}/{T2:.4f} · 同号 {SAME}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
