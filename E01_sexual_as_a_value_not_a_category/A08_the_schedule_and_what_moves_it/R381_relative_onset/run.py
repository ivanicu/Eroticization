import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A123 R381 -- `S` 与「来得早」重叠的那一段,是不是人群时间表

`#335b`:`S × EARLY` 的次可加不是天花板。`#333c` 指出候选机制:
**`#130` 的人群时间表**(内容类早、关系类晚 -> 早来的兴趣本来就更冷门,
`corr(起始均值, S) = −0.1952`)。

ESTIMAND        `EARLY_rel` = 每人「他的兴趣**相对于人群时间表**来得多早」
                (= 该类别的**人群中位起始** − 他的实际起始,再跨类别平均)——
                **按构造已扣掉人群时间表**。用它替换 `EARLY` 重跑 `#334` 的全族检验。
KILL            **若 `S × EARLY_rel` 的交互塌掉 -> 重叠那一段**就是**人群时间表,机制说清了;
                若仍在 -> 重叠不是时间表,而是别的东西。**
POSITIVE CTRL   种入 `S × EARLY_rel` = 0.15 -> 必须过全族阈。
NEGATIVE CTRL   置换(全族阈的构造)。
⚠ 先报         **`corr(EARLY_rel, S)`** —— 若它仍与 `S` 强相关,这一步没扣干净,后面读不了。
⚠ 族           仍是**三个**:`S × c3⁻` · `S × EARLY_rel` · `c3⁻ × EARLY_rel`。
IMPOSSIBLE      人群中位起始是**在同一批人身上**估的,所以 `EARLY_rel` 含一点自我参照;
                n=6,717 下这点偏差很小,但不是零。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
nc=np.isfinite(ONS).sum(1)
MED=np.nanmedian(ONS,0)                          # ⚠ 人群时间表:每类别的中位起始
REL=MED[None,:]-ONS                              # 正 = 比人群典型**更早**
MO=np.where(nc>=5,np.nanmean(ONS,1),np.nan)
RELm=np.where(nc>=5,np.nanmean(REL,1),np.nan)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(RELm)&np.isfinite(MO)&np.isfinite(sh)&ok
n=int(m0.sum()); zz=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zs,zc,ncz=zz(S),zz(C3),zz(nc.astype(float)); y=zz(sh)
ze_abs,ze_rel=zz(-MO),zz(RELm)
print(f"n={n:,} · 人群时间表(类别中位起始)范围 [{MED.min():.1f}, {MED.max():.1f}]")
print(f"⚠ 先报:`corr(EARLY_rel, S)` = **{np.corrcoef(ze_rel,zs)[0,1]:+.4f}** "
      f"(对比 `corr(EARLY_abs, S)` = **{np.corrcoef(ze_abs,zs)[0,1]:+.4f}**)")
print(f"   `corr(EARLY_rel, EARLY_abs)` = **{np.corrcoef(ze_rel,ze_abs)[0,1]:+.4f}**")
def famtest(ze,tag,NP=1000,yy=None):
    X=np.column_stack([np.ones(n),zs,zc,ze,ncz,zs*zc,zs*ze,zc*ze])
    XtXi=np.linalg.pinv(X.T@X)
    def ts(v):
        b,*_=np.linalg.lstsq(X,v,rcond=None); r=v-X@b
        s2=float(r@r)/(n-X.shape[1]); se=np.sqrt(np.diag(s2*XtXi)); return b[[5,6,7]]/se[[5,6,7]],b[[5,6,7]]
    t0,b0=ts(y if yy is None else yy)
    rg=np.random.default_rng(20260804)
    T_=np.array([ts((y if yy is None else yy)[rg.permutation(n)])[0] for _ in range(NP)])
    mx=np.abs(T_).max(1); thr=float(np.percentile(mx,95))
    pfw=[float((mx>=abs(t0[i])).mean()) for i in range(3)]
    L=['S × c3⁻',f'S × {tag}',f'c3⁻ × {tag}']
    print(f"\n【{tag}】全族 95% 阈 **{thr:.3f}**")
    for i in range(3):
        print(f"   {L[i]:<16} 系数 **{b0[i]:+.4f}** · |t| **{abs(t0[i]):.3f}** · "
              f"全族 p **{pfw[i]:.3f}** · 过阈 **{'是' if abs(t0[i])>thr else '否'}**")
    return b0,t0,thr,pfw
bA,tA,thrA,pA=famtest(ze_abs,'EARLY_abs')
bR,tR,thrR,pR=famtest(ze_rel,'EARLY_rel')
rg=np.random.default_rng(19)
yp=0.12*zs+0.12*zc+0.12*ze_rel+0.15*zs*ze_rel+rg.standard_normal(n); yp=(yp-yp.mean())/yp.std()
bP,tP,thrP,pP=famtest(ze_rel,'EARLY_rel(正对照 种入0.15)',NP=300,yy=yp)
T=pd.DataFrame([dict(v_arm='EARLY_abs',v_b=float(bA[1]),v_t=float(tA[1]),v_pfw=pA[1]),
                dict(v_arm='EARLY_rel',v_b=float(bR[1]),v_t=float(tR[1]),v_pfw=pR[1])])
check_columns(T,'R381'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rel.csv',index=False)
gg=Gate('重叠那一段是不是人群时间表')
gg.asserted('★ 正对照:种入 `S × EARLY_rel` = 0.15 必须过全族阈',abs(tP[1])>thrP,
            f"|t| {abs(tP[1]):.2f} vs 阈 {thrP:.3f}")
gg.asserted('⚠ 先报的控制:`EARLY_rel` 与 `S` 的相关是否明显小于 `EARLY_abs` 的',
            abs(np.corrcoef(ze_rel,zs)[0,1])<abs(np.corrcoef(ze_abs,zs)[0,1]),
            f"rel {np.corrcoef(ze_rel,zs)[0,1]:+.4f} vs abs {np.corrcoef(ze_abs,zs)[0,1]:+.4f}")
gg.asserted('★ 注册的 kill:`S × EARLY_rel` 的交互塌不塌',
            abs(tR[1])<thrR,
            f"abs 版 |t| **{abs(tA[1]):.3f}**(过阈 {'是' if abs(tA[1])>thrA else '否'})-> "
            f"rel 版 |t| **{abs(tR[1]):.3f}**(过阈 {'是' if abs(tR[1])>thrR else '否'})"
            f" —— **塌了则重叠那段就是人群时间表**")
gg.asserted('⚠ 边界:人群中位起始在同一批人身上估',True,'`EARLY_rel` 含一点自我参照,n 大时很小但不是零')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
