import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A127 R390 -- `c3⁻` 的两端与羞耻的关系,对不对称

羞耻的三条路里 `c3⁻` 的增量最大(**0.949pp**,是 `S` 的两倍),
而这个项目对它做的全部工作是**命名失败五次**与**描述两次**。
`S` 与 `EARLY` 都被追到了机制层,**`c3⁻` 没有**。

⚠ **格层做不了**:两套仪器题目不相交,类别没有 `c3⁻` 载荷(`#384` 已如实报过不可算)。
⚠ `#334` 已答过「与 `S` / `EARLY` 的交互」——**都是零**。
**所以真正没做过的是:`c3⁻` 的**两端**分别与羞耻的关系是否对称。**

ESTIMAND        **斜率**(不是相关 —— 半样本内的相关会被**全距受限**压低,斜率不会):
                以中位为节点的**分段线性**模型 `羞耻 ~ z + z·1[z>中位]`,
                报两端斜率与它们的**差**;并对 `perm_finite` 置换零。
KILL            **若两端斜率明显不同 -> 「广度型」这条路是**单向**的,那是一个新的、可上页面的形状;
                若相同 -> 它是一条真正的**双向**维度,那也值得写。**
POSITIVE CTRL   合成一个**已知单向**的关系 -> 差必须被抓到;合成一个**线性**的 -> 差必须落零。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ KNOB         节点位置(中位 / 三分位 / 四分位)-> 报规格曲线。
IMPOSSIBLE      分段线性只测**一个**节点处的折点;若关系是平滑弯曲的,它会低估。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]; S=Q[0]
m0=np.isfinite(C3)&np.isfinite(sh)&np.isfinite(S)&ok
n=int(m0.sum()); zz=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zc,zs,y=zz(C3),zz(S),zz(sh)
def seg(yv,q=0.5,x=None):
    x=zc if x is None else x
    k=np.quantile(x,q); hi=(x>k).astype(float)
    X=np.column_stack([np.ones(n),x,(x-k)*hi,hi])
    b,*_=np.linalg.lstsq(X,yv,rcond=None); r=yv-X@b
    s2=float(r@r)/(n-X.shape[1]); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(b[1]+b[2]),float(b[2]),float(se[2])
lo,hi,dif,sed=seg(y)
print(f"n={n:,} · 节点 = `c3⁻` 中位")
print(f"★ 下半斜率 **{lo:+.4f}** · 上半斜率 **{hi:+.4f}** · **差 {dif:+.4f}**(se {sed:.4f},|t| {abs(dif/sed):.2f})")
def perm_finite(v,seed):
    z2=v.copy(); z2=z2[np.random.default_rng(seed).permutation(len(z2))]; return z2
nul=np.array([seg(perm_finite(y,700+i))[2] for i in range(200)])
q=float((np.abs(nul)>=abs(dif)).mean())
print(f"负对照(打乱人 200 次):差 **{nul.mean():+.4f} ± {nul.std():.4f}** · |零| ≥ |观测| **{q:.3f}**")
print(f"\n⚠ 旋钮 · 节点位置:")
SPEC=[]
for qq in (0.33,0.5,0.67):
    a,b_,dd,ss=seg(y,q=qq); SPEC.append((qq,dd,ss))
    print(f"   分位 {qq:.2f}: 下 **{a:+.4f}** · 上 **{b_:+.4f}** · 差 **{dd:+.4f}**(|t| {abs(dd/ss):.2f})")
rg=np.random.default_rng(23)
print(f"\n正对照:")
for tag,yv in (('已知**单向**(只有上半有斜率)',np.where(zc>0,0.25*zc,0)+rg.standard_normal(n)),
               ('已知**线性**',0.15*zc+rg.standard_normal(n))):
    yv=(yv-yv.mean())/yv.std(); a,b_,dd,ss=seg(yv)
    print(f"   {tag:<26} 下 **{a:+.4f}** · 上 **{b_:+.4f}** · 差 **{dd:+.4f}**(|t| {abs(dd/ss):.2f})")
    if '单向' in tag: p1=(dd,ss)
    else: p2=(dd,ss)
mde=2.8*sed
T=pd.DataFrame([dict(v_arm=f'分位{q_:.2f}',v_diff=dd_,v_se=ss_) for q_,dd_,ss_ in SPEC])
check_columns(T,'R390'); T.to_csv(pathlib.Path(__file__).parent/'results'/'ends.csv',index=False)
gg=Gate('`c3⁻` 的两端对不对称')
gg.asserted('★ 正对照:已知单向必须被抓到、已知线性必须落零',
            abs(p1[0]/p1[1])>3 and abs(p2[0]/p2[1])<2,
            f"单向 差 {p1[0]:+.4f}(|t| {abs(p1[0]/p1[1]):.2f})· 线性 差 {p2[0]:+.4f}(|t| {abs(p2[0]/p2[1]):.2f})")
gg.negative_control('★ 负对照:打乱人后的两端斜率差',float(nul.mean()),dif,
    null_spread=float(nul.std()),null_kind='打乱人 —— 保持 `c3⁻` 的分布,只打散羞耻')
gg.asserted('★ 注册的 kill:两端斜率是否明显不同',abs(dif)>2*sed,
            f"差 **{dif:+.4f}** vs 2×se **{2*sed:.4f}**(|零| ≥ |观测| {q:.3f})")
gg.asserted('⚠ 规格曲线:三个节点位置下差的符号是否一致',
            len(set(np.sign([x[1] for x in SPEC])))==1,
            ' · '.join(f"{q_:.2f} {dd_:+.4f}" for q_,dd_,_ in SPEC))
gg.null_claim_uses_null_criteria('★ guard 21:若判对称(零),三件套在不在',
    'NULL' if abs(dif)<=2*sed else 'EFFECT',perm_quantile=q,mde=mde,
    sensitivity_shown=f"单向合成 |t| {abs(p1[0]/p1[1]):.1f}",meaningful=0.10)
gg.asserted('⚠ 边界:分段线性只测一个节点处的折点',True,'若关系是平滑弯曲的,它会低估')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
