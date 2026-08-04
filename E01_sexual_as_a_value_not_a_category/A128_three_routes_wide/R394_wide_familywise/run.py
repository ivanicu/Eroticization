import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A128 R394 -- 三条路与那条次可加,在宽口径的更大人群里还在不在

`#334` 是**已上页面**的:三条路增量相当,`S × c3⁻` 与 `c3⁻ × EARLY` 相加,
**`S × EARLY` 次可加(全族 p = 0.004)**。它是在 `cov>=8`(n=6,717)上做的。
`#346a`:宽口径(`cov>=4`,n=12,720)对 `c3` 几乎免费。

⚠ **人群不同**(`#349c`):宽口径多了 6,003 个低覆盖的人,而 `#346b` 说他们**系统上不一样**。
**所以这不是「同一个问题、更多数据」** —— 读法里必须带着这句。
⚠ 门槛与对比**全部**由**同一次运行的窄口径参照臂**给(`#329b`)。

ESTIMAND        `cov>=4` 与 `cov>=8` 同一次运行:三条路的增量 ΔR² ·
                三个预先声明的两两交互的 **全族 max-|t| 阈**(置换 1000 次)。
KILL            **若 `S × EARLY` 在宽口径下仍过全族阈 -> `#334` 在更大人群里复现,页面更稳;
                若不过 -> 它是窄口径人群的性质,页面上那句要加人群限定。**
POSITIVE CTRL   两个口径各自种入 `S × EARLY` = 0.15 -> 都必须过阈。
NEGATIVE CTRL   置换(全族阈的构造)。
IMPOSSIBLE      两个口径的人群不同 -> 差异可以是**人群**的,不是**功效**的;本轮分不开这两者。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in onsc])
ncat=np.isfinite(ONS).sum(1); MO=np.where(ncat>=5,np.nanmean(ONS,1),np.nan)
cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
def run(TH,NP=1000,plant=0.0):
    mk=cov>=TH
    R_=np.flatnonzero(mk&np.isfinite(sh))
    v=load_of(R_); c3=np.where(mk,score_of(v),np.nan)
    k=np.isfinite(c3)&np.isfinite(sh)
    if float(np.corrcoef(c3[k],sh[k])[0,1])<0: c3=-c3
    S=Spos(mk); EARLY=-MO
    m=mk&np.isfinite(S)&np.isfinite(c3)&np.isfinite(EARLY)&np.isfinite(sh)
    n=int(m.sum()); z=lambda x:(x[m]-x[m].mean())/max(x[m].std(),1e-12)
    zs,zc,ze,ncz=z(S),z(c3),z(EARLY),z(ncat.astype(float)); y=z(sh)
    if plant: y=(y+plant*zs*ze); y=(y-y.mean())/y.std()
    MAIN=[zs,zc,ze,ncz]
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,y,rcond=None)
        r=y-X@b; return 1-float(r@r)/float(((y-y.mean())**2).sum())
    full=r2(MAIN)
    inc=[full-r2([c for j,c in enumerate(MAIN) if j!=i]) for i in range(3)]
    X=np.column_stack([np.ones(n)]+MAIN+[zs*zc,zs*ze,zc*ze]); XtXi=np.linalg.pinv(X.T@X)
    def ts(vv):
        b,*_=np.linalg.lstsq(X,vv,rcond=None); r=vv-X@b
        s2=float(r@r)/(n-X.shape[1]); se=np.sqrt(np.diag(s2*XtXi)); return b[[5,6,7]]/se[[5,6,7]],b[[5,6,7]]
    t0,b0=ts(y)
    rg=np.random.default_rng(20260804)
    T_=np.array([ts(y[rg.permutation(n)])[0] for _ in range(NP)])
    mx=np.abs(T_).max(1); thr=float(np.percentile(mx,95))
    pfw=[float((mx>=abs(t0[i])).mean()) for i in range(3)]
    return dict(n=n,inc=inc,t=t0,b=b0,thr=thr,pfw=pfw)
L=['S × c3⁻','S × EARLY','c3⁻ × EARLY']; LI=['S','c3⁻','EARLY']
for TH in (8,4):
    r=run(TH)
    print(f"\n【cov>={TH}】n={r['n']:,} · 全族 95% 阈 **{r['thr']:.3f}**")
    print(f"   增量 ΔR²:" + ' · '.join(f"{LI[i]} **{100*r['inc'][i]:.3f}pp**" for i in range(3)))
    for i in range(3):
        print(f"   {L[i]:<12} 系数 **{r['b'][i]:+.4f}** · |t| **{abs(r['t'][i]):.3f}** · "
              f"全族 p **{r['pfw'][i]:.3f}** · 过阈 **{'是' if abs(r['t'][i])>r['thr'] else '否'}**")
    if TH==8: r8=r
    else: r4=r
print(f"\n正对照(两口径各自种入 `S × EARLY` = 0.15,300 次置换):")
for TH in (8,4):
    p=run(TH,NP=300,plant=0.15)
    print(f"   cov>={TH}: |t| **{abs(p['t'][1]):.2f}** vs 阈 {p['thr']:.3f} -> "
          f"**{'过' if abs(p['t'][1])>p['thr'] else '不过'}**")
    if TH==4: pc4=p
T=pd.DataFrame([dict(v_th=8,t=float(r8['t'][1]),pfw=r8['pfw'][1],n=r8['n']),
                dict(v_th=4,t=float(r4['t'][1]),pfw=r4['pfw'][1],n=r4['n'])])
check_columns(T,'R394'); T.to_csv(pathlib.Path(__file__).parent/'results'/'wide_fw.csv',index=False)
gg=Gate('三条路与那条次可加,在宽口径里还在不在')
gg.asserted('★ 正对照:两口径各自种入 0.15 都必须过阈',abs(pc4['t'][1])>pc4['thr'],
            f"宽 |t| {abs(pc4['t'][1]):.2f} vs 阈 {pc4['thr']:.3f}")
gg.asserted('★ 注册的 kill:`S × EARLY` 在宽口径下是否仍过全族阈',
            abs(r4['t'][1])>r4['thr'],
            f"宽 |t| **{abs(r4['t'][1]):.3f}** vs 阈 **{r4['thr']:.3f}**(全族 p **{r4['pfw'][1]:.3f}**);"
            f"参照臂 窄 |t| {abs(r8['t'][1]):.3f} vs 阈 {r8['thr']:.3f}")
gg.asserted('⚠ 另两个交互在宽口径下是否仍为零',
            abs(r4['t'][0])<r4['thr'] and abs(r4['t'][2])<r4['thr'],
            f"S×c3⁻ p {r4['pfw'][0]:.3f} · c3⁻×EARLY p {r4['pfw'][2]:.3f}")
gg.asserted('⚠⚠ 人群不同,不是「同一个问题、更多数据」',True,
            f"宽口径多了 {r4['n']-r8['n']:,} 个低覆盖的人,而 `#346b` 说他们系统上不一样 —— "
            f"差异可以是**人群**的,不是**功效**的;本轮分不开这两者")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
