import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A111 R365 -- 去衰减数字的区间:`rel` 在分母里

`#317a` 证明了去衰减在这份数据上是一次观察。**那就用它** —— 但用之前先修一个形状问题:
**去衰减是 `r / √rel`,而 `rel` 在分母里。** `rel` 小的时候,`rel` 自己的抽样误差会被
**非线性地放大**,而且放大是**不对称**的(rel 偏小 -> 读数暴涨)。
**报一个点估计是误导。**

ESTIMAND        对 `form ↔ 羞耻` 这个可完整重建的例子:**人层自助 ≥400 次**,
                每次**同时**重估 ① 两个指标的相关(-> SB 信度)② 原始相关 ③ 去衰减值;
                报 **2.5/50/97.5 分位**,并与「把 rel 当常数」的朴素区间对比。
KILL            **若自助区间明显比朴素区间宽或不对称 -> 点估计不可发布,页面上要换成区间;
                若两者几乎一样 -> `rel≈0.38` 还没小到让非线性咬人,点估计可以留。**
POSITIVE CTRL   在 `#363` 的合成网格里取一格(真相关已知)-> **自助区间必须覆盖真值**,
                且覆盖率在 20 次重复里接近 95%。
NEGATIVE CTRL   `perm_finite` 打乱人 -> 区间必须覆盖 0。
IMPOSSIBLE      自助只传播**抽样**不确定性,不传播「两个指标是不是真的平行」这个模型假设;
                后者由 `#317a` 的不变性检验单独支持。
"""
import numpy as np, pandas as pd, hashlib
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ani=pd.to_numeric(d['animated'],errors='coerce').values.astype(float)
wri=pd.to_numeric(d['written'],errors='coerce').values.astype(float)
m=np.isfinite(ani)&np.isfinite(wri)&np.isfinite(sh)
A0,W0,Y=ani[m],wri[m],sh[m]; n=len(Y)
def dis_of(a,w,y):
    za=(a-a.mean())/max(a.std(),1e-12); zw=(w-w.mean())/max(w.std(),1e-12)
    r=float(np.corrcoef(za,zw)[0,1]); rel=min(2*r/(1+r),1.0) if r>0 else np.nan
    f=za+zw; raw=float(np.corrcoef(f,(y-y.mean())/max(y.std(),1e-12))[0,1])
    return raw,rel,(raw/np.sqrt(rel) if np.isfinite(rel) and rel>0 else np.nan)
raw0,rel0,dis0=dis_of(A0,W0,Y)
print(f"n={n:,} · 原始 r **{raw0:+.4f}** · SB 信度 **{rel0:.4f}** · 点估计去衰减 **{dis0:+.4f}**")
B=600; rg=np.random.default_rng(808)
bs=np.array([dis_of(*(v[i] for v in (A0,W0,Y))) for i in
             (rg.integers(0,n,n) for _ in range(B))])
q=np.nanpercentile(bs[:,2],[2.5,50,97.5]); qrel=np.nanpercentile(bs[:,1],[2.5,50,97.5])
qraw=np.nanpercentile(bs[:,0],[2.5,50,97.5])
se_raw=float(np.nanstd(bs[:,0]))
naive=(raw0-1.96*se_raw)/np.sqrt(rel0),(raw0+1.96*se_raw)/np.sqrt(rel0)
print(f"\n自助 {B} 次(同时重估信度与相关):")
print(f"   信度      2.5% **{qrel[0]:.4f}** · 中位 **{qrel[1]:.4f}** · 97.5% **{qrel[2]:.4f}**")
print(f"   原始 r    2.5% **{qraw[0]:+.4f}** · 中位 **{qraw[1]:+.4f}** · 97.5% **{qraw[2]:+.4f}**")
print(f"   **去衰减  2.5% {q[0]:+.4f} · 中位 {q[1]:+.4f} · 97.5% {q[2]:+.4f}**  宽度 **{q[2]-q[0]:.4f}**")
print(f"   朴素区间(把 rel 当常数):[{naive[0]:+.4f}, {naive[1]:+.4f}] 宽度 **{naive[1]-naive[0]:.4f}**")
asym=abs((q[2]-q[1])-(q[1]-q[0]))/max(q[2]-q[0],1e-9)
print(f"   -> 自助/朴素 宽度比 **{(q[2]-q[0])/max(naive[1]-naive[0],1e-9):.2f}×** · "
      f"不对称度 **{100*asym:.1f}%**")
# ⚠ 页面上报的是 **#316 的口径**(还要求块覆盖 `cov>=8`),n 与本轮不同 ->
#    页面的那个数必须在**它自己的口径**上算区间,否则是拿一套的区间配另一套的点估计。
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
cov=np.zeros(NN)
for _,qq in keep.iterrows():
    ss=lg[lg.qi==qq.qi]; vc=ss.option.value_counts(); ss=ss[ss.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(ss.person.unique()))
    if len(ppl)>=1200: cov[ppl]+=1
m316=np.isfinite(ani)&np.isfinite(wri)&np.isfinite(sh)&(cov>=8)
A1,W1,Y1=ani[m316],wri[m316],sh[m316]; n1=len(Y1)
r1,rel1,d1=dis_of(A1,W1,Y1)
rg3=np.random.default_rng(4242)
b1=np.array([dis_of(A1[i],W1[i],Y1[i]) for i in (rg3.integers(0,n1,n1) for _ in range(600))])
q1=np.nanpercentile(b1[:,2],[2.5,50,97.5])
print(f"\n★ **页面口径**(`#316`,加 `cov>=8`)n={n1:,}:原始 **{r1:+.4f}** · 信度 **{rel1:.4f}** · "
      f"去衰减 **{d1:+.4f}**  **95% 区间 [{q1[0]:+.4f}, {q1[2]:+.4f}]**")

print(f"\n正对照:合成一格(真相关 0.20,目标信度 0.40),20 次重复的覆盖率")
cov=0; TRUE=0.20; R=0.40
for t in range(20):
    rg2=np.random.default_rng(3000+t); L=rg2.standard_normal(n)
    y=TRUE*L+np.sqrt(1-TRUE**2)*rg2.standard_normal(n)
    ri=R/(2-R); s=np.sqrt(1/ri-1)
    a=L+s*rg2.standard_normal(n); w=L+s*rg2.standard_normal(n)
    bb=np.array([dis_of(a[i],w[i],y[i]) for i in (rg2.integers(0,n,n) for _ in range(200))])
    lo,hi=np.nanpercentile(bb[:,2],[2.5,97.5]); cov+= (lo<=TRUE<=hi)
print(f"   覆盖真值 **{cov}/20**({100*cov/20:.0f}%)")
rgn=np.random.default_rng(99)
Yp=Y[rgn.permutation(n)]
bn=np.array([dis_of(A0[i],W0[i],Yp[i]) for i in (rgn.integers(0,n,n) for _ in range(300))])
ln,hn=np.nanpercentile(bn[:,2],[2.5,97.5])
print(f"负对照(打乱人):去衰减区间 [{ln:+.4f}, {hn:+.4f}] · 覆盖 0 **{'是' if ln<=0<=hn else '否'}**")
T=pd.DataFrame([dict(v_q='2.5%',dis=q[0],rel=qrel[0],raw=qraw[0]),
                dict(v_q='50%',dis=q[1],rel=qrel[1],raw=qraw[1]),
                dict(v_q='97.5%',dis=q[2],rel=qrel[2],raw=qraw[2])])
check_columns(T,'R365'); T.to_csv(pathlib.Path(__file__).parent/'results'/'intervals.csv',index=False)
gg=Gate('去衰减数字的区间')
gg.asserted('★ 正对照:自助区间对已知真值的覆盖率(20 次,期望 ~19/20)',cov>=17,f"{cov}/20")
gg.asserted('★ 负对照:打乱人后区间必须覆盖 0',ln<=0<=hn,f"[{ln:+.4f}, {hn:+.4f}]")
gg.asserted('★ 注册的 kill:自助区间是否明显比朴素区间宽或不对称',
            (q[2]-q[0])/max(naive[1]-naive[0],1e-9)>1.10 or asym>0.10,
            f"宽度比 {(q[2]-q[0])/max(naive[1]-naive[0],1e-9):.2f}× · 不对称度 {100*asym:.1f}%")
gg.asserted('⚠ 边界:自助只传播抽样不确定性,不传播「两个指标是否真平行」',True,
            '后者由 `#317a` 的不变性检验单独支持')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
