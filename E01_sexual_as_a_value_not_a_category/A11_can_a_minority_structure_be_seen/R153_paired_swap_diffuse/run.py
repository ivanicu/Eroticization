import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R19 -- 保边际的弥散对照。#121e 指定的修法,也是这个问题的最后一次尝试(#111c)。

#121d:我的弥散种植改变列和,而真实-vs-零的对比是在精确保边际下做的,两者不可比。
#121e 的修法:**成对交换** —— A 得到稀有项、失去常见项,B 反向。行和与列和都精确不变。
这样"弥散"和"集中"两个种植世界与 curveball 零在**同一个尺度**上,可以直接比。

  弥散(保边际): 很多人各做少量成对交换 -> 人层面轻微异质,集中度上尾不凸
  集中(保边际): 少数人各做大量成对交换 -> 上尾凸起
  两者交换**总数相同**,只是分配给多少人不同 —— 这是"同等弥散量下还有没有多余集中"的精确形式。

⚠ #121a 的教训写进设计:**只条件于勾选数 k,绝不条件于平均意外度 S**(它与被测结构共变)。
⚠ #121b 的教训:每个正对照都必须先证明它**产生了预测的签名**,再读真实数据。

ESTIMAND        每人稀有项计数的(p95 − p50),在勾选数分层内标准化;真实 vs 零 vs 两个保边际种植。
KILL            threshold-free。前置:集中种植必须产生正的尾部超出且显著大于弥散种植的,
                否则整轮 UNVERIFIED 并停止(#111c)。
POSITIVE CTRL   集中种植;必须凸起。
NEGATIVE CTRL   第二次独立 curveball;以及交换总数为 0 的退化种植。
IMPOSSIBLE      "强烈"与"更愿意勾罕见项"分不开。
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns, check_coverage
warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
def curveball(M,rng,per_row=5.):
    A=[set(np.flatnonzero(r).tolist()) for r in M]; n=len(A)
    for _ in range(int(per_row*n)):
        i,j=int(rng.integers(n)),int(rng.integers(n))
        if i==j: continue
        ai,aj=A[i],A[j]; inter=ai&aj
        di=list(ai-inter); dj=list(aj-inter); L=di+dj
        if not L: continue
        rng.shuffle(L); k=len(di)
        A[i]=inter|set(L[:k]); A[j]=inter|set(L[k:])
    out=np.zeros_like(M)
    for i,s in enumerate(A): out[i,list(s)]=1.
    return out
def paired_swap(M,rng,n_swaps,frac):
    """保边际:A 得稀有失常见,B 得常见失稀有。行和列和都精确不变。
    n_swaps = 总交换次数;frac = 承担这些交换的人的比例(小 -> 集中,大 -> 弥散)。"""
    Mw=M.copy(); n,m=M.shape
    o=np.argsort(M.mean(0)); rare=o[:max(3,m//5)]; common=o[-max(3,m//5):]
    pool=np.flatnonzero(rng.random(n)<frac)
    if len(pool)<4: return Mw
    done=0; tries=0
    while done<n_swaps and tries<n_swaps*40:
        tries+=1
        i,j=pool[rng.integers(len(pool))],pool[rng.integers(len(pool))]
        if i==j: continue
        r=rare[rng.integers(len(rare))]; c=common[rng.integers(len(common))]
        if Mw[i,r]==0 and Mw[i,c]==1 and Mw[j,r]==1 and Mw[j,c]==0:
            Mw[i,r]=1.; Mw[i,c]=0.; Mw[j,r]=0.; Mw[j,c]=1.; done+=1
    return Mw
def profile(build,seed=1):
    cnt=np.zeros(len(ALLP)); rare=np.zeros(len(ALLP)); used=0
    for t in IDENT:
        M0=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        M=build(M0,seed)
        assert np.allclose(M.sum(0),M0.sum(0)) and np.allclose(M.sum(1),M0.sum(1)),"边际未保住"
        thr=np.quantile(M0.mean(0),0.20); israre=(M0.mean(0)<=thr).astype(float)
        cnt[idx]+=M.sum(1); rare[idx]+=M@israre; used+=1
    check_coverage(used,len(IDENT),'A11R19')
    ok=cnt>=15
    return check_columns(pd.DataFrame(dict(v_k=cnt[ok],v_rare=rare[ok])),'A11R19')
def tailminus(D):
    d=D.copy(); d['cell']=np.digitize(d.v_k,np.quantile(d.v_k,[.25,.5,.75]))
    z=d.groupby('cell').v_rare.transform(lambda s:(s-s.mean())/(s.std()+1e-9))
    q=np.percentile(z,[50,95]); return q[1]-q[0]
NSW=3000
B={'real':  lambda M,s: M,
   'null':  lambda M,s: curveball(M,np.random.default_rng(7000+s)),
   'null2': lambda M,s: curveball(M,np.random.default_rng(8000+s)),
   'deg':   lambda M,s: paired_swap(curveball(M,np.random.default_rng(7000+s)),np.random.default_rng(11),0,0.5),
   'conc':  lambda M,s: paired_swap(curveball(M,np.random.default_rng(7000+s)),np.random.default_rng(12),NSW,0.05),
   'diff':  lambda M,s: paired_swap(curveball(M,np.random.default_rng(7000+s)),np.random.default_rng(13),NSW,0.60)}
out={}
for n,f in B.items():
    D=profile(f); out[n]=(tailminus(D),float(D.v_rare.mean()),len(D))
    print(f"  {n:6s} 尾−中 {out[n][0]:+.3f}  平均稀有项 {out[n][1]:.3f}  n={out[n][2]:,}",flush=True)
base=out['null'][0]; jit=abs(out['null2'][0]-base)
print(f"\n=== 保边际,交换总数相同({NSW}/块),只改分配给多少人 ===")
print(f"  {'臂':7s} {'尾−中':>8s} {'相对零':>8s}   {'平均稀有项':>10s}")
for n in ['null','null2','deg','conc','diff','real']:
    print(f"  {n:7s} {out[n][0]:+8.3f} {out[n][0]-base:+8.3f}   {out[n][1]:10.3f}")
print(f"\n  零的抖动 {jit:.3f}")
ec=out['conc'][0]-base; ed=out['diff'][0]-base; er=out['real'][0]-base
g=Gate("#121e 保边际下:同等交换量,集中 vs 弥散")
g.asserted("退化种植(0 次交换)等于零", abs(out['deg'][0]-base)<2*jit,
           f"{out['deg'][0]-base:+.3f} vs 2*抖动 {2*jit:.3f}")
g.asserted("集中种植产生了预测的签名(正的尾部超出)(#121b)", ec>2*jit,
           f"集中 {ec:+.3f} vs 2*抖动 {2*jit:.3f}")
g.asserted("集中与弥散在同等交换量下可分", abs(ec-ed)>2*jit,
           f"集中 {ec:+.3f} vs 弥散 {ed:+.3f},差 {abs(ec-ed):.3f}")
g.require_resolvable_first("真实相对零", effect=er, spread=jit)
g.artifact_cannot_explain("弥散能解释多少", artifact=ed, effect=er, spread=jit)
print(); print(g)
if g.verdict():
    print(f"\n  真实 {er:+.3f}   集中 {ec:+.3f}   弥散 {ed:+.3f}")
    print(f"  -> 真实更接近{'集中:存在对少数选项异常强烈的人' if abs(er-ec)<abs(er-ed) else '弥散:只有连续维度'}")
else:
    print(f"\n  -> UNVERIFIED。按 #111c,这是第二次尝试,停止,不追第三轮。")
print(f"\nartifact sha1 {hashlib.sha1(str(out).encode()).hexdigest()[:12]}")
