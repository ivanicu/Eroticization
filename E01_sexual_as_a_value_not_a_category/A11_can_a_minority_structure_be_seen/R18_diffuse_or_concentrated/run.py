import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R18 -- 「整体偏好罕见」还是「对某一样异常强烈」?这是恋物的定义,而它从没被测过。

#120 的 NEXT 指向 #104 的"勾选数分层内的意外度秩"。设计时重读 #99 发现:那条路线要回答的
"有没有子群",#99 已经用 S 回答了 —— 对称加宽,是一条连续维度不是一个子群。

但 #99 测的是 S 的**水平**,没测它的**形状**。而恋物在现象学上不是"整体偏好罕见",是
"**对某一样东西异常强烈**"。这两件事在 S 的均值上无法区分,在**人内集中度**上可以:

  弥散(连续维度): 一个人的稀有度均匀铺在他所有勾选上 -> 集中度与零一致
  集中(恋物子群): 少数人的稀有度堆在极少数选项上   -> 集中度分布出现右尾

统计量必须同时:逐人 · 条件于勾选数 · 条件于总意外度(否则 S 高的人机械地有更多集中空间)。
curveball 精确保留每人勾选数与每项基率,只毁掉"谁勾了什么",所以它是这个问题精确匹配的零。

ESTIMAND        每人勾选项的意外度**集中度**(最大值、变异系数、稀有项计数),
                在 (勾选数 x 平均意外度) 分层内与固定边际零比较,报**整个分位数曲线**(#99 的教训)。
IDENTIFICATION  identified;零保留两个条件量,只毁掉分配。
WORLDS          diffuse      集中度分布与零一致 -> #99 的连续维度是全部
                concentrated 上尾超出零 -> 存在"对某几样异常强烈"的人,而 S 的均值看不见他们
KILL            threshold-free,逐分位数对自身自助 SE;gate 顺序按 #120d。
POSITIVE CTRL   (1) 集中种植:5% 的人各自把 3 个最稀有选项勾满 -> 上尾必须凸起,中位数不动
                (2) 弥散种植:所有人整体偏向稀有 -> 中位数移动,上尾按比例移动
                两者必须彼此可分,否则读数不被许可(#120b 的教训:用它们自己的 SE 判,不用常数)
NEGATIVE CTRL   固定边际零本身;并报第二次独立抽样以确认零的稳定性。
IMPOSSIBLE      "强烈"无法与"回答问卷时更愿意勾罕见项"分开;本轮测的是分布形状,不是体验强度。
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
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
print(f"块 {len(IDENT)}  人 {len(ALLP):,}",flush=True)
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
def plant_conc(M,who,rng,k=3):
    """集中种植:指定的人把这个块里最稀有的 k 个选项勾满(换掉常见的,保持勾选数)。"""
    Mw=M.copy(); o=np.argsort(M.mean(0)); rare=o[:k]; common=o[-6:]
    for i in who:
        for r in rare:
            if Mw[i,r]==0:
                c=[c for c in common if Mw[i,c]==1]
                if c: Mw[i,c[0]]=0.; Mw[i,r]=1.
    return Mw
def plant_diff(M,rng,g=0.25):
    """弥散种植:所有人以概率 g 把一个常见勾选换成一个稀有的。"""
    Mw=M.copy(); o=np.argsort(M.mean(0)); rare=o[:6]; med=o[len(o)//2:]
    for i in range(M.shape[0]):
        if rng.random()<g:
            c=[c for c in med if Mw[i,c]==1]; r=[r for r in rare if Mw[i,r]==0]
            if c and r: Mw[i,c[0]]=0.; Mw[i,r[0]]=1.
    return Mw
QS=[10,25,50,75,90,95,99]
def profile(build,seed):
    """每人:勾选数 k、平均意外度 S、集中度(最大意外度、变异系数、稀有项数)。"""
    tot=np.zeros(len(ALLP)); cnt=np.zeros(len(ALLP))
    mx=np.zeros(len(ALLP)); sq=np.zeros(len(ALLP)); rare=np.zeros(len(ALLP))
    used=0
    for t in IDENT:
        M0=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        M=build(M0,t,seed)
        ref=-np.log(np.clip(M0.mean(0),1e-4,1.))          # 基率永远取自真实矩阵,零不移动尺子
        thr=np.quantile(M0.mean(0),0.20)                   # 该块最稀有的 20% 选项
        israre=(M0.mean(0)<=thr).astype(float)
        s=M@ref; k=M.sum(1)
        tot[idx]+=s; cnt[idx]+=k
        mx[idx]=np.maximum(mx[idx],(M*ref[None,:]).max(1))
        sq[idx]+=(M*(ref[None,:]**2)).sum(1)
        rare[idx]+=M@israre
        used+=1
    check_coverage(used,len(IDENT),'A11R18 profile')
    ok=cnt>=15
    S=tot[ok]/cnt[ok]
    var=sq[ok]/cnt[ok]-S**2
    return pd.DataFrame(dict(v_k=cnt[ok],v_S=S,v_max=mx[ok],
                             v_cv=np.sqrt(np.maximum(var,0))/np.maximum(S,1e-9),
                             v_rare=rare[ok]))
def matched_quantiles(D,col):
    """在 (勾选数 x 平均意外度) 分层内标准化,再取分位数 —— 条件于两个量。"""
    d=D.copy()
    kb=np.digitize(d.v_k,np.quantile(d.v_k,[.25,.5,.75]))
    sb=np.digitize(d.v_S,np.quantile(d.v_S,[.25,.5,.75]))
    d['cell']=kb*4+sb
    z=d.groupby('cell')[col].transform(lambda s:(s-s.mean())/(s.std()+1e-9))
    return np.percentile(z,QS)
BUILDS={
 'real':   lambda M,t,s: M,
 'null':   lambda M,t,s: curveball(M,np.random.default_rng(7000+s)),
 'null2':  lambda M,t,s: curveball(M,np.random.default_rng(8000+s)),
 'p_conc': lambda M,t,s: plant_conc(curveball(M,np.random.default_rng(7000+s)),
             np.flatnonzero(np.random.default_rng(900).random(M.shape[0])<0.05),
             np.random.default_rng(901)),
 'p_diff': lambda M,t,s: plant_diff(curveball(M,np.random.default_rng(7000+s)),
             np.random.default_rng(902)),
}
res={}
for name,fn in BUILDS.items():
    D=check_columns(profile(fn,1),'A11R18')
    res[name]={c:matched_quantiles(D,c) for c in ['v_max','v_cv','v_rare']}
    res[name]['_n']=len(D)
    print(f"  {name} done  n={len(D):,}",flush=True)
OUT=pathlib.Path(__file__).parent/'results'
for c in ['v_max','v_cv','v_rare']:
    T=pd.DataFrame({a:res[a][c] for a in BUILDS},index=[f'p{q}' for q in QS])
    print(f"\n=== {c}(在勾选数 x 平均意外度分层内标准化) ===")
    print(T.round(3).to_string())
    T.to_csv(OUT/f'{c}.csv')
def tail(a,c): return res[a][c][QS.index(95)]
def med(a,c):  return res[a][c][QS.index(50)]
c='v_rare'
print(f"\n=== 判别:上尾 vs 中位数,{c} ===")
for a in BUILDS:
    print(f"  {a:8s} p50 {med(a,c):+.3f}   p95 {tail(a,c):+.3f}   尾−中 {tail(a,c)-med(a,c):+.3f}")
dconc=tail('p_conc',c)-med('p_conc',c); ddiff=tail('p_diff',c)-med('p_diff',c)
dnull=tail('null',c)-med('null',c); dreal=tail('real',c)-med('real',c)
g=Gate("弥散的连续维度,还是集中的恋物子群?")
g.require_resolvable_first("真实的尾−中 相对零", effect=dreal-dnull, spread=abs(tail('null2',c)-med('null2',c)-dnull)+1e-4)
g.asserted("两个种植彼此可分(用它们的差,不用常数阈值)",
           abs(dconc-ddiff)>2*abs(dnull-(tail('null2',c)-med('null2',c))),
           f"集中 {dconc:+.3f} vs 弥散 {ddiff:+.3f},差 {abs(dconc-ddiff):.3f};零的抖动 "
           f"{abs(dnull-(tail('null2',c)-med('null2',c))):.3f}")
g.artifact_cannot_explain("第二次独立零", artifact=tail('null2',c)-med('null2',c)-dnull,
                          effect=dreal-dnull, spread=1e-3)
print(); print(g)
print(f"\n  真实 尾−中 {dreal:+.3f}   零 {dnull:+.3f}   集中种植 {dconc:+.3f}   弥散种植 {ddiff:+.3f}")
if g.verdict():
    print(f"  -> 真实更接近{'集中(存在恋物式子群)' if abs(dreal-dconc)<abs(dreal-ddiff) else '弥散(只有连续维度)'}")
print(f"\nartifact sha1 {hashlib.sha1(str(res).encode()).hexdigest()[:12]}")
