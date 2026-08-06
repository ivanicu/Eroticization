import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A20 R02 -- 「它唯一挂得住的外部锚是性别」这句话的力量,全部来自零

#153e:`#101`/`#102` 的零是 `curveball(M, np.random.default_rng(8100))` —— **一次**保边际实现,
而且**每个块都用同一个种子**。

这条比 `#114` 更值得查,原因有二:
① 它支撑 README 上「**唯一挂得住的外部锚是性别 +0.093;人格五因素全部 |r| ≤ 0.056**」,
   而**那句话的力量全部来自零** —— 「五因素全都很小」是一个**否定**,
   而否定的强度完全取决于零有多干净;
② curveball 是**保边际**的,实现方差比一次简单打乱更难直觉,
   而 `#105c` 已证明保边际零对某些统计量**结构上盲**。

ESTIMAND        零臂(curveball 亲和)与 8 个外部变量的相关,在 **10 次独立重抽**下的分布;
                特别是 max|r| 在人格五因素上的分布。
IDENTIFICATION  与 `A11/R148` **逐字相同**的构造,只把 `default_rng(8100)` 换成 10 个种子,
                并且**每个块用不同的种子**(原轮每块同种子,这本身是一个要记录的差异)。
SCOPE           A11/R148 的原口径(>=6 个已识别块的人)。
WORLDS          CLEAN  零的相关全部远小于 0.056 -> 那个否定是硬的,声明不动
                SOFT   零的相关能达到 0.05 上下 -> 「五因素全都 ≤0.056」与零分不开,
                       「唯一的外部锚」这句话要改成「唯一**能与零区分**的外部锚」
KILL            条件式:真实臂必须复现原轮的 +0.093 与 ≤0.056,才读零的重抽。
POSITIVE CTRL   真实臂复现。
NEGATIVE CTRL   —— 本轮**就是**在造零对照。
NOISE FLOOR     10 次重抽(curveball 每次 23 块 x 5n 次交换,比 `#114` 贵一个量级)。
MULTIPLICITY    8 个外部变量 x 10 抽,整格发表。
IMPOSSIBLE      重抽只检验**实现方差**;若 curveball 对这个统计量结构上盲(#105c),
                重抽看不见 —— 那需要一个种植对照,本轮不做,写进 NEXT。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
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

def affinity(seed=None):
    """seed=None -> 真实臂;否则 curveball 零臂(每块一个不同种子)。"""
    tot=np.zeros(len(ALLP)); cnt=np.zeros(len(ALLP)); nblk=np.zeros(len(ALLP))
    cntr=np.zeros(len(ALLP))
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        Mw=M if seed is None else curveball(M,np.random.default_rng(seed*1000+bi))
        tot[idx]+=Mw@ref; cnt[idx]+=Mw.sum(1); cntr[idx]+=M.sum(1); nblk[idx]+=1
    ok=nblk>=6
    Sv=tot[ok]/np.maximum(cnt[ok],1); picks=cntr[ok]
    X=np.c_[np.ones(ok.sum()),picks,np.log(np.maximum(picks,1))]
    return Sv-X@np.linalg.lstsq(X,Sv,rcond=None)[0], ok, picks

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
onset=pd.DataFrame({c:df[c].map(BIN) for c in ons}).mean(axis=1)
EXT={'biomale':pd.to_numeric(df.get('biomale'),errors='coerce'),'age':df['age'].map(AGEMAP),
     'openness':pd.to_numeric(df.get('opennessvariable'),errors='coerce'),
     'neuroticism':pd.to_numeric(df.get('neuroticismvariable'),errors='coerce'),
     'extroversion':pd.to_numeric(df.get('extroversionvariable'),errors='coerce'),
     'conscientious':pd.to_numeric(df.get('consciensiousnessvariable'),errors='coerce'),
     'agreeable':pd.to_numeric(df.get('agreeablenessvariable'),errors='coerce'),
     'powerlessness':pd.to_numeric(df.get('powerlessnessvariable'),errors='coerce')}
PERS=['openness','neuroticism','extroversion','conscientious','agreeable','powerlessness']

AR,ok,picks=affinity(None); ids=ALLP[ok]
def cors(a):
    out={}
    for nm,v in EXT.items():
        y=pd.to_numeric(v,errors='coerce').reindex(ids).values; m=np.isfinite(y)
        if m.sum()<500: continue
        out[nm]=float(np.corrcoef(a[m],y[m])[0,1])
    return out
real=cors(AR)
print(f"\n真实臂(复现 A11/R148):" + "  ".join(f"{k} {v:+.4f}" for k,v in real.items()),flush=True)
print(f"  人格五因素 max|r| = {max(abs(real[p]) for p in PERS if p in real):.4f}"
      f"(#101/#102 报 <=0.056);性别 {real['biomale']:+.4f}(报 +0.093)",flush=True)

NDRAW=10
rows=[]
for s in range(1,NDRAW+1):
    AN,_,_=affinity(s); c=cors(AN)
    rows.append({'draw':s,**c}); print(f"  抽 {s}: 人格 max|r| "
        f"{max(abs(c[p]) for p in PERS if p in c):.4f}  性别 {c['biomale']:+.4f}",flush=True)
T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'redraw.csv',index=False)
pm=np.array([max(abs(r[p]) for p in PERS if p in r) for r in rows])
gm=np.abs(T['biomale'].values)
print(f"\n零臂 10 抽:人格 max|r| 均值 {pm.mean():.4f}  sd {pm.std():.4f}  范围 [{pm.min():.4f}, {pm.max():.4f}]")
print(f"          性别 |r| 均值 {gm.mean():.4f}  范围 [{gm.min():.4f}, {gm.max():.4f}]")
real_pm=max(abs(real[p]) for p in PERS if p in real)
g=Gate('「唯一的外部锚是性别」这句话的力量,零撑得住吗')
g.asserted('真实臂复现原轮(正对照)',abs(real['biomale']-0.093)<0.02 and real_pm<0.08,
           f"性别 {real['biomale']:+.4f}(报 +0.093);人格 max|r| {real_pm:.4f}(报 ≤0.056)")
g.require_resolvable_first('真实的人格 max|r| 是否离得开零臂的分布',
                           abs(real_pm-pm.mean()),float(pm.std()))
g.offset_control('真实人格 max|r| vs 零臂 10 抽',float(real_pm),float(pm.mean()),float(pm.std()),
                 null_kind='curveball 保边际零的 10 次独立实现(每块不同种子)')
g.require_resolvable_first('性别是否离得开零臂',abs(abs(real['biomale'])-gm.mean()),float(gm.std()),
                           family='gender')
g.offset_control('性别 |r| vs 零臂 10 抽',float(abs(real['biomale'])),float(gm.mean()),float(gm.std()),
                 null_kind='同上')
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
