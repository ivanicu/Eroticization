import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R22 -- 把 #126c 的守卫回溯到两个现存声明上。

#126c 是本项目第一个**在设计时被漏掉**的混淆(前四十一个都是写下了但写错了)。
守卫 check_disjoint_items 已落地。#126 的 NEXT:回头查同样组合 item 派生量的两个声明。

审计一 · #100 的信度(+0.432)
  分半是按**不相交的块**分的,这一层是干净的。但残差化用的 `picks` 是**全部块**的勾选总数,
  两半的 item 都在里面 -> 两个半分被同一个含有彼此的量残差化。
  修法:A 半用 A 半自己的勾选数残差化,B 半用 B 半自己的。

审计二 · #116 的顺序效应(+0.0159)
  预测器是完整的 68 列偏好矩阵 P,而结局 y = (起始A < 起始B)。
  **该对自己的两个评分列就在 P 里**,而 #114 已证明评分经由回忆偏差扭曲起始年龄。
  #116 把它们作为**协变量**加了进去(控制线性效应),但它们**仍在预测器里**,模型可以非线性地用。
  修法:把 best[a]、best[b] 两列从 P 里剔掉。

两个审计都只改一件事,其余完全照抄原轮,以便可比。

KILL            threshold-free;若某个声明在剔除共享 item 后掉到不可分辨,该声明必须改写。
POSITIVE CTRL   #100 用种植特质;#116 用种植顺序信号。都必须在剔除后仍被测出。
NEGATIVE CTRL   各自原轮的零。
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools, re
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_disjoint_items, check_columns
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
# ---------------- 审计一:#100 的信度 ----------------
print("=== 审计一 · #100 的分半信度 ===",flush=True)
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
def plant_trait(M,aff,rng):
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in np.flatnonzero(aff>0):
        d=0
        for _ in range(8*int(aff[i])):
            if d>=aff[i]: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; d+=1
    return Mw
sb=lambda r: 2*r/(1+r) if r>-1 else np.nan
def reliability(build,mode,seed=1,reps=6):
    blk={}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        Mw=build(M,seed); ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        tot=Mw@ref; k=Mw.sum(1)
        for j,gi in enumerate(idx):
            if k[j]>0: blk.setdefault(gi,[]).append((bi,tot[j],k[j]))
    ks=[i for i,v in blk.items() if len(v)>=6]
    out=[]
    for rep in range(reps):
        rr=np.random.default_rng(700+rep); a=[];b=[];ka=[];kb=[]
        for i in ks:
            v=blk[i]; o=rr.permutation(len(v)); h=len(v)//2
            A=[v[j] for j in o[:h]]; B=[v[j] for j in o[h:2*h]]
            sa=sum(x[2] for x in A); sbk=sum(x[2] for x in B)
            a.append(sum(x[1] for x in A)/max(sa,1)); b.append(sum(x[1] for x in B)/max(sbk,1))
            ka.append(sa); kb.append(sbk)
        a=np.array(a);b=np.array(b);ka=np.array(ka);kb=np.array(kb)
        if mode=='shared':            # #100 原样:两半都用**全部**勾选数
            tot_k=ka+kb
            X=np.c_[np.ones(len(tot_k)),tot_k,np.log(np.maximum(tot_k,1))]
            ra=a-X@np.linalg.lstsq(X,a,rcond=None)[0]; rb=b-X@np.linalg.lstsq(X,b,rcond=None)[0]
        else:                          # 修正:各半用**自己**的勾选数
            Xa=np.c_[np.ones(len(ka)),ka,np.log(np.maximum(ka,1))]
            Xb=np.c_[np.ones(len(kb)),kb,np.log(np.maximum(kb,1))]
            ra=a-Xa@np.linalg.lstsq(Xa,a,rcond=None)[0]; rb=b-Xb@np.linalg.lstsq(Xb,b,rcond=None)[0]
        out.append(np.corrcoef(ra,rb)[0,1])
    return sb(float(np.mean(out))),len(ks)
rgc=np.random.default_rng(6101); aff=np.floor(rgc.exponential(0.8,size=len(ALLP))).astype(float)
B={'real':lambda M,s: M,
   'null':lambda M,s: curveball(M,np.random.default_rng(6200+s)),
   'plant':lambda M,s: plant_trait(curveball(M,np.random.default_rng(6200+s)),
                                   aff[[PM[p] for p in RAW[list(RAW)[0]]['ppl']]][:M.shape[0]]
                                   if False else aff[:M.shape[0]],np.random.default_rng(6300+s))}
res100={}
for mode in ['shared','own']:
    for w,f in B.items():
        r,n=reliability(f,mode); res100[(mode,w)]=r
        print(f"  {mode:6s} {w:6s} 残差化信度 {r:+.4f}   n={n:,}",flush=True)
# ---------------- 审计二:#116 的顺序效应 ----------------
print("\n=== 审计二 · #116 的顺序效应 ===",flush=True)
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
CONCRETE={1,2,4,13,14,16,17,18,20,26,29,30}; RELATIONAL={3,5,6,8,9,10,11,12,15,21,22,23,24,27}
KIND={**{i:'C' for i in CONCRETE},**{i:'R' for i in RELATIONAL}}
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
Pm=(np.nan_to_num(R)>0).astype(float); breadth=Pm.sum(1)
meanrating=np.nanmean(np.where(np.isfinite(R),R,np.nan),axis=1)
meanrating=np.where(np.isfinite(meanrating),meanrating,np.nanmean(meanrating))
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
mc=np.nanmean(V[:,sorted(CONCRETE)],axis=1); mr=np.nanmean(V[:,sorted(RELATIONAL)],axis=1)
mc=np.where(np.isfinite(mc),mc,np.nanmean(mc)); mr=np.where(np.isfinite(mr),mr,np.nanmean(mr))
zs=lambda X:(X-X.mean(0))/(X.std(0)+1e-9)
BASE=zs(np.c_[male,agev,breadth,prec,mc,mr,mc-mr,meanrating])
def norm(s): return re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<10 or n0<10: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,seed,alpha=50.,reps=8):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<10 or (1-y[tr]).sum()<10: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def eff(idx,y,C,Pred,seed=1,ndraw=3):
    inc=lambda yy: ridge_auc(np.c_[C,Pred],yy,seed)-ridge_auc(C,yy,seed)
    X=np.c_[np.ones(len(idx)),C]; w=np.linalg.lstsq(X,y,rcond=None)[0]; lin=np.clip(X@w,0.02,0.98)
    off=np.nanmean([inc((np.random.default_rng(seed+300+d).random(len(idx))<lin).astype(float))
                    for d in range(ndraw)])
    return inc(y)-off
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)
       if KIND[a]==KIND[b] and a in best and b in best]
np.random.default_rng(3).shuffle(pairs)
rows=[]
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)
    ra=np.nan_to_num(R[idx,best[a]]); rb=np.nan_to_num(R[idx,best[b]])
    C=np.c_[BASE[idx],zs(np.c_[ra,rb,ra-rb])]
    keepcols=[i for i in range(Pm.shape[1]) if i not in (best[a],best[b])]
    check_disjoint_items(keepcols,[best[a],best[b]],'A11R22 #116')     # 机械断言,不是散文
    rows.append(dict(a=a,b=b,e_with=eff(idx,y,C,Pm[idx]),
                     e_without=eff(idx,y,C,Pm[idx][:,keepcols]),n=len(idx)))
    if len(rows)>=45: break
D=check_columns(pd.DataFrame(rows),'A11R22')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
print(f"  {len(D)} 对   含该对评分列 {D.e_with.mean():+.4f} ± {se(D.e_with.values):.4f}   "
      f"剔除后 {D.e_without.mean():+.4f} ± {se(D.e_without.values):.4f}   "
      f"保留 {100*D.e_without.mean()/D.e_with.mean():.0f}%")
print("\n=== 判定 ===")
g=Gate("#100 的信度:共享的勾选数协变量有没有制造它?")
g.asserted("修正后零仍与零无法区分", abs(res100[('own','null')])<0.10,
           f"零 {res100[('own','null')]:+.4f}")
g.positive_control("种植特质在修正后仍被测出", planted=res100[('own','plant')],
                   floor=res100[('own','null')], spread=0.02)
g.artifact_cannot_explain("共享协变量的影响",
                          artifact=res100[('shared','real')]-res100[('own','real')],
                          effect=res100[('own','real')], spread=0.02)
print(g)
g2=Gate("#116 的顺序效应:该对自己的评分列有没有在做工?")
g2.require_resolvable_first("剔除后的效应", effect=D.e_without.mean(), spread=se(D.e_without.values))
g2.artifact_cannot_explain("剔除带来的变化", artifact=D.e_with.mean()-D.e_without.mean(),
                           effect=D.e_without.mean(), spread=se(D.e_without.values))
print(); print(g2)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
