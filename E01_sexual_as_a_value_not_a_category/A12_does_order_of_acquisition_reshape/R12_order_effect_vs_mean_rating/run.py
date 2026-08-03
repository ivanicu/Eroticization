import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R12 -- #107/#110 的顺序效应,是不是同一条通路?

#115 撤回了 #112a:"早=中心"其实大半是"什么都给高分的人"。#115e 指出 #107/#110 的顺序对设计
可能走同一条通路,而且从没测过:

  结局 y = (起始A < 起始B)。#114a 已量出回忆偏差把心爱的往前拉,
  所以 y 经由 (评分A − 评分B) 依赖评分;而偏好轮廓预测评分。
  #107/#110 的协变量里有广度、早熟度、具体/关系倾向 —— 没有人均评分,也没有这两项各自的评分。

四臂,同一批对、同一个估计量,只改协变量:

  base       #110 原协变量(性别·年龄·广度·早熟·具体均值·关系均值·差)
  +mean      加人均评分
  +pairrat   加这一对两项各自的评分,以及它们的差            <- 直接堵 (评分A − 评分B) 那条通路
  +both      两个都加                                        <- 最严

判定写在前面:+both 存活 -> #107/#110 站住;崩塌 -> 整条顺序线归约为作答水平。

零用 #109e 认定的**合成无信号世界**(不是置换),并在 gate 里命名 (#109c)。

ESTIMAND        同类对的 per-pair 效应,四种协变量集。
KILL            threshold-free;+both 臂对自身 SE。
POSITIVE CTRL   种植的顺序信号必须在 +both 臂仍被测出。
NEGATIVE CTRL   纯置换 y(#110f 已知有小残差,照报)。
IMPOSSIBLE      与 #110 相同。
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools, re
sys.path.insert(0,str(ROOT))
from lib.gates import Gate
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
CONCRETE={1,2,4,13,14,16,17,18,20,26,29,30}; RELATIONAL={3,5,6,8,9,10,11,12,15,21,22,23,24,27}
KIND={**{i:'C' for i in CONCRETE},**{i:'R' for i in RELATIONAL}}
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
meanrating=np.nanmean(np.where(np.isfinite(R),R,np.nan),axis=1)
meanrating=np.where(np.isfinite(meanrating),meanrating,np.nanmean(meanrating))
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
mc=np.nanmean(V[:,sorted(CONCRETE)],axis=1); mr=np.nanmean(V[:,sorted(RELATIONAL)],axis=1)
mc=np.where(np.isfinite(mc),mc,np.nanmean(mc)); mr=np.where(np.isfinite(mr),mr,np.nanmean(mr))
BASE=np.c_[male,agev,breadth,prec,mc,mr,mc-mr]
def zs(X): return (X-X.mean(0))/(X.std(0)+1e-9)
BASE=zs(BASE); MEANR=zs(np.c_[meanrating])
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
print(f"能对上评分列的类别 {len(best)}",flush=True)
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
def eff(idx,y,C,seed=1,ndraw=3):
    inc=lambda yy: ridge_auc(np.c_[C,P[idx]],yy,seed)-ridge_auc(C,yy,seed)
    X=np.c_[np.ones(len(idx)),C]
    w=np.linalg.lstsq(X,y,rcond=None)[0]; lin=np.clip(X@w,0.02,0.98)
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
    PR=zs(np.c_[ra,rb,ra-rb])
    C={'base':BASE[idx],'+mean':np.c_[BASE[idx],MEANR[idx]],
       '+pairrat':np.c_[BASE[idx],PR],'+both':np.c_[BASE[idx],MEANR[idx],PR]}
    rp=np.random.default_rng(88)
    for tag,CC in C.items():
        rows.append(dict(a=a,b=b,arm=tag,n=len(idx),e=eff(idx,y,CC)))
    rows.append(dict(a=a,b=b,arm='null_both',n=len(idx),
                     e=eff(idx,y[rp.permutation(len(y))],C['+both'])))
    w=np.random.default_rng(66).normal(size=P.shape[1]); sig=P[idx]@w
    y0=y[rp.permutation(len(y))]
    yp=np.where(np.random.default_rng(67).random(len(idx))<0.20,
                (sig>np.median(sig)).astype(float),y0)
    rows.append(dict(a=a,b=b,arm='plant_both',n=len(idx),e=eff(idx,yp,C['+both'])))
    if len(rows)%60==0: print(f"  {len(rows)//6} pairs",flush=True)
    if len(rows)>=6*45: break
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
S=D.groupby('arm').agg(k=('e','size'),n=('n','median'),e=('e','mean'),
                       se=('e',lambda s:s.std()/np.sqrt(len(s))))
order=['base','+mean','+pairrat','+both','null_both','plant_both']
print("\n=== 同类对顺序效应,按协变量集 ===")
print(S.reindex(order).round(4).to_string())
g=Gate("#115e 顺序效应是不是同一条通路?")
g.negative_control("纯置换 y(+both 臂)", null=S.loc['null_both','e'], effect=S.loc['+both','e'])
g.positive_control("种植在 +both 臂仍被测出", planted=S.loc['plant_both','e'],
                   floor=S.loc['null_both','e'], spread=S.loc['plant_both','se'])
g.offset_control("+both 臂的效应", effect=S.loc['+both','e'], offset=0.0,
                 spread=S.loc['+both','se'], null_kind='合成无信号世界(已在 eff 内减掉)')
print(); print(g)
kr=S.loc['+both','e']/max(S.loc['base','e'],1e-9)
print(f"\n  base {S.loc['base','e']:+.4f} -> +人均评分 {S.loc['+mean','e']:+.4f}"
      f" -> +对内两项评分 {S.loc['+pairrat','e']:+.4f} -> +both {S.loc['+both','e']:+.4f}")
print(f"  +both 保留 {100*kr:.0f}%")
print(f"  -> {'#107/#110 站住,顺序效应不是作答水平' if kr>0.5 and S.loc['+both','e']>2*S.loc['+both','se'] else '同一条通路,顺序线也要撤回或改写'}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
