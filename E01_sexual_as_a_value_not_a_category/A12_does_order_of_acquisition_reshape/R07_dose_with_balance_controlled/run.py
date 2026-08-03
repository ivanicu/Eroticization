import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R07 -- 剂量,把 #110d 测出来的伪影机制控制掉。

#110c 剂量再次读不出来,而 #110d 把原因直接测了出来:间隔大的对,顺序标签更一边倒
(corr(平均间隔, 类别平衡度) = -0.597;小间隔对 0.379 vs 大间隔对 0.291),AUC 在类别不平衡时
的估计 regime 不同 —— 平种植因此也能产生 +0.0079 的斜率,和剂量无关。

#110e:把类别平衡度放进对间回归。共线性 -0.597 < 0.8,可以同时估。

而且这一轮把判决条件写死在前面:**只有当平种植的假斜率降到观测斜率之下,才允许读剂量。**
上一轮两者是 +0.0079 vs +0.0080 —— 一模一样,所以什么都不能读。

ESTIMAND        per-pair 效应(合成无信号世界偏移)对 (a) 平均间隔 (b) 平均首获年龄 的斜率,
                回归里同时含 log n 和 **类别平衡度**。
KILL            threshold-free。前置条件:平种植的假斜率 < 观测斜率的一半,否则整轮 UNVERIFIED。
POSITIVE CTRL   (1) 种间隔剂量 -> 正斜率;(2) 平种植 -> 零斜率(现在必须真的接近零)
NEGATIVE CTRL   纯置换 y 的效应(#110f 已知有 +0.0079 残差,报出来不掩盖)
IMPOSSIBLE      因果方向;依赖年龄的共因。
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
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
P=(df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float).fillna(0.).values
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values; breadth=P.sum(1)
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
mc=np.nanmean(V[:,sorted(CONCRETE)],axis=1); mr=np.nanmean(V[:,sorted(RELATIONAL)],axis=1)
mc=np.where(np.isfinite(mc),mc,np.nanmean(mc)); mr=np.where(np.isfinite(mr),mr,np.nanmean(mr))
COV=np.c_[male,agev,breadth,prec,mc,mr,mc-mr]; COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
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
def inc(idx,y,seed=1):
    return ridge_auc(np.c_[COV[idx],P[idx]],y,seed)-ridge_auc(COV[idx],y,seed)
def synth_offset(idx,y,ndraw=3,seed=1):
    """#109e 的真值偏移:y 只由协变量生成,P 无额外信息,直接量 (full-base)。"""
    X=np.c_[np.ones(len(idx)),COV[idx]]
    w=np.linalg.lstsq(X,y,rcond=None)[0]; lin=np.clip(X@w,0.02,0.98)
    return np.nanmean([inc(idx,(np.random.default_rng(seed+300+d).random(len(idx))<lin).astype(float),seed)
                       for d in range(ndraw)])
def eff(idx,y,seed=1): return inc(idx,y,seed)-synth_offset(idx,y,seed=seed)
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)]
np.random.default_rng(3).shuffle(pairs)
rows=[]
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)
    rp=np.random.default_rng(88)
    rows.append(dict(a=a,b=b,kind='same' if KIND[a]==KIND[b] else 'cross',n=len(idx),
        gap=float(np.mean(np.abs(V[idx,a]-V[idx,b]))),
        first=float(np.mean(np.minimum(V[idx,a],V[idx,b]))),
        bal=float(min(y.mean(),1-y.mean())),
        e=eff(idx,y), e_null=eff(idx,y[rp.permutation(len(y))])))
    if len(rows)%20==0: print(f"  {len(rows)} pairs",flush=True)
    if len(rows)>=130: break
ctrl=[]
for mode in ['dose','flat']:
    cs=0
    for a,b in pairs:
        if KIND[a]!=KIND[b]: continue
        m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
        if len(idx)<400: continue
        gp=float(np.mean(np.abs(V[idx,a]-V[idx,b])))
        y_=(V[idx,a]<V[idx,b]).astype(float); bl=float(min(y_.mean(),1-y_.mean()))
        rp=np.random.default_rng(31); wv=rp.normal(size=P.shape[1])
        sig=(P[idx]@wv>np.median(P[idx]@wv)).astype(float)
        y0=(V[idx,a]<V[idx,b]).astype(float); y0=y0[rp.permutation(len(y0))]
        g=(0.05+0.045*gp) if mode=='dose' else 0.22
        yp=np.where(rp.random(len(idx))<g,sig,y0)
        ctrl.append(dict(mode=mode,gap=gp,bal=bl,n=len(idx),e=eff(idx,yp,seed=5))); cs+=1
        if cs>=30: break
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def slope(x,y,*ctrls,B=2000,seed=0):
    """#110e:控制变量任意个,类别平衡度必须在其中。"""
    cols=[np.ones(len(x)),(x-x.mean())/(x.std()+1e-9)]
    for c in ctrls: cols.append((c-c.mean())/(c.std()+1e-9))
    X=np.column_stack(cols)
    b=np.linalg.lstsq(X,y,rcond=None)[0][1]
    rb=np.random.default_rng(seed)
    bs=[np.linalg.lstsq(X[i],y[i],rcond=None)[0][1] for i in
        (rb.integers(0,len(x),len(x)) for _ in range(B))]
    return float(b),float(np.std(bs))
print(f"\n=== {len(D)} 对   同类 {(D.kind=='same').sum()}  跨类 {(D.kind=='cross').sum()} ===")
print(D.groupby('kind').agg(n=('n','median'),gap=('gap','mean'),first=('first','mean'),
      e=('e','mean'),se=('e',lambda s:s.std()/np.sqrt(len(s))),
      e_null=('e_null','mean')).round(4).to_string())
S=D[D.kind=='same']; X_=D[D.kind=='cross']
out=[]
for nm,sub in [('same',S),('cross',X_)]:
    for dose in ['gap','first']:
        b,sd=slope(sub[dose].values,sub.e.values,np.log(sub.n.values),sub.bal.values)
        out.append(dict(kind=nm,dose=dose,slope=b,sd=sd))
R=pd.DataFrame(out); print("\n=== 对间剂量斜率(含 log n) ===")
print(R.round(4).to_string(index=False))
cd=C[C['mode']=='dose']; cf=C[C['mode']=='flat']
bd,sdd=slope(cd.gap.values,cd.e.values,np.log(cd.n.values),cd.bal.values)
bf,sdf=slope(cf.gap.values,cf.e.values,np.log(cf.n.values),cf.bal.values)
print(f"\n=== 正对照 ===\n  种剂量 斜率 {bd:+.4f} ± {sdd:.4f}   平种植 斜率 {bf:+.4f} ± {sdf:.4f}")
se=lambda s: s.std()/np.sqrt(len(s))
sg=R[(R.kind=='same')&(R.dose=='gap')].iloc[0]; sf=R[(R.kind=='same')&(R.dose=='first')].iloc[0]
print(f"  (上一轮:假斜率 +0.0079 vs 首获年龄斜率 +0.0080 —— 一模一样,什么都读不了)")
xg=R[(R.kind=='cross')&(R.dose=='gap')].iloc[0]
g=Gate("#109f 剂量,在对的零上")
g.negative_control("纯置换 y 的效应(合成偏移是否有偏)", null=abs(S.e_null.mean()),
                   effect=abs(S.e.mean()))
g.resolvable("种植的对间剂量被测出", effect=bd, spread=sdd)
g.asserted("平种植的假斜率 < 观测斜率的一半 (#110c 的前置条件)",
           abs(bf)<0.5*max(abs(sg.slope),abs(sf.slope)),
           f"假斜率 |{bf:+.4f}| vs 观测斜率一半 {0.5*max(abs(sg.slope),abs(sf.slope)):.4f}")
g.asserted("同类效应本轮为正(功率前提)", S.e.mean()>2*se(S.e),
           f"{S.e.mean():+.4f} vs 2*SE {2*se(S.e):.4f}")
g.resolvable("同类·间隔斜率", effect=sg.slope, spread=sg.sd)
g.resolvable("同类·首获年龄斜率", effect=sf.slope, spread=sf.sd)
print(); print(g)
print(f"\n  伪影检查:同类间隔斜率 {sg.slope:+.4f}   跨类间隔斜率 {xg.slope:+.4f}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
