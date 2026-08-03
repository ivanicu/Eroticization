import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R04 -- 剂量,这次不牺牲功率。

#108 证明对内分层把每格砍到 ~430 人,而在那个 n 上效应保留 -6% —— 剂量格子里的 ~0 是失明,不是零。
#108e 给了保功率的修法:不在对内分,在对间做。每个对贡献一个用全样本(~1300)算出的效应,
和它自己的典型间隔、典型首获年龄;然后把效应回归到剂量上。n 变成对数,而每个效应保留全样本功率。

同时修 #108d:每对效应用更多重复和更多置换抽样,把跨轮 2 倍的估计量方差压下去。

ESTIMAND        同类对内,per-pair 效应 对 (a) 该对的平均间隔 (b) 该对的平均首获年龄 的回归斜率。
IDENTIFICATION  每个效应用全样本,剂量是对的属性,不是人的属性 -> 不牺牲任何 n。
CONFOUNDS(跑之前):
                (a) 间隔大的对更可能跨类 -> 只在同类臂内回归,跨类臂单独报作对照
                (b) 间隔大 -> 顺序标签测得更可靠 -> 伪影。若为伪影,同类与跨类斜率应相同
                (c) 对的样本量与间隔相关 -> 把 log n 放进回归
WORLDS          重塑  同类臂斜率为正且大于跨类臂
                共因  斜率为零
                伪影  两臂斜率相同且为正
KILL            threshold-free;斜率对自身自助 SE。
POSITIVE CTRL   两个:(1) 种一个强度随对间隔递增的效应 -> 必须测出正斜率
                     (2) 种一个跟间隔无关的平效应     -> 必须测出零斜率
NEGATIVE CTRL   零已在 eff() 内部作为偏移减掉 (#106c)。
IMPOSSIBLE      因果方向;以及一个本身依赖年龄的共因。
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
def ridge_auc(X,y,rng,alpha=50.,reps=8):
    out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<10 or (1-y[tr]).sum()<10: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def eff(idx,y,seed=1,ndraw=3):
    base=ridge_auc(COV[idx],y,np.random.default_rng(seed))
    full=ridge_auc(np.c_[COV[idx],P[idx]],y,np.random.default_rng(seed))
    pm=[]
    for d in range(ndraw):
        rp=np.random.default_rng(seed+70+d); yp=y[rp.permutation(len(y))]
        pm.append(ridge_auc(np.c_[COV[idx],P[idx]],yp,np.random.default_rng(seed))-
                  ridge_auc(COV[idx],yp,np.random.default_rng(seed)))
    return (full-base)-np.nanmean(pm)
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)]
np.random.default_rng(3).shuffle(pairs)
rows=[]
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)
    rows.append(dict(a=a,b=b,kind='same' if KIND[a]==KIND[b] else 'cross',n=len(idx),
                     gap=float(np.mean(np.abs(V[idx,a]-V[idx,b]))),
                     first=float(np.mean(np.minimum(V[idx,a],V[idx,b]))),
                     bal=float(min(y.mean(),1-y.mean())),e=eff(idx,y)))
    if len(rows)%25==0: print(f"  {len(rows)} pairs",flush=True)
    if len(rows)>=150: break
# 正对照:对间剂量种植 / 平种植
ctrl=[]
for mode in ['dose','flat']:
    cs=[]
    for a,b in pairs:
        if KIND[a]!=KIND[b]: continue
        m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
        if len(idx)<400: continue
        gp=float(np.mean(np.abs(V[idx,a]-V[idx,b])))
        rp=np.random.default_rng(31); wv=rp.normal(size=P.shape[1])
        sig=(P[idx]@wv>np.median(P[idx]@wv)).astype(float)
        y0=(V[idx,a]<V[idx,b]).astype(float); y0=y0[rp.permutation(len(y0))]
        g=(0.05+0.045*gp) if mode=='dose' else 0.22
        yp=np.where(rp.random(len(idx))<g,sig,y0)
        cs.append(dict(mode=mode,gap=gp,e=eff(idx,yp,seed=5)))
        if len(cs)>=30: break
    ctrl+=cs
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def slope(x,y,z=None,B=2000,seed=0):
    X=np.c_[np.ones(len(x)),(x-x.mean())/(x.std()+1e-9)]
    if z is not None: X=np.c_[X,(z-z.mean())/(z.std()+1e-9)]
    b=np.linalg.lstsq(X,y,rcond=None)[0][1]
    rb=np.random.default_rng(seed)
    bs=[np.linalg.lstsq(X[i],y[i],rcond=None)[0][1] for i in
        (rb.integers(0,len(x),len(x)) for _ in range(B))]
    return float(b),float(np.std(bs))
print(f"\n=== {len(D)} 对   同类 {(D.kind=='same').sum()}  跨类 {(D.kind=='cross').sum()} ===")
print(D.groupby('kind').agg(n=('n','median'),gap=('gap','mean'),first=('first','mean'),
                            e=('e','mean'),se=('e',lambda s:s.std()/np.sqrt(len(s)))).round(4).to_string())
S=D[D.kind=='same']; X=D[D.kind=='cross']
out=[]
for nm,sub in [('same',S),('cross',X)]:
    for dose in ['gap','first']:
        b,sd=slope(sub[dose].values,sub.e.values,np.log(sub.n.values))
        out.append(dict(kind=nm,dose=dose,slope=b,sd=sd,ratio=b/max(sd,1e-9)))
R=pd.DataFrame(out); print("\n=== 对间剂量斜率(已含 log n 作协变量) ===")
print(R.round(4).to_string(index=False))
cd=C[C['mode']=='dose']; cf=C[C['mode']=='flat']
bd,sdd=slope(cd.gap.values,cd.e.values); bf,sdf=slope(cf.gap.values,cf.e.values)
print(f"\n=== 正对照 ===\n  种剂量: 斜率 {bd:+.4f} ± {sdd:.4f} ({bd/sdd:.1f}x)"
      f"   平种植: 斜率 {bf:+.4f} ± {sdf:.4f} ({bf/sdf:.1f}x)")
g=Gate("#108e 对间剂量:顺序效应随间隔/早发而增强吗?")
g.resolvable("种植的对间剂量被测出", effect=bd, spread=sdd)
g.asserted("平种植不产生假斜率", abs(bf)<2*sdf, f"|{bf:+.4f}| < {2*sdf:.4f}")
g.asserted("同类效应在本轮仍为正(功率前提)",
           S.e.mean()>2*S.e.std()/np.sqrt(len(S)),
           f"{S.e.mean():+.4f} vs 2*SE {2*S.e.std()/np.sqrt(len(S)):.4f}")
sg=R[(R.kind=='same')&(R.dose=='gap')].iloc[0]
sf=R[(R.kind=='same')&(R.dose=='first')].iloc[0]
g.resolvable("同类·间隔斜率", effect=sg.slope, spread=sg.sd)
g.resolvable("同类·首获年龄斜率", effect=sf.slope, spread=sf.sd)
print(); print(g)
xg=R[(R.kind=='cross')&(R.dose=='gap')].iloc[0]
print(f"\n  伪影检查:同类间隔斜率 {sg.slope:+.4f}   跨类间隔斜率 {xg.slope:+.4f}")
print(f"  (可靠性伪影两臂应相同;重塑只应同类臂更强)")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
