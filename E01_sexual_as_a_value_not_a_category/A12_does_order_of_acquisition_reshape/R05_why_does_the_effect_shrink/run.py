import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R05 -- 我两轮前提交的效应,在我自己的三次重估里一路缩小。先查这个。

  #107 (R02)  +0.0236    n>=400,reps=5,2 次 *分层* 置换
  #108 (R03)  +0.0121    n>=900,reps=4,1 次 *平* 置换
  #109 (R04)  +0.0058    n>=400,reps=8,3 次 *平* 置换

置换抽样数增加只该降方差,不该系统性移动位置。所以漂移一定来自别处,而唯一的结构性差异是
置换的**种类**:分层置换在协变量分层内打乱,保留协变量对结果的预测力;平置换把它一起毁掉。

  分层置换 -> 基线仍在 0.56 附近,过拟合代价在**正确的回归 regime** 里测
  平置换   -> 基线塌到 0.50,过拟合代价在**另一个 regime** 里测

两者测的不是同一个偏移。哪个是对的偏移,决定 #106/#107 是发现还是伪影。

ESTIMAND        同一批对、同一个 ridge、同一批人,只改 (置换种类 x 置换抽样数 x reps),
                看效应在哪个轴上移动、在哪个轴上收敛。
IDENTIFICATION  完全受控的仪器实验;唯一变的是估计量的三个旋钮。
WORLDS          artifact   效应随抽样数单调趋零 -> #106/#107 是偏移估计不足
                regime     效应随置换种类跳变而对抽样数收敛 -> 两个偏移之一是错的 regime
                stable     两轴都收敛到同一个正值 -> 漂移是我记错了配置
KILL            threshold-free:效应对每个旋钮的曲线,以及它是否收敛。
POSITIVE CTRL   把一个已知强度的信号种进去,同样扫旋钮 —— 真信号在所有配置下都该被测出,
                且收敛到同一个值。不收敛的话是估计量坏了,不是效应假。
NEGATIVE CTRL   纯置换数据扫同样的旋钮 —— 必须处处为零。
IMPOSSIBLE      无。这是纯仪器测量。
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
def strata(idx):
    q=lambda v:np.digitize(v,np.quantile(v,[.33,.66]))
    return (male[idx]>0).astype(int)*9+q(breadth[idx])*3+q(prec[idx])
def eff(idx,y,ptype,ndraw,reps,seed=1):
    base=ridge_auc(COV[idx],y,np.random.default_rng(seed),reps=reps)
    full=ridge_auc(np.c_[COV[idx],P[idx]],y,np.random.default_rng(seed),reps=reps)
    st=strata(idx) if ptype=='strat' else None
    pm=[]
    for d in range(ndraw):
        rp=np.random.default_rng(seed+90+d)
        if ptype=='strat':
            yp=y.copy()
            for s in np.unique(st):
                w=np.flatnonzero(st==s)
                if len(w)>1: yp[w]=y[w][rp.permutation(len(w))]
        else: yp=y[rp.permutation(len(y))]
        pm.append(ridge_auc(np.c_[COV[idx],P[idx]],yp,np.random.default_rng(seed),reps=reps)-
                  ridge_auc(COV[idx],yp,np.random.default_rng(seed),reps=reps))
    return (full-base)-np.nanmean(pm)
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2) if KIND[a]==KIND[b]]
np.random.default_rng(3).shuffle(pairs)
SEL=[]
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)>=400: SEL.append((a,b,idx))
    if len(SEL)>=40: break
print(f"同类对 {len(SEL)}  中位 n {int(np.median([len(i) for _,_,i in SEL]))}",flush=True)
rows=[]
for ptype in ['strat','plain']:
    for nd in [1,2,4,8]:
        for reps in [5,8]:
            es=[];  ps=[];  ns=[]
            for a,b,idx in SEL:
                y=(V[idx,a]<V[idx,b]).astype(float)
                es.append(eff(idx,y,ptype,nd,reps))
                rp=np.random.default_rng(77); wv=rp.normal(size=P.shape[1])
                sig=(P[idx]@wv>np.median(P[idx]@wv)).astype(float)
                y0=y[rp.permutation(len(y))]
                yp=np.where(rp.random(len(idx))<0.20,sig,y0)          # 已知强度的种植
                ps.append(eff(idx,yp,ptype,nd,reps))
                ns.append(eff(idx,y0,ptype,nd,reps))                  # 纯置换,必须处处为零
            rows.append(dict(ptype=ptype,ndraw=nd,reps=reps,
                             real=np.nanmean(es),real_se=np.nanstd(es)/np.sqrt(len(es)),
                             plant=np.nanmean(ps),null=np.nanmean(ns)))
            print(f"  {ptype} nd={nd} reps={reps}  real {rows[-1]['real']:+.4f}  "
                  f"plant {rows[-1]['plant']:+.4f}  null {rows[-1]['null']:+.4f}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print("\n=== 真实效应,按 (置换种类 x 抽样数 x reps) ===")
print(D.pivot_table(index=['ptype','reps'],columns='ndraw',values='real').round(4).to_string())
print("\n=== 已知种植(强度 0.20),同样的旋钮 ===")
print(D.pivot_table(index=['ptype','reps'],columns='ndraw',values='plant').round(4).to_string())
print("\n=== 纯置换数据(必须处处为零) ===")
print(D.pivot_table(index=['ptype','reps'],columns='ndraw',values='null').round(4).to_string())
S=D[D.ptype=='strat']; PL=D[D.ptype=='plain']
g=Gate("效应为什么缩小?")
g.asserted("纯置换数据处处为零", D.null.abs().max()<0.012, f"max |null| = {D.null.abs().max():.4f}")
g.asserted("已知种植在所有配置下都被测出", D.plant.min()>0.02,
           f"min plant = {D.plant.min():.4f} over {len(D)} 配置")
g.asserted("已知种植跨配置稳定(极差 < 自身一半)",
           (D.plant.max()-D.plant.min())<0.5*D.plant.mean(),
           f"极差 {D.plant.max()-D.plant.min():.4f} vs 均值一半 {0.5*D.plant.mean():.4f}")
g.asserted("真实效应随抽样数收敛(nd=4 与 nd=8 之差 < nd=1 与 nd=2 之差)",
           abs(S[S.ndraw==8].real.mean()-S[S.ndraw==4].real.mean())<
           abs(S[S.ndraw==2].real.mean()-S[S.ndraw==1].real.mean())+1e-9,
           f"nd 4->8 {abs(S[S.ndraw==8].real.mean()-S[S.ndraw==4].real.mean()):.4f} vs "
           f"1->2 {abs(S[S.ndraw==2].real.mean()-S[S.ndraw==1].real.mean()):.4f}")
print(); print(g)
print(f"\n  分层置换下收敛值 {S[S.ndraw==8].real.mean():+.4f}   平置换下 {PL[PL.ndraw==8].real.mean():+.4f}")
print(f"  两种偏移之差 {S[S.ndraw==8].real.mean()-PL[PL.ndraw==8].real.mean():+.4f}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
