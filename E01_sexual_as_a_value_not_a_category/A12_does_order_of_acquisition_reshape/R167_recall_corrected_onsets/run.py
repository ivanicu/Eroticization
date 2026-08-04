import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R11 -- 用 #114 量出来的回忆偏差,反向校正起始年龄,再重跑 #112。

#114a:心爱的兴趣被报告得更早,斜率 -0.2000 年/评分标准差(19.8x SE),约合 0.74 年。
#114c:这威胁 #112a,通路比直觉的更微妙 —— #112 的"早"是 rel = 本项 − 其余项均值,
       格内本项的拉前量恒定,但**其余项均值的拉前量不恒定**(心爱之物多的人其余项整体被拉前),
       于是"早"与"你有多少心爱之物"挂钩,而后者与轮廓相关。
#114e:畸变的大小和形状已知,可以直接反着加回去。

  校正:  onset_corrected = onset + 0.2 * z(该项评分)
          (斜率是负的,所以加回 +0.2*z 就是把心爱的推回它本来的位置)
  同时:  把**人均评分**放进协变量 —— 这正是 114c 那条通路走的量,而 #112 的协变量里只有广度。

四个臂,同一批人同一个统计量,只改这两件事:

  raw        限定到能校正的类别集,未校正,原协变量           <- 与 #112 可比的基线
  corr       校正起始年龄
  meanrating 未校正,但人均评分进协变量
  both       两个都做                                        <- 最严的一臂

判定写在前面:both 臂存活 -> 回忆偏差被排除,#112a 站住;both 臂崩塌 -> #112a 撤回。

ESTIMAND        早-A 组减晚-A 组的留出 R2 差,四臂。
KILL            threshold-free;both 臂对自身自助 SE,零的种类必须命名。
POSITIVE CTRL   种植的重塑必须在 both 臂下仍被测出(否则 both 臂的零不可读)。
NEGATIVE CTRL   早/晚标签在匹配格内打乱,both 臂上跑。
IMPOSSIBLE      校正用的是**人群平均**的畸变斜率;若畸变强度因人而异,残差仍在。
                以及 #114 的 IMPOSSIBLE 原样继承:拉前 vs 早获得后被珍视,分不开。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
sys.path.insert(0,str(ROOT))
from lib.gates import Gate
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V0=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
meanrating=np.nanmean(np.where(np.isfinite(R),R,np.nan),axis=1)
meanrating=np.where(np.isfinite(meanrating),meanrating,np.nanmean(meanrating))
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
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
CATS=sorted(best); RIS=[best[j] for j in CATS]
V=V0[:,CATS]                                   # 只用能校正的类别集,raw 与 corr 因此完全可比
BETA=0.2                                        # #114a 量出来的畸变斜率,单位 年/评分标准差
Z=np.zeros_like(V)
for k,ri in enumerate(RIS):
    r=R[:,ri]; mu=np.nanmean(r); sd=np.nanstd(r)
    Z[:,k]=np.where(np.isfinite(r),(r-mu)/(sd+1e-9),0.)
VC=V+BETA*Z                                     # 校正:把心爱的推回去
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
COV_A=np.c_[male,agev,breadth,prec]
COV_B=np.c_[male,agev,breadth,prec,meanrating]  # 114c 的通路
COV_A=(COV_A-COV_A.mean(0))/(COV_A.std(0)+1e-9)
COV_B=(COV_B-COV_B.mean(0))/(COV_B.std(0)+1e-9)
print(f"可校正类别 {len(CATS)}   校正量 {BETA} 年/评分SD   "
      f"corr(V, VC) = {np.corrcoef(V[np.isfinite(V)],VC[np.isfinite(V)])[0,1]:.4f}",flush=True)
def r2(y,p):
    v=np.var(y); return 1-np.mean((y-p)**2)/v if v>0 else np.nan
def ridge_r2(X,y,seed,alpha=50.,reps=6):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(r2(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out)
def run_cat(k,Vm,COV,plant=0.,shuffle=False,seed=1):
    j=CATS[k]; ri=RIS[k]
    own=Vm[:,k]; others=np.delete(Vm,k,axis=1)
    rel=own-np.nanmean(others,axis=1)
    y=R[:,ri].copy()
    m=np.isfinite(own)&np.isfinite(rel)&np.isfinite(y)&(np.isfinite(others).sum(1)>=5)
    idx=np.flatnonzero(m)
    if len(idx)<600: return None
    early=(rel[idx]<np.median(rel[idx]))
    if shuffle: early=np.random.default_rng(seed+400).permutation(early)
    Xp=np.delete(P[idx],ri,axis=1)
    if plant>0:
        w=np.random.default_rng(55).normal(size=Xp.shape[1]); sh=Xp@w
        sh=(sh-sh.mean())/(sh.std()+1e-9)
        y=y.copy(); y[idx]=np.where(early,y[idx]+plant*sh*np.nanstd(y[idx]),y[idx])
    yy=y[idx]; Xf=np.c_[Xp,COV[idx]]
    grid=np.digitize(yy,[0.5,1.5,2.5,3.5,4.5])
    keep=[]
    for gq in np.unique(grid):
        w=np.flatnonzero(grid==gq); e=w[early[w]]; l=w[~early[w]]
        kk=min(len(e),len(l))
        if kk<25: continue
        rp=np.random.default_rng(seed+7)
        keep+=list(rp.permutation(e)[:kk])+list(rp.permutation(l)[:kk])
    if len(keep)<400: return None
    keep=np.array(keep); ee=early[keep]
    dy=float(yy[keep][ee].mean()-yy[keep][~ee].mean())
    if abs(dy)>0.05: return None
    a=ridge_r2(Xf[keep][ee],yy[keep][ee],seed); b=ridge_r2(Xf[keep][~ee],yy[keep][~ee],seed)
    return dict(cat=j,n=len(keep),diff=a-b,early=a,late=b,dy=dy)
ARMS=[('raw',V,COV_A,0.,False),('corr',VC,COV_A,0.,False),
      ('meanrating',V,COV_B,0.,False),('both',VC,COV_B,0.,False),
      ('shuf_both',VC,COV_B,0.,True),('plant_both',VC,COV_B,.25,False)]
rows=[]
for k in range(len(CATS)):
    for tag,Vm,CV,pl,sh in ARMS:
        r=run_cat(k,Vm,CV,plant=pl,shuffle=sh)
        if r: rows.append(dict(arm=tag,**r))
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
S=D.groupby('arm').agg(k=('diff','size'),n=('n','median'),diff=('diff','mean'),
    se=('diff',lambda s:s.std()/np.sqrt(len(s))),early=('early','mean'),late=('late','mean'),dy=('dy','mean'))
print("\n=== 早-A 组 减 晚-A 组,四臂 ===")
print(S.reindex([a[0] for a in ARMS]).round(4).to_string())
g=Gate("#114e 回忆偏差校正后,#112a 还在吗?")
g.asserted("评分匹配在 both 臂成立", abs(S.loc['both','dy'])<0.05, f"|Δ评分| {abs(S.loc['both','dy']):.4f}")
g.negative_control("打乱标签(both 臂)", null=S.loc['shuf_both','diff'], effect=S.loc['both','diff'])
g.positive_control("种植在 both 臂仍被测出", planted=S.loc['plant_both','diff'],
                   floor=S.loc['shuf_both','diff'], spread=S.loc['plant_both','se'])
g.resolvable("both 臂的效应", effect=S.loc['both','diff'], spread=S.loc['both','se'])
print(); print(g)
kr=S.loc['both','diff']/max(S.loc['raw','diff'],1e-9)
print(f"\n  raw {S.loc['raw','diff']:+.4f} -> 仅校正起始 {S.loc['corr','diff']:+.4f}"
      f" -> 仅加人均评分 {S.loc['meanrating','diff']:+.4f} -> both {S.loc['both','diff']:+.4f}")
print(f"  both 保留 {100*kr:.0f}%")
print(f"  -> {'#112a 站住:回忆偏差不解释它' if kr>0.5 and S.loc['both','diff']>2*S.loc['both','se'] else '#112a 必须撤回或改写'}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
