import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R08 -- #107 自己写下的那个 NEXT,我绕了三轮才回来。

问题一句话:一个兴趣**来得早**的时候,它是否与这个人其余的偏好**结合得更紧**?

  重塑:  A 在早期到来 -> 它重塑了后续读出的空间 -> 早-A 的人其余轮廓更"A 形"
          -> 用其余偏好预测 A 的评分,在早-A 组里更准
  纯读出:A 只是被读出来的一个值,什么时候读出来不改变空间 -> 两组预测精度相同

这个设计里**没有间隔,没有类别平衡度** —— #110d/#111a 的伪影结构上不存在。

⚠ 跑之前写下的最强混淆:**早-A 的人可能就是更喜欢 A**(兴趣强 -> 来得早)。
  评分是直接观测的,所以用 **A 的评分分层匹配**即可,而且这一步必须在代码里断言,不能写在散文里 (#96a)。
⚠ 第二个:早-A 的人可能整体早熟。早熟度进协变量。
⚠ 第三个:早-A 组和晚-A 组的**广度**可能不同 -> 广度进协变量。

ESTIMAND        用其余偏好预测 A 的评分的留出 R2,早-A 组减晚-A 组,在 A 的评分分层内匹配之后。
IDENTIFICATION  两组在 A 的评分上按格匹配;A 本身及其起始年龄不进预测器。
WORLDS          重塑    早-A 组更高,且在匹配后存活
        读出    两组相同
        强度    差异被 A 的评分匹配吸收掉 -> 是"喜欢得更强",不是"来得更早"
KILL            threshold-free;差值对自身自助 SE,零的种类必须命名 (#109c)。
POSITIVE CTRL   造一个世界:早-A 的人其余轮廓真的更 A 形(按已知强度)-> 必须测出,且单调。
NEGATIVE CTRL   把"早/晚"标签在匹配格内打乱 -> 差值必须为零。
IMPOSSIBLE      因果方向;以及"来得早"与"记得来得早"分不开。
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values; breadth=P.sum(1)
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
COV=np.c_[male,agev,breadth,prec]; COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
# 把每个 onset 类别对到最相关的评分列(只用名字里的关键词,失败就跳过该类别)
import re
def norm(s): return re.sub(r'[^a-z]',' ',s.lower())
key={}
for j,c in enumerate(ons):
    m=re.search(r'interest in ([a-z /-]+)',norm(c)) or re.search(r'interest in ([a-z /-]+)',norm(c))
    if m: key[j]=set(w for w in m.group(1).split() if len(w)>4)
best={}
for j,ws in key.items():
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
print(f"onset 类别 {len(ons)}  能对上评分列的 {len(best)}",flush=True)
def r2(y,p): 
    v=np.var(y)
    return 1-np.mean((y-p)**2)/v if v>0 else np.nan
def ridge_r2(X,y,seed,alpha=50.,reps=6):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(r2(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out)
def run_cat(j,ri,plant=0.0,shuffle=False,seed=1):
    """早-A 组 vs 晚-A 组,用其余偏好预测 A 的评分的留出 R2 之差。"""
    own=V[:,j]; others=np.delete(V,j,axis=1)
    rel=own-np.nanmean(others,axis=1)                    # A 在这个人自己序列里的相对位置
    y=R[:,ri]
    m=np.isfinite(own)&np.isfinite(rel)&np.isfinite(y)&(np.isfinite(others).sum(1)>=5)
    idx=np.flatnonzero(m)
    if len(idx)<600: return None
    early=(rel[idx]<np.median(rel[idx]))
    if shuffle: early=np.random.default_rng(seed+400).permutation(early)
    Xp=np.delete(P[idx],ri,axis=1)                        # A 的评分列本身剔除
    if plant>0:                                           # 早组的轮廓真的更 A 形
        w=np.random.default_rng(55).normal(size=Xp.shape[1])
        sh=(Xp@w); sh=(sh-sh.mean())/(sh.std()+1e-9)
        y=y.copy(); y[idx]=np.where(early,y[idx]+plant*sh*np.nanstd(y[idx]),y[idx])
    yy=y[idx]; Xf=np.c_[Xp,COV[idx]]
    # #96a:按 A 的评分分层匹配,断言在代码里
    grid=np.digitize(yy,[0.5,1.5,2.5,3.5,4.5])
    keep=[]
    for gq in np.unique(grid):
        w=np.flatnonzero(grid==gq); e=w[early[w]]; l=w[~early[w]]
        k=min(len(e),len(l))
        if k<25: continue
        rp=np.random.default_rng(seed+7)
        keep+=list(rp.permutation(e)[:k])+list(rp.permutation(l)[:k])
    if len(keep)<400: return None
    keep=np.array(keep); ee=early[keep]
    matched=abs(yy[keep][ee].mean()-yy[keep][~ee].mean())<0.05
    if not matched: return None
    a=ridge_r2(Xf[keep][ee],yy[keep][ee],seed); b=ridge_r2(Xf[keep][~ee],yy[keep][~ee],seed)
    return dict(cat=j,n=len(keep),diff=a-b,early=a,late=b,
                dy=float(yy[keep][ee].mean()-yy[keep][~ee].mean()))
rows=[]
for j,ri in best.items():
    for tag,pl,sh in [('real',0.,False),('shuf',0.,True),
                      ('plant10',0.10,False),('plant25',0.25,False)]:
        r=run_cat(j,ri,plant=pl,shuffle=sh)
        if r: rows.append(dict(arm=tag,**r))
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
S=D.groupby('arm').agg(k=('diff','size'),n=('n','median'),diff=('diff','mean'),
    se=('diff',lambda s:s.std()/np.sqrt(len(s))),early=('early','mean'),late=('late','mean'),
    dy=('dy','mean'))
print("\n=== 早-A 组 减 晚-A 组:用其余偏好预测 A 的评分的留出 R2 ===")
print(S.round(4).to_string())
g=Gate("兴趣来得早,是否与其余偏好结合得更紧?")
g.asserted("A 的评分在两组间已匹配 (#96a,断言在代码里)",
           abs(S.loc['real','dy'])<0.05, f"|Δ评分| = {abs(S.loc['real','dy']):.4f} < 0.05")
g.negative_control("早/晚标签在匹配格内打乱", null=S.loc['shuf','diff'],
                   effect=S.loc['real','diff'])
g.positive_control("种植的重塑被测出", planted=S.loc['plant25','diff'],
                   floor=S.loc['shuf','diff'], spread=S.loc['plant25','se'])
g.no_sign_crossing("种植阶梯", [S.loc['shuf','diff']+1e-6,S.loc['plant10','diff'],S.loc['plant25','diff']])
g.resolvable("早减晚的差值", effect=S.loc['real','diff'], spread=S.loc['real','se'])
print(); print(g)
if g.verdict():
    print(f"\n  -> 早到来的兴趣与其余偏好结合更紧:R2 差 {S.loc['real','diff']:+.4f}")
else:
    print(f"\n  -> 未确立。早 {S.loc['real','early']:+.4f}  晚 {S.loc['real','late']:+.4f}  "
          f"差 {S.loc['real','diff']:+.4f} ± {S.loc['real','se']:.4f}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
