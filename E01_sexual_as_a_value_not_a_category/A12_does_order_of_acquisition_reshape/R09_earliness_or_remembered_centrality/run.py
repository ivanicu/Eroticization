import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R09 -- "来得早" 还是 "记得它是核心"?#112e 唯一没排除的对手。

#112 立住:一个兴趣来得越早,其余偏好越能预测它(R2 0.0437 -> 0.0675,+54%),而且两组对它的
评分逐格匹配(Δ=0.0000),所以不是"喜欢得更强"。#112e 写下唯一没排除的对手:

  回忆偏差:一个人把 X 当作自己的核心 -> 既报告更早的起始年龄,也让其余轮廓围绕 X 组织。
            "早" 不是真实时间,是"核心性"穿着时间的衣服。

两个分离器,同一次运行,而且是不同种类的证据:

  (1) 交互签名。 回忆偏差应在评分最高的兴趣上最强 -> 早晚差应随评分升高而增大。
                 重塑不预测这个交互(匹配已在每格内固定评分)。
                 做法:只用低评分格(0-2)跑一遍,只用高评分格(3-5)跑一遍。
  (2) 直接替代。 若"早"只是核心性,那么用**直接测到的核心性**匹配之后,早的效应该消失。
                 核心性 = A 的评分在这个人自己所有评分里的百分位(不是绝对评分 —— #112 匹配的是绝对值)。
                 这是比 #112 更锋利的匹配:它固定的是"A 对这个人有多突出",而不是"A 有多高分"。

ESTIMAND        早-A 组减晚-A 组的留出 R2 差,在 (a) 绝对评分匹配 (b) 相对核心性匹配
                (c) 只低评分 (d) 只高评分 四种条件下。
KILL            threshold-free。判定规则写在前面:
                核心性匹配后效应消失 -> "早"就是核心性,#112a 必须改写
                交互存在(高评分格远大于低评分格) -> 回忆偏差
                两者都不 -> #112a 在这个对手面前站住
POSITIVE CTRL   种植的重塑必须在每种匹配下都被测出。
NEGATIVE CTRL   早/晚标签在匹配格内打乱。
IMPOSSIBLE      因果方向。以及一个既改变真实时间又改变核心性的第三共因。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
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
# 核心性:A 的评分在这个人自己所有评分里的百分位
Rz=np.where(np.isfinite(R),R,np.nan)
rk=np.argsort(np.argsort(np.nan_to_num(Rz,nan=-1),axis=1),axis=1).astype(float)
CENT=rk/max(Rz.shape[1]-1,1)
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
def r2(y,p):
    v=np.var(y); return 1-np.mean((y-p)**2)/v if v>0 else np.nan
def ridge_r2(X,y,seed,alpha=50.,reps=6):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(r2(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out)
def run_cat(j,ri,match='abs',band=None,plant=0.0,shuffle=False,seed=1):
    own=V[:,j]; others=np.delete(V,j,axis=1)
    rel=own-np.nanmean(others,axis=1)
    y=R[:,ri].copy(); cent=CENT[:,ri]
    m=np.isfinite(own)&np.isfinite(rel)&np.isfinite(y)&(np.isfinite(others).sum(1)>=5)
    if band=='low':  m&= (y<=2)
    if band=='high': m&= (y>=3)
    idx=np.flatnonzero(m)
    if len(idx)<600: return None
    early=(rel[idx]<np.median(rel[idx]))
    if shuffle: early=np.random.default_rng(seed+400).permutation(early)
    Xp=np.delete(P[idx],ri,axis=1)
    if plant>0:
        w=np.random.default_rng(55).normal(size=Xp.shape[1]); sh=Xp@w
        sh=(sh-sh.mean())/(sh.std()+1e-9)
        y=y.copy(); y[idx]=np.where(early,y[idx]+plant*sh*np.nanstd(y[idx]),y[idx])
    yy=y[idx]; cc=cent[idx]; Xf=np.c_[Xp,COV[idx]]
    grid=np.digitize(yy,[0.5,1.5,2.5,3.5,4.5]) if match=='abs' \
         else np.digitize(cc,np.quantile(cc,[.2,.4,.6,.8]))
    keep=[]
    for gq in np.unique(grid):
        w=np.flatnonzero(grid==gq); e=w[early[w]]; l=w[~early[w]]
        k=min(len(e),len(l))
        if k<25: continue
        rp=np.random.default_rng(seed+7)
        keep+=list(rp.permutation(e)[:k])+list(rp.permutation(l)[:k])
    if len(keep)<300: return None
    keep=np.array(keep); ee=early[keep]
    dy=float(yy[keep][ee].mean()-yy[keep][~ee].mean())
    dc=float(cc[keep][ee].mean()-cc[keep][~ee].mean())
    if match=='abs'  and abs(dy)>0.05: return None
    if match=='cent' and abs(dc)>0.03: return None
    a=ridge_r2(Xf[keep][ee],yy[keep][ee],seed); b=ridge_r2(Xf[keep][~ee],yy[keep][~ee],seed)
    return dict(cat=j,n=len(keep),diff=a-b,early=a,late=b,dy=dy,dc=dc)
ARMS=[('abs_all','abs',None,0.,False),('cent_all','cent',None,0.,False),
      ('abs_low','abs','low',0.,False),('abs_high','abs','high',0.,False),
      ('shuf_abs','abs',None,0.,True),('shuf_cent','cent',None,0.,True),
      ('plant_abs','abs',None,.25,False),('plant_cent','cent',None,.25,False)]
rows=[]
for j,ri in best.items():
    for tag,mt,bd,pl,sh in ARMS:
        r=run_cat(j,ri,match=mt,band=bd,plant=pl,shuffle=sh)
        if r: rows.append(dict(arm=tag,**r))
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
S=D.groupby('arm').agg(k=('diff','size'),n=('n','median'),diff=('diff','mean'),
    se=('diff',lambda s:s.std()/np.sqrt(len(s))),early=('early','mean'),late=('late','mean'),
    dy=('dy','mean'),dc=('dc','mean'))
print("\n=== 早-A 组 减 晚-A 组,四种匹配/分带 ===")
print(S.reindex([a[0] for a in ARMS]).round(4).to_string())
g=Gate("「来得早」还是「记得它是核心」?")
g.asserted("绝对评分匹配成立", abs(S.loc['abs_all','dy'])<0.05, f"|Δ评分| {abs(S.loc['abs_all','dy']):.4f}")
g.asserted("核心性匹配成立", abs(S.loc['cent_all','dc'])<0.03, f"|Δ核心性| {abs(S.loc['cent_all','dc']):.4f}")
g.negative_control("打乱标签(绝对匹配)", null=S.loc['shuf_abs','diff'], effect=S.loc['abs_all','diff'])
g.negative_control("打乱标签(核心性匹配)", null=S.loc['shuf_cent','diff'], effect=S.loc['cent_all','diff'])
g.positive_control("种植在核心性匹配下仍被测出", planted=S.loc['plant_cent','diff'],
                   floor=S.loc['shuf_cent','diff'], spread=S.loc['plant_cent','se'])
g.resolvable("核心性匹配后的效应", effect=S.loc['cent_all','diff'], spread=S.loc['cent_all','se'])
print(); print(g)
lo,hi=S.loc['abs_low','diff'],S.loc['abs_high','diff']
slo,shi=S.loc['abs_low','se'],S.loc['abs_high','se']
print(f"\n  交互签名:低评分带 {lo:+.4f} ± {slo:.4f}   高评分带 {hi:+.4f} ± {shi:.4f}")
print(f"  回忆偏差预测「高远大于低」;实测 高−低 = {hi-lo:+.4f} ± {np.hypot(slo,shi):.4f}"
      f"  -> {'有交互,回忆偏差侧' if hi-lo>2*np.hypot(slo,shi) else '无交互,不支持回忆偏差'}")
keep_frac=S.loc['cent_all','diff']/max(S.loc['abs_all','diff'],1e-9)
print(f"\n  核心性替代检验:绝对匹配 {S.loc['abs_all','diff']:+.4f} -> 核心性匹配 "
      f"{S.loc['cent_all','diff']:+.4f}   保留 {100*keep_frac:.0f}%")
print(f"  -> {'「早」不只是核心性,#112a 站住' if keep_frac>0.5 else '「早」大部分是核心性,#112a 必须改写'}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
