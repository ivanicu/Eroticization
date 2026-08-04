import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A06 R03 -- 用强制单选打 #26 那条 UNVERIFIED。

#61:13,530 个有恋物的人里 82.7% 说色情诱导了本来不会有的兴趣,而这句自述没有时序签名、
     没有结构签名,只追踪"这个人整体勾了多少"。
#26:但"其中 85% 是作答风格"被降级为 UNVERIFIED —— **全部题项都是情欲内容且无反向计分**,
     "泛泛同意"与"泛泛认可情欲事物"分不开。
#123b:强制单选**按构造消除作答水平**(必须且只能选一个),实测各选项人群的平均给分极差
      0.324–0.506 < 0.73 sd。**它正是 #26 缺的那个仪器,而 166 轮里只用过 2 轮。**

问题一句话:**自称"色情诱导了我"的人,他们选的东西和别人不一样吗?**

  携带内容: 强制单选的分布不同 -> 那句自述毕竟说了关于偏好的什么
  只是水平: 在广度匹配之后分布相同 -> #61 的"只追踪你勾了多少"第一次拿到免于作答风格的确认

⚠ 跑之前写下的混淆:
  (a) 广度 —— 诱导自述与广度相关 rho +0.2922 (#26)。必须匹配,且断言匹配成立(#96a)
  (b) 性别·年龄 -> 协变量
  (c) 进入某个强制单选题本身是门控的 -> 只在**都答了该题**的人之间比

ESTIMAND        用强制单选(one-hot)预测诱导自述的留出 AUC 增量,零用合成无信号世界 (#109e)。
KILL            threshold-free;零的种类必须命名 (#109c);gate 顺序按 #120d。
POSITIVE CTRL   把诱导标签**部分替换**为一个由选择决定的标签 -> 必须单调开火;
                g=0 必须等于真实臂(退化不开火,#118d 的教训:不要在投影层面种)。
NEGATIVE CTRL   诱导标签在(广度 x 性别)分层内打乱。
IMPOSSIBLE      因果方向;以及"诱导"自述本身的效度。
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns, check_coverage
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
IND=[c for c in df.columns if c.startswith('Do you feel as though you\'ve "induced"')][0]
raw=df[IND]
ind=raw.map(lambda s: 1.0 if isinstance(s,str) and s.startswith('Yes') else (0.0 if s=='No' else np.nan))
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
P=(df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float).fillna(0.).values
breadth=P.sum(1)
lvl=df[rate].apply(pd.to_numeric,errors='coerce').mean(axis=1).values
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
zs=lambda X:(X-X.mean(0))/(X.std(0)+1e-9)
COV=zs(np.c_[male,agev,breadth,np.nan_to_num(lvl)])
FC=inv[inv['kind']=='FORCED_CHOICE_MOST']['col'].tolist()
print(f"诱导自述: Yes {int((ind==1).sum()):,}  No {int((ind==0).sum()):,}  "
      f"(不适用/缺失 {int(ind.isna().sum()):,})")
print(f"强制单选列 {len(FC)}",flush=True)
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<15 or n0<15: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,seed,alpha=30.,reps=8):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<15 or (1-y[tr]).sum()<15: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def eff(H,C,y,seed=1,ndraw=3):
    """强制单选 one-hot 在协变量之上的 AUC 增量,减去合成无信号世界的偏移。"""
    inc=lambda yy: ridge_auc(np.c_[C,H],yy,seed)-ridge_auc(C,yy,seed)
    X=np.c_[np.ones(len(y)),C]; w=np.linalg.lstsq(X,y,rcond=None)[0]; lin=np.clip(X@w,0.02,0.98)
    off=np.nanmean([inc((np.random.default_rng(seed+400+d).random(len(y))<lin).astype(float))
                    for d in range(ndraw)])
    return inc(y)-off
rows=[]; used=0
for c in FC:
    s=df[c]
    m=np.isfinite(ind)&s.notna()
    idx=np.flatnonzero(m)
    if len(idx)<800: continue
    used+=1
    y=ind.values[idx]; C=COV[idx]
    cats=s.iloc[idx].astype('category')
    if len(cats.cat.categories)<3: continue
    H=pd.get_dummies(cats,drop_first=True).values.astype(float)
    # #96a:广度匹配必须断言。五分层残差 +0.58~+0.72(诊断于 R03 第一次跑),
    # 因为诱导自述与广度 rho +0.29 而分层太粗。改用**卡尺 1:1 匹配**:
    # 每个 No 配一个广度最接近的 Yes,超出卡尺则丢弃。
    b_=breadth[idx]; SD=float(np.std(b_))
    CAL=0.25*SD                                   # 卡尺 = 广度自身 sd 的四分之一
    order=np.argsort(b_)                          # 用排序做最近邻,避免 O(n^2)
    yes=np.flatnonzero(y==1); no=np.flatnonzero(y==0)
    yb=b_[yes]; ys=np.argsort(yb); yes_s=yes[ys]; yb_s=yb[ys]
    taken=np.zeros(len(yes_s),bool); keep=[]
    for j in no:
        k=np.searchsorted(yb_s,b_[j])
        best=None
        for d in range(0,60):
            for kk in (k-d,k+d):
                if 0<=kk<len(yes_s) and not taken[kk] and abs(yb_s[kk]-b_[j])<=CAL:
                    best=kk; break
            if best is not None: break
        if best is not None:
            taken[best]=True; keep+= [j,yes_s[best]]
    if len(keep)<400: continue
    keep=np.array(keep)
    db=float(b_[keep][y[keep]==1].mean()-b_[keep][y[keep]==0].mean())
    # #102a:容差不写成常数,写成广度自身 sd 的比例
    if abs(db)>0.05*SD: continue
    rp=np.random.default_rng(88)
    e=eff(H[keep],C[keep],y[keep])
    yshuf=y[keep].copy()
    st=np.digitize(b_[keep],np.quantile(b_[keep],[.33,.66]))*2+(male[idx][keep]>0).astype(int)
    for q in np.unique(st):
        w_=np.flatnonzero(st==q)
        if len(w_)>1: yshuf[w_]=y[keep][w_][rp.permutation(len(w_))]
    en=eff(H[keep],C[keep],yshuf,seed=2)
    # 正对照:把一部分标签替换成由选择决定的标签(在**标签**层面,不在投影层面 #118d)
    ps={}
    lab=(H[keep]@np.random.default_rng(3).normal(size=H.shape[1])>0).astype(float)
    for g in [0.0,0.10,0.25]:
        yp=np.where(np.random.default_rng(4).random(len(keep))<g,lab,y[keep])
        ps[g]=eff(H[keep],C[keep],yp,seed=5)
    rows.append(dict(v_col=c[:24],n=len(keep),k=H.shape[1]+1,db=db,e=e,e_null=en,
                     p0=ps[0.0],p10=ps[0.10],p25=ps[0.25]))
    print(f"  {c[:22]:<22} n={len(keep):5,} 选项{H.shape[1]+1:3d}  e={e:+.4f} 零={en:+.4f} "
          f"种植 {ps[0.0]:+.4f}/{ps[0.10]:+.4f}/{ps[0.25]:+.4f}",flush=True)
check_coverage(len(rows),len(FC),'A06R03',tol=0.35)
print(f"  纳入 {len(rows)}/{len(FC)} 个强制单选题(其余 n 不足或匹配不成立)")
D=check_columns(pd.DataFrame(rows),'A06R03')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
print(f"\n=== 用强制单选预测「色情诱导了我」的自述,广度已匹配 ===")
print(D.round(4).to_string(index=False))
E=D.e.values; N=D.e_null.values
print(f"\n  效应 {E.mean():+.4f} ± {se(E):.4f}   零 {N.mean():+.4f} ± {se(N):.4f}   {len(D)} 题")
print(f"  种植阶梯 {D.p0.mean():+.4f} / {D.p10.mean():+.4f} / {D.p25.mean():+.4f}")
g=Gate("自称被色情诱导的人,选的东西不一样吗?")
g.asserted("广度匹配成立(#96a,卡尺 1:1,容差 = 0.05 sd)", abs(D.db.mean())<0.05*float(np.std(breadth)),
           f"|Δ广度| = {abs(D.db.mean()):.3f} < {0.05*float(np.std(breadth)):.3f}")
g.asserted("g=0 等于真实臂(退化种植不开火)", abs(D.p0.mean()-E.mean())<1e-9,
           f"{D.p0.mean():+.4f} vs {E.mean():+.4f}")
g.no_sign_crossing("种植阶梯(相对 g=0)", [D.p10.mean()-D.p0.mean(), D.p25.mean()-D.p0.mean()])
g.positive_control("种植 0.25 被测出", planted=D.p25.mean()-D.p0.mean(), floor=0.0, spread=se(E))
g.require_resolvable_first("真实效应", effect=E.mean(), spread=se(E))
g.negative_control("诱导标签在(广度x性别)分层内打乱", null=N.mean(), effect=E.mean())
print(); print(g)
if g.verdict():
    print(f"\n  -> 携带内容:那句自述不只是「我勾得多」")
else:
    print(f"\n  -> 强制单选看不出差别。#61 的「只追踪你勾了多少」拿到免于作答风格的确认" 
          if abs(E.mean())<2*se(E) else "\n  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
