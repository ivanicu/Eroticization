import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A05 R15 -- 广度是内容还是作答风格?用强制单选绕开 #26 的中介问题。

#26 把"广度里 9-13% 是作答风格"降级为 UNVERIFIED,理由不是测错了,而是**无法解释**:
作答风格指标与广度相关 +0.385(而广度就是结局),所以把它当干扰减掉是第四次中介问题;
而"泛泛同意"与"泛泛认可情欲事物"在本 release 分不开,因为全部题项都是情欲内容、无反向计分。

#123b/#124b:强制单选**按构造消除作答水平**,`#26` 的反对在它上面失效。
而且这个设计**根本不需要那个作答风格指标** —— 中介问题因此不存在:

  广度 = 纯作答风格: "什么都说是"在强制单选上无处施力
                    -> 广度不该预测你选哪一个
  广度 = 真实内容:   口味真的广的人,被迫只选一个时会选得**不一样**
                    -> 广度应当预测选择

⚠ 跑之前写下的混淆:
  (a) 性别·年龄 -> 协变量
  (b) 进入某个强制单选题本身是门控的 -> 只在都答了该题的人之间比
  (c) 广度与"答了多少题"共变 -> 把该人答过的强制单选题数放进协变量
  (d) #124f:退化臂必须复用参照臂的种子,并用 degenerate_matches_reference 断言

ESTIMAND        用强制单选 one-hot 预测广度(高 vs 低三分位)的留出 AUC 增量,
                零用合成无信号世界 (#109e)。
KILL            threshold-free;gate 顺序按 #120d;零的种类必须命名 (#109c)。
POSITIVE CTRL   把广度标签部分替换为由选择决定的标签,g=0/0.06/0.15;必须单调,g=0 精确复现真实臂。
NEGATIVE CTRL   广度标签在(性别 x 年龄)分层内打乱。
IMPOSSIBLE      "内容"与"更愿意在问卷上做出区分"分不开;本轮测的是选择是否携带广度信息。
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns, check_coverage
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
P=(R>0).astype(float).fillna(0.).values; breadth=P.sum(1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
FC=inv[inv['kind']=='FORCED_CHOICE_MOST']['col'].tolist()
nfc=df[FC].notna().sum(axis=1).values          # 混淆 (c):这个人答了多少个强制单选题
zs=lambda X:(X-X.mean(0))/(X.std(0)+1e-9)
COV=zs(np.c_[male,agev,nfc])
print(f"广度 均值 {breadth.mean():.1f} sd {breadth.std():.1f}   强制单选列 {len(FC)}",flush=True)
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
    inc=lambda yy: ridge_auc(np.c_[C,H],yy,seed)-ridge_auc(C,yy,seed)
    X=np.c_[np.ones(len(y)),C]; w=np.linalg.lstsq(X,y,rcond=None)[0]; lin=np.clip(X@w,0.02,0.98)
    off=np.nanmean([inc((np.random.default_rng(seed+400+d).random(len(y))<lin).astype(float))
                    for d in range(ndraw)])
    return inc(y)-off
rows=[]
for c in FC:
    s=df[c]; m=s.notna()&np.isfinite(breadth); idx=np.flatnonzero(m)
    if len(idx)<800: continue
    b=breadth[idx]; lo,hi=np.quantile(b,[1/3,2/3])
    sel=np.flatnonzero((b<=lo)|(b>=hi))
    if len(sel)<600: continue
    gi=idx[sel]; y=(breadth[gi]>=hi).astype(float)
    cats=s.iloc[gi].astype('category')
    if len(cats.cat.categories)<3: continue
    H=pd.get_dummies(cats,drop_first=True).values.astype(float)
    C=COV[gi]
    SEED=1                                        # #124f:退化臂与真实臂共用这个种子
    e=eff(H,C,y,seed=SEED)
    lab=(H@np.random.default_rng(3).normal(size=H.shape[1])>0).astype(float)
    ps={}
    for g in [0.0,0.06,0.15]:
        yp=np.where(np.random.default_rng(4).random(len(y))<g,lab,y)
        ps[g]=eff(H,C,yp,seed=SEED)               # 同一个种子
    rp=np.random.default_rng(88); ysh=y.copy()
    st=(male[gi]>0).astype(int)*5+agev[gi].astype(int)
    for q in np.unique(st):
        w_=np.flatnonzero(st==q)
        if len(w_)>1: ysh[w_]=y[w_][rp.permutation(len(w_))]
    en=eff(H,C,ysh,seed=SEED)
    rows.append(dict(v_col=c[:24],n=len(gi),k=H.shape[1]+1,e=e,e_null=en,
                     p0=ps[0.0],p06=ps[0.06],p15=ps[0.15]))
    print(f"  {c[:22]:<22} n={len(gi):5,} 选项{H.shape[1]+1:3d}  e={e:+.4f} 零={en:+.4f} "
          f"种植 {ps[0.0]:+.4f}/{ps[0.06]:+.4f}/{ps[0.15]:+.4f}",flush=True)
check_coverage(len(rows),len(FC),'A05R15',tol=0.35)
print(f"  纳入 {len(rows)}/{len(FC)} 个强制单选题")
D=check_columns(pd.DataFrame(rows),'A05R15')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
E=D.e.values; N=D.e_null.values
print(f"\n=== 用强制单选预测广度(高 vs 低三分位) ===")
print(D.round(4).to_string(index=False))
print(f"\n  效应 {E.mean():+.4f} ± {se(E):.4f}   零 {N.mean():+.4f}   {len(D)} 题")
print(f"  种植阶梯 {D.p0.mean():+.4f} / {D.p06.mean():+.4f} / {D.p15.mean():+.4f}")
g=Gate("广度是内容,还是作答风格?")
g.degenerate_matches_reference("g=0 种植精确复现真实臂 (#124f)",
                               degenerate=D.p0.mean(), reference=E.mean())
g.no_sign_crossing("种植阶梯(相对 g=0)", [D.p06.mean()-D.p0.mean(), D.p15.mean()-D.p0.mean()])
g.positive_control("种植 0.15 被测出", planted=D.p15.mean()-D.p0.mean(), floor=0.0, spread=se(E))
g.require_resolvable_first("真实效应", effect=E.mean(), spread=se(E))
g.negative_control("广度标签在(性别x年龄)分层内打乱", null=N.mean(), effect=E.mean())
print(); print(g)
if g.verdict():
    print(f"\n  -> 广度携带内容:口味广的人被迫只选一个时,选得不一样。")
    print(f"     「广度是作答风格」在一个「什么都说是」无处施力的仪器上被否定。")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
