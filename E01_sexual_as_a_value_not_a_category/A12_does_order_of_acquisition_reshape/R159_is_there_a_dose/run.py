import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R03 -- 剂量:共因不做的那个预测。

#107 证明同类对里的顺序效应不是"你是哪一类人"(吸收类型后保留 82%,6.8x SE),并在 107e 留下最后
一个对手:一个作用在比"具体/关系"更细层次上的共因。更细的分类学不是解决办法 —— 更锋利的是剂量。

  静态共因:  "我本来就是 A 型" -> 我先得到 A,也长得像 A 型。
              这个解释里,先后之间隔了多久、发生在几岁,都不该有影响。
  表征重塑:  A 在 B 到来之前的那段窗口里重塑表征。
              -> 间隔越大,重塑越多,轮廓差异越大。
              -> 发生得越早(可塑性越高),重塑越多。
  两条都是共因不做的剂量预测,而且不需要任何分类学。

⚠ 跑之前写下的最强混淆:顺序标签在间隔大时测得更可靠(分箱造成的并列更少),所以任何真实效应
  在大间隔处都会看起来更大 —— 这是测量伪影,不是剂量。
  控制:同一个剂量在跨类对里也测。#107 已确定跨类臂是类型驱动的。伪影在两臂应等量出现;
  重塑只应在同类臂更强。而年龄剂量不受这个伪影影响(年龄不是顺序标签)。

ESTIMAND        偏移校正后的 AUC 增量,按 (a) 对内间隔 三分位 (b) 先获得者的年龄 三分位 分层,
                同类臂与跨类臂各自计算。
KILL            threshold-free;零是过拟合偏移 -> offset_control (#106c)。
POSITIVE CTRL   两个,必须同时过:
                (1) 种一个真的随间隔递增的效应 -> 必须测出递增
                (2) 种一个跟间隔无关的平效应   -> 必须测不出剂量(否则剂量是分层伪影)
NEGATIVE CTRL   分层置换,每格各自的偏移。
IMPOSSIBLE      因果方向。剂量让共因更难解释,但不排除一个本身依赖年龄的共因。
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
COV=np.c_[male,agev,breadth,prec,mc,mr,mc-mr]      # 类型协变量已在基线里 (#107)
COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<8 or n0<8: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,rng,alpha=50.,reps=4):
    out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<8 or (1-y[tr]).sum()<8: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def eff(idx,y,seed=1):
    """偏移校正后的增量:满模型 - 协变量模型,减去分层置换的同一量。"""
    base=ridge_auc(COV[idx],y,np.random.default_rng(seed))
    full=ridge_auc(np.c_[COV[idx],P[idx]],y,np.random.default_rng(seed))
    rp=np.random.default_rng(seed+50); yp=y[rp.permutation(len(y))]
    pb=ridge_auc(COV[idx],yp,np.random.default_rng(seed))
    pf=ridge_auc(np.c_[COV[idx],P[idx]],yp,np.random.default_rng(seed))
    return (full-base)-(pf-pb)
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)]
np.random.default_rng(3).shuffle(pairs)
rows=[]; NP=0
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<900: continue
    kind='same' if KIND[a]==KIND[b] else 'cross'
    y=(V[idx,a]<V[idx,b]).astype(float)
    gap=np.abs(V[idx,a]-V[idx,b]); first=np.minimum(V[idx,a],V[idx,b])
    for dose,vals in [('gap',gap),('age_first',first)]:
        t=np.digitize(vals,np.quantile(vals,[1/3,2/3]))
        for lv in [0,1,2]:
            w=np.flatnonzero(t==lv)
            if len(w)<250: continue
            rows.append(dict(a=a,b=b,kind=kind,dose=dose,level=lv,n=len(w),
                             mval=float(np.mean(vals[w])),e=eff(idx[w],y[w])))
    NP+=1
    if NP%15==0: print(f"  {NP} pairs",flush=True)
    if NP>=60: break
# ---- 两个正对照 ----
ctrl=[]
for mode in ['dose_dependent','flat']:
    cs=[]
    for a,b in pairs[:14]:
        m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
        if len(idx)<900: continue
        gap=np.abs(V[idx,a]-V[idx,b]); t=np.digitize(gap,np.quantile(gap,[1/3,2/3]))
        rp=np.random.default_rng(21); wv=rp.normal(size=P.shape[1]); sig=(P[idx]@wv)
        sig=(sig>np.median(sig)).astype(float)
        y0=(V[idx,a]<V[idx,b]).astype(float); y0=y0[rp.permutation(len(y0))]   # 地板已毁
        for lv in [0,1,2]:
            w=np.flatnonzero(t==lv)
            if len(w)<250: continue
            g=(0.10+0.15*lv) if mode=='dose_dependent' else 0.25
            yp=np.where(rp.random(len(w))<g,sig[w],y0[w])
            cs.append(dict(mode=mode,level=lv,e=eff(idx[w],yp,seed=3)))
        if len([c for c in cs if c['mode']==mode])>=27: break
    ctrl+=cs
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print(f"\n=== {D.groupby(['a','b']).ngroups} 对 ===")
for dose in ['gap','age_first']:
    print(f"\n--- 剂量 = {dose} ---")
    t=D[D.dose==dose].groupby(['kind','level']).agg(mval=('mval','mean'),e=('e','mean'),
        se=('e',lambda s:s.std()/np.sqrt(len(s))),k=('e','size'))
    print(t.round(4).to_string())
print("\n=== 正对照 ===")
print(C.groupby(['mode','level']).e.mean().unstack('level').round(4).to_string())
cd=C[C['mode']=='dose_dependent'].groupby('level').e.mean()
cf=C[C['mode']=='flat'].groupby('level').e.mean()
G=D[(D.dose=='gap')&(D.kind=='same')].groupby('level').e
A_=D[(D.dose=='age_first')&(D.kind=='same')].groupby('level').e
Gx=D[(D.dose=='gap')&(D.kind=='cross')].groupby('level').e
g=Gate("顺序效应有剂量吗?共因不预测剂量。")
g.asserted("种植的剂量效应被测出(递增)", cd.iloc[0]<cd.iloc[1]<cd.iloc[2], f"{cd.round(4).tolist()}")
g.asserted("平种植不产生假剂量", abs(cf.iloc[2]-cf.iloc[0])<abs(cd.iloc[2]-cd.iloc[0])/2,
           f"flat {cf.round(4).tolist()} vs dose {cd.round(4).tolist()}")
g.offset_control("同类·间隔剂量 (高-低)", effect=G.mean().iloc[2], offset=G.mean().iloc[0],
                 spread=float(np.hypot(G.sem().iloc[2],G.sem().iloc[0])))
g.offset_control("同类·年龄剂量 (低龄-高龄)", effect=A_.mean().iloc[0], offset=A_.mean().iloc[2],
                 spread=float(np.hypot(A_.sem().iloc[0],A_.sem().iloc[2])))
print(); print(g)
print(f"\n  伪影检查:同类间隔剂量 {G.mean().iloc[2]-G.mean().iloc[0]:+.4f}   "
      f"跨类间隔剂量 {Gx.mean().iloc[2]-Gx.mean().iloc[0]:+.4f}")
print(f"  (测量可靠性伪影应在两臂等量;重塑只应在同类臂更强)")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
