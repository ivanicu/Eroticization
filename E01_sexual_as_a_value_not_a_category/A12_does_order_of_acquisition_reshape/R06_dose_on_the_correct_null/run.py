import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R06 -- 剂量,这次用对的零。

#109f:#108 和 #109 的剂量结论全部作废,因为两轮都用了平置换 —— 一个把偏移低估 +0.0207 的错零。
#109e 给了更准的偏移来源:合成无信号世界(误差 +0.0040 的分层置换之外,它是直接量出来的真值)。

于是剂量问题回到完全未测,而它仍然是共因唯一不做的预测:

  静态共因: "我本来就是 A 型" -> 先得到 A 也长得像 A 型。隔多久、几岁,不该有影响。
  表征重塑: A 在 B 到来前的窗口里重塑表征 -> 间隔越大重塑越多;越早发生可塑性越高。

功率按 #108e 保:每个对贡献一个用全样本算出的效应,剂量是对的属性,回归在对间做。

ESTIMAND        per-pair 效应(用合成无信号世界的偏移)对该对的 (a) 平均间隔 (b) 平均首获年龄 的斜率,
                同类臂为估计目标,跨类臂作可靠性伪影对照。
IDENTIFICATION  偏移不是置换而是直接测量:把 y 回归到协变量上生成 y_synth(P 无额外信息),
                在 y_synth 上量 (full-base)。这是 #109e 认定的真值来源。
CONFOUNDS(跑之前):
                (a) 间隔大 -> 顺序标签测得更可靠 -> 伪影。若为伪影,同类与跨类斜率应相同
                (b) 对的 n 与间隔相关 -> log n 进回归
                (c) 首获年龄与早熟度相关 -> 早熟度已在个人层协变量里
WORLDS          重塑  同类臂斜率显著且强于跨类臂
                共因  斜率为零
                伪影  两臂斜率相同
KILL            threshold-free;斜率对自身自助 SE。零的种类在 gate 里必须命名 (#109c)。
POSITIVE CTRL   (1) 种一个强度随对间隔递增的效应 -> 必须测出正斜率
                (2) 种一个跟间隔无关的平效应     -> 必须测出零斜率
                (3) 纯置换 y -> 效应必须为零(检验合成偏移本身没有偏)
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
        rp=np.random.default_rng(31); wv=rp.normal(size=P.shape[1])
        sig=(P[idx]@wv>np.median(P[idx]@wv)).astype(float)
        y0=(V[idx,a]<V[idx,b]).astype(float); y0=y0[rp.permutation(len(y0))]
        g=(0.05+0.045*gp) if mode=='dose' else 0.22
        yp=np.where(rp.random(len(idx))<g,sig,y0)
        ctrl.append(dict(mode=mode,gap=gp,e=eff(idx,yp,seed=5))); cs+=1
        if cs>=30: break
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def slope(x,y,z,B=2000,seed=0):
    X=np.c_[np.ones(len(x)),(x-x.mean())/(x.std()+1e-9),(z-z.mean())/(z.std()+1e-9)]
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
        b,sd=slope(sub[dose].values,sub.e.values,np.log(sub.n.values))
        out.append(dict(kind=nm,dose=dose,slope=b,sd=sd))
R=pd.DataFrame(out); print("\n=== 对间剂量斜率(含 log n) ===")
print(R.round(4).to_string(index=False))
cd=C[C['mode']=='dose']; cf=C[C['mode']=='flat']
bd,sdd=slope(cd.gap.values,cd.e.values,np.ones(len(cd))+np.arange(len(cd))*1e-9)
bf,sdf=slope(cf.gap.values,cf.e.values,np.ones(len(cf))+np.arange(len(cf))*1e-9)
print(f"\n=== 正对照 ===\n  种剂量 斜率 {bd:+.4f} ± {sdd:.4f}   平种植 斜率 {bf:+.4f} ± {sdf:.4f}")
se=lambda s: s.std()/np.sqrt(len(s))
sg=R[(R.kind=='same')&(R.dose=='gap')].iloc[0]; sf=R[(R.kind=='same')&(R.dose=='first')].iloc[0]
xg=R[(R.kind=='cross')&(R.dose=='gap')].iloc[0]
g=Gate("#109f 剂量,在对的零上")
g.negative_control("纯置换 y 的效应(合成偏移是否有偏)", null=abs(S.e_null.mean()),
                   effect=abs(S.e.mean()))
g.resolvable("种植的对间剂量被测出", effect=bd, spread=sdd)
g.asserted("平种植不产生假斜率", abs(bf)<2*sdf, f"|{bf:+.4f}| < {2*sdf:.4f}")
g.asserted("同类效应本轮为正(功率前提)", S.e.mean()>2*se(S.e),
           f"{S.e.mean():+.4f} vs 2*SE {2*se(S.e):.4f}")
g.resolvable("同类·间隔斜率", effect=sg.slope, spread=sg.sd)
g.resolvable("同类·首获年龄斜率", effect=sf.slope, spread=sf.sd)
print(); print(g)
print(f"\n  伪影检查:同类间隔斜率 {sg.slope:+.4f}   跨类间隔斜率 {xg.slope:+.4f}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
