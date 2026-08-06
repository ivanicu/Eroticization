import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R02 -- 稳定类型,还是表征重塑?#106e 的分离器。

#106 发现顺序携带信息(+0.0254, 11.1x SE, 104/120 对为正),并在跑之前就写下了最强替代解释:
一个稳定的"类型"同时驱动顺序和内容。偏关系型的人既更晚获得关系性兴趣,也有关系性偏好 ——
顺序通过类型预测轮廓,而不是通过重塑。

两个分离器,同一次运行:

  (1) 对的种类。 具体 = 物体/身体属性/物质,不需要建模他人意图。
                关系 = 需要建模他人的意图、地位、同意或关系本身。
      类型驱动 -> 效应集中在跨类对;重塑驱动 -> 同类对里也在。
  (2) 直接吸收类型。 把每个人自己的"具体类平均起始年龄"和"关系类平均起始年龄"放进基线协变量。
      这直接吸收"你是哪一类"。效应在这之后还活着 -> 类型被排除得更彻底。

ESTIMAND        顺序的留出 AUC 增量,按对的种类分层,并在加入类型协变量前后各测一次。
IDENTIFICATION  A、B 本身及其派生量不进预测器。种类划分在看结果之前写死在代码里。
CONFOUNDS(跑之前写下):
                (a) 跨类对的起始年龄间隔更大,顺序更"确定",y 的方差更小 -> 报告每类的 y 平衡度
                (b) 同类对可能功率更低 -> 种植正对照必须在每一类里分别跑,否则同类的零无意义
WORLDS          type          效应几乎全在跨类对,且加入类型协变量后崩塌
                restructuring 同类对里也在,且在类型协变量之后存活
KILL            threshold-free;零是过拟合偏移,用 offset_control(#106c)。
POSITIVE CTRL   种植阶梯在跨类和同类各跑一遍 —— 一类里检测不到的话,那一类的零不可读。
NEGATIVE CTRL   分层置换,每类各自的偏移。
IMPOSSIBLE      排除"第三共因作用在比具体/关系更细的层次上" —— 本设计只排除这一个层次的类型。
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
# 种类划分,写死在代码里,在看任何结果之前
CONCRETE={1,2,4,13,14,16,17,18,20,26,29,30}
RELATIONAL={3,5,6,8,9,10,11,12,15,21,22,23,24,27}
KIND={**{i:'C' for i in CONCRETE},**{i:'R' for i in RELATIONAL}}
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
P=(df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float).fillna(0.).values
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
breadth=P.sum(1)
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
cidx=sorted(CONCRETE); ridx=sorted(RELATIONAL)
mc=np.nanmean(V[:,cidx],axis=1); mr=np.nanmean(V[:,ridx],axis=1)
mc=np.where(np.isfinite(mc),mc,np.nanmean(mc)); mr=np.where(np.isfinite(mr),mr,np.nanmean(mr))
BASE=np.c_[male,agev,breadth,prec]
TYPE=np.c_[male,agev,breadth,prec,mc,mr,mc-mr]     # 直接吸收"你是哪一类"
BASE=(BASE-BASE.mean(0))/(BASE.std(0)+1e-9); TYPE=(TYPE-TYPE.mean(0))/(TYPE.std(0)+1e-9)
print(f"具体类 {len(CONCRETE)}  关系类 {len(RELATIONAL)}  偏好项 {P.shape[1]}",flush=True)
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<10 or n0<10: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,rng,alpha=50.,reps=5):
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
rows=[]
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)]
np.random.default_rng(3).shuffle(pairs)
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)
    kind='same' if KIND[a]==KIND[b] else 'cross'
    st=strata(idx); rec={}
    for tag,COVm in [('base',BASE),('type',TYPE)]:
        base=ridge_auc(COVm[idx],y,np.random.default_rng(1))
        full=ridge_auc(np.c_[COVm[idx],P[idx]],y,np.random.default_rng(1))
        pm=[]
        for d in range(2):
            rp=np.random.default_rng(200+d); yp=y.copy()
            for s in np.unique(st):
                w=np.flatnonzero(st==s)
                if len(w)>1: yp[w]=y[w][rp.permutation(len(w))]
            pm.append(ridge_auc(np.c_[COVm[idx],P[idx]],yp,np.random.default_rng(1))-
                      ridge_auc(COVm[idx],yp,np.random.default_rng(1)))
        rec[tag]=(full-base)-np.nanmean(pm)
    rows.append(dict(a=a,b=b,kind=kind,n=len(idx),bal=min(y.mean(),1-y.mean()),
                     gap=float(np.nanmean(np.abs(V[idx,a]-V[idx,b]))),
                     eff_base=rec['base'],eff_type=rec['type']))
    if len(rows)%25==0: print(f"  {len(rows)} pairs",flush=True)
    if len(rows)>=140: break
# 正对照:两类各跑一遍
ctrl=[]
for kind in ['cross','same']:
    for g in [0.0,0.25]:
        cs=[]
        for a,b in pairs:
            k='same' if KIND[a]==KIND[b] else 'cross'
            if k!=kind: continue
            m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
            if len(idx)<400: continue
            y=(V[idx,a]<V[idx,b]).astype(float)
            rp=np.random.default_rng(7); w=rp.normal(size=P.shape[1]); sig=P[idx]@w
            yp=(rp.random(len(idx))<(y*(1-g)+(sig>np.median(sig)).astype(float)*g)).astype(float)
            cs.append(ridge_auc(np.c_[BASE[idx],P[idx]],yp,np.random.default_rng(1))-
                      ridge_auc(BASE[idx],yp,np.random.default_rng(1)))
            if len(cs)>=10: break
        ctrl.append(dict(kind=kind,g=g,inc=np.nanmean(cs)))
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print(f"\n=== {len(D)} 对   跨类 {(D.kind=='cross').sum()}  同类 {(D.kind=='same').sum()} ===")
S=D.groupby('kind').agg(n=('n','median'),bal=('bal','mean'),gap=('gap','mean'),
                        eff_base=('eff_base','mean'),eff_type=('eff_type','mean'),
                        pos=('eff_base',lambda s:(s>0).mean()),k=('n','size'))
print(S.round(4).to_string())
print("\n=== 正对照:两类各自的种植阶梯 ===")
print(C.pivot_table(index='kind',columns='g',values='inc').round(4).to_string())
se=lambda s: s.std()/np.sqrt(len(s))
cb=D[D.kind=='cross'].eff_base; sb=D[D.kind=='same'].eff_base
ct=D[D.kind=='cross'].eff_type; st_=D[D.kind=='same'].eff_type
g=Gate("是稳定类型,还是表征重塑?")
pw=C.pivot_table(index='kind',columns='g',values='inc')
g.positive_control("种植在同类对里能被检测", planted=pw.loc['same',0.25],
                   floor=pw.loc['same',0.0], spread=se(sb))
g.positive_control("种植在跨类对里能被检测", planted=pw.loc['cross',0.25],
                   floor=pw.loc['cross',0.0], spread=se(cb))
g.resolvable("同类对的效应(基线协变量)", effect=sb.mean(), spread=se(sb))
g.resolvable("同类对的效应(加入类型协变量后)", effect=st_.mean(), spread=se(st_))
print(); print(g)
if g.verdict():
    print(f"\n  跨类  基线 {cb.mean():+.4f}  加类型协变量后 {ct.mean():+.4f}  "
          f"保留 {100*ct.mean()/max(cb.mean(),1e-9):.0f}%")
    print(f"  同类  基线 {sb.mean():+.4f}  加类型协变量后 {st_.mean():+.4f}  "
          f"保留 {100*st_.mean()/max(sb.mean(),1e-9):.0f}%")
    if st_.mean()>2*se(st_) and st_.mean()>0.5*ct.mean():
        print("\n  -> 重塑侧。同类对里效应仍在,且在直接吸收了「你是哪一类」之后存活。")
        print("     顺序携带的信息不是「你是关系型还是具体型」。")
    elif st_.mean()<2*se(st_):
        print("\n  -> 类型侧。同类对里效应消失,#106 的顺序痕迹是稳定类型在两处的同一投影。")
    else:
        print("\n  -> 两者都有,按大小分不开。")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
