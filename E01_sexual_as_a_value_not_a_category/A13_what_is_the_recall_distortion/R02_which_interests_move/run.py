import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A13 R02 -- 哪些兴趣在移动?记忆衰退 vs 叙事固化。

#119c 的残余:"随时间累积"对两个机制都成立。#119 的 NEXT 提议"低评分兴趣单独看",但设计时发现
一个硬约束:dev 的基准是**样本内**的(留一人),所以在任一分层里 mean(dev) 被钉在 0 附近,
高评分组下移就机械地逼着低评分组上移。**两组之差是被构造强制的,不能读。**

能读的是**形状**:mean(dev) 沿评分 0..5 的曲线,在不同年龄层之间怎么变。加权均值被钉住,
但曲线的**形状**没有被钉住。两个机制预测不同的形状变化:

  叙事固化: 故事是关于"现在什么重要" -> 只有最高分的项被往前拉
            -> 年龄带来的变化**集中在评分顶端**(凸的)
  记忆衰退: 回忆整体变不准 -> 沿评分的梯度整体变陡
            -> 年龄带来的变化**沿评分线性铺开**

⚠ 跑之前:两个机制都能产生 #119a 的斜率增长,所以斜率不判别(#113c 的教训)。形状才判别。

ESTIMAND        mean(dev) 在 (评分等级 x 年龄层) 上的表;年龄带来的变化沿评分的**曲率**
                (二次项系数) 对**线性项**的比。
IDENTIFICATION  identified 相对于两个正对照;曲率本身没有绝对含义,只有与已知形状的对照才有。
KILL            threshold-free:真实曲率落在两个已知形状之间的哪一侧,以及自助 SE。
POSITIVE CTRL   (1) 顶端种植:只把评分>=4 的项按年龄拉前 -> 必须测出**凸**形状
                (2) 线性种植:按 z(评分) 线性拉前,同样按年龄 -> 必须测出**线性**形状
                两者必须**彼此可分**,否则形状读数不被许可。
NEGATIVE CTRL   评分在同类别同年龄层内打乱 -> 曲线必须平。
IMPOSSIBLE      "现在重要"与"当时就重要"分不开;本轮只测形状,不定因果。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns, check_coverage
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float)
AGEMID={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=df['age'].map(AGEMID).values
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
def build(Vm):
    prec=np.nanmean(Vm,axis=1); out=[]; avail=0; used=0
    for j,ri in best.items():
        col=Vm[:,j]; y=R[:,ri]
        m=np.isfinite(col)&np.isfinite(y)&(np.isfinite(Vm).sum(1)>=6)&np.isfinite(prec)&np.isfinite(age)
        avail+=1; idx=np.flatnonzero(m)
        if len(idx)<400: continue
        used+=1
        s=col[idx].sum(); n=len(idx)
        loo=(s-col[idx])/(n-1)
        dev=(col[idx]-loo)-(prec[idx]-np.nanmean(prec[idx]))
        # #119b:层内标准化,消掉 sd(dev) 随年龄扩大的尺度伪影
        d=pd.DataFrame(dict(v_cat=j,person=idx,dev=dev,rat=y[idx],v_age=age[idx]))
        d['dev']=d.groupby('v_age').dev.transform(lambda s:(s-s.mean())/(s.std()+1e-9))
        out.append(d)
    check_coverage(used,avail,'A13R02 build',tol=0.10)
    if used<avail: print(f"  纳入 {used}/{avail} 类别(其余 n<400)",flush=True)
    return check_columns(pd.concat(out,ignore_index=True),'A13R02')
def shape(L):
    """返回 (曲率, 线性) —— 年龄带来的 mean(dev) 变化沿评分的二次/一次分解。"""
    T=L.groupby(['v_age','rat']).dev.mean().unstack('rat')
    rr=np.array([c for c in T.columns if np.isfinite(T[c]).all()],dtype=float)
    if len(rr)<4: return np.nan,np.nan,T
    M=T[rr].values                                # 年龄 x 评分
    ages=np.array(T.index,dtype=float)
    slope_per_rating=np.array([np.polyfit(ages,M[:,k],1)[0] for k in range(M.shape[1])])
    z=(rr-rr.mean())/(rr.std()+1e-9)
    c=np.polyfit(z,slope_per_rating,2)            # [二次, 一次, 常数]
    return float(c[0]),float(c[1]),T
def boot_shape(L,B=250,seed=0):
    rb=np.random.default_rng(seed); ppl=L.person.unique()
    gi={p:np.flatnonzero(L.person.values==p) for p in ppl}
    q,l,_=shape(L); qs=[];ls=[]
    for _ in range(B):
        ix=np.concatenate([gi[p] for p in rb.choice(ppl,len(ppl))])
        a,b,_=shape(L.iloc[ix]); qs.append(a); ls.append(b)
    return q,float(np.nanstd(qs)),l,float(np.nanstd(ls))
def plant(kind,g=0.10):
    Vp=V.copy()
    for j,ri in best.items():
        r=R[:,ri]; z=(r-np.nanmean(r))/(np.nanstd(r)+1e-9)
        k=g*(age-np.nanmean(age))
        pull = (r>=4).astype(float) if kind=='top' else np.nan_to_num(z)
        Vp[:,j]=Vp[:,j]-np.nan_to_num(pull)*np.nan_to_num(k)
    return Vp
res={}
L=build(V); print(f"长表 {len(L):,} 行  {L.person.nunique():,} 人",flush=True)
res['real']=boot_shape(L,seed=1)
Ls=L.copy(); rp=np.random.default_rng(9)
Ls['rat']=Ls.groupby(['v_cat','v_age']).rat.transform(lambda s: rp.permutation(s.values))
res['shuf']=boot_shape(Ls,B=150,seed=2)
for kind in ['top','linear']:
    res[f'plant_{kind}']=boot_shape(build(plant(kind)),B=150,seed=3)
_,_,T=shape(L)
print("\n=== mean(dev) 表:年龄层 x 评分等级(层内标准化,负 = 早于时间表) ===")
print(T.round(3).to_string())
print("\n=== 年龄带来的变化沿评分的形状分解 ===")
print(f"  {'臂':12s} {'曲率(二次)':>14s} {'线性':>14s}   曲率/线性")
for a,(q,sq,l,sl) in res.items():
    print(f"  {a:12s} {q:+8.5f}±{sq:.5f} {l:+8.5f}±{sl:.5f}   {q/l if abs(l)>1e-9 else float('nan'):+.2f}")
q,sq,l,sl=res['real']; qt,_,lt,_=res['plant_top']; ql,_,ll,_=res['plant_linear']
g=Gate("哪些兴趣在移动?叙事固化(顶端) vs 记忆衰退(线性)")
g.negative_control("评分在同类别同年龄层内打乱(曲率)", null=res['shuf'][0], effect=q)
g.asserted("两个种植形状彼此可分",
           abs(qt/max(abs(lt),1e-9) - ql/max(abs(ll),1e-9))>0.5,
           f"顶端种植 曲率/线性 {qt/lt:+.2f}  vs  线性种植 {ql/ll:+.2f}")
g.resolvable("真实曲率", effect=q, spread=sq)
g.resolvable("真实线性项", effect=l, spread=sl)
print(); print(g)
if g.verdict():
    rr=q/l if abs(l)>1e-9 else float('nan')
    dt=abs(rr-qt/lt); dl=abs(rr-ql/ll)
    print(f"\n  真实 曲率/线性 {rr:+.2f}   顶端种植 {qt/lt:+.2f}   线性种植 {ql/ll:+.2f}")
    print(f"  -> 更接近{'顶端(叙事固化)' if dt<dl else '线性(记忆衰退)'}")
print(f"\nartifact sha1 {hashlib.sha1(str(res).encode()).hexdigest()[:12]}")
