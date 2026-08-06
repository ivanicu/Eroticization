import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R10 -- 记忆会不会把心爱的东西往前拉?一个真正判别的分离器。

#113c 承认我上一轮写的"判别式"不判别(两个世界做同一个预测)。这一轮的分离器不碰绑定统计量,
测的是**起始年龄报告本身**:

  回忆偏差: 记忆把心爱的兴趣往前拉 -> 一个人评分高的兴趣,应当系统性地**早于**人群时间表
            对它的预测;评分低的不应当。
  重塑:     起始年龄报告是准的,重塑发生在报告之后的因果链上 -> 偏离与评分无关。

这两条做**相反**的预测,所以它是一个真的判别式。而且它用的是 A03 花了 22 轮建立的人群时间表,
完全不依赖 #112/#113 的绑定统计量。

ESTIMAND        每个(人 x 类别)的**时间表偏离** = 报告起始年龄 − 人群对该类别的中位起始年龄,
                再减掉这个人自己的整体早熟度(否则早熟的人每一项都偏早)。
                然后:偏离对该人对该类别的评分的斜率。
IDENTIFICATION  identified;人群时间表用**留一人**计算,所以一个人不参与自己的基准。
CONFOUNDS(跑之前):
                (a) 早熟度 -> 已在偏离定义里减掉,并断言残差与早熟度不相关
                (b) 广度 -> 广度大的人评分高的项多,进协变量
                (c) 类别本身的评分-年龄关系(有些类别既晚又受欢迎)-> 类别固定效应(按类别中心化)
                (d) 天花板:评分 5 的项无法再高 -> 报告评分的分布
WORLDS          recall     斜率显著为负(评分越高,报告越早)
                accurate   斜率为零
KILL            threshold-free;斜率对自身自助 SE,零的种类必须命名。
POSITIVE CTRL   人为把高评分项的起始年龄往前拉已知的量 -> 必须测出,且单调。
NEGATIVE CTRL   把评分在**同一类别内**打乱 -> 斜率必须为零(保留类别效应,毁掉人-项配对)。
IMPOSSIBLE      "记忆拉前" 与 "早获得的东西后来更被珍视" 分不开 —— 两者都产生负斜率。
                本轮只能测这个负斜率在不在,不能定方向。
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
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
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
    """返回长表:每个 (人, 类别) 一行,含时间表偏离(留一人基准、去早熟)和评分。"""
    prec=np.nanmean(Vm,axis=1)
    rows=[]
    for j,ri in best.items():
        col=Vm[:,j]; y=R[:,ri]
        m=np.isfinite(col)&np.isfinite(y)&(np.isfinite(Vm).sum(1)>=6)&np.isfinite(prec)
        idx=np.flatnonzero(m)
        if len(idx)<400: continue
        s=col[idx].sum(); n=len(idx)
        loo=(s-col[idx])/(n-1)                       # 留一人的人群基准
        dev=(col[idx]-loo)-(prec[idx]-np.nanmean(prec[idx]))   # 去掉这个人的整体早熟
        dev=dev-dev.mean()                            # 类别固定效应
        rows.append(pd.DataFrame(dict(cat=j,person=idx,dev=dev,rating=y[idx],
                                      breadth=breadth[idx],prec=prec[idx])))
    return pd.concat(rows,ignore_index=True)
def slope(L,B=1500,seed=0):
    x=(L.rating.values-L.rating.mean())/L.rating.std()
    Z=np.c_[np.ones(len(L)),x,(L.breadth.values-L.breadth.mean())/L.breadth.std(),
            (L.prec.values-L.prec.mean())/L.prec.std()]
    y=L.dev.values
    b=np.linalg.lstsq(Z,y,rcond=None)[0][1]
    rb=np.random.default_rng(seed)
    # 按人自助,因为同一个人贡献多行
    ppl=L.person.unique(); gi={p:np.flatnonzero(L.person.values==p) for p in ppl}
    bs=[]
    for _ in range(B):
        pick=rb.choice(ppl,len(ppl))
        ix=np.concatenate([gi[p] for p in pick])
        bs.append(np.linalg.lstsq(Z[ix],y[ix],rcond=None)[0][1])
    return float(b),float(np.std(bs))
L=build(V)
print(f"长表 {len(L):,} 行  {L.person.nunique():,} 人  {L.cat.nunique()} 类别",flush=True)
b_real,se_real=slope(L)
# 负对照:评分在同一类别内打乱(保留类别效应,毁掉人-项配对)
Ls=L.copy(); rp=np.random.default_rng(9)
Ls['rating']=Ls.groupby('cat').rating.transform(lambda s: rp.permutation(s.values))
b_shuf,se_shuf=slope(Ls,seed=1)
# 正对照:把高评分项的起始年龄人为往前拉
res=[]
for g in [0.0,0.3,0.8]:
    Vp=V.copy()
    for j,ri in best.items():
        hi=np.isfinite(R[:,ri])&(R[:,ri]>=3)
        Vp[hi,j]=Vp[hi,j]-g
    Lp=build(Vp); bp,sp=slope(Lp,B=600,seed=2)
    res.append(dict(g=g,slope=bp,se=sp))
C=pd.DataFrame(res)
out=pd.DataFrame([dict(arm='real',slope=b_real,se=se_real),
                  dict(arm='rating shuffled within category',slope=b_shuf,se=se_shuf)])
OUT=pathlib.Path(__file__).parent/'results'
pd.concat([out,C.assign(arm=lambda d:'plant g='+d.g.astype(str))]).to_csv(OUT/'grid.csv',index=False)
print("\n=== 时间表偏离 对 评分 的斜率(负 = 心爱的被报告得更早) ===")
print(out.round(4).to_string(index=False))
print("\n=== 正对照:人为把高评分项往前拉 g 年 ===")
print(C.round(4).to_string(index=False))
corr_prec=np.corrcoef(L.dev,L.prec)[0,1]
g=Gate("记忆会不会把心爱的东西往前拉?")
g.asserted("偏离已去掉早熟度(残差与早熟度不相关)", abs(corr_prec)<0.05,
           f"corr(dev, precocity) = {corr_prec:+.4f}")
g.negative_control("评分在同类别内打乱", null=b_shuf, effect=b_real)
g.no_sign_crossing("种植阶梯(拉前 g 年 -> 斜率应更负)",
                   [-(C.slope.iloc[0]-C.slope.iloc[0])-1e-9,
                    -(C.slope.iloc[1]-C.slope.iloc[0]),-(C.slope.iloc[2]-C.slope.iloc[0])])
g.positive_control("种植 0.8 年被测出", planted=-(C.slope.iloc[2]-C.slope.iloc[0]),
                   floor=0.0, spread=C.se.iloc[2])
g.resolvable("真实斜率", effect=b_real, spread=se_real)
print(); print(g)
print(f"\n  真实斜率 {b_real:+.4f} 年/评分标准差   种植 0.8 年给 {C.slope.iloc[2]-C.slope.iloc[0]:+.4f}")
if abs(b_real)>2*se_real:
    eq=0.8*abs(b_real)/max(abs(C.slope.iloc[2]-C.slope.iloc[0]),1e-9)
    print(f"  -> 相当于把高评分项往前拉了约 {eq:.2f} 年")
    print(f"  -> {'回忆偏差侧:心爱的确实被报告得更早' if b_real<0 else '反向:心爱的被报告得更晚'}")
else:
    print(f"  -> 斜率不可分辨,起始年龄报告没有随评分系统性偏移")
print(f"\nartifact sha1 {hashlib.sha1(out.to_csv(index=False).encode()).hexdigest()[:12]}")
