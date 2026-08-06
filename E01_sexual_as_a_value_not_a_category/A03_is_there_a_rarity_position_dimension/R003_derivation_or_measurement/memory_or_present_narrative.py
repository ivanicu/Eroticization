import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A13 R01 -- 那个回忆偏差本身是什么?一个真的判别式。

#114a 量出:心爱的性兴趣被报告得更早,-0.2000 年/评分标准差(19.8x SE),约 0.74 年。
我一直只把它当干扰用(#115 #116 #117)。它本身是这个项目里最有意思的未开采对象,
而且两个心理学机制做**相反**的预测:

  记忆重构:  畸变来自回忆过程 -> 回忆间隔越长,重构越多 -> **受访者年龄越大,斜率越负**
  当下叙事:  畸变来自报告时刻的自我叙事整合(把现在最爱的说成"从小就有")
              -> 与流逝时间无关 -> **斜率与受访者年龄无关**

这不是干扰控制,这是问"人是怎么构造自己的性欲自传的"。

⚠ 跑之前写下的最强混淆:
  (a) 年龄大的人记录的起始项更多 -> 样本不同。控制:每层内报出 n 与项数
  (b) 回忆间隔 = 当前年龄 − 起始年龄,与起始年龄纠缠 -> 用**受访者年龄**分层,
      它对"哪个类别"是外生的;类别固定效应已在偏离定义里
  (c) 年龄大的人评分分布可能不同 -> 每层内评分标准化

ESTIMAND        #114 的斜率(时间表偏离 对 评分),按受访者年龄分层估计;以及形式化的交互项。
IDENTIFICATION  identified;人群时间表用留一人计算,受访者年龄不进入偏离的定义。
WORLDS          memory     斜率随年龄单调更负
                narrative  斜率与年龄无关
                mixed      有趋势但不单调
KILL            threshold-free;交互斜率对自身自助 SE(按人自助)。
POSITIVE CTRL   (1) 种一个真的随年龄增强的畸变 -> 必须测出趋势
                (2) 种一个平的畸变 -> **不得**产生假趋势(#110d 的教训)
NEGATIVE CTRL   评分在同类别内打乱 -> 每层斜率必须为零。
IMPOSSIBLE      年龄是分箱的(14-17 ... 29-32),跨度只有约 15 年;若畸变在更长尺度上才增长,
                本设计看不见。这是范围限制,不是零。
"""
import pandas as pd, numpy as np, warnings, hashlib
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
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
AGEMID={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=df['age'].map(AGEMID).values
import re
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
print(f"能对上评分列的类别 {len(best)}   年龄有效的人 {np.isfinite(age).sum():,}",flush=True)
def build(Vm,agev):
    prec=np.nanmean(Vm,axis=1)
    out=[]; avail=0; used=0
    for j,ri in best.items():
        col=Vm[:,j]; y=R[:,ri]
        m=np.isfinite(col)&np.isfinite(y)&(np.isfinite(Vm).sum(1)>=6)&np.isfinite(prec)&np.isfinite(agev)
        avail+=1
        idx=np.flatnonzero(m)
        if len(idx)<400: continue
        used+=1
        s=col[idx].sum(); n=len(idx)
        loo=(s-col[idx])/(n-1)
        dev=(col[idx]-loo)-(prec[idx]-np.nanmean(prec[idx]))
        dev=dev-dev.mean()
        yz=(y[idx]-np.nanmean(y[idx]))/(np.nanstd(y[idx])+1e-9)   # 每类别内标准化评分
        out.append(pd.DataFrame(dict(v_cat=j,person=idx,dev=dev,rz=yz,
                                     v_age=agev[idx],breadth=breadth[idx],prec=prec[idx])))
    # #118c 要求:截断必须显式声明并打印。2 个类别因 n<400 被排除,这是纳入标准不是成本控制。
    check_coverage(used,avail,'A13R01 build',tol=0.10)
    if used<avail: print(f"  纳入 {used}/{avail} 类别(其余 n<400)",flush=True)
    return check_columns(pd.concat(out,ignore_index=True),'A13R01')
def slope_by(L,col=None,B=800,seed=0):
    """按人自助的斜率;col 给定时按该列分层。"""
    def fit(Lx):
        Z=np.c_[np.ones(len(Lx)),Lx.rz.values,
                (Lx.breadth.values-Lx.breadth.mean())/(Lx.breadth.std()+1e-9),
                (Lx.prec.values-Lx.prec.mean())/(Lx.prec.std()+1e-9)]
        return np.linalg.lstsq(Z,Lx.dev.values,rcond=None)[0][1]
    ppl=L.person.unique(); gi={p:np.flatnonzero(L.person.values==p) for p in ppl}
    rb=np.random.default_rng(seed)
    def boot(Lx):
        b=fit(Lx); pp=Lx.person.unique()
        g2={p:np.flatnonzero(Lx.person.values==p) for p in pp}
        bs=[]
        for _ in range(B):
            ix=np.concatenate([g2[p] for p in rb.choice(pp,len(pp))])
            bs.append(fit(Lx.iloc[ix]))
        return b,float(np.std(bs))
    if col is None: return boot(L)
    return {v:boot(L[L[col]==v]) for v in sorted(L[col].unique())}
L=build(V,age)
print(f"长表 {len(L):,} 行  {L.person.nunique():,} 人",flush=True)
b_all,se_all=slope_by(L)
by=slope_by(L,'v_age',B=500,seed=1)
rows=[dict(arm='real',v_age=a,slope=v[0],se=v[1],
           n=int((L.v_age==a).sum()),ppl=int(L[L.v_age==a].person.nunique())) for a,v in by.items()]
# 负对照:评分在同类别内打乱
Ls=L.copy(); rp=np.random.default_rng(9)
Ls['rz']=Ls.groupby('v_cat').rz.transform(lambda s: rp.permutation(s.values))
bys=slope_by(Ls,'v_age',B=300,seed=2)
rows+=[dict(arm='shuf',v_age=a,slope=v[0],se=v[1],n=0,ppl=0) for a,v in bys.items()]
# 正对照:(1) 随年龄增强的畸变 (2) 平的畸变
for tag,fn in [('plant_age',lambda ag: 0.06*(ag-np.nanmean(ag))),
               ('plant_flat',lambda ag: np.full_like(ag,0.5))]:
    Vp=V.copy()
    k=fn(age)
    for j,ri in best.items():
        z=(R[:,ri]-np.nanmean(R[:,ri]))/(np.nanstd(R[:,ri])+1e-9)
        Vp[:,j]=Vp[:,j]-np.nan_to_num(z)*np.nan_to_num(k)
    Lp=build(Vp,age); byp=slope_by(Lp,'v_age',B=250,seed=3)
    rows+=[dict(arm=tag,v_age=a,slope=v[0],se=v[1],n=0,ppl=0) for a,v in byp.items()]
D=check_columns(pd.DataFrame(rows),'A13R01 out')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print(f"\n=== #114 的斜率,按受访者年龄分层(负 = 心爱的被报告得更早) ===")
print(f"  全样本 {b_all:+.4f} ± {se_all:.4f}\n")
piv=D.pivot_table(index='v_age',columns='arm',values='slope')
n=D[D.arm=='real'].set_index('v_age')[['n','ppl']]
print(piv[['real','shuf','plant_age','plant_flat']].join(n).round(4).to_string())
def trend(arm):
    d=D[D.arm==arm]; x=d.v_age.values; y=d.slope.values; w=1/np.maximum(d.se.values,1e-9)**2
    Z=np.c_[np.ones(len(x)),(x-x.mean())/x.std()]
    W=np.diag(w)
    b=np.linalg.solve(Z.T@W@Z,Z.T@W@y)[1]
    cov=np.linalg.inv(Z.T@W@Z)
    return float(b),float(np.sqrt(cov[1,1]))
tr={a:trend(a) for a in ['real','shuf','plant_age','plant_flat']}
print("\n=== 年龄趋势(斜率随年龄的变化率;负 = 年龄越大畸变越强) ===")
for a,(b,s) in tr.items(): print(f"  {a:11s} {b:+.4f} ± {s:.4f}   {abs(b)/max(s,1e-9):.1f}x")
g=Gate("回忆畸变:记忆重构,还是当下叙事?")
g.negative_control("评分同类别内打乱(全样本斜率)",
                   null=float(np.mean([v[0] for v in bys.values()])), effect=b_all)
g.resolvable("种植的年龄趋势被测出", effect=tr['plant_age'][0], spread=tr['plant_age'][1])
g.asserted("平种植不产生假年龄趋势",
           abs(tr['plant_flat'][0])<2*tr['plant_flat'][1],
           f"|{tr['plant_flat'][0]:+.4f}| < {2*tr['plant_flat'][1]:.4f}")
g.resolvable("真实的年龄趋势", effect=tr['real'][0], spread=tr['real'][1])
print(); print(g)
b,s=tr['real']
if abs(b)<2*s:
    print(f"\n  -> 当下叙事侧:畸变**不随年龄增强**({b:+.4f} ± {s:.4f})。")
    print(f"     它不像是回忆间隔越长丢失越多,更像是报告时刻把现在最爱的整合进自传。")
elif b<0:
    print(f"\n  -> 记忆重构侧:年龄越大畸变越强({b:+.4f} ± {s:.4f})。")
else:
    print(f"\n  -> 反向趋势({b:+.4f}),两个机制都不预测。")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
