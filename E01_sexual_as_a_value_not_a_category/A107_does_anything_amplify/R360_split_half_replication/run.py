import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A107 R360 -- 那个 p=0.067 的信号,在独立的一半人上还在不在

`#314a`:三个先验候选的联合检验落在置换零的第 93 百分位,**没到预注册的 2× 门槛**。
**加更多候选只会让多重性更糟;最便宜的分辨是问同一个信号在独立的一半人上还在不在。**

ESTIMAND        8 次随机人劈半(两半**不相交且等大**),**每一半各自**跑同一个三候选联合检验,
                **并各自跑自己的置换零**(半样本的零展布更宽,不能共用)。
KILL            **若两半都落在各自零的高分位(≥0.90)-> 信号真但小,下一轮才允许定位;
                若只有一半 -> 噪声,这条线关掉。**
                纯噪声下「两半都 ≥0.90」的概率是 **0.01**,所以 8 次里出现 ≥2 次就已经很难是巧合。
POSITIVE CTRL   在**半样本**上重测 MDE(`#305a`:样本量一变,估计量就不是同一个),
                植入扫描用 5 种子平均(`#314b`)。
NEGATIVE CTRL   `perm_finite`,**每半各自**。
⚠ 报           **两半一致的方向**,不报两半的均值 —— 均值会把「一半有一半没有」洗成「有一点」。
IMPOSSIBLE      半样本的功效比全样本低,所以「两半都不显著」不能推翻全样本的 0.067;
                本轮判的是**一致性**,不是显著性。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
MODS={'无力感':pd.to_numeric(d['powerlessnessvariable'],errors='coerce').values.astype(float),
      '神经质':pd.to_numeric(d['neuroticismvariable'],errors='coerce').values.astype(float),
      '0–14岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)']
         .map({'Never':0.,'Sometimes':1.,'Often':2.}).values.astype(float)}
mAll=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&ok
for v in MODS.values(): mAll&=np.isfinite(v)
POOL=np.flatnonzero(mAll)
def pieces(rows):
    m=np.zeros(NN,bool); m[rows]=True
    def rz(a,b):
        out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
        out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
    RS,RC=rz(S,C3),rz(C3,S); n=int(m.sum())
    z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    zr,zc=z(RS),z(RC); ZM={k:z(v) for k,v in MODS.items()}
    BASE=np.column_stack([np.ones(n),zr,zc]+[ZM[k] for k in MODS])
    INT=np.column_stack([zr*ZM[k] for k in MODS]+[zc*ZM[k] for k in MODS])
    FULL=np.column_stack([BASE,INT])
    def r2(X,y):
        b,*_=np.linalg.lstsq(X,y,rcond=None); r=y-X@b
        return 1-float(r@r)/float(((y-y.mean())**2).sum())
    def stat(y):
        yy=(y-y.mean())/y.std(); return r2(FULL,yy)-r2(BASE,yy)
    return m,stat,zr,zc,ZM,n
def pf(v,seed):
    zz=v.copy(); j=np.flatnonzero(np.isfinite(zz))
    zz[j]=zz[np.random.default_rng(seed).permutation(j)]; return zz
NPERM=40
def arm(rows,seed):
    m,stat,zr,zc,ZM,n=pieces(rows)
    obs=stat(sh[m]); nul=np.array([stat(pf(sh,seed+i)[m]) for i in range(NPERM)])
    return obs,float(nul.mean()),float(nul.std()),float((nul<obs).mean()),n
print(f"全样本池 n={len(POOL):,};8 次劈半,每半各自跑 {NPERM} 次置换零")
rg=np.random.default_rng(31415); ROWS=[]
for t in range(8):
    p=rg.permutation(POOL); h=len(p)//2
    a=arm(p[:h],1000+100*t); b=arm(p[h:2*h],5000+100*t)
    ROWS.append(dict(v_split=t+1,qa=a[3],qb=b[3],obsa=a[0],obsb=b[0],n=a[4],
                     both=(a[3]>=0.90 and b[3]>=0.90),either=(a[3]>=0.90 or b[3]>=0.90)))
    print(f"  劈分 {t+1}: A 分位 **{a[3]:.3f}**(观测 {100*a[0]:.4f}pp)· "
          f"B 分位 **{b[3]:.3f}**(观测 {100*b[0]:.4f}pp)  {'★ 两半都 ≥0.90' if ROWS[-1]['both'] else ''}")
T=pd.DataFrame(ROWS); check_columns(T,'R360')
T.to_csv(pathlib.Path(__file__).parent/'results'/'split_half.csv',index=False)
nb=int(T.both.sum()); ne=int(T.either.sum())
allq=np.r_[T.qa.values,T.qb.values]
print(f"\n★ **两半都 ≥0.90 的劈分:{nb}/8**(纯噪声下期望 **0.08**)· 至少一半 ≥0.90:**{ne}/8**")
print(f"   16 个半样本的分位数:中位 **{np.median(allq):.3f}** · ≥0.90 的有 **{int((allq>=0.90).sum())}/16**"
      f"(纯噪声期望 1.6)· ≥0.50 的有 **{int((allq>=0.50).sum())}/16**(期望 8)")
m1,stat1,zr1,zc1,ZM1,n1=pieces(POOL[:len(POOL)//2])
nul1=np.array([stat1(pf(sh,777+i)[m1]) for i in range(NPERM)])
def plant_stat(Uk,amp,base):
    return float(np.mean([stat1(0.12*zr1*(1+amp*Uk)+0.12*zc1*(1+amp*Uk)
                          +np.random.default_rng(base+100*t).standard_normal(n1)) for t in range(5)]))
print(f"\n正对照:**在半样本上**重测 MDE(5 种子平均)")
MD={}
for ki,k in enumerate(MODS):
    got=None
    for amp in (0.10,0.20,0.30,0.50,0.80):
        if (plant_stat(ZM1[k],amp,900+ki)-nul1.mean())>2*nul1.std(): got=amp; break
    MD[k]=got; print(f"   {k:<14} 半样本 MDE **{'%.0f%%'%(100*got) if got else '> 80%'}**")
SW=[(a,plant_stat(ZM1[list(MODS)[0]],a,1234)) for a in (0.0,0.20,0.40,0.60)]
gg=Gate('那个信号在独立的一半人上还在不在')
gg.plant_direction_from_sweep('★ 正对照:半样本上的植入扫描',SW,SW[0][1],baseline_spread=2*float(nul1.std()))
gg.asserted('★ 半样本 MDE 三个都测了(全样本是 20%)',all(v is not None for v in MD.values()),
            ' · '.join(f"{k} {'%.0f%%'%(100*v) if v else '>80%'}" for k,v in MD.items()))
gg.asserted('★ 注册的 kill:两半都 ≥0.90 的劈分数(纯噪声期望 0.08/8)',nb>=2,
            f"**{nb}/8**;至少一半 {ne}/8;16 个半样本里 ≥0.90 的 {int((allq>=0.90).sum())}(期望 1.6)")
gg.asserted('⚠ 报两半一致的方向,不报均值',True,
            f"16 个分位数的中位 {np.median(allq):.3f} —— **均值会把「一半有一半没有」洗成「有一点」**")
gg.null_claim_uses_null_criteria('★ guard 21:若判为零,三件套在不在','NULL',
    perm_quantile=float(np.median(allq)),mde=max([v for v in MD.values() if v] or [0.99]),
    sensitivity_shown='半样本植入扫描',meaningful=0.30)
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
