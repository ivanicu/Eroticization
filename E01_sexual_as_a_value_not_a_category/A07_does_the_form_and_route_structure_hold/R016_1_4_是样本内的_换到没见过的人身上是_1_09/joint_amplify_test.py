import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A107 R359 -- 有没有东西**放大**这份羞耻

⚠ **候选集与判据已在跑之前单独提交**(见 `PREREGISTRATION.md` 的 R359 段)。
三个**先验**候选:**无力感 · 神经质 · 0–14 岁被打屁股**。
⚠ **污染声明**:我看过 `#309c` 里神经质那两格的数值 -> **神经质的单变量读数按污染处理**。
⚠ **本轮不报「哪一个」** —— 联合检验说不出,而事后挑最大的那个正是 `#309c` 的陷阱。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
NEW={'无力感':pd.to_numeric(d['powerlessnessvariable'],errors='coerce').values.astype(float),
     '神经质':pd.to_numeric(d['neuroticismvariable'],errors='coerce').values.astype(float),
     '0–14岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)']
        .map({'Never':0.,'Sometimes':1.,'Often':2.}).values.astype(float)}
OLD={'成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map(
        {'Repressed':-1.,'Neutral':0.,'Liberated':1.}).values.astype(float),
     '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map(
        {'Monogamous':0.,'Not monogamous':1.}).values.astype(float),
     '年龄':d['age'].map(AGE).values.astype(float),
     '开放性':pd.to_numeric(d['opennessvariable'],errors='coerce').values.astype(float)}
def run_set(MODS,tag):
    m=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&ok
    for v in MODS.values(): m&=np.isfinite(v)
    def resid(a,b):
        out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
        out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
    RS=resid(S,C3); RC=resid(C3,S); n=int(m.sum())
    z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    zr,zc=z(RS),z(RC); ZM={k:z(v) for k,v in MODS.items()}
    BASE=np.column_stack([np.ones(n),zr,zc]+[ZM[k] for k in MODS])
    INT=np.column_stack([zr*ZM[k] for k in MODS]+[zc*ZM[k] for k in MODS])
    def r2(X,y):
        b,*_=np.linalg.lstsq(X,y,rcond=None); r=y-X@b
        return 1-float(r@r)/float(((y-y.mean())**2).sum())
    def stat(y):
        yy=(y-y.mean())/y.std(); return r2(np.column_stack([BASE,INT]),yy)-r2(BASE,yy)
    def pf(v,seed):
        zz=v.copy(); j=np.flatnonzero(np.isfinite(zz))
        zz[j]=zz[np.random.default_rng(seed).permutation(j)]; return zz
    obs=stat(sh[m]); nul=np.array([stat(pf(sh,600+i)[m]) for i in range(60)])
    q=float((nul>=obs).mean())
    print(f"\n【{tag}】 n={n:,} · {INT.shape[1]} 个交互项")
    print(f"   观测 R² 下降 **{100*obs:.4f}pp** · 置换零 **{100*nul.mean():.4f} ± {100*nul.std():.4f}pp** "
          f"-> **{(obs-nul.mean())/max(2*nul.std(),1e-12):+.2f}× 的 2×展布** · 零里 ≥ 观测 **{100*q:.1f}%**")
    # ⚠ 第一版每个强度只用**一个**噪声种子 -> 扫描非单调 -> 正对照 FAIL -> MDE 不可读。
    #    非单调的来源是噪声,不是种入 —— 所以每个强度**平均 5 个种子**。
    NSEED=5
    def plant_stat(Uk,amp,base_seed):
        return float(np.mean([stat(0.12*zr*(1+amp*Uk)+0.12*zc*(1+amp*Uk)
                              +np.random.default_rng(base_seed+100*t).standard_normal(n))
                              for t in range(NSEED)]))
    MDEV={}
    for ki,k in enumerate(MODS):
        Uk=ZM[k]; got=None
        for amp in (0.10,0.15,0.20,0.30,0.50):
            if (plant_stat(Uk,amp,55+ki)-nul.mean())>2*nul.std(): got=amp; break
        MDEV[k]=got
        print(f"     逐变量 MDE · {k:<14} **{'%.0f%%'%(100*got) if got else '> 50%'}**({NSEED} 种子均值)")
    SW=[]
    Uk=ZM[list(MODS)[0]]
    for amp in (0.0,0.10,0.20,0.30,0.50): SW.append((amp,plant_stat(Uk,amp,66)))
    return dict(obs=obs,nmean=float(nul.mean()),nsd=float(nul.std()),q=q,MDEV=MDEV,SW=SW,n=n,k=INT.shape[1])
A3=run_set(NEW,'★ 三个先验候选:无力感 · 神经质 · 被打屁股')
A7=run_set({**OLD,**NEW},'次要:七个变量全放(4 旧 + 3 新)')
T=pd.DataFrame([dict(v_set='三个新候选',**{k:v for k,v in A3.items() if k in('obs','nmean','nsd','q','n','k')}),
                dict(v_set='七个全放',**{k:v for k,v in A7.items() if k in('obs','nmean','nsd','q','n','k')})])
check_columns(T,'R359'); T.to_csv(pathlib.Path(__file__).parent/'results'/'amplify.csv',index=False)
worst=max([v for v in A3['MDEV'].values() if v] or [0.99])
gg=Gate('有没有东西放大这份羞耻')
gg.plant_direction_from_sweep('★ 正对照:单变量放大强度扫描',A3['SW'],A3['SW'][0][1],
                              baseline_spread=2*A3['nsd'])
gg.asserted('★ 逐变量 MDE:三个都测过了吗',all(v is not None for v in A3['MDEV'].values()),
            ' · '.join(f"{k} {'%.0f%%'%(100*v) if v else '>50%'}" for k,v in A3['MDEV'].items()))
gg.asserted('★ 注册的 kill:联合 R² 下降是否超过置换零的 2×展布',
            (A3['obs']-A3['nmean'])>2*A3['nsd'],
            f"观测 {100*A3['obs']:.4f}pp vs 零 {100*A3['nmean']:.4f} ± {100*A3['nsd']:.4f}pp "
            f"(零里 ≥ 观测 {100*A3['q']:.1f}%)")
gg.null_claim_uses_null_criteria('★ guard 21:这个零可不可发布','NULL',
    perm_quantile=A3['q'],mde=worst,sensitivity_shown=f"逐变量植入 {int(100*worst)}% 抓到",
    meaningful=0.30)
gg.asserted('⚠ 污染:神经质的单变量读数不可作为独立证据(`PREREGISTRATION.md` 跑前声明)',True,
            '我看过 `#309c` 里神经质那两格;它进候选集的理由是先验的,单独读数按污染处理')
gg.asserted('⚠ 本轮不报「哪一个」',True,'联合检验说不出;事后挑最大的那个正是 `#309c` 的陷阱')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
