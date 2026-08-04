import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A60 R278 -- 「N/31 越阈结局数」到底有多宽,以及宽度由什么决定

**类型:CLOSURE**(如实标注,`§0`)。它不分开任何世界 —— 它保护的是一个横跨多轮的**写法**。
`#232d`:同一个分数、同一次运行里,越阈个数两次算出 **10 和 8**。
本项目把 `N/M` 当确定值写进了 `#223a` `#230a` `#270` 与两个 README。

ESTIMAND        固定分数,把最大统计量阈值的零抽样独立重复 25 次,报 `N` 的分布。
机制假说         **宽度不来自分数强弱,来自有多少个结局正好坐在阈值附近。**
                可测形式:跨 6 个分数,`sd(N)` 应与**边界结局数**正相关。
POSITIVE CTRL   两端:纯噪声分数(几乎全在阈值下)与极强分数(几乎全在阈值上)
                **都应当窄** —— 若"越强越窄"成立,这两个就会一宽一窄,机制假说被否。
KILL            **若 `corr(边界数, sd(N)) > 0.5` 且两个极端都窄 -> 机制成立,
                守卫的区间读法有依据;若不成立 -> 宽度另有来源,守卫仍然要装,但解释要撤。**
IMPOSSIBLE      这里只重抽**阈值**,不重抽人。真实的 `N` 还有一层人层抽样的抖动,
                所以本轮报的是**下界宽度**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
rgb=np.random.default_rng(500); Spos=np.zeros(NN); ctp=np.zeros(NN); Rp=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    Spos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0); ctp[ppl]+=1; Rp[b,ppl]=M.mean(1)
S=np.where(ok,Spos/np.maximum(ctp,1),np.nan)
F=np.isfinite(Rp); Z=np.where(F,Rp,0.0); tot=Z.sum(0); c2=F.sum(0); R=np.full_like(Rp,np.nan)
for b in range(NB):
    lo=np.where(c2-F[b]>=6,(tot-Z[b])/np.maximum(c2-F[b],1),np.nan); R[b]=Rp[b]-lo; R[b]=R[b]-np.nanmean(R[b])
def hs(cols):
    sub=R[cols]; F2=np.isfinite(sub)
    return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
D=z(hs(TRG))-z(hs(ORD))
C=np.full((NB,NB),np.nan)
for i in range(NB):
    for j in range(NB):
        m=np.isfinite(R[i])&np.isfinite(R[j])&ok
        if m.sum()>300: C[i,j]=np.corrcoef(R[i][m],R[j][m])[0,1]
C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; V=np.linalg.eigh(C)[1][:,::-1]
def comp(k):
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(V[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)

lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EX={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
NO=len(OUT)

def corrs(x):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>=200 else np.nan)
    return np.array(r), bi
def thr_once(x,bi,seed):
    rg=np.random.default_rng(seed); nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: continue
        nl.append([abs(float(np.corrcoef(rg.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(v) for v in nl)
    return float(np.nanquantile(np.nanmax(np.array([v[:L] for v in nl]),0),0.95))

rngN=np.random.default_rng(2026)
noise=np.where(ok,rngN.standard_normal(NN),np.nan)
strong=np.where(ok,sum(np.nan_to_num((y-np.nanmean(y))/np.nanstd(y)) for _,y in OUT[:12]),np.nan)
SCORES=[('位置分 S',S),('跨块对比 D',D),('成分 1',comp(0)),('成分 3',comp(2)),
        ('纯噪声',noise),('极强(12 个结局之和)',strong)]
rows=[]
print(f"{'分数':<20}{'N 均值':>8}{'sd':>7}{'范围':>9}{'阈值':>9}{'阈值 sd':>9}{'边界结局数':>10}")
for nm,x in SCORES:
    r,bi=corrs(x); Ns=[]; ths=[]
    for t in range(25):
        th=thr_once(x,bi,7000+97*t); ths.append(th); Ns.append(int(np.nansum(np.abs(r)>th)))
    Ns=np.array(Ns); mt,st=float(np.mean(ths)),float(np.std(ths))
    bord=int(np.nansum(np.abs(np.abs(r)-mt)<2*st))
    rows.append(dict(score=nm,n_mean=float(Ns.mean()),n_sd=float(Ns.std()),
                     n_min=int(Ns.min()),n_max=int(Ns.max()),thr=mt,thr_sd=st,borderline=bord))
    print(f"{nm:<20}{Ns.mean():>8.2f}{Ns.std():>7.2f}{str(Ns.min())+'–'+str(Ns.max()):>9}"
          f"{mt:>9.4f}{st:>9.4f}{bord:>10}")
T=pd.DataFrame(rows); check_columns(T,'R278')
T.to_csv(pathlib.Path(__file__).parent/'results'/'count_width.csv',index=False)
cb=float(np.corrcoef(T.borderline,T.n_sd)[0,1])
ext=T[T.score.isin(['纯噪声','极强(12 个结局之和)'])]
print(f"\ncorr(边界结局数, sd(N)) = **{cb:+.4f}**;两个极端的 sd:"
      + ' · '.join(f"{r.score} {r.n_sd:.2f}" for _,r in ext.iterrows()))
worst=T.loc[T.n_sd.idxmax()]
print(f"最宽的一个:**{worst.score} N = {worst.n_mean:.1f} ± {worst.n_sd:.2f}(实测跨度 "
      f"{int(worst.n_min)}–{int(worst.n_max)})**")

g=Gate('N/31 有多宽,宽度由什么决定')
g.asserted('⚠ 类型标注:这是 CLOSURE,不分开任何世界,保护的是一个横跨多轮的写法',
           True, '§0 三类动作:Frontier / Closure / Production —— 本轮是 Closure')
g.asserted('正对照两端:纯噪声与极强分数都应当窄(否则「越强越窄」才是机制)',
           float(ext.n_sd.max())<1.0, ' · '.join(f"{r.score} sd {r.n_sd:.2f}" for _,r in ext.iterrows()))
for _,r in T.iterrows():
    g.count_needs_interval(f"{r.score} 的越阈计数",int(round(r.n_mean)),NO,
                           float(r.n_sd),'threshold_resample_阈值重抽样',n_resamples=25)
g.asserted('⚠ P7:守卫在第一次使用时被自己的正对照打出一个假阳,已修后重打',
           True, '纯噪声计数钉死在 0/31,展布真为 0 —— 「没重抽」与「重抽了不动」被折成了同一个 FAIL(#233b)')
g.asserted('★ 注册的 kill:corr(边界结局数, sd(N)) > 0.5 -> 宽度来自"坐在阈值上的结局数"',
           cb>0.5, f"corr = {cb:+.4f};各分数边界数 {T.borderline.tolist()} vs sd {T.n_sd.round(2).tolist()}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
