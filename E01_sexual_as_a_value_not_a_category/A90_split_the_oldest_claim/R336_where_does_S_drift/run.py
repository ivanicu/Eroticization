import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A90 R336 -- S 的非不变性到底落在哪一格

`#289a`:去掉位置分 S,非不变性完全消失 ⇒ **S 是来源**。
`#290a`:S ↔ 羞耻在两组间没有可分辨差别 ⇒ **来源不在羞耻这一格**。
**⇒ 它在别的格上,而保留集 22 道里还有 21 格没查。**

ESTIMAND        对保留集每一格,算两组各自的 `corr(S, 结局)` 与它们的差;
                **每格各自的 offset**(同一格上「随机劈两组之差」的分布,因为每格 n 与方差不同);
                **全族阈值用最大统计量**(`#184` 同款),不逐格看 p。
KILL            **若集中在少数几格 -> 非不变性是局部的,可以指名道姓;
                若分散在很多格 -> 是 S 这个量本身在两组里的含义漂移,那更难处理也更重要。**
POSITIVE CTRL   沿用 `#290b` 修好的**交互项驱动**版本(`SEX × S`),**不能再用主效应**。
NEGATIVE CTRL   跨人置换 S(**只在有限值内**)。
⚠ 计数           按两层区间读(`#255c`)。
IMPOSSIBLE      「哪一格」是这份问卷的 22 个结局里的哪一格;S 的含义漂移可能落在没被问的地方。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S=np.where(ok,pos/np.maximum(cov,1),np.nan)
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
EX={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
GEN=[c for c,_ in OUT if any(k in str(c) for k in
     ('opposite gender','biological *female*','biological *male*','biological female','biological male'))]
SEEN=[c for c,_ in OUT if str(c).startswith('"I sometimes find people') or c=='animated']
HOLD=[(c,y) for c,y in OUT if c not in GEN and c not in SEEN]
print(f"保留集 **{len(HOLD)}** 道(`#333` 的同一套:剔除性别指涉 {len(GEN)} 道与已看过的 {len(SEEN)} 道)")
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
print(f"两组 n = {len(g0):,} / {len(g1):,}")
def corr_on(x,y,rows):
    m=np.zeros(NN,bool); m[rows]=True; m&=np.isfinite(x)&np.isfinite(y)
    return float(np.corrcoef(x[m],y[m])[0,1]) if m.sum()>250 else np.nan
def gaps(x,r0,r1,cells):
    return np.array([corr_on(x,y,r1)-corr_on(x,y,r0) for _,y in cells])
obs=gaps(S,g0,g1,HOLD)
rngS=np.random.default_rng(20260804)
NUL=[]
for t in range(80):
    p=rngS.permutation(np.flatnonzero(have))
    NUL.append(gaps(S,p[:len(g0)],p[len(g0):len(g0)+len(g1)],HOLD))
NUL=np.array(NUL)
per_sd=np.nanstd(NUL,axis=0)
fam=[float(np.nanquantile(np.nanmax(np.abs(NUL[rngS.choice(len(NUL),len(NUL),True)]),axis=1),0.95))
     for _ in range(12)]
THR=float(np.mean(fam)); THR_SD=float(np.std(fam))
hit=[(HOLD[i][0],obs[i],per_sd[i]) for i in range(len(HOLD)) if np.isfinite(obs[i]) and abs(obs[i])>THR]
print(f"\n**全族阈值(最大统计量,80 次随机劈)= {THR:.4f} ± {THR_SD:.4f}**")
print(f"**越阈的格数:{len(hit)}/{len(HOLD)}**")
for nm,v,sd in sorted(hit,key=lambda t:-abs(t[1])):
    print(f"   {nm[:46]:<48} 差 **{v:+.4f}**(该格 offset ±{sd:.4f})")
oi=np.argsort(-np.abs(np.nan_to_num(obs)))
print(f"   最大的三格(无论是否越阈):"+' · '.join(
    f"{HOLD[i][0][:20]} {obs[i]:+.3f}/±{per_sd[i]:.3f}" for i in oi[:3]))
# 越阈那一格的两组数值 —— 只有它们说得出是关于人的哪句话
for nm,v,sd in hit:
    y=dict(HOLD)[nm]
    a=corr_on(S,y,g0); b=corr_on(S,y,g1)
    print(f"   ⇒ `{nm[:44]}`:组 0 **{a:+.4f}** vs 组 1 **{b:+.4f}**(差 {b-a:+.4f})")
    # 顺带:同一格上「支配」的对偶「臣服」是否也漂移
sub=[c for c,_ in HOLD if 'submissive' in str(c)]
if sub:
    y=dict(HOLD)[sub[0]]
    print(f"   对偶格 `{sub[0][:44]}`:组 0 {corr_on(S,y,g0):+.4f} vs 组 1 {corr_on(S,y,g1):+.4f}"
          f"(差 {corr_on(S,y,g1)-corr_on(S,y,g0):+.4f},阈值 {THR:.4f})")

zS=np.where(np.isfinite(S),(S-np.nanmean(S))/np.nanstd(S),0.0)
n_=rngS.standard_normal(NN)
FK=np.where(have,0.35*np.nan_to_num(SEX)*zS+n_,np.nan)
HOLD_P=HOLD[:-1]+[('__plant__',FK)]
obs_p=gaps(S,g0,g1,HOLD_P)
print(f"\n正对照(交互项驱动的假结局放进保留集,`#290b` 修好的版本):"
      f"该格差 **{obs_p[-1]:+.4f}** vs 全族阈值 {THR:.4f} -> "
      f"**{'被抓到 ✅' if abs(obs_p[-1])>THR else '没抓到 ❌'}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul_cnt=[int(np.nansum(np.abs(gaps(perm_finite(S,400+i),g0,g1,HOLD))>THR)) for i in range(10)]
print(f"负对照(置换 S,只在有限值内):越阈格数 {np.mean(nul_cnt):.1f} ± {np.std(nul_cnt):.1f}")
T=pd.DataFrame([dict(outcome=HOLD[i][0][:46],gap=float(obs[i]),cell_sd=float(per_sd[i]),
                     over=bool(np.isfinite(obs[i]) and abs(obs[i])>THR)) for i in range(len(HOLD))])
check_columns(T,'R336'); T.to_csv(pathlib.Path(__file__).parent/'results'/'per_cell_gaps.csv',index=False)

g=Gate('S 的非不变性落在哪一格')
g.asserted('正对照:交互项驱动的假结局必须被全族阈值抓到',
           abs(obs_p[-1])>THR, f"假结局差 {obs_p[-1]:+.4f} vs 阈值 {THR:.4f}")
g.negative_control('置换 S 后的越阈格数',float(np.mean(nul_cnt)),float(len(hit)),
                   null_spread=float(np.std(nul_cnt)),
                   null_kind='跨人置换 S(只在有限值内)—— 只打掉配对,保留每格的边际')
g.count_needs_interval('越阈的格数',len(hit),len(HOLD),float(np.std(nul_cnt)),
                       'threshold_resample_阈值重抽样',n_resamples=12,seed_spread=THR_SD*20)
g.asserted('★ 注册的 kill:集中在少数几格 -> 局部的,可指名道姓;分散在很多格 -> S 的含义漂移',
           0<len(hit)<=4, f"{len(hit)}/{len(HOLD)} 越阈;"
           +(' · '.join(f"{n[:18]} {v:+.3f}" for n,v,_ in hit) if hit else '无'))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
