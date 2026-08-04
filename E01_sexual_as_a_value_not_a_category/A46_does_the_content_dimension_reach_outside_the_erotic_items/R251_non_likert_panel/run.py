import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A46 R251 -- 那个无名维度,够不够得着问卷情色部分之外

`#205` 的三块都齐了(是一个维度 · 载荷可复现 · 有具体选项表),**唯独名字还欠着**。
转到一个**不需要命名**的用途:**它能不能预测它自己之外的行为?**
`#188` 的六道结局**全都是同一份问卷里的情色 Likert 题**。

ESTIMAND        Cres(以及 S · rho_i 作对照)对 **11 个非 Likert、非情色**字段的相关,
                最大统计量零给全族阈值。
KILL            **若 Cres 在这些字段上一个都不越阈值 -> 它的作用域是"这份问卷的情色题",
                而这必须写进它的每一次陈述。**
POSITIVE CTRL   合成一个由 Cres 构造的结局塞进同一条面板 -> 必须强测到。
NEGATIVE CTRL   每个字段在分析样本内打乱(`#184b` 的教训)。
PRIOR           `#101`/`#102` 已测过**大五对 S**(全部 |r| ≤ 0.056)——
                **但从没测过大五对内容维度**。这是本轮的新东西。
编码(跑之前写死)
                age 5 段 -> 段中点 · 大五与无力感:已是数值
                关系风格 Monogamous=0 / Not monogamous=1
                0–14 岁被打屁股 Never=0 / Sometimes=1 / Often=2
                自评吸引力 7 点 -> −3..+3 · 成长期性开放度 Repressed=−1 / Neutral=0 / Liberated=+1
                ⚠ 「Which describes you best?」是**支配/被支配自我认同** -> **属情色,排除**。
IMPOSSIBLE      11 个字段的编码由我写(不是独立编码者);
                但它们是**人口学与既有量表**,不是我为这条轴发明的分类 —— 这与 `#203c` 的污染问题不同。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
_,RHO=betas(V)
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
FIELDS={
 'age(段中点)':d['age'].map(AGE),
 'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格(非一对一=1)':d['Personally, your preferred relationship style is: (4jib23m)']
        .map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)']
        .map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)']
        .map({'Repressed':-1,'Neutral':0,'Liberated':1}),
}
print(f"字段 {len(FIELDS)} 个;覆盖 " + ' · '.join(f"{k.split('(')[0]}:{int(v.notna().sum())//1000}k"
                                                 for k,v in list(FIELDS.items())[:5]))

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(d); con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    con[ppl]+=Z@v[:,-1]; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
Cb=np.where(ok,con/np.maximum(cnt,1),np.nan); Sb=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
bi=np.flatnonzero(np.isfinite(Sb)&np.isfinite(Cb)&np.isfinite(KB))
X0=np.c_[np.ones(len(bi)),Sb[bi]]
Cr=np.full(NN,np.nan); Cr[bi]=Cb[bi]-X0@np.linalg.lstsq(X0,Cb[bi],rcond=None)[0]
check_residualized(Cr[bi],Sb[bi],'R251 内容残差')
rng=np.random.default_rng(20260803)
FIELDS['【正对照】由 Cres 造的合成结局']=pd.Series(np.where(np.isfinite(Cr),
        Cr+rng.standard_normal(NN)*np.nanstd(Cr[bi])*2,np.nan))
P={'Cres(内容)':Cr,'S(位置)':Sb,'rho_i(何时)':RHO}

def cr(y,x,ii):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    X=np.c_[np.ones(len(jj)),KB[jj]]
    ry=y[jj]-X@np.linalg.lstsq(X,y[jj],rcond=None)[0]
    rx=x[jj]-X@np.linalg.lstsq(X,x[jj],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), len(jj)

rows=[]; nulls={k:[] for k in P}
for name,ser in FIELDS.items():
    y=ser.values.astype(float)
    rec=dict(field=name)
    for k,x in P.items():
        r,n=cr(y,x,bi); rec[k]=r; rec['n']=n
        ps=[]
        for _ in range(40):
            yp=y.copy(); yp[bi]=rng.permutation(y[bi]); v,_=cr(yp,x,bi)
            if np.isfinite(v): ps.append(abs(v))
        if len(ps)>=20 and '正对照' not in name: nulls[k].append(ps)
    rows.append(rec)
T=pd.DataFrame(rows); check_columns(T,'R251')
check_coverage(len(T),len(FIELDS),'R251 面板',tol=0.0)
thr={}
for k in P:
    L=min(len(x) for x in nulls[k])
    thr[k]=float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nulls[k]]),axis=0),0.95))
T.to_csv(pathlib.Path(__file__).parent/'results'/'nonlikert.csv',index=False)
print(f"\n全族阈值:" + ' · '.join(f"{k} {thr[k]:.4f}" for k in P) + "\n")
print(f"{'Cres':>10}{'S':>10}{'rho_i':>10}{'n':>8}  字段")
for _,r in T.iterrows():
    m=lambda k:'★' if abs(r[k])>thr[k] else ' '
    print(f"{r['Cres(内容)']:>+9.4f}{m('Cres(内容)')}{r['S(位置)']:>+9.4f}{m('S(位置)')}"
          f"{r['rho_i(何时)']:>+9.4f}{m('rho_i(何时)')}{int(r.n):>8,}  {r.field}")
real=T[~T.field.str.contains('正对照')]
nC=int((real['Cres(内容)'].abs()>thr['Cres(内容)']).sum())
nS=int((real['S(位置)'].abs()>thr['S(位置)']).sum())
nR=int((real['rho_i(何时)'].abs()>thr['rho_i(何时)']).sum())
print(f"\n越阈值:Cres {nC}/{len(real)} · S {nS}/{len(real)} · rho_i {nR}/{len(real)}")
pc=T[T.field.str.contains('正对照')].iloc[0]
g=Gate('内容维度够不够得着情色题之外')
g.asserted('正对照:由 Cres 造的合成结局必须强测到',abs(pc['Cres(内容)'])>0.3,f"{pc['Cres(内容)']:+.4f}")
g.asserted('先验对照:大五对 S 应当很小(`#101`/`#102` 报的是 |r| ≤ 0.056)',
           float(real[real.field.isin(['openness','conscientiousness','extroversion','neuroticism','agreeableness'])]
                 ['S(位置)'].abs().max())<0.09,
           f"实测最大 {float(real[real.field.isin(['openness','conscientiousness','extroversion','neuroticism','agreeableness'])]['S(位置)'].abs().max()):.4f}")
g.asserted('注册的 kill:Cres 一个都不越阈值 -> 作用域只是这份问卷的情色题',nC==0,
           f"Cres 越阈值 {nC}/{len(real)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
