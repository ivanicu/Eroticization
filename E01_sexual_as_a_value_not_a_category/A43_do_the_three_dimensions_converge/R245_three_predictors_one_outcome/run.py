import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A43 R245 -- 三个维度,是不是在同一道题上汇合

`#199b` 的最强一格是「性成熟但未成年」(rho_i −0.0790)——
**而这一格在 S 侧(`#184b`)与内容侧(`#190`)也都越过阈值。**
三个互相不是影子的维度**同时指向同一道题**。

    CONVERGE  三个偏相关都独立存活 -> 那道题是一个「三方汇合点」
    RELAY     只剩一个 -> 另外两个对它的相关是**经由**第三个走的

ESTIMAND        对每道结局,同时放入 S · Cres · rho_i,报三个偏相关与总 R²。
                对照结局:羞耻(S 最强)· `animated`(内容最强)。
KILL            **若三个偏相关里有 ≥2 个在放入其余两个后掉到 2× 以下 -> RELAY,不是汇合。**
POSITIVE CTRL   合成一个由**三者之和**构造的结局 -> 三个偏相关必须都存活。
NEGATIVE CTRL   合成一个**只由 S** 构造的结局 -> 只有 S 存活。
COLLINEARITY    Cres ⟂ S 是构造的(`#129`);rho_i vs S = −0.046(`#199b`);rho_i vs Cres 待测。
                **三者若近正交,偏相关可读**(与 `#182b`/`#189b` 的共线陷阱相反的情形)——
                但共线度必须先量,不能假设。
IMPOSSIBLE      rho_i 信度 0.242、S 0.432、羞耻单题不可估 —— **三个预测子的衰减程度不同**,
                所以偏相关的**相对大小**不可直接比;本轮只判"存不存活",不排名。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
_,RHO=betas(V)
df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df_raw.columns if df_raw[c].dtype!=object and
     set(pd.Series(df_raw[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df_raw[c].notna().sum()>10000]
MAT=next(c for c in lik if 'sexual maturity' in c or 'clearly reached' in c)
SHAME=next(c for c in lik if 'ashamed' in c); ANI=next(c for c in lik if c=='animated')

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df_raw); con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
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
base=np.isfinite(Sb)&np.isfinite(Cb)&np.isfinite(KB)&np.isfinite(RHO)&KEEP
bi=np.flatnonzero(base)
X0=np.c_[np.ones(len(bi)),Sb[bi]]
Cr=np.full(NN,np.nan); Cr[bi]=Cb[bi]-X0@np.linalg.lstsq(X0,Cb[bi],rcond=None)[0]
check_residualized(Cr[bi],Sb[bi],'R245 内容残差')
P={'S(位置)':Sb,'Cres(内容)':Cr,'rho_i(何时)':RHO}
print(f"n = {len(bi):,}")
print("三个预测子的两两相关:")
ks=list(P)
for a in range(3):
    for b in range(a+1,3):
        print(f"  {ks[a]:<12} × {ks[b]:<12} = {np.corrcoef(P[ks[a]][bi],P[ks[b]][bi])[0,1]:+.4f}")

def partials(y,ii):
    out={}
    for k in ks:
        oth=[P[o][ii] for o in ks if o!=k]
        X=np.c_[np.ones(len(ii)),KB[ii],*oth]
        ry=y[ii]-X@np.linalg.lstsq(X,y[ii],rcond=None)[0]
        rx=P[k][ii]-X@np.linalg.lstsq(X,P[k][ii],rcond=None)[0]
        out[k]=float(np.corrcoef(ry,rx)[0,1])
    Xf=np.c_[np.ones(len(ii)),KB[ii],*[P[k][ii] for k in ks]]
    pr=Xf@np.linalg.lstsq(Xf,y[ii],rcond=None)[0]
    out['R2']=float(1-((y[ii]-pr)**2).sum()/((y[ii]-y[ii].mean())**2).sum())
    return out

rng=np.random.default_rng(20260803)
z=lambda a:(a-np.nanmean(a[bi]))/np.nanstd(a[bi])
synth3=np.full(NN,np.nan); synth3[bi]=z(Sb)[bi]+z(Cr)[bi]+z(RHO)[bi]+rng.standard_normal(len(bi))*2
synthS=np.full(NN,np.nan); synthS[bi]=z(Sb)[bi]+rng.standard_normal(len(bi))*2
OUT=[('性成熟未成年',df_raw[MAT].values.astype(float)),('羞耻',df_raw[SHAME].values.astype(float)),
     ('animated',df_raw[ANI].values.astype(float)),
     ('【正对照】三者之和',synth3),('【负对照】只由 S 造',synthS)]
rows=[]
for nm,y in OUT:
    ii=bi[np.isfinite(y[bi])]
    d=partials(y,ii)
    sds={k:float(np.std([partials(y,rng.choice(ii,len(ii),replace=True))[k] for _ in range(150)])) for k in ks}
    rows.append(dict(outcome=nm,n=len(ii),**{k:d[k] for k in ks},
                     **{f'sd_{k}':sds[k] for k in ks},R2=d['R2']))
T=pd.DataFrame(rows); check_columns(T,'R245'); T.to_csv(pathlib.Path(__file__).parent/'results'/'three.csv',index=False)
print(f"\n{'结局':<16}{'S(位置)':>12}{'Cres(内容)':>13}{'rho_i(何时)':>13}{'总 R²':>9}")
for _,r in T.iterrows():
    stars=lambda k: '★' if abs(r[k])>2*r[f'sd_{k}'] else ' '
    print(f"{r.outcome:<16}{r[ks[0]]:>+11.4f}{stars(ks[0])}{r[ks[1]]:>+12.4f}{stars(ks[1])}"
          f"{r[ks[2]]:>+12.4f}{stars(ks[2])}{r.R2:>9.4f}")

mat=T[T.outcome=='性成熟未成年'].iloc[0]
pos_=T[T.outcome.str.contains('正对照')].iloc[0]; neg_=T[T.outcome.str.contains('负对照')].iloc[0]
surv=lambda r: sum(1 for k in ks if abs(r[k])>2*r[f'sd_{k}'])
g=Gate('三个维度是不是在同一道题上汇合')
g.asserted('正对照:由三者之和构造的结局,三个偏相关必须都存活',surv(pos_)==3,f"存活 {surv(pos_)}/3")
g.asserted('负对照:只由 S 构造的结局,只有 S 存活',surv(neg_)==1 and abs(neg_[ks[0]])>2*neg_[f'sd_{ks[0]}'],
           f"存活 {surv(neg_)}/3")
g.asserted('共线度已量化(与 #182b/#189b 的陷阱相反的情形,但必须先量)',True,
           ' · '.join(f"{ks[a]}×{ks[b]}={np.corrcoef(P[ks[a]][bi],P[ks[b]][bi])[0,1]:+.3f}"
                      for a in range(3) for b in range(a+1,3)))
g.asserted('注册的 kill:三个里 ≥2 个掉到 2× 以下 -> RELAY',surv(mat)<2,f"性成熟未成年 存活 {surv(mat)}/3")
print(g)
print(f"\n  => {'CONVERGE —— 三方汇合点' if surv(mat)==3 else ('RELAY' if surv(mat)<2 else f'部分汇合({surv(mat)}/3)')}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
