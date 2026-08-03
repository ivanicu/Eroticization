import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The agent/patient indicators sit near each other in the survey, so their intercorrelation
could be common-method inflation. Decisive test: do the axes predict WHICH OPTIONS a person
endorses inside blocks they were never measured on -- held out, cross-domain, and against
the propensity baseline that must be beaten.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(31)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
A=pd.read_csv('data/derived/agent_patient.csv')

# two candidate axes, built ONLY from the direct questions
power_cols=[c for c in A.columns if any(k in c for k in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
gaze_cols =[c for c in A.columns if 'exhibition' in c or 'voyeur' in c]
def zscore(s): return (s-s.mean())/(s.std()+1e-9)
sign={c:(-1 if 'worship' in c else 1) for c in power_cols}
POWER=pd.concat([zscore(pd.to_numeric(A[c],errors='coerce'))*sign[c] for c in power_cols],axis=1).mean(axis=1)
GAZE =zscore(pd.to_numeric(A[[c for c in gaze_cols if 'exhibition' in c][0]],errors='coerce')) \
     -zscore(pd.to_numeric(A[[c for c in gaze_cols if 'voyeur' in c][0]],errors='coerce'))
print(f"POWER axis from {len(power_cols)} indicators, n={POWER.notna().sum()}")
print(f"GAZE  axis from 2 indicators,            n={GAZE.notna().sum()}")
print(f"r(POWER,GAZE) = {POWER.corr(GAZE):+.3f}   <- near zero => two axes, not one\n")

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
cc=['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable']
COV=df[cc].apply(pd.to_numeric,errors='coerce'); COV=COV.fillna(COV.median()); COV=((COV-COV.mean())/(COV.std()+1e-9)).fillna(0.)

keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
res=[]
for _,q in keep.iterrows():
    sub=lg[lg.qi==q.qi]; ppl=np.array(sorted(sub.person.unique())); opt=np.array(sorted(sub.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    if len(ppl)<1200: continue
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)     # profile shape only
    p_=POWER.reindex(ppl).values; g_=GAZE.reindex(ppl).values
    prop=(M.sum(1)/M.shape[1])
    ok=~(np.isnan(p_)|np.isnan(g_))
    if ok.sum()<800: continue
    R,p_,g_,prop,ppl_ok=R[ok],p_[ok],g_[ok],prop[ok],ppl[ok]
    Xc=np.c_[np.ones(ok.sum()),COV.loc[ppl_ok].values,prop]
    idx=rng.permutation(ok.sum()); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
    def held_out_r2(X):
        b,*_=lstsq(X[tr],R[tr],rcond=None); P=X[te]@b
        return 1-((R[te]-P)**2).sum()/ (R[te]**2).sum()
    base=held_out_r2(Xc)
    full=held_out_r2(np.c_[Xc,p_,g_])
    pwr =held_out_r2(np.c_[Xc,p_])
    gz  =held_out_r2(np.c_[Xc,g_])
    nul =np.mean([held_out_r2(np.c_[Xc,rng.permutation(p_),rng.permutation(g_)]) for _ in range(10)])
    res.append(dict(qi=q.qi,col=q.col[:46],n=int(ok.sum()),
                    base=base,d_power=pwr-base,d_gaze=gz-base,d_both=full-base,d_null=nul-base))
T=pd.DataFrame(res).sort_values('d_both',ascending=False)
T.to_csv('data/derived/axis_prediction.csv',index=False)
print("held-out R2 GAIN over (demographics + propensity) baseline, per block:\n")
print(T.assign(**{k:T[k].round(4) for k in ['base','d_power','d_gaze','d_both','d_null']}).to_string(index=False))
print(f"\nblocks: {len(T)}")
print(f"median gain from POWER axis : {T.d_power.median():+.4f}")
print(f"median gain from GAZE  axis : {T.d_gaze.median():+.4f}")
print(f"median gain from BOTH       : {T.d_both.median():+.4f}")
print(f"median PERMUTED-axis gain   : {T.d_null.median():+.4f}   <- the floor")
print(f"blocks where both-axis gain exceeds 5x the permuted floor: {int((T.d_both> 5*abs(T.d_null)).sum())}/{len(T)}")
