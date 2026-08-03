import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Confound: the 4 role-varying blocks are all bodily fluids. "Fluids are special" predicts the
same result as "roles are what the axis reads".
Separator: inside the SAME block, same people, same n, build a MATCHED PLACEBO contrast from
the non-role options. If POWER predicts the role contrast but not the placebo contrast, the
fluid explanation dies. If it predicts both, the role explanation dies.
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(101)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
A=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in A.columns if any(k in c for k in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=pd.concat([z(pd.to_numeric(A[c],errors='coerce'))*sg[c] for c in pc],axis=1).mean(axis=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
COV=df[['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable']].apply(pd.to_numeric,errors='coerce')
COV=COV.fillna(COV.median()); COV=((COV-COV.mean())/(COV.std()+1e-9)).fillna(0.)

FLUID={8:'saliva',7:'precum',9:'squirt',83:'male ejaculate',11:'urine',6:'breast milk',10:'sweat'}
SELF=r'(myself|\bmy\b)'; OTH=r'(others|other )'
out=[]
for qi,name in FLUID.items():
    sub=lg[lg.qi==qi]
    if len(sub)==0: continue
    ppl=np.array(sorted(sub.person.unique())); opt=np.array(sorted(sub.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower()
    s_i=np.flatnonzero(lo.str.contains(SELF,regex=True).values)
    o_i=np.flatnonzero(lo.str.contains(OTH,regex=True).values)
    n_i=np.flatnonzero(~(lo.str.contains(SELF,regex=True)|lo.str.contains(OTH,regex=True)).values)
    if len(s_i)<2 or len(o_i)<2 or len(n_i)<2: continue
    p_=POWER.reindex(ppl).values; ok=~np.isnan(p_)
    if ok.sum()<400: continue
    M,p_,ppl_ok=M[ok],p_[ok],ppl[ok]
    role = M[:,s_i].mean(1)-M[:,o_i].mean(1)                       # self-minus-other
    # matched placebo: split the NON-role options into two halves, same arity as role contrast
    plac=[]
    for _ in range(200):
        pm=rng.permutation(n_i); h=len(pm)//2
        plac.append(M[:,pm[:h]].mean(1)-M[:,pm[h:2*h]].mean(1))
    prop=M.mean(1)
    D=np.c_[np.ones(ok.sum()),COV.loc[ppl_ok].values,prop]
    def part(y):
        b,*_=lstsq(D,y,rcond=None); yr=y-D@b
        b2,*_=lstsq(D,p_,rcond=None); pr=p_-D@b2
        return np.corrcoef(yr,pr)[0,1]
    r_role=part(role); r_plac=np.array([part(v) for v in plac])
    out.append(dict(block=name,n=int(ok.sum()),n_self=len(s_i),n_other=len(o_i),n_neutral=len(n_i),
                    r_role=round(r_role,3), plac_mean=round(float(np.mean(r_plac)),3),
                    plac_p95=round(float(np.percentile(np.abs(r_plac),95)),3),
                    ratio=round(abs(r_role)/max(float(np.percentile(np.abs(r_plac),95)),1e-6),2)))
R=pd.DataFrame(out); R.to_csv('data/derived/placebo_contrast.csv',index=False)
print(R.to_string(index=False))
print(f"\n  r_role  = corr(POWER, self-minus-other endorsement), demographics+propensity removed")
print(f"  plac_p95= 95th pct of |r| for 200 MATCHED contrasts built from the SAME block's non-role options")
print(f"\n  blocks where |r_role| exceeds its own placebo 95th pct: {int((R.ratio>1).sum())}/{len(R)}")
print(f"  median r_role = {R.r_role.median():+.3f}   median placebo p95 = {R.plac_p95.median():.3f}   median ratio = {R.ratio.median():.2f}x")
