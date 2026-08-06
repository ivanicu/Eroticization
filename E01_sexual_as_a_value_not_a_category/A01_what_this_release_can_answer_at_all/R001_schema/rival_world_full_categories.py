"""
E01 A03 R11 -- re-run the CONFIRMED rival-world result on the categories it never chose.

#31 found that A03's rounds obtain their data by exec-ing the RSA round's loader, which filters
onset columns down to the 27 having a matching AROUSAL-RATING column. That filter exists because
the RSA needed to PAIR onset with preference. Mapping the exec graph shows six A03 rounds inherit
it, including R09 -- the synthetic-rival-world round that produced one of only two CONFIRMED claims
in the whole project ("onset carries structure preference does not", observed 0.959 vs a rival
world at 0.441).

R09 needs the pairing for its residualisation, so the filter is defensible there -- but it was never
CHOSEN, and #31 showed the same accidental filter moved an unrelated result by 0.06. A CONFIRMED
claim resting on an unexamined inherited filter is exactly what should be re-run.

ESTIMAND        the top residual eigenvalue against its purpose-built rival world, computed on
                every category set the question admits.
IDENTIFICATION  identified at each cell; the rival world is regenerated inside each cell so the
                comparison never crosses category sets.
WORLDS          A  the CONFIRMED result is about the phenomenon: it holds on 27 and on the wider set
                B  it was about the filter: it weakens or vanishes when the filter is lifted
KILL (CONDITIONAL -- threshold evaluated ONLY if the gate passes)
      gate: in EVERY cell the injected positive control must be detected and the zero-injection
            control must NOT be
      then: observed above the rival world's 97.5th pct in every cell -> CONFIRMED, filter-free
            observed inside the rival world in any cell               -> the confirmation was the filter
            otherwise                                                  -> UNVERIFIED
POSITIVE CTRL   person-level rank-2 injection at amplitude 0.3 must be detected in every cell;
                amplitude 0.0 must not.
NEGATIVE CTRL   the rival world itself -- onset generated AS preference-covariance plus noise.
NOISE FLOOR     sd of the rival-world distribution per cell.
MULTIPLICITY    3 category sets x 3 noise levels x 3 seeds, all reported.
SEEDS           3.
IMPOSSIBLE      categories with no arousal-rating twin cannot enter the RESIDUALISATION, so the
                widest admissible set is bounded by the question, not by my choice. That bound is
                itself reported.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, eigvalsh
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,'15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
MIDS=np.array([2,5.5,7.5,9.5,11.5,13.5,15.5,17.5,22,28]); EDGES=[0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,25.5,99]
binify=lambda x: MIDS[np.clip(np.digitize(x,EDGES)-1,0,len(MIDS)-1)]
KEY={'bondage':'bondage','humiliation':'humiliation','nonconsent':'nonconsent','sadomasochism':'sadomasochism',
 'sensory':'sensory','transformation':'transform','specific roles':'roles','mental alteration':'mentalalteration',
 'pregnancy':'pregnancy','genderplay':'genderplay','exhibitionism':'exhibitionself','multiple partners':'multiplepartners',
 'incest':'incest','bestiality':'bestiality','abnormal body':'abnormalbody','bodily-secretions':'secretions',
 'mythical':'mythical','creepy':'creepy','brutality':'brutality','vore':'vore','clothing':'clothing',
 'body parts':'appearance','gentleness':'gentleness','power dynamics':'powerdynamic','dirtiness':'dirty',
 'eagerness':'eagerness','objects':'objects','toys':'toys'}
EXTRA={'age-related':'agegap','older':'older','progression':'progression','regression':'regression',
       'multiple partners':'multiplepartners','preferred':'appearance'}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
def pairs_for(mode):
    K=dict(KEY) if mode!='matched27' else dict(KEY)
    if mode=='matched_plus': K.update(EXTRA)
    out=[]
    for c in ons:
        lc=c.lower()
        for k,v in K.items():
            if k in lc and v in df.columns: out.append((c,v)); break
    return list(dict.fromkeys(out))
SETS={'matched27 (as published)':pairs_for('matched27'),'matched_plus':pairs_for('matched_plus')}
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; age=df['age'].map(AGEMAP)
def build(pairs):
    O=pd.DataFrame({v:df[c].map(BIN) for c,v in pairs})
    P=pd.DataFrame({v:pd.to_numeric(df[v],errors='coerce') for _,v in pairs})
    def dc(M):
        Z=M.copy(); Z=Z.sub(Z.mean(axis=1),axis=0); return Z.sub(Z.mean(axis=0),axis=1)
    Or,Pr=dc(O),dc(P)
    for c in Or.columns:
        m=Or[c].notna()&age.notna()
        X=np.c_[np.ones(m.sum()),age[m].values]; b,*_=lstsq(X,Or.loc[m,c].values,rcond=None)
        Or.loc[m,c]=Or.loc[m,c].values-X@b
    return O,Or,Pr
def top_eig(Ovals,O,Pr,mask,k,CP,iu):
    Om=pd.DataFrame(np.where(mask,Ovals,np.nan),columns=O.columns)
    Z=Om.sub(Om.mean(axis=1),axis=0); Z=Z.sub(Z.mean(axis=0),axis=1)
    CO=np.full((k,k),np.nan)
    for i in range(k):
        for j in range(i+1,k):
            m=Z.iloc[:,i].notna()&Z.iloc[:,j].notna()&Pr.iloc[:,i].notna()&Pr.iloc[:,j].notna()
            if m.sum()<150: continue
            CO[i,j]=CO[j,i]=np.corrcoef(Z.iloc[:,i][m],Z.iloc[:,j][m])[0,1]
    ok=~np.isnan(CO[iu])&~np.isnan(CP[iu])
    co,cp=CO[iu][ok],CP[iu][ok]
    b,*_=lstsq(np.c_[np.ones(len(cp)),cp],co,rcond=None)
    r=co-np.c_[np.ones(len(cp)),cp]@b
    M=np.zeros((k,k)); M[iu[0][ok],iu[1][ok]]=r; M=M+M.T
    return float(np.abs(eigvalsh(M))[::-1][0])
rows=[]
for name,pairs in SETS.items():
    O,Or,Pr=build(pairs); k=len(pairs); mask=O.notna().values; iu=np.triu_indices(k,1)
    CP=np.full((k,k),np.nan)
    for i in range(k):
        for j in range(i+1,k):
            m=Pr.iloc[:,i].notna()&Pr.iloc[:,j].notna()&O.iloc[:,i].notna()&O.iloc[:,j].notna()
            if m.sum()<150: continue
            CP[i,j]=CP[j,i]=np.corrcoef(Pr.iloc[:,i][m],Pr.iloc[:,j][m])[0,1]
    obs=top_eig(O.values,O,Pr,mask,k,CP,iu)
    S=np.nan_to_num(CP,nan=0.0); np.fill_diagonal(S,1.0)
    w,V=np.linalg.eigh(S); w=np.clip(w,1e-6,None); L=V@np.diag(np.sqrt(w))
    mu=np.nanmean(O.values); sd=np.nanstd(O.values)
    def draw(rng,noise,inject=0.0):
        Zs=rng.normal(size=(len(O),k))@L.T
        if inject>0:
            U=rng.normal(size=(k,2)); Zs=Zs+inject*(rng.normal(size=(len(O),2))@U.T)
        Zs=Zs/Zs.std(0,keepdims=True)
        return binify(mu+sd*(np.sqrt(1-noise**2)*Zs+noise*rng.normal(size=Zs.shape)))
    for noise in [0.3,0.5,0.7]:
        vals=[]
        for seed in (2,12,22):
            rng=np.random.default_rng(seed)
            vals += [top_eig(draw(rng,noise),O,Pr,mask,k,CP,iu) for _ in range(40)]
        v=np.array(vals); hi=np.percentile(v,97.5)
        rows.append(dict(catset=name,k=k,noise=noise,observed=round(obs,3),rival_mean=round(v.mean(),3),
                         rival_p975=round(hi,3),above=bool(obs>hi)))
    # FIX (#32): the zero-injection control must be compared to draws from the SAME stream, or
    # seed noise alone decides it. Draw the reference ONCE and test both injections against it.
    rng_ref=np.random.default_rng(101)
    base=np.array([top_eig(draw(rng_ref,0.5),O,Pr,mask,k,CP,iu) for _ in range(60)])
    hi=np.percentile(base,97.5)
    for inj in [0.0,0.3]:
        rng=np.random.default_rng(202)
        vi=np.array([top_eig(draw(rng,0.5,inject=inj),O,Pr,mask,k,CP,iu) for _ in range(60)])
        rows.append(dict(catset=name,k=k,noise=f"POSCTRL inj={inj}",observed=round(vi.mean(),3),
                         rival_mean=round(base.mean(),3),rival_p975=round(hi,3),
                         above=bool(vi.mean()>hi)))
T=pd.DataFrame(rows); T.to_csv(OUT/'full_categories.csv',index=False)
print(T.to_string(index=False))
real=T[~T.noise.astype(str).str.startswith('POSCTRL')]
pos=T[T.noise.astype(str)=='POSCTRL inj=0.3']; zero=T[T.noise.astype(str)=='POSCTRL inj=0.0']
gate=bool(pos.above.all() and (~zero.above).all())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control detected in every cell (inj=0.3): {'PASS' if pos.above.all() else 'FAIL'}")
print(f"  zero-injection NOT detected in every cell        : {'PASS' if (~zero.above).all() else 'FAIL'}")
if not gate: print("  -> gate FAILED : UNVERIFIED, threshold not evaluated")
elif real.above.all(): print(f"  -> observed above the rival world in ALL {len(real)} cells : CONFIRMED, and filter-free")
else: print(f"  -> observed inside the rival world in {int((~real.above).sum())} cells : the confirmation was the filter")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
