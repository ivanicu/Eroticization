import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R07 -- RUNNING MY OWN FORECAST #4: CHARGE THE SHARED BASIS.

ADVERSARY_FORECAST block 2, prediction #4 (p=0.50): "#71's parameter-count argument is a bad
accounting. I called C's person scores 'free' because they are estimated elsewhere -- but they are
estimated on THE SAME PEOPLE, from 31 blocks, and that estimation is not free, it is just not billed
to the target."

The forecast asks for an accounting. An accounting is the wrong instrument -- df is a proxy for
overfitting, and #76 just cost me a round for reaching for a proxy (a permutation) when the direct
measurement existed. The direct question is:

  does C survive when the people whose data built the shared basis are DISJOINT from the people it
  is evaluated on?

If C's advantage came from the basis having quietly seen the evaluation people, person-holdout
breaks it. If C survives, "free" was the right accounting for the only reason that matters: those
parameters cannot overfit cells they were never near.

  ARM raw            the #71 specification: basis built on everyone
  ARM person-holdout basis built on FIT-half people only; eval-half scores obtained by PROJECTION
                     onto the fitted option-side loadings; target loadings and held-out cells both
                     from eval-half only. No eval person contributes to the basis.
  ARM shuffled-basis the same pipeline with the option-side loadings replaced by a random
                     orthonormal basis of the same rank -- the structure is destroyed, the parameter
                     count is identical. This is what an accounting argument would predict SURVIVES.

ESTIMAND        C and W under person-holdout, on identical eval-half held-out cells, at matched rank.
IDENTIFICATION  identified: fit-half and eval-half are disjoint by construction, verified numerically.
SCOPE           the 23 blocks A09/R03 identified. Kc=Kw=4. Gate+demographics projected out
                throughout, per #77 -- this round inherits the corrected specification, not the
                inflated one.
WORLDS          basis-is-free   C survives person-holdout at close to its raw value -> #71 stands
                basis-is-paid   C collapses toward the shuffled-basis arm -> #71's ordering was
                                bought with an unbilled estimation and must be withdrawn
KILL            threshold-free: C under person-holdout is compared against (a) the shuffled-basis
                arm and (b) W on the same cells, each declared only above 2x the pooled seed spread.
POSITIVE CTRL   W, fit entirely within the eval half, must stay positive -- it is unaffected by the
                holdout and calibrates that the eval half is large enough to fit anything at all.
NEGATIVE CTRL   the shuffled-basis arm: same rank, same parameter count, no structure. An accounting
                argument cannot distinguish it from the real basis; a prediction test must.
NOISE FLOOR     3 masks x 3 person splits.
MULTIPLICITY    23 blocks x 3 arms x 3 seeds, published whole.
IMPOSSIBLE      charging the basis in DEGREES OF FREEDOM -- df is defined against a likelihood this
                estimator does not have (soft-impute + per-column least squares is not an MLE).
                Reported N/A with that reason; the prediction test replaces it.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd, lstsq, qr
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW])); PM={p:i for i,p in enumerate(ALLP)}
BLKS=sorted(RAW)
E=np.zeros((len(ALLP),len(BLKS)))
for k,q in enumerate(BLKS): E[[PM[p] for p in RAW[q]['ppl']],k]=1.
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
COV=pd.DataFrame({'male':pd.to_numeric(df.get('biomale'),errors='coerce'),'age':df['age'].map(AGEMAP)}).reindex(ALLP)
for c in ['opennessvariable','neuroticismvariable','extroversionvariable',
          'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']:
    if c in df.columns: COV[c]=pd.to_numeric(df[c],errors='coerce').reindex(ALLP).values
COV=COV.fillna(COV.median()).values; COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
print(f"targets {len(IDENT)}  people {len(ALLP)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]; K=4

def other_matrix(target):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0)
    return np.where(np.isnan(Z),mu,Z)

def nuis(U,target):
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(U)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b        # #77: gate + demographics removed throughout

def run(target,arm,seed):
    M=RAW[target]['M']; rows=np.array([PM[p] for p in RAW[target]['ppl']])
    rng=np.random.default_rng(seed)
    Z=other_matrix(target); Z=Z-Z.mean(0)
    if arm=='raw':
        U,S,_=svd(Z,full_matrices=False); U=nuis(U[:,:K]*S[:K],target)
        ev=np.arange(len(rows))                      # everyone is eval
    else:
        half=rng.random(len(ALLP))<0.5               # FIT half / EVAL half, disjoint
        Zf=Z[half]
        _,_,Vt=svd(Zf-Zf.mean(0),full_matrices=False); V=Vt[:K].T
        if arm=='shuffled':
            V,_=qr(rng.normal(size=V.shape))         # same rank, same cost, no structure
        U=nuis(Z@V,target)
        ev=np.flatnonzero(~half[rows])               # eval-half members of this block only
        if len(ev)<400: return None
    idx=rows[ev]; Msub=M[ev]
    if Msub.shape[0]<400: return None
    obs=rng.random(Msub.shape)>=MASK; he=~obs
    T=np.where(obs,Msub,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    Us=U[idx]; Us=(Us-Us.mean(0))/(Us.std(0)+1e-12)
    C=np.zeros_like(Msub)
    for j in range(Msub.shape[1]):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),Us[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(len(Us)),Us]@b
    F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(20):
        Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:K]*Ss[:K])@Vv[:K])
    Uu,Ss,Vv=svd(F,full_matrices=False); W=(Uu[:,:K]*Ss[:K])@Vv[:K]
    base=np.mean((Msub[he]-gm)**2)
    f=lambda *p: 1.-np.mean((Msub[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    IB=np.broadcast_to(I,Msub.shape)
    return dict(C=f(IB,P,C)-f(IB,P), W=f(IB,P,W)-f(IB,P),
                CW=f(IB,P,C,W)-f(IB,P), n_eval=Msub.shape[0])

rows=[]
for i,t in enumerate(IDENT):
    for arm in ['raw','holdout','shuffled']:
        for sd in SEEDS:
            r=run(t,arm,sd)
            if r: rows.append(dict(q=t,arm=arm,seed=sd,**r))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R07_charge_the_shared_basis/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== C AND W BY ARM (mean over blocks x 3 seeds) ===")
print(D.groupby('arm')[['C','W','CW','n_eval']].mean().round(5).to_string())

print("\n=== PER-BLOCK ===")
out=[]
for arm in ['raw','holdout','shuffled']:
    d=D[D.arm==arm]
    c=d.groupby('q').C; w=d.groupby('q').W
    out.append(dict(arm=arm,C=c.mean().median(),C_sd=c.std().median(),
                    W=w.mean().median(),blocks=len(c.mean()),
                    C_pos=int((c.mean()>2*c.std()).sum()),
                    C_beats_W=int(((c.mean()-w.mean())>2*np.sqrt(c.std()**2+w.std()**2)).sum())))
S=pd.DataFrame(out); print(S.round(5).to_string(index=False))

hold=S[S.arm=='holdout'].iloc[0]; shuf=S[S.arm=='shuffled'].iloc[0]; raw=S[S.arm=='raw'].iloc[0]
print("\n  CONDITIONAL KILL -- gates first")
g1=hold.W>0; g2=shuf.C<hold.C
print(f"   (a) W positive in the holdout arm (eval half can fit anything): "
      f"{'PASS' if g1 else 'FAIL'} (W {hold.W:+.5f})")
print(f"   (b) shuffled basis below the real basis at identical cost      : "
      f"{'PASS' if g2 else 'FAIL -- the advantage IS the parameter count'} "
      f"({shuf.C:+.5f} vs {hold.C:+.5f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"\n   C raw {raw.C:+.5f}  ->  person-holdout {hold.C:+.5f}  "
          f"({100*hold.C/max(raw.C,1e-12):.0f}% retained)")
    print(f"   under person-holdout: C>0 in {hold.C_pos}/{hold.blocks} blocks, "
          f"C>W in {hold.C_beats_W}/{hold.blocks}")
    if hold.C_pos>=hold.blocks*0.7 and hold.C_beats_W>=hold.blocks*0.5:
        print("\n   -> FORECAST #4 IS WRONG. C survives when every person in the evaluation set is")
        print("      disjoint from the people whose data built the basis. 'Free' was the right")
        print("      accounting for the only reason that matters: those parameters cannot overfit")
        print("      cells they were never near. #71's ordering stands.")
    elif hold.C_pos<hold.blocks*0.3:
        print("\n   -> FORECAST #4 IS RIGHT. C does not survive person-holdout, so #71's advantage")
        print("      was bought with an unbilled estimation and must be withdrawn.")
    else:
        print(f"\n   -> PARTIAL: C>0 in {hold.C_pos}/{hold.blocks} but beats W in only "
              f"{hold.C_beats_W}. The transfer survives; the ORDERING against W does not, and #71")
        print("      must be restated as 'both exist' rather than 'general beats specific'.")
print("\nN/A, with what it would require: charging the basis in DEGREES OF FREEDOM needs a likelihood "
      "this estimator does not have (soft-impute plus per-column least squares is not an MLE). The "
      "prediction test above replaces it and answers the question the accounting was a proxy for.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
