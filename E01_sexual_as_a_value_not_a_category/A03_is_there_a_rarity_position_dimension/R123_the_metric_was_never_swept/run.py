import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R08 -- THE ONE AXIS THE SPECIFICATION CURVE NEVER SWEPT: THE LOSS.

A09 and A10 swept rank, estimand, null, block, ordering, projection and score type. Every single
number in both arcs is SQUARED ERROR on a binary cell. realstat G4 lists metric-robustness among the
axes a specification curve must cover, and this project has none.

ADVERSARY_FORECAST block 2, prediction #3 (p=0.55): "the item/interaction comparison is
scale-dependent... under log-loss the interaction's share rises, because squared error on a 0/1
outcome under-weights confident-and-wrong predictions, which is exactly where a person-specific
readout would differ from a base rate."

That is a sharp, mechanistic prediction and it is testable two ways, because the metric enters twice:

  ARM eval   the SAME least-squares components, evaluated under Brier / log-loss / L1.
             Tests whether the ORDERING is a property of the evaluation scale.
  ARM fit    components FITTED on the logit scale -- item as column log-odds, person as a log-odds
             offset, cross-block as per-column logistic regression -- evaluated under log-loss.
             Tests whether the ordering is a property of the ESTIMATOR, which is a different claim.

ESTIMAND        Shapley skill of {I, P, C, W} under three losses, and of {I, P, C} under a logistic
                estimator, on identical held-out cells; and the item:interaction ratio under each.
IDENTIFICATION  identified. Skill is 1 - loss(model)/loss(grand-mean baseline) under each loss, so
                the three are on comparable [0,1]-ish scales without being the same number.
SCOPE           the 23 blocks A09/R114 identified. Kc=Kw=4. Gate + demographics projected out (#77).
WORLDS          scale-free   the ordering item > interaction holds under all three losses and both
                             estimators -> A09/A10's headline is not an artefact of squared error
                metric-bound the ordering flips under log-loss -> #67/#68/#70 are statements about
                             Brier score and must be reworded as such
KILL            threshold-free: the ordering is declared per (metric, estimator) cell only above 2x
                that cell's own seed spread, and all six cells are published including disagreements.
POSITIVE CTRL   a synthetic world with a KNOWN rank-5 shared person factor: every metric must
                recover a large interaction there. A metric that cannot see a planted interaction
                cannot be used to say the real one is small.
NEGATIVE CTRL   person-permutation, run under every metric: must be ~0 in all of them.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 3 losses x 2 estimators x 3 seeds, published whole.
IMPOSSIBLE      a loss that is neither proper nor rank-based (e.g. raw accuracy) -- it is not a
                scoring rule and a decomposition under it is not interpretable. Excluded by design,
                not omitted.
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
from numpy.linalg import svd, lstsq
from math import factorial
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
               'R114_fixed_margin_null/results/grid.csv')
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
print(f"targets {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]; K=4; EPS=1e-3
def logit(p): p=np.clip(p,EPS,1-EPS); return np.log(p/(1-p))
def sig(z): return 1/(1+np.exp(-np.clip(z,-30,30)))

LOSSES={'brier': lambda y,p:(y-p)**2,
        'logloss': lambda y,p:-(y*np.log(np.clip(p,EPS,1))+(1-y)*np.log(np.clip(1-p,EPS,1))),
        'l1':    lambda y,p:np.abs(y-p)}

def other_scores(target,blocks=None):
    B=blocks if blocks is not None else RAW
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=B[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in B[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:K]*S[:K]
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b

def shapley(v,names):
    n=len(names); out={}
    for c in names:
        o=[x for x in names if x!=c]; tot=0.
        for r in range(n):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(n-len(S)-1)/factorial(n)*(
                    v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def decompose(M,U_all,rows,seed,estimator,permute=False):
    n,m=M.shape
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm)
    U=U_all[rows]
    if permute: U=U[np.random.default_rng(seed+555).permutation(n)]
    U=(U-U.mean(0))/(U.std(0)+1e-12)
    if estimator=='ls':
        I=(cm-gm)[None,:]; T1=T-gm-I
        rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
        Rres=T1-P
        C=np.zeros_like(M)
        for j in range(m):
            k=obs[:,j]
            if k.sum()<50: continue
            b,*_=lstsq(np.c_[np.ones(k.sum()),U[k]],Rres[k,j],rcond=None)
            C[:,j]=np.c_[np.ones(n),U]@b
        F=np.where(np.isnan(Rres),0.,Rres)
        for _ in range(20):
            Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:K]*Ss[:K])@Vv[:K])
        Uu,Ss,Vv=svd(F,full_matrices=False); W=(Uu[:,:K]*Ss[:K])@Vv[:K]
        comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':W}
        link=lambda z: np.clip(gm+z,EPS,1-EPS); names='IPCW'
    else:                                                   # logistic: everything on the logit scale
        g0=logit(gm); I=(logit(cm)-g0)[None,:]
        base=np.clip(sig(g0+I),EPS,1-EPS)
        num=np.nansum(np.where(obs,T-base,np.nan),axis=1)
        den=np.nansum(np.where(obs,base*(1-base),np.nan),axis=1)+1e-9
        P=(num/den)[:,None]                                 # one Newton step for the person offset
        eta=g0+I+P; mu_=np.clip(sig(eta),EPS,1-EPS)
        wgt=mu_*(1-mu_); z=np.where(obs,(T-mu_)/np.maximum(wgt,1e-6),0.)
        C=np.zeros_like(M)
        for j in range(m):                                   # one IRLS step per column
            k=obs[:,j]
            if k.sum()<50: continue
            X=np.c_[np.ones(k.sum()),U[k]]; sw=np.sqrt(wgt[k,j])
            b,*_=lstsq(X*sw[:,None],z[k,j]*sw,rcond=None)
            C[:,j]=np.c_[np.ones(n),U]@b
        comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C}
        link=lambda zz: np.clip(sig(g0+zz),EPS,1-EPS); names='IPC'
    out={}
    for lname,L in LOSSES.items():
        b0=np.mean(L(M[he],np.full(he.sum(),gm)))
        v={}
        for bits in range(1<<len(names)):
            S=frozenset([c for j,c in enumerate(names) if bits>>j&1])
            p=link(sum(comp[c] for c in S)) if S else np.full(M.shape,gm)
            v[S]=1.-np.mean(L(M[he],np.asarray(p)[he]))/b0
        sh=shapley(v,names)
        for c in names: out[f"{lname}_{c}"]=sh[c]
        out[f"{lname}_full"]=v[frozenset(names)]
    return out

def synth_rank5(seed):
    rng=np.random.default_rng(9000+seed); F=rng.normal(size=(len(ALLP),5)); out={}
    for q in BLKS:
        M=RAW[q]['M']; n,m=M.shape; rows=[PM[p] for p in RAW[q]['ppl']]
        L=rng.normal(size=(5,m))*0.30
        p=np.clip(M.mean(0)[None,:]+F[rows]@L,0.02,0.98)
        out[q]=dict(M=(rng.random((n,m))<p).astype(float),ppl=RAW[q]['ppl'])
    return out

rows=[]
for i,t in enumerate(IDENT):
    U=other_scores(t); r_=np.array([PM[p] for p in RAW[t]['ppl']])
    for est in ['ls','logit']:
        for sd in SEEDS:
            rows.append(dict(q=t,world='real',est=est,arm='real',seed=sd,
                             **decompose(RAW[t]['M'],U,r_,sd,est)))
            rows.append(dict(q=t,world='real',est=est,arm='perm',seed=sd,
                             **decompose(RAW[t]['M'],U,r_,sd,est,permute=True)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
CT=IDENT[:8]
for i,t in enumerate(CT):
    S5=synth_rank5(1); U=other_scores(t,S5); r_=np.array([PM[p] for p in S5[t]['ppl']])
    for sd in SEEDS:
        rows.append(dict(q=t,world='rank5',est='ls',arm='real',seed=sd,
                         **decompose(S5[t]['M'],U,r_,sd,'ls')))
    print(f"  [ctrl] {i+1}/{len(CT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R123_the_metric_was_never_swept/results/'
D.to_csv(OUT+'grid.csv',index=False)

R=D[(D.world=='real')&(D.arm=='real')]; N=D[(D.world=='real')&(D.arm=='perm')]
print("\n=== POSITIVE CONTROL: a planted rank-5 interaction, seen by each metric ===")
c5=D[D.world=='rank5']
print(c5.groupby('est')[[f"{l}_{c}" for l in LOSSES for c in 'IPCW']].mean().round(4).T.to_string())

print("\n=== NEGATIVE CONTROL: person-permutation, C under each metric ===")
print(N.groupby('est')[[f"{l}_C" for l in LOSSES]].mean().round(5).to_string())

print("\n=== THE DECOMPOSITION UNDER EACH LOSS (mean over 23 blocks x 3 seeds) ===")
for est in ['ls','logit']:
    d=R[R.est==est]; names='IPCW' if est=='ls' else 'IPC'
    t=pd.DataFrame({l:[d[f"{l}_{c}"].mean() for c in names]+[d[f"{l}_full"].mean()]
                    for l in LOSSES},index=list(names)+['full'])
    print(f"\n  estimator = {est}")
    print(t.round(4).to_string())
    inter=(t.loc['C']+t.loc['W']) if est=='ls' else t.loc['C']
    print("  item : interaction  " + "   ".join(
        f"{l} {t.loc['I',l]/max(inter[l],1e-9):6.1f}x" for l in LOSSES))

print("\n=== THE ORDERING, per (estimator, loss), declared above 2x its own seed spread ===")
out=[]
for est in ['ls','logit']:
    d=R[R.est==est]
    for l in LOSSES:
        I=d.groupby('q')[f"{l}_I"].mean()
        X=d.groupby('q')[f"{l}_C"].mean()+(d.groupby('q')[f"{l}_W"].mean() if est=='ls' else 0)
        sp=np.sqrt(d.groupby('q')[f"{l}_I"].std()**2+d.groupby('q')[f"{l}_C"].std()**2)
        gap=I-X
        out.append(dict(est=est,loss=l,I=I.median(),X=X.median(),gap=gap.median(),
                        spread2=2*sp.median(),item_wins=int((gap>2*sp).sum()),
                        inter_wins=int((-gap>2*sp).sum()),tied=int((gap.abs()<=2*sp).sum())))
S=pd.DataFrame(out); print(S.round(4).to_string(index=False))

print("\n  CONDITIONAL KILL -- gates first")
g_pos=all(c5[f"{l}_C"].mean()+c5[f"{l}_W"].mean()>0.01 for l in LOSSES)
g_neg=all(abs(N[N.est=='ls'][f"{l}_C"].mean())<0.005 for l in LOSSES)
print(f"   (a) every metric sees a planted rank-5 interaction : {'PASS' if g_pos else 'FAIL'}")
print(f"   (b) person-permutation ~0 under every metric        : {'PASS' if g_neg else 'FAIL'}")
if not(g_pos and g_neg): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    flips=S[S.inter_wins>S.item_wins]
    print(f"\n   cells where the ITEM wins: {S.item_wins.tolist()}  (of 23 blocks each)")
    print(f"   cells where the INTERACTION wins: {S.inter_wins.tolist()}")
    if len(flips)==0:
        print("\n   -> SCALE-FREE. The ordering holds under Brier, log-loss and L1, and under both a")
        print("      least-squares and a logistic estimator. FORECAST #3 IS WRONG: A09/A10's headline")
        print("      is not an artefact of squared error.")
    else:
        print(f"\n   -> METRIC-BOUND in {len(flips)} of {len(S)} cells: {flips[['est','loss']].values.tolist()}")
        print("      #67/#68/#70 are statements about the loss they were computed under and must be")
        print("      reworded. FORECAST #3 IS RIGHT.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
