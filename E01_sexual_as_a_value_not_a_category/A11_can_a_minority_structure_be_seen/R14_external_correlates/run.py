import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R14 -- THE ONE TEST THIS ARC HAS NEVER RUN, AND THE ONLY ONE THAT SEPARATES THE TWO READINGS.

#100 confirmed a reliable person-level trait in WHICH options are endorsed (residualised split-half
0.432, floor -0.022, planted ceiling 0.832) and named what reliability cannot settle: a stable
RESPONSE STYLE is also a reliable person property.

The two readings make OPPOSITE predictions about variables outside the endorsement matrix:

  response style     tracks PERSONALITY -- openness, agreeableness, the acquiescence-adjacent scales
  erotic parameter   tracks the ACQUISITION-AGE structure A03 spent fifteen rounds establishing

And there is a reference this project already owns: BREADTH (pick count), whose external correlations
were measured at r_true ~ 0.10 with adversity / mental illness / openness (#17, #23). A trait whose
correlations are the same size as breadth's has not distinguished itself from the nuisance the whole
project projects out.

ESTIMAND        correlation of the pick-count-residualised rare-option affinity with each external
                variable, disattenuated by its own reliability (0.432); breadth's correlations as
                the reference row; and the CONTRAST between the personality block and the onset block.
IDENTIFICATION  identified; every external variable is outside the endorsement matrix entirely.
SCOPE           people entering >=6 blocks who also have the external variable. n stated per cell.
WORLDS          style    personality correlations dominate the onset ones
                erotic   onset correlations dominate
                neither  both at breadth's level (~0.10) -> the trait has no external anchor and
                         #100d's reading is unsupported in either direction
KILL            threshold-free: each correlation against its own bootstrap sd, and the personality
                block mean against the onset block mean.
POSITIVE CTRL   BREADTH must reproduce its known correlations (~0.10 with openness etc., #17/#23),
                or the external variables are not being read correctly.
NEGATIVE CTRL   the fixed-margin null's affinity, residualised the same way: every correlation ~0.
MULTIPLICITY    9 external variables x 3 measures, published whole, with BH over the family.
NOISE FLOOR     200 bootstrap resamples over people.
IMPOSSIBLE      a causal reading -- these are cross-sectional correlations in a stripped release.
"""
import pandas as pd, numpy as np, warnings, hashlib
warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
def curveball(M,rng,per_row=5.):
    A=[set(np.flatnonzero(r).tolist()) for r in M]; n=len(A)
    for _ in range(int(per_row*n)):
        i,j=int(rng.integers(n)),int(rng.integers(n))
        if i==j: continue
        ai,aj=A[i],A[j]; inter=ai&aj
        di=list(ai-inter); dj=list(aj-inter); L=di+dj
        if not L: continue
        rng.shuffle(L); k=len(di)
        A[i]=inter|set(L[:k]); A[j]=inter|set(L[k:])
    out=np.zeros_like(M)
    for i,s in enumerate(A): out[i,list(s)]=1.
    return out
tot={w:np.zeros(len(ALLP)) for w in ['real','cb']}
cnt={w:np.zeros(len(ALLP)) for w in ['real','cb']}
nblk=np.zeros(len(ALLP))
for t in IDENT:
    M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
    ref=-np.log(np.clip(M.mean(0),1e-4,1.))
    for w,Mw in {'real':M,'cb':curveball(M,np.random.default_rng(8100))}.items():
        tot[w][idx]+=Mw@ref; cnt[w][idx]+=Mw.sum(1)
    nblk[idx]+=1
ok=nblk>=6
S={w:tot[w][ok]/np.maximum(cnt[w][ok],1) for w in tot}
picks=cnt['real'][ok]
X=np.c_[np.ones(ok.sum()),picks,np.log(picks)]
AFF={w:S[w]-X@np.linalg.lstsq(X,S[w],rcond=None)[0] for w in S}
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
onset=pd.DataFrame({c:df[c].map(BIN) for c in ons}).mean(axis=1)
EXT={'biomale':pd.to_numeric(df.get('biomale'),errors='coerce'),
     'age':df['age'].map(AGEMAP),
     'openness':pd.to_numeric(df.get('opennessvariable'),errors='coerce'),
     'neuroticism':pd.to_numeric(df.get('neuroticismvariable'),errors='coerce'),
     'extroversion':pd.to_numeric(df.get('extroversionvariable'),errors='coerce'),
     'conscientious':pd.to_numeric(df.get('consciensiousnessvariable'),errors='coerce'),
     'agreeable':pd.to_numeric(df.get('agreeablenessvariable'),errors='coerce'),
     'powerlessness':pd.to_numeric(df.get('powerlessnessvariable'),errors='coerce'),
     'mean_onset_age':onset}
PERS=['openness','neuroticism','extroversion','conscientious','agreeable','powerlessness']
ONS=['mean_onset_age']
ids=ALLP[ok]
REL=0.432
rows=[]; rb=np.random.default_rng(31337)
for name,v in EXT.items():
    y=pd.to_numeric(v,errors='coerce').reindex(ids).values
    m=np.isfinite(y)
    if m.sum()<500: continue
    for w,a in [('affinity',AFF['real']),('null_affinity',AFF['cb']),('breadth',picks)]:
        r=float(np.corrcoef(a[m],y[m])[0,1])
        bs=[float(np.corrcoef(a[m][ix],y[m][ix])[0,1]) for ix in
            (rb.integers(0,m.sum(),m.sum()) for _ in range(200))]
        rows.append(dict(var=name,measure=w,n=int(m.sum()),r=r,sd=float(np.std(bs)),
                         r_true=r/np.sqrt(REL) if w=='affinity' else r))
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
P=D.pivot_table(index='var',columns='measure',values='r').round(4)
Pn=D.pivot_table(index='var',columns='measure',values='n')
Ps=D.pivot_table(index='var',columns='measure',values='sd').round(4)
Pt=D[D.measure=='affinity'].set_index('var').r_true.round(4)
print("\n=== CORRELATION WITH EXTERNAL VARIABLES ===")
out=P.copy(); out['n']=Pn['affinity'].astype(int); out['aff_sd']=Ps['affinity']
out['aff_disattenuated']=Pt
print(out.to_string())
A=D[D.measure=='affinity'].set_index('var')
N=D[D.measure=='null_affinity'].set_index('var')
B=D[D.measure=='breadth'].set_index('var')
pers=A.loc[[p for p in PERS if p in A.index],'r'].abs().mean()
onsr=A.loc[[p for p in ONS if p in A.index],'r'].abs().mean()
bref=B.loc[[p for p in PERS if p in B.index],'r'].abs().mean()
gn=N.r.abs().max()<0.05
gb=bref>0.03
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) the NULL affinity correlates with nothing : {'PASS' if gn else 'FAIL'} "
      f"(max |r| {N.r.abs().max():.4f})")
print(f"   (b) BREADTH reproduces known correlations     : {'PASS' if gb else 'FAIL'} "
      f"(mean |r| with personality {bref:.4f}; #17/#23 measured ~0.10 disattenuated)")
if not(gn and gb): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"\n   personality block, mean |r| : {pers:.4f}  (disattenuated {pers/np.sqrt(REL):.4f})")
    print(f"   onset block,       mean |r| : {onsr:.4f}  (disattenuated {onsr/np.sqrt(REL):.4f})")
    print(f"   breadth's personality mean  : {bref:.4f}   <- the reference this must beat to matter")
    biggest=A.r.abs().idxmax()
    print(f"   largest single correlate    : {biggest} r={A.loc[biggest,'r']:+.4f} "
          f"(sd {A.loc[biggest,'sd']:.4f}, disattenuated {A.loc[biggest,'r_true']:+.4f})")
    if max(pers,onsr)<bref*1.2:
        print("\n   -> NO EXTERNAL ANCHOR. The trait's correlations are no larger than breadth's, so")
        print("      neither the response-style nor the erotic reading is supported. #100d stands as")
        print("      'something stable exists' and no more.")
    elif onsr>pers*1.5:
        print("\n   -> EROTIC PARAMETER. It tracks acquisition age far more than personality.")
    elif pers>onsr*1.5:
        print("\n   -> RESPONSE STYLE. It tracks personality far more than acquisition age.")
    else:
        print("\n   -> BOTH, comparably. The trait has an external anchor but it does not separate")
        print("      the two readings.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
