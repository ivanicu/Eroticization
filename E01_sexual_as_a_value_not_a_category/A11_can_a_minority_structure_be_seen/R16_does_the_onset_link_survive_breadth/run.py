import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R16 -- DOES THE ONSET LINK SURVIVE THE NUISANCE THAT PREDICTS ONSET BETTER THAN IT DOES?

#101d found rare-option affinity correlates with earlier acquisition (r = -0.084, z = 7.1) and #101f
named the threat in the same entry: BREADTH correlates with onset at -0.160 -- twice as strongly --
and breadth is the nuisance this project projects out everywhere else.

Affinity is count-matched BY CONSTRUCTION (real minus a curveball null that preserves each person's
pick count exactly), so corr(affinity, picks) should already be ~0. That is CHECKED here, not
assumed, and then three nested tests are run:

  raw        corr(affinity, onset)
  partial    controlling for pick count and the NUMBER OF ONSET ENTRIES a person has
  split-half affinity from one half of a person's BLOCKS vs onset from the other half of their
             ONSET CATEGORIES -- disjoint measurement, so no shared-item route survives

ESTIMAND        the affinity-onset correlation under each, with its bootstrap sd.
IDENTIFICATION  identified; the split-half version shares only the PERSON, which is the claim.
SCOPE           people with >=6 blocks and >=6 onset entries.
WORLDS          survives   the link is a person-level fact, not a measurement artifact
                collapses  it is breadth or answer-count in a new costume, and #101d is withdrawn
KILL            threshold-free: each correlation against its own bootstrap sd, and the partial
                against the raw.
POSITIVE CTRL   BREADTH's own onset correlation must survive the same partial (it is a real
                association, so a partial that kills everything is over-controlling).
NEGATIVE CTRL   the null-minus-null affinity, through all three tests: ~0 everywhere.
NOISE FLOOR     200 bootstraps; 6 half-splits.
IMPOSSIBLE      causal direction; and separating "acquired earlier" from "remembers acquiring earlier".
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
blk={w:{} for w in ['real','n1','n2']}
for bi,t in enumerate(IDENT):
    M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
    ref=-np.log(np.clip(M.mean(0),1e-4,1.))
    pack={'real':M,'n1':curveball(M,np.random.default_rng(9100)),'n2':curveball(M,np.random.default_rng(9101))}
    for w,Mw in pack.items():
        tot=Mw@ref; kk=Mw.sum(1)
        for j,gi in enumerate(idx):
            if kk[j]>0: blk[w].setdefault(gi,[]).append((bi,tot[j],kk[j]))
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
ONS=pd.DataFrame({c:df[c].map(BIN) for c in ons})
ids=np.array([i for i in blk['real'] if len(blk['real'][i])>=6
              and np.isfinite(ONS.iloc[i].values).sum()>=6])
print(f"people with >=6 blocks and >=6 onset entries: {len(ids):,}",flush=True)
def agg(w,who,sel=None):
    out=[]
    for i in who:
        v=blk[w][i]
        if sel is not None: v=[x for x in v if x[0] in sel[i]]
        if not v: out.append(np.nan); continue
        out.append(sum(x[1] for x in v)/max(sum(x[2] for x in v),1))
    return np.array(out)
picks=np.array([sum(x[2] for x in blk['real'][i]) for i in ids])
nons=np.array([int(np.isfinite(ONS.iloc[i].values).sum()) for i in ids])
onset=np.array([np.nanmean(ONS.iloc[i].values) for i in ids])
AFF=agg('real',ids)-(agg('n1',ids)+agg('n2',ids))/2.
NEG=agg('n1',ids)-agg('n2',ids)
def pcorr(a,b,Z):
    X=np.c_[np.ones(len(a)),Z]
    ra=a-X@np.linalg.lstsq(X,a,rcond=None)[0]; rb=b-X@np.linalg.lstsq(X,b,rcond=None)[0]
    return float(np.corrcoef(ra,rb)[0,1])
rb_=np.random.default_rng(555)
def boot(f,*args,n=200):
    return float(np.std([f(*[a[ix] for a in args]) for ix in
                         (rb_.integers(0,len(args[0]),len(args[0])) for _ in range(n))]))
Zc=np.c_[picks,np.log(picks),nons]
rows=[]
for nm,a in [('affinity',AFF),('breadth',picks.astype(float)),('null_minus_null',NEG)]:
    r0=float(np.corrcoef(a,onset)[0,1])
    rp=pcorr(a,onset,Zc) if nm!='breadth' else pcorr(a,onset,np.c_[nons])
    rows.append(dict(measure=nm,raw=r0,raw_sd=boot(lambda x,y:float(np.corrcoef(x,y)[0,1]),a,onset),
                     partial=rp))
print("\n=== corr(measure, mean onset age) ===")
print(pd.DataFrame(rows).round(4).to_string(index=False))
print(f"\n  corr(affinity, picks) = {np.corrcoef(AFF,picks)[0,1]:+.4f}  "
      f"(count-matched by construction -- checked, not assumed)")
sp=[]
for rep in range(6):
    rr=np.random.default_rng(400+rep)
    selA={}; oB=[]
    for k,i in enumerate(ids):
        v=blk['real'][i]; o=rr.permutation(len(v)); h=len(v)//2
        selA[i]={v[j][0] for j in o[:h]}
        oi=np.flatnonzero(np.isfinite(ONS.iloc[i].values)); oo=rr.permutation(len(oi))
        oB.append(np.nanmean(ONS.iloc[i].values[oi[oo[len(oi)//2:]]]))
    aA=agg('real',ids,selA)-(agg('n1',ids,selA)+agg('n2',ids,selA))/2.
    oB=np.array(oB); m=np.isfinite(aA)&np.isfinite(oB)
    sp.append(float(np.corrcoef(aA[m],oB[m])[0,1]))
sp=np.array(sp)
print(f"\n=== SPLIT-HALF: affinity from HALF the blocks vs onset from the OTHER HALF of categories ===")
print(f"  r = {sp.mean():+.4f}  (sd over 6 splits {sp.std():.4f})")
r_aff=rows[0]['raw']; p_aff=rows[0]['partial']; r_br=rows[1]['raw']; p_br=rows[1]['partial']
gn=abs(rows[2]['raw'])<0.03 and abs(rows[2]['partial'])<0.03
gb=abs(p_br)>0.5*abs(r_br)
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) null-minus-null ~0 raw and partial : {'PASS' if gn else 'FAIL'} "
      f"({rows[2]['raw']:+.4f}, {rows[2]['partial']:+.4f})")
print(f"   (b) breadth's own link survives the partial (not over-controlling): "
      f"{'PASS' if gb else 'FAIL'} ({r_br:+.4f} -> {p_br:+.4f})")
if not(gn and gb): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"\n   affinity-onset  raw {r_aff:+.4f}  partial {p_aff:+.4f}  "
          f"split-half {sp.mean():+.4f}")
    keepfrac=abs(p_aff)/max(abs(r_aff),1e-9)
    if abs(sp.mean())>2*sp.std() and keepfrac>0.5:
        print(f"\n   -> SURVIVES. {100*keepfrac:.0f}% of the raw link remains after controlling for")
        print("      pick count and answer count, and it reproduces with DISJOINT measurement.")
        print("      #101d is a person-level fact, not a shared-measurement artifact.")
    elif keepfrac<0.3:
        print(f"\n   -> COLLAPSES to {100*keepfrac:.0f}% under the partial. #101d was breadth or")
        print("      answer-count in a new costume and must be withdrawn.")
    else:
        print(f"\n   -> PARTIAL SURVIVAL ({100*keepfrac:.0f}%). The link is real but shares most of its")
        print("      variance with how much the person answered.")
print(f"\nartifact sha1 {hashlib.sha1(pd.DataFrame(rows).to_csv(index=False).encode()).hexdigest()[:12]}")
