"""
E01 A03 R16 -- the rarity ordering: acquisition, report, or CENSORING?

#53 found prevalence ordering predicts held-out pairwise acquisition order at 60.75%, largely
independent of the developmental ordering, and I flagged it as an unreported finding. Before
interpreting it, the alternative that would manufacture it for free:

  CENSORING. The sample is 18-32. A category acquired late is HELD BY FEWER PEOPLE at survey time
  simply because some have not reached it yet. So "rare" and "late" are mechanically linked by the
  age window, with no acquisition or reporting mechanism involved.

DISCRIMINATOR: censoring is an age effect. It must be STRONGER IN YOUNGER RESPONDENTS, where more
acquisition is still outstanding, and WEAKER IN THE OLDEST BAND, where more of it has completed. A
real acquisition or report mechanism has no reason to vary with the respondent's age.

ESTIMAND        prevalence-ordering ranking accuracy, computed within age bands, using each band's
                OWN prevalence so the predictor is not imported across bands.
IDENTIFICATION  identified; age is observed in 5 bands.
WORLDS          A  censoring: accuracy falls monotonically with band age
                B  acquisition/report: accuracy flat across bands
KILL (CONDITIONAL) gate -- ceiling first (#50): the ONSET ordering must work in every band, else a
                   band cannot support the comparison and is dropped rather than averaged in.
                   then: oldest-band prevalence accuracy < youngest by >4 points -> CENSORING
                         difference < 2 points                                    -> NOT censoring
                         otherwise                                                -> partial
POSITIVE CTRL   the onset ordering within each band.
NEGATIVE CTRL   random orderings, 200 per band (fixing #53's own complaint that one permutation per
                seed made chance the least precise row in the table).
SEEDS           5.
MULTIPLICITY    5 bands x 3 orderings x 5 seeds, all reported.
IMPOSSIBLE      an older sample. The release is 14-32, so censoring can be tested for a gradient but
                never removed.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
O=pd.DataFrame({c:df[c].map(BIN) for c in ons}); V=O.values; mask=~np.isnan(V)
BANDS=['14-17','18-20','21-24','25-28','29-32']
age=df['age'].values
def acc(order_vals,people,rng,cap=25000):
    right=0;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(12,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b] or order_vals[a]==order_vals[b]: continue
            right+=((order_vals[a]<order_vals[b])==(V[i,a]<V[i,b])); tot+=1
            if tot>=cap: break
        if tot>=cap: break
    return 100*right/max(tot,1), tot
rows=[]
for band in BANDS:
    idx=np.flatnonzero((age==band)&(mask.sum(1)>=6))
    if len(idx)<400: print(f"  band {band}: only {len(idx)} people, dropped"); continue
    prev_band=np.array([mask[idx][:,j].mean() for j in range(V.shape[1])])   # band's OWN prevalence
    for seed in (1,2,3,4,5):
        rng=np.random.default_rng(seed)
        p=rng.permutation(idx); tr,te=p[:len(p)//2],p[len(p)//2:]
        onset_tr=np.nanmean(V[tr],axis=0)
        a_on,_=acc(onset_tr,te,rng); a_pr,n=acc(-prev_band,te,rng)
        rnd=[acc(rng.permutation(onset_tr),te,rng)[0] for _ in range(40)]
        rows.append(dict(band=band,seed=seed,n_people=len(idx),onset=a_on,prevalence=a_pr,
                         random=float(np.mean(rnd)),random_sd=float(np.std(rnd)),pairs=n))
G=pd.DataFrame(rows); G.to_csv(OUT/'censoring.csv',index=False)
S=G.groupby('band')[['onset','prevalence','random','random_sd','n_people']].median().round(2)
S=S.reindex([b for b in BANDS if b in S.index])
print("\n=== ranking accuracy by respondent age band (each band's own prevalence) ===")
print(S.to_string())
print(f"\n  random baseline now averaged over 40 permutations per seed (fixing #53's own complaint);")
print(f"  its sd within a band is {G.random_sd.median():.2f}, against #53's 7.68 seed spread")
ok=S.index.tolist()
gate=bool((S['onset']>55).all())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  onset ordering works in every band : {'PASS' if gate else 'FAIL'} (min {S['onset'].min():.2f}%)")
if not gate: print("  -> gate FAILED : UNVERIFIED")
else:
    young=float(S.loc[ok[0],'prevalence']); old=float(S.loc[ok[-1],'prevalence'])
    print(f"  prevalence accuracy: youngest band {young:.2f}%  ->  oldest {old:.2f}%   change {old-young:+.2f}")
    if young-old>4: print("  -> CENSORING : the rarity ordering is an artifact of the age window")
    elif abs(young-old)<2: print("  -> NOT CENSORING : rarity predicts acquisition order independently of the age window")
    else: print(f"  -> PARTIAL : {young-old:+.2f} points of the rarity effect is censoring")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
