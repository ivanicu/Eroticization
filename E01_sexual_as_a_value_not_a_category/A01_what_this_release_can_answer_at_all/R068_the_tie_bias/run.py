"""
E01 A03 R19 -- 66.7% is accuracy on the pairs this release can order. What is it on all of them?

#52 reported the schedule at 66.71% held-out pairwise accuracy and flagged, honestly but without
quantifying it, that 36.3% of within-person pairs are TIED by the 2-year binning and were EXCLUDED.
Excluded pairs are not a random third: a tie means the two onsets fell in the same bin, so ties are
the pairs with the SMALLEST true gaps -- the hard ones. Scoring only the resolvable pairs inflates
accuracy by exactly the amount the hard pairs would have dragged it down.

  quantify: accuracy as a function of the observed onset gap. If accuracy rises steeply with gap,
  66.7% is "accuracy on easy pairs" and the all-pairs figure -- ties scored at 0.5, which is what a
  predictor that cannot break them earns -- is lower.

ESTIMAND        held-out pairwise accuracy by |onset gap|, and the all-pairs accuracy with ties at
                half credit.
IDENTIFICATION  identified; the gap is observed and the tie rule is arithmetic.
WORLDS          A  gap-independent: accuracy flat across gaps, and 66.7% generalises
                B  gap-driven: accuracy rises with gap, and the headline is an easy-pairs number
KILL (CONDITIONAL) gate: accuracy must exceed chance in EVERY non-tied gap bin -- if the ordering
                   fails on small gaps entirely, the bins are not comparable.
                   then: accuracy at the smallest gap bin within 5 points of the largest -> FLAT
                         spread across bins > 15 points -> GAP-DRIVEN, and the headline is restated
POSITIVE CTRL   the largest-gap bin, where the ordering should do best.
NEGATIVE CTRL   a random ordering, per gap bin -- must sit at chance in all of them.
SEEDS           4.
MULTIPLICITY    5 gap bins x 2 orderings x 4 seeds, all reported, plus the all-pairs figure.
IMPOSSIBLE      un-binning the onsets. The tie rate is a property of the release, so the all-pairs
                figure is the honest ceiling for THIS data, not an estimate of the true one.
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
keep=np.flatnonzero(mask.sum(1)>=6)
GAPS=[(0,0),(0.1,2.5),(2.5,4.5),(4.5,8.5),(8.5,99)]
rows=[]
for seed in (1,2,3,4):
    rng=np.random.default_rng(seed)
    p=rng.permutation(keep); tr,te=p[:len(p)//2],p[len(p)//2:]
    order=np.nanmean(V[tr],axis=0); rnd=rng.permutation(order)
    tally={g:[0,0] for g in GAPS}; rtally={g:[0,0] for g in GAPS}
    for i in te:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(12,len(j))):
            a,b=rng.choice(j,2,replace=False)
            gap=abs(V[i,a]-V[i,b])
            g=next(gg for gg in GAPS if gg[0]<=gap<=gg[1] or (gg==(0,0) and gap==0))
            if gap==0:
                tally[(0,0)][1]+=1; rtally[(0,0)][1]+=1; continue
            for ov,t in ((order,tally),(rnd,rtally)):
                if ov[a]==ov[b]: continue
                t[g][0]+=((ov[a]<ov[b])==(V[i,a]<V[i,b])); t[g][1]+=1
    for g in GAPS:
        r,n=tally[g]; rr,rn=rtally[g]
        rows.append(dict(seed=seed,gap_lo=g[0],gap_hi=g[1],n=n,
                         acc=(100*r/n if n and g!=(0,0) else np.nan),
                         rnd=(100*rr/rn if rn and g!=(0,0) else np.nan)))
G=pd.DataFrame(rows); G.to_csv(OUT/'tie_bias.csv',index=False)
S=G.groupby(['gap_lo','gap_hi']).agg(acc=('acc','median'),rnd=('rnd','median'),n=('n','median')).reset_index()
S['share']=100*S['n']/S['n'].sum()
print("\n=== accuracy by true onset gap (years) ===")
print(S.round(2).to_string(index=False))
nz=S[S.gap_lo>0]
tie_share=float(S[S.gap_lo==0]['share'].iloc[0])
resolved=float((nz.acc*nz.n).sum()/nz.n.sum())
allpairs=(resolved*nz.n.sum()+50.0*S[S.gap_lo==0]['n'].iloc[0])/S['n'].sum()
print(f"\n  tied pairs: {tie_share:.1f}% of all pairs")
print(f"  accuracy on resolvable pairs (the published figure): {resolved:.2f}%")
print(f"  ALL-PAIRS accuracy with ties at half credit         : {allpairs:.2f}%")
ga=bool((nz.acc>52).all()); gb=bool(nz.rnd.between(45,55).all())
spread=float(nz.acc.max()-nz.acc.min())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) accuracy above chance in every non-tied bin : {'PASS' if ga else 'FAIL'} (min {nz.acc.min():.1f}%)")
print(f"  (b) random at chance in every bin               : {'PASS' if gb else 'FAIL'} ({nz.rnd.min():.1f}-{nz.rnd.max():.1f})")
if not (ga and gb): print("  -> gate FAILED : UNVERIFIED")
elif spread<5: print(f"  -> FLAT across gaps ({spread:.1f} points): 66.7% generalises")
elif spread>15: print(f"  -> GAP-DRIVEN ({spread:.1f} points across bins): the headline is an easy-pairs number")
else: print(f"  -> partial: {spread:.1f} points across gap bins")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
