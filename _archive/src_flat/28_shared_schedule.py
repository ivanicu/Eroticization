"""
Pairwise co-variation says categories cluster in time. It does not say there is a shared ORDER.
Direct test: within each person, rank their own onsets and correlate with the POPULATION mean
onset ranking. A shared maturational schedule => positive per-person rank correlation.
Null: shuffle that person's own onsets across the categories THEY answered (preserves how many
they answered and their own value distribution exactly).
Then: does an individual's DEVIATION from the schedule mean anything, or is it noise?
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(4711)
exec(open('src/24_attack_rsa.py').read().split('print("=== is onset a proxy')[0])
catcols=[v for _,v in pairs]
V=O.values; mask=~np.isnan(V)
popmean=np.nanmean(V,axis=0)
keep=mask.sum(1)>=6
print(f"people with >=6 onsets: {keep.sum():,}")

obs=[];nul=[]
for i in np.flatnonzero(keep):
    j=np.flatnonzero(mask[i]); v=V[i,j]
    if len(set(v))<3: continue
    obs.append(stats.spearmanr(v,popmean[j]).statistic)
    nul.append(stats.spearmanr(rng.permutation(v),popmean[j]).statistic)
obs=np.array(obs); nul=np.array(nul)
print(f"\nper-person rank agreement with the population schedule")
print(f"  observed : mean {np.nanmean(obs):+.3f}   median {np.nanmedian(obs):+.3f}   sd {np.nanstd(obs):.3f}")
print(f"  null     : mean {np.nanmean(nul):+.3f}   sd {np.nanstd(nul):.3f}")
print(f"  share of people with positive agreement : {np.nanmean(obs>0):.1%}  (null {np.nanmean(nul>0):.1%})")
print(f"  Cohen d vs null : {(np.nanmean(obs)-np.nanmean(nul))/np.nanstd(nul):.2f}")

# is a person's DEVIATION from the schedule a stable trait or noise?
dev=[];ids=[]
for i in np.flatnonzero(keep):
    j=np.flatnonzero(mask[i]); v=V[i,j]
    if len(j)<8 or len(set(v))<3: continue
    a,b=j[::2],j[1::2]
    ra=stats.spearmanr(V[i,a],popmean[a]).statistic; rb=stats.spearmanr(V[i,b],popmean[b]).statistic
    if not (np.isnan(ra) or np.isnan(rb)): dev.append((ra,rb)); ids.append(i)
dev=np.array(dev)
r=np.corrcoef(dev[:,0],dev[:,1])[0,1]; sb=2*r/(1+r)
print(f"\nis 'how closely you follow the schedule' a stable individual trait?")
print(f"  split-half over categories: r={r:+.3f}  Spearman-Brown={sb:+.3f}  n={len(dev):,}")
print(f"  -> {'a real individual difference' if sb>0.3 else 'mostly measurement noise; do not treat as a trait'}")

# does schedule-adherence relate to anything?
sched=pd.Series(dev.mean(1),index=np.array(ids))
for c in ['biomale','opennessvariable','neuroticismvariable']:
    y=pd.to_numeric(df[c],errors='coerce').reindex(sched.index)
    m=y.notna()&sched.notna()
    if m.sum()>500: print(f"  corr(schedule adherence, {c:22s}) = {np.corrcoef(sched[m],y[m])[0,1]:+.3f}")
tot=pd.Series(np.nanmean(V,axis=1)).reindex(sched.index)
m=tot.notna()&sched.notna()
print(f"  corr(schedule adherence, own mean onset  ) = {np.corrcoef(sched[m],tot[m])[0,1]:+.3f}")
