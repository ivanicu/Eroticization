"""
Coord 4 reads as receiving-vs-giving, yet |r| with the POWER (submissive/dominant) composite
is only 0.112. If 'who am I in the scene' were one folk axis these would be near-collinear.
Test the three candidate role axes against each other on the same people.
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
exec(open('src/16_dimensionality.py').read().split("allq=list(B)")[0])
allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index))
G=np.load('data/derived/gcca_G.npy')
A=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in A.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=pd.concat([z(pd.to_numeric(A[c],errors='coerce'))*sg[c] for c in pc],axis=1).mean(axis=1)
ex=[c for c in A.columns if 'exhibition' in c][0]; vo=[c for c in A.columns if 'voyeur' in c][0]
GAZE=z(pd.to_numeric(A[ex],errors='coerce'))-z(pd.to_numeric(A[vo],errors='coerce'))
X=pd.DataFrame({'POWER submit/dominate':POWER.reindex(pool).values,
                'GAZE  seen/seeing':GAZE.reindex(pool).values,
                'SUBSTANCE receive/give (coord4)':G[:,3]})
C=X.corr(min_periods=200)
print("=== three candidate 'who am I in the scene' axes, same people ===\n")
print(C.round(3).to_string())
n=X.dropna().shape[0]
print(f"\n  pairwise-complete n (all three) = {n:,}")
off=C.values[np.triu_indices(3,1)]
print(f"  mean |off-diagonal r| = {np.abs(off).mean():.3f}")
print(f"  largest |r| between any two = {np.abs(off).max():.3f}")
print("\n  A single folk 'top/bottom' axis predicts |r| near 1. Observed max:", round(float(np.abs(off).max()),3))
ev=np.linalg.eigvalsh(C.fillna(0).values)[::-1]
print(f"  eigenvalues {np.round(ev,2)}  -> effective dimensionality {(ev.sum()**2)/(ev**2).sum():.2f} of 3")
