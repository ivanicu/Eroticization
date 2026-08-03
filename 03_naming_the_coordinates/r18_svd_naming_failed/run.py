import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Name the 4 surviving coordinates. A name is only admissible if it is NOT block-local:
for each factor report (a) how many distinct blocks its top-25 loadings span, and
(b) leave-top-block-out stability -- refit with the single most-loaded block deleted and
correlate person scores. A factor whose identity dies when one block is removed was that
block's local structure wearing a general name.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore'); rng=np.random.default_rng(31415)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])   # reuse loaders

allq=list(B); K=6
pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
ppl,F,L,blkid,optname=factors(allq,K,pool)
print(f"persons {len(ppl):,}  items {L.shape[1]}\n")

for k in range(4):
    l=L[k]; order=np.argsort(l)
    top=np.argsort(-np.abs(l))[:25]
    span=len(set(blkid[top]))
    # leave-top-block-out
    tb=pd.Series(blkid[top]).value_counts().index[0]
    ppl2,F2,_,_,_=factors([q for q in allq if q!=tb],K,pool)
    stab=max(abs(np.corrcoef(F[:,k],F2[:,j])[0,1]) for j in range(K))
    print(f"=== COORDINATE {k+1}   top-25 spans {span} blocks   drop-top-block stability r={stab:.2f} (dropped qi={tb}) ===")
    print("   +  "+" | ".join(f"{optname[i][:30]}" for i in order[::-1][:8]))
    print("   -  "+" | ".join(f"{optname[i][:30]}" for i in order[:8]))
    print()
np.save('data/derived/named_F.npy',F); np.save('data/derived/named_L.npy',L)
pd.DataFrame({'item':optname,'block':blkid,**{f'f{j+1}':L[j] for j in range(K)}}).to_csv('data/derived/named_loadings.csv',index=False)
