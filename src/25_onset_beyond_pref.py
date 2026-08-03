"""
Two remaining alternatives to 'acquired as a unit':
 (T) topical near-synonymy -- secretions/dirtiness, vore/bestiality are the same thing twice.
     Control: partial out lexical overlap of the category NAMES (no LLM, no embedding).
 (R) onset structure is redundant with preference structure and adds nothing.
     Control: residualise the onset-similarity matrix on the preference-similarity matrix and
     test whether the RESIDUAL still has non-random structure. Redundant => nothing left.
Already noticed: secretions+abnormal-body has r_onset=+0.149 with r_pref=-0.019, so the two
matrices are demonstrably not collinear. Quantify it.
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore'); rng=np.random.default_rng(999)
exec(open('src/24_attack_rsa.py').read().split('print("=== is onset a proxy')[0])
Ostr=Ores.copy()
for c in Ores.columns:
    m=Ores[c].notna()&Pres[c].notna()
    X=np.c_[np.ones(m.sum()),Pres.loc[m,c].values]; b,*_=lstsq(X,Ores.loc[m,c].values,rcond=None)
    Ostr.loc[m,c]=Ores.loc[m,c].values-X@b
CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan); NN=np.zeros((k,k))
for i in range(k):
    for j in range(i+1,k):
        m=Ostr.iloc[:,i].notna()&Ostr.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
        if m.sum()<150: continue
        CO[i,j]=CO[j,i]=np.corrcoef(Ostr.iloc[:,i][m],Ostr.iloc[:,j][m])[0,1]
        CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
        NN[i,j]=NN[j,i]=m.sum()
iu=np.triu_indices(k,1); ok=~np.isnan(CO[iu])&~np.isnan(CP[iu])
co,cp=CO[iu][ok],CP[iu][ok]

def toks(s): return set(re.findall(r'[a-z]{4,}',s.lower()))
T=[toks(n) for n in names]
lex=np.array([len(T[i]&T[j])/max(1,len(T[i]|T[j])) for i,j in zip(iu[0][ok],iu[1][ok])])
def partial(x,y,Z):
    Z=np.c_[np.ones(len(x)),Z]
    rx=x-Z@lstsq(Z,x,rcond=None)[0]; ry=y-Z@lstsq(Z,y,rcond=None)[0]
    return np.corrcoef(rx,ry)[0,1]
print("=== (T) topical near-synonymy ===")
print(f"  pairs with any shared content word : {int((lex>0).sum())}/{len(lex)}")
print(f"  RSA                                : {np.corrcoef(co,cp)[0,1]:+.3f}")
print(f"  RSA | lexical overlap partialled   : {partial(co,cp,lex[:,None]):+.3f}")
hi=lex>0; lo=~hi
print(f"  RSA within name-overlapping pairs  : {np.corrcoef(co[hi],cp[hi])[0,1]:+.3f}  (n={hi.sum()})")
print(f"  RSA within NON-overlapping pairs   : {np.corrcoef(co[lo],cp[lo])[0,1]:+.3f}  (n={lo.sum()})")

print("\n=== (R) does onset carry structure preference does NOT? ===")
b,*_=lstsq(np.c_[np.ones(len(cp)),cp],co,rcond=None)
res=co-np.c_[np.ones(len(cp)),cp]@b
Rm=np.zeros((k,k)); Rm[iu[0][ok],iu[1][ok]]=res; Rm=Rm+Rm.T
ev=np.abs(np.linalg.eigvalsh(Rm))[::-1]
nulls=[]
for _ in range(400):
    v=rng.permutation(res); M=np.zeros((k,k)); M[iu[0][ok],iu[1][ok]]=v; M=M+M.T
    nulls.append(np.abs(np.linalg.eigvalsh(M))[::-1][:3])
nulls=np.array(nulls)
print(f"  residual matrix top-3 |eigenvalues| : {np.round(ev[:3],3)}")
print(f"  permuted-residual null (mean)       : {np.round(nulls.mean(0),3)}")
print(f"  z                                   : {np.round((ev[:3]-nulls.mean(0))/nulls.std(0),1)}")
print(f"  sd of residual = {res.std():.3f}  vs sd of raw onset-corr = {co.std():.3f}"
      f"  -> {res.std()/co.std():.0%} of onset structure is NOT explained by preference")
pairs_sorted=np.argsort(-np.abs(res))
print("\n  largest onset-similarity NOT predicted by preference:")
for t in pairs_sorted[:6]:
    i,j=iu[0][ok][t],iu[1][ok][t]
    print(f"    resid={res[t]:+.3f}  r_onset={CO[i,j]:+.3f}  r_pref={CP[i,j]:+.3f}  n={int(NN[i,j])}  {names[i]} + {names[j]}")
