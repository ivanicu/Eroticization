import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Claim: the POWER readout only shows up where the OPTION SET itself varies in role.
Make it a number: score each block's role-variance from its own option text (pure string
rule, no LLM, no hand-coding of the outcome), then correlate with the measured d_power.
Pre-registered direction: positive. A null here kills the interpretation.
"""
import pandas as pd, numpy as np, re
o=pd.read_csv('data/derived/options.csv'); T=pd.read_csv('data/derived/axis_prediction.csv')
SELF=r'\b(myself|my|me|i )\b'; OTH=r'\b(others|other|someone else|them|their)\b'
rows=[]
for qi,g in o.groupby('qi'):
    txt=g.option.str.lower()
    s=txt.str.contains(SELF,regex=True); t=txt.str.contains(OTH,regex=True)
    # role variance = both poles present, scaled by how balanced they are
    ns,nt,k=int(s.sum()),int(t.sum()),len(g)
    bal=(min(ns,nt)/max(ns,nt)) if max(ns,nt)>0 else 0.0
    rows.append(dict(qi=qi, k=k, n_self=ns, n_other=nt,
                     role_var=round((ns+nt)/k*bal,3)))
RV=pd.DataFrame(rows).merge(T[['qi','col','n','d_power','d_both']],on='qi')
RV=RV.sort_values('role_var',ascending=False)
print(RV[['qi','col','k','n_self','n_other','role_var','d_power']].assign(
      col=RV.col.str[:44], d_power=RV.d_power.round(4)).to_string(index=False))
r=np.corrcoef(RV.role_var,RV.d_power)[0,1]
from scipy import stats
rho,p=stats.spearmanr(RV.role_var,RV.d_power)
perm=[np.corrcoef(RV.role_var,np.random.permutation(RV.d_power))[0,1] for _ in range(5000)]
print(f"\nblocks={len(RV)}")
print(f"corr(option-set role variance, POWER-axis predictive gain) : r={r:+.3f}  rho={rho:+.3f}  p_spearman={p:.4f}")
print(f"permutation p (|r|) = {np.mean(np.abs(perm)>=abs(r)):.4f}")
hi=RV[RV.role_var>0]; lo=RV[RV.role_var==0]
print(f"\nblocks whose options DO contrast self/other (n={len(hi)}): median d_power = {hi.d_power.median():+.4f}")
print(f"blocks whose options do NOT             (n={len(lo)}): median d_power = {lo.d_power.median():+.4f}")
print(f"ratio = {hi.d_power.median()/max(abs(lo.d_power.median()),1e-6):.1f}x")
RV.to_csv('data/derived/role_variance.csv',index=False)
