import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #424c found 0 of 134 anchor correlations negative and read it as "no specialist taste".
   But the 67 columns are highly intercorrelated -- one shared factor landing positively
   under chance would make ALL of them positive at once. Is 0/134 a finding or a shadow?

Worlds
  A  finding : under a null that PRESERVES the 67-column correlation structure but breaks
               the relation to the anchors, all-positive is rare -> #424c stands.
  B  shadow  : all-positive is common in that null -> #424c must be DOWNGRADED to
               "one direction (the family may be averaged)" and "no specialist taste"
               is withdrawn.

Null = row permutation (lib/nulls.row_perm): whole PERSON rows are moved together, so
every column-column correlation and every missingness pattern is preserved exactly, and
only the pairing with the anchors is destroyed. (Per-column independent shuffling would
break the inter-column correlation -- a null biased toward PASS, #390's lesson.)

Pre-registered BEFORE running (threshold, not chosen after):
  DOWNGRADE if all-positive occurs in > 5% of the null draws.
  CONTROL-1  the null must preserve inter-column correlation: mean |corr| of the 67x67
             matrix must match the observed to < 0.01.
  CONTROL-2  the null must destroy the anchor relation: mean null r must sit at ~0.
This is FRONTIER and its POSITIVE outcome is one I do not want (frontier §3).
"""
import pandas as pd, numpy as np, json
from lib.gates import Gate
from lib.nulls import row_perm

g=Gate("R469 shared-factor null")
inv=pd.read_csv('data/derived/inventory.csv')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda c: pd.to_numeric(df[c],errors='coerce')
A1,A2='Totalsexacts','totalfetishcategory'

# rebuild R468's column set EXACTLY (same precondition), not from memory
ALL=[c for c in inv[inv.kind=='RATING_0_5'].col if c in df.columns]
COLS=[c for c in ALL if pd.to_numeric(df[c],errors='coerce').nunique()>2]
X=np.column_stack([num(c).values for c in COLS])
print(f"重建 R468 的列集:{len(COLS)} / {len(ALL)}  (R468 记录为 67/68)")
assert len(COLS)==67, f"column set drifted: {len(COLS)}"

anchors={A1:num(A1).values, A2:num(A2).values}

def sign_count(Xm):
    """number of (column, anchor) pairs with r<0, and the total tested"""
    neg=0; tot=0
    for a in anchors.values():
        for j in range(Xm.shape[1]):
            m=np.isfinite(Xm[:,j])&np.isfinite(a)
            if m.sum()<200: continue
            r=np.corrcoef(Xm[m,j],a[m])[0,1]
            if np.isfinite(r): tot+=1; neg+= (r<0)
    return neg,tot

obs_neg,obs_tot=sign_count(X)
print(f"观测:负号 {obs_neg} / {obs_tot}")

# ---- the null
NP=200
full=np.ones(X.shape[0],dtype=bool)   # permute over ALL people; row_perm keeps NaN with the row
negs=[]; allpos=0
corr_obs=np.corrcoef(np.nan_to_num(X - np.nanmean(X,0),nan=0.0),rowvar=False)
mac_obs=float(np.nanmean(np.abs(corr_obs[np.triu_indices(len(COLS),1)])))
mac_null=[]
mean_r_null=[]
for i in range(NP):
    cols=row_perm([X[:,j] for j in range(X.shape[1])], full, seed=31000+i)
    Xn=np.column_stack(cols)
    n,t=sign_count(Xn); negs.append(n); allpos += (n==0)
    if i<10:
        cn=np.corrcoef(np.nan_to_num(Xn-np.nanmean(Xn,0),nan=0.0),rowvar=False)
        mac_null.append(float(np.nanmean(np.abs(cn[np.triu_indices(len(COLS),1)]))))
    rs=[]
    for a in anchors.values():
        for j in range(Xn.shape[1]):
            m=np.isfinite(Xn[:,j])&np.isfinite(a)
            if m.sum()>=200:
                r=np.corrcoef(Xn[m,j],a[m])[0,1]
                if np.isfinite(r): rs.append(r)
    mean_r_null.append(float(np.mean(rs)))

negs=np.array(negs)
rate=allpos/NP
print(f"\n零里「全为正」的出现率 = {allpos}/{NP} = **{rate:.1%}**  (预注册阈 5%)")
print(f"零里负号数:min {negs.min()} · 中位 {np.median(negs):.0f} · max {negs.max()}  (观测 {obs_neg})")
print(f"CONTROL-1 列间平均|corr|:观测 {mac_obs:.4f} · 零 {np.mean(mac_null):.4f} · 差 {abs(mac_obs-np.mean(mac_null)):.5f}")
print(f"CONTROL-2 零里平均 r = {np.mean(mean_r_null):+.5f}  (应 ~0)")

g.asserted("CONTROL-1 the null preserves inter-column correlation",
           abs(mac_obs-np.mean(mac_null))<0.01,
           f"|{mac_obs:.4f} - {np.mean(mac_null):.4f}| = {abs(mac_obs-np.mean(mac_null)):.5f}", kind="control")
g.asserted("CONTROL-2 the null destroys the anchor relation",
           abs(np.mean(mean_r_null))<0.01, f"mean null r = {np.mean(mean_r_null):+.5f}", kind="control")

stands = rate<=0.05
g.asserted("KILL the all-positive pattern is rare under the shared-factor null",
           stands, f"all-positive rate = {rate:.1%} (pre-registered <= 5%)")

verdict = "STANDS" if stands else "DOWNGRADE"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict, obs_neg=int(obs_neg), obs_tot=int(obs_tot), nperm=NP,
               allpos_rate=rate, null_neg_min=int(negs.min()),
               null_neg_median=float(np.median(negs)), null_neg_max=int(negs.max()),
               mac_obs=mac_obs, mac_null=float(np.mean(mac_null)),
               mean_r_null=float(np.mean(mean_r_null))),
          open(HERE/'results/verdict.json','w'), indent=1)
pd.DataFrame(dict(perm=range(NP), n_negative=negs)).to_csv(HERE/'results/null_sign_counts.csv',index=False)
print(g.verdict())
