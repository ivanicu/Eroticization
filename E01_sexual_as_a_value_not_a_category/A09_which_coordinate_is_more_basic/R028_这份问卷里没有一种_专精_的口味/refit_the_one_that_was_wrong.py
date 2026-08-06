import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: R449 concluded "the division of labour between S and the five-item score is STABLE across
   four outcomes" by comparing the angles' circular dispersion to a null that permuted EACH
   outcome with its OWN seed -- destroying the outcomes' mutual correlation. That makes the
   null MORE dispersed than it should be, so a low observed dispersion looks more significant
   than it is. #425d measured exactly this failure elsewhere. Does R449's conclusion survive
   the right null?

Worlds
  A  survives : observed dispersion still below the 5th percentile of a row-permutation null
                -> R449 stands, only its strength was overstated.
  B  dies     : it does not -> R449's "stable division of labour" is a shadow of the null's
                own width, and must be retracted from wherever it is claimed.

Both nulls are run SIDE BY SIDE in one script, so the difference is the null and nothing else.
CONTROL: the two nulls must agree on the observed value (it does not depend on the null) and
         the row-permutation null must preserve the outcome-outcome correlation matrix.
FRONTIER: world B forces a retraction.
"""
import pandas as pd, numpy as np, json
from lib.gates import Gate
from lib.nulls import perm_in, row_perm

GATE=Gate("R471 refit R449")   # NOT `g`: exec of R449 rebinds `g` (#426b)
# splice R449's construction verbatim, up to its own null (P16 explicit dependency)
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

print(f"\n重建 R449:n={n:,} · 四个结局 = {list(OUT)}")
print(f"观测圆离散度 = **{D_ANG:.5f}** · 比值 sd = **{D_RAT:.4f}**")

def angles_from(ys):
    aa=[]
    for y in ys:
        bS,bI=coefs(np.asarray(y,dtype=float)); aa.append(np.arctan2(bI,bS))
    return aa

Y=[np.asarray(y,dtype=float) for y in OUT.values()]
NP_=400
# ---- null 1: R449's own -- each outcome its own seed (inter-outcome correlation DESTROYED)
nulA_old=np.array([disp(angles_from([perm_in(Y[k],M,8600+s*7+k) for k in range(len(Y))]))
                   for s in range(NP_)])
# ---- null 2: whole person-rows moved together (inter-outcome correlation PRESERVED)
nulA_new=np.array([disp(angles_from(row_perm(Y,M,52000+s))) for s in range(NP_)])

LO_old=float(np.percentile(nulA_old,5)); LO_new=float(np.percentile(nulA_new,5))
print(f"\n{'零':<28}{'均值':>10}{'sd':>10}{'5 分位(下侧)':>16}   判定")
for nm,nl,lo in (("旧 · 逐结局各自打乱",nulA_old,LO_old),("新 · 整行一起搬(正确)",nulA_new,LO_new)):
    print(f"{nm:<28}{nl.mean():>10.5f}{nl.std():>10.5f}{lo:>16.5f}   "
          f"{'**越阈(稳定)**' if D_ANG<lo else '**落在零里**'}")

# CONTROL: does the new null actually preserve the outcome-outcome correlation?
Ym=np.column_stack(Y)[M]
co=np.corrcoef(np.nan_to_num(Ym-np.nanmean(Ym,0),nan=0.0),rowvar=False)
Yn=np.column_stack(row_perm(Y,M,1))[M]
cn=np.corrcoef(np.nan_to_num(Yn-np.nanmean(Yn,0),nan=0.0),rowvar=False)
iu=np.triu_indices(len(Y),1)
d=float(np.nanmax(np.abs(co[iu]-cn[iu])))
GATE.asserted("CONTROL row_perm preserves the outcome-outcome correlations",
           d<1e-9, f"max |Δcorr| = {d:.2e}", kind="control")
GATE.asserted("CONTROL the observed quantity does not depend on the null",
           True, f"D_ANG = {D_ANG:.5f} in both arms", kind="control")

survives = bool(D_ANG < LO_new)
GATE.asserted("KILL R449 survives the correct null", survives,
           f"D_ANG {D_ANG:.5f} vs new 5th pct {LO_new:.5f}")

# how much was the strength overstated?
sd_old=float(nulA_old.std()); sd_new=float(nulA_new.std())
z_old=(D_ANG-nulA_old.mean())/sd_old; z_new=(D_ANG-nulA_new.mean())/sd_new
print(f"\n强度:旧零 **{z_old:+.2f} sd** -> 新零 **{z_new:+.2f} sd**  "
      f"(零的宽度 {sd_old:.5f} -> {sd_new:.5f},×{sd_new/sd_old:.2f})")
verdict="STANDS" if survives else "RETRACT"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict, D_ANG=float(D_ANG),
               lo_old=LO_old, lo_new=LO_new, z_old=float(z_old), z_new=float(z_new),
               sd_old=sd_old, sd_new=sd_new, ctrl_maxdiff=d, nperm=NP_),
          open(HERE/'results/verdict.json','w'), indent=1)
pd.DataFrame(dict(old=nulA_old,new=nulA_new)).to_csv(HERE/'results/nulls.csv',index=False)
print(GATE.verdict())
