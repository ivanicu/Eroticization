import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: `pornhabit` is the one variable in #428's list whose direction is now MEASURED
   (#428c: +0.181 / +0.196 on both count anchors, higher = more use). It has never been
   asked a psychological question. Is how much porn a person uses its own dimension, or
   just another reading of how broadly they are involved?

Worlds
  A  a re-reading of involvement : once S, c3-, category count and age are controlled, the
     coefficients collapse toward zero. Then it is not an independent quantity and the page
     gains nothing.
  B  its own dimension          : coefficients survive the same controls -- and the SIGN
     PATTERN across the four outcomes is the sentence about people.

Outcomes and controls are R449's, spliced verbatim, so nothing was chosen for this question:
  outcomes = shame · "could I stop" · therapeutic · how much acted on
  controls = z(S) · -z(five-item) · c3- · category count   (+ age, added here)

Pre-registered BEFORE running:
  KILL     : if all four |z| fall under the family-wise threshold, write "no independent
             contribution" and STOP -- do not chase a third round (frontier ss3).
  CONTROL  : the same fit with `pornhabit` REPLACED by a permutation of itself must give
             |z| under threshold on all four (the null has to be able to say no).
  CONTROL2 : the raw (uncontrolled) association must be non-zero, else there is nothing to
             explain away and "collapse" would be vacuous.
  MULTIPLICITY: 4 outcomes, family-wise threshold from the null's own max.
FRONTIER.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show

GATE=Gate("R477 pornhabit")   # NOT `g` -- exec of R449 rebinds it (#427e)
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values
AGE=pd.to_numeric(raw['age'],errors='coerce').values if raw['age'].dtype!=object else None
if AGE is None or np.isfinite(AGE).sum()<1000:
    AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
    AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)

# #392e precondition, printed before it enters anything
fin=np.isfinite(PH)
print(f"\n`pornhabit`:n={int(fin.sum()):,} · 值集 {sorted(set(PH[fin].tolist()))} · "
      f"众数 **{float(pd.Series(PH[fin]).mode().iloc[0]):g}** · 方向已确认(`#428c`:高 = 用得多)")

MM = M & np.isfinite(PH) & np.isfinite(AGE)
print(f"共同可用 n = **{int(MM.sum()):,}**")

def fit(y, focal):
    X=np.column_stack([np.ones(int(MM.sum())), z(focal,MM), z(A,MM), z(Bv,MM),
                       z(C3,MM), z(ncat,MM), z(AGE,MM)])
    yy=z(np.asarray(y,dtype=float),MM)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); return float(b[1])

NP_=400
rows=[]; nulls={}
for nm,y in OUT.items():
    obs=fit(y,PH)
    nl=np.array([fit(y, perm_in(PH,MM,seed=8800+i)) for i in range(NP_)])
    nulls[nm]=nl
    # raw (uncontrolled) association, for CONTROL2
    yy=np.asarray(y,dtype=float); m2=MM&np.isfinite(yy)
    raw_r=float(np.corrcoef(PH[m2],yy[m2])[0,1])
    rows.append(dict(outcome=nm, b=obs, null_sd=float(nl.std()),
                     z=obs/max(nl.std(),1e-12), raw_r=raw_r))
T=pd.DataFrame(rows)

# family-wise threshold from the null's own max across the four
# ⚠ #433a: the null draws are COEFFICIENTS; T.z is a Z-SCORE. Comparing the two is a unit
# mismatch that passes everything. Standardise the null by each outcome's own null sd FIRST.
NM=np.column_stack([nulls[nm]/max(nulls[nm].std(),1e-12) for nm in OUT])
THR=float(np.percentile(np.abs(NM).max(1),95))
T['sig']=T.z.abs()>THR
print(f"\n族内阈(四个结局的零的最大值的 95 分位)= **{THR:.2f}**")
show(T, HERE/'results/pornhabit.csv', n=8, label="pornhabit")

nsig=int(T.sig.sum())
GATE.asserted("CONTROL2 the uncontrolled association is non-zero",
              bool((T.raw_r.abs()>0.03).any()),
              f"max |raw r| = {T.raw_r.abs().max():.4f}", kind="control")
GATE.asserted("CONTROL the permuted focal cannot clear the threshold",
              bool(np.mean(np.abs(NM).max(1)>THR)<=0.06),
              f"null exceedance = {np.mean(np.abs(NM).max(1)>THR):.3f}", kind="control")
GATE.asserted("KILL pornhabit has an independent contribution", nsig>0,
              f"{nsig}/4 outcomes clear the family-wise threshold")

verdict = "OWN_DIMENSION" if nsig>0 else "NO_INDEPENDENT_CONTRIBUTION"
print(f"\n判决 = {verdict}  ({nsig}/4 越阈)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),thr=THR,nsig=nsig,
               rows=T.to_dict('records')), open(HERE/'results/verdict.json','w'),
          indent=1, default=str)
print(GATE.verdict())
