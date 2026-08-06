import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
Q: RATING_BINNED_FIB {0,1,2,3,5,8} is the MIRROR of RATING_NEG_FIB {0,-1,-2,-3,-5,-8},
   whose -8 = most erotic (#392c, read backwards for four rounds).
   Is |value| = intensity in BOTH -- i.e. do the two families disagree in NUMERIC sign
   about the same construct?

Method = #419b's: anchor against two counts whose direction is fixed BY BEING A COUNT
(Totalsexacts, totalfetishcategory) and which play no part in either family's definition.

Worlds
  A  mirror-coded: |v| = intensity in both -> BINNED_FIB corr with counts is POSITIVE,
     NEG_FIB is NEGATIVE. Numeric use of the two families gives OPPOSITE signs.
  B  both ascend in the same numeric direction -> both same sign vs the counts.
  C  BINNED_FIB is not an arousal-intensity scale at all -> no consistent sign.

Pre-registered (branch-labelled per #379c):
  KILL-A  : if BINNED_FIB's mean does NOT correlate positively with BOTH counts -> A dead.
  KILL-B  : if the two families' anchor correlations have the SAME sign -> B lives, A dead.
  CONTROL : NEG_FIB must reproduce #392c's negative sign on the same two anchors, on this
            data, in this script -- the positive control that the anchors work at all.
This is FRONTIER: outcome changes whether a whole 7-column family can enter a quantity.
"""
import pandas as pd, numpy as np, json
from lib.gates import Gate
from lib.nulls import perm_in

g=Gate("R467 mirror FIB")
inv=pd.read_csv('data/derived/inventory.csv')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda c: pd.to_numeric(df[c],errors='coerce')

BIN=[c for c in inv[inv.kind=='RATING_BINNED_FIB'].col if c in df.columns]
NEG=[c for c in inv[inv.kind=='RATING_NEG_FIB'].col if c in df.columns]
A1,A2='Totalsexacts','totalfetishcategory'
a1,a2=num(A1),num(A2)

# #392e: print value set / mode BEFORE anything enters a quantity
rows=[]
for fam,cols in (('BINNED_FIB',BIN),('NEG_FIB',NEG)):
    for c in cols:
        v=num(c); ok=v.notna()
        rows.append(dict(family=fam,col=c[:70],n=int(ok.sum()),
                         vals=str(sorted(v.dropna().unique())),
                         mode=float(v.mode().iloc[0]) if ok.sum() else np.nan))
pd.DataFrame(rows).to_csv(f'{pathlib.Path(__file__).parent}/results/value_sets.csv',index=False)

def famscore(cols):
    M=pd.concat([num(c) for c in cols],axis=1)
    return M.mean(axis=1), M.notna().sum(axis=1)

res=[]
for fam,cols in (('BINNED_FIB',BIN),('NEG_FIB',NEG)):
    s,k=famscore(cols)
    for an,a in ((A1,a1),(A2,a2)):
        m=s.notna()&a.notna()&(k>=2)
        r=float(np.corrcoef(s[m],a[m])[0,1])
        # offset null: permute the family score within the SAME mask (#385c/#438)
        nl=[float(np.corrcoef(perm_in(s.values,m.values,seed=1000+i)[m.values],a[m])[0,1])
            for i in range(200)]
        sd=float(np.std(nl))
        res.append(dict(family=fam,anchor=an,n=int(m.sum()),r=r,
                        null_sd=sd,z=r/sd if sd>0 else np.nan))
R=pd.DataFrame(res); R.to_csv(f'{pathlib.Path(__file__).parent}/results/anchors.csv',index=False)
print(R.to_string(index=False))

b=R[R.family=='BINNED_FIB']; n=R[R.family=='NEG_FIB']
bpos=bool((b.r>0).all()); nneg=bool((n.r<0).all())
same=bool(np.sign(b.r.mean())==np.sign(n.r.mean()))

g.asserted("CONTROL NEG_FIB reproduces #392c's negative sign on both anchors",
           nneg, f"NEG_FIB r = {list(np.round(n.r,4))}", kind="control")
g.asserted("KILL-A BINNED_FIB positive on BOTH counts",
           bpos, f"BINNED_FIB r = {list(np.round(b.r,4))}")
g.asserted("KILL-B the two families do NOT share a sign",
           not same, f"mean signs: BIN {np.sign(b.r.mean()):+.0f} / NEG {np.sign(n.r.mean()):+.0f}")
g.asserted("every |z| > 3 (the anchors are not noise)",
           bool((R.z.abs()>3).all()), f"|z| = {list(np.round(R.z.abs(),1))}", kind="control")

# is the family USED anywhere? (the #422c move: a decidable existence check)
blob="\n".join(p.read_text(errors='ignore')
               for p in list(pathlib.Path('.').rglob('run.py'))+list(pathlib.Path('lib').rglob('*.py'))
               if '.git' not in str(p))
used=[c for c in BIN if c in blob]
pg="\n".join(pathlib.Path(f).read_text() for f in ('README.md','README_zh.md'))
onpage=[c for c in BIN if c in pg]
print(f"\n用过吗:代码 {len(used)}/{len(BIN)} · 页面 {len(onpage)}/{len(BIN)}")
verdict = ("MIRROR-CODED" if (bpos and nneg and not same) else "UNVERIFIED")
print(f"\n判决 = {verdict}   (代码/页面用量 = {len(used)}/{len(onpage)})")
json.dump(dict(verdict=verdict,bin_pos=bpos,neg_neg=nneg,same_sign=same,
               used_in_code=len(used),used_on_page=len(onpage),n_bin=len(BIN),n_neg=len(NEG)),
          open(f'{pathlib.Path(__file__).parent}/results/verdict.json','w'),indent=1)
print(g.verdict())

# the content of these 7 items is BODY IMAGE ("when I picture my own nude ...").
# now that the direction is anchored, they say something about people:
s,k=famscore(BIN); m=s.notna()&(k>=2)
print(f"\n身体意象族:n={int(m.sum())} · 与性行为计数 +{res[0]['r']:.3f} · 与恋物类别计数 +{res[1]['r']:.3f}")
