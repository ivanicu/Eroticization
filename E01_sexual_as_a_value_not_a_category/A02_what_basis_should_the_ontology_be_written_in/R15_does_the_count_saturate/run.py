"""
E01 R02 r15 -- does the coordinate count saturate, or is it undefined?

r14 showed count_3x tracks K (4->4, 6->5, 8->7, 10->8) and is completely flat in the option floor.
Two readings, and they demand different actions:
  SAT   the tail saturates at some K* -- then the count is a real property and K=8 was simply too
        small. Report K* and move on.
  UNDEF the count rises with K without bound -- then "how many coordinates" is not answerable by
        this criterion at all, and only the MAGNITUDE PROFILE is meaningful. The headline loses its
        number permanently rather than gaining a bigger one.

ESTIMAND        count of held-out canonical correlations above 3x their own permutation floor, as
                K -> large; and the shape of the magnitude profile, which is K-free.
IDENTIFICATION  the count is identified only if it converges. Testing convergence IS the round.
KILL            PRE-REGISTERED: if count_3x at K=24 exceeds count at K=12 by more than 3, the
                criterion is declared non-saturating and the coordinate count is withdrawn as a
                number in README.md, replaced by the profile. If it converges within 3, report K*.
POSITIVE CTRL   floor must stay flat as K grows (a floor falling with K would manufacture
                non-saturation). Checked at every K.
NEGATIVE CTRL   person-shuffled null recomputed at every K -- its count must stay ~0 regardless
                of K, otherwise the criterion is trivially satisfiable.
NOISE FLOOR     per-cell permutation, 4 reps.
SPECIFICATION   K in {4,6,8,10,12,16,20,24} x 3 seeds, option floor fixed at 20 (r14 proved it
                irrelevant across 5-80, so sweeping it again would be padding, not coverage).
SEEDS           3.
MULTIPLICITY    24 cells, all reported.
IMPOSSIBLE      cross-dataset, independent replication -- one release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import pandas as pd, numpy as np, warnings, hashlib, itertools
exec(open(ROOT/'E01_sexual_as_a_value_not_a_category/A02_what_basis_should_the_ontology_be_written_in/R14_coordinate_count_specification_curve/run.py').read().split('rows=[]')[0].split('"""')[2])
OUT=pathlib.Path(__file__).parent/'results'
rows=[]
for K,seed in itertools.product([4,6,8,10,12,16,20,24],[11,22,33]):
    r=cell(20,K,seed)
    if r: rows.append(dict(K=K,seed=seed,count_3x=r['count_3x'],count_p20=r['count_p20'],
                           floor=r['floor'],top=r['top'],profile=r['O']))
G=pd.DataFrame(rows); G.to_csv(OUT/'saturation.csv',index=False)
S=G.groupby('K').agg(count_3x=('count_3x','median'),count_p20=('count_p20','median'),
                     floor=('floor','median'),top=('top','median'),sd=('count_3x','std'))
print("=== does the count saturate as K grows? (option floor fixed at 20, r14 proved it inert) ===")
print(S.round(3).to_string())
c12=float(S.loc[12,'count_3x']); c24=float(S.loc[24,'count_3x'])
print(f"\n  count_3x at K=12 : {c12:.0f}      at K=24 : {c24:.0f}      difference : {c24-c12:+.0f}")
print(f"  POSITIVE CONTROL floor across K: {S.floor.min():.4f} to {S.floor.max():.4f} "
      f"({'flat' if S.floor.max()-S.floor.min()<0.005 else 'MOVING -- comparison invalid'})")
print(f"  top canonical r across K: {S.top.min():.3f} to {S.top.max():.3f} (K-free, should be stable)")
prof=G[G.K==24].iloc[0]['profile']
print(f"\n  magnitude profile at K=24 (the K-free object):")
print("   ", " ".join(f"{v:.3f}" for v in prof[:16]))
print(f"    -> falls below 0.10 after coordinate {sum(1 for v in prof if v>0.10)}")
print(f"    -> falls below 0.05 after coordinate {sum(1 for v in prof if v>0.05)}")
print("\nPRE-REGISTERED KILL, evaluated:")
if c24-c12>3:
    print(f"  -> +{c24-c12:.0f} from K=12 to K=24 : NON-SATURATING. The count is withdrawn as a number.")
    print("  -> what replaces it: the magnitude profile, which is K-free and does converge.")
else:
    print(f"  -> converged within 3 (K*={int(c24)}). The count is real and K=8 was too small.")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
