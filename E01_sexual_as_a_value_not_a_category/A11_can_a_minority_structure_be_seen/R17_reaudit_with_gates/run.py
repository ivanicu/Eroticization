import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R17 -- RE-READING THE ARC'S LOAD-BEARING CLAIMS THROUGH lib/gates.py.

#103 built the comparison rules as code and validated them by replaying eight historical failures.
The obvious next question is whether the claims the arc RESTS on survive the instrument that was
built because so many gates did not.

Three claims, re-read from their own committed artifacts -- no recomputation, so nothing can drift:

  #95   the rare-option signal is LICENSED           (A11/R05 grid)
  #99   the distribution is symmetrically WIDER      (A11/R12 grid)
  #100  the trait is RELIABLE and is not pick count  (A11/R13 grid)

This is Closure, labelled as such (P0): it protects existing conclusions rather than separating
worlds. Its value is that the protection is now performed by something that has been shown to catch
my specific errors.
"""
import pandas as pd, numpy as np, sys as _s
_s.path.insert(0,str(ROOT))
from lib.gates import Gate
A=pathlib.Path('E01_sexual_as_a_value_not_a_category/A11_can_a_minority_structure_be_seen')

# ---------------- #95 : the licensing of the rare-option signal ----------------
D5=pd.read_csv(A/'R05_control_at_the_derived_magnitude/results/grid.csv',keep_default_na=False)
D5['S']=pd.to_numeric(D5['S']); D5['qq']=pd.to_numeric(D5['qq'])
def gap(a,b,q):
    r=D5[(D5.world==a)&(D5.qq==q)]['S']; n=D5[(D5.world==b)&(D5.qq==q)]['S']
    return r.mean()-n.mean(), 2*np.sqrt(r.std()**2+n.std()**2)
e0,s0=gap('n0','n0cb',95); e1,s1=gap('n1','n1cb',95); e2,s2=gap('n2','n2cb',95)
er,sr=gap('real','cb',95)
g95=Gate("#95  is the rare-option signal licensed?")
g95.asserted("no-op arm == the real comparison", abs(e0-er)<max(s0,sr),
             f"|{e0:+.4f} - {er:+.4f}| = {abs(e0-er):.4f} < {max(s0,sr):.4f}")
g95.positive_control("plant at 2 swaps/block", planted=e2, floor=e0, spread=s0/2)
g95.no_sign_crossing("plant ladder monotone & one-signed", [e0,e1,e2])
g95.resolvable("real effect at p95", effect=er, spread=sr/2)
print(g95); print()

# ---------------- #99 : symmetric widening ----------------
D12=pd.read_csv(A/'R12_is_it_elevation_or_just_width/results/grid.csv',keep_default_na=False)
for c in D12.columns:
    if c.startswith('p'): D12[c]=pd.to_numeric(D12[c])
G12=D12.groupby('world')[[c for c in D12.columns if c.startswith('p')]].mean()
B12=np.load(A/'R12_is_it_elevation_or_just_width/results/boot.npy')
bs=dict(zip(['p1','p5','p10','p25','p50','p75','p90','p95','p99'],B12.std(0)))
up=(G12.loc['real','p95']-G12.loc['cb','p95'])+(G12.loc['real','p90']-G12.loc['cb','p90'])
dn=(G12.loc['cb','p5']-G12.loc['real','p5'])+(G12.loc['cb','p10']-G12.loc['real','p10'])
g99=Gate("#99  is the widening symmetric?")
g99.resolvable("upper elevation (p90+p95)", effect=up, spread=np.hypot(bs['p90'],bs['p95']))
g99.resolvable("lower depression (p5+p10)", effect=dn, spread=np.hypot(bs['p5'],bs['p10']))
g99.asserted("median unmoved", abs(G12.loc['real','p50']-G12.loc['cb','p50'])<2*bs['p50'],
             f"|{G12.loc['real','p50']-G12.loc['cb','p50']:+.4f}| < {2*bs['p50']:.4f}")
g99.asserted("ratio within 2x of symmetry", 0.5<up/dn<2.0, f"up/dn = {up/dn:.2f}")
print(g99); print()

# ---------------- #100 : the trait is reliable and is not pick count ----------------
D13=pd.read_csv(A/'R13_is_it_a_trait/results/grid.csv',keep_default_na=False)
for c in ['rel','rel_resid','r_half_sd','r_vs_picks']: D13[c]=pd.to_numeric(D13[c])
G13=D13.groupby('world')[['rel','rel_resid','r_half_sd','r_vs_picks']].mean()
g100=Gate("#100  is the trait reliable, and is it pick count?")
g100.negative_control("null reliability (residualised)", null=G13.loc['cb','rel_resid'],
                      effect=G13.loc['real','rel_resid'])
g100.positive_control("planted trait recovered", planted=G13.loc['plant','rel_resid'],
                      floor=G13.loc['cb','rel_resid'], spread=G13.loc['plant','r_half_sd'])
g100.resolvable("real residualised reliability", effect=G13.loc['real','rel_resid'],
                spread=G13.loc['real','r_half_sd'])
g100.covers_every_arm("every arm read on the SAME column", checked=['real','cb','plant'],
                      arms=['real','cb','plant'])
g100.asserted("not pick count: >50% of raw reliability survives",
              G13.loc['real','rel_resid']>0.5*G13.loc['real','rel'],
              f"{G13.loc['real','rel_resid']:.4f} > 0.5*{G13.loc['real','rel']:.4f}")
print(g100); print()
print("="*78)
for nm,g in [('#95  licensed signal',g95),('#99  symmetric widening',g99),('#100 reliable trait',g100)]:
    print(f"  {nm:26s} {'SURVIVES' if g.verdict() else 'DOES NOT SURVIVE'} the instrument built from my own failures")
