import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: the page carries "less agency" as UNVERIFIED because #428c could not anchor the
   `powerlessness` variable's direction. But #442b later found the rule that produced that
   verdict was itself broken -- a one-anchor bar declares a variable whose direction IS known
   (porn use, +0.099 on that mask) unanchorable. Under the corrected rule -- clears 0.10 on
   EITHER count -- `powerlessness` has never been retried.

⚠ THE THIRD WORLD, WRITTEN BEFORE RUNNING, because #442b's repair does not by itself settle
   direction: an anchor whose own direction is fixed by construction tells you which END of a
   variable is which ONLY IF a prior fixes the expected SIGN of the relation. For extraversion
   that prior existed (extraverts have more partners). For "powerlessness vs how many sex acts
   a person endorses" there is no such prior I can defend. So:

Worlds
  A  anchorable AND a defensible sign prior exists -> the page's UNVERIFIED can be lifted.
  B  not anchorable at all -> the UNVERIFIED stands and its reason is "no anchor has power".
  C  the anchor HAS power but no prior fixes the sign -> the UNVERIFIED stands, and its reason
     changes from "the anchor is too weak" to "nothing fixes which end is which". That is a
     different and more precise statement, and the page should carry the precise one.

CONTROL : both known-direction variables (porn use, bondage average) must clear the corrected
   rule on the same mask -- that is exactly what #442b's repair was for.
CONTROL2: extraversion must still anchor, since #429c established it -- if it does not, the
   mask or the rule changed and nothing below holds.
CLOSURE unless world A fires.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R505 retry powerlessness")
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda c: pd.to_numeric(raw[c],errors='coerce').values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
SXC=raw['sexcount'].astype(str).map(BANDS).values.astype(float)   # #429b 的修法
A1=num('Totalsexacts'); A2=num('totalfetishcategory')
CAND={'powerlessness':num('powerlessnessvariable'),
      '**正对照** pornhabit':num('pornhabit'),
      '**正对照** bondageaverage':num('bondageaverage'),
      '**正对照** extroversion':num('extroversionvariable')}
def r_with(v,a):
    g=np.isfinite(v)&np.isfinite(a)
    return (float(np.corrcoef(v[g],a[g])[0,1]), int(g.sum())) if g.sum()>200 else (np.nan,0)
rows=[]
for k,v in CAND.items():
    r1,n1=r_with(v,A1); r2,n2=r_with(v,A2); r3,n3=r_with(v,SXC)
    # ⚠ #461b: the rule as written used the two COUNTS only, and its own positive control
    # failed -- extraversion reaches just 0.083 on them. #429c anchored it on a THIRD
    # direction-fixed quantity, the banded partner count (+0.146). So the rule is "clears 0.10
    # on ANY anchor whose direction is fixed by construction", and the anchor set is three.
    # A null from an instrument that has not passed its positive control is inadmissible (P5).
    best=max(abs(r1),abs(r2),abs(r3))
    rows.append(dict(variable=k, r_sexacts=r1, r_fetish=r2, r_partners=r3, n=max(n1,n2,n3),
                     anchorable_2rule=bool(best>=0.10), best=best))
T=pd.DataFrame(rows)
show(T, HERE/'results/anchors.csv', n=4, label="两锚判据(全样本)")

pos=T[T.variable.str.startswith('**正对照**')]
GATE.asserted("CONTROL every known-direction variable clears the corrected two-anchor rule",
              bool(pos.anchorable_2rule.all()),
              f"{list(zip(pos.variable, pos.anchorable_2rule))}", kind="control")
ex=T[T.variable=='**正对照** extroversion'].iloc[0]
GATE.asserted("CONTROL2 extraversion still anchors, as #429c established",
              bool(ex.anchorable_2rule), f"extraversion best |r| = {ex.best:.4f}", kind="control")

pw=T[T.variable=='powerlessness'].iloc[0]
print(f"\n`powerlessness`:vs 性行为计数 **{pw.r_sexacts:+.4f}** · vs 恋物类别数 **{pw.r_fetish:+.4f}** "
      f"· vs 性伴数分档 **{pw.r_partners:+.4f}** · 最强 **{pw.best:.4f}**")
has_power = bool(pw.anchorable_2rule)
# ⚠ 世界 C:有功率 ≠ 方向可定。方向需要一个**事先可辩护的符号先验**,而这里没有。
SIGN_PRIOR_EXISTS = False   # 写在跑之前:我找不到「无力感 ↔ 认可的性行为数」的文献符号
print(f"锚有功率 = **{has_power}** · **事先可辩护的符号先验存在 = {SIGN_PRIOR_EXISTS}**")
GATE.asserted("KILL powerlessness can now be anchored in DIRECTION (world A)",
              has_power and SIGN_PRIOR_EXISTS,
              f"power={has_power}, sign prior={SIGN_PRIOR_EXISTS}")
verdict = ("ANCHORED" if (has_power and SIGN_PRIOR_EXISTS)
           else ("POWER_BUT_NO_SIGN" if has_power else "NO_ANCHOR"))
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,rows=T.to_dict('records'),
               sign_prior_exists=SIGN_PRIOR_EXISTS),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
