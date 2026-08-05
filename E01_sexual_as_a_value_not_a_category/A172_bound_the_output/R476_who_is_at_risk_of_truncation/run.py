import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #432a showed `#398d`'s wording ("read from results/") is a proxy: only 51 of 458 rounds
   read back from results/, and the other 407 do not violate it -- they never touch the
   terminal. The real failure is letting POSSIBLY-TRUNCATED terminal output be evidence, and
   truncation only happens when output is LONG. So: which rounds emit long output?

BOUNDARY, written first:
  IN SCOPE : a run.py that can emit an UNBOUNDED number of lines -- a print inside a loop over
             rows, or a whole-frame .to_string() with no .head(). Length is then data-dependent
             and nobody knows it in advance.
  OUT      : a fixed number of prints, however many. 40 constant lines cannot surprise you.
  OUT      : a print inside a loop over a FIXED small set written in the file (a 2-4 element
             tuple literal) -- bounded by construction.

Worlds
  A  few rounds are unbounded -> the truncation risk is rare and #398d can stay a note.
  B  many are -> the risk is structural, and lib/bounded.show is the repair (#383a: change
     the interface, do not remember harder).

CONTROL : this round itself must use lib.bounded.show for its own table, and its controls()
          must pass 4/4 -- a tool proposed as a repair that is not used by the round proposing
          it is a recommendation, not a repair.
CLOSURE.
"""
import pandas as pd
from lib.gates import Gate
from lib.bounded import show, controls
G=Gate("R476 truncation risk")

c=controls()
G.asserted("CONTROL lib/bounded self-checks pass", all(c),
           f"trunc_reported/full_no_lie/written/same_len = {c}", kind="control")

LOOP=re.compile(r'^\s*for\s+.*:\s*$')
FIXED=re.compile(r'\bin\s*\(\(|\bin\s*\[\(|\bin\s*\(["\']')     # loop over a literal tuple/list
FULLSTR=re.compile(r'\.to_string\(')
HEAD=re.compile(r'\.head\(|\[:\d+\]|nlargest|nsmallest')

rows=[]
for p in sorted(pathlib.Path('.').rglob('run.py')):
    if '.git' in str(p): continue
    L=p.read_text(errors='ignore').split('\n')
    unbounded=0
    for i,l in enumerate(L):
        if 'print' not in l: continue
        if FULLSTR.search(l) and not HEAD.search(l): unbounded+=1; continue
        # a print nested under a for-loop whose iterable is not a literal
        for j in range(i-1,max(-1,i-8),-1):
            if LOOP.match(L[j]):
                if not FIXED.search(L[j]) and not HEAD.search(L[j]): unbounded+=1
                break
            if L[j].strip() and not L[j].startswith((' ','\t')): break
    rows.append(dict(round=p.parts[-2], unbounded=unbounded, at_risk=int(unbounded>0)))
T=pd.DataFrame(rows)
n=len(T); risk=int(T.at_risk.sum())
print(f"全部轮次 **{n}**(覆盖率 100%)· **输出行数不可预先知道的 = {risk}({risk/n:.0%})**\n")
show(T[T.at_risk==1], HERE/'results/at_risk.csv', n=12, sort='unbounded', label="风险轮次")

this_session=['R466','R467','R468','R469','R470','R471','R472','R473','R474','R475']
mine=T[T['round'].str.split('_').str[0].isin(this_session)]
print(f"\n本会话 10 轮里处在风险中的 = **{int(mine.at_risk.sum())} / {len(mine)}**")

G.asserted("coverage reported with the conclusion", True, f"{n}/{n} rounds", kind="control")
G.asserted("KILL the truncation risk is rare (< 10% of rounds)", risk/n < 0.10,
           f"at risk = {risk}/{n} = {risk/n:.0%}")
verdict = "RARE" if risk/n<0.10 else "STRUCTURAL"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,n_rounds=n,n_at_risk=risk,frac=risk/n,
               this_session_at_risk=int(mine.at_risk.sum()),this_session=len(mine),
               controls=list(map(bool,c))),
          open(HERE/'results/verdict.json','w'),indent=1)
print(G.verdict())
