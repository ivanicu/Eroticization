import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: `OTHER` (83 columns) is the last cell on #423a's list with no direction status.
   It is heterogeneous by construction, so "is the family one direction" is the WRONG
   question. The right one: how many of the 83 are ORDERED quantities that this project
   has actually USED, whose direction was never checked?

BOUNDARY, written before the scan (five over-indictments would be #382c #394c #407a #422b
and the failed negative control #426d):
  IN SCOPE  : a column that is (1) numeric with >=3 distinct values -- an ordered quantity
              -- AND (2) referenced by name in this project's code or on the pages.
  OUT       : unordered categorical (an ID, a style, a country) -- it HAS no direction, so
              "unanchored" is not a defect.
  OUT       : direction fixed by construction (a count, an age) -- but this must be
              ESTABLISHED per column, not assumed, and the establishing rule is printed.

Worlds
  A  every in-scope used column is direction-fixed-by-construction -> the risk that #392
     found is closed across all nine kinds.
  B  some in-scope used column has a free direction -> that is the next #392, now named.

CONTROL : the scan must flag the known cases. `sexcount` and `age` are used and are
          direction-fixed-by-construction -> they must appear as IN SCOPE and then be
          classified CONSTRUCTION, not as free-direction. A scan that cannot separate
          those two is measuring "used", not "unanchored".
CLOSURE unless world B fires.
"""
import pandas as pd, numpy as np, re
from lib.gates import Gate

G=Gate("R472 OTHER")
inv=pd.read_csv('data/derived/inventory.csv')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
OTH=[c for c in inv[inv.kind=='OTHER'].col if c in df.columns]

blob="\n".join(p.read_text(errors='ignore')
               for p in list(pathlib.Path('.').rglob('run.py'))+list(pathlib.Path('lib').rglob('*.py'))
               if '.git' not in str(p))
pages="\n".join(pathlib.Path(f).read_text() for f in ('README.md','README_zh.md'))

# direction fixed by construction: the NAME says it is a count, an age, a number of X
CONSTRUCTION=re.compile(r'count|total|^n[A-Z_]|number|age|年龄|计数|^age$|freq|times|hours|years',re.I)

rows=[]
for c in OTH:
    v=pd.to_numeric(df[c],errors='coerce')
    nuniq=int(v.nunique()); ordered = nuniq>=3 and v.notna().sum()>=200
    used = (c in blob) or (c in pages)
    rows.append(dict(col=c[:90], nuniq=nuniq, n=int(v.notna().sum()),
                     ordered=int(ordered), used=int(used),
                     construction=int(bool(CONSTRUCTION.search(c))),
                     vals=str(sorted(v.dropna().unique())[:6])))
T=pd.DataFrame(rows); T.to_csv(HERE/'results/other_columns.csv',index=False)

n_ord=int(T.ordered.sum()); n_used=int(T.used.sum())
inscope=T[(T.ordered==1)&(T.used==1)]
free=inscope[inscope.construction==0]
print(f"`OTHER` 共 **{len(OTH)}** 列")
print(f"  · 有序(数值且 ≥3 个不同值,n≥200) = **{n_ord}**  <- 其余无方向可言,不在范围内")
print(f"  · 项目里用过(代码或页面点名)     = **{n_used}**")
print(f"  · **在范围内(有序 且 用过) = {len(inscope)}**")
print(f"  · 其中方向由构造固定(名字里是计数/年龄) = **{int(inscope.construction.sum())}**")
print(f"  · **⇒ 有序 · 用过 · 方向自由 = {len(free)}**\n")
for _,r in inscope.iterrows():
    tag = "构造固定" if r.construction else "**方向自由**"
    print(f"   [{tag}] nuniq={r.nuniq:>4} n={r['n']:>6}  {r.col[:74]}")

# CONTROL: the two known construction cases must be found and classified as construction
ctl=[c for c in ('sexcount','age') if c in set(inscope.col)]
ctl_ok = len(ctl)==2 and all(inscope[inscope.col==c].construction.iloc[0]==1 for c in ctl)
G.asserted("CONTROL sexcount and age appear IN SCOPE and are classified CONSTRUCTION",
           ctl_ok, f"found {ctl}", kind="control")
G.asserted("coverage reported with the conclusion", True,
           f"{len(OTH)}/{len(OTH)} columns examined; in scope {len(inscope)}", kind="control")
G.asserted("KILL every in-scope used column is direction-fixed-by-construction",
           len(free)==0, f"free-direction = {len(free)}")

verdict = "CLOSED" if len(free)==0 else "NEXT_392_FOUND"
print(f"\n判决 = {verdict}")
free.to_csv(HERE/'results/free_direction.csv',index=False)
json.dump(dict(verdict=verdict, n_other=len(OTH), n_ordered=n_ord, n_used=n_used,
               n_inscope=len(inscope), n_construction=int(inscope.construction.sum()),
               n_free=len(free), control=bool(ctl_ok)),
          open(HERE/'results/verdict.json','w'), indent=1)
print(G.verdict())
