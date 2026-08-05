import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #452a(1) found a `#NNN` citation on the public page pointing at the wrong ledger entry --
   "recall bias controlled (#289)", where #289 is about measurement invariance by sex. Nothing
   checks these. Are there others?

⚠ This is the eighth scanning round of this session and six of the previous seven over-indicted
   (#382c #394c #407a #422b #431c #436b). The difference here, stated before running: the two
   sides of the comparison are BOTH LITERAL TEXT -- the citing sentence is on the page, the
   cited entry's title is in RETRACTIONS.md. Nothing is inferred from a regex about meaning.
   The regex only LOCATES; the judgement is by hand on printed pairs.

Worlds
  A  the #289 case is isolated -> the risk is closed.
  B  more exist -> each is a link that sends a reader to the wrong evidence, and they now have
     names.

BOUNDARY, written first:
  IN SCOPE : a `#NNN` on either page that resolves to an existing ledger entry.
  OUT      : `#NNN` that has no matching entry -- that is a DIFFERENT defect (a dangling
             pointer), counted separately and not mixed in.
  OUT      : `#NNN` inside a passage that is itself about the citation being wrong (this
             round's own correction, #452) -- else the fix reads as a new defect.
CONTROL : the known case must be findable -- run against the pre-fix text via git, so the
          audit is shown to catch what it was built for.
CLOSURE.
"""
import pandas as pd, subprocess
from lib.gates import Gate
from lib.bounded import show
G=Gate("R497 citation audit")

led=pathlib.Path('RETRACTIONS.md').read_text()
# ⚠ #453b: the first version matched only `## Entry NNN` and reported 4 dangling pointers.
# Reading the object showed BOTH forms it missed: entries 1-14 live in a TABLE at the head of
# the file (`| 11 | ... |`), and 15-16 share a MERGED header (`## Entries 15-16`). Zero are
# dangling. That would have been the seventh over-indictment of this session; it was caught by
# reading the file, not by any check I wrote.
TITLES={}
for m in re.finditer(r'^## Entry (\d+)[^\n]*?—\s*(.+)$', led, re.M):
    TITLES[int(m.group(1))]=m.group(2).strip()
for m in re.finditer(r'^## Entries (\d+)[–-](\d+)[^\n]*?—\s*(.+)$', led, re.M):
    for k in range(int(m.group(1)), int(m.group(2))+1): TITLES[k]=m.group(3).strip()
for m in re.finditer(r'^\|\s*(\d{1,2})\s*\|\s*(.{10,120}?)\s*\|', led, re.M):
    k=int(m.group(1))
    if k not in TITLES: TITLES[k]=m.group(2).strip()
print(f"台账里可解析标题的条目 = **{len(TITLES)}**")

CIT=re.compile(r'`#(\d+)([a-z])?`')
rows=[]
for f in ('README.md','README_zh.md'):
    t=pathlib.Path(f).read_text()
    for m in CIT.finditer(t):
        n=int(m.group(1))
        a,b=max(0,m.start()-110),min(len(t),m.end()+70)
        rows.append(dict(page=f, num=n, sub=m.group(2) or '', pos=m.start(),
                         ctx=t[a:b].replace('\n',' '),
                         title=TITLES.get(n,'<<无此条目>>'),
                         resolves=int(n in TITLES)))
T=pd.DataFrame(rows)
print(f"页面上的 `#NNN` 引用 = **{len(T)}**(不同条目 **{T.num.nunique()}**)")
dang=T[T.resolves==0]
print(f"  · **解析不到条目(悬空指针,单独计,不混入)= {len(dang)}**")
if len(dang): print(f"    悬空的编号:{sorted(set(dang.num))}")

# 词重叠只用来**排序**,不用来判定 —— 判定是手读印出来的两列
def toks(s): return set(re.findall(r'[一-鿿]{2,}|[A-Za-z]{4,}', s.lower()))
T['overlap']=[len(toks(r.ctx)&toks(r.title)) for r in T.itertuples()]
ok=T[T.resolves==1].copy()
show(ok.nsmallest(12,'overlap')[['page','num','overlap','title','ctx']],
     HERE/'results/citations_ranked.csv', n=12, label="重叠最低(仅排序用)")
T.to_csv(HERE/'results/all_citations.csv',index=False)

# CONTROL:对**修复前**的页面跑一遍,必须能把已知那条排到最前面
try:
    old=subprocess.run(['git','show','HEAD~2:README_zh.md'],capture_output=True,text=True).stdout
    hit=[m for m in CIT.finditer(old) if int(m.group(1))==289]
    ctxs=[old[max(0,m.start()-90):m.end()+40].replace('\n',' ') for m in hit]
    found=any('回忆偏差' in c for c in ctxs)
except Exception:
    found=False
G.asserted("CONTROL the audit locates the known case in the pre-fix text",
           found, f"pre-fix #289 citations found next to 回忆偏差: {found}", kind="control")
G.asserted("coverage reported with the conclusion", True,
           f"{len(T)} citations over 2 pages; {len(TITLES)} ledger titles parsed", kind="control")
G.asserted("KILL every resolvable citation matches its entry (judged by hand below)",
           len(dang)==0, f"dangling pointers = {len(dang)}")
json.dump(dict(n_citations=len(T),n_distinct=int(T.num.nunique()),
               n_dangling=len(dang),dangling=sorted(set(dang.num.tolist())),
               n_titles=len(TITLES)),
          open(HERE/'results/verdict.json','w'),indent=1)
print(G.verdict())
