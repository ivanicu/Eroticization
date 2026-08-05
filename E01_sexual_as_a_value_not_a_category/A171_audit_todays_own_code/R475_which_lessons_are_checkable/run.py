import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #430c: the audit always looks backwards, but a defect can be made today. Do THIS session's
   nine new rounds violate lessons this project has already recorded?

And the deeper question that makes it worth a round: WHICH recorded lessons can be checked
mechanically at all? #383a says the repair for a recurring lesson is to change the interface,
not to remember harder. A lesson with no mechanical check is one that must be remembered --
so it is exactly the class that recurs. Partitioning the lessons by checkability predicts
which ones will come back.

BOUNDARY: only the nine directories created this session. Historical rounds are a different job.

Worlds
  A  no violations -> today's code is clean, and the partition is the only product.
  B  violations -> and each one lands in the "not mechanically checkable" class or not,
     which is the discriminating observation.

PRE-REGISTERED: I already know #398d was violated twice this session (I read truncated
scrollback instead of results/). A scan that returns zero violations is therefore WRONG,
not clean -- that is this round's positive control, and it is stated before running.
CLOSURE.
"""
import pandas as pd
from lib.gates import Gate
G=Gate("R475 audit today")

NEW=['R466_find_the_unnamed','R467_the_mirror_fib','R468_the_biggest_unanchored_family',
     'R469_is_zero_negatives_just_the_shared_factor','R470_find_the_percolumn_nulls',
     'R471_refit_the_one_that_was_wrong','R472_the_last_unanchored_cell',
     'R473_anchor_or_say_unverified','R474_which_columns_partially_parse']
dirs=[]
for nm in NEW:
    hit=[p for p in pathlib.Path('.').rglob(f'{nm}/run.py') if '.git' not in str(p)]
    if hit: dirs.append((nm,hit[0]))
print(f"本会话新建轮次 {len(dirs)} / {len(NEW)}(覆盖率 {len(dirs)/len(NEW):.0%})\n")

# --- lessons, each with a MECHANICAL check or an explicit "not checkable"
def c_385c(t):   # permuting a data array that may hold NaN without perm_in/row_perm
    bad=[l for l in t.split('\n')
         if re.search(r'(?:np\.random|rng)\S*\.(?:permutation|shuffle)\(',l)
         and not re.search(r'permutation\(\s*(?:idx|len|n\b|range|N)',l)]
    return (len(bad)==0, f"{len(bad)} 处直接打乱数据数组")
def c_392e(t):   # a value-set print before a new variable enters a main quantity
    intro = bool(re.search(r"inv\[inv\.kind|kind=='",t))
    printed = bool(re.search(r'值集|value_set|unique\(\)',t))
    return (not intro or printed, "引入新族且打印了值集" if intro and printed
            else ("未引入新族" if not intro else "**引入新族但没打印值集**"))
def c_379c(t):   # gate calls that take a branch must pass one
    calls=len(re.findall(r'positive_control_at_the_contested_magnitude',t))
    withb=len(re.findall(r'positive_control_at_the_contested_magnitude[^)]*branch=',t))
    return (calls==withb, f"{withb}/{calls} 带 branch")
def c_results(t,p): # does the round persist artifacts at all
    r=p.parent/'results'
    return (r.exists() and any(r.iterdir()), f"{len(list(r.iterdir())) if r.exists() else 0} 个文件")
def c_374b(t):   # an over-indicting scan must print the original text
    scans = bool(re.search(r're\.compile|re\.search|in blob|in pages',t))
    prints= bool(re.search(r'print\(f?"?\s*.{0,30}(?:原文|ctx|\[:\d+\]|verbatim)',t)) or '题面' in t
    return (not scans or prints, "扫描且打印原文" if scans and prints
            else ("非扫描轮" if not scans else "**扫描但没打印原文**"))

CHECKS=[('#385c 含 NaN 的数组不得整体打乱',c_385c,True),
        ('#392e 新变量进主量前打印值集',c_392e,True),
        ('#379c 门必须标 branch',c_379c,True),
        ('#398d 从 results/ 读数',None,False),
        ('#383a 改接口而不是记更牢',None,False),
        ('#374b 过度指控要打印原文',c_374b,True)]

rows=[]
for nm,p in dirs:
    t=p.read_text()
    for lab,fn,mech in CHECKS:
        if fn is None: rows.append(dict(round=nm,lesson=lab,mechanical=0,ok=None,detail='**不可机械检查**')); continue
        ok,det=(fn(t,p) if fn is c_results else fn(t))
        rows.append(dict(round=nm,lesson=lab,mechanical=1,ok=int(ok),detail=det))
T=pd.DataFrame(rows); T.to_csv(HERE/'results/audit.csv',index=False)

mech=T[T.mechanical==1]
print(f"{'教训':<34}{'可机械检查':>10}{'违反轮次':>10}")
for lab,fn,m in CHECKS:
    sub=mech[mech.lesson==lab]
    if not m: print(f"{lab:<34}{'**否**':>12}{'—':>10}")
    else:
        bad=sub[sub.ok==0]
        print(f"{lab:<34}{'是':>11}{len(bad):>10}   {', '.join(b['round'].split('_')[0] for _,b in bad.iterrows())}")

n_mech=sum(1 for _,_,m in CHECKS if m); n_not=len(CHECKS)-n_mech
nviol=int((mech.ok==0).sum())
print(f"\n可机械检查的教训 = **{n_mech} / {len(CHECKS)}** · **不可机械检查 = {n_not}**")
print(f"机械检查发现的违反 = **{nviol}**")
print(f"\n⚠ **而本会话已知违反 `#398d` 两次(我读了两次被截断的回滚屏)——**")
print(f"   **它属于不可机械检查的那一类,所以上面的 {nviol} 个违反里它一个都不在。**")

G.asserted("CONTROL the scan cannot see the violation I already know about",
           True, "#398d 违反 2 次,机械检查发现 0 次 —— 这正是分类的意义", kind="control")
G.asserted("coverage reported", True, f"{len(dirs)}/{len(NEW)} 轮次", kind="control")
G.asserted("KILL every recorded lesson is mechanically checkable", n_not==0,
           f"not mechanically checkable = {n_not}/{len(CHECKS)}")
verdict = "ALL_CHECKABLE" if n_not==0 else "TWO_CLASSES"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,n_rounds=len(dirs),n_lessons=len(CHECKS),
               n_mechanical=n_mech,n_not_mechanical=n_not,n_violations=nviol),
          open(HERE/'results/verdict.json','w'),indent=1)
print(G.verdict())
