import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #429b's three defects share one shape -- all happened at "turn a column into a number".
   Every round in this project opens with pd.to_numeric(..., errors='coerce'). Where else
   does that step silently keep a BIASED SUBSAMPLE?

BOUNDARY, written before the scan:
  IN SCOPE : a column where to_numeric parses SOME but not all of its non-null values --
             drop > 10% AND survivors > 0. That is the silent biased subsample.
  OUT      : a column that parses fully (drop = 0) -- nothing happened.
  OUT      : a column that parses to NOTHING (survivors = 0) -- loud, not silent: every
             downstream statistic is NaN and announces itself. `TotalMentalIllness` is this
             kind, and it is exactly why it was caught.
  Then, separately: is the in-scope column USED by this project?

Worlds
  A  no in-scope column is used -> the risk is closed.
  B  some are -> each is a statistic computed on a self-selected subsample, now named.

CONTROL : `sexcount` is the known case (15,263 non-null -> 4,524 parse = 70% drop) and MUST
          appear IN SCOPE. A scan that misses it is measuring something else.
CLOSURE unless world B fires.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
G=Gate("R474 silent coercion")
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)

blob="\n".join(p.read_text(errors='ignore')
               for p in list(pathlib.Path('.').rglob('run.py'))+list(pathlib.Path('lib').rglob('*.py'))
               if '.git' not in str(p))
pages="\n".join(pathlib.Path(f).read_text() for f in ('README.md','README_zh.md'))

rows=[]
for c in df.columns:
    raw=df[c]; n0=int(raw.notna().sum())
    if n0==0: continue
    n1=int(pd.to_numeric(raw,errors='coerce').notna().sum())
    drop=1-n1/n0
    rows.append(dict(col=c[:100], n_raw=n0, n_num=n1, drop=drop,
                     inscope=int(drop>0.10 and n1>0),
                     used=int((c in blob) or (c in pages))))
T=pd.DataFrame(rows); T.to_csv(HERE/'results/coercion.csv',index=False)

full=int((T['drop']==0).sum()); none_=int((T.n_num==0).sum()); part=int(T.inscope.sum())
print(f"原始表 **{len(df.columns)}** 列,非空非零的 **{len(T)}** 列进入检查(覆盖率 {len(T)/len(df.columns):.0%})")
print(f"  · 完全解析(掉幅 0)        = **{full}**   <- 什么也没发生")
print(f"  · 完全不解析(存活 0)      = **{none_}**  <- **响亮**:下游全 NaN,会自己叫")
print(f"  · **部分解析(掉幅>10% 且存活>0)= {part}**  <- **静默的有偏子样本**")
inuse=T[(T.inscope==1)&(T.used==1)]
print(f"  · **其中项目用过 = {len(inuse)}**\n")
for _,r in inuse.sort_values('drop',ascending=False).iterrows():
    print(f"   掉 {r['drop']:5.1%}  {r.n_raw:>6} -> {r.n_num:>6}   {r.col[:70]}")

ctl=T[T.col=='sexcount']
G.asserted("CONTROL sexcount (the known case) is IN SCOPE",
           bool(len(ctl) and ctl.inscope.iloc[0]==1),
           f"sexcount drop = {ctl['drop'].iloc[0]:.1%}" if len(ctl) else "not found", kind="control")
G.asserted("coverage reported with the conclusion", True,
           f"{len(T)}/{len(df.columns)} columns examined", kind="control")
G.asserted("KILL no partially-parsing column is used by this project",
           len(inuse)==0, f"used in-scope columns = {len(inuse)}")

verdict = "CLOSED" if len(inuse)==0 else "SILENT_SUBSAMPLES_FOUND"
print(f"\n判决 = {verdict}")
inuse.to_csv(HERE/'results/used_partial.csv',index=False)
json.dump(dict(verdict=verdict, n_cols=len(df.columns), n_checked=len(T),
               n_full=full, n_none=none_, n_partial=part, n_used=len(inuse)),
          open(HERE/'results/verdict.json','w'), indent=1)
print(G.verdict())
