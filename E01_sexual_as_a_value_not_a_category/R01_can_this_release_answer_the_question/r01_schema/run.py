import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

import pandas as pd, numpy as np, json, re, sys
pd.set_option('display.width', 200)
df = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
print("SHAPE:", df.shape)
rows=[]
for c in df.columns:
    s = df[c]
    nn = s.notna().sum()
    u  = s.dropna().unique()
    rows.append(dict(col=c, dtype=str(s.dtype), n=int(nn), pct=round(100*nn/len(df),1),
                     nuniq=int(len(u)), sample=str(sorted(map(str,u))[:6])[:150]))
sch = pd.DataFrame(rows)
sch.to_csv('data/derived/schema.csv', index=False)
print("\n--- value-set families (how many cols share each unique value-set) ---")
fam = {}
for c in df.columns:
    u = tuple(sorted(map(str, df[c].dropna().unique())))
    if len(u) <= 12:
        fam.setdefault(u, []).append(c)
for u, cols in sorted(fam.items(), key=lambda kv: -len(kv[1]))[:12]:
    print(f"\n[{len(cols)} cols] values={list(u)[:12]}")
    print("   e.g.:", "; ".join(cols[:4])[:240])
