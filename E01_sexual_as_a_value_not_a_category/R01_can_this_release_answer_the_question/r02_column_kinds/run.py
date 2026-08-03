import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

import pandas as pd, numpy as np, re
df = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
cols = list(df.columns)

def kind(c):
    lc = c.lower()
    u = df[c].dropna().unique()
    su = set(map(str,u))
    if lc.startswith('total'): return 'TOTAL'
    if lc.endswith('most'):    return 'FORCED_CHOICE_MOST'
    if 'how old were you' in lc or 'at what age' in lc: return 'AGE_ONSET'
    if su and su <= {'0.0','1.0','2.0','3.0','4.0','5.0'}: return 'RATING_0_5'
    if su and su <= {'0.0','1.0','2.0','3.0','5.0','8.0'}: return 'RATING_BINNED_FIB'
    if su and su <= {'0.0','-1.0','-2.0','-3.0','-5.0','-8.0'}: return 'RATING_NEG_FIB'
    if su and su <= {'-3.0','-2.0','-1.0','0.0','1.0','2.0','3.0'}: return 'LIKERT_PM3'
    if 'which of the following' in lc or 'which one' in lc: return 'MULTISELECT'
    return 'OTHER'

inv = pd.DataFrame([dict(col=c, kind=kind(c), n=int(df[c].notna().sum()),
                         nuniq=int(df[c].nunique())) for c in cols])
inv.to_csv('data/derived/inventory.csv', index=False)
print(inv.groupby('kind').agg(cols=('col','size'), median_n=('n','median')).sort_values('cols', ascending=False))
print()
for k in ['RATING_0_5','RATING_BINNED_FIB','RATING_NEG_FIB','FORCED_CHOICE_MOST','TOTAL','MULTISELECT']:
    sub = inv[inv['kind']==k]['col'].tolist()
    print(f"\n########## {k}  ({len(sub)}) ##########")
    for c in sub: print("  ", c[:150].replace('\n',' '))
