import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

import pandas as pd, numpy as np, re, json
df = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
inv = pd.read_csv('data/derived/inventory.csv')
ms = inv[inv['kind']=='MULTISELECT']['col'].tolist()

# A multiselect answer is a comma-joined option list. Options themselves contain commas
# inside parentheses -> split only on commas NOT inside parens.
def split_opts(s):
    if not isinstance(s,str): return []
    out, depth, cur = [], 0, []
    for ch in s:
        if ch=='(': depth+=1
        elif ch==')': depth-=1
        if ch==',' and depth==0:
            out.append(''.join(cur).strip()); cur=[]
        else: cur.append(ch)
    if cur: out.append(''.join(cur).strip())
    return [o for o in out if o]

records=[]; qmeta=[]
for qi,c in enumerate(ms):
    vals = df[c].dropna()
    if len(vals)==0: continue
    opts = {}
    for idx,v in vals.items():
        for o in split_opts(v):
            opts.setdefault(o,0); opts[o]+=1
            records.append((idx, qi, o))
    # "most erotic" single-pick questions produce 1 option per person -> flag them
    mean_picks = np.mean([len(split_opts(v)) for v in vals])
    qmeta.append(dict(qi=qi, col=c, n_respondents=int(len(vals)),
                      n_options=len(opts), mean_picks=round(float(mean_picks),2),
                      single_pick=bool(mean_picks<1.02)))
qm = pd.DataFrame(qmeta)
long = pd.DataFrame(records, columns=['person','qi','option'])
long['item'] = long['qi'].astype(str)+'::'+long['option']

qm.to_csv('data/derived/multiselect_questions.csv', index=False)
long.to_parquet('data/derived/endorsements_long.parquet')

print("multiselect questions kept :", len(qm))
print("  of which SINGLE-PICK     :", int(qm.single_pick.sum()), "(forced choice)")
print("  of which MULTI-PICK      :", int((~qm.single_pick).sum()), "(check-all-that-apply)")
print("distinct options (items)   :", long['item'].nunique())
print("endorsement rows           :", len(long))
print("persons appearing          :", long['person'].nunique())
print()
print("respondents per question  -> median %d  min %d  max %d" % (qm.n_respondents.median(), qm.n_respondents.min(), qm.n_respondents.max()))
print("options per question      -> median %d  min %d  max %d" % (qm.n_options.median(), qm.n_options.min(), qm.n_options.max()))
print("\n--- 12 biggest multi-pick blocks ---")
print(qm[~qm.single_pick].nlargest(12,'n_options')[['col','n_respondents','n_options','mean_picks']].to_string(max_colwidth=78))
