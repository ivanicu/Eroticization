import pandas as pd, numpy as np
df  = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
qm  = pd.read_csv('data/derived/multiselect_questions.csv')
lg  = pd.read_parquet('data/derived/endorsements_long.parquet')
inv = pd.read_csv('data/derived/inventory.csv')
gate_cols = inv[inv['kind'].isin(['RATING_0_5','RATING_BINNED_FIB'])]['col'].tolist()

seen = lg.groupby('qi')['person'].apply(set)
print("Q: is block entry gated on a parent arousal rating?  (top corr gate per block)\n")
rows=[]
for qi,persons in seen.items():
    if len(persons) > 14000 or len(persons) < 300: continue
    ind = pd.Series(0, index=df.index); ind.loc[list(persons)] = 1
    best=(None,0,0)
    for g in gate_cols:
        v = df[g]
        m = v.notna()
        if m.sum() < 1000: continue
        # P(enter block | gate>0) vs P(enter block | gate==0)
        hi = ind[m & (v>0)].mean(); lo = ind[m & (v==0)].mean()
        if pd.isna(hi) or pd.isna(lo): continue
        if hi-lo > best[1]: best=(g, hi-lo, hi)
    if best[0]: rows.append(dict(qi=qi, n_in_block=len(persons), gate=best[0][:60],
                                 lift=round(best[1],3), p_enter_if_gate_pos=round(best[2],3)))
r = pd.DataFrame(rows).sort_values('lift', ascending=False)
r.to_csv('data/derived/branching.csv', index=False)
print(r.head(18).to_string(index=False))
print("\nblocks with lift > 0.5 :", int((r.lift>0.5).sum()), "of", len(r))
print("median lift            :", round(float(r.lift.median()),3))
