import pandas as pd, numpy as np
rng = np.random.default_rng(20260802)
qm = pd.read_csv('data/derived/multiselect_questions.csv')
lg = pd.read_parquet('data/derived/endorsements_long.parquet')

multi = qm[(~qm.single_pick) & (qm.n_options>=8) & (qm.n_respondents>=800)]
print("usable multi-pick blocks:", len(multi), " options:", int(multi.n_options.sum()))

def block_matrix(qi):
    sub = lg[lg.qi==qi]
    people = sorted(sub.person.unique()); opts = sorted(sub.option.unique())
    pi = {p:i for i,p in enumerate(people)}; oi = {o:i for i,o in enumerate(opts)}
    M = np.zeros((len(people), len(opts)), dtype=np.int8)
    M[sub.person.map(pi).values, sub.option.map(oi).values] = 1
    return M, people, opts

def splithalf_person(M, reps=40):
    """Reliability of a PERSON's option profile: split options in half, correlate
    the person's half-profile deviations. Spearman-Brown corrected."""
    n,k = M.shape
    if k < 8: return np.nan
    dev = M - M.mean(0, keepdims=True)          # remove item base rate
    out=[]
    for _ in range(reps):
        idx = rng.permutation(k); a,b = idx[:k//2], idx[k//2:2*(k//2)]
        ra, rb = dev[:,a].mean(1), dev[:,b].mean(1)
        if ra.std()==0 or rb.std()==0: continue
        r = np.corrcoef(ra,rb)[0,1]
        out.append(2*r/(1+r) if r>-1 else np.nan)
    return float(np.nanmedian(out))

def splithalf_item(M, reps=40):
    """Reliability of an ITEM's base rate: split PEOPLE in half, correlate item means."""
    n,k = M.shape
    out=[]
    for _ in range(reps):
        idx = rng.permutation(n); a,b = idx[:n//2], idx[n//2:]
        ma, mb = M[a].mean(0), M[b].mean(0)
        r = np.corrcoef(ma,mb)[0,1]
        out.append(2*r/(1+r))
    return float(np.nanmedian(out))

rows=[]
for _,q in multi.iterrows():
    M,people,opts = block_matrix(q.qi)
    rows.append(dict(qi=q.qi, col=q.col[:58], n=M.shape[0], k=M.shape[1],
                     density=round(float(M.mean()),3),
                     rel_person=round(splithalf_person(M),3),
                     rel_item=round(splithalf_item(M),3)))
res = pd.DataFrame(rows).sort_values('rel_person', ascending=False)
res.to_csv('data/derived/reliability.csv', index=False)
print(res.to_string(index=False))
print()
print("PERSON-profile reliability : median %.3f   [this is the CEILING for individual-level claims]" % res.rel_person.median())
print("ITEM-baserate reliability  : median %.3f   [ceiling for population-level claims]" % res.rel_item.median())
