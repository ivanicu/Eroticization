import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: RATING_0_5 is 68 columns, the largest never-anchored family and the most-used one
   outside the coordinates. Is it ONE direction?

Worlds
  A  one direction: every column agrees in sign against the two count anchors
     -> the family can be averaged, and a family score is meaningful.
  B  mixed sign: it is NOT one direction -> any quantity that averages them is wrong.
     Sub-worlds inside B (must be distinguished, #423e NEXT):
       B1 coding-reversed  (same content, opposite numeric direction)
       B2 content-reversed (the item asks about DISLIKE -- sign is correct and meaningful)

Method = #419b / #423c: anchor against Totalsexacts and totalfetishcategory, whose
direction is fixed by being counts and which define neither family.
Null = permute-within-mask (#385c: permuting an array with NaN contaminates the null).

Pre-registered (branch-labelled per #379c):
  PRECONDITION (#392e) : print value set + mode for EVERY column BEFORE it enters anything.
                         Columns whose value set is {0,1} are NOT 0-5 ratings -> excluded,
                         and the exclusion is reported, not silent.
  CONTROL              : NEG_FIB must reproduce -0.4502 / -0.1863 in this script.
  KILL-A               : any column significant in the MINORITY sign -> A dead, B lives.
  MULTIPLICITY (#379)  : 2 anchors x k columns, reported over the WHOLE family with a
                         family-wise threshold, not per-column.
FRONTIER: the outcome decides whether a 68-column family may ever be averaged.
"""
import pandas as pd, numpy as np, json
from lib.gates import Gate
from lib.nulls import perm_in

g=Gate("R468 anchor RATING_0_5")
inv=pd.read_csv('data/derived/inventory.csv')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda c: pd.to_numeric(df[c],errors='coerce')
A1,A2='Totalsexacts','totalfetishcategory'; a1,a2=num(A1),num(A2)

ALL=[c for c in inv[inv.kind=='RATING_0_5'].col if c in df.columns]
# ---- PRECONDITION #392e : value set + mode for every column, BEFORE use
vs=[]
for c in ALL:
    v=num(c); u=sorted(v.dropna().unique())
    vs.append(dict(col=c[:80], n=int(v.notna().sum()), nuniq=len(u),
                   vmin=(u[0] if u else np.nan), vmax=(u[-1] if u else np.nan),
                   vals=str(u[:8]), mode=float(v.mode().iloc[0]) if v.notna().any() else np.nan))
V=pd.DataFrame(vs); V.to_csv(HERE/'results/value_sets.csv',index=False)
BINARY=[r['col'] for r in vs if r['nuniq']<=2]
EMPTY =[r['col'] for r in vs if r['n']==0]
print(f"68 列的值集已写盘。⚠ 值集只有 <=2 个值(不是 0-5 评分)= {len(BINARY)} 列 -> 排除")
for b in BINARY[:12]: print(f"    排除 · {b}")
if EMPTY: print(f"⚠ 全空 = {len(EMPTY)} 列 -> 排除")
COLS=[c for c in ALL if c[:80] not in set(BINARY)|set(EMPTY)]
print(f"进入标定的列数 = {len(COLS)} / {len(ALL)}   (覆盖率 {len(COLS)/len(ALL):.0%})")

def anchor(s, a, seed0, nperm=200):
    m=s.notna()&a.notna()
    if m.sum()<200: return None
    r=float(np.corrcoef(s[m],a[m])[0,1])
    nl=[float(np.corrcoef(perm_in(s.values,m.values,seed=seed0+i)[m.values],a[m])[0,1])
        for i in range(nperm)]
    sd=float(np.std(nl))
    return dict(n=int(m.sum()), r=r, null_sd=sd, z=(r/sd if sd>0 else np.nan))

rows=[]
for i,c in enumerate(COLS):
    s=num(c)
    for an,a in ((A1,a1),(A2,a2)):
        d=anchor(s,a,seed0=7000+37*i)
        if d: rows.append(dict(col=c[:80], anchor=an, **d))
R=pd.DataFrame(rows); R.to_csv(HERE/'results/anchors.csv',index=False)

# ---- CONTROL: the anchors must work, in this script
NEG=[c for c in inv[inv.kind=='RATING_NEG_FIB'].col if c in df.columns]
ns=pd.concat([num(c) for c in NEG],axis=1); nk=ns.notna().sum(1); nm=ns.mean(1).where(nk>=2)
ctl=[anchor(nm,a,seed0=999) for a in (a1,a2)]
ok_ctl=all(d is not None and d['r']<0 and abs(d['z'])>10 for d in ctl)
g.asserted("CONTROL NEG_FIB reproduces #392c in this script",
           ok_ctl, f"r = {[round(d['r'],4) for d in ctl]}", kind="control")

# ---- MULTIPLICITY over the WHOLE family (#379)
K=len(R); thr=float(np.abs(np.percentile(np.random.default_rng(5).normal(size=200000),
                                         [100*(1-0.05/(2*K)),100*(0.05/(2*K))])).max())
R['sig']=R.z.abs()>thr
pos=R[(R.sig)&(R.r>0)]; neg=R[(R.sig)&(R.r<0)]
maj,mino=(pos,neg) if len(pos)>=len(neg) else (neg,pos)
print(f"\n族内阈(Bonferroni, K={K} 个检验) |z| > {thr:.2f}")
print(f"显著且为正 = {len(pos)} · 显著且为负 = {len(neg)} · 不显著 = {int((~R.sig).sum())}")
print(f"多数号 = {'正' if maj is pos else '负'} · **少数号显著列 = {len(mino)}**")

one_direction = len(mino)==0
g.asserted("KILL-A the family is ONE direction (no significant minority-sign column)",
           one_direction, f"minority-sign significant = {len(mino)} of {K}")
g.asserted("coverage reported with the conclusion",
           True, f"{len(COLS)}/{len(ALL)} columns entered ({len(COLS)/len(ALL):.0%})", kind="control")

# ---- B: print the item text of the minority-sign columns (#423e NEXT: B1 vs B2)
if len(mino):
    both=mino.groupby('col').size(); both=both[both>=1].index.tolist()
    print(f"\n⚠ 少数号显著的题面({len(both)} 题)—— 编码反 vs 内容反,必须看题面:")
    out=[]
    for c in both:
        rr=mino[mino.col==c]
        print(f"  · r = {list(np.round(rr.r,3))}  |  {c}")
        out.append(dict(col=c, r=list(np.round(rr.r,4))))
    json.dump(out, open(HERE/'results/minority_sign_items.json','w'), indent=1, ensure_ascii=False)

verdict = "ONE_DIRECTION" if one_direction else "NOT_ONE_DIRECTION"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict, n_all=len(ALL), n_binary_excluded=len(BINARY),
               n_entered=len(COLS), K=K, thr=thr, n_pos=len(pos), n_neg=len(neg),
               n_minority=len(mino), control=ok_ctl),
          open(HERE/'results/verdict.json','w'), indent=1)
print(g.verdict())
