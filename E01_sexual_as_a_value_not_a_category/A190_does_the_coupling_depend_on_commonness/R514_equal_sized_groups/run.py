import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #469c's difference grew monotonically with where the split was placed. As the cut moved, the
   "uncommon" group went from 11 blocks to 21 -- so the trend may be a GROUP-SIZE artifact
   rather than a commonness effect: more blocks means a less noisy group score.

Worlds
  A  the difference is stable across equal-sized groups -> the monotone trend WAS group size,
     and the effect's SIZE becomes stable too.
  B  it still grows with how extreme the groups are -> then it is tracking how far apart the
     two ends are, which IS a dose of commonness, and the conclusion strengthens.

Design: take the k MOST-entered and the k LEAST-entered blocks and DISCARD THE MIDDLE, so the
two groups always have the same number of blocks. Sweep k = 6, 8, 10, 12, 14.
⚠ #469d's noise argument does not follow automatically once block counts are equal -- so the
   mean picks per person in each group is MEASURED here rather than asserted.
⚠ Discarding the middle shrinks the usable sample, so n and the MDE are reported PER k before
   the差 is read (#413b).
CONTROL : at the largest k the design must approach #469's middle cut in magnitude.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R514 equal-sized groups")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
rate=np.full(NB,np.nan); RAR=np.full((NB,NN),np.nan); NPICK=np.full((NB,NN),np.nan)
for b,(Mb,ppl) in enumerate(MB):
    rr=-np.log(np.clip(Mb.mean(0),1e-4,1.)); nb=Mb.sum(1)
    rate[b]=float(Mb.mean())
    RAR[b,ppl]=np.where(nb>0,(Mb@rr)/np.maximum(nb,1),np.nan)
    NPICK[b,ppl]=nb
_RAR=RAR.copy(); _rate=rate.copy(); _NPICK=NPICK.copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float)
order=np.argsort(_rate)                      # 低 -> 高
rows=[]
for k in (6,8,10,12,14):
    lo_b=order[:k]; hi_b=order[-k:]
    def agg(sel,need=3):
        sub=_RAR[sel]; cnt=np.isfinite(sub).sum(0)
        return np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    L=agg(lo_b); H=agg(hi_b)
    mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L); n=int(mm.sum())
    if n<800: rows.append(dict(k=k,n=n,note='样本不足')); continue
    X=np.column_stack([np.ones(n), z(H,mm), z(L,mm), z(ncat,mm), z(AGE,mm)])
    b=np.linalg.lstsq(X,z(SH,mm),rcond=None)[0]
    rg=np.random.default_rng(200+k); idx=np.flatnonzero(mm); bs=[]
    for _ in range(300):
        take=rg.choice(idx,len(idx),replace=True)
        m2=np.zeros(len(mm),bool); m2[np.unique(take)]=True; kk=int(m2.sum())
        X2=np.column_stack([np.ones(kk), z(H,m2), z(L,m2), z(ncat,m2), z(AGE,m2)])
        bb=np.linalg.lstsq(X2,z(SH,m2),rcond=None)[0]; bs.append(float(bb[2]-bb[1]))
    bs=np.array(bs); lo_ci,hi_ci=np.percentile(bs,[2.5,97.5])
    pk_h=float(np.nanmean(_NPICK[hi_b][:,mm])); pk_l=float(np.nanmean(_NPICK[lo_b][:,mm]))
    rows.append(dict(k=k, n=n, b_common=float(b[1]), b_uncommon=float(b[2]),
                     diff=float(b[2]-b[1]), lo=lo_ci, hi=hi_ci, mde=float(1.96*bs.std()),
                     picks_common=pk_h, picks_uncommon=pk_l,
                     sig=bool((lo_ci>0)==(hi_ci>0))))
T=pd.DataFrame(rows)
show(T[['k','n','b_common','b_uncommon','diff','lo','hi','mde','sig']],
     HERE/'results/equal_k.csv', n=6, label="等块数 k 扫描")
print("\n⚠ **`#469d` 的噪声论证在等块数下要重新测,不能沿用**:")
for _,r in T.iterrows():
    if 'picks_common' in r and r.picks_common==r.picks_common:
        print(f"   k={int(r.k)}:每人每块平均勾选项 常见 **{r.picks_common:.2f}** vs "
              f"罕见 **{r.picks_uncommon:.2f}**")
ok=T[T['diff'].notna()]
GATE.asserted("CONTROL the largest k approaches #469's middle-cut magnitude",
              abs(float(ok.iloc[-1]['diff'])-0.0534)<0.06,
              f"k={int(ok.iloc[-1].k)} diff = {float(ok.iloc[-1]['diff']):+.4f} vs #469 +0.0534",
              kind="control")
allpos=bool((ok['diff']>0).all()); nsig=int(ok.sig.sum())
sp=float(ok['diff'].max()-ok['diff'].min())
print(f"\n所有 k 的差都为正 = **{allpos}** · 区间不含零的 k = **{nsig}/{len(ok)}** · "
      f"差的极差 = **{sp:.4f}**")
stable = allpos and sp<0.04
GATE.asserted("KILL the difference is stable across equal-sized groups (world A)",
              stable, f"all positive={allpos}, spread={sp:.4f}")
verdict = "STABLE" if stable else ("DOSE_OF_EXTREMITY" if allpos else "UNSTABLE")
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,rows=T.to_dict('records'),spread=sp,all_positive=allpos),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
