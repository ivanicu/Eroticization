import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: in #472's shame arm the COMMON-domain slope is negative at small k (-0.0395 at k=8, -0.0365
   at k=10) and only turns positive at k=14 (+0.0332). Nobody has asked what that minus sign
   is. If it is real, then in the domains everyone enters, being drawn to the rarer options
   goes with LESS shame -- the opposite sign to this page's main axis.

Worlds
  A  really negative -> a reversal inside the most-entered domains, which would also explain
     why the pooled slope is diluted.
  B  inside the null -> small-sample noise; #472's reading is unchanged and the page should not
     mention it.

⚠ THE NULL IS A DIFFERENT KIND FROM THE LAST FOUR ROUNDS, and it is named here rather than
   inherited: those tested a BETWEEN-GROUP DIFFERENCE, where zero is not the natural
   expectation. This is a SINGLE-GROUP SLOPE, and the question "should this zero be zero?"
   answers yes -- if the rarity score carried no shame signal in these blocks, permuting it
   gives exactly zero. So this is a **negative_control**, not an offset_control.
⚠ MDE FIRST (#413b): n is 1,665 at k = 8, so a null here is unreadable without it.
Spec curve over k = 6, 8, 10, 12 for the common group alone.
CONTROL : the uncommon group at the same k must come back POSITIVE and clear -- the same
   instrument on the same people must see the effect it is known to see.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R519 the negative slope in common domains")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
rate=np.full(NB,np.nan); RAR=np.full((NB,NN),np.nan)
for b,(Mb,ppl) in enumerate(MB):
    rr=-np.log(np.clip(Mb.mean(0),1e-4,1.)); nb=Mb.sum(1)
    rate[b]=float(Mb.mean()); RAR[b,ppl]=np.where(nb>0,(Mb@rr)/np.maximum(nb,1),np.nan)
_RAR=RAR.copy(); _rate=rate.copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float)
order=np.argsort(_rate)

def slope(sel, need=3, nperm=400, seed=0):
    sub=_RAR[list(sel)]; cnt=np.isfinite(sub).sum(0)
    V=np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    mm = M & np.isfinite(AGE) & np.isfinite(V); n=int(mm.sum())
    if n<600: return None
    def fit(v,idx):
        X=np.column_stack([np.ones(int(idx.sum())), z(v,idx), z(ncat,idx), z(AGE,idx)])
        return float(np.linalg.lstsq(X,z(SH,idx),rcond=None)[0][1])
    b=fit(V,mm)
    nul=np.array([fit(perm_in(V,mm,seed=seed+i),mm) for i in range(nperm)])
    sd=float(nul.std()); thr=float(np.percentile(np.abs(nul),95))
    return dict(n=n, b=b, null_mean=float(nul.mean()), null_sd=sd,
                thr=thr, mde=float(1.96*sd), sig=bool(abs(b)>thr))

rows=[]
for k in (6,8,10,12):
    c=slope(order[-k:], seed=4000+k); u=slope(order[:k], seed=5000+k)
    if c: rows.append(dict(group='常见领域', k=k, **c))
    if u: rows.append(dict(group='罕见领域(对照)', k=k, **u))
T=pd.DataFrame(rows)
show(T[['group','k','n','b','mde','thr','sig']], HERE/'results/single_group.csv',
     n=8, label="单组斜率(零 = negative_control)")
com=T[T['group']=='常见领域']; unc=T[T['group'].str.startswith('罕见')]
print(f"\n⚠ **MDE 先报**:常见组各 k 的 MDE = {[round(float(x),4) for x in com['mde']]}")
print(f"常见组斜率 = {[round(float(x),4) for x in com['b']]} · 越阈 {int(com['sig'].sum())}/{len(com)}")
print(f"罕见组(对照)斜率 = {[round(float(x),4) for x in unc['b']]} · 越阈 {int(unc['sig'].sum())}/{len(unc)}")
GATE.asserted("CONTROL the uncommon group shows the known positive effect",
              bool((unc['b']>0).all()) and int(unc['sig'].sum())>=len(unc)//2,
              f"uncommon slopes {[round(float(x),4) for x in unc['b']]}, "
              f"{int(unc['sig'].sum())}/{len(unc)} clear", kind="control")
GATE.asserted("CONTROL the null is centred at zero (it is a negative_control)",
              bool(com['null_mean'].abs().max()<0.01),
              f"max |null mean| = {com['null_mean'].abs().max():.5f}", kind="control")
neg_real = bool((com['b']<0).any() and com[com['b']<0]['sig'].any())
GATE.asserted("KILL the common-domain slope is really negative", neg_real,
              f"negative and clearing: "
              f"{[(int(r.k),round(float(r.b),4)) for _,r in com.iterrows() if r.b<0 and r.sig]}")
verdict = "REAL_REVERSAL" if neg_real else "NOISE"
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
