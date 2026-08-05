import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #469/#470 established a MODERATOR -- the rarity-to-shame slope is steeper in domains fewer
   people enter. The page's other outcome, the sense of being able to change what arouses you,
   is uncorrelated with shame (-0.002) and has its own antecedents (#458b). Is the same
   moderator there too?

Worlds
  A  moderated as well -> the moderator is a property of the DOMAIN, not of the shame path.
  B  not moderated -> the moderator belongs to the shame path alone, which adds a moderator-
     level layer to the fork #458b/#463 found. That is the stronger statement about people.

PRE-REGISTERED PREDICTION: **B**. (Ten predictions this session, four right -- written to be
killed.)
⚠ THE DESIGN IS REUSED VERBATIM from #469b/#470a: rarity computed WITHIN block, groups split
   by the block's population pick rate, both group scores in ONE model, and the equal-block-count
   sweep over k. **Only the outcome changes**, so a difference cannot be a difference of method.
⚠ Per #470b, what is reported is DIRECTION STABILITY across the sweep, not one cell's size.
CONTROL : rerunning the shame arm here must reproduce #470a's k-sweep -- if it does not, the
   pipeline drifted and the comparison is void.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R516 is the moderator cross-outcome")
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
YS={'羞耻(对照臂)':np.asarray(OUT['羞耻'],dtype=float),
    '**能不能改(新结局)**':np.asarray(OUT['能不能改'],dtype=float)}
order=np.argsort(_rate)

def run(y,k,need=3,nboot=300,seed=0):
    lo_b=order[:k]; hi_b=order[-k:]
    def agg(sel):
        sub=_RAR[list(sel)]; cnt=np.isfinite(sub).sum(0)
        return np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    L=agg(lo_b); H=agg(hi_b)
    mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L); n=int(mm.sum())
    if n<800: return None
    X=np.column_stack([np.ones(n), z(H,mm), z(L,mm), z(ncat,mm), z(AGE,mm)])
    b=np.linalg.lstsq(X,z(y,mm),rcond=None)[0]
    rg=np.random.default_rng(seed+k); idx=np.flatnonzero(mm); bs=[]
    for _ in range(nboot):
        take=rg.choice(idx,len(idx),replace=True)
        m2=np.zeros(len(mm),bool); m2[np.unique(take)]=True; kk=int(m2.sum())
        X2=np.column_stack([np.ones(kk), z(H,m2), z(L,m2), z(ncat,m2), z(AGE,m2)])
        bb=np.linalg.lstsq(X2,z(y,m2),rcond=None)[0]; bs.append(float(bb[2]-bb[1]))
    bs=np.array(bs); lo,hi=np.percentile(bs,[2.5,97.5])
    return dict(k=k,n=n,b_common=float(b[1]),b_uncommon=float(b[2]),
                diff=float(b[2]-b[1]),lo=lo,hi=hi,sig=bool((lo>0)==(hi>0)))

rows=[]
for nm,y in YS.items():
    for k in (8,10,12,14):
        r=run(y,k,seed=900 if '羞耻' in nm else 1900)
        if r: rows.append(dict(outcome=nm,**r))
T=pd.DataFrame(rows)
show(T[['outcome','k','n','b_common','b_uncommon','diff','lo','hi','sig']],
     HERE/'results/cross_outcome.csv', n=10, label="两个结局 × 等块数 k")

sh=T[T['outcome'].str.startswith('羞耻')]
GATE.asserted("CONTROL the shame arm reproduces #470a's sweep",
              all(abs(float(sh[sh['k']==kk]['diff'].iloc[0])-v)<0.012
                  for kk,v in ((8,0.0765),(10,0.0581),(12,0.0202),(14,0.0610))),
              f"shame diffs = {[round(float(x),4) for x in sh['diff']]}", kind="control")
be=T[T['outcome'].str.contains('能不能改')]
pos_sh=int((sh['diff']>0).sum()); pos_be=int((be['diff']>0).sum())
sig_sh=int(sh['sig'].sum());     sig_be=int(be['sig'].sum())
print(f"\n羞耻臂:差为正 **{pos_sh}/{len(sh)}** · 区间不含零 **{sig_sh}/{len(sh)}** · "
      f"中位差 **{sh['diff'].median():+.4f}**")
print(f"能不能改臂:差为正 **{pos_be}/{len(be)}** · 区间不含零 **{sig_be}/{len(be)}** · "
      f"中位差 **{be['diff'].median():+.4f}**")
moderated = sig_be>=len(be)//2 and pos_be==len(be)
GATE.asserted("KILL the moderator is cross-outcome (world A)", moderated,
              f"changeability: {pos_be}/{len(be)} positive, {sig_be}/{len(be)} significant")
verdict = "CROSS_OUTCOME" if moderated else "SHAME_PATH_ONLY"
print(f"\n判决 = **{verdict}**  (预注册预测 = B / SHAME_PATH_ONLY)")
json.dump(dict(verdict=verdict,rows=T.to_dict('records'),
               shame=dict(pos=pos_sh,sig=sig_sh,med=float(sh['diff'].median())),
               belief=dict(pos=pos_be,sig=sig_be,med=float(be['diff'].median())),
               prediction="B"),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
