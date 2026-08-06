import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #470b left one named gap: the difference dips to +0.020 at k = 12 while k = 10 and k = 14
   give +0.058 and +0.061. Is that one block entering, or a random walk?

Worlds
  A  a specific block -> its content is the explanation and can be pointed at (R372's rule
     still binds: pointing is DESCRIPTION, never a construct name).
  B  a random walk -> the dip is noise, already covered by #470b's "size not estimable to
     better than a factor of four", and nothing more is owed.

Two decidable things, not one:
  (1) step sweep k = 10..14 with the NAME of each newly entering block printed;
  (2) leave-one-block-out at k = 12 -- remove each of the 24 blocks in turn and see whether any
      single removal restores the difference. ⚠ That is a MAX over 24, so the bar is the null
      OF THE MAXIMUM (#440b), not a single bootstrap sd.
⚠ THIS IS THE THIRD ROUND ON THIS QUESTION. If it lands UNVERIFIED, frontier ss3 says switch
   direction rather than open a fourth -- written before running, not decided after.
CONTROL : k = 10, 12, 14 must reproduce #470a (+0.0581, +0.0202, +0.0610).
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R515 what happens at k12")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
rate=np.full(NB,np.nan); RAR=np.full((NB,NN),np.nan)
for b,(Mb,ppl) in enumerate(MB):
    rr=-np.log(np.clip(Mb.mean(0),1e-4,1.)); nb=Mb.sum(1)
    rate[b]=float(Mb.mean()); RAR[b,ppl]=np.where(nb>0,(Mb@rr)/np.maximum(nb,1),np.nan)
_RAR=RAR.copy(); _rate=rate.copy()
NAMES=[]
qm=pd.read_csv('data/derived/multiselect_questions.csv')
try:
    keepq=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)].reset_index(drop=True)
    NAMES=[str(keepq.iloc[i]['col'])[:60] for i in range(min(len(keepq),NB))]
except Exception: pass
if len(NAMES)<NB: NAMES=[f'block{i}' for i in range(NB)]
_NAMES=list(NAMES)
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float)
order=np.argsort(_rate)

def diff_of(lo_b, hi_b, need=3):
    def agg(sel):
        sub=_RAR[list(sel)]; cnt=np.isfinite(sub).sum(0)
        return np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    L=agg(lo_b); H=agg(hi_b)
    mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L); n=int(mm.sum())
    if n<600: return np.nan, n
    X=np.column_stack([np.ones(n), z(H,mm), z(L,mm), z(ncat,mm), z(AGE,mm)])
    b=np.linalg.lstsq(X,z(SH,mm),rcond=None)[0]
    return float(b[2]-b[1]), n

rows=[]
for k in range(10,15):
    lo_b=list(order[:k]); hi_b=list(order[-k:])
    d,n=diff_of(lo_b,hi_b)
    new_lo=_NAMES[order[k-1]] if k>10 else '—'
    new_hi=_NAMES[order[-k]] if k>10 else '—'
    rows.append(dict(k=k,n=n,diff=d,new_uncommon=new_lo[:46],new_common=new_hi[:46]))
S=pd.DataFrame(rows); show(S, HERE/'results/step_sweep.csv', n=6, label="逐步进入(点名)")
GATE.asserted("CONTROL k=10/12/14 reproduce #470a",
              all(abs(float(S[S['k']==kk]['diff'].iloc[0])-v)<0.01
                  for kk,v in ((10,0.0581),(12,0.0202),(14,0.0610))),
              f"{[round(float(x),4) for x in S['diff']]} vs #470a 0.0581/…/0.0202/…/0.0610",
              kind="control")

k=12; lo12=list(order[:k]); hi12=list(order[-k:])
base,_=diff_of(lo12,hi12)
loo=[]
for grp,blocks in (('罕见组',lo12),('常见组',hi12)):
    for b in blocks:
        L=[x for x in lo12 if x!=b]; H=[x for x in hi12 if x!=b]
        d,nn=diff_of(L,H)
        loo.append(dict(group=grp, block=b, name=_NAMES[b][:44], diff=d,
                        shift=abs(d-base) if d==d else np.nan))
L=pd.DataFrame(loo)
rg=np.random.default_rng(313); idx=np.flatnonzero(M & np.isfinite(AGE)); nul=[]
for _ in range(200):
    take=rg.choice(idx,len(idx),replace=True)
    m2=np.zeros(len(M),bool); m2[np.unique(take)]=True
    sh=[]
    for b in lo12[:6]+hi12[:6]:
        Lb=[x for x in lo12 if x!=b]; Hb=[x for x in hi12 if x!=b]
        def agg(sel):
            sub=_RAR[list(sel)]; cnt=np.isfinite(sub).sum(0)
            return np.where(cnt>=3, np.nanmean(sub,0), np.nan)
        LL=agg(Lb); HH=agg(Hb); mm=m2&np.isfinite(HH)&np.isfinite(LL); nn=int(mm.sum())
        if nn<600: continue
        X=np.column_stack([np.ones(nn), z(HH,mm), z(LL,mm), z(ncat,mm), z(AGE,mm)])
        bb=np.linalg.lstsq(X,z(SH,mm),rcond=None)[0]; sh.append(abs(float(bb[2]-bb[1])-base))
    if sh: nul.append(max(sh))
nul=np.array(nul); thr=float(np.percentile(nul,95))
show(L.nlargest(6,'shift'), HERE/'results/leave_one_block.csv', n=6, label="k=12 留一块(位移最大)")
print(f"\nk=12 基线差 = **{base:+.4f}** · **留一块最大位移 = {L['shift'].max():.4f}** · "
      f"**最大值零的 95 分位 = {thr:.4f}**")
one_block = bool(L['shift'].max()>thr)
GATE.asserted("KILL one specific block explains the dip", one_block,
              f"max shift {L['shift'].max():.4f} vs null-of-max 95th pct {thr:.4f}")
verdict = "ONE_BLOCK" if one_block else "RANDOM_WALK"
print(f"\n判决 = **{verdict}**")
print(f"⚠ **第三轮**:若为 RANDOM_WALK,按 frontier §3 **换方向,不开第四轮**(已写在跑之前)。")
json.dump(dict(verdict=verdict,base=base,max_shift=float(L['shift'].max()),thr=thr,
               steps=S.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
