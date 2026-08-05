import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: at the most-entered end of the map the rarity-to-shame slope is near zero (and negative at
   the narrowest cuts). Is that a real absence, or two groups of people paying opposite costs
   and cancelling in the average?

Worlds
  A  a real near-zero -> both halves sit near zero; "no shame cost in the crowded domains" is
     a statement about everyone there.
  B  cancellation -> the halves take opposite signs, and the average is an artefact: two kinds
     of people, opposite costs. That is much the stronger sentence about people.

The split is `c3-`, which is already on this page with its direction established -- nothing new
is invented to make the cut.
⚠ POWER FIRST, AND IT MAY KILL THE QUESTION: halving n doubles the resolution needed. The MDE
   of each half is computed and printed BEFORE any slope is read; if it exceeds the slopes this
   design has ever produced at that end, the question is unanswerable there and that is the
   honest output (#413b).
⚠ The null is a **negative_control**: a single-group slope, and if rarity carried no shame
   signal in those blocks, permuting it gives zero. Zero is the right expectation.
⚠ #477d is now standing practice: any KILL that PASSES states what would have made it fail.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R523 zero or cancellation")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
keep=list(range(NB)); Ra,Rb=prof_(A,keep),prof_(B,keep)
C=np.zeros((NB,NB))
for i in keep:
    for j in keep:
        g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        if g.sum()>200: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
Rr=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
Rr=np.where(np.isfinite(Ra)|np.isfinite(Rb),Rr,np.nan)
Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
nu=(V[:,2][:,None]*Zm).sum(0); de=(Fm*np.abs(V[:,2])[:,None]).sum(0)
_C3RAW=np.where(de>1e-9,nu/np.maximum(de,1e-9),np.nan)
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
C3M = C3 if np.corrcoef(*[x[M&np.isfinite(C3)&np.isfinite(_C3RAW)] for x in (C3,_C3RAW)])[0,1]>0 else C3
order=np.argsort(_rate)

def slope(V, mask, nperm=400, seed=0):
    mm=mask&np.isfinite(V); n=int(mm.sum())
    if n<400: return None
    def fit(v,idx):
        X=np.column_stack([np.ones(int(idx.sum())), z(v,idx), z(ncat,idx), z(AGE,idx)])
        return float(np.linalg.lstsq(X,z(SH,idx),rcond=None)[0][1])
    b=fit(V,mm)
    nul=np.array([fit(perm_in(V,mm,seed=seed+i),mm) for i in range(nperm)])
    thr=float(np.percentile(np.abs(nul),95))
    return dict(n=n,b=b,thr=thr,mde=float(1.96*nul.std()),sig=bool(abs(b)>thr))

rows=[]
for k in (6,8,10):
    sub=_RAR[list(order[-k:])]; cnt=np.isfinite(sub).sum(0)
    Vk=np.where(cnt>=3, np.nanmean(sub,0), np.nan)
    base=M & np.isfinite(AGE) & np.isfinite(Vk) & np.isfinite(C3M)
    med=float(np.nanmedian(C3M[base]))
    for lab,msk in (('全体',base),('`c3⁻` 高半',base&(C3M>=med)),('`c3⁻` 低半',base&(C3M<med))):
        r=slope(Vk,msk,seed=7000+k*7+len(lab))
        if r: rows.append(dict(k=k, half=lab, **r))
T=pd.DataFrame(rows)
print("⚠ **功率先报**(劈半会把分辨率变差,可能直接杀掉这个问题):")
for k in (6,8,10):
    s=T[(T['k']==k)&(T['half']!='全体')]
    print(f"   k={k}:两半 n={list(s['n'])} · **MDE={[round(float(x),4) for x in s['mde']]}** · "
          f"全体 MDE={float(T[(T['k']==k)&(T['half']=='全体')]['mde'].iloc[0]):.4f}")
show(T, HERE/'results/halves.csv', n=12, label="常见端 × `c3⁻` 劈半")

GATE.asserted("CONTROL the pooled arm reproduces #474 (k=6 −0.055, k=8 −0.039, k=10 +0.040)",
              all(abs(float(T[(T['k']==k)&(T['half']=='全体')]['b'].iloc[0])-v)<0.012
                  for k,v in ((6,-0.0548),(8,-0.0386),(10,0.0398))),
              f"pooled = {[round(float(x),4) for x in T[T['half']=='全体']['b']]}", kind="control")
opp=[]
for k in (6,8,10):
    hi=float(T[(T['k']==k)&(T['half']=='`c3⁻` 高半')]['b'].iloc[0])
    lo=float(T[(T['k']==k)&(T['half']=='`c3⁻` 低半')]['b'].iloc[0])
    hs=bool(T[(T['k']==k)&(T['half']=='`c3⁻` 高半')]['sig'].iloc[0])
    ls=bool(T[(T['k']==k)&(T['half']=='`c3⁻` 低半')]['sig'].iloc[0])
    opp.append((k,hi,lo,hi*lo<0,hs,ls))
    print(f"   k={k}:高半 **{hi:+.4f}**{'✅' if hs else ''} · 低半 **{lo:+.4f}**{'✅' if ls else ''} · "
          f"{'**符号相反**' if hi*lo<0 else '同号'}")
cancel = any(o[3] and (o[4] or o[5]) for o in opp)
GATE.asserted("KILL the near-zero is two halves cancelling", cancel,
              f"opposite-signed with at least one half clearing: "
              f"{[(o[0],round(o[1],3),round(o[2],3)) for o in opp if o[3] and (o[4] or o[5])]}")
verdict = "CANCELLATION" if cancel else "REAL_NEAR_ZERO_OR_UNRESOLVED"
print(f"\n判决 = **{verdict}**")
npk,nmiss=GATE.passing_kill_audit({
 "CONTROL the pooled arm reproduces #474 (k=6 −0.055, k=8 −0.039, k=10 +0.040)": None})
json.dump(dict(verdict=verdict,rows=T.to_dict('records'),opposite=opp,
               passing_kills=npk,without_floor=nmiss),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
