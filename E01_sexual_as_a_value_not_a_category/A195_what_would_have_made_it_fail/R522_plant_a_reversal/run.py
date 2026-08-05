import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #477d shipped the rule that a KILL which PASSES must say what result would have made it
   fail. The most recent passing kill on the page is #474's -- the common-domain slope is
   positive and clears at k = 10 and k = 12. It has never been asked that question.

The answer has to be measured, not asserted: plant a REVERSAL of known size into the outcome
and see at what size the design reports one.

Worlds
  A  the design recovers planted reversals and reports them -> #474's positive result is a
     resolving pass, and the size it could have detected can be stated on the page.
  B  it does not -> #474 "passed" partly because the design is insensitive in that direction,
     and its positive half must be narrowed.

⚠ A LADDER, NOT ONE PLANT (#372c(2)): delta = 0, -0.02, -0.05, -0.10 in standardised units. If
   the recovered slope does not move monotonically with delta, the plant is what is being
   measured, not the instrument.
⚠ delta = 0 must reproduce #474's numbers exactly, or the harness is not the same one.
This round is the FIRST real use of `Gate.passing_kill_audit`, and its output is persisted.
CLOSURE (it puts a floor under an existing claim; it opens nothing).
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R522 plant a reversal")
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

def arm(k, delta, nperm=400, seed=0, need=3):
    sub=_RAR[list(order[-k:])]; cnt=np.isfinite(sub).sum(0)
    V=np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    mm = M & np.isfinite(AGE) & np.isfinite(V); n=int(mm.sum())
    y=z(SH,mm) + delta*z(V,mm)                       # 种入已知大小的**反向**效应
    X=np.column_stack([np.ones(n), z(V,mm), z(ncat,mm), z(AGE,mm)])
    b=float(np.linalg.lstsq(X,y,rcond=None)[0][1])
    def fitp(v):
        X2=np.column_stack([np.ones(n), z(v,mm), z(ncat,mm), z(AGE,mm)])
        return float(np.linalg.lstsq(X2,y,rcond=None)[0][1])
    nul=np.array([fitp(perm_in(V,mm,seed=seed+i)) for i in range(nperm)])
    thr=float(np.percentile(np.abs(nul),95))
    return dict(k=k, delta=delta, n=n, slope=b, thr=thr,
                clears=bool(abs(b)>thr), reports_negative=bool(b<0 and abs(b)>thr))

rows=[arm(k,d,seed=6000+int(abs(d)*1000)+k) for k in (10,12) for d in (0.0,-0.02,-0.05,-0.10)]
T=pd.DataFrame(rows)
show(T, HERE/'results/ladder.csv', n=8, label="反向种植阶梯")

z0=T[(T['delta']==0.0)]
GATE.asserted("CONTROL delta=0 reproduces #474 (k=10 +0.0398, k=12 +0.0755)",
              abs(float(z0[z0['k']==10]['slope'].iloc[0])-0.0398)<0.008 and
              abs(float(z0[z0['k']==12]['slope'].iloc[0])-0.0755)<0.008,
              f"delta=0 slopes = {[round(float(x),4) for x in z0['slope']]}", kind="control")
mono=all(list(T[T['k']==k]['slope'])==sorted(T[T['k']==k]['slope'],reverse=True) for k in (10,12))
GATE.asserted("CONTROL the recovered slope moves monotonically with the plant",
              mono, f"k=10 {[round(float(x),4) for x in T[T['k']==10]['slope']]}; "
                    f"k=12 {[round(float(x),4) for x in T[T['k']==12]['slope']]}", kind="control")
det={}
for k in (10,12):
    sub=T[(T['k']==k)&(T['delta']<0)&(T['reports_negative'])]
    det[k]=float(sub['delta'].max()) if len(sub) else None
print(f"\n**能被报成「负」的最小反向种植**:k=10 -> {det[10]} · k=12 -> {det[12]}")
ok=all(v is not None for v in det.values())
GATE.asserted("KILL the design detects a planted reversal", ok,
              f"smallest detected reversal: {det}")
verdict = "RESOLVING_PASS" if ok else "INSENSITIVE"
print(f"\n判决 = **{verdict}**")
npk,nmiss=GATE.passing_kill_audit({
 "KILL the design detects a planted reversal":
   f"若种入 {det[10]} / {det[12]} 大小的反向效应仍报不出负号,这一条就失败",
 "CONTROL delta=0 reproduces #474 (k=10 +0.0398, k=12 +0.0755)": None,
})
json.dump(dict(verdict=verdict,detected=det,rows=T.to_dict('records'),
               passing_kills=npk,without_floor=nmiss),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
