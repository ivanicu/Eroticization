import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #459c found that counting DOMAINS ENTERED predicts feeling able to change while counting
   PICKS does not. The remaining worry is that both are answering behaviour. This round splits
   answering behaviour against itself.

The separator: **width and depth are both produced by answering**, and they can be made
orthogonal-ish by construction --
    width = `ncat`            how many separate domains a person entered
    depth = `PICKS / ncat`    how many options they ticked per domain they entered
If answering behaviour were the common explanation, it should show in BOTH, because it
determines both. If only width shows, answering behaviour is not what is being measured.

Worlds
  A  width only  -> answering behaviour is excluded as the common cause, and #459c becomes a
     statement about people rather than about the survey.
  B  both, or neither -> they cannot be separated here, and #459's reading narrows back to
     "something about how much was answered".

⚠ COLLINEARITY FIRST (#459's NEXT): width and depth are built from the same two numbers, so
   their correlation is reported BEFORE the fits -- if they are near-collinear, "neither
   clears" would just be the two eating each other, not an answer.
⚠ depth is a RATIO whose denominator is a count, positive and far from zero by construction,
   so `share()`'s refusal does not apply -- but the denominator's distribution is printed, as
   that is what makes the ratio admissible.
MULTIPLICITY: 2 quantities x 2 outcomes -> the null of the maximum (#440b).
CONTROL : `ncat`'s coefficients must stay near #459b's (-0.0353 on changeability) when depth
   joins the model, or the model changed more than intended.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R504 wide or deep")
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy(); _COVB=np.array(COVB,dtype=float).copy()
PICKS=np.zeros(NN); _seen=np.zeros(NN,bool)
for _Mb,_ppl in MB:
    PICKS[_ppl]+=_Mb.sum(1); _seen[_ppl]=True
PICKS=np.where(_seen,PICKS,np.nan)
_PICKS=PICKS.copy(); _NCAT=np.array(ncat,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float); BE=np.asarray(OUT['能不能改'],dtype=float)

WIDE=_NCAT
with np.errstate(divide='ignore',invalid='ignore'):
    DEEP=np.where(_NCAT>0, _PICKS/np.maximum(_NCAT,1e-9), np.nan)
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY) & np.isfinite(PH) & np.isfinite(WIDE) & np.isfinite(DEEP)
n=int(MM.sum())
print(f"n = **{n:,}**")
print(f"分母 `ncat` 分布(比值可采信的前提):min **{np.nanmin(WIDE[MM]):.0f}** · "
      f"中位 **{np.nanmedian(WIDE[MM]):.0f}** · max **{np.nanmax(WIDE[MM]):.0f}** —— "
      f"**由构造为正且远离零**")
r_wd=float(np.corrcoef(WIDE[MM],DEEP[MM])[0,1])
print(f"⚠ **共线性先报**:corr(进得广, 进得深) = **{r_wd:+.4f}**")
GATE.asserted("CONTROL width and depth are not near-collinear",
              abs(r_wd)<0.8, f"corr = {r_wd:+.4f}", kind="control")

Q={'进得广 ncat':WIDE, '进得深 PICKS/ncat':DEEP}
OTHERS={'冷门程度 S':A,'广度型 c3⁻':C3,'常规也管用(−五题)':Bv,'色情使用量':PH,'起始年龄':_EARLY}
def coef(focal_key, y, over=None):
    cols=[np.ones(n)]
    src=Q if over is None else over
    cols.append(z(np.asarray(src[focal_key],dtype=float),MM))
    cols += [z(np.asarray(v,dtype=float),MM) for k,v in src.items() if k!=focal_key]
    cols += [z(np.asarray(v,dtype=float),MM) for v in OTHERS.values()]+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])
rows=[dict(quantity=k, b_羞耻=coef(k,SH), b_能不能改=coef(k,BE)) for k in Q]
NP_=400; nul=np.zeros((NP_,4))
for i in range(NP_):
    j=0
    for k in Q:
        pq={kk:(perm_in(np.asarray(v,dtype=float),MM,seed=23000+i) if kk==k else v)
            for kk,v in Q.items()}
        for y in (SH,BE): nul[i,j]=abs(coef(k,y,over=pq)); j+=1
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows)
T['sig_羞耻']=T.b_羞耻.abs()>thr; T['sig_能不能改']=T.b_能不能改.abs()>thr
show(T, HERE/'results/wide_deep.csv', n=4, label="广 vs 深")
print(f"   **族内阈(4 格里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

w=T[T.quantity=='进得广 ncat'].iloc[0]; d=T[T.quantity=='进得深 PICKS/ncat'].iloc[0]
GATE.asserted("CONTROL ncat's coefficient survives depth joining the model",
              abs(w.b_能不能改-(-0.0353))<0.02,
              f"ncat on changeability = {w.b_能不能改:+.4f} vs #459b -0.0353", kind="control")
only_wide = bool(w.sig_能不能改) and not bool(d.sig_能不能改) and not bool(d.sig_羞耻)
GATE.asserted("KILL only width predicts (answering behaviour excluded as common cause)",
              only_wide,
              f"wide {w.b_能不能改:+.4f} sig={bool(w.sig_能不能改)}; "
              f"deep {d.b_能不能改:+.4f} sig={bool(d.sig_能不能改)}")
verdict = "WIDTH_ONLY" if only_wide else "CANNOT_SEPARATE"
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,thr=thr,corr_wide_deep=r_wd,
               rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
