import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #440c says "I cannot change this" is carried by the SHAPE of a person's pattern rather than
   by any content. `c3-` is the THIRD eigenvector of the block x block residual correlation.
   What are the first two, and do they predict the belief too?

Worlds
  A  one contrast          : only the third predicts it -> the "shape" is a specific, unique
     contrast, and asking what it is remains a live question.
  B  the whole structure   : the first two predict it as well -> then the belief tracks the
     STRENGTH of the residual structure, not one contrast, and the page's current wording
     ("the shape of a person's pattern") is too specific and must be narrowed.

⚠ SIGN (#368a -- this project's fifth time): an eigenvector's sign is arbitrary, so a SIGNED
   coefficient on an unanchored one is a coin flip. Each of the three is anchored SEPARATELY
   against a quantity whose direction is fixed by construction (a count). If the anchor has no
   power (|r| < 0.10), the eigenvector is declared NOT ANCHORABLE and only |coefficient| is
   reported -- absolute value is sign-immune.
⚠ MULTIPLICITY: three coefficients -> the bar is the NULL OF THE MAXIMUM over the three
   (#440b, the first time that lesson is used as a PRECONDITION rather than a repair).
CONTROL : the third eigenvector reproduced here must match `#440`'s c3 (r ~ 1).
CONTROL2: the anchor must work on a variable whose direction IS known (`pornhabit`, #428c).
FRONTIER.
"""
import numpy as np, pandas as pd, warnings, json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R485 one contrast or the structure")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])

# same code path as c3_drop, generalised to any component k (P4: reuse, do not rewrite)
keep=list(range(NB))
Ra,Rb=prof_(A,keep),prof_(B,keep)
C=np.zeros((NB,NB))
for i in keep:
    for j in keep:
        mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); w=w[o]; V=V[:,o]
Rr=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
Rr=np.where(np.isfinite(Ra)|np.isfinite(Rb),Rr,np.nan)
Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
def score(k):
    num=(V[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
COMP={f'c{k+1}':score(k) for k in range(3)}
print(f"NB = **{NB}** · 前三个特征值 = {np.round(w[:3],4).tolist()} "
      f"(占迹 {100*w[:3].sum()/w.sum():.1f}%)")

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float); ACT=np.asarray(OUT['实践了多少'],dtype=float)
BEL=np.asarray(OUT['能不能改'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(SHAME) & np.isfinite(BEL) & np.isfinite(ACT)

g3=MM & np.isfinite(COMP['c3']) & np.isfinite(C3)
r3=float(np.corrcoef(-COMP['c3'][g3], C3[g3])[0,1])
GATE.asserted("CONTROL the third component reproduces #440's c3",
              abs(r3)>0.95, f"r = {r3:+.4f}", kind="control")

# ---- anchor each component separately (#368a)
ANC=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)
PH =pd.to_numeric(raw['pornhabit'],errors='coerce').values.astype(float)
def anchor(v):
    mm=MM&np.isfinite(v)&np.isfinite(ANC)
    return float(np.corrcoef(v[mm],ANC[mm])[0,1]), int(mm.sum())
mmp=MM&np.isfinite(PH)&np.isfinite(ANC)
r_ctl=float(np.corrcoef(PH[mmp],ANC[mmp])[0,1])
GATE.asserted("CONTROL2 the anchor works on a variable whose direction is known",
              abs(r_ctl)>=0.10, f"pornhabit vs count anchor r = {r_ctl:+.4f}", kind="control")

def coefabs(v, idx):
    X=np.column_stack([np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
                       z(ncat,idx), z(AGE,idx), z(SHAME,idx), z(ACT,idx)])
    return float(np.linalg.lstsq(X,z(BEL,idx),rcond=None)[0][1])

rows=[]
for nm,v in COMP.items():
    ra,na=anchor(v); ok=abs(ra)>=0.10
    gd=MM&np.isfinite(v); b=coefabs(v,gd)
    rows.append(dict(comp=nm, anchor_r=ra, anchorable=int(ok), n=na,
                     coef=b, abs_coef=abs(b),
                     signed=(f"{b:+.4f}" if ok else "**符号不可用**")))
T=pd.DataFrame(rows)

NP_=400
nul=np.zeros((NP_,3))
for i in range(NP_):
    for j,(nm,v) in enumerate(COMP.items()):
        gd=MM&np.isfinite(v)
        nul[i,j]=abs(coefabs(perm_in(v,gd,seed=13000+7*i+j), gd))
mx=nul.max(1); thr=float(np.percentile(mx,95))
T['sig']=T.abs_coef>thr
show(T[['comp','anchor_r','anchorable','coef','abs_coef','sig','signed']],
     HERE/'results/components.csv', n=6, label="三个成分")
print(f"\n**多重性:三个成分中最大 |系数| 的零分布 95 分位 = {thr:.5f}**(`#440b` 首次前置使用)")

nsig=int(T.sig.sum())
GATE.asserted("KILL only the third component predicts the belief",
              nsig==1 and bool(T[T.comp=='c3'].sig.iloc[0]),
              f"{nsig}/3 clear the null-of-the-maximum: {list(T[T.sig].comp)}")
verdict = ("ONE_CONTRAST" if (nsig==1 and bool(T[T.comp=='c3'].sig.iloc[0]))
           else ("WHOLE_STRUCTURE" if nsig>1 else "NONE_CLEARS"))
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,thr=thr,nsig=nsig,r_c3=r3,anchor_control=r_ctl,
               eig=[float(x) for x in w[:3]], rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
