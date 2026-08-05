import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #441d(2) is the sharpest knife on the board: `c2` correlates -0.38 with the count of
   endorsed sex acts, so "the leading residual structure predicts I-cannot-change-this" may
   just be "being more involved predicts it", restated.

Worlds
  A  it IS involvement : all three coefficients collapse once the count is controlled ->
     #440 and #441 both narrow sharply, and "nothing to do with content" becomes
     "synonymous with involvement".
  B  something else too : coefficients survive -> the residual structure carries something
     beyond involvement, and WHAT IT IS becomes the only large open question on the page.

PRE-REGISTERED PREDICTION (before running): **B in part** -- `c3` survives (it correlates
only +0.46 with the count, not 1) but **`c2` collapses**.
⚠ My predictions were wrong three times in four this session (#433b #434 #437 vs #439 right).
This is written to be killed.

⚠ Adding a control absorbs real signal as well as confound, so BOTH arms are reported --
   before and after -- and the comparison, not the after-number, is the result.
Multiplicity stays the NULL OF THE MAXIMUM over the three components (#440b).
Sign handling unchanged: c1 is not anchorable, so only |coefficient| is used for it.
CONTROL : the "before" arm must reproduce #441's numbers exactly (nothing changed but the
          control set).
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

SXA=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)
MM = MM & np.isfinite(SXA)
def coefabs(v, idx, with_count=False):
    cols=[np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
          z(ncat,idx), z(AGE,idx), z(SHAME,idx), z(ACT,idx)]
    if with_count: cols.append(z(SXA,idx))
    return float(np.linalg.lstsq(np.column_stack(cols),z(BEL,idx),rcond=None)[0][1])

rows=[]
for nm,v in COMP.items():
    ra,na=anchor(v); ok=abs(ra)>=0.10
    gd=MM&np.isfinite(v)
    b0=coefabs(v,gd,False); b1=coefabs(v,gd,True)
    rows.append(dict(comp=nm, anchor_r=ra, anchorable=int(ok), n=int(gd.sum()),
                     before=b0, after=b1, abs_before=abs(b0), abs_after=abs(b1),
                     kept=abs(b1)/max(abs(b0),1e-12),
                     signed=(f"{b1:+.4f}" if ok else "**符号不可用**")))
T=pd.DataFrame(rows)

NP_=400
THR={}
for wc in (False,True):
    nul=np.zeros((NP_,3))
    for i in range(NP_):
        for j,(nm,v) in enumerate(COMP.items()):
            gd=MM&np.isfinite(v)
            nul[i,j]=abs(coefabs(perm_in(v,gd,seed=14000+7*i+j), gd, wc))
    THR[wc]=float(np.percentile(nul.max(1),95))
thr=THR[True]
T['thr_before']=THR[False]; T['thr_after']=THR[True]
T['sig_before']=T.abs_before>THR[False]; T['sig']=T.abs_after>THR[True]
show(T[['comp','anchorable','before','after','kept','sig_before','sig','signed']],
     HERE/'results/components.csv', n=6, label="加计数控制前后")
print(f"   阈:不控计数 **{THR[False]:.5f}** · 控计数 **{THR[True]:.5f}**")
print(f"\n**多重性:三个成分中最大 |系数| 的零分布 95 分位 = {thr:.5f}**(`#440b` 首次前置使用)")

nsig=int(T.sig.sum())
GATE.asserted("CONTROL the before-arm reproduces #441 (c3 = 0.0355)",
              abs(abs(T[T.comp=='c3'].before.iloc[0])-0.035535)<2e-3,
              f"c3 before = {T[T.comp=='c3'].before.iloc[0]:+.5f} vs #441 +0.03554", kind="control")
GATE.asserted("KILL the structure survives the involvement control", nsig>0,
              f"{nsig}/3 still clear after controlling the sex-act count: {list(T[T.sig].comp)}")
verdict = ("BEYOND_INVOLVEMENT" if nsig>0 else "IT_IS_INVOLVEMENT")
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,thr=thr,nsig=nsig,r_c3=r3,anchor_control=r_ctl,
               eig=[float(x) for x in w[:3]], rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
