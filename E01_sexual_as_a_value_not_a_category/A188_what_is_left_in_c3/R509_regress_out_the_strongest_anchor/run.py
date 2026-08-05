import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #464a found `c3`'s strongest countable correlate is the endorsed-sex-act count (+0.460),
   the strongest pair in the whole table -- and #421 had already established that this count is
   the ODD ONE OUT among the five breadth measures. So: take that count out of `c3` and see
   what is left. Does the residual still push only shame?

Worlds
  A  it survives -> `c3-` carries something the sex-act count does not, and THAT is what five
     naming attempts should have been aiming at.
  B  it collapses -> `c3-` is largely "how many sex acts a person has endorsed", inverted --
     and that IS a name. The naming attempts may have failed because the answer was too plain.

⚠ WRITTEN BEFORE RUNNING (#464's NEXT demanded it, and #456 is why): r = 0.460 means the count
   explains **0.460^2 = 21%** of `c3`'s variance. So removing it leaves **79%** -- "it will take
   a lot" is wrong as an intuition; the honest pre-registration is that MOST of `c3` survives
   the removal by construction, and a collapse of the COEFFICIENT would therefore be
   informative rather than expected.
Both arms are reported (#442a). Same seven-quantity model and mask as #463, so only the focal
   changes. Null of the maximum over the family.
CONTROL : the raw arm must reproduce #463b (+0.0784 / -0.0188 for `c3-`).
CONTROL2: the residual must actually be orthogonal to the count (corr ~ 0), else the removal
   did not happen.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R509 what is left in c3")
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
def sc(k):
    nu=(V[:,k][:,None]*Zm).sum(0); de=(Fm*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(de>1e-9,nu/np.maximum(de,1e-9),np.nan)
_C1=sc(0).copy(); _C3=sc(2).copy(); _m=m.copy()
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
PH=num('pornhabit'); SXA=num('Totalsexacts'); FET=num('totalfetishcategory')
SH=np.asarray(OUT['羞耻'],dtype=float); BE=np.asarray(OUT['能不能改'],dtype=float)
C1P=-_C1 if np.corrcoef(*[x[_m&np.isfinite(_C1)&np.isfinite(FET)] for x in (_C1,FET)])[0,1]<0 else _C1
C3M=-_C3   # `c3⁻`,页面的约定

MM = M & np.isfinite(AGE) & np.isfinite(PH) & np.isfinite(_EARLY) & np.isfinite(C1P) \
     & np.isfinite(C3M) & np.isfinite(SXA)
n=int(MM.sum())
r0=float(np.corrcoef(C3M[MM],SXA[MM])[0,1])
print(f"n = **{n:,}** · corr(`c3⁻`, 性行为计数) = **{r0:+.4f}** -> "
      f"它解释 `c3⁻` 方差的 **{r0**2:.1%}**,**剩下 {1-r0**2:.1%}**")
print(f"⚠ **事前预期(写在跑之前)**:大部分 `c3⁻` 由构造存活,"
      f"**所以系数若塌掉,那是有信息的,不是意料之中的**")

X0=np.column_stack([np.ones(n), (SXA[MM]-SXA[MM].mean())/SXA[MM].std()])
yv=(C3M[MM]-C3M[MM].mean())/C3M[MM].std()
C3R=np.full(len(MM),np.nan); C3R[MM]=yv-X0@np.linalg.lstsq(X0,yv,rcond=None)[0]
r_check=float(np.corrcoef(C3R[MM],SXA[MM])[0,1])
GATE.asserted("CONTROL2 the removal actually happened", abs(r_check)<1e-8,
              f"corr(residual, count) = {r_check:.2e}", kind="control")

def model(focal_key, focal_vec, y, over=None):
    Q={focal_key:focal_vec, '`c1⁺` 卷入广度':C1P, '冷门程度 S':A, '常规也管用(−五题)':Bv,
       '答题类别数':ncat, '色情使用量':PH, '起始年龄(晚)':_EARLY}
    src=Q if over is None else over
    cols=[np.ones(n), z(np.asarray(src[focal_key],dtype=float),MM)]+[
        z(np.asarray(v,dtype=float),MM) for k,v in src.items() if k!=focal_key]+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])

ARMS={'原始 `c3⁻`':C3M, '**扣掉性行为计数后的 `c3⁻`**':C3R}
rows=[dict(arm=k, b_羞耻=model('f',v,SH), b_能不能改=model('f',v,BE)) for k,v in ARMS.items()]
NP_=400; nul=np.zeros((NP_,4)); j=0
for i in range(NP_):
    j=0
    for k,v in ARMS.items():
        pv=perm_in(np.asarray(v,dtype=float),MM,seed=25000+i)
        for y in (SH,BE): nul[i,j]=abs(model('f',pv,y)); j+=1
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows)
T['sig_羞耻']=T.b_羞耻.abs()>thr; T['sig_能不能改']=T.b_能不能改.abs()>thr
T['keep']=[abs(r.b_羞耻)/max(abs(T.b_羞耻.iloc[0]),1e-12) for _,r in T.iterrows()]
show(T, HERE/'results/two_arms.csv', n=4, label="扣计数前后")
print(f"   **族内阈 = {thr:.5f}**")

raw_arm=T.iloc[0]; res_arm=T.iloc[1]
GATE.asserted("CONTROL the raw arm reproduces #463b",
              abs(raw_arm.b_羞耻-0.0784)<0.01 and abs(raw_arm.b_能不能改-(-0.0188))<0.01,
              f"raw = {raw_arm.b_羞耻:+.4f} / {raw_arm.b_能不能改:+.4f} vs #463b +0.0784 / -0.0188",
              kind="control")
survives = bool(res_arm.sig_羞耻) and not bool(res_arm.sig_能不能改)
print(f"\n扣掉之后:羞耻 **{res_arm.b_羞耻:+.4f}**(保留 **{res_arm.keep:.3f}**)· "
      f"能不能改 **{res_arm.b_能不能改:+.4f}**")
GATE.asserted("KILL the residual still pushes only shame (world A)", survives,
              f"residual: shame {res_arm.b_羞耻:+.4f} sig={bool(res_arm.sig_羞耻)}, "
              f"changeability {res_arm.b_能不能改:+.4f} sig={bool(res_arm.sig_能不能改)}")
verdict = "CARRIES_MORE" if survives else "MOSTLY_THE_COUNT"
print(f"\n判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,n=n,thr=thr,r_c3_count=r0,var_explained=r0**2,
                rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
