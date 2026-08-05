import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #465c shrank the naming target to "the 87% that survives removing the sex-act count".
   #464a's table has five more countable anchors. Remove ALL SIX at once and see what is left.

Worlds
  A  the residual still predicts shame -> there is a piece of the breadth-type coordinate that
     NO countable quantity on this page touches, and that is why five naming attempts kept
     grasping at air. That sentence can go on the page.
  B  it collapses -> `c3-` is the negative of some weighted blend of those counts. Its name is
     plain, and the naming attempts failed because I kept looking for something richer.

⚠ TWO PRE-COMPUTATIONS, BEFORE THE OUTCOME FITS (#465a's lesson: I read a correlation as if it
   were a variance share -- 0.46 sounds large, 0.46^2 = 21% is what it can remove):
   (1) the six quantities' mutual correlations, since over-removal is the failure mode here;
   (2) their JOINT R-squared on `c3-` -- that, not any single r, is what the removal can take.
Both arms reported (#442a). Same model and mask as #509 so only the focal changes.
CONTROL : the raw arm must reproduce #509 (+0.0784 on shame).
CONTROL2: the residual must be orthogonal to all six.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R510 remove every countable")
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
_C1=sc(0).copy(); _C3=sc(2).copy()
PICKS=np.zeros(NN); _seen=np.zeros(NN,bool)
for _Mb,_ppl in MB:
    PICKS[_ppl]+=_Mb.sum(1); _seen[_ppl]=True
_PICKS=np.where(_seen,PICKS,np.nan)
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy(); _COVB=np.array(COVB,dtype=float).copy()
_NCAT=np.array(ncat,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
PH=num('pornhabit'); FET=num('totalfetishcategory')
COUNTS={'性行为计数':num('Totalsexacts'),'恋物类别数':FET,
        '性伴数分档':raw['sexcount'].astype(str).map(BANDS).values.astype(float),
        '块覆盖数':_COVB,'总勾选项数':_PICKS,'起始类别数':_NCAT}
SH=np.asarray(OUT['羞耻'],dtype=float); BE=np.asarray(OUT['能不能改'],dtype=float)
C1P=-_C1 if np.corrcoef(*[x[m&np.isfinite(_C1)&np.isfinite(FET)] for x in (_C1,FET)])[0,1]<0 else _C1
C3M=-_C3
MM = M & np.isfinite(AGE) & np.isfinite(PH) & np.isfinite(_EARLY) & np.isfinite(C1P) & np.isfinite(C3M)
for v in COUNTS.values(): MM = MM & np.isfinite(v)
n=int(MM.sum()); print(f"n = **{n:,}**")

KS=list(COUNTS)
Rm=np.array([[float(np.corrcoef(COUNTS[a][MM],COUNTS[b][MM])[0,1]) for b in KS] for a in KS])
iu=np.triu_indices(len(KS),1)
print(f"⚠ **前置 1 — 共线性**:六个量两两 |r| 中位 **{np.median(np.abs(Rm[iu])):.3f}** · "
      f"最大 **{np.max(np.abs(Rm[iu])):.3f}**")
Xc=np.column_stack([np.ones(n)]+[z(v,MM) for v in COUNTS.values()])
yv=z(C3M,MM); bb=np.linalg.lstsq(Xc,yv,rcond=None)[0]; res=yv-Xc@bb
R2=1-float(res.var()/yv.var())
print(f"⚠ **前置 2 — 六个量联合能解释 `c3⁻` 的 {R2:.1%}**(单看最强那个是 21.2%)-> "
      f"**剩下 {1-R2:.1%} 是任何计数都碰不到的上界**")
C3R=np.full(len(MM),np.nan); C3R[MM]=res
mx=max(abs(float(np.corrcoef(C3R[MM],v[MM])[0,1])) for v in COUNTS.values())
GATE.asserted("CONTROL2 the residual is orthogonal to all six", mx<1e-8,
              f"max |corr| = {mx:.2e}", kind="control")

def model(fv,y,over=None):
    Q={'f':fv,'`c1⁺` 卷入广度':C1P,'冷门程度 S':A,'常规也管用(−五题)':Bv,
       '色情使用量':PH,'起始年龄(晚)':_EARLY}
    src=Q if over is None else over
    cols=[np.ones(n), z(np.asarray(src['f'],dtype=float),MM)]+[
        z(np.asarray(v,dtype=float),MM) for k,v in src.items() if k!='f']+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])
ARMS={'原始 `c3⁻`':C3M,'**扣掉全部六个计数后**':C3R}
rows=[dict(arm=k,b_羞耻=model(v,SH),b_能不能改=model(v,BE)) for k,v in ARMS.items()]
NP_=400; nul=np.zeros((NP_,4))
for i in range(NP_):
    j=0
    for k,v in ARMS.items():
        pv=perm_in(np.asarray(v,dtype=float),MM,seed=26000+i)
        for y in (SH,BE): nul[i,j]=abs(model(pv,y)); j+=1
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows); T['sig_羞耻']=T.b_羞耻.abs()>thr; T['sig_能不能改']=T.b_能不能改.abs()>thr
T['keep']=T.b_羞耻.abs()/max(abs(T.b_羞耻.iloc[0]),1e-12)
show(T, HERE/'results/arms.csv', n=4, label="扣掉全部计数前后")
print(f"   **族内阈 = {thr:.5f}**")
GATE.asserted("CONTROL the raw arm reproduces #509", abs(T.b_羞耻.iloc[0]-0.0784)<0.01,
              f"raw shame = {T.b_羞耻.iloc[0]:+.4f} vs #509 +0.0784", kind="control")
r=T.iloc[1]; survives=bool(r.sig_羞耻)
print(f"\n扣掉全部六个之后:羞耻 **{r.b_羞耻:+.4f}**(保留 **{r.keep:.3f}**)· "
      f"能不能改 **{r.b_能不能改:+.4f}**")
GATE.asserted("KILL a piece survives that NO countable quantity touches", survives,
              f"residual shame {r.b_羞耻:+.4f}, sig={survives}")
verdict = "UNCOUNTABLE_PIECE_SURVIVES" if survives else "IT_IS_A_BLEND_OF_COUNTS"
print(f"\n判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,n=n,thr=thr,joint_R2=R2,
                collin_med=float(np.median(np.abs(Rm[iu]))),
                rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
