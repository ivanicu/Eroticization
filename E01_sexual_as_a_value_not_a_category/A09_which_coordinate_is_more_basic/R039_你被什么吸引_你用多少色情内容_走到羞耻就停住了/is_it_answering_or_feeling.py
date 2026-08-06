import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #458c(1) named the weak link: "how many categories a person answered" pushes only the sense
   of changeability, and that could be about ANSWERING rather than about feeling. Is it?

⚠ TWO CORRECTIONS FROM READING, BEFORE ANY COMPUTATION (#452a's rule paid off again):
  (1) my own NEXT said to reuse "`#100`'s caliper matching". The project's record says the
      caliper DID NOT WORK (`#226a`: medians 49/62/69, still increasing) and what worked was
      `#234a`'s DECILE-BINNED RESAMPLING to a common distribution. The signpost I wrote cited
      the wrong entry AND the broken tool.
  (2) more decisive: **matching on the pick count would remove the focal itself.** Matching is
      the wrong instrument for this question, whichever version. Written down, not swapped
      silently.

The separator that does fit: `#421` already enumerated FIVE measures of how much a person
engaged -- onset-category count, fetish-category count, endorsed sex acts, block coverage,
total picks -- and established they are NOT interchangeable (six of ten pairs reach 0.6; the
sex-act count is the outlier at 0.26-0.51). So:

Worlds
  A  it is answering volume : the other volume measures behave like `ncat` -- pushing the
     changeability outcome and not shame. Then `#458b`'s single-sidedness is about the survey,
     and the page must say so.
  B  `ncat` is not just volume : the volume measures disagree with each other. Then "answered
     more categories" is carrying something the others do not, and calling it an answering
     artifact was wrong.

MULTIPLICITY: measures x 2 outcomes -> the null of the maximum (#440b).
CONTROL : `ncat`'s own two coefficients must reproduce #458b (+0.0102 / -0.0354) when it is run
   in the same model, or the arm is not comparable.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R503 answering or feeling")
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
_COVB=np.array(COVB,dtype=float).copy()
# PICKS is not in R416's spliced portion; R465 builds it from the block matrices. Reused
# verbatim rather than re-invented (P4), then PRIVATISED before the next splice (#449a).
PICKS=np.zeros(NN); _seen=np.zeros(NN,bool)
for _Mb,_ppl in MB:
    PICKS[_ppl]+=_Mb.sum(1); _seen[_ppl]=True
PICKS=np.where(_seen,PICKS,np.nan)
_PICKS=np.array(PICKS,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values.astype(float)
FET=pd.to_numeric(raw['totalfetishcategory'],errors='coerce').values.astype(float)
SXA=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float); BE=np.asarray(OUT['能不能改'],dtype=float)
VOL={'起始类别数 ncat':np.asarray(ncat,dtype=float), '块覆盖数 COVB':_COVB,
     '总勾选项数 PICKS':_PICKS, '恋物类别数':FET, '认可的性行为数':SXA}
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY) & np.isfinite(PH)
for v in VOL.values(): MM = MM & np.isfinite(v)
n=int(MM.sum()); print(f"n = **{n:,}**(五个广度量都要有值)")
print("五个量彼此的相关(`#421` 说它们不可互换):")
KS=list(VOL)
Rm=np.array([[float(np.corrcoef(VOL[a][MM],VOL[b][MM])[0,1]) for b in KS] for a in KS])
for i,a in enumerate(KS):
    print(f"   {a:<16} " + " ".join(f"{Rm[i,j]:+.2f}" for j in range(len(KS))))

# 控制集 = #458b 的其余五个量(不含被检的那个 volume 量)+ 年龄
OTHERS={'冷门程度 S':A,'广度型 c3⁻':C3,'常规也管用(−五题)':Bv,'色情使用量':PH,'起始年龄':_EARLY}
def coef(v,y):
    cols=[np.ones(n), z(np.asarray(v,dtype=float),MM)]+[z(np.asarray(u,dtype=float),MM)
          for u in OTHERS.values()]+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])
rows=[dict(measure=k, b_羞耻=coef(v,SH), b_能不能改=coef(v,BE)) for k,v in VOL.items()]
NP_=400; nul=np.zeros((NP_,len(VOL)*2))
for i in range(NP_):
    j=0
    for k,v in VOL.items():
        pv=perm_in(np.asarray(v,dtype=float),MM,seed=22000+i)
        for y in (SH,BE): nul[i,j]=abs(coef(pv,y)); j+=1
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows)
T['sig_羞耻']=T.b_羞耻.abs()>thr; T['sig_能不能改']=T.b_能不能改.abs()>thr
T['pushes']=[('两个' if r.sig_羞耻 and r.sig_能不能改 else
              ('只推羞耻' if r.sig_羞耻 else ('只推能不能改' if r.sig_能不能改 else '都不推')))
             for _,r in T.iterrows()]
show(T, HERE/'results/volume_measures.csv', n=6, label="五个广度量")
print(f"   **族内阈(10 格里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

nc=T[T.measure=='起始类别数 ncat'].iloc[0]
GATE.asserted("CONTROL ncat reproduces #458b in this model",
              abs(nc.b_羞耻-0.0102)<0.02 and abs(nc.b_能不能改-(-0.0354))<0.02,
              f"ncat = {nc.b_羞耻:+.4f} / {nc.b_能不能改:+.4f} vs #458b +0.0102 / -0.0354",
              kind="control")
pat=list(T.pushes)
same_as_ncat = sum(1 for p in pat if p=='只推能不能改')
print(f"\n**与 `ncat` 同型(只推能不能改)的量 = {same_as_ncat}/5** · 全部模式:{pat}")
GATE.asserted("KILL it is answering volume (the other measures behave like ncat)",
              same_as_ncat>=4, f"{same_as_ncat}/5 measures push only changeability")
verdict = "ANSWERING_VOLUME" if same_as_ncat>=4 else "NCAT_IS_NOT_JUST_VOLUME"
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,thr=thr,same=same_as_ncat,
               rows=T.to_dict('records'),corr=Rm.tolist(),keys=KS),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
