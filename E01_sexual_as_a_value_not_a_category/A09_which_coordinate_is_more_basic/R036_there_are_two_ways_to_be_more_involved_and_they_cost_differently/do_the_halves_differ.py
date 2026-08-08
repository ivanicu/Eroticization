import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #442b found c2 and c3 are BOTH anchored by the same count -- endorsed sex acts -- but with
   OPPOSITE signs (-0.376 and +0.460). So the leading residual structure cuts "being involved"
   into two halves. Turned to face the same way, do the two halves have the same consequences?

Worlds
  A  one thing twice : the two profiles across the four outcomes agree -> they are two readings
     of involvement, and #441b's "the whole structure" narrows again.
  B  two kinds       : the profiles differ -> there are two ways of being more involved and
     they land differently on people. That is the sentence, and it is the first concrete shape
     of "what involvement does not contain" (#442d(3)).

Orientation (#442b): c2 is FLIPPED so that both point toward MORE involvement --
   c2+ = -c2 (correlates +0.376 with the count), c3 as-is (+0.460).
   ⚠ Signs are stated in the ANCHORED frame, never in the frame numpy happened to return.
Outcomes: the four R449 carries. Controls: R449's set, which never controls the outcome itself.
MULTIPLICITY: 2 components x 4 outcomes -> the bar is the NULL OF THE MAXIMUM over all 8
   (#440b, now standing practice).
PROFILE COMPARISON: correlation between the two 4-vectors of coefficients, with a person-level
   bootstrap interval. Four points is few, so the interval is the result, not the point estimate.
⚠ BARS ARE SET BY A POSITIVE CONTROL, NOT BY ME (#442b): the anchor rule used here is
   "clears 0.10 on either count", which is exactly the rule that passes both known-direction
   variables on this mask.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R487 two ways of being involved")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
keep=list(range(NB)); Ra,Rb=prof_(A,keep),prof_(B,keep)
C=np.zeros((NB,NB))
for i in keep:
    for j in keep:
        mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
Rr=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
Rr=np.where(np.isfinite(Ra)|np.isfinite(Rb),Rr,np.nan)
Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
def score(k):
    num=(V[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
c2_raw, c3_raw = score(1), score(2)

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
SXA=pd.to_numeric(raw['Totalsexacts'],errors='coerce').values.astype(float)
# ⚠ #443a: `anc` collides -- R449's splice rebinds it to an ARRAY. Same class as #427e, and the
# failure is a TypeError only because the name is later CALLED; had it been read, it would have
# been silent. Cross-round splicing imports a namespace, not just numbers.
def anchor_r(v):
    g=m&np.isfinite(v)&np.isfinite(SXA); return float(np.corrcoef(v[g],SXA[g])[0,1])
c2 = -c2_raw if anchor_r(c2_raw)<0 else c2_raw
c3 = -c3_raw if anchor_r(c3_raw)<0 else c3_raw
ANC2, ANC3 = anchor_r(c2), anchor_r(c3)      # captured BEFORE the R449 splice runs
print(f"锚定后(两者都朝「更卷入」):c2⁺ vs 计数 **{ANC2:+.4f}** · c3 vs 计数 **{ANC3:+.4f}**")
GATE.asserted("CONTROL both components now face the same way",
              ANC2>0 and ANC3>0, f"{ANC2:+.4f} / {ANC3:+.4f}", kind="control")

_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
MM = M & np.isfinite(AGE) & np.isfinite(c2) & np.isfinite(c3)
print(f"n = **{int(MM.sum()):,}** · 四个结局 = {list(OUT)}")

def coef(v, y, idx):
    X=np.column_stack([np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
                       z(ncat,idx), z(AGE,idx)])
    return float(np.linalg.lstsq(X,z(np.asarray(y,dtype=float),idx),rcond=None)[0][1])

COMP={'c2⁺ 卷入之一':c2, 'c3 卷入之二':c3}
rows=[]
for cn,cv in COMP.items():
    for on,y in OUT.items():
        rows.append(dict(comp=cn, outcome=on, b=coef(cv,y,MM)))
T=pd.DataFrame(rows)

NP_=400
nul=np.zeros((NP_,len(COMP)*len(OUT)))
for i in range(NP_):
    k=0
    for cn,cv in COMP.items():
        pv=perm_in(cv,MM,seed=15000+i)
        for on,y in OUT.items():
            nul[i,k]=abs(coef(pv,y,MM)); k+=1
thr=float(np.percentile(nul.max(1),95))
T['sig']=T.b.abs()>thr
show(T, HERE/'results/profiles.csv', n=10, label="两半 × 四结局")
print(f"   **多重性阈(8 个里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

p2=T[T.comp=='c2⁺ 卷入之一'].set_index('outcome').b
p3=T[T.comp=='c3 卷入之二'].set_index('outcome').b.reindex(p2.index)
r_prof=float(np.corrcoef(p2,p3)[0,1])
rng=np.random.default_rng(53); idxall=np.flatnonzero(MM); bs=[]
for _ in range(300):
    take=rng.choice(idxall,len(idxall),replace=True)
    mm2=np.zeros(len(MM),bool); mm2[np.unique(take)]=True
    a2=[coef(c2,y,mm2) for y in OUT.values()]; a3=[coef(c3,y,mm2) for y in OUT.values()]
    bs.append(float(np.corrcoef(a2,a3)[0,1]))
bs=np.array(bs); plo,phi=np.percentile(bs,[2.5,97.5])
print(f"\n两条剖面的相关 = **{r_prof:+.4f}** · 自助区间 **[{plo:+.4f}, {phi:+.4f}]**")
print(f"   (⚠ 只有 4 个点 -> **区间才是结果,点估不是**)")

same = plo>0.5
GATE.asserted("KILL the two halves are one thing (profiles agree)", same,
              f"profile r = {r_prof:+.4f}, boot [{plo:+.4f}, {phi:+.4f}]")
verdict = "ONE_THING_TWICE" if same else "TWO_KINDS"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,n=int(MM.sum()),thr=thr,profile_r=r_prof,
               boot=[plo,phi],anchors=[ANC2,ANC3],rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
