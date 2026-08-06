import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #442b gave `c1` a direction after the one-anchor rule was repaired -- "more fetish
   categories, less sense that this can be changed" -- and put it on the page. **That sentence
   has never been tested on its own.** #458b showed `c3-` pushes only shame. Where does `c1`
   land?

Worlds
  A  c1 pushes only changeability -> the two leading components fall on the TWO DIFFERENT
     paths, and "breadth type" versus "fetish breadth" is the concrete content of the fork
     #458b found. The page could then name both sides of it for the first time.
  B  c1 pushes both, or neither -> the fork is not carried by the components, and #458b's
     reading narrows.

⚠ `c1`'s direction is taken from #442b's anchor (fetish-category count, r = -0.131) and is NOT
   re-anchored here; it is oriented so that higher = MORE fetish categories, which is the frame
   the page's sentence is written in.
⚠ The block mask is smaller than #458b's sample -- reported before the result, because a null
   on a smaller n is not the same null (#413b).
Model: `c1+` joined to #458b's six quantities + current age. Both outcomes. Null of the maximum.
CONTROL : `c3-` must still push only shame in this model -- if it does not, the model changed
   more than intended and `c1`'s row is not comparable to #458b's.
CONTROL2: `c1+`'s anchor correlation must reproduce #442b's -0.131 in magnitude.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R507 where does c1 land")
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
_C1=sc(0).copy()

raw0=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
FET0=pd.to_numeric(raw0['totalfetishcategory'],errors='coerce').values.astype(float)
gg=m&np.isfinite(_C1)&np.isfinite(FET0)
r_raw=float(np.corrcoef(_C1[gg],FET0[gg])[0,1])
C1P = -_C1 if r_raw<0 else _C1          # 朝「恋物类别更多」(#442b 的框架),不重新锚定
r_anch=float(np.corrcoef(C1P[gg],FET0[gg])[0,1])
print(f"`c1` 原始与恋物类别数 **{r_raw:+.4f}**(`#442b` 报 −0.1307)· "
      f"定向后 **{r_anch:+.4f}**(高 = 恋物类别更多)")
GATE.asserted("CONTROL2 c1's anchor reproduces #442b in magnitude",
              abs(abs(r_raw)-0.1307)<0.02, f"|r| = {abs(r_raw):.4f} vs 0.1307", kind="control")
_C1P=C1P.copy()

_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float); BE=np.asarray(OUT['能不能改'],dtype=float)
Q={'**`c1⁺` 恋物广度**':_C1P, '冷门程度 S':A, '广度型 c3⁻':C3, '常规也管用(−五题)':Bv,
   '答题类别数':ncat, '色情使用量':PH, '起始年龄(晚)':_EARLY}
MM = M & np.isfinite(AGE) & np.isfinite(PH) & np.isfinite(_EARLY) & np.isfinite(_C1P)
n=int(MM.sum())
print(f"n = **{n:,}**(块掩码,`#458b` 是 6,715 -> **先报,零不是同一个零**)")

def coef(key,y,over=None):
    src=Q if over is None else over
    cols=[np.ones(n), z(np.asarray(src[key],dtype=float),MM)]+[
        z(np.asarray(v,dtype=float),MM) for k,v in src.items() if k!=key]+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])
rows=[dict(quantity=k, b_羞耻=coef(k,SH), b_能不能改=coef(k,BE)) for k in Q]
NP_=400; nul=np.zeros((NP_,len(Q)*2))
for i in range(NP_):
    j=0
    for k in Q:
        pq={kk:(perm_in(np.asarray(v,dtype=float),MM,seed=24000+i) if kk==k else v)
            for kk,v in Q.items()}
        for y in (SH,BE): nul[i,j]=abs(coef(k,y,over=pq)); j+=1
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows)
T['sig_羞耻']=T.b_羞耻.abs()>thr; T['sig_能不能改']=T.b_能不能改.abs()>thr
T['pushes']=[('两个' if r.sig_羞耻 and r.sig_能不能改 else
              ('只推羞耻' if r.sig_羞耻 else ('只推能不能改' if r.sig_能不能改 else '都不推')))
             for _,r in T.iterrows()]
show(T[['quantity','b_羞耻','b_能不能改','pushes']], HERE/'results/c1_landing.csv', n=8,
     label="七个量 × 两结局")
print(f"   **族内阈(14 格里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

c3r=T[T.quantity=='广度型 c3⁻'].iloc[0]
GATE.asserted("CONTROL c3- still pushes only shame in this model",
              c3r.pushes=='只推羞耻', f"c3- pushes: {c3r.pushes}", kind="control")
c1r=T[T.quantity=='**`c1⁺` 恋物广度**'].iloc[0]
print(f"\n★ **`c1⁺`**:羞耻 **{c1r.b_羞耻:+.4f}** · 能不能改 **{c1r.b_能不能改:+.4f}** -> "
      f"**{c1r.pushes}**")
onlyBE = c1r.pushes=='只推能不能改'
GATE.asserted("KILL c1 lands on the OTHER path (pushes only changeability)",
              onlyBE, f"c1+ pushes {c1r.pushes}")
verdict = ("TWO_COMPONENTS_TWO_PATHS" if onlyBE else
           ("BOTH" if c1r.pushes=='两个' else ("SHAME_ONLY" if c1r.pushes=='只推羞耻' else "NEITHER")))
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,thr=thr,r_anchor=r_raw,
               rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
