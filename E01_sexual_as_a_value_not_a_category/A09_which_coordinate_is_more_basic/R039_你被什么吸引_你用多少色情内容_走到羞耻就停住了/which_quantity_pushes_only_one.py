import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: shame and the sense of being able to change correlate **-0.0022** with each other, yet both
   are moved by when a person's tastes arrived. What else moves each -- and is there a quantity
   that moves only one? That would be the fork between the two paths.

⚠ THE ARITHMETIC, WRITTEN BEFORE RUNNING (#457's NEXT demanded this, because #456 was a round
   whose emptiness was only noticed afterwards):
   if the two outcomes are near-orthogonal, then summed over any set of predictors that
   accounts for them, the products r(X,Y1)*r(X,Y2) must very nearly CANCEL. **So it is
   arithmetically near-impossible for every predictor to push both the same way.** Sign
   disagreement somewhere is not a discovery -- it is forced.
   ⇒ world "A: everything pushes both" is essentially excluded before the round starts, and
   the round is therefore NOT "does something push only one" (something must) but:
   **WHICH ones agree, which disagree, and does the cancellation actually show up in the data.**
   That reframing is the pre-registration.

Worlds
  B1 a clean fork : one or two quantities load on exactly one outcome, the rest on both ->
     the fork can be named.
  B2 diffuse cancellation : the agreement pattern is spread thin with no quantity clearly
     single-sided -> the two paths separate through many small opposing contributions, and no
     single fork exists to name.

Quantities (ALL already on the page, none built for this question, each with its direction
   established elsewhere): S · c3- · five-item score · category count · porn use · onset age.
Both raw correlations AND partial coefficients (the page's control set minus the focal) are
   reported, because a fork that exists only in one of the two is a different claim.
MULTIPLICITY: 6 quantities x 2 outcomes -> the null of the maximum (#440b).
CONTROL : the cancellation sum must actually be near zero, or my pre-registered arithmetic is
   wrong and nothing below follows from it.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R502 which quantity pushes only one")
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
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY) & np.isfinite(PH)
n=int(MM.sum())
Q={'冷门程度 S':A, '广度型 c3⁻':C3, '常规也管用(−五题)':Bv,
   '答题类别数':ncat, '色情使用量':PH, '起始年龄(晚)':_EARLY}
YS={'羞耻':SH, '能不能改':BE}
r_out=float(np.corrcoef(SH[MM],BE[MM])[0,1])
print(f"n = **{n:,}** · 两个结局彼此相关 **{r_out:+.4f}**")

def rr(u,v): return float(np.corrcoef(np.asarray(u)[MM],np.asarray(v)[MM])[0,1])
prods=[rr(v,SH)*rr(v,BE) for v in Q.values()]
S_=float(np.sum(prods))
print(f"⚠ **事前算术的检查**:Σ r(X,羞耻)·r(X,能不能改) = **{S_:+.5f}** "
      f"(逐项 {[round(p,4) for p in prods]})")
GATE.asserted("CONTROL the cancellation my pre-registration assumed actually happens",
              abs(S_)<0.02, f"sum of products = {S_:+.5f}", kind="control")

def coef(v,y,extra_drop):
    cols=[np.ones(n), z(np.asarray(v,dtype=float),MM)]+[z(np.asarray(u,dtype=float),MM)
          for k,u in Q.items() if k!=extra_drop]+[z(AGE,MM)]
    return float(np.linalg.lstsq(np.column_stack(cols),
                                 z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])
rows=[]
for k,v in Q.items():
    d=dict(quantity=k)
    for yn,y in YS.items():
        d[f'r_{yn}']=rr(v,y); d[f'b_{yn}']=coef(v,y,k)
    rows.append(d)
T=pd.DataFrame(rows)
NP_=400; nul=np.zeros((NP_,len(Q)*len(YS)))
for i in range(NP_):
    j=0
    for k,v in Q.items():
        pv=perm_in(np.asarray(v,dtype=float),MM,seed=21000+i)
        for yn,y in YS.items(): nul[i,j]=abs(coef(pv,y,k)); j+=1
thr=float(np.percentile(nul.max(1),95))
for yn in YS: T[f'sig_{yn}']=T[f'b_{yn}'].abs()>thr
T['pushes']=[('两个' if r[f'sig_羞耻'] and r[f'sig_能不能改']
              else ('只推羞耻' if r[f'sig_羞耻'] else
                    ('只推能不能改' if r[f'sig_能不能改'] else '都不推')))
             for _,r in T.iterrows()]
show(T[['quantity','b_羞耻','b_能不能改','sig_羞耻','sig_能不能改','pushes']],
     HERE/'results/profiles.csv', n=6, label="六个量 × 两结局")
print(f"   **族内阈(12 格里最大 |b| 的零分布 95 分位)= {thr:.5f}**")
single=T[T.pushes.str.startswith('只推')]
print(f"\n**只推其中一个的量 = {len(single)}**:{list(single.quantity)}")
GATE.asserted("KILL there is a clean fork (at least one quantity pushes exactly one outcome)",
              len(single)>0, f"single-sided quantities = {list(single.quantity)}")
verdict = "CLEAN_FORK" if len(single)>0 else "DIFFUSE"
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,thr=thr,r_outcomes=r_out,sum_products=S_,
               rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
