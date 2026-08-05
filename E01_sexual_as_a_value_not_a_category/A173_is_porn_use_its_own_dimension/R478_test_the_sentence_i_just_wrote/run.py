import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #433b closed with a reading I invented in the same sentence: "use -> shame -> behaviour,
   and the middle segment is the only one carrying information." That was never measured.
   #431c's lesson is that a smooth-reading chain is exactly what must be tested as an
   assertion. So: is there a mediation at all?

PRE-REGISTERED PREDICTION, written before running (#433's NEXT):
   I expect WORLD A -- no mediation -- because the second segment's uncontrolled coefficient
   is only z=1.02. And world A RETRACTS the sentence I just wrote, which is why it is worth
   running (frontier ss3).

Worlds
  A  no mediation      : the indirect effect a*b sits inside its own null -> retract.
  B  suppression       : controlling shame makes pornhabit's coefficient on "acted" LARGER,
                         i.e. the indirect path is NEGATIVE while the total is positive. That
                         is a different sentence, not the one I wrote.
  C  the chain as read : indirect positive and a real share of the total.

⚠ THE IDENTITY TRAP (#431/#440): in OLS, c = c' + a*b EXACTLY. So "the coefficient moved by
   a*b" is an ALGEBRAIC IDENTITY, not evidence. The measurement is whether a*b exceeds its
   OWN sampling null, and what its sign is.
Nulls: (1) permutation of pornhabit -- kills a, so indirect -> 0; (2) person-level bootstrap
   for the interval. Both reported.
CONTROL: the identity must hold numerically (c - c' - a*b ~ 0) -- if it does not, the model
   is mis-specified and nothing below means anything.
FRONTIER.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show

GATE=Gate("R478 mediation")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float); ACT=np.asarray(OUT['实践了多少'],dtype=float)
MM = M & np.isfinite(PH) & np.isfinite(AGE) & np.isfinite(SHAME) & np.isfinite(ACT)
print(f"n = **{int(MM.sum()):,}** · 中介 = 羞耻 · 结局 = 实践了多少 · focal = `pornhabit`")

BASE=lambda idx: np.column_stack([np.ones(idx.sum()), z(A,idx), z(Bv,idx),
                                  z(C3,idx), z(ncat,idx), z(AGE,idx)])
def paths(ph, idx):
    Xb=BASE(idx); p=z(ph,idx); sh=z(SHAME,idx); ac=z(ACT,idx)
    a  = np.linalg.lstsq(np.column_stack([Xb,p]),      sh, rcond=None)[0][-1]
    c  = np.linalg.lstsq(np.column_stack([Xb,p]),      ac, rcond=None)[0][-1]
    bb = np.linalg.lstsq(np.column_stack([Xb,p,sh]),   ac, rcond=None)[0]
    cp, b = float(bb[-2]), float(bb[-1])
    return float(a), float(b), float(c), cp

a,b,c,cp = paths(PH, MM)
ind = a*b
print(f"\na(用量->羞耻) **{a:+.4f}** · b(羞耻->实践,控用量) **{b:+.4f}** · "
      f"c(总) **{c:+.4f}** · c'(控羞耻) **{cp:+.4f}**")
print(f"间接 a*b = **{ind:+.5f}** · 直接 c' = **{cp:+.4f}** · 恒等式残差 c−c'−ab = **{c-cp-ind:+.2e}**")

GATE.asserted("CONTROL the OLS identity c = c' + a*b holds numerically",
              abs(c-cp-ind)<1e-9, f"residual = {c-cp-ind:.2e}", kind="control")

NP_=400
nul=np.array([ (lambda t: t[0]*t[1])(paths(perm_in(PH,MM,seed=9100+i), MM)) for i in range(NP_)])
rng=np.random.default_rng(7)
idxall=np.flatnonzero(MM)
bs=[]
for i in range(400):
    take=rng.choice(idxall,len(idxall),replace=True)
    mm=np.zeros(len(MM),bool); mm[np.unique(take)]=True
    t=paths(PH,mm); bs.append(t[0]*t[1])
bs=np.array(bs)
lo,hi=np.percentile(bs,[2.5,97.5])
thr=float(np.percentile(np.abs(nul),95))
T=pd.DataFrame([dict(quantity='a 用量->羞耻',value=a),dict(quantity="b 羞耻->实践",value=b),
                dict(quantity='c 总',value=c),dict(quantity="c' 直接",value=cp),
                dict(quantity='a*b 间接',value=ind),dict(quantity='间接零 95 分位|.|',value=thr),
                dict(quantity='自助 2.5%',value=lo),dict(quantity='自助 97.5%',value=hi)])
show(T,HERE/'results/mediation.csv',n=8,label="中介")

sig = abs(ind)>thr and (lo>0)==(hi>0)
share = ind/c if c!=0 else np.nan
print(f"\n间接 |{ind:+.5f}| vs 零阈 {thr:.5f} -> {'**越阈**' if abs(ind)>thr else '**落在零里**'} · "
      f"自助区间 [{lo:+.5f}, {hi:+.5f}] {'不含 0' if (lo>0)==(hi>0) else '**含 0**'}")
print(f"间接 / 总 = **{share:+.1%}** · 而 c' **{cp:+.4f}** {'>' if abs(cp)>abs(c) else '<'} c **{c:+.4f}**"
      f" -> {'**抑制(世界 B)**' if abs(cp)>abs(c) and ind*c<0 else ''}")

GATE.asserted("CONTROL the permutation null kills the indirect path",
              abs(nul.mean())<abs(ind) or abs(nul.mean())<1e-3,
              f"null mean = {nul.mean():+.5f}", kind="control")
GATE.asserted("KILL the mediation I wrote exists (indirect POSITIVE and beyond its null)",
              bool(sig and ind>0), f"indirect {ind:+.5f}, sig={sig}")
verdict = ("CHAIN_AS_READ" if (sig and ind>0) else
           ("SUPPRESSION" if (sig and ind<0) else "NO_MEDIATION"))
print(f"\n判决 = {verdict}   (预注册预测 = NO_MEDIATION)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),a=a,b=b,c=c,cp=cp,indirect=ind,
               null_thr=thr,boot_lo=lo,boot_hi=hi,share=share,prediction="NO_MEDIATION"),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
