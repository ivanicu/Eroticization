import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: this page carries two strong arcs that have never been joined. One: things that arrived
   EARLIER go with MORE shame (replicated across calibers, ~2x in non-biomale respondents).
   Two (new, #437-#442): breadth type predicts feeling UNABLE TO CHANGE what arouses you, and
   that is not explained by shame, by experience, by involvement, or by any single block.
   So: **does what arrived early feel more fixed?**

Worlds
  A  early = fixed : later mean onset goes with feeling MORE able to change (positive
     coefficient on onset age). Then "I can't change this" is partly a fact about WHEN it
     arrived, and the two arcs join.
  B  no relation   : onset says nothing about the belief once current age is held. Then the
     sense of being unable to change is not about timing, and the arcs stay separate -- which
     is itself worth knowing, because timing is the most obvious folk explanation.
  C  the reverse   : later onset goes with feeling LESS able to change.

PRE-REGISTERED PREDICTION: **A**. ⚠ Five of my seven predictions this session were wrong, so
this is written to be killed.

focal    = `EARLY`, the person's mean onset age. Direction fixed BY CONSTRUCTION (higher =
           later), so no anchoring is needed.
⚠ CURRENT AGE IS A HARD CONTROL: onset age and current age are mechanically linked, and the
   page's own onset results all carry it.
controls = R449's set (S, -(five-item), c3-, category count) + current age.
POSITIVE CONTROL, and it is the design's spine: the SHAME cell must reproduce the page's own
   established result (earlier onset -> more shame). If it does not, the instrument is not
   measuring what the page already knows, and the belief cell is unreadable.
NULL     = permute the focal within the analysis mask -- a single-column statistic, so
   `perm_in` is the right family here (#426c), and the world being excluded is "this quantity
   is unrelated to the outcome", which is exactly what the person-permutation excludes (#450b).
MULTIPLICITY: four outcomes -> the null of the maximum (#440b).
⚠ DEVIATION, named: #450's NEXT was a page-level null-family audit -- pure Closure, and the
   seventh scanning round after six over-indictments. The standing instruction is psychology
   first, and #436's own tally found method rounds outnumbering object rounds. The audit is
   deferred, not dropped.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R495 does what came early feel fixed")
# ⚠ #451a: R449 cuts R416 before EARLY is defined, so EARLY is spliced from R416 FIRST and
# PRIVATISED before R449 runs (#449a's structural rule -- a cross-round exec imports a
# namespace, so anything that must survive it cannot share a name with the other round).
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
print(f"`EARLY` 取自 `R416`,已私有化(n 有值 = {int(np.isfinite(_EARLY).sum()):,})")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
E=_EARLY
MM = M & np.isfinite(AGE) & np.isfinite(E)
print(f"n = **{int(MM.sum()):,}** · `EARLY` 值域 [{np.nanmin(E[MM]):.1f}, {np.nanmax(E[MM]):.1f}] · "
      f"中位 **{np.nanmedian(E[MM]):.1f}** · **方向由构造固定(大 = 晚),不需要锚定**")
print(f"⚠ `EARLY` 与当前年龄相关 **{np.corrcoef(E[MM],AGE[MM])[0,1]:+.4f}** -> 当前年龄是硬控制")

def coef(v,y,idx):
    X=np.column_stack([np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
                       z(C3,idx), z(ncat,idx), z(AGE,idx)])
    return float(np.linalg.lstsq(X,z(np.asarray(y,dtype=float),idx),rcond=None)[0][1])

rows=[dict(outcome=nm, b=coef(E,y,MM)) for nm,y in OUT.items()]
NP_=400
nul=np.zeros((NP_,len(OUT)))
for i in range(NP_):
    pv=perm_in(E,MM,seed=18000+i)
    for j,(nm,y) in enumerate(OUT.items()): nul[i,j]=abs(coef(pv,y,MM))
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows); T['sig']=T.b.abs()>thr
show(T, HERE/'results/onset.csv', n=6, label="起始年龄 × 四结局")
print(f"   **多重性阈(四个里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

sh=T[T.outcome=='羞耻'].iloc[0]
# ⚠ #451b: the control's SIGN as first written was wrong -- I asserted b>0 while its own label
# says "later onset -> LESS shame", which is b<0. Fifth gate this session that tested something
# other than its own sentence (#433a #439d #440b #444a). Both versions are kept: the mis-written
# one FAILS, the sentence it meant to test PASSES, and the difference is a typo in a direction,
# which is exactly the class of error #392 cost four rounds.
GATE.asserted("POSITIVE CONTROL as first written (b>0) -- MIS-SPECIFIED, kept on the record",
              bool(sh.sig and sh.b>0), f"asserted b>0; observed {sh.b:+.4f}", kind="control")
GATE.asserted("POSITIVE CONTROL as meant: later onset -> LESS shame (b<0), reproducing the page",
              bool(sh.sig and sh.b<0),
              f"onset -> shame b = {sh.b:+.4f} (page: earlier onset -> more shame)", kind="control")
bel=T[T.outcome=='能不能改'].iloc[0]
print(f"\n★ **关键一格 —— 能不能改**:b = **{bel.b:+.4f}** · 阈 {thr:.5f} -> "
      f"{'**越阈**' if bool(bel.sig) else '**未越阈**'}")
GATE.asserted("KILL later onset goes with feeling MORE able to change (world A)",
              bool(bel.sig and bel.b>0), f"b = {bel.b:+.4f}, sig = {bool(bel.sig)}")
verdict = ("EARLY_FEELS_FIXED" if (bel.sig and bel.b>0)
           else ("REVERSED" if (bel.sig and bel.b<0) else "TIMING_SAYS_NOTHING"))
print(f"\n判决 = {verdict}   (预注册预测 = A / EARLY_FEELS_FIXED)")
# MDE so a null is readable
mde=1.96*float(nul[:, [i for i,(nm,_) in enumerate(OUT.items()) if nm=='能不能改'][0]].std())
print(f"⚠ 该格的 MDE ≈ **{mde:.4f}**(观测 {abs(bel.b):.4f})")
_json.dump(dict(verdict=verdict,n=int(MM.sum()),thr=thr,mde=mde,
                rows=T.to_dict('records'),prediction="A"),
           open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
