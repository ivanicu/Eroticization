import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #455e(4) named the weak spot: "one antecedent" is currently just another way of saying "two
   gradients agree". That is falsifiable. If onset really reaches both outcomes through one
   thing, then holding one outcome should flatten the other's gradient.

Worlds
  A  one antecedent : both gradients collapse when the other outcome is controlled -> the page
     can merge its two arcs into one.
  B  two paths from one start : at least one gradient survives -> onset has a route to that
     outcome that does not pass through the other, and "one antecedent" must be withdrawn in
     favour of "two paths sharing a starting point".

PRE-REGISTERED PREDICTION: **B**, because #437b already showed the changeability cell does not
move when shame enters the model -- same control, different focal (`c3-` there, onset here).
⚠ Five of my eight predictions this session were wrong. Written to be killed.

Everything identical to #454/#455: the page's own four strata, the same control set, the same
monotonicity statistic. The ONLY change is which outcome joins the controls, and BOTH ARMS are
reported (#442a: a control absorbs real signal as well as confound).
⚠ DIRECTIONS ARE PRE-SPECIFIED and not re-chosen after seeing the data (#455c): changeability
   rises with later onset (+1), shame falls with later onset (-1). One-sided, because a
   two-sided ordering test cannot fire at four strata at all.
CONTROL : the uncontrolled arms must reproduce #455 exactly -- nothing changed but the control.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R500 one antecedent or two paths")
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
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY)
n=int(MM.sum()); q=np.quantile(_EARLY[MM],[0.25,0.5,0.75]); sub=np.digitize(_EARLY[MM],q)
BASE=[np.ones(n), z(A,MM), z(Bv,MM), z(C3,MM), z(ncat,MM), z(AGE,MM)]
PRED={'能不能改':+1, '羞耻':-1}                 # 事先给定,不在看数之后改(#455c)
OUTS={'能不能改':OUT['能不能改'], '羞耻':OUT['羞耻']}

def prof(y, extra=None, s=None):
    X=np.column_stack(BASE+([z(np.asarray(extra,dtype=float),MM)] if extra is not None else []))
    yy=z(np.asarray(y,dtype=float),MM)
    r=yy-X@np.linalg.lstsq(X,yy,rcond=None)[0]
    ss=sub if s is None else s
    return np.array([r[ss==k].mean() for k in range(4)])

rg=np.random.default_rng(101); rows=[]
for k,y in OUTS.items():
    other=[v for kk,v in OUTS.items() if kk!=k][0]
    for lab,ex in (('不控另一结局',None), ('**控另一结局**',other)):
        p=prof(y,extra=ex); rho=float(spearmanr(np.arange(4),p).statistic)
        nul=np.array([float(spearmanr(np.arange(4),prof(y,extra=ex,s=pm)).statistic)
                      for pm in (rg.permutation(sub) for _ in range(3000))])
        one=float(np.mean(np.sign(PRED[k])*nul >= np.sign(PRED[k])*rho))
        rows.append(dict(outcome=k, arm=lab, s1=p[0], s4=p[3], steep=float(p[3]-p[0]),
                         rho=rho, p_one_sided=one,
                         mono=bool(abs(rho)>=0.8 and np.sign(rho)==PRED[k] and one<0.05)))
T=pd.DataFrame(rows)
show(T, HERE/'results/arms.csv', n=6, label="两个结局 × 两臂")

c=T[(T.outcome=='能不能改')&(T.arm=='不控另一结局')].iloc[0]
GATE.asserted("CONTROL the uncontrolled arm reproduces #455 exactly",
              abs(c.steep-0.127643)<1e-4 and abs(c.rho-1.0)<1e-9,
              f"steepness {c.steep:+.6f} vs #455 +0.127643; rho {c.rho:+.3f}", kind="control")
sh=T[(T.outcome=='羞耻')&(T.arm=='不控另一结局')].iloc[0]
GATE.asserted("CONTROL2 the shame arm reproduces #455 too",
              abs(sh.steep-(-0.146519))<1e-4, f"steepness {sh.steep:+.6f} vs #455 -0.146519",
              kind="control")

print("\n保留比(控另一结局后的陡度 ÷ 不控时的陡度):")
keep={}
for k in OUTS:
    a=T[(T.outcome==k)&(T.arm=='不控另一结局')].steep.iloc[0]
    b=T[(T.outcome==k)&(T.arm=='**控另一结局**')].steep.iloc[0]
    keep[k]=abs(b)/max(abs(a),1e-12)
    print(f"   {k:<8} {a:+.4f} -> {b:+.4f}  **{keep[k]:.3f}**  "
          f"{'**仍单调越阈**' if bool(T[(T.outcome==k)&(T.arm=='**控另一结局**')].mono.iloc[0]) else '**塌了**'}")
survivors=[k for k in OUTS if bool(T[(T.outcome==k)&(T.arm=='**控另一结局**')].mono.iloc[0])]
GATE.asserted("KILL one antecedent (BOTH gradients collapse when the other outcome is held)",
              len(survivors)==0, f"survivors = {survivors}")
verdict = "ONE_ANTECEDENT" if len(survivors)==0 else "TWO_PATHS_ONE_START"
print(f"\n判决 = **{verdict}**  (预注册预测 = B / TWO_PATHS_ONE_START)")
json.dump(dict(verdict=verdict,n=n,keep=keep,survivors=survivors,
               rows=T.to_dict('records'),prediction="B"),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())

# ---------------------------------------------------------------- #456d(发布前追加)
# 近乎零的中介,两个方向都是 —— 这需要一个机制,而它可能是页面已知的:
# 若两个结局本身几乎不相关,那么互相控制**本来就**动不了对方。
_sh=np.asarray(OUT['羞耻'],dtype=float); _be=np.asarray(OUT['能不能改'],dtype=float)
_g=MM&np.isfinite(_sh)&np.isfinite(_be)
r_out=float(np.corrcoef(_sh[_g],_be[_g])[0,1])
print(f"\n**两个结局本身的相关 = {r_out:+.4f}**(n={int(_g.sum()):,})")
print(f"⇒ 近乎零的中介**不是意外** —— 互相控制动不了对方,因为两者本来就几乎不相关。")
GATE.asserted("CONTROL3 the near-zero mediation has a mechanism: the outcomes are near-orthogonal",
              abs(r_out)<0.15, f"corr(shame, changeability) = {r_out:+.4f}", kind="control")
json.dump(dict(r_outcomes=r_out), open(HERE/'results/outcome_corr.json','w'), indent=1)
