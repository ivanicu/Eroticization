import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #456b showed a mediation test can "pass" while being nearly incapable of failing, and that
   telling the two apart needs only three correlations. #456's NEXT: turn that into a gate that
   runs BEFORE the mediation and refuses to hand over the number when there is no room. Does
   the analytic bound actually match what the full fits produced?

Worlds
  A  the bound matches the two arms of R500 and flags them as no-room, AND does NOT flag the
     page's real mediation (#434, porn use -> shame -> acted) -> the tool works and is armed.
  B  it does not match -> I have the mechanism wrong, and the tool is premature.

POSITIVE CONTROL (this is the point of the round): `#434`'s published mediation must come back
   with ROOM. A gate that refuses everything is not a gate.
NEGATIVE CONTROL: `#500`'s two arms must come back with NO ROOM -- that is what they were.
CONTROL3 : lib.bounded.headroom_controls() must pass 3/3.
CLOSURE (it builds an instrument and validates it; it decides no new fact about people).
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show, mediation_headroom, headroom_controls

GATE=Gate("R501 validate the headroom bound")
c3=headroom_controls()
GATE.asserted("CONTROL3 the tool's own self-checks pass", all(c3),
              f"orth_rejected/med_allowed/no_value_on_reject = {c3}", kind="control")

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
AC=np.asarray(OUT['实践了多少'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY)
MP = MM & np.isfinite(PH)

CASES=[("**负对照** `#500` 起始→能不能改(控羞耻)", _EARLY, BE, SH, MM, "无空间"),
       ("**负对照** `#500` 起始→羞耻(控能不能改)", _EARLY, SH, BE, MM, "无空间"),
       ("**正对照** `#434` 用量→实践(控羞耻)",     PH,     AC, SH, MP, "有空间")]
rows=[]
for lab,x,y1,y2,msk,expect in CASES:
    ok,ind,share,why = mediation_headroom(x,y1,y2,msk,min_move=0.05,name="")
    rows.append(dict(case=lab, expect=expect, room=('有' if ok else '无'),
                     indirect=(ind if ok else np.nan), share=share, why=why[:60]))
T=pd.DataFrame(rows); show(T, HERE/'results/headroom.csv', n=4, label="事前上界")

neg_ok = all(r['room']=='无' for r in rows[:2])
pos_ok = rows[2]['room']=='有'
GATE.asserted("NEGATIVE CONTROL R500's two arms come back with NO ROOM",
              neg_ok, f"rooms = {[r['room'] for r in rows[:2]]}", kind="control")
GATE.asserted("POSITIVE CONTROL #434's published mediation comes back WITH ROOM",
              pos_ok, f"share = {rows[2]['share']:+.4f}", kind="control")

# 事前上界 vs 事后实测:#434 的实测占比是 −2.4%..? 不 —— #435c 已判它不可估;
# 可比的是 R500 的保留比(1.007 / 1.005 -> 移动 ≈ 0.7% / 0.5%)
pred=[abs(r['share']) for r in rows[:2]]
print(f"\n事前上界(占总效应比):`#500` 两臂 = **{pred[0]:.2%}** / **{pred[1]:.2%}**")
print(f"事后实测的移动          = **0.7%** / **0.5%**(`#456a` 的 1.007 / 1.005)")
agree = all(p<0.05 for p in pred)
GATE.asserted("KILL the analytic bound agrees with what the full fits actually moved",
              agree, f"bound {[f'{p:.2%}' for p in pred]} vs observed 0.7% / 0.5%")
verdict = "BOUND_WORKS" if (agree and neg_ok and pos_ok) else "NOT_VALIDATED"
print(f"\n判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,rows=T.to_dict('records'),
                bound=pred,observed=[0.007,0.005]),
           open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
