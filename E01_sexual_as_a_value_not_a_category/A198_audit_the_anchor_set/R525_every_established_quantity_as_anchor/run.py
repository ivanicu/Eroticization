import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #480b: the bar of my anchoring rule had been audited twice, but the ANCHOR SET never had --
   it only ever held counts, which is why `powerlessness` looked unanchorable for two rounds.
   #480's NEXT: put every direction-established quantity into the anchor set and rerun every
   variable whose direction is still only inferred.

Candidate set = **4** (openness, conscientiousness, neuroticism, agreeableness). Small enough
to hand-read in full, which is the condition #474a set for trusting this kind of sweep.

⚠ PRE-CHECK THAT NARROWS THE ANCHOR SET AGAIN, and it must be applied before any number is
   read: all four candidates share a scoring pipeline with `extroversion` (#429d: identical n,
   range and step). A shared coding convention would reproduce the expected sign whatever the
   truth is, so **`extroversion` is NOT a valid anchor for them** -- it is inside the same
   instrument. Valid anchors are the ones outside that pipeline.

⚠ SIGN PRIORS, WRITTEN HERE BEFORE ANY CORRELATION IS COMPUTED (this is the repair of #480c,
   where the prior arrived after the number):
     openness           x count-of-interests : weakly POSITIVE (openness to experience covers
                        breadth of interests) -- defensible but weak.
     conscientiousness  x any anchor here     : **NO DEFENSIBLE PRIOR** -- stated now, so a
                        later "it correlates with X" cannot be dressed as confirmation.
     neuroticism        x any anchor here     : **NO DEFENSIBLE PRIOR** against counts of
                        sexual interests.
     agreeableness      x any anchor here     : **NO DEFENSIBLE PRIOR**.
   => at most ONE of the four can have its direction read even if all four clear the bar. The
   rest can reach world C from #461a: power without a sign.

Worlds
  A  some clear AND have a prior -> their direction moves from inferred to readable.
  B  they clear without a prior -> "power but no sign", which is a sharper UNVERIFIED than the
     current one.
  C  none clears -> the inherited-direction inference (#429d) remains the only thing available.
CONTROL : `powerlessness` must clear via `extroversion` here, reproducing #480b -- it is the
   known-positive case for the widened anchor set.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show, anchor_rule, PositiveControlFailed

GATE=Gate("R525 audit the anchor set")
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
PIPE=['opennessvariable','consciensiousnessvariable','neuroticismvariable','agreeablenessvariable']
CAND={k.replace('variable',''):num(k) for k in PIPE}
ANCH_OUT={'性行为计数':num('Totalsexacts'),'恋物类别数':num('totalfetishcategory'),
          '性伴数分档':raw['sexcount'].astype(str).map(BANDS).values.astype(float),
          'pornhabit':num('pornhabit'),'bondageaverage':num('bondageaverage')}
PC={'pornhabit':num('pornhabit'),'bondageaverage':num('bondageaverage')}
FULL=np.ones(len(raw),dtype=bool)
print("⚠ **前置**:四个候选与 `extroversion` **同一条计分流水线**(`#429d`)->")
print("   **`extroversion` 对它们不是合法的锚**(共享编码约定会复现期待的符号,无论真相如何)。")
print(f"   合法锚集(流水线之外)= {list(ANCH_OUT)}")

res,pc = anchor_rule(CAND, ANCH_OUT, PC, FULL)
rows=[]
PRIOR={'openness':'弱正(开放性涵盖兴趣广度)—— 可辩护但弱',
       'consciensiousness':'**无可辩护先验**(跑之前已写下)',
       'neuroticism':'**无可辩护先验**(跑之前已写下)',
       'agreeableness':'**无可辩护先验**(跑之前已写下)'}
for k,(ok,best,rs) in res.items():
    top=max(rs.items(), key=lambda kv: abs(kv[1]) if kv[1]==kv[1] else -1)
    rows.append(dict(variable=k, clears=ok, best=best, best_anchor=top[0],
                     prior=PRIOR.get(k,'—')))
T=pd.DataFrame(rows)
show(T[['variable','clears','best','best_anchor']], HERE/'results/candidates.csv', n=6,
     label="四个候选 × 流水线外的锚")
T.to_csv(HERE/'results/with_priors.csv',index=False)
for _,r in T.iterrows():
    print(f"   {r['variable']:<18} 最强 **{r['best']:.3f}**({r['best_anchor']}) · "
          f"{'过杠' if r['clears'] else '不过'} · 先验:{r['prior']}")

# CONTROL:powerlessness 在加宽后的锚集里必须过杠(复现 #480b)
res2,_ = anchor_rule({'powerlessness':num('powerlessnessvariable')},
                     dict(ANCH_OUT, **{'extroversion(方向已确立)':num('extroversionvariable')}),
                     PC, FULL)
pw_ok, pw_best, _ = res2['powerlessness']
GATE.asserted("CONTROL powerlessness clears via the widened anchor set (reproduces #480b)",
              pw_ok and abs(pw_best-0.174)<0.02,
              f"powerlessness best = {pw_best:.3f}, clears = {pw_ok}", kind="control")
GATE.asserted("CONTROL extroversion is excluded from these four's anchor set",
              'extroversion' not in str(list(ANCH_OUT)),
              "same scoring pipeline -> not an independent anchor", kind="control")

cleared=[r['variable'] for _,r in T.iterrows() if r['clears']]
readable=[r['variable'] for _,r in T.iterrows() if r['clears'] and '无可辩护' not in r['prior']]
print(f"\n过杠的 = **{cleared or '无'}** · 其中**有事先可辩护符号先验**的 = **{readable or '无'}**")
GATE.asserted("KILL at least one of the four gets a readable direction", bool(readable),
              f"cleared {cleared}, readable {readable}")
verdict = ("READABLE" if readable else ("POWER_NO_SIGN" if cleared else "NOTHING_CLEARS"))
print(f"\n判决 = **{verdict}**")
npk,nmiss=GATE.passing_kill_audit({})
json.dump(dict(verdict=verdict,cleared=cleared,readable=readable,
               rows=T.to_dict('records'),powerlessness_best=pw_best,
               passing_kills=npk,without_floor=nmiss),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
