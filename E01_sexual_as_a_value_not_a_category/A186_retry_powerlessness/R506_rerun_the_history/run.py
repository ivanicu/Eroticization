import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #461b was the SECOND time I wrote an anchoring rule too narrow and only found out when its
   own positive control failed. #461's NEXT: make the positive control a REQUIREMENT of the
   rule, then RERUN THE HISTORY -- because a tool that is never pointed at the past is a
   promise, not a repair.

Worlds
  A  all three historical verdicts stand under the enforced rule -> the tool is a safety net
     for the future and the page needs no correction.
  B  one changes -> that verdict was a mis-call produced by a too-narrow rule, and it now has
     a name and a page fix.

The three, run on THEIR OWN masks (not mixed -- #461's NEXT says so explicitly):
  1. `#428c` full sample: the eight ordered `OTHER` variables, where the original rule was a
     SINGLE count anchor.
  2. `#442b` block mask: the three components c1/c2/c3, where the rule was TWO count anchors.
  3. `#461`  full sample: `powerlessness`, where the rule was two counts and then three.
CONTROL : `anchor_rule_controls()` 3/3, including "raises when a positive control fails".
CONTROL2: pointing the tool at the ORIGINAL narrow anchor set must RAISE -- if it does not,
   the enforcement is not doing anything.
CLOSURE unless world B fires.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show, anchor_rule, anchor_rule_controls, PositiveControlFailed

GATE=Gate("R506 rerun the history")
c=anchor_rule_controls()
GATE.asserted("CONTROL the rule's own self-checks pass", all(c),
              f"raises/passes/both = {c}", kind="control")

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
PARTNERS=raw['sexcount'].astype(str).map(BANDS).values.astype(float)
ANCH={'性行为计数':num('Totalsexacts'),'恋物类别数':num('totalfetishcategory'),
      '性伴数分档':PARTNERS}
PC={'pornhabit':num('pornhabit'),'bondageaverage':num('bondageaverage'),
    'extroversion':num('extroversionvariable')}
FULL=np.ones(len(raw),dtype=bool)

# CONTROL2:把工具指向**原来那套窄锚**,必须抛错
raised=False
try:
    anchor_rule({'x':num('powerlessnessvariable')},{'性行为计数':ANCH['性行为计数']},PC,FULL)
except PositiveControlFailed as e:
    raised=True; why=str(e)[:110]
GATE.asserted("CONTROL2 the original narrow anchor set now RAISES",
              raised, why if raised else "did not raise -- enforcement is inert", kind="control")
print(f"⚠ 指向原来的单锚集 -> **抛错**:{why if raised else '没抛(执行是空的)'}")

rows=[]
# ---- 历史 1:#428c 全样本
C1={k:num(v) for k,v in
    {'openness':'opennessvariable','conscientiousness':'consciensiousnessvariable',
     'extroversion':'extroversionvariable','neuroticism':'neuroticismvariable',
     'agreeableness':'agreeablenessvariable','powerlessness':'powerlessnessvariable',
     'pornhabit':'pornhabit','bondageaverage':'bondageaverage'}.items()}
v1,_=anchor_rule(C1,ANCH,PC,FULL)
for k,(ok,best,rs) in v1.items():
    rows.append(dict(case='#428c 全样本',variable=k,anchorable=ok,best=best))
# ---- 历史 3:#461 全样本(同一次调用即可,单列出来对照台账)
rows.append(dict(case='#461 全样本',variable='powerlessness',
                 anchorable=v1['powerlessness'][0],best=v1['powerlessness'][1]))
# ---- 历史 2:#442b 块掩码上的三个成分
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
C2={f'c{k+1}':sc(k) for k in range(3)}
v2,pc2=anchor_rule(C2,ANCH,PC,m)
for k,(ok,best,rs) in v2.items():
    rows.append(dict(case='#442b 块掩码',variable=k,anchorable=ok,best=best))
T=pd.DataFrame(rows); show(T, HERE/'results/history.csv', n=12, label="历史三次")
print("\n正对照在块掩码上的 best |r|:", {k:round(v[1],3) for k,v in pc2.items()})

# ⚠ #462b: the first version of this dict was TYPED FROM MEMORY and marked the four other Big
# Five traits as anchorable. The record says the opposite, verbatim: they "inherit it at one
# remove ... an inference about the pipeline, NOT a measurement of each trait" (#429d). They
# were never claimed anchorable, so they are NOT in the comparison set. The check's EXPECTED
# VALUES must be read from the record, exactly as its METHOD must be (#452a, #459a) -- and this
# is the first time the memory-typed part was the expectation rather than the method.
LEDGER={'#428c 全样本':{'powerlessness':False,'pornhabit':True,'bondageaverage':True,
                       'extroversion':True},
        '#442b 块掩码':{'c1':True,'c2':True,'c3':True},
        '#461 全样本':{'powerlessness':False}}
diff=[(r['case'],r['variable'],bool(r['anchorable']),LEDGER[r['case']].get(r['variable']))
      for _,r in T.iterrows()
      if r['variable'] in LEDGER[r['case']] and bool(r['anchorable'])!=LEDGER[r['case']][r['variable']]]
print(f"\n**与台账记录不一致的判定 = {len(diff)}**")
for d in diff: print(f"   ⚠ {d[0]} · {d[1]}:现在 {d[2]},台账 {d[3]}")
GATE.asserted("KILL every historical verdict stands under the enforced rule",
              len(diff)==0, f"{len(diff)} verdicts changed: {diff}")
verdict = "HISTORY_STANDS" if len(diff)==0 else "A_VERDICT_CHANGED"
print(f"\n判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,n_diff=len(diff),diff=[list(map(str,d)) for d in diff],
                rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
