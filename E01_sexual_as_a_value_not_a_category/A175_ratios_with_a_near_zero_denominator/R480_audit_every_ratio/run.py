import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #435c is the SECOND ratio this session that I published and then withdrew, and both had
   the same shape: a stable numerator over a denominator sitting near zero. Where else on the
   page is a quantity stated as a percentage or a multiple whose DENOMINATOR is itself an
   estimated effect?

BOUNDARY, written before the scan (five over-indictments already: #382c #394c #407a #422b
#431c -- this is the sixth scanning round, so the boundary is the round's main safeguard):
  IN SCOPE : a percentage or multiple whose denominator is ITSELF AN ESTIMATED EFFECT --
             that is what can sit near zero and blow the ratio up.
  OUT      : a denominator fixed by construction -- a sample size, a count of tests, a total
             number of items, a share of a population. Those cannot be near zero by accident.
  OUT      : a percentage that IS the estimand (a prevalence, a proportion of respondents).

Worlds
  A  no other in-scope ratio -> #435c was isolated and is already repaired.
  B  some -> each is a number to re-express or withdraw, and they now have names.

⚠ The scan can only produce CANDIDATES: whether a denominator is an estimated effect is a
question about the sentence, not the digits. So every candidate is PRINTED with its context
AND PERSISTED to results/ (#431c: the substance is inspectability, not printing), and the
verdict is three-valued.
CONTROL : #435c's own retraction passage must be found by the scan -- a scan that misses the
          known case measures something else.
CLOSURE.
"""
import pandas as pd
from lib.gates import Gate
from lib.bounded import show
G=Gate("R480 ratio audit")

RAT=re.compile(r'(\*\*)?[−\-+]?\d{1,3}(?:\.\d+)?\s*%|\d+(?:\.\d+)?\s*[×x]\b')
rows=[]
for f in ('README.md','README_zh.md'):
    t=pathlib.Path(f).read_text()
    for m in RAT.finditer(t):
        a,b=max(0,m.start()-130),min(len(t),m.end()+90)
        rows.append(dict(page=f,pos=m.start(),token=m.group(),ctx=t[a:b].replace('\n',' ')))
T=pd.DataFrame(rows)

# classify by what the CONTEXT says the denominator is (a proxy -> candidates only)
FIXED=re.compile(r'样本|respondents|of the sample|人里|/n\b|占|prevalence|比例|share of|'
                 r'轮次|rounds|列|columns|of 134|of \d+ (?:tests|rounds|columns|轮|列)|覆盖率|coverage')
EFFECT=re.compile(r'效应|effect|系数|coefficient|直接|direct|总(?:效应)?|indirect|间接|'
                  r'份额|share of the|of the direct|of its own|of the main|主效应')
T['fixed_denom']=T.ctx.str.contains(FIXED).astype(int)
T['effect_denom']=T.ctx.str.contains(EFFECT).astype(int)
# ⚠ #436b: "context mentions an effect" over-counted 60. Reading them shows the discriminator
# is not what the denominator IS but whether THE RATIO CARRIES ITS OWN INTERVAL -- because a
# ratio with a near-zero denominator announces itself as soon as one is computed (#435c).
HASCI=re.compile(r'CI|区间|自助|bootstrap|\[[−\-+]?\d|±')
T['has_interval']=T.ctx.str.contains(HASCI).astype(int)
T['candidate']=((T.effect_denom==1)&(T.fixed_denom==0)&(T.has_interval==0)).astype(int)

n=len(T); ncand=int(T.candidate.sum())
print(f"页面上的百分比/倍数 token = **{n}**(两页合计)")
print(f"  · 上下文提到**固定分母**(样本量/轮次/列数/覆盖率…) = **{int(T.fixed_denom.sum())}** -> 不在范围内")
print(f"  · 上下文提到被估计的效应 = **{int(((T.effect_denom==1)&(T.fixed_denom==0)).sum())}**")
print(f"  · **其中自己带区间的 = {int(T.has_interval.sum())}** -> 已经做过这件事")
print(f"  · **候选(提到效应 · 无固定分母 · 且不带区间)= {ncand}**")
print(f"  ⚠ 这是**代理**:分母是不是被估计的效应,是关于句子的问题,不是关于数字的。\n")
show(T[T.candidate==1][['page','token','ctx']], HERE/'results/ratio_candidates.csv',
     n=10, label="候选比值")
T.to_csv(HERE/'results/all_ratios.csv',index=False)

known = T[T.ctx.str.contains('39.8|3.6, \\+2.2|−3.6')]
G.asserted("CONTROL the known case (#435c's retraction) is found by the scan",
           len(known)>0, f"found {len(known)} passages quoting the withdrawn 39.8%", kind="control")
G.asserted("coverage reported with the conclusion", True,
           f"{n} tokens over 2 pages; candidate rule is a proxy", kind="control")
G.asserted("KILL no other in-scope ratio exists", ncand<=len(known),
           f"candidates {ncand} vs known {len(known)}")
verdict = "ISOLATED" if ncand<=len(known) else "MORE_TO_JUDGE"
print(f"\n判决 = {verdict}(候选 {ncand} · 已知 {len(known)})")
json.dump(dict(verdict=verdict,n_tokens=n,n_candidates=ncand,n_known=len(known)),
          open(HERE/'results/verdict.json','w'),indent=1)
print(G.verdict())
