import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #473c: a cross-outcome reference arm answers two questions at once -- "is this moderated"
   and "is this a measurement artifact". Which of this project's block-splitting conclusions
   ever ran one?

BOUNDARY, written before the scan (this is the NINTH scanning round; seven of the previous
eight over-indicted):
  IN SCOPE : the round's MAIN QUANTITY comes from comparing two groups of BLOCKS.
  OUT      : block grouping used only as a control or a matching device.
  OUT      : a round whose conclusion is itself a statement about blocks rather than about a
             relationship measured within them.
⚠ AND THE THING THAT MAKES THIS ONE DIFFERENT, stated up front: the candidate set is EIGHT.
   Every one is read by hand, in full. The eight earlier scans had candidate sets of 24 to 631
   and were judged by a text proxy -- which is exactly when they over-indicted.

Worlds
  A  every in-scope conclusion has a reference arm (or is itself cross-outcome) -> closed.
  B  some do not -> each is a #473b-shaped risk: a gate that passes while pointing at the
     wrong thing.
CONTROL : the two known-good cases (#472's design, #473's design) must come back as HAVING a
   reference arm -- a scan that cannot see the ones built with it is measuring nothing.
CLOSURE.
"""
import pandas as pd
from lib.gates import Gate
from lib.bounded import show
G=Gate("R518 reference-arm audit")

SPLIT=re.compile(r'order\[:\s*k\]|order\[-\s*k\s*:\]|hi_b|lo_b|块分成两组|按块.*分组')
cands=[]
tot=0
for p in sorted(pathlib.Path('.').rglob('run.py')):
    if '.git' in str(p): continue
    # ⚠ #474b: the first run found NINE candidates, and the ninth was THIS FILE -- the regex
    # literals above are themselves a match. A scan that walks the directory it lives in will
    # find itself, and the ninth scanning round's artifact is a new one: self-inclusion.
    if p.resolve()==pathlib.Path(__file__).resolve(): continue
    tot+=1; t=p.read_text(errors='ignore')
    if SPLIT.search(t): cands.append((p.parts[-2], t))
print(f"全部轮次 **{tot}** · 候选 **{len(cands)}** -> **小到可以整批手读**(这正是本轮与前八次的差别)")

# 手判,依据写在每一行里,并持久化(#431c)
VERDICT={
 'R445_shame_and_healing':      ('不在范围内','主量是**人层**羞耻×治疗性的联合分布,不是块分两组'),
 'R174_is_the_trait_link_just_coverage':('不在范围内','块数只用作**卡钳匹配**的协变量,不是被比较的两组'),
 'R453_split_S_in_two':         ('在范围内 · 已覆盖','把块拆成越轨/普通两半,而主量**本身跨 29 个结局** -> 参照臂问题不适用'),
 'R513_rarity_shame_by_domain_commonness':('在范围内 · 事后覆盖','当时只跑羞耻;**`R516` 用同一设计跑了第二个结局**,并给出分离'),
 'R514_equal_sized_groups':     ('在范围内 · 事后覆盖','同上,`R516` 的 k 扫描覆盖了它'),
 'R515_what_happens_at_k12':    ('不在范围内','诊断轮:问的是 k=12 那一格的来源,不作结局主张'),
 'R516_same_design_other_outcome':('在范围内 · **自带**','它**就是**参照臂那一轮'),
 'R517_is_the_other_path_moderated':('在范围内 · **自带**','参照臂正是它否掉自己的门的那个东西'),
}
rows=[]
for nm,t in cands:
    v,why=VERDICT.get(nm,('**未判定**','不在手判清单里 -> 本轮遗漏,必须补'))
    outs=sorted(set(re.findall(r'OUT\[.(羞耻|能不能改|治疗性|实践了多少).\]',t)))
    rows.append(dict(round=nm[:44], verdict=v, outcomes_in_code=len(outs), why=why))
T=pd.DataFrame(rows); show(T[['round','verdict','outcomes_in_code']],
                           HERE/'results/reference_arm_audit.csv', n=8, label="八轮手判")
T.to_csv(HERE/'results/full_with_reasons.csv',index=False)
for _,r in T.iterrows(): print(f"   · {r['round']:<42} {r['verdict']:<16} {r['why']}")

undone=[r for _,r in T.iterrows() if '未判定' in r['verdict']]
G.asserted("coverage: every candidate was hand-judged",
           len(undone)==0, f"unjudged = {[r['round'] for r in undone]}", kind="control")
known=[r for _,r in T.iterrows() if r['round'].startswith(('R516','R517'))]
G.asserted("CONTROL the two rounds built WITH a reference arm are seen as having one",
           all('自带' in r['verdict'] for r in known),
           f"{[(r['round'],r['verdict']) for r in known]}", kind="control")
gaps=[r for _,r in T.iterrows() if r['verdict'].startswith('在范围内') and '自带' not in r['verdict']
      and '覆盖' not in r['verdict']]
G.asserted("KILL every in-scope conclusion has a reference arm or is cross-outcome",
           len(gaps)==0, f"gaps = {[r['round'] for r in gaps]}")
verdict = "CLOSED" if len(gaps)==0 else "GAPS_FOUND"
print(f"\n在范围内 **{int(T.verdict.str.startswith('在范围内').sum())}** · 假阳 **{int((T.verdict=='不在范围内').sum())}** · "
      f"**缺参照臂 = {len(gaps)}** -> 判决 = **{verdict}**")
json.dump(dict(verdict=verdict,total_rounds=tot,candidates=len(cands),
               gaps=[r['round'] for r in gaps],rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(G.verdict())
