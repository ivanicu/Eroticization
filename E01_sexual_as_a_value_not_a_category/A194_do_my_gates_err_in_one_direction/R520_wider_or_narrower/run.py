import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #475c claimed all seven mis-specified gates this session erred in the SAME direction --
   the criterion always easier to satisfy than the gate's own name -- and concluded that is a
   directional bias rather than random wording noise. #475's NEXT: read the seven verbatim
   from the ledger and judge each.

Worlds
  A  all seven wider -> a directional bias; the repair must be structural.
  B  some narrower -> it is not uniformly directional, and #475c must be narrowed or withdrawn.

⚠ THE SOURCE IS THE LEDGER, NOT MY MEMORY (#452a #459a #462b): each passage is extracted from
   RETRACTIONS.md and printed, and the judgement is made against the printed text.
⚠ The judgement (wider / narrower / equivalent) is MINE ALONE and is labelled as such -- there
   is no second reader here, so this is a description of my own gates, not an audited one.
CONTROL : each of the seven must actually be found in the ledger; a missing one is not silently
   dropped.
CLOSURE.
"""
import pandas as pd
from lib.gates import Gate
from lib.bounded import show
G=Gate("R520 wider or narrower")

led=pathlib.Path('RETRACTIONS.md').read_text()
IDS=['433a','439d','440b','444a','451b','455b','475c']
found={}
for i in IDS:
    m=re.search(r'\*\*`#'+i+r'`[^\n]*\n(?:[^\n]*\n){0,4}', led)
    found[i]=m.group(0).strip() if m else None
missing=[i for i,v in found.items() if v is None]
G.asserted("CONTROL every one of the seven is found in the ledger",
           not missing, f"missing = {missing}", kind="control")
print(f"从台账取到 **{len(IDS)-len(missing)}/{len(IDS)}** 条原文")

# 手判,依据是上面打印的原文;kind 取自该门在原轮里的角色
J={
 '433a':('KILL','更宽','名字=「系数越族内阈」;判据把 **z 分数**比到**系数尺度**的阈 -> 四个全过,应为 1/4'),
 '439d':('KILL','更宽','名字=「这条路经由经验」;判据只测「间接 ≠ 0」,不要求同向与实质'),
 '440b':('KILL','更宽','名字=「某一块的移除超过它自己的展布」;判据把 **32 取最大**比到 1 个 sd'),
 '444a':('KILL','更宽','名字=「两半由不同的块构成」;判据 |r|<0.5 **按构造必然成立** -> 不可能失败'),
 '451b':('对照','**更严**','名字=「晚起始 → 更少羞耻」= b<0;判据写成 **b>0** -> **在正确的数据上失败**'),
 '455b':('KILL','**更严**','名字=「两条都单调」;判据 `rho>=0.8` **排除完美负单调** -> 在真的单调上失败'),
 '475c':('KILL','更宽','名字=「常见组斜率真的为负」;判据「**任一**负的 k 越阈」-> 一格即判,忽略反向的三格'),
}
rows=[dict(id=f'#{i}', kind=J[i][0], verdict=J[i][1], why=J[i][2],
           excerpt=(found[i] or '')[:90].replace('\n',' ')) for i in IDS]
T=pd.DataFrame(rows)
show(T[['id','kind','verdict']], HERE/'results/gate_direction.csv', n=8, label="七条门的方向")
T.to_csv(HERE/'results/full_with_excerpts.csv',index=False)
for _,r in T.iterrows(): print(f"   {r['id']:<7}{r['kind']:<6}{r['verdict']:<8}{r['why']}")

wide=int((T['verdict']=='更宽').sum()); narrow=int((T['verdict']=='**更严**').sum())
kills=T[T['kind']=='KILL']; kw=int((kills['verdict']=='更宽').sum())
print(f"\n**更宽 {wide} · 更严 {narrow}**(共 {len(T)})")
print(f"⇒ **`#475c` 说「七次全是更宽、没有一次反过来」—— 那是**假的**。**")
print(f"而**在 KILL 里**:更宽 **{kw}/{len(kills)}**;两条「更严」的是**对照**与一个**方向判据**。")
G.asserted("KILL all seven err in the same direction (world A)", narrow==0,
           f"wider {wide}, narrower {narrow} -> #475c's claim is false as stated")
verdict = "DIRECTIONAL" if narrow==0 else "NOT_UNIFORM"
print(f"\n判决 = **{verdict}**  (`#475c` 主张 DIRECTIONAL)")
json.dump(dict(verdict=verdict,wide=wide,narrow=narrow,kill_wide=kw,n_kill=len(kills),
               rows=T.drop(columns=['excerpt']).to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(G.verdict())
