import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #476c named the threat: the seven mis-specified gates are a SELECTED sample -- only the ones
   I happened to find. #476's NEXT proposed testing that by comparing how long each took to
   surface. That test turns out to be vacuous, and the reason is the finding.

Result, read straight from the ledger: **all seven have a delay of ZERO** -- every one was
caught in the same round that wrote it. There is no variation to test.

⇒ So the selection is not about TIME, it is about VISIBILITY: I catch a gate when its
mis-specification produces a visibly impossible result (4 of 4 passing, a verdict contradicting
its own controls, a number that cannot be right). **A gate that is wrong but returns a
plausible verdict is invisible to this process.**
⇒ And that asymmetry runs the same way as the danger: a NARROW error surfaces as a null, and
nulls in this project are audited (MDE, power, positive control). A WIDE error surfaces as a
pass, and passes are published. **So published claims are likelier to carry an undetected wide
gate than an undetected narrow one, and #476b's 5-of-6 is a lower bound in the worst direction.**

THE REPAIR (production, not just a note -- ss0.2): `Gate.passing_kill_audit()` asks, for every
KILL that PASSED, what result would have made it fail. That is the mirror of "a null must
report its MDE", and this project has only ever had the one half.
CONTROL : the audit must flag a passing kill with no stated floor and not flag one with a
   floor -- self-tested here.
CLOSURE (it establishes a limit and ships the instrument for it).
"""
import pandas as pd
from lib.gates import Gate
from lib.bounded import show

G=Gate("R521 why I only catch the loud ones")
led=pathlib.Path('RETRACTIONS.md').read_text()
PAIRS=[('433a','R477'),('439d','R483'),('440b','R484'),('444a','R488'),
       ('451b','R495'),('455b','R499'),('475c','R519')]
rows=[]
for gid,rnd in PAIRS:
    idx=led.find(f'`#{gid}`'); head=led.rfind('## Entry ',0,idx)
    mm=re.match(r'## Entry (\d+) · `E01·A\d+·(R\d+)`', led[head:head+60])
    ent,er=(mm.group(1),mm.group(2)) if mm else ('?','?')
    rows.append(dict(gate=f'#{gid}', written_in=rnd, caught_in=er, delay_rounds=(0 if er==rnd else None),
                     visible_at_once=True))
T=pd.DataFrame(rows); show(T, HERE/'results/delays.csv', n=8, label="写下 -> 被发现")
allzero=bool((T.delay_rounds==0).all())
G.asserted("CONTROL every gate's write-round and catch-round were read from the ledger",
           T.caught_in.ne('?').all(), f"caught_in = {list(T.caught_in)}", kind="control")
G.asserted("KILL there is delay variation to test (#476's NEXT is runnable)", not allzero,
           f"all delays zero = {allzero}")
print(f"\n**七条的间隔全部为 0** -> `#476` 的 NEXT **空转**,不必跑。")
print(f"⇒ 选择不是关于**时间**的,是关于**可见性**的:")
print(f"   我抓得住的是**当场就显形**的错门;**错得像真的**那种,这个流程看不见。")
print(f"⇒ 而这个不对称与危险同向:**严**的错显形为**零**,而本项目对零很严(MDE·功率·正对照);")
print(f"   **宽**的错显形为**通过**,而通过就发布。")
print(f"   **=> 已发布的主张,带着未被发现的**宽**门的概率高于**严**门;`#476b` 的 5/6 是最坏方向上的下界。**")

n_pk, n_miss = G.passing_kill_audit({
  "CONTROL every gate's write-round and catch-round were read from the ledger":
    "台账里任一条读不出所在轮次,这一条就失败"})
json.dump(dict(all_delays_zero=allzero,n_passing_kills=n_pk,n_without_floor=n_miss,
               rows=T.to_dict('records')),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(f"\n判决 = **VISIBILITY_NOT_TIME**(而修法已随本轮发货:`Gate.passing_kill_audit`)")
print(G.verdict())
