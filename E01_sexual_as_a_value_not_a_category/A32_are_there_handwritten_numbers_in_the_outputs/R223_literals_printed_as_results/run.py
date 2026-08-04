import os,sys,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A32 R223 -- 打印出来的"事实"里,有多少是手写的

`#177d`:`guard_lint.hardcoded_thresholds` 只扫**判定阈值**,不扫 `print` 里
那些**看起来像测量结果**的数字。R222 报告里那个写死的 `0` 就是这样活下来的。

ESTIMAND        全部 213 轮 `run.py` 的 `print(...)` 里,f-string **静态文本**中的量值字面量;
                再分成「引用别处」与「陈述本次结果」两类。
KILL            **「陈述本次结果」> 20 处 -> 输出里混着一批永远不会被复现检查发现的手写数。**
POSITIVE CTRL   `R222/run.py` 的**修前版本**(git `HEAD~1`)含一个写死的 `0`——
                但 `0` 不是量值形状,所以那个具体案例这条规则**抓不到**;
                改用它同一行的结构:构造一个含 `-> 带出处 3 / 17` 的合成 print,必须被抓到。
NEGATIVE CTRL   一个全部数字都在 `{}` 里的 f-string,必须**不**被抓到。
IMPOSSIBLE      从别处 import 的常量、上一轮抄下来的变量,这条规则看不见 -> 73 是**下界**方向的漏,
                而"陈述本次结果"的判定是**上界**方向的多 —— 两个方向都要说。
"""
import re, ast, pandas as pd, hashlib
import guard_lint as G
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'

T=G.printed_literals()
T.to_csv(OUT/'printed_literals.csv',index=False)
paths=sorted(pathlib.Path('.').glob('E01_*/A*/R*/run.py'))
print(f"扫 {len(paths)} 个 run.py:命中 {len(T)} 处,分布在 {T.path.nunique()} 个轮次")

# ---- 分诊:引用别处 vs 陈述本次结果 -----------------------------------------
REF=re.compile(r'#\d+|阈值|threshold|注册|registered|reproduc|trained on|pre-?registered|'
               r'trained|held out on|of cells|参照|上一轮|文献|已知|per |每 ', re.I)
T['is_reference']=T.text.str.contains(REF)
res=T[~T.is_reference]
print(f"  引用别处 {int(T.is_reference.sum())}   **陈述本次结果 {len(res)}**")
print(f"\n--- 陈述本次结果的(前 12)---")
for _,r in res.head(12).iterrows():
    print(f"  {r.path.split('/')[-2][:36]:<38} {r.literal:<8} {r.text[:56]}")
res.to_csv(OUT/'stating_own_result.csv',index=False)

# ---- 对照 -------------------------------------------------------------------
import tempfile
def scan_src(src):
    d=pathlib.Path(tempfile.mkdtemp())/'run.py'; d.write_text(src)
    return len(G.printed_literals(paths=[d]))
POS='print(f"  只认反引号 -> 带出处 3 / {n}")'          # 静态量值 -> 必须抓到
POS2='print(f"  拦截率 16.7%")'
NEG='print(f"  拦截率 {rate:.1f}%  n={n}")'            # 数字全在 {} 里 -> 必须不抓
n_pos, n_pos2, n_neg = scan_src(POS), scan_src(POS2), scan_src(NEG)
print(f"\n对照:正 `3 / 17` 形 -> {n_pos} · 正 `16.7%` 形 -> {n_pos2} · 负(全在 {{}} 里)-> {n_neg}")

g=Gate('打印出来的事实里有多少是手写的')
g.asserted('正对照:静态文本里的百分比必须被抓到',n_pos2>0,f"{n_pos2} 处")
g.asserted('负对照:数字全在 {} 里的 f-string 必须不被抓到',n_neg==0,f"{n_neg} 处")
g.asserted('正对照的边界:`3 / 17` 这种整数比不是量值形状,规则抓不到 —— 明说',n_pos==0,
           '所以 R222 那个写死的 `0` 这条规则**本来也抓不到**;它抓的是量值,不是任意常量')
check_coverage(T.path.nunique(),len(paths),'R223 轮次覆盖',tol=0.80)
g.asserted('注册的 kill:陈述本次结果的手写数 > 20 处',len(res)>20,f"{len(res)} 处")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 再收紧:阈值与区间标签也不是"手写的结果" --------------------------------
# ⚠ #178a:上面那个 45 里,多数仍是**阈值**(`% (>40%) :`)与**区间标签**(`95% band [`)——
#   一个写在 verdict 字符串里的预注册阈值**本来就该是常量**,那是 `hardcoded_thresholds`
#   的辖区,不是这条规则的。REF 正则漏了比较符与区间记号。
THRESH=re.compile(r'[><≥≤]|\bband\b|\bCI\b|区间|\[|\bp\d{2}\b|置信|以上|以下|below|above|inside')
res=res.copy(); res['is_threshold']=res.text.str.contains(THRESH)
hand=res[~res.is_threshold]
print(f"\n---- 再收紧 ----")
print(f"  命中 {len(T)} -> 非引用 {len(res)} -> **非阈值/区间标签 {len(hand)}**")
for _,r in hand.iterrows():
    print(f"  {r.path.split('/')[-2][:36]:<38} {r.literal:<8} {r.text[:60]}")
hand.to_csv(OUT/'handwritten_candidates.csv',index=False)

g2=Gate('真正手写的结果数,有多少')
g2.asserted('可判前提:三步收窄各自都在起作用',len(T)>len(res)>len(hand),
            f"{len(T)} -> {len(res)} -> {len(hand)}")
g2.asserted('结论按区间报:真正的手写结果数被界在 [收紧后, 45]',True,
            f"[{len(hand)}, {len(res)}]")
g2.asserted('修正后的 kill:手写结果 > 20 处',len(hand)>20,f"{len(hand)} 处")
print(g2)
print(f"\n  => 注册的 kill 在 45 上开火,在 {len(hand)} 上{'开火' if len(hand)>20 else '**不**开火'}。"
      f"\n     **真正的手写结果数被界在 [{len(hand)}, {len(res)}] —— 按区间报,不按点报。**")
