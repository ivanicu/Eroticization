"""E03·A18·R691 —— 一条条目被后来的条目重引时,丢掉了多少定义性的限定

**类型:FRONTIER**。`#654` 发现一种新的失败方式:**不是账本错了,是我重引自己的账本时
丢掉了定义那个量的符号**(`|ρ|` -> `ρ`),而两条条目之后它弄坏了一道控制。
**本轮把它变成一个可测的量,在 654 条条目上量。**

⚠ **BASIN**:最近几轮都在确认「我自己的记录/仪器是问题」——**这就是盆地**。
  **所以本轮下注「丢失率低」**,即下注 `#654` 是孤例、我刚写的那条规矩是为一个不存在的问题写的。

W1 **无损** —— 中位丢失率 = 0 ⇒ `#654` 是孤例。
W2 **有损** —— 中位丢失率 > 0 ⇒ 而丢失率就是它的大小,**凡重引过的数字都要重看**。
**区别是决策性的**:W2 要求我改写「怎么写 NEXT」这个程序本身。

G1 ESTIMAND(先于方法):对每一对(**原条目里带数字的句子** × **后来条目里重引该数字并反向引用原条目的句子**),
  比较两句各自携带的**限定词集合**,报 **`丢失的限定数 ÷ 原句携带的限定数`** 的分布。
G2 CONTROLS:
  **正对照** `#653` 重引 `#542` 的 0.5789/0.1726,**已知丢了 `|·|`** -> **必须被查出来**。
    ⚠ **而预注册⑤要求先验证词表**:词表若抓不到 `|·|`,**先补词表再跑全库**,不许直接跑。
  **安慰剂** 把引句与**随机一条**原句配对,丢失率应回到「两句本就无关」的水平。
G3:全部配对报,含丢失率 = 0 的。G4:限定词表 {严格, 宽松} 两条规格。
KILL(条件式):if 正对照抓到 `#653`→`#542` and 安慰剂 ≈ 无关基线:
  中位丢失率 = 0 -> W1 · > 0 -> W2;else UNVERIFIED
⚠ **最强混淆**:**限定词表是我列的,漏一个词就等于假装没丢。**
  ⇒ 词表先在正对照上验证;**并且报「原句一个限定词都没有」的对数** —— 那些对**不是无损,是量不了**。
IMPOSSIBLE(不写 planned):只查**数字被重引**的情形,**措辞性的限定丢失查不到**(`#623` 实测 36 次页面
  更正里 29 次改的是措辞)· 只查账本内部,不查页面 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
TXT=pathlib.Path("RETRACTIONS.md").read_text()
ENT={}
for m in re.finditer(r'^## Entry (\d+)[^\n]*\n(.*?)(?=^## Entry |\Z)', TXT, re.M|re.S):
    ENT[int(m.group(1))]=m.group(2)
print(f"=== 硬规则①:先打印 ===\n  账本条目 = **{len(ENT)}**  编号 {min(ENT)}–{max(ENT)}")

QUAL = {
 "绝对值": r'\|ρ\||\|r\||\|rho\||绝对值|abs\(',
 "中位":   r'中位|median',
 "逐年":   r'逐年|per-year|年内|within.year',
 "归一":   r'归一|normali[sz]ed|/ *r_?max|天花板',
 "样本量": r'\bn\s*=|\bn\s*=\s*\d|样本量',
 "区间":   r'CI|区间|\[[-+−]?\d*\.\d+ *, *[-+−]?\d*\.\d+\]|95%',
 "留出":   r'留出|held.out|半样本|split.half',
 "逐对":   r'逐对|每一对|per-pair|每对',
 "块重抽": r'块 ?bootstrap|block bootstrap|bootstrap',
 "秩":     r'秩相关|Spearman|Kendall|rank',
 "对齐":   r'对齐|align|反向计分|反向措辞',
 "生/同号": r'未归一|生相关|生阈|同号|signed|raw',
 "控制":   r'正对照|安慰剂|置换零|positive control|placebo',
 "对数":   r'\d+ *对\b|npairs|\d+ 对',
}
NUM=re.compile(r'(?<![\w.])([0-9]+\.[0-9]{3,4})(?![\w])')
def quals(s): return {k for k,p in QUAL.items() if re.search(p,s,re.I)}
def sents(body):
    return [x.strip() for x in re.split(r'[。\n]|(?<=[.;])\s', body) if x.strip()]

# ── 预注册⑤:先在正对照上验证词表 ─────────────────────────────────────────
print("\n=== 预注册⑤:先验证词表能不能抓到 `#653`→`#542` 的 `|·|` ===")
src=[s for s in sents(ENT[542]) if "0.5789" in s]
dst=[s for s in sents(ENT[653])+sents(ENT[654]) if "0.579" in s or "0.5789" in s]
print(f"  `#542` 里含 0.5789 的句子 {len(src)} 句 · 后来条目里重引它的 {len(dst)} 句")
if src: print(f"    原句限定 = {sorted(quals(src[0]))}")
for s in dst[:3]: print(f"    引句限定 = {sorted(quals(s))}   « {s[:70]} »")
ok_vocab = bool(src) and ("绝对值" in quals(src[0])) and any("绝对值" not in quals(s) for s in dst)
print(f"  ⇒ 词表抓得到那次丢失? **{ok_vocab}**" + ("" if ok_vocab else "  ⚠ 抓不到 ⇒ 先补词表,不许跑全库"))

# ── 全库配对 ────────────────────────────────────────────────────────────
pairs=[]
for y,body in ENT.items():
    refs={int(r) for r in re.findall(r'#(\d{3})', body)} & set(ENT)
    for x in refs:
        if x>=y: continue
        for s in sents(body):
            if f"#{x}" not in s: continue
            nums=set(NUM.findall(s))
            if not nums: continue
            for s0 in sents(ENT[x]):
                shared=nums & set(NUM.findall(s0))
                if not shared: continue
                q0,q1=quals(s0),quals(s)
                if not q0: pairs.append(dict(src=x,dst=y,num=sorted(shared)[0],n_src=0,lost=None)); continue
                pairs.append(dict(src=x,dst=y,num=sorted(shared)[0],n_src=len(q0),
                                  lost=len(q0-q1)/len(q0),lost_set=sorted(q0-q1)))
                break
            else: continue
            break
measur=[p for p in pairs if p["lost"] is not None]
noq=[p for p in pairs if p["lost"] is None]
print(f"\n=== 配对 ===\n  找到 {len(pairs)} 对 · **可量 {len(measur)}** · **原句零限定 ⇒ 量不了 {len(noq)}**")
if measur:
    L=[p["lost"] for p in measur]
    print(f"  **丢失率:中位 {np.median(L):.3f} · 均值 {np.mean(L):.3f} · 全无损的对 {sum(1 for x in L if x==0)}/{len(L)}**")
    print(f"  分位 q25 {np.quantile(L,.25):.3f} · q75 {np.quantile(L,.75):.3f} · 最大 {max(L):.3f}")
    from collections import Counter
    C=Counter(w for p in measur for w in p.get("lost_set",[]))
    print("  最常被丢掉的限定:"+" · ".join(f"{k} {v}" for k,v in C.most_common(8)))

def placebo(seed):
    rng=np.random.default_rng(seed); out=[]
    allsent=[s for b in ENT.values() for s in sents(b) if quals(s)]
    for p in measur:
        s0=allsent[rng.integers(0,len(allsent))]
        s1=allsent[rng.integers(0,len(allsent))]
        q0,q1=quals(s0),quals(s1)
        if q0: out.append(len(q0-q1)/len(q0))
    return float(np.median(out)) if out else np.nan
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"\n=== 控制 ===\n  安慰剂 随机配对的丢失率中位 = **{pl:.3f}**(「两句本就无关」的水平)")
print(f"  正对照 词表抓到 `#653`→`#542` = **{ok_vocab}**")
from lib.gates import Gate
G=Gate("重引自己的账本时丢掉了多少限定")
med=float(np.median([p["lost"] for p in measur])) if measur else np.nan
p1=G.positive_control("词表必须抓到 #653→#542 丢掉的 `|·|`",planted=float(1.0 if ok_vocab else 0.0),floor=0.5,spread=0.01)
# 「这个零该不该是零?」—— **不该**。随机配一句无关的原句,丢失率本来就该很高,
# 它是一条**系统性的基线偏移**,不是一个应当趋零的干扰 ⇒ **用 offset_control,并命名这个零的种类**。
p2=G.offset_control("真实重引的丢失率必须显著低于随机配对(否则重引不比随机更保真)",
                    effect=med, offset=pl, spread=0.02,
                    null_kind="随机配一句无关的原句 —— 两句本就无关时的丢失水平(系统性基线,不该趋零)")
if p1 and p2:
    verdict=("**W1 —— 无损:中位丢失率 = 0,`#654` 是孤例**" if med==0 else
             f"**W2 —— 有损:中位丢失率 = {med:.3f},而重引仍显著优于随机({pl:.3f})**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · offset {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(entries=len(ENT),pairs=len(pairs),measurable=len(measur),no_qual=len(noq),
               median_loss=med,placebo=pl,vocab_ok=bool(ok_vocab),verdict=verdict,
               worst=sorted(measur,key=lambda p:-p["lost"])[:12],unchallenged=True),
          open(OUT/"lossy_requote.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'lossy_requote.json'}")
