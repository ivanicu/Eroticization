"""E02·A12·R668 —— 那个配对启发式,错配率是多少?

`#631` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。

⚠ **一处标注的偏离,写在最前面。**
   `#631` 的 NEXT 写的是「在临时分支上伪造 N=12 组」。**本轮走纯文本支。**
   理由不是省事,是 `#173c` 立下的规则本身:*「正对照就跑这一支 —— 否则对照要穿过
   git/cwd/文件名三层,任何一层出错都会让对照静默失败」*。
   同一个估计量,少穿三层。**偏离是刻意的,并且被标注。**

⚠ **§3 梯度检查:「无关段落」必须是个真陷阱。**
   它必须**携带一个新来的量** —— 否则配错进去也报不出 `替换`,错配率会假地是 0。
   ⇒ 设计:改写段**丢掉** `+0.7373` 且不增新量;无关段**增加** `+0.9999`。
   若配对挑中改写段 -> `拿走`(对);挑中无关段 -> `替换 +0.9999`(**错配**)。

G1 ESTIMAND(先于方法):`错配率(档, 版)` = 报出 `替换`/`多候选` 且其候选来自**无关段**的比例。
  三档改写强度:词重合约 **0.9 / 0.7 / 0.4**(实测每组真实重合度,一并报,不凭标称)。
KILL(条件式,预注册于 `#631`):
  if g=0(两段都只改措辞、无量变动)一行都不报:
      0.9 档与 0.7 档错配 **必须为 0**;0.4 档允许 `UNMATCHED`,**但不许错配**
      任一档出现错配 -> **记为缺陷并修**
  else: UNVERIFIED
G3:三档 × 两版 6 格 + 每组的真实词重合度。
IMPOSSIBLE(不写 planned):合成文本的词分布不等于真实页面 ·
  只验「一走零来 vs 无关段一来」这一种陷阱形状 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
norm = lambda m: m.strip().replace(' ', '').replace('倍', '×').replace('−', '-')
paras = lambda s: [x for x in re.split(r'\n\s*\n', s) if x.strip()]
words = lambda s: set(re.findall(r'\w+', s))
ov = lambda a, b: len(words(a) & words(b)) / max(len(words(a)), 1)


def judge(old, new):
    """与已接入的 `numbers_that_left` 同一套判定,跑在纯文本上(`#173c`)。"""
    o = {norm(m) for m in A._MAGNUM.findall(old)}; n = {norm(m) for m in A._MAGNUM.findall(new)}
    NP = paras(new); newcomers = n - o; out = []
    for tk in sorted(o - n):
        src = [x for x in paras(old) if tk in norm(x) or tk in x]
        rep, kind, best = None, "UNMATCHED", None
        if src and NP:
            best = max(NP, key=lambda q: ov(src[0], q))
            if ov(src[0], best) >= 0.50:
                cand = sorted({norm(m) for m in A._MAGNUM.findall(best) if norm(m) in newcomers})
                if not cand: rep, kind = None, "拿走"
                elif len(cand) == 1: rep, kind = cand[0], "替换"
                else: rep, kind = " | ".join(cand), "多候选"
        out.append((tk, kind, rep, best))
    return out


BODY = ("这一段讨论的是社会层的谴责与它的对象,句子较长以便词重合度可控,"
        "其中包含若干可以被替换掉的词语来调节相似度")
FILLER = ["替换甲", "替换乙", "替换丙", "替换丁", "替换戊", "替换己", "替换庚", "替换辛",
          "替换壬", "替换癸", "替换子", "替换丑", "替换寅", "替换卯", "替换辰", "替换巳"]
UNREL = "另一段完全无关的文字,谈的是别的题目,与上面那一段没有共同的实质内容"


def make(strength, seed):
    """strength 越低,改写越狠(词重合越低)。返回 (old, new, 标称档)。"""
    rng = np.random.default_rng(seed)
    toks = list(BODY)
    keep = int(len(toks) * strength)
    idx = sorted(rng.choice(len(toks), size=len(toks) - keep, replace=False)) if keep < len(toks) else []
    new_toks = list(toks)
    for i, j in enumerate(idx): new_toks[j] = FILLER[i % len(FILLER)][0]
    old = f"{''.join(toks)},其中一个量是 +0.7373。\n\n{UNREL}。"
    new = f"{''.join(new_toks)}。\n\n{UNREL},而这里出现了 +0.9999。"
    return old, new


rows = []
for tier, s in (("0.9", 0.92), ("0.7", 0.72), ("0.4", 0.42)):
    for g_ in range(12):
        old, new = make(s, 20260806 + g_)
        src = [x for x in paras(old) if "+0.7373" in x][0]
        NP = paras(new)
        real_ov = max(ov(src, q) for q in NP)
        best_is_unrel = max(NP, key=lambda q: ov(src, q)).startswith("另一段")
        res = judge(old, new)
        tk, kind, rep, best = res[0]
        mis = (kind in ("替换", "多候选")) and (best or "").startswith("另一段")
        rows.append(dict(tier=tier, grp=g_, real_ov=round(real_ov, 3), kind=kind,
                         rep=rep, matched_unrelated=bool(best_is_unrel), mismatch=bool(mis)))
T = pd.DataFrame(rows)
print("=== G3:三档 × 12 组 ===")
for tier in ("0.9", "0.7", "0.4"):
    s = T[T.tier == tier]
    print(f"  标称 {tier} 档 · 实测词重合中位 {s.real_ov.median():.3f} · "
          f"kind 分布 {dict(s.kind.value_counts())} · **错配 {int(s.mismatch.sum())}/{len(s)}**")

# g=0:两段都只改措辞、无量变动
old0 = f"{BODY},其中一个量是 +0.7373。\n\n{UNREL}。"
new0 = f"{BODY}(措辞略改),其中一个量是 +0.7373。\n\n{UNREL}(措辞略改)。"
g0 = judge(old0, new0)
print(f"\n  g=0(只改措辞、量不动)-> 报出 **{len(g0)}** 行(须 0)")

G = Gate("那个配对启发式,错配率是多少?")
pos_hi = int(T[T.tier.isin(("0.9", "0.7"))].mismatch.sum())
pos_lo = int(T[T.tier == "0.4"].mismatch.sum())
pla_ok = G.negative_control("g=0:只改措辞不得报出任何行", null=float(len(g0)), effect=12.0,
                            null_spread=0.4, null_kind="两段都只改措辞、没有量变动")
# 正对照:这个陷阱必须**能**被踩中 —— 造一个必然错配的极端组
old_x = f"{BODY},其中一个量是 +0.7373。\n\n{UNREL}。"
new_x = f"完全不同的文字甲乙丙丁。\n\n{UNREL},而这里出现了 +0.9999。"
x = judge(old_x, new_x)
trap_works = any(k in ("替换", "多候选") and (b or "").startswith("另一段") for _, k, _, b in x)
print(f"  正对照(把改写段整段换掉,只剩无关段可配)-> 踩中陷阱?**{trap_works}**")
pos_ok = G.positive_control("正对照:这个陷阱必须能被踩中", planted=float(trap_works), floor=0.0, spread=0.4)
if pos_ok and pla_ok:
    verdict = (f"0.9/0.7 档错配 **{pos_hi}** · 0.4 档错配 **{pos_lo}** ⇒ "
               + ("**合格**" if pos_hi == 0 and pos_lo == 0 else "**记为缺陷并修**"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(rows=T.to_dict("records"), mismatch_hi=pos_hi, mismatch_lo=pos_lo,
               g0_rows=len(g0), trap_works=bool(trap_works), verdict=verdict,
               deviation="预注册写「临时分支伪造」,本轮走纯文本支(`#173c`),刻意且标注",
               unchallenged=True),
          open(OUT/"mismatch_rate.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'mismatch_rate.json'}")
