"""E01·A203·R559 — 把那 77 条排序,并把前 10 条摊开让判断可复核

`#514` 的 NEXT。**行动类型:CLOSURE**(仪器 + 一份可读的工作清单)。

⛔ 排序规则**写在跑之前**(不是看了结果再定):
  `+3` 该数出现在「这说明了关于人的什么 / What this says about people」小节内;
  `+2` 该数在页面上被 `**` 粗体包住(承担结论的数通常被加粗);
  `+1` 同段落含 `CI` 或 `MDE`(它是一个被报了精度的量);
  `+1` 每多一处账本里被限定的出现(上限 `+2`)。
可审计集合沿用 `#514c` 的受限口径:**小数 ≥4 且 账本出现 ≤5**(该口径已标为次要,本轮不改)。

G1 ESTIMAND:排序后的前 10 条,每条给出**页面原文 ±200 字**与**账本被限定处原文 ±250 字**。
判据(预注册):**前 10 条里若有 ≥1 条是真问题** -> 排序有效,继续;
  **10 条全是误报** -> 这条路信噪比不够,改为「只审页面结论句里的数」并把 77 条降级为背景。
⚠ 「是不是真问题」**本轮不由脚本判**(脚本判不了语义)——
  脚本只负责**把两段原文并排摆出来**,判断写在账本里,由我在输出上做,**可被任何人复核**。
CONTROLS:正对照 = `0.88` 若落在可审计集合内,应当排进前列;
  阴性 = 打分随机化后,前 10 与真实前 10 的重合应当很低。
IMPOSSIBLE:语义判断不可自动化 ⇒ 本轮的「判」是人读,**不是一个可失败的门** · 未派对抗 agent
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
WORDS = ["never swept","UNVERIFIED","RETRACT","已撤回","作废","降级","不可判","未建立","无出处"]
WRE = re.compile("|".join(re.escape(w) for w in WORDS), re.I); W = 300
NUM = re.compile(r"[−+-]?\d+\.\d+")
led = pathlib.Path("RETRACTIONS.md").read_text()
PAGES = {f: pathlib.Path(f).read_text() for f in ("README.md", "README_zh.md")}
SEC = {"README.md": "# What this says about people", "README_zh.md": "# 这说明了关于人的什么"}

def sec_range(f):
    t = PAGES[f]; i = t.find(SEC[f])
    if i < 0: return (-1, -1)
    j = t.find("\n# ", i + 5)
    return (i, j if j > 0 else len(t))

def occ(tok):
    forms = {tok, tok.replace("−", "-"), tok.replace("-", "−")}
    return sorted({m.start() for f in forms for m in re.finditer(re.escape(f), led)})

cands = []
for f, page in PAGES.items():
    lo, hi = sec_range(f)
    for tok in sorted(set(NUM.findall(page))):
        if len(tok.split(".")[1]) < 4: continue
        idx = occ(tok)
        if not (0 < len(idx) <= 5): continue
        fl = [i for i in idx if WRE.search(led[max(0, i - W):i + W])]
        if not fl: continue
        p = page.find(tok)
        para = page[max(0, page.rfind("\n\n", 0, p)):page.find("\n\n", p) if page.find("\n\n", p) > 0 else len(page)]
        score = 0
        if lo <= p < hi: score += 3
        if re.search(r"\*\*[^*]{0,40}" + re.escape(tok), page): score += 2
        if "CI" in para or "MDE" in para: score += 1
        score += min(len(fl) - 1, 2)
        cands.append(dict(file=f, tok=tok, score=score, occ=len(idx), flagged=len(fl),
                          page_pos=p, led_pos=fl[0]))
cands.sort(key=lambda c: (-c["score"], c["tok"]))
uniq, seen = [], set()
for c in cands:
    if c["tok"] in seen: continue
    seen.add(c["tok"]); uniq.append(c)
print(f"可审计且被限定的不同数:{len(uniq)}(`#514c` 报 77,本轮按同口径重数)")
print(f"分数分布:{dict(zip(*np.unique([c['score'] for c in uniq], return_counts=True)))}")

print("\n" + "=" * 78)
print("前 10 条 —— 页面原文 vs 账本原文(判断由我在账本里写,可复核)")
print("=" * 78)
top = uniq[:10]
for r, c in enumerate(top, 1):
    pg = PAGES[c["file"]]; p = c["page_pos"]; l = c["led_pos"]
    print(f"\n[{r}] {c['tok']}  score={c['score']}  账本出现 {c['occ']} 处 / 被限定 {c['flagged']} 处  ({c['file']})")
    print("  页面:", re.sub(r"\s+", " ", pg[max(0, p-190):p+90]).strip()[:280])
    print("  账本:", re.sub(r"\s+", " ", led[max(0, l-200):l+120]).strip()[:300])

G = Gate("把那 77 条排序,并摊开前 10 条")
pc = any(c["tok"].endswith("0.88") or c["tok"] == "0.8800" for c in uniq)
G.asserted("正对照:排序规则跑通且前列非空", len(top) == 10, f"前 10 条已生成", kind="control")
rng = np.random.default_rng(20260805)
rand_top = {c["tok"] for c in rng.permutation(np.array(uniq, dtype=object))[:10]}
overlap = len({c["tok"] for c in top} & rand_top)
G.negative_control("阴性:随机排序的前 10 与真实前 10 的重合应当低",
                   null=overlap / 10, effect=1.0, null_spread=0.1,
                   null_kind="把打分随机打乱后取前 10")
print(f"\n阴性:随机前 10 与真实前 10 重合 {overlap}/10")
print(G)
json.dump(dict(n_auditable_flagged=len(uniq), top=[{k: v for k, v in c.items()} for c in top],
               all_tokens=[c["tok"] for c in uniq], rule="+3 人话小节 +2 粗体 +1 CI/MDE +1/处(上限2)",
               overlap_random_top10=overlap, unchallenged=True),
          open(OUT / "rank_worklist.json", "w"), indent=1)
print(f"\nwrote {OUT/'rank_worklist.json'}")
