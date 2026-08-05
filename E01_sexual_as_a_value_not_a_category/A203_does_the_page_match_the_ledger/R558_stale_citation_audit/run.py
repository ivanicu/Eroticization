"""E01·A203·R558 — 页面上那个数,账本里它旁边的话是不是「已经撤回」

`#513` 的 NEXT。**行动类型:CLOSURE**(仪器)。
`#512` 的审计只查**数在不在**,而 `#513b` 证明那不够:`0.88` 在账本里,
**但它旁边的话是「Computed at K=6 and never swept」**,页面照抄了数、没照抄限定。

⛔ **STRONGEST CONFOUND,写在跑之前:`RETRACTIONS.md` 本身就是一本撤回账本。**
  「撤回 / UNVERIFIED / 降级」这些词在里面**到处都是** -> 词表的**基率**可能接近 1,
  那样这个检查就是**空的**(realstat 的「不会失败的检查」)。
  ⇒ **必须同时量基率**:在账本里随机取同样宽度的窗口,看词表命中率。
     **基率 ≥ 0.8 -> 本仪器无分辨力,当场判 UNVERIFIED,不报清单。**

G1 ESTIMAND:对页面上每一个在账本里出现过的数,取它在账本里**每一处**的 ±300 字窗口,
  查是否命中**预先声明的强限定词表**;**任一处命中即入清单**,并报该数的总出现次数。
词表(写在跑之前,只收**限定一个数**的强标记,不收泛泛的叙述词):
  `never swept` · `UNVERIFIED` · `RETRACT` · `已撤回` · `作废` · `降级` · `不可判` · `未建立` · `无出处`
CONTROLS:
  正对照 **`0.88` 必须被命中**(`#513b` 已知它的账本语境是 `never swept`)—— 不中则词表没用;
  阴性   **基率**:随机窗口的命中率;它同时是「这个检查会不会到处开火」的度量。
KILL(条件式,预注册):
  if 正对照命中 and 基率 < 0.8:
      清单为空 -> 页面无过期引用;非空 -> **逐条列出,它们是下一批**
  else: UNVERIFIED
IMPOSSIBLE:一个数在账本里出现多次时,**无法判定页面引的是哪一处** ->
  本仪器采「任一处被限定即入清单」的**保守**口径,因而**会高报** · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
WORDS = ["never swept", "UNVERIFIED", "RETRACT", "已撤回", "作废", "降级", "不可判", "未建立", "无出处"]
W = 300
NUM = re.compile(r"[−+-]?\d+\.\d+")
led = pathlib.Path("RETRACTIONS.md").read_text()
WRE = re.compile("|".join(re.escape(w) for w in WORDS), re.I)

def ctx_flagged(tok):
    """返回 (出现次数, 被限定的处数)。tok 用页面上的原样字串,并做 −/- 双形匹配。"""
    forms = {tok, tok.replace("−", "-"), tok.replace("-", "−")}
    idx = []
    for f in forms:
        idx += [m.start() for m in re.finditer(re.escape(f), led)]
    idx = sorted(set(idx))
    flagged = sum(1 for i in idx if WRE.search(led[max(0, i - W):i + W]))
    return len(idx), flagged

# 阴性:基率
rng = np.random.default_rng(20260805)
pos = rng.integers(W, len(led) - W, 2000)
base = float(np.mean([bool(WRE.search(led[i - W:i + W])) for i in pos]))
print(f"账本长度 {len(led):,};词表 {len(WORDS)} 个")
print(f"⛔ 阴性/基率:随机 ±{W} 字窗口的命中率 = **{base:.3f}**  "
      f"{'-> 有分辨力' if base < 0.8 else '-> ⛔ 无分辨力,本仪器空转'}")

n88, f88 = ctx_flagged("0.88")
print(f"正对照 `0.88`:账本出现 {n88} 处,其中被限定 {f88} 处 -> "
      f"{'命中 ✅' if f88 > 0 else '未命中 ⛔ 词表没用'}")

rows, stale = [], []
for f in ("README.md", "README_zh.md"):
    page = pathlib.Path(f).read_text()
    toks = sorted({t for t in NUM.findall(page)})
    for t in toks:
        n, fl = ctx_flagged(t)
        if n == 0: continue
        rows.append(dict(file=f, tok=t, occ=n, flagged=fl))
        if fl > 0: stale.append((f, t, n, fl))
uniq = {t for _, t, _, _ in stale}
print(f"\n页面上出现在账本里的**不同数** {len({r['tok'] for r in rows})} 个;"
      f"**其中至少有一处旁边带限定词的:{len(uniq)} 个**")
top = sorted({(t, n, fl) for _, t, n, fl in stale}, key=lambda x: -x[2])[:12]
for t, n, fl in top: print(f"   {t:>10s}  出现 {n:3d} 处,被限定 {fl:3d} 处")

G = Gate("页面上那个数,账本里它旁边的话是不是已经撤回")
G.asserted("正对照:`0.88` 必须被命中", f88 > 0, f"{n88} 处出现,{f88} 处被限定", kind="control")
G.negative_control("阴性:随机窗口的基率(它决定这个检查会不会到处开火)",
                   null=base, effect=len(uniq) / max(len({r['tok'] for r in rows}), 1),
                   null_spread=0.02, null_kind="账本内随机 ±300 字窗口")
print("\n" + "=" * 70)
if f88 > 0 and base < 0.8:
    verdict = (f"**{len(uniq)} 个数在账本里至少有一处被限定 -> 逐条处理是下一批**"
               if uniq else "页面无过期引用")
    print(f"控制齐备 ⇒ 评判。基率 {base:.3f} · {verdict}")
    print("⚠ 通过的 KILL 会怎样失败:一个数在账本里出现多次时,**无法判定页面引的是哪一处** ——"
          "本仪器采「任一处被限定即入清单」的保守口径,**因而必然高报**;"
          f"而基率 {base:.3f} 意味着**随机一段话里就有这么高的概率出现限定词**,"
          "所以清单里的多数条目**很可能只是撞上了账本的叙述风格**。")
else:
    verdict = f"UNVERIFIED —— 正对照={f88>0} 基率={base:.3f}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(words=WORDS, window=W, base_rate=base, pc_0_88=dict(occ=n88, flagged=f88),
               n_distinct=len({r["tok"] for r in rows}), n_stale=len(uniq),
               stale=sorted(uniq), verdict=verdict, unchallenged=True),
          open(OUT / "stale_citation_audit.json", "w"), indent=1)
print(f"\nwrote {OUT/'stale_citation_audit.json'}")
