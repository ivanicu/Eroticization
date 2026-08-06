"""E02·A12·R660 —— 那 29 次措辞更正里,有多少次是「限定语被删掉,而数字留了下来」?

`#623` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。

⚠ **P4 当场救了一次,而且救得很彻底 —— 这一条写在结果前面。**
   我在 `#623` 的 NEXT 里预注册的判据是「**删掉的句子含限定词且不含小数**」。
   而 `readme_ledger_audit.py` 里**已经有** `qualifiers_stripped`,它的判据是:
   **一条限定语行被删掉,而它限定的那个数字仍然活在页面上、且新家里没有任何限定语。**
   ⇒ **已有的严格更好**:我的版本会对「几乎每一次删除」开火(中文散文里 `而/但/只` 无处不在);
     它的版本要求**「幸存的裸数字」**,那才是真正的缺陷形状。
   ⇒ **本轮改用已有的那个,并如实记:我预注册的判据不如仓库里已经躺着的代码。**(`L21`/`P4`)

⚠ 而它还有一件事必须记:**这个函数从来没有被接进 `readme_gate`** —— `#562`/`#563` 的
   「造好而没被调用的检查等于没有」,第三次。

G1 ESTIMAND(先于方法):`可检测率` = 在 `#623` 认定的**措辞更正**提交里,
  `qualifiers_stripped(父提交 -> 该提交)` **开火**的比例。
CONTROLS:
  正对照:走**纯文本支** `qualifiers_stripped_texts`(它的 docstring 明写「正对照就跑这一支」——
    否则对照要穿过 git/cwd/文件名三层,任何一层出错都会让对照**静默失败**)。
    造一对:删掉一条带限定语的行,而它的数字在新文本里**裸着**活下来 -> **必须开火**。
  **g=0 两支**:①删掉一条**不含限定词**的行(数字仍活)-> **不得开火**;
    ②删掉一条**含限定词但数字没活下来**的行 -> **不得开火**。
KILL(条件式,预注册于 `#623`):
  if 正对照开火 and 两支 g=0 都不开火:
      可检测率 > 50% -> **接进 `readme_gate` 作第七条规则**
      否则 -> **不接**,并如实报可检测率
  else: UNVERIFIED
G3:29 条逐条发布(开火/不开火/管道失败)。G4:`rev_to` 取 {该提交 / 当前工作区} 两档。
IMPOSSIBLE(不写 planned):**它只抓「限定语被删而数字留下」这一种措辞缺陷** ——
  其余措辞缺陷(改写、语气、范围词替换)仍然没有单位可数 ·
  条目→提交映射靠 `git log -S`(`#623` 已记 6 次失配)· `[unchallenged]`
"""
import os, sys, pathlib, json, re, subprocess, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r'(?<![\w.])\d+\.\d{2,4}(?![\w])')
sh = lambda *a: subprocess.run(a, capture_output=True, text=True).stdout

# ── 控制先跑(纯文本支,按它自己的 docstring)──────────────────
G = Gate("那 29 次措辞更正里,有多少次是「限定语被删掉,而数字留了下来」?")
# ⚠ 第一版正对照失败,而原因在代码里,不在逻辑里:`_MAGNUM` 只认
#   **带符号的小数 / 百分数 / 倍数**(`[+−-]\s?\d\.\d{2,4}` | `%` | `×倍`),
#   而我的合成用的是**裸** `0.432` —— 它压根不在这个函数关心的数字种类里。
#   ⇒ **正对照要用它认得的那种数**。这也说明:主测那一片零,在修好对照之前不可采信(P5★)。
OLD_P = "⚠ 这个 +0.432 只在同一份问卷内部成立,不能外推。\n位置分是 +0.432。\n别的句子。"
NEW_P = "位置分是 +0.432。\n别的句子。"
pos = A.qualifiers_stripped_texts(OLD_P, NEW_P)
print(f"=== 正对照(纯文本支):删掉限定语行,数字 0.432 裸着活下来 -> 开火 **{len(pos)}** 条 ===")
if len(pos): print(f"    {pos.iloc[0].to_dict()}")
OLD_N1 = "这一句陈述 +0.432 并且到此为止。\n位置分是 +0.432。\n别的句子。"
g01 = A.qualifiers_stripped_texts(OLD_N1, "位置分是 +0.432。\n别的句子。")
OLD_N2 = "⚠ 这个 +0.999 只在同一份问卷内部成立,不能外推。\n别的句子。"
g02 = A.qualifiers_stripped_texts(OLD_N2, "别的句子。")
print(f"  g=0①:删掉**不含限定词**的行(数字仍活)-> 开火 **{len(g01)}**(须 0)")
print(f"  g=0②:删掉含限定词但**数字没活下来**的行 -> 开火 **{len(g02)}**(须 0)")
pos_ok = G.positive_control("正对照:限定语被删且数字裸着幸存,必须开火",
                            planted=float(len(pos)), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:两支都不得开火", null=float(len(g01) + len(g02)),
                            effect=float(max(len(pos), 1)), null_spread=0.4,
                            null_kind="删除不含限定词的行 / 限定的数字并未幸存")

# ── 主测:`#623` 认定的措辞更正提交 ───────────────────────────
led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
PAT = re.compile(r'两份 README 已改|两版 README|页面已改|页面.{0,6}更正|已更正.{0,10}页面|正确写法')
CAND = sorted({ent[j] for j, l in enumerate(led) if PAT.search(l) and ent[j] > 0})
rows = []
for N in CAND:
    log = sh("git", "log", "--format=%H", "-S", f"## Entry {N} ", "--", "RETRACTIONS.md").split()
    if not log: rows.append(dict(entry=N, status="找不到提交")); continue
    c = log[-1]; par = sh("git", "rev-parse", f"{c}^").strip()
    if not par: rows.append(dict(entry=N, status="无父提交")); continue
    d = sh("git", "show", c, "--", "README.md")
    rem = {n for l in d.splitlines() if l.startswith("-") and not l.startswith("---") for n in NUM.findall(l)}
    add = {n for l in d.splitlines() if l.startswith("+") and not l.startswith("+++") for n in NUM.findall(l)}
    kind = "数字更正" if (rem - add) else "措辞更正"
    try:
        D, plumb = A.qualifiers_stripped(par, c)
        ok = all(v[0] > 0 for v in plumb.values())
        rows.append(dict(entry=N, status=kind, fired=int(len(D)), plumbing_ok=bool(ok)))
    except Exception as e:
        rows.append(dict(entry=N, status=kind, fired=-1, err=type(e).__name__))
T = pd.DataFrame(rows)
W = T[T.status == "措辞更正"]
fired = W[W.fired > 0] if "fired" in W else W
rate = len(fired) / len(W) if len(W) else 0.0
print(f"\n=== G3:{len(T)} 条候选 ===")
print(T.status.value_counts().to_string())
print(f"\n**措辞更正 {len(W)} 条 · 规则开火 {len(fired)} 条 · 可检测率 **{rate*100:.1f}%**")
if len(W):
    for r in W.itertuples():
        print(f"    #{int(r.entry):4d} 开火 {getattr(r,'fired','?')} · 管道 {getattr(r,'plumbing_ok','?')}")

if pos_ok and pla_ok:
    verdict = (f"可检测率 **{rate*100:.1f}%** ⇒ "
               + ("**接进 `readme_gate` 作第七条规则**" if rate > 0.50 else "**不接**,如实报"))
    print(f"\n控制齐备 ⇒ 评判(门槛 50%)。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(table=T.to_dict("records"), n_wording=len(W), n_fired=len(fired), rate=rate,
               positive=len(pos), g0=[len(g01), len(g02)], verdict=verdict,
               note="预注册判据不如仓库里已有的 qualifiers_stripped;已改用后者", unchallenged=True),
          open(OUT/"qualifier_gate.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'qualifier_gate.json'}")
