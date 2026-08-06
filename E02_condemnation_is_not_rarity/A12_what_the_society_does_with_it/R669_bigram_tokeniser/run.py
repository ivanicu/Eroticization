"""E02·A12·R669 —— 换成字符 2-gram 之后,中文那一半真的看见了吗?

`#632` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。`#620` 的回测是硬门槛。

⚠ **§3 梯度检查:2-gram 会让「什么都相似」,所以 g=0 必须是真正无关的两段中文。**
   一个把阈值降到人人都过的分词,和一个谁都不过的分词,同样没用 ——
   **只报正对照升上去了,是半个结果。**

G1 ESTIMAND(先于方法,三个):
  E1 `ov` 在 `#632d` 那两段中文上的值:**旧 0.000 -> 新必须 > 0.50**;
  E2 g=0:两段**真正无关**的中文,新 `ov` 必须 **< 0.50**;
  E3 **中文页 `UNMATCHED` 比例**在 32 个历史版本上的**前后变化** ——
     **不降就如实报**,不许只报 E1。
分词(先写死):**CJK 字符取 2-gram;非 CJK 仍按 `\\w+` 取词;两者取并集。**
KILL(条件式,预注册于 `#632`):
  if E1 > 0.50 and E2 < 0.50:
      **回测**:32 个历史版本上 `numbers_that_left` 的 `token` 列**一个不少** -> 接受
      少一个 -> **不改**
  else: UNVERIFIED
G3:三个估计量 + 回测表。G4:2-gram / 3-gram 两档。
IMPOSSIBLE(不写 planned):**2-gram 是形状相似,不是语义相似** ——
  两段讲不同事情但用词高度重叠的中文,它仍会判成相似 · `[unchallenged]`
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
CJK = re.compile(r'[一-鿿]')
run = lambda *a: subprocess.run(a, capture_output=True, text=True)


def toks_old(s): return set(re.findall(r'\w+', s))


def toks_new(s, n=2):
    """CJK 取 n-gram;非 CJK 仍按 \\w+ 取词;并集。"""
    cj = "".join(CJK.findall(s))
    grams = {cj[i:i+n] for i in range(max(len(cj)-n+1, 0))}
    lat = set(re.findall(r'[A-Za-z0-9_+\-.]+', s))
    return grams | lat


def ov(a, b, f): 
    A_, B_ = f(a), f(b)
    return len(A_ & B_) / max(len(A_), 1)


# ── E1 / E2 ─────────────────────────────────────────────────
P1 = "这一段讨论社会层的谴责与它的对象句子较长以便词重合度可控"
P2 = "这一段讨论个人层的偏好与它的来源句子较短以便词重合度可控"
U1 = "另一段完全无关的文字谈的是农业耕作与季节安排跟上面那一段没有共同的实质内容"
print("=== E1 / E2:分词换成 2-gram 之后 ===")
for nm, a, b in (("E1 相似的两段中文", P1, P2), ("E2 真正无关的两段中文", P1, U1)):
    o_, n_ = ov(a, b, toks_old), ov(a, b, toks_new)
    print(f"  {nm:20s} 旧 `\\w+` **{o_:.3f}** -> 新 2-gram **{n_:.3f}**")
E1 = ov(P1, P2, toks_new); E2 = ov(P1, U1, toks_new)
print(f"  判据:E1 > 0.50 ? **{E1>0.50}** · E2 < 0.50 ? **{E2<0.50}**")

# ── E3:32 个历史版本上的 UNMATCHED 比例(新旧两版分词)────────
norm = lambda m: m.strip().replace(' ', '').replace('倍', '×').replace('−', '-')
paras = lambda s: [x for x in re.split(r'\n\s*\n', s) if x.strip()]


def judge(old, new, f):
    o = {norm(m) for m in A._MAGNUM.findall(old)}; n = {norm(m) for m in A._MAGNUM.findall(new)}
    NP = paras(new); newc = n - o; out = []
    for tk in sorted(o - n):
        src = [x for x in paras(old) if tk in norm(x) or tk in x]
        kind = "UNMATCHED"
        if src and NP:
            best = max(NP, key=lambda q: ov(src[0], q, f))
            if ov(src[0], best, f) >= 0.50:
                cand = sorted({norm(m) for m in A._MAGNUM.findall(best) if norm(m) in newc})
                kind = "拿走" if not cand else ("替换" if len(cand) == 1 else "多候选")
        out.append((tk, kind))
    return out


led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
PAT = re.compile(r'两份 README 已改|两版 README|页面已改|页面.{0,6}更正|已更正.{0,10}页面|正确写法')
CAND = sorted({ent[j] for j, l in enumerate(led) if PAT.search(l) and ent[j] > 0})
stat = {("old", f): {"n": 0, "unm": 0} for f in ("README.md", "README_zh.md")}
stat.update({("new", f): {"n": 0, "unm": 0} for f in ("README.md", "README_zh.md")})
miss = 0; nver = 0
for N in CAND:
    log = run("git", "log", "--format=%H", "-S", f"## Entry {N} ", "--", "RETRACTIONS.md").stdout.split()
    if not log: continue
    par = run("git", "rev-parse", f"{log[-1]}^").stdout.strip()
    if not par: continue
    nver += 1
    for f in ("README.md", "README_zh.md"):
        r = run("git", "show", f"{par}:{f}")
        if r.returncode != 0: continue
        old, new = r.stdout, pathlib.Path(f).read_text()
        ro = judge(old, new, toks_old); rn = judge(old, new, toks_new)
        if {t for t, _ in ro} - {t for t, _ in rn}: miss += 1
        for tag, res in (("old", ro), ("new", rn)):
            stat[(tag, f)]["n"] += len(res)
            stat[(tag, f)]["unm"] += sum(1 for _, k in res if k == "UNMATCHED")
print(f"\n=== E3 · 回测:{nver} 个历史版本 · **token 列丢失 {miss} 次** ===")
rows = []
for f in ("README.md", "README_zh.md"):
    o_, n_ = stat[("old", f)], stat[("new", f)]
    ro = o_["unm"]/o_["n"] if o_["n"] else float("nan")
    rn = n_["unm"]/n_["n"] if n_["n"] else float("nan")
    rows.append(dict(file=f, n=o_["n"], unmatched_old=o_["unm"], unmatched_new=n_["unm"],
                     rate_old=ro, rate_new=rn))
    print(f"  {f:14s} 共 {o_['n']:4d} 行 · UNMATCHED 旧 {o_['unm']:4d}({ro*100:5.1f}%) "
          f"-> 新 {n_['unm']:4d}({rn*100:5.1f}%)")
T = pd.DataFrame(rows)

G = Gate("换成字符 2-gram 之后,中文那一半真的看见了吗?")
pos_ok = G.positive_control("正对照:相似的两段中文 ov 必须 > 0.50", planted=float(E1), floor=0.50, spread=0.02)
pla_ok = G.negative_control("g=0:真正无关的两段中文 ov 必须 < 0.50", null=float(E2), effect=float(E1),
                            null_spread=0.02, null_kind="两段讲完全不同题目的中文")
zh = T[T.file == "README_zh.md"].iloc[0]
dropped = zh.rate_new < zh.rate_old
if pos_ok and pla_ok:
    verdict = (f"回测丢失 {miss} ⇒ " + ("**接受并接入**" if miss == 0 else "**不改**")
               + f";中文页 UNMATCHED {zh.rate_old*100:.1f}% -> {zh.rate_new*100:.1f}%"
               + ("(**下降**)" if dropped else "(**未下降 —— 换分词没解决问题,如实报**)"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4:2-gram vs 3-gram ===")
for n in (2, 3):
    f = lambda s, _n=n: toks_new(s, _n)
    print(f"  {n}-gram: E1 {ov(P1,P2,f):.3f} · E2 {ov(P1,U1,f):.3f}")
json.dump(dict(E1=E1, E2=E2, backtest_versions=nver, token_loss=miss,
               unmatched=T.to_dict("records"), verdict=verdict, unchallenged=True),
          open(OUT/"bigram.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'bigram.json'}")
