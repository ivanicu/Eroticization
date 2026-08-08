"""#783 · E03·A43·R222 —— 页面上那些「虔诚/非虔诚」比,有几个是裸的点估计?

`#782`① 预注册的第一句是 **「先数」**:
  「⚠ 按 `#623` 先回测:**先数出页面上有多少个比值型点估计**,再决定是逐个挂还是整体改成
   『只在 0/30 那四个世代上报点值,其余只报区间』。」
⇒ **不许先定规则再数。** 本轮先数,再按数出来的东西定规则。

G1 估计量(方法之前先命名,两个量):
  (a) **N_裸** = 页面(中英两版)上发表的「虔诚/非虔诚」比里,**不带自己的不确定区间**的个数
  (b) 对每一个,它的**年份自助 95% 区间**,以及**区间是否含 1.0**

识别:
  (a) 是对一个文本语料的**普查** —— 识别没问题,但**抽取器本身是一具仪器**
      (`realstat §4`:「一个 search 就是一具仪器,而它没有正对照」),所以必须回测。
  (b) 世代格已在 `#782` 的 44 格里;**而合并版(`#778` 的 0.409–0.431 / 0.328–0.341)不在**
      —— 那正是页面上被引用得最多的两个数,本轮**当场把它们的区间算出来**。

⚠⚠ **三个世界,而第三个是元分离(`frontier §3`):**
  A **区间是装饰**:每个已发表比值的区间都窄且排除 1.0 ⇒ 挂上去不改任何句子,纯记账。
  B **页面超发**:≥1 个已发表比值的区间含 1.0,**或**根本算不出区间(合并版)
     ⇒ 承载它的那句话必须缩。
  C **我的分解本身错了**:页面上已发表的一半根本不是点估计,而是**规格区间**
     (`0.409–0.431` 是两种分层的跨度,`[0.32, 0.50]` 是两种切法的跨度)——
     **那是另一种不确定性**,给它挂自助区间是类型错误。
     ⇒ C 由一个 A/B 都产生不了的观测量分离:**同一批格的「自助宽度」与「规格跨度」谁大。**

预测矩阵(粗数,形状才是重点):
  | 世界 | 现在 | 若区间都窄且排除 1.0 | 若 ≥1 含 1.0 或算不出 | 若自助宽 ≫ 规格跨度 |
  | A   | 0.35 | 0.85 | 0.05 | 0.10 |
  | B   | 0.45 | 0.05 | 0.85 | 0.30 |
  | C   | 0.20 | 0.10 | 0.10 | **0.60** |

⚠ 跑之前写下的最强混淆:**抽取器会把「不是这个量」的小数一起捞进来**
  —— `#781` 那行的自由度杠杆 `0.147 · 0.146 · 0.075 · 0.022`、`#780` 的 `1.49×/1.74×`
  都是三位小数,**长得和比值一模一样**。⇒ 同一轮里就放对照:
  **手标名单在抽取器跑之前写死在下面**,回测同时报 **recall 与 precision**,
  假阳性逐个列出来 —— `#623` 要的正是这个,而不是只报 recall。
⚠ 手标名单是我自己造的 ⇒ **它是我自己的仪器**,所以**有信息的是两者的差**:
  机器找到而我没标的 = 我漏了;我标了而机器找不到的 = 抽取器盲。两个方向都报。

预注册判词(条件式,不是阈值 —— `#P16`):
  if 抽取器正对照开火(recall ≥ 0.8)and 自助真的在变:
      if 任一已发表比值的区间含 1.0 **或** 算不出区间:
          -> B:那句话必须改写 / 那个数必须换成区间
      elif 自助宽度中位 > 2 × 规格跨度中位:
          -> C:页面公布的是规格不确定性,**藏起了更大的抽样不确定性**,两个都要写
      else:
          -> A:纯加法,挂上去即可
  else:
      -> UNVERIFIED(抽取器或自助不合格,不许下判)
  ⚠ 2× 这个阈值不是新造的:`#778` 判「两个统计量是否一致」用的就是 2×,沿用同一个。

本轮改不了的(`realstat §2` 登记,每条附「需要什么」):
  · **合并版比值的因果方向** —— 需要面板,GSS 是重复截面。**结构性不可能。**
  · **英文页与中文页是否逐字对应** —— 需要一具对齐仪器;本轮只分别数,**不声称两页等价**。
  · **自助是对年份点重抽,不是对人** —— 需要个体层 replicate weights;GSS 有 `vpsu/vstrat`
    但仅 1975 年后,**跨全 28 年不可用**,所以年份自助是这里能做的最好的,且它**低估**了人层面的抖动。
"""
import pandas as pd, numpy as np, json, pathlib, sys, re
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate

RNG = np.random.default_rng(222)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"; OUT.mkdir(exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# ① 手标名单 —— 跑之前写死。这是 ground truth,不是抽取器的输出。
#    规则:只收「虔诚层某统计量 ÷ 非虔诚层同一统计量」这一个量,其余一律不收。
# ─────────────────────────────────────────────────────────────────────────────
# ⚠⚠ 第一版这张表没有 `unit` 一列,于是靠**数值**去对 `#782` 的 44 格 —— 而那是**单位盲**:
#    合并版的 0.328 对上了一个**世代格**(纯属数值巧合),0.226 对上了 0.222 那一格(差 0.004),
#    0.375 是**整个网格的中位数**,它根本不是一个格。⇒ `realstat §4`「一个 search 就是一具仪器」
#    那一条的补救原文:**先把「仪器的单位」与「主张的单位」写成两个字符串,并要求它们相等**,
#    再谈对照。加上 `unit` 之后,可对格的从 14 个掉到 **9 个**,而**纯数值会误配 5 个**。
HAND = {  # value -> (哪一轮, 单位, 它是什么, 是否已带区间)
    "0.409": ("#778", "pooled", "合并 · 水平比 · 分层(b)三题三分位", False),
    "0.431": ("#778", "pooled", "合并 · 水平比 · 分层(a)attend三档", False),
    "0.328": ("#778", "pooled", "合并 · 端点(相对基线)比 · 分层(b)", False),
    "0.341": ("#778", "pooled", "合并 · 端点(相对基线)比 · 分层(a)", False),
    "0.321": ("#779", "cell", "世代1930–49 · 水平比 · k3", False),
    "0.222": ("#779", "cell", "世代1930–49 · 端点相对基线比 · k3", False),
    "0.435": ("#779", "cell", "世代1950–64 · 水平比 · k3", False),
    "0.230": ("#779", "cell", "世代1950–64 · 端点相对基线比 · k3", False),
    "0.226": ("#779", "derived", "两世代端点比的合并值(⚠ 推导,不是独立测量)", False),
    "0.501": ("#780", "cell", "世代1930–49 · 水平比 · k2", False),
    "0.185": ("#781", "gridstat", "50格网格 · 比值的 5% 分位", False),
    "0.375": ("#781", "gridstat", "50格网格 · 比值的中位", False),
    "0.692": ("#781", "gridstat", "50格网格 · 比值的 95% 分位", False),
    "1.128": ("#781", "cell", "触发撤回那一格 · Seg-B k3 n≥60 1975–94 水平", True),  # #782 给了区间
    "0.438": ("#781", "cell", "同世代另一格", False),
    "0.623": ("#781", "cell", "同世代另一格", False),
    "0.821": ("#781", "cell", "同世代另一格", False),
}
# 手标的「规格跨度」形式(不是点估计,而是一段区间 —— 世界 C 的对象)
HAND_SPAN = {
    "0.409–0.431": ("#778", "水平比 · 跨两种分层"),
    "0.328–0.341": ("#778", "端点相对基线比 · 跨两种分层"),
    "[0.32, 0.50]": ("#780", "水平比 · 跨两种切法"),
    "[0.22, 0.37]": ("#780", "端点相对基线比 · 跨两种切法"),
}
# 手标的「长得像但不是这个量」—— 抽取器若捞到它们就是假阳性
DECOYS = ["0.147", "0.146", "0.075", "0.022",     # #781 自由度杠杆
          "0.181", "0.539", "1.26", "1.49", "1.74", "2.38", "18.5"]

# ─────────────────────────────────────────────────────────────────────────────
# ② 抽取器 —— 机械规则,跑在 #778–#782 的锚窗内
# ─────────────────────────────────────────────────────────────────────────────
NUM = re.compile(r"(?<![\d.×倍%])([01]\.\d{3})(?![\d]|\s*[×倍%点])")
def windows(path):
    txt = open(path, encoding="utf-8").read()
    ancs = list(re.finditer(r"\[#(\d+)「[^」]*」\]", txt))
    out = {}
    for i, m in enumerate(ancs):
        if m.group(1) in ("778", "779", "780", "781", "782"):
            out[m.group(1)] = txt[(ancs[i-1].end() if i else 0):m.end()]
    return out

print("=== ① 抽取器的回测(`#623`:先回测再用;recall 与 precision 都报)===")
found = {}
for lang, path in (("zh", ROOT/"README_zh.md"), ("en", ROOT/"README.md")):
    w = windows(path)
    hits = {}
    for r, t in w.items():
        for v in NUM.findall(t):
            hits.setdefault(v, []).append(r)
    found[lang] = hits
    print(f"  {lang} 页:锚窗 {len(w)} 个 · 抽到 {len(hits)} 个不同数值")

zh = set(found["zh"]); hand = set(HAND)
tp = zh & hand; fn = hand - zh; fp_all = zh - hand
fp_decoy = fp_all & set(DECOYS); fp_other = fp_all - set(DECOYS)
recall = len(tp)/len(hand); precision = len(tp)/max(1, len(zh))
print(f"  手标 {len(hand)} 个 · 抽到 {len(zh)} 个 · 命中 {len(tp)}")
print(f"  **recall {recall:.2f}({len(tp)}/{len(hand)}) · precision {precision:.2f}({len(tp)}/{len(zh)})**")
print(f"  抽取器漏掉(我标了它没找到,= 仪器盲):{sorted(fn) if fn else '无'}")
print(f"  抽取器多捞的诱饵(我预先列为「长得像但不是」):{sorted(fp_decoy) if fp_decoy else '无'}")
print(f"  ⚠ 两边都没预料到的(机器找到、我既没标也没列进诱饵):{sorted(fp_other) if fp_other else '无'}")

# ─────────────────────────────────────────────────────────────────────────────
# ③ 把已发表的比值对到 `#782` 的 44 格,看谁有区间、谁没有
# ─────────────────────────────────────────────────────────────────────────────
grid = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A43_两条已发表的话从没并排读过/"
                      "R221_a_real_sign_flip_or_my_criterion_was_loose/results/ratio_interval.json"))["cells"]
def match(v, tol=0.004):
    return [c for c in grid if abs(c["ratio"]-v) <= tol]

print(f"\n=== ② 已发表的比值 × `#782` 的 {len(grid)} 格:谁有区间? ===")
print("  ⚠ **单位必须相等才允许配对**:只有 unit=='cell' 的数才是一个格;"
      "pooled/gridstat/derived 三类**不是格**,数值上对得再准也是假配对。")
pub, false_pairs = [], []
for v, (rnd, unit, what, has_iv) in sorted(HAND.items()):
    ms = match(float(v))
    rec = dict(value=float(v), round=rnd, unit=unit, what=what, already_has_interval=has_iv,
               matched=len(ms) if unit == "cell" else 0, lo=None, hi=None, covers1=None)
    if ms and unit == "cell":
        c = min(ms, key=lambda c: abs(c["ratio"]-float(v)))
        rec.update(lo=c["lo"], hi=c["hi"], covers1=c["covers1"], width=c["hi"]-c["lo"])
        tag = f"[{rec['lo']:.3f}, {rec['hi']:.3f}]" + ("  ⚠含1.0" if rec["covers1"] else "")
    elif ms:
        c = min(ms, key=lambda c: abs(c["ratio"]-float(v)))
        false_pairs.append(dict(value=float(v), unit=unit, would_match=c["cohort"], stat=c["stat"]))
        tag = f"**单位不符 ⇒ 无区间**(数值上会误配到 {c['cohort']}·{c['stat']})"
    else:
        tag = "**无区间可挂**"
    pub.append(rec)
    print(f"  {v}  {rnd}  {unit:8s} {what[:34]:36s} {tag}")

no_iv = [p for p in pub if p["matched"] == 0]
cov = [p for p in pub if p["covers1"]]
print(f"\n  已发表 {len(pub)} 个 · **单位对得上、因而有区间的 {len(pub)-len(no_iv)} 个** · "
      f"**没有区间的 {len(no_iv)} 个**")
print(f"  ⚠ 其中 **{len(false_pairs)} 个是纯数值配对会误配的**(第一版就是这么配的):"
      f"{[f['value'] for f in false_pairs]}")
print(f"  区间含 1.0 的:{len(cov)} 个 —— {[p['value'] for p in cov] if cov else '无'}")

# ─────────────────────────────────────────────────────────────────────────────
# ④ 把对不上的合并版比值,当场算出它自己的区间(#778 的四个数)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n=== ③ 合并版(`#778`)的四个数从没有过区间 —— 当场算 ===")
VALID = {"homosex": (1, 4), "attend": (0, 8), "reliten": (1, 4), "fund": (1, 3)}
for c, rg in VALID.items():
    dr, _ = check_kept_codes(ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta", c, rg)
    if dr: print(f"  #766 前瞻:{c} 删 " + " · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a, b, n, sh in dr[:2]))
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
d = pd.read_stata(gp, columns=["year"]+list(VALID), convert_categoricals=False)
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, lo=VALID[c][0], hi=VALID[c][1]: (v >= lo) & (v <= hi)) for c in VALID})
M["year"] = d.year
cat = pd.read_stata(gp, columns=["homosex"], convert_categoricals=True)
for c in aligned({"homosex": list(cat["homosex"].cat.categories)[:4]}, "strict"): M[c] = -M[c]+5
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
z = lambda s: (s-s.mean())/s.std(ddof=1)
sub = M.dropna(subset=["homosex", "attend", "reliten", "fund", "year"]).copy()
sub["REL"] = z(sub[["attend", "reliten", "fund"]]).mean(axis=1)
sub["a"] = pd.cut(sub.attend, [-1, 1, 5, 8], labels=[0, 1, 2]).astype(float)
sub["b"] = sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def series(col, kk, nmin=120):
    rows = []
    for y, gy in sub[sub[col] == kk].groupby("year"):
        if len(gy) < nmin: continue
        rows.append((int(y), float(gy.homosex.mean()), float((gy.homosex == 4).mean())))
    return rows
def ratio_boot(rowsA, rowsB, stat, B=4000):
    yA = np.array([r[0] for r in rowsA], float); yB = np.array([r[0] for r in rowsB], float)
    j = 1 if stat == "水平" else 2
    vA = np.array([r[j] for r in rowsA]); vB = np.array([r[j] for r in rowsB])
    f0A, f0B = rowsA[0][2], rowsB[0][2]
    def rat(ia, ib):
        sa, sb = slope(yA[ia], vA[ia]), slope(yB[ib], vB[ib])
        if stat != "水平": sa, sb = sa/f0A, sb/f0B
        return sa/sb if abs(sb) > 1e-12 else np.nan
    obs = rat(np.arange(len(yA)), np.arange(len(yB)))
    bs = np.array([r for r in (rat(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB)))
                               for _ in range(B)) if np.isfinite(r)])
    return obs, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), len(bs)

pooled = []
for cname, col in (("(a) attend三档", "a"), ("(b) 三题三分位", "b")):
    rA, rB = series(col, 2), series(col, 0)
    for st in ("水平", "端点相对基线"):
        o, l, h, nb = ratio_boot(rA, rB, st)
        pooled.append(dict(cut=cname, stat=st, ratio=o, lo=l, hi=h, width=h-l,
                           covers1=bool(l <= 1.0 <= h), nA=len(rA), nB=len(rB), nboot=nb))
        print(f"  {cname:14s} {st:8s} 比值 {o:.3f} · **95% 区间 [{l:.3f}, {h:.3f}]** 宽 {h-l:.3f}"
              f" · 年 {len(rA)}/{len(rB)}{'  ⚠含1.0' if l <= 1.0 <= h else ''}")
PLO, PHI = min(p["lo"] for p in pooled), max(p["hi"] for p in pooled)
PSPAN = (min(p["ratio"] for p in pooled), max(p["ratio"] for p in pooled))
print(f"  ⇒ 四个合并数的**点估计跨度 [{PSPAN[0]:.3f}, {PSPAN[1]:.3f}]**(页面写的「三分之一到五分之二」),"
      f"而**把各自的抽样不确定性并进来是 [{PLO:.3f}, {PHI:.3f}]**")
print(f"    —— 宽了 {(PHI-PLO)/(PSPAN[1]-PSPAN[0]):.1f} 倍。**页面那句话的精度紧了这么多。**")

# ─────────────────────────────────────────────────────────────────────────────
# ⑤ 世界 C 的分离量:自助宽度 vs 规格跨度
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n=== ④ 世界 C 的分离量:页面公布的「规格跨度」对上「自助宽度」 ===")
SPANS = {"0.409–0.431": ("合并·水平", [p for p in pooled if p["stat"] == "水平"]),
         "0.328–0.341": ("合并·端点相对基线", [p for p in pooled if p["stat"] != "水平"]),
         "[0.32, 0.50]": ("世代1930–49·水平·跨切法", [c for c in grid if c["cohort"] == "1930–1949" and c["stat"] == "水平"]),
         "[0.22, 0.37]": ("世代1930–49·端点·跨切法", [c for c in grid if c["cohort"] == "1930–1949" and c["stat"] != "水平"])}
cmp_rows = []
for span, (label, cells) in SPANS.items():
    rs = [c["ratio"] for c in cells]
    span_w = max(rs)-min(rs)
    boot_w = float(np.median([c["width"] if "width" in c else c["hi"]-c["lo"] for c in cells]))
    cmp_rows.append(dict(span=span, label=label, n_cells=len(cells), span_width=span_w,
                         boot_width_median=boot_w, ratio=boot_w/span_w if span_w > 1e-9 else np.inf))
    print(f"  {span:14s} {label:24s} 规格跨度 {span_w:.3f} · 自助宽度中位 {boot_w:.3f}"
          f" · **自助/规格 = {boot_w/span_w:.1f}×**")
med_ratio = float(np.median([r["ratio"] for r in cmp_rows]))
print(f"  ⇒ 四处的中位:**自助宽度是规格跨度的 {med_ratio:.1f} 倍**(预注册阈值 2×)")

# ─────────────────────────────────────────────────────────────────────────────
# ⑥ 闸 + 判词(条件式,不是阈值)
# ─────────────────────────────────────────────────────────────────────────────
G = Gate("#783 · 页面上那些比值有几个是裸的点估计")
G.asserted("① 正控:抽取器必须能看见已知在页上的数(recall ≥ 0.8)",
           bool(recall >= 0.8), f"recall {recall:.2f}({len(tp)}/{len(hand)}) · precision {precision:.2f}", kind="control")
G.asserted("② 正控:自助真的在变(每个合并格的有效重抽 >1000 且区间非零宽)",
           bool(all(p["nboot"] > 1000 and p["width"] > 1e-6 for p in pooled)),
           f"有效重抽 {[p['nboot'] for p in pooled]} · 宽 {[round(p['width'],3) for p in pooled]}", kind="control")
G.asserted("③ 负控:抽取器必须**不是**把窗内所有三位小数都捞进来(诱饵不得全中)",
           bool(len(fp_decoy) < len(DECOYS)),
           f"预列诱饵 {len(DECOYS)} 个 · 抽取器捞到 {len(fp_decoy)} 个", kind="control")
G.asserted("④ 负控:单位检查必须真的在挡东西(否则它是装饰;第一版没有它,当场误配)",
           bool(len(false_pairs) > 0),
           f"纯数值会误配 {len(false_pairs)} 个:{[f['value'] for f in false_pairs]}", kind="control")
ctrl_ok = bool(recall >= 0.8 and all(p["nboot"] > 1000 for p in pooled)
               and len(fp_decoy) < len(DECOYS) and len(false_pairs) > 0)
B_fires = bool(len(cov) > 0 or len(no_iv) > 0)
C_fires = bool(med_ratio > 2.0)
# ⚠ 判词行按本库的约定写:**PASS = 页面现在的写法站得住**,FAIL = 被推翻(`#782` 同一写法)
G.asserted("⑤ kill(预注册):页面现写法要站住,需每个已发表比值都能对上格**且**区间排除 1.0",
           not B_fires, f"含1.0 {len(cov)} 个 · 无区间可挂 {len(no_iv)} 个", kind="kill")
G.asserted("⑥ kill(预注册):页面公布的「规格跨度」要能代表不确定性,需自助宽度中位 ≤ 2× 规格跨度",
           not C_fires, f"自助/规格中位 {med_ratio:.1f}× vs 阈值 2×", kind="kill")
print(); print(G)

print("\n" + "="*92)
if not ctrl_ok:
    verdict = "**UNVERIFIED:抽取器或自助不合格,本轮不下判。**"
else:
    parts = []
    if B_fires:
        parts.append(f"**B 开火**:页面发表 {len(pub)} 个「虔诚/非虔诚」比,**{len(no_iv)} 个对不上任何已算过的格**"
                     f"(合并版 `#778` 的四个数 + 网格分位数),**{len(cov)} 个的区间含 1.0**。"
                     f"⇒ 承载它们的句子必须挂区间或改写。")
    if C_fires:
        parts.append(f"**C 也开火,而它比 B 更重**:页面上已发表的四处「区间」全是**规格跨度**,"
                     f"而同一批格的**自助宽度是它的 {med_ratio:.1f} 倍** —— "
                     f"⇒ **页面一直在公布小的那种不确定性,而把大的那种留在了幕后。**"
                     f"我的 A/B 分解假定「点估计缺区间」,而真正的问题是**已发表的区间是另一种东西。**")
    if not parts:
        parts.append("**A**:每个已发表比值都有窄区间且排除 1.0 ⇒ 纯加法,挂上即可。")
    verdict = "\n  ".join(parts)
print(verdict)

json.dump(dict(recall=recall, precision=precision, hand=len(hand), extracted=len(zh),
               false_pairs=false_pairs, pooled_point_span=list(PSPAN), pooled_union=[PLO, PHI],
               missed=sorted(fn), decoys_caught=sorted(fp_decoy), unexpected=sorted(fp_other),
               published=pub, n_no_interval=len(no_iv), n_covers1=len(cov),
               pooled=pooled, span_vs_boot=cmp_rows, med_boot_over_span=med_ratio,
               B_fires=B_fires, C_fires=C_fires, verdict=verdict,
               gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"naked_point_estimates.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'naked_point_estimates.json'}")
