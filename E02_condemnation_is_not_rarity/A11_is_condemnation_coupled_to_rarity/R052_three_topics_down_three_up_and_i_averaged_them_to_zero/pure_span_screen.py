"""E02·A198·R542 — 先只问「量表挪得动吗」,一个 lnOR 都不算

`#496` 的 NEXT。**行动类型:PRODUCTION**(造一份合格话题清单,不分离任何世界)。
⛔ **本轮不计算任何 `lnOR`** —— 只算三个切点上的**谴责占比跨度**。
   因此它**不产生关于「可见痕迹」猜想的任何证据**,**不消耗样本外点**。脚本末尾对此有 assert。

`#496c` 的前置判据(已冻结):阈值斜率检验要求**谴责占比跨度 ≥ 0.30**。
  GSS `homosex` 只有 0.10–0.12(分布集中在两端)-> 结构上不合格;
  NSFG 那六个话题是 0.55–0.78 -> 合格。

G1 ESTIMAND:对每一道态度题,`span = max(share) − min(share)`,
  share 取三个嵌套切点(最严 / 中 / 最宽)上的谴责者占比。**只此一个量。**

预注册的决策(写在跑之前):
  **≥3 个新话题** 同时满足 (a) `span ≥ 0.30` (b) 在同一份问卷里配得上一个行为变量
  -> 猜想**可被检验**,下一轮做多点样本外;
  **< 3 个** -> **直说「这个猜想在现有仪器上无法被检验」**,写进页面「做不到什么」,留在冻结文件里。

⚠ 这不是一个 kill(没有世界被分离),所以**没有条件式 kill**,只有一个**预注册的决策规则**。
  诚实标注:把筛选叫成检验,就是 realstat 的「验证串不是计算」那一类。
IMPOSSIBLE:span 只说量表挪不挪得动,**不说斜率是否存在** · 行为变量的配对仍需人判断,
  而我已看过六条斜率 ⇒ 配对**只列候选,不下结论** · 未派对抗 agent ⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

SPAN_MIN = 0.30
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NS = ROOT / "data/external/nsfg"
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
rows = []


def span_of(v, levels):
    """v: 原始码;levels: 有序的可用等级(从「最谴责」到「最不谴责」)。
    三个嵌套切点 = 前 1 / 前 2 / 前 3 个等级算谴责。"""
    v = v[np.isin(v, levels)]
    if len(v) < 500 or len(levels) < 4: return None
    sh = [float(np.isin(v, levels[:k]).mean()) for k in (1, 2, 3)]
    return dict(shares=[round(x, 4) for x in sh], span=round(max(sh) - min(sh), 4), n=int(len(v)))


# ---------------------------------------------------------------- GSS
print("=== GSS 态度题筛选(规则①:逐题打印码与 n)===")
GSS_ITEMS = {  # 变量 -> 从「最谴责」到「最不谴责」的等级序
    "homosex": [1, 2, 3, 4], "premarsx": [1, 2, 3, 4], "xmarsex": [1, 2, 3, 4],
    "teensex": [1, 2, 3, 4], "pornlaw": [1, 2, 3], "sexeduc": [1, 2, 3],
    "abany": [1, 2], "divlaw": [1, 2, 3], "spanking": [1, 2, 3, 4],
    "suicide1": [1, 2], "grass": [1, 2], "cappun": [1, 2],
    "letdie1": [1, 2], "fepol": [1, 2],
}
it = pd.read_stata(DTA, iterator=True); vl = it.variable_labels()
have = [c for c in GSS_ITEMS if c in vl]
g = pd.read_stata(DTA, columns=have, convert_categoricals=False)
for c in have:
    s = span_of(g[c].dropna().values, GSS_ITEMS[c])
    if s is None:
        print(f"  {c:10s} 等级<4 或 n<500 -> 不适用(需要至少 4 级才有三个嵌套切点)")
        rows.append(dict(src="GSS", var=c, label=str(vl[c])[:52], span=None, note="等级不足"))
        continue
    rows.append(dict(src="GSS", var=c, label=str(vl[c])[:52], **s))
    print(f"  {c:10s} n={s['n']:6d} shares={s['shares']} **span={s['span']:.4f}** "
          f"{'✅' if s['span'] >= SPAN_MIN else '⛔'}  {str(vl[c])[:40]}")

# ---------------------------------------------------------------- NSFG
print("\n=== NSFG 2011-2013 IH 段筛选 ===")
def parse_dct(p):
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out


LAY = parse_dct(NS / "setup" / "2011_2013_FemRespSetup.dct")
NIT = ["staytog", "samesex", "sxok18", "sxok16", "chunless", "chsuppor",
       "gayadopt", "okcohab", "marrfail", "chcohab", "prvntdiv"]
cols = {n: LAY[n] for n in NIT if n in LAY}
buf = {n: [] for n in cols}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n, (s, w, _) in cols.items():
        v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
REV = {"okcohab", "chunless", "marrfail", "prvntdiv"}   # 禁止式/反向题干
for n in cols:
    v = np.array(buf[n]); v = v[np.isfinite(v)]
    lv = [1, 2, 3, 4] if n in REV else [4, 3, 2, 1]     # 从最谴责到最不谴责
    s = span_of(v, lv)
    if s is None: continue
    rows.append(dict(src="NSFG", var=n, label=cols[n][2][:52], **s))
    print(f"  {n:10s} n={s['n']:6d} shares={s['shares']} **span={s['span']:.4f}** "
          f"{'✅' if s['span'] >= SPAN_MIN else '⛔'}  {cols[n][2][:40]}")

# ---------------------------------------------------------------- 决策
ok = [r for r in rows if r.get("span") is not None and r["span"] >= SPAN_MIN]
used = {"samesex", "sxok18", "sxok16", "staytog", "okcohab", "chsuppor"}   # `#494a` 已算过斜率
new = [r for r in ok if r["var"] not in used]
print("\n" + "=" * 70)
print(f"合格(span ≥ {SPAN_MIN}):{len(ok)} 道;其中**我从没算过斜率的新话题**:{len(new)} 道")
for r in new: print(f"   ✅ {r['src']:5s} {r['var']:10s} span={r['span']:.4f}  {r['label']}")
print(f"\n不合格:{[r['var'] for r in rows if r.get('span') is not None and r['span'] < SPAN_MIN]}")
print(f"等级不足(<4 级,无法取三个嵌套切点):"
      f"{[r['var'] for r in rows if r.get('span') is None]}")

# ⛔ 本轮没算过任何 lnOR —— 用一个可失败的断言把它钉住
# ⚠ 第一版这行**自指**了:它在自己的源码里匹配到了自己写的模式,于是必然失败。
#    一个检查若把自己算进被检查的总体,它就不是检查。只看这行之前的源码。
_src = open(__file__).read().split("# ⛔ 本轮没算过任何 lnOR")[0]
assert not any(k in _src for k in ["np.log" + "(", "math.log" + "("]), "本轮不得计算 lnOR"
decision = ("猜想**可被检验**:下一轮做多点样本外" if len(new) >= 3 else
            f"**只有 {len(new)} 个新合格话题(<3)-> 这个猜想在现有仪器上无法被多点检验**,"
            f"写进页面「做不到什么」,留在冻结文件里")
print(f"\n预注册的决策 ⇒ {decision}")
print("⚠ 这不是一个 kill:span 只说**量表挪不挪得动**,不说斜率是否存在。"
      "把筛选叫成检验就是「验证串不是计算」那一类。")
json.dump(dict(span_min=SPAN_MIN, rows=rows, eligible=[r["var"] for r in ok],
               eligible_new=[r["var"] for r in new], decision=decision,
               computed_any_lnor=False, unchallenged=True),
          open(OUT / "pure_span_screen.json", "w"), indent=1)
print(f"\nwrote {OUT/'pure_span_screen.json'}")
