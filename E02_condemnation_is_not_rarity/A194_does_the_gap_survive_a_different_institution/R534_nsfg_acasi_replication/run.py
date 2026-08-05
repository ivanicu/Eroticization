"""E02·A194·R534 — 换一个机构、换一种施测模式,那道差距还在不在

`#488` 的 NEXT:换回跨仪器(规则 ④),先按规则 ① 数 NSFG 有没有态度变量。

⚠ **规则 ① 差点被我自己用坏,记下来。** 我先用关键词搜 `wrong/approv/moral/should/attitude`,
**全部返回 0**,而正确答案是**有** —— 题目写作
`IH-1 Sexual relations between two same-sex adults is all right`。
**若我停在那个 0 上,就会宣布「NSFG 结构上没有态度变量」并关掉这条线,而那是假的。**
⇒ **一个关键词搜出的 0,和一个变量表读出的 0,不是同一种 0。**

⚠ 而 `.dat` 是**无表头定宽 ASCII** —— 我在 `#487` 之前下了 545 MB 却**一个布局文件都没下**,
本轮补下 `.dct`/`.sas`。**「下载成功」不等于「可读」。**

NSFG 与 GSS 同构的一对(不同机构 · 不同题面 · 不同抽样 · **ACASI 自填**):
  `samesex`    IH-1 谴责题,AGDGFMT:1 强同意 · 2 同意 · 3 不同意 · 4 强烈不同意
  ⛔ **`5 = 既不同意也不反对`,被放在数值量程的顶端 —— 非序数等级,必须剔除。**
     这是本项目第三次撞到同型陷阱(`SCCS176 code 2`、`SCCS743 「婚姻头几年」`)。8/9 = 拒答/不知道。
  `samesexany` 曾有过同性性接触(行为)· `samyearnum` 过去 12 个月同性伴侣数 · `attract` 吸引取向

G1 ESTIMAND(先于方法):`lnOR(谴责, 行为阳性)`,谴责 := `samesex ∈ {3,4}`。
  ⚠ 这**不是**复现 `#486b` 的**那个数**(话题不同:同性 vs 色情),
     而是复现**那个设计** —— 差距在不同机构、不同模式下**存不存在、多大**。
  ⚠ NSFG 只有 2 波可用 ⇒ **无时代杠杆,结构上无法复现「五十一年没动」**。先写下。

⭐ 而 ACASI 正是 `#486b` 留下的那个未解混淆要的东西:
  GSS 的这些题是**访员施测**,NSFG 的敏感题是**自填**。
WORLDS:
  W-INTERVIEWER 差距有一部分是访员在场造的 -> 自填模式下差距**明显更小**
  W-REAL        差距反映真实关联           -> 量级相近
  | World         | now | 明显更小 | 量级相近 |
  | W-INTERVIEWER | 0.4 | 0.85     | 0.10     |
  | W-REAL        | 0.6 | 0.15     | 0.90     |
⛔ **STRONGEST CONFOUND,写在跑之前:话题也变了**(同性 vs 色情)。
  量级差**不可单独归因于模式**。这是跨调查比较的结构性限制,在此不可拆。

CONTROLS:
  正对照 `samesex` × 宗教出席 `attndnow` —— 门槛 = **NSFG 内部实测参照分布 q95**(`#485a`)
  阴性   参照分布本身(`samesex` × 一批与性无关的 NSFG 变量),**测量,非挑选**
  精度   人层 bootstrap
KILL(条件式,预注册):
  if 正对照触发 and 参照分布中位 < 0.5·|效应|:
      |lnOR| 明显小于 GSS 的 −1.42(比如 |lnOR| < 0.7)-> W-INTERVIEWER 得到支持
      |lnOR| 与之同量级(0.7 ~ 2.2)                  -> W-REAL
      两者都不是                                      -> 未决
  else: UNVERIFIED
IMPOSSIBLE:话题与模式不可拆 · 只有 2 波 ⇒ 无时代分辨 · 无干预 ⇒ 非因果 ·
  未派对抗 agent(会话约束)⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, math, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
GSS_LNOR = -1.4183          # `#487c` 的 xmovie 臂,作为量级参照


def parse_dct(path):
    """Stata 字典 -> {name: (start0, width, label)}。定宽 ASCII 的唯一入口。"""
    out = {}
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(path, errors="replace"):
        m = pat.search(line)
        if m:
            col, name, w, lab = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
            out[name.lower()] = (col - 1, w, lab)
    return out


def read_fixed(dat, layout, names):
    cols = {n: layout[n] for n in names if n in layout}
    missing = [n for n in names if n not in layout]
    if missing: print(f"  ⚠ 字典里没有: {missing}")
    rows = {n: [] for n in cols}
    with open(dat, errors="replace") as f:
        for line in f:
            for n, (s, w, _) in cols.items():
                v = line[s:s + w].strip()
                rows[n].append(float(v) if v not in ("", ".") else np.nan)
    return pd.DataFrame(rows)


WAVES = [("2017_2019_Fem", "2017_2019_FemRespSetup.dct", "2017_2019_FemRespData.dat"),
         ("2011_2013_Fem", "2011_2013_FemRespSetup.dct", "2011_2013_FemRespData.dat")]
NEED = ["samesex", "samesexany", "samyearnum", "attract", "attndnow", "age_r",
        "educat", "hieduc", "religion", "poverty", "hisp", "rmarital", "chsuppor"]

frames = {}
for nm, dct, dat in WAVES:
    lp = NS / "setup" / dct
    if not lp.exists(): print(f"跳过 {nm}: 无字典"); continue
    lay = parse_dct(lp)
    print(f"\n=== {nm}: 字典 {len(lay)} 个变量 ===")
    df = read_fixed(NS / dat, lay, NEED)
    frames[nm] = df
    print(f"  行数 {len(df)}")
    for c in ["samesex", "samesexany", "samyearnum"]:
        if c in df:
            vc = df[c].value_counts().head(7).to_dict()
            print(f"  {c:12s} n={df[c].notna().sum():6d}  值分布(前7)={vc}")

W = "2017_2019_Fem"
d = frames[W].copy()
# ⛔ 剔除非序数与缺失码
d = d[d.samesex.isin([1, 2, 3, 4])].copy()
d["condemn"] = d.samesex.isin([3, 4]).astype(float)
# ⛔ NSFG 惯例:1=Yes, 5=No, 7=Not ascertained。不是 {0,1}。
# 我第一版写了 isin([0,1]),只留下 1,基率变成恰好 1.0000 —— `#296b` 那一类,
# 一个落在边界上的比率是仪器失败,不是极端成功。下面的 assert 把它固定住。
d["ever"] = np.where(d.samesexany == 1, 1.0, np.where(d.samesexany == 5, 0.0, np.nan))
assert 0.01 < np.nanmean(d.ever) < 0.99, f"ever 基率 {np.nanmean(d.ever)} 落在边界 -> 码读错了"
d["yr"] = np.where(d.samyearnum.between(0, 90), (d.samyearnum > 0).astype(float), np.nan)
print(f"\n=== 主格 {W}(已剔除 samesex ∈ {{5,8,9}})===")
print(f"n={len(d)}  谴责率={d.condemn.mean():.4f}  曾有过={d.ever.mean():.4f}  "
      f"过去一年有={d.yr.mean():.4f}")


def lnor(dd, ycol):
    a, b = dd[dd.condemn == 1], dd[dd.condemn == 0]
    a, b = a[ycol].dropna(), b[ycol].dropna()
    if len(a) < 30 or len(b) < 30: return np.nan, 0
    p1, p0 = a.mean(), b.mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan, 0
    return math.log((p1 / (1 - p1)) / (p0 / (1 - p0))), len(a) + len(b)


main, n_main = lnor(d, "ever")
print(f"\n主:lnOR(谴责, 曾有过同性性接触) = {main:+.4f}  n={n_main}")
print(f"   GSS 参照(色情,访员施测)= {GSS_LNOR:+.4f}")


def boot(dd, ycol, B=600, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(B):
        s = dd.iloc[rng.integers(0, len(dd), len(dd))]
        v, _ = lnor(s, ycol)
        if np.isfinite(v): out.append(v)
    return np.array(out)


bs = np.concatenate([boot(d, "ever", 600, s) for s in SEEDS])
lo, hi = np.quantile(bs, [.025, .975])
print(f"   95% CI [{lo:+.4f}, {hi:+.4f}]  sd={bs.std():.4f}  "
      f"seed_spread={np.std([boot(d,'ever',600,s).mean() for s in SEEDS]):.5f}")

# ---------------------------------------------------------------- G4
print("\n=== G4 规格曲线(全格公布)===")
spec = []
for wv in frames:
    dd0 = frames[wv]
    dd0 = dd0[dd0.samesex.isin([1, 2, 3, 4])].copy()
    dd0["condemn"] = dd0.samesex.isin([3, 4]).astype(float)
    dd0["ever"] = np.where(dd0.samesexany == 1, 1.0, np.where(dd0.samesexany == 5, 0.0, np.nan))
    dd0["yr"] = np.where(dd0.samyearnum.between(0, 90), (dd0.samyearnum > 0).astype(float), np.nan)
    dd0["strict"] = (dd0.samesex == 4).astype(float)
    for ycol in ("ever", "yr"):
        v, n = lnor(dd0, ycol)
        if np.isfinite(v): spec.append(dict(wave=wv, cond="disagree", beh=ycol, lnor=v, n=n))
    d2 = dd0.copy(); d2["condemn"] = d2["strict"]
    for ycol in ("ever", "yr"):
        v, n = lnor(d2, ycol)
        if np.isfinite(v): spec.append(dict(wave=wv, cond="strongly", beh=ycol, lnor=v, n=n))
for s in spec:
    print(f"  {s['wave']:14s} {s['cond']:9s} {s['beh']:5s} lnOR={s['lnor']:+.4f}  n={s['n']}")
ls_ = [s["lnor"] for s in spec]
sg = [np.sign(v) for v in ls_]; dom = max(set(sg), key=sg.count)
print(f"\nspec_survival: {sg.count(dom)}/{len(sg)} = {sg.count(dom)/len(sg):.0%} 同号 ({dom:+.0f});"
      f" 范围 [{min(ls_):+.4f}, {max(ls_):+.4f}]")

# ---------------------------------------------------------------- 控制
G = Gate("换一个机构、换一种施测模式,那道差距还在不在?(NSFG,ACASI)")
ref = []
for c in ["age_r", "educat", "hieduc", "poverty", "hisp", "rmarital", "chsuppor", "religion"]:
    if c not in d.columns: continue
    s = d.dropna(subset=[c])
    s = s[s[c] < 90]
    if len(s) < 500 or s[c].nunique() < 3: continue
    r = float(np.corrcoef(s.condemn, s[c])[0, 1])
    ref.append(dict(var=c, r=r, n=len(s)))
ar = np.array([abs(x["r"]) for x in ref])
print(f"\n参照分布(NSFG 内 {len(ref)} 个变量)|r| 中位={np.median(ar):.4f} q95={np.quantile(ar,.95):.4f}")
for x in sorted(ref, key=lambda z: -abs(z["r"]))[:4]:
    print(f"   |{x['r']:+.4f}| n={x['n']:5d}  {x['var']}")

pc = d.dropna(subset=["attndnow"]); pc = pc[pc.attndnow < 90]
pc_r = float(np.corrcoef(pc.condemn, pc.attndnow)[0, 1])
T = float(np.quantile(ar, .95))
print(f"正对照 谴责×宗教出席 r={pc_r:+.4f} n={len(pc)}  门槛(参照 q95)={T:.4f}")
pc_ok = G.positive_control("正对照:谴责×宗教出席(门槛=NSFG 内部参照 q95)",
                           planted=abs(pc_r), floor=T, spread=1e-9)
nc_ok = G.negative_control("阴性:NSFG 内部参照分布中位(测量,非挑选)",
                           null=float(np.median(ar)), effect=abs(main),
                           null_spread=float(ar.std()), null_kind="同问卷无关变量参照分布")
G.has_error_bar("lnOR", value=main, spread=float(bs.std()), spread_source="bootstrap_人层")

print("\n" + "=" * 70)
if pc_ok and nc_ok:
    if abs(main) < 0.7:
        verdict = f"|lnOR|={abs(main):.4f} < 0.7 -> W-INTERVIEWER 得到支持"
    elif 0.7 <= abs(main) <= 2.2:
        verdict = (f"|lnOR|={abs(main):.4f} 与 GSS 的 {abs(GSS_LNOR):.2f} 同量级 -> "
                   f"**W-REAL:差距在自填模式下依然存在且同量级**")
    else:
        verdict = f"|lnOR|={abs(main):.4f} 超出预注册区间 -> 未决"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:**话题也变了**(同性 vs 色情),"
          "量级相近可能是两个不同话题恰好给出相近的数;要拆开需要**同一话题两种模式**,"
          "而 GSS 与 NSFG 都不提供。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(wave=W, n=int(len(d)), condemn_rate=float(d.condemn.mean()),
               ever_rate=float(d.ever.mean()), lnor=main, ci=[float(lo), float(hi)],
               gss_reference=GSS_LNOR, spec=spec, reference=ref,
               positive=dict(r=pc_r, threshold=T, ok=bool(pc_ok)),
               verdict=verdict, seeds=SEEDS, unchallenged=True,
               excluded_codes="samesex 5=neither / 8=refused / 9=DK"),
          open(OUT / "nsfg_acasi_replication.json", "w"), indent=1)
print(f"\nwrote {OUT/'nsfg_acasi_replication.json'}")
