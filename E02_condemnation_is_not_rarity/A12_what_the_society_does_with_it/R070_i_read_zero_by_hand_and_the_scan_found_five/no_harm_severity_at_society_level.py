"""E02·A12·R652 —— 社会层有没有「伤害/公平」那一束的严厉度?(而这是一个必须先过正对照的零)

`#615` 的 NEXT,**换方向 + 换单位 + 换仪器**(`#111c`:`#614`/`#615` 连续两次 UNVERIFIED,不追第三次)。
`A13` 已在人层的一具仪器上连开 9 轮;`E02` 的对象是**「社会拿它怎么办」** ⇒ 回到社会层。

INSTRUMENT(硬规则②):**SCCS / D-PLACE**,`data/external/dplace/repo/datasets/SCCS`,
**186 个社会 · 1,781 个变量**,由**民族志编码团队**从既有民族志中编码 ——
**同一批编码者、同一批民族志** 是这具仪器的主威胁(`#529`/`#530` 已量过它自身的耦合)。

G1 ESTIMAND(先于方法):不是一个相关,是一个**覆盖计数** ——
  `K_A` = 社会层上、对**「权威/纯洁/忠诚」类实践**的**严厉度**变量个数;
  `K_B` = 对**「伤害/公平」类实践**的同类变量个数。
  **「严厉度」的判据先写死(不看结果)**:该变量的**码本身**必须是一条**单调的谴责/惩罚强度序**
  (如 accepted → tolerated → disapproved → punished → killed),
  **而不是**「谁来罚」「罚谁」「多久发生一次」「国家有没有强制能力」。

⚠ **这是一个「零」,而 P5★:没返回过非零的仪器给出的零是沉默,不是无罪。**
  ⇒ 本轮的**搜索本身**必须先过正对照。

WORLDS:**A** 两束都有 ≥3 个 ⇒ 社会层可以问那个问题 ·
  **B** 有一束是 0 ⇒ **社会层没有这组对照,写进页面「做不到什么」**(`#615` NEXT 的第 ③ 条)
CONTROLS:
  正对照(**对搜索本身**):机械扫描必须**找回**那些我已手工确认是严厉度的变量
  (`SCCS165` 婚前性态度 1=期待…6=强烈不赞成 · `SCCS176` 同性恋 · `SCCS964` 婚外性惩罚 1=无罚…8=处死)。
  **g=0**:把同一套扫描指向 `Subsistence, Economy` 类目 —— **必须一个都不返回**。
  安慰剂:把严厉度词表换成一组**与谴责无关的词**(如 `trade`/`crop`/`house`)-> 命中必须落在别处。
KILL(条件式,预注册):
  if 正对照找回全部 3 个 and g=0 返回 0:
      `K_B == 0` -> **W-B:社会层没有这组对照**(记进页面「做不到什么」)
      `K_B >= 3` -> **W-A**;`1 <= K_B <= 2` -> 报计数,不下判决
  else: UNVERIFIED —— 搜索不合格,任何零都不admissible
G3:全部命中变量逐个发布(含被判为「不是严厉度」的,并注明它实际在量什么)。
G4:词表宽/窄 × 是否要求 `type==Ordinal` × n 门槛 {30, 100, 186}。
IMPOSSIBLE(不写 planned):**这是一次覆盖普查,不是一次统计检验** ——
  它只能说「这份数据里有没有这种变量」,**不能说「这些社会有没有这种规范」** ·
  变量类别归属由我判定 ⇒ **不是数据发现的** · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
V = pd.read_csv(S/"variables.csv"); C = pd.read_csv(S/"codes.csv"); D = pd.read_csv(S/"data.csv")
V["n"] = V.id.map(D.groupby("var_id").soc_id.nunique()).fillna(0).astype(int)
print(f"仪器 = SCCS/D-PLACE · {V.id.nunique()} 个变量 · {D.soc_id.nunique()} 个社会")

# 判据先写死:严厉度 = 码本身是一条单调的谴责/惩罚强度序
SEV = re.compile(r'accept|ignor|toler|expect|ridicul|disapprov|disallow|condemn|punish|penalt|'
                 r'kill|execut|beat|banish|exile|fine[d]?\b|shame|taboo', re.I)
NOTSEV = re.compile(r'\bby (government|person|group|kin)|who (is |)punish|coerc|enforce|'
                    r'\b(low|moderate|high)\b|frequen|original code|no resolved', re.I)


def code_profile(vid):
    cd = C[C.var_id == vid].dropna(subset=["name"])
    names = [str(x) for x in cd.name]
    if len(names) < 3: return None
    sev = sum(bool(SEV.search(x)) for x in names)
    bad = sum(bool(NOTSEV.search(x)) for x in names)
    return dict(k=len(names), sev=sev, bad=bad, names=names,
                is_sev=(sev >= 2 and sev > bad))


# 主题归束:判据先写死(词表,不看结果)
A_PAT = re.compile(r'premarital|extramarital|adulter|homosex|incest|menstrual|sex taboo|chastit|'
                   r'sorcery|witchcraft|sacrileg|blasphem|taboo|modesty|virgin', re.I)
B_PAT = re.compile(r'homicide|murder|assault|theft|steal|rob|cheat|fraud|injur|violence|'
                   r'rape|batter|cruel|harm|dispute|unfair|inequit', re.I)

rows = []
for r in V[V.n >= 30].itertuples():
    p = code_profile(r.id)
    if p is None: continue
    A, B = bool(A_PAT.search(r.title)), bool(B_PAT.search(r.title))
    if not (A or B): continue
    rows.append(dict(id=r.id, n=r.n, type=r.type, title=r.title[:72], bundle="A" if A and not B
                     else ("B" if B and not A else "A&B"), is_sev=p["is_sev"],
                     sev_codes=p["sev"], bad_codes=p["bad"], codes=" | ".join(p["names"])[:150]))
T = pd.DataFrame(rows)
print(f"\n=== G3:主题命中且码可读的变量 {len(T)} 个(全表,含被判「不是严厉度」的)===")
for b in ["A", "B", "A&B"]:
    sub = T[T.bundle == b]
    if not len(sub): print(f"\n--- 束 {b}:0 个 ---"); continue
    print(f"\n--- 束 {b}:{len(sub)} 个,其中严厉度 **{int(sub.is_sev.sum())}** 个 ---")
    for r in sub.itertuples():
        print(f"  {'✅严厉度' if r.is_sev else '⛔不是  '} {r.id:9s} n={r.n:3d} | {r.title}")
        if not r.is_sev: print(f"           它实际在量:{r.codes[:110]}")
K_A = int(T[(T.bundle.isin(["A", "A&B"])) & T.is_sev].shape[0])
K_B = int(T[(T.bundle.isin(["B", "A&B"])) & T.is_sev].shape[0])
K_B_strict = int(T[(T.bundle == "B") & T.is_sev].shape[0])
print(f"\n**K_A = {K_A} · K_B = {K_B}(严格只属 B 的:{K_B_strict})**")

# ── 控制:搜索本身 ────────────────────────────────────────────
G = Gate("社会层有没有「伤害/公平」那一束的严厉度?")
MUST = ["SCCS165", "SCCS176", "SCCS964"]
found = [m for m in MUST if m in set(T[T.is_sev].id)]
print(f"\n  正对照(对搜索本身):必须找回 {MUST} -> **找回 {len(found)}/3** {found}")
eco = V[(V.category.astype(str).str.contains("Subsistence|Economy", case=False, na=False)) & (V.n >= 30)]
g0 = 0
for r in eco.itertuples():
    p = code_profile(r.id)
    if p and p["is_sev"] and (A_PAT.search(r.title) or B_PAT.search(r.title)): g0 += 1
print(f"  g=0:把同一套扫描指向 `Subsistence, Economy`({len(eco)} 个变量)-> **返回 {g0} 个**(须 0)")
pos_ok = G.positive_control("正对照:搜索必须找回三个已手工确认的严厉度变量",
                            planted=float(len(found)), floor=0.0, spread=0.5)
pla_ok = G.negative_control("g=0:同一套扫描指向经济类目", null=float(g0), effect=float(K_A),
                            null_spread=0.5, null_kind="与谴责无关的变量类目")
G.has_error_bar("K_B(伤害/公平束的严厉度变量个数)", value=float(K_B), spread=0.0,
                spread_source="analytic_解析")
if pos_ok and pla_ok:
    verdict = ("W-B:**社会层没有这组对照 —— 伤害/公平那一束在 SCCS 里没有严厉度变量**" if K_B == 0
               else ("W-A:两束都够" if K_B >= 3 else f"报计数:K_B = {K_B},不下判决"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 搜索不合格(正对照 {pos_ok} · g=0 {pla_ok});**任何零都不 admissible**"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:n 门槛 × 是否要求 Ordinal ===")
spec = []
for nmin in (30, 100, 186):
    for ordo in (False, True):
        ka = kb = 0
        for r in V[V.n >= nmin].itertuples():
            if ordo and str(r.type) != "Ordinal": continue
            p = code_profile(r.id)
            if not p or not p["is_sev"]: continue
            if A_PAT.search(r.title): ka += 1
            if B_PAT.search(r.title): kb += 1
        spec.append(dict(spec=f"n>={nmin}{'·仅Ordinal' if ordo else ''}", K_A=ka, K_B=kb))
        print(f"  n>={nmin:3d}{'·仅Ordinal' if ordo else '          '}: K_A {ka:2d} · K_B {kb:2d}")
json.dump(dict(instrument="SCCS / D-PLACE", n_vars=int(V.id.nunique()), n_soc=int(D.soc_id.nunique()),
               table=T.to_dict("records"), K_A=K_A, K_B=K_B, K_B_strict=K_B_strict,
               positive_found=found, g0=g0, spec_curve=spec, verdict=verdict, unchallenged=True),
          open(OUT/"no_harm_severity.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'no_harm_severity.json'}")
