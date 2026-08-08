"""E02·A12·R654 —— 性那一侧,有没有同型的「对谁」维度?

`#617` 的 NEXT。⚠ **不许把「暴力有梯度」直接推广到性** —— 本轮就是去查它到底有没有。

⚠ **`#616` 刚教过我一次:手读四个标题判了零,机器扫描找回五个。**
   所以本轮**先手读七个,再用机器全扫复核** —— 而 **P5★:零必须先过正对照。**

INSTRUMENT(硬规则②):SCCS/D-PLACE,1,781 个变量,民族志编码。

G1 ESTIMAND(先于方法):**覆盖普查**,不是相关。
  `K_sex_target` = 满足两条的 SCCS 变量个数:
  (a) 关于一个**性行为**;(b) 把**同一个行为**按**对方是谁**分档
      (亲属/非亲属 · 本群/外群 · 已婚/未婚 · 同族/异族)。
  **判据先写死,(b) 明确排除**:「哪个性别的**行为者**被罚」「**谁承受**惩罚」「什么**情境**下」——
  这三类都不是「对谁做」。

WORLDS:**A** `K_sex_target >= 2` ⇒ 可以跑同型分离器 ·
  **B** `K_sex_target == 0` ⇒ **性那一侧在 SCCS 里没有距离维度,写进页面「做不到什么」**
CONTROLS:
  正对照(**对搜索本身**):同一套扫描指向**暴力**,必须找回 `781/782/783` 与 `1768/1769/1770`
  —— 它们已知是按「对谁」分档的。**这是「零」可采信的前提(P5★)。**
  g=0:指向 `Subsistence, Economy` 类目,必须返回 0。
  安慰剂:把「对谁」词表换成一组与对象无关的词(`season`/`crop`/`tool`)-> 命中必须落在别处。
KILL(条件式,预注册):
  if 正对照找回 >= 4 个 and g=0 返回 0:
      `K_sex_target == 0` -> **W-B**;`>= 2` -> **W-A**;`== 1` -> 报计数
  else: UNVERIFIED —— 搜索不合格,任何零都不 admissible
G3:全部性相关变量的**真实非缺失覆盖**表 + 扫描判定,含被判「不是对谁」的并注明它实际在量什么。
G4:词表宽/窄 × 是否要求码里 >=2 个目标类别 × n 门槛 {0, 30, 60}。
IMPOSSIBLE(不写 planned):**普查不是统计检验** ——
  只能说「这份数据里有没有这种变量」,**不能说「这些社会有没有这种规范」** ·
  「是不是性行为」「是不是按对谁分档」**由我判定**,不是数据发现的 · `[unchallenged]`
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
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
COV = {c: int(W[c].notna().sum()) for c in W.columns}
print(f"仪器 = SCCS/D-PLACE · {V.id.nunique()} 变量 · 186 社会")

# 判据先写死(在看结果之前)
SEXP = re.compile(r'sex|premarital|extramarital|adulter|incest|homosex|coitus|intercourse|'
                  r'virgin|chastit|concubin|polygyn|mating|rape', re.I)
VIOP = re.compile(r'violence|aggress|homicide|assault|feud|warfare|raid', re.I)
# 「对谁」= 码里区分 >=2 个**对象**类别(而不是行为者性别 / 谁受罚 / 情境)
TGT = re.compile(r'\b(kin|kinsm|relative|clan|lineage|cousin|sibling|in-?law|affine|'
                 r'own (group|community|society|ethnic)|other (group|societ|ethnic|people|tribe)|'
                 r'outsider|stranger|foreign|local community|same society|same ethnic|'
                 r'married|unmarried|spouse of another|another man\'s wife)\b', re.I)
NOTTGT = re.compile(r'who (is |)punish|punishment ot|by (the |)(husband|family|government)|'
                    r'(male|female|men|women|boys|girls) (are |is |)punish|double standard|'
                    r'single standard|unless pregnan|limited context', re.I)


# ⚠ 第一版扫描**只读码**,于是正对照只找回 2/6 —— `SCCS781/782/783` 的「对谁」
#   (local community / same society / other societies)**写在标题里,不在码里**。
#   **一个只读码的扫描,对我要找的那种设计恰好是瞎的。** 修法:标题与码一起读。
def profile(vid, title=""):
    cd = C[C.var_id == vid].dropna(subset=["name"])
    names = [str(x) for x in cd.name]
    if len(names) < 2: return None
    names = names + [str(title)]          # 标题作为第 n+1 条「码」参与判定
    tgt = sum(bool(TGT.search(x)) for x in names)
    bad = sum(bool(NOTTGT.search(x)) for x in names)
    return dict(names=names, tgt=tgt, bad=bad, is_tgt=(tgt >= 2 and tgt > bad))


# ⚠ 第一版 g=0 自身构造错了:传 `cat=` 时它**跳过了标题过滤**,于是变成「另一套扫描」——
#   realstat 的「对照因它自己的原因而失败」:**对照的两侧不是同一个对象**。
#   修法:g=0 必须是**同一套两级过滤**(标题匹配主题 AND 码里有 >=2 个目标类别),只换主题词表。
def sweep(pat):
    hits = []
    for r in V.itertuples():
        if not pat.search(str(r.title)): continue
        p = profile(r.id, r.title)
        if p is None: continue
        hits.append(dict(id=r.id, cov=COV.get(r.id, 0), title=str(r.title)[:66],
                         is_tgt=p["is_tgt"], tgt=p["tgt"], bad=p["bad"],
                         codes=" | ".join(p["names"])[:130]))
    return pd.DataFrame(hits)


# ── 先把 `#617` 点名的七个的真实覆盖钉住(G3 的第一部分)──────────
SEV7 = ["SCCS165", "SCCS169", "SCCS176", "SCCS960", "SCCS961", "SCCS962", "SCCS964"]
print("\n=== G3-a:`#616` 把这七个都报成 n=186 —— 真实非缺失覆盖 ===")
for x in SEV7:
    rows = int(D[D.var_id == x].soc_id.nunique())
    print(f"  {x:9s} 行数 {rows:3d} · **真实覆盖 {COV.get(x,0):3d}** · "
          f"按「对谁」分档? **{profile(x, V.set_index('id').loc[x,'title'])['is_tgt']}**  | {V.set_index('id').loc[x,'title'][:52]}")

# ── 正对照:同一套扫描指向暴力 ───────────────────────────────
vio = sweep(VIOP)
vio_hit = vio[vio.is_tgt]
MUST = {"SCCS781","SCCS782","SCCS783","SCCS1768","SCCS1769","SCCS1770"}
found = sorted(MUST & set(vio_hit.id))
print(f"\n=== 正对照:扫描指向**暴力** -> 命中 {len(vio)} 个,判为「按对谁分档」**{len(vio_hit)}** 个 ===")
for r in vio_hit.itertuples(): print(f"  ✅ {r.id:9s} cov={r.cov:3d} | {r.title}")
print(f"  **必须找回的 6 个里找回 {len(found)}**:{found}")

# ── 主问:扫描指向性 ─────────────────────────────────────────
sex = sweep(SEXP)
sex_hit = sex[sex.is_tgt]
print(f"\n=== G3-b:扫描指向**性** -> 命中 {len(sex)} 个,判为「按对谁分档」**{len(sex_hit)}** 个 ===")
for r in sex_hit.itertuples(): print(f"  ✅ {r.id:9s} cov={r.cov:3d} | {r.title}\n        码: {r.codes}")
print("\n  --- 被判「不是按对谁分档」的,取覆盖最高的 8 个,并注明它实际在量什么 ---")
for r in sex[~sex.is_tgt].sort_values("cov", ascending=False).head(8).itertuples():
    print(f"  ⛔ {r.id:9s} cov={r.cov:3d} | {r.title}\n        它实际在量:{r.codes[:110]}")
K = int(len(sex_hit))

SUBP = re.compile(r'crop|cultivat|fishing|herding|harvest|granary|plough|irrigat|milking|gathering', re.I)
g0 = sweep(SUBP)
g0n = int(g0.is_tgt.sum())
print(f"\n  g=0:**同一套两级过滤**,主题词换成耕作/渔猎/畜牧({len(g0)} 个标题命中)"
      f"-> 判为「按对谁分档」**{g0n}** 个(须 0)")

G = Gate("性那一侧,有没有同型的「对谁」维度?")
# ⚠ 第一版这里 floor=0.0 spread=0.5 -> 门槛 1.0,**比我在 docstring 里写下的「>=4」松**。
#   判据与代码不一致就是没有判据。对齐:floor=3.0, spread=0.5 -> 门槛 4.0。
pos_ok = G.positive_control("正对照:扫描指向暴力必须找回 >=4 个已知按对谁分档的变量",
                            planted=float(len(found)), floor=3.0, spread=0.5)
pla_ok = G.negative_control("g=0:同一套两级过滤,主题换成耕作/渔猎/畜牧", null=float(g0n), effect=float(max(len(found), 1)),
                            null_spread=0.5, null_kind="与「对谁」无关的变量类目")
if pos_ok and pla_ok:
    verdict = ("W-B:**性那一侧在 SCCS 里没有「对谁」维度**(写进页面「做不到什么」)" if K == 0 else
               ("W-A:**有 %d 个,可以跑同型分离器**" % K if K >= 2 else f"报计数:K = {K},不下判决"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 搜索不合格(正对照 {pos_ok} · g=0 {pla_ok});**任何零都不 admissible**"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:目标类别数门槛 × 覆盖门槛 ===")
spec = []
for need in (1, 2, 3):
    for nmin in (0, 30, 60):
        k = sum(1 for r in sex.itertuples() if r.tgt >= need and r.tgt > r.bad and r.cov >= nmin)
        kv = sum(1 for r in vio.itertuples() if r.tgt >= need and r.tgt > r.bad and r.cov >= nmin)
        spec.append(dict(need=need, nmin=nmin, K_sex=k, K_vio=kv))
        print(f"  需 >={need} 个目标类别 · 覆盖 >={nmin:2d}: **性 {k}** · 暴力 {kv}")

# ══════════════════════════════════════════════════════════════════
# ⛔ 第三次修,而这一次修的是**单位**,不是词表 —— 前两次都没修对地方。
#   正对照两次都只找回 2/6,因为 **Ross 的每个变量只带一个目标**(写在标题里),
#   **三个变量才构成那个距离维度**;而我的规则要求「一个变量里有 >=2 个目标类别」。
#   realstat:**先把仪器的单位与主张的单位写成两个字符串,并要求它们相等,再去设计对照。**
#   仪器的单位 = 「变量」· 主张的单位 = 「三元组」 ⇒ 不等 ⇒ 前两版从设计上就看不见它。
#   修法:**按「同一词干 + 只有目标短语不同」的变量族来搜。**
print("\n" + "="*66)
print("=== 修正的单位:按**变量族**搜(同一词干,只有目标短语不同)===")
TGTPH = re.compile(r'(members of the local community|members of same ethnic group[^,]*|'
                   r'members of other ethnic groups?|people in other societies|'
                   r'members of the same society[^,]*|own (group|community|kin)|'
                   r'other (group|societies|tribes?|peoples?)|kin|non-?kin|'
                   r'unmarried|married|stranger|outsider|foreigner)', re.I)

def stem(title):
    s = TGTPH.sub("<TARGET>", str(title))
    s = re.sub(r'[^a-z<>]+', ' ', s.lower()).strip()
    return s

def families(pat):
    g = {}
    for r in V.itertuples():
        if not pat.search(str(r.title)): continue
        if not TGTPH.search(str(r.title)): continue
        g.setdefault(stem(r.title), []).append((r.id, COV.get(r.id, 0), str(r.title)[:64]))
    return {k: v for k, v in g.items() if len(v) >= 2}

fam_vio = families(VIOP); fam_sex = families(SEXP)
print(f"\n  正对照(暴力):**{len(fam_vio)} 个族**")
for k, v in fam_vio.items():
    print(f"    族「{k[:44]}」{len(v)} 个:")
    for i, c, ti in v: print(f"       {i:9s} cov={c:3d} | {ti}")
print(f"\n  主问(性):**{len(fam_sex)} 个族**")
if not fam_sex: print("    (一个都没有)")
for k, v in fam_sex.items():
    print(f"    族「{k[:44]}」{len(v)} 个:")
    for i, c, ti in v: print(f"       {i:9s} cov={c:3d} | {ti}")
found2 = sorted(MUST & {i for v in fam_vio.values() for i, _, _ in v})
K2 = sum(len(v) for v in fam_sex.values())
print(f"\n  **正对照按族搜找回 {len(found2)}/6:{found2}**")
print(f"  **性那一侧的族成员总数 K2 = {K2}**")
G2 = Gate("修正单位之后:性那一侧有没有「对谁」的变量族?")
pos2 = G2.positive_control("正对照(按族):必须找回 >=4 个已知按对谁分档的变量",
                           planted=float(len(found2)), floor=3.0, spread=0.5)
g0_fam = families(SUBP)
g0_2 = sum(len(v) for v in g0_fam.values())
print(f"  g=0(耕作/渔猎/畜牧,同一套按族搜)-> {g0_2} 个族成员(须 0)")
pla2 = G2.negative_control("g=0(按族):主题换成耕作/渔猎/畜牧", null=float(g0_2),
                           effect=float(max(len(found2), 1)), null_spread=0.5,
                           null_kind="与「对谁」无关的主题")
if pos2 and pla2:
    v2 = ("W-B:**性那一侧在 SCCS 里没有「对谁」的变量族**(写进页面「做不到什么」)" if K2 == 0
          else f"W-A:**性那一侧有 {K2} 个族成员**")
    print(f"\n  控制齐备 ⇒ 评判。**{v2}**")
else:
    v2 = f"UNVERIFIED(按族)—— 正对照 {pos2} · g=0 {pla2}"
    print(f"\n  ⚠ {v2}")
print(G2)
verdict = verdict + " ‖ 修正单位后:" + v2

json.dump(dict(instrument="SCCS / D-PLACE", coverage_of_7={x: COV.get(x, 0) for x in SEV7},
               sex_hits=sex.to_dict("records"), vio_hits=vio.to_dict("records"),
               K_sex_target=K, positive_found=found, g0=g0n, spec_curve=spec,
               verdict=verdict, family_positive=found2, K_family_sex=K2, g0_family=g0_2,
               families_vio={k: [x[0] for x in v] for k, v in fam_vio.items()},
               families_sex={k: [x[0] for x in v] for k, v in fam_sex.items()}, unchallenged=True),
          open(OUT/"does_sex_have_a_distance_axis.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'does_sex_have_a_distance_axis.json'}")
