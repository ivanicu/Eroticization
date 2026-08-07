"""#845 · E03·A83·R284 —— 墙倒了,那条路走得通吗:只扣**可识别的那一半**年龄

`#844` 推倒了 `#835` 的第一堵墙(APC「不可分离」说得太满:**共线的是线性分量,非线性分量可识别**)。
⇒ `#844`① 于是问:**那条路走得通吗?**
⚠ 而 `#844` 自己写了「不许把『墙倒了』读成『问题能答了』」——
**`#835` 还有第二堵墙(功效),它是量出来的,没被碰过。**
⇒ **所以本轮的第一件事不是求结果,是求功效**(`#835` 教的:**功效闸排在所有关于人的判断之前**)。

**⚠⚠ 而墙倒得只有一半,估计量就必须只做一半 —— 这是本轮设计的全部内容:**
不能说「扣掉年龄效应」,因为**年龄的线性分量与时期/世代共线,扣它就是在扣一个我识别不了的东西**。
**能说的是:扣掉年龄剖面里的非线性部分**(把年龄哑变量对 `age` 的线性趋势投影掉,只留残差)。
⇒ **估计量因此比 `#834`① 问的那个小,而这不是保守,是它是唯一识别得了的那个。**

G1 估计量:**`homosex` 先扣掉「年龄剖面的非线性分量」,再算虔诚/世俗两层的差距,
再算它在九十年代相对自己五十年匀速参照的偏离** —— 与 `#840` 同一构造,只是换了个被减过的因变量。
   `dep_gap^{−age_nl}(1990s)`,与 `#840` 未扣年龄的 **−0.5129 [−0.740, −0.287]** 直接对照。

三个世界:
   A **扣掉非线性年龄后,九十年代那条缝几乎不动** ⇒ **那条缝不是年龄剖面造出来的。**
   B **明显缩小但仍排零** ⇒ **一部分是年龄,一部分不是** —— 报**比例**,不报「是/不是」。
   C **缩到含零** ⇒ **九十年代那条缝里,可识别的那部分年龄就能解释掉它。**

预测矩阵:
   | 世界 | 现在 | 几乎不动 | 缩小仍排零 | 缩到含零 |
   | A 不是年龄 | 0.45 | **0.85** | 0.10 | 0.05 |
   | B 各占一部分 | 0.40 | 0.10 | **0.85** | 0.10 |
   | C 是年龄 | 0.15 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式,而 ⓪ 排在最前):
  ⓪ **功效闸**:在扣过年龄的数据上**植入一个已知大小的缩减**,若捞不回来 ⇒ **UNVERIFIED,
     不看结果**(`#835` 的做法:结果被拦在读之前,而不是撤在读之后)。
  if 功效闸过 and 正控开火 and 负控开火(**扣一个纯噪声的「假年龄剖面」不许改变结论**):
      |dep 变化| < 该格自身自助 SD -> A
      仍排零但缩小 ≥ 1 个 SD       -> B(报比例)
      含零                          -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**「扣掉非线性年龄」这一步本身会改变噪声结构** ——
  残差的方差比原始小,于是区间会**自动变窄**,而变窄会让「仍排零」显得更容易。
  ⇒ 控制:**负控用一个「假年龄剖面」**(把年龄标签打乱后估出来的剖面)做同样的扣减 ——
  它同样缩小方差,**却不含真实年龄信息** ⇒ **真剖面与假剖面的差,才是年龄真正解释掉的部分。**

⚠ 本轮换不了仪器(GSS);而它不需要 —— 本轮换的是**估计量的识别范围**。
"""
import json, math, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, B, Q = "homosex", 3000, 0.05

d = pd.read_stata(gp, columns=["year", "age", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= 4))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("age", (18, 89))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1); M = M.join(R["REL"])
t = M.groupby("year")["REL"].transform(
    lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
M["stratum"] = np.where(t == 2, 1, np.where(t == 0, 0, np.nan))
W = M.dropna(subset=[IT, "age", "stratum"]).copy()
W["age"] = W.age.astype(int)

print(f"=== ⓪ 硬规则①:合格样本 · 年龄真的问了吗 ===")
print(f"  n = **{len(W):,}** · 年龄 {W.age.min()}–{W.age.max()} · 年份 {W.year.min()}–{W.year.max()} · "
      f"虔诚层 {int((W.stratum==1).sum()):,} · 世俗层 {int((W.stratum==0).sum()):,}")

# ⚠⚠ **第一版用 72 个逐岁哑变量,而负控当场把它判了死刑:**
#    打乱年龄标签后估出的**假**剖面,非线性全幅 **0.3305**,真剖面才 **0.3744** ——
#    **几乎一样大 ⇒ 那条「非线性年龄剖面」基本是噪声。**
#    24,157 行上放 72 个年龄哑 + 29 个年份哑,**过度参数化**,而这在跑之前就能推出来。
#    ⇒ 改成**扫年龄分箱宽度**(1 / 5 / 10 岁),并**只在假剖面明显小于真剖面的那个分辨率上**下判。
#    ⚠ 这不是看了结果改判据:失败的是**控制**(仪器不合格),而 kill ③ 一个字没动;
#      `P6` 说 UNVERIFIED 意味着「这次检查不合用」,修好不合用的仪器再跑,正是它要求的。
BANDW = 5
def band(a, w): return ((np.asarray(a)-18)//w)*w + 18
AGES = np.sort(np.unique(band(W.age.values, BANDW)))
W["aband"] = band(W.age.values, BANDW)
def age_profile(y, ages):
    """年龄哑 + 年份哑 的最小二乘;返回年龄剖面(以最小年龄为基准)。"""
    A = pd.get_dummies(pd.Categorical(ages, categories=AGES), drop_first=True).astype(float).values
    P = pd.get_dummies(W.year.values, drop_first=True).astype(float).values
    X = np.column_stack([np.ones(len(y)), A, P])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.concatenate([[0.0], b[1:1+A.shape[1]]])

def nonlinear_part(prof):
    """⚠ 只留**非线性**分量 —— 线性分量与时期/世代共线,识别不了(`#844`)。"""
    x = AGES.astype(float)
    return prof - np.polyval(np.polyfit(x, prof, 1), x)

# 先扫分辨率:每个宽度上比较「真剖面 vs 假剖面」的非线性全幅
print(f"\n=== ①a 先扫年龄分箱宽度 —— 只在**假剖面明显小于真剖面**的分辨率上下判 ===")
rg0 = np.random.default_rng(7); SWEEP = {}
for w in (1, 5, 10):
    ag = np.sort(np.unique(band(W.age.values, w)))
    def _prof(y, a, ag=ag):
        A = pd.get_dummies(pd.Categorical(a, categories=ag), drop_first=True).astype(float).values
        P = pd.get_dummies(W.year.values, drop_first=True).astype(float).values
        X = np.column_stack([np.ones(len(y)), A, P])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return np.concatenate([[0.0], b[1:1+A.shape[1]]])
    def _nl(pr, ag=ag):
        x = ag.astype(float); return pr - np.polyval(np.polyfit(x, pr, 1), x)
    tr = _nl(_prof(W[IT].to_numpy(float), band(W.age.values, w)))
    sh = _nl(_prof(W[IT].to_numpy(float), band(rg0.permutation(W.age.values), w)))
    SWEEP[w] = dict(k=len(ag), true=float(tr.max()-tr.min()), sham=float(sh.max()-sh.min()))
    SWEEP[w]["ratio"] = SWEEP[w]["sham"]/SWEEP[w]["true"]
    print(f"  {w:>2} 岁箱({len(ag)} 个哑变量):真剖面非线性全幅 **{SWEEP[w]['true']:.4f}** · "
          f"假剖面 **{SWEEP[w]['sham']:.4f}** ⇒ 假/真 = **{SWEEP[w]['ratio']:.3f}**"
          f"{'  ⇒ **不可用(假的几乎一样大)**' if SWEEP[w]['ratio'] > 0.5 else '  ⇒ **可用**'}")
usable = [w for w in SWEEP if SWEEP[w]["ratio"] <= 0.5]
print(f"  ⇒ 可用的分辨率:**{usable or '无'}** ⇒ 本轮采用 **{BANDW} 岁箱**")
prof = age_profile(W[IT].to_numpy(float), W.aband.values)
nl = nonlinear_part(prof)
adj = dict(zip(AGES, nl))
W["y_raw"] = W[IT].astype(float)
W["y_adj"] = W.y_raw - W.aband.map(adj).astype(float)
print(f"\n=== ① 年龄剖面:只扣**可识别的那一半** ===")
print(f"  年龄剖面全幅 **{prof.max()-prof.min():+.4f}** · 其中**非线性分量全幅 {nl.max()-nl.min():+.4f}** "
      f"(占 {100*(nl.max()-nl.min())/(prof.max()-prof.min()):.0f}%)")
print(f"  ⚠ **线性那一半没有被扣,而这不是保守,是它识别不了**(`#844`:共线的正是线性分量)。")

def cells(col):
    ys = {}
    for y, g in W.groupby("year"):
        a = g[g.stratum == 1][col].to_numpy(float); b = g[g.stratum == 0][col].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    return ys
def dep(ys, dc=1990, Bv=B, seed=845):
    S = sorted(ys); span = S[-1]-S[0]
    grp = {}
    for y in S: grp.setdefault((y//10)*10, []).append(y)
    yy = grp[dc]; g = lambda Z, y: Z[y][0].mean()-Z[y][1].mean()
    ref = (g(ys, S[-1])-g(ys, S[0]))*(yy[-1]-yy[0])/span
    obs = (g(ys, yy[-1])-g(ys, yy[0])) - ref
    rg = np.random.default_rng(seed); out = np.empty(Bv)
    r = lambda a: a[rg.integers(0, len(a), len(a))]
    for i in range(Bv):
        Z = {y: (r(ys[y][0]), r(ys[y][1])) for y in (S[0], S[-1], yy[0], yy[-1])}
        out[i] = (g(Z, yy[-1])-g(Z, yy[0])) - (g(Z, S[-1])-g(Z, S[0]))*(yy[-1]-yy[0])/span
    lo, hi = np.quantile(out, [.025, .975])
    return dict(obs=float(obs), lo=float(lo), hi=float(hi), sd=float(np.std(out)),
                excl=bool(lo > 0 or hi < 0))

raw, adjr = dep(cells("y_raw")), dep(cells("y_adj"))
print(f"\n=== ② 九十年代那一格:扣年龄前 vs 扣掉可识别的那一半后 ===")
print(f"  未扣年龄:**{raw['obs']:+.4f}** [{raw['lo']:+.4f},{raw['hi']:+.4f}] "
      f"{'**排零**' if raw['excl'] else '含零'}(SD {raw['sd']:.4f})")
print(f"  扣非线性年龄:**{adjr['obs']:+.4f}** [{adjr['lo']:+.4f},{adjr['hi']:+.4f}] "
      f"{'**排零**' if adjr['excl'] else '含零'}(SD {adjr['sd']:.4f})")
delta = adjr["obs"]-raw["obs"]
print(f"  ⇒ 变化 **{delta:+.4f}**,而该格自身自助 SD = **{raw['sd']:.4f}** ⇒ "
      f"**{abs(delta)/raw['sd']:.2f} 个 SD**")

print(f"\n=== ③ 控制(⓪ 功效闸排在所有关于人的判断之前)===")
# ⓪ 功效闸:植入一个已知大小的缩减,必须捞得回来
PLANT = 0.5*abs(raw["obs"])
ysp = cells("y_raw"); Sp = sorted(ysp)
grp = {}
for y in Sp: grp.setdefault((y//10)*10, []).append(y)
y1 = grp[1990][-1]
ysp2 = {y: (a.copy(), b.copy()) for y, (a, b) in ysp.items()}
ysp2[y1] = (ysp2[y1][0] + PLANT, ysp2[y1][1])          # 只动虔诚层尾年 ⇒ 差距上移
pl = dep(ysp2)
got = pl["obs"]-raw["obs"]
print(f"  ⓪ **功效闸**:在九十年代尾年给虔诚层植入 +{PLANT:.4f} ⇒ `dep_gap` "
      f"{raw['obs']:+.4f} → **{pl['obs']:+.4f}**(动 **{got:+.4f}**,预期 **+{PLANT:.4f}**;"
      f"该年非全程端点 ⇒ 参照不变,`#838`③ 已验过整张网格)")
print(f"     ⇒ 捞回比例 **{got/PLANT:.3f}** · 相对该格 SD **{abs(got)/raw['sd']:.2f} 个 SD** ⇒ "
      f"{'**有功效**' if abs(got)/raw['sd'] > 2 else '**没功效 ⇒ 不看结果**'}")
# 负控:假年龄剖面(打乱年龄标签后估出的剖面)——同样缩小方差,却不含真实年龄信息
rg = np.random.default_rng(7)
sham = W.copy(); sham["aband"] = band(rg.permutation(W.age.values), BANDW)
sprof = nonlinear_part(age_profile(sham[IT].to_numpy(float), sham.aband.values))
sadj = dict(zip(AGES, sprof))
W["y_sham"] = W.y_raw - W.aband.map(sadj).astype(float)
shamr = dep(cells("y_sham"))
print(f"  负控(**假年龄剖面**:打乱年龄标签后估出的剖面,同样缩方差却不含真实年龄信息)")
print(f"     ⇒ **{shamr['obs']:+.4f}** [{shamr['lo']:+.4f},{shamr['hi']:+.4f}] · "
      f"相对未扣的变化 **{shamr['obs']-raw['obs']:+.4f}**")
print(f"     ⚠ **真剖面与假剖面的差 = 年龄真正解释掉的部分 = "
      f"{(adjr['obs']-raw['obs'])-(shamr['obs']-raw['obs']):+.4f}**")
print(f"     ⚠ 假剖面非线性全幅 {sprof.max()-sprof.min():.4f} vs 真剖面 {nl.max()-nl.min():.4f}"
      f" —— **假的必须明显更小,否则「真剖面」也只是噪声**")

G = Gate("#845 · 只扣可识别的那一半年龄")
G.asserted("⓪ **功效闸(排在所有关于人的判断之前,`#835` 的做法)**:植入一个等于观测量一半的缩减,"
           "必须捞回 ≥2 个自助 SD",
           bool(abs(got)/raw["sd"] > 2 and abs(got/PLANT - 1) < 0.15),
           f"捞回 {got:+.4f}/预期 +{PLANT:.4f} = {got/PLANT:.3f} · {abs(got)/raw['sd']:.2f} 个 SD",
           kind="control")
G.asserted("① 识别控制(`#844` 的直接后果):**只扣年龄剖面的非线性分量** —— "
           "线性分量与时期/世代共线,**扣它就是在扣一个识别不了的东西**;非线性全幅必须 >0",
           bool(nl.max()-nl.min() > 1e-6),
           f"剖面全幅 {prof.max()-prof.min():.4f} · 非线性全幅 {nl.max()-nl.min():.4f}", kind="control")
G.asserted("② 负控(跑前写下的最强混淆):**扣减本身会缩小方差,让「仍排零」显得更容易** ⇒ "
           "用**假年龄剖面**(打乱年龄标签)做同样的扣减,它同样缩方差却不含真实年龄信息;"
           "**假剖面的非线性全幅必须明显小于真剖面**,否则真剖面也只是噪声",
           bool(sprof.max()-sprof.min() < 0.5*(nl.max()-nl.min())),
           f"假 {sprof.max()-sprof.min():.4f} vs 真 {nl.max()-nl.min():.4f}", kind="control")
G.asserted("③ kill(预注册):「九十年代那条缝不是年龄剖面造出来的」要成立,"
           "需扣掉非线性年龄后**仍排除零**",
           bool(adjr["excl"]), f"扣后 {adjr['obs']:+.4f} [{adjr['lo']:+.4f},{adjr['hi']:+.4f}]",
           kind="kill", yardstick="扣年龄后 `dep_gap(1990s)` 自己的 95% 自助区间",
           yardstick_noise=float(adjr["sd"]))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
true_share = ((adjr["obs"]-raw["obs"])-(shamr["obs"]-raw["obs"]))/abs(raw["obs"])
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif abs(delta) < raw["sd"]:
    VERD = (f"**A 九十年代那条缝不是年龄剖面造出来的。** 扣掉**可识别的那一半**年龄之后,"
            f"`dep_gap` 从 {raw['obs']:+.4f} 动到 {adjr['obs']:+.4f} —— "
            f"**只动了 {abs(delta)/raw['sd']:.2f} 个自助 SD**,而且仍然排除零。\n"
            f"  ⚠ 扣掉**假**年龄剖面也会动 {shamr['obs']-raw['obs']:+.4f} ⇒ "
            f"**年龄真正解释掉的只有 {true_share:+.1%}**。\n"
            f"  ⇒ **一句关于人的话:九十年代信教与不信教之间那条张开的缝,不是「老人少了」——\n"
            f"  把年龄剖面里能识别的那一半扣干净,缝还在原地。**")
else:
    VERD = (f"**B 一部分是年龄,一部分不是 —— 报比例,不报是/不是。** "
            f"扣后 {adjr['obs']:+.4f} [{adjr['lo']:+.4f},{adjr['hi']:+.4f}]"
            f"{'(仍排零)' if adjr['excl'] else '(**含零**)'},"
            f"动了 {abs(delta)/raw['sd']:.2f} 个 SD;而**假**剖面也动 {shamr['obs']-raw['obs']:+.4f}"
            f" ⇒ **年龄真正解释掉的是 {true_share:+.1%}**。")
print(VERD)
print(f"\n⚠ **这一轮能问的比 `#834`① 问的小,而这不是保守**:线性那一半年龄与时期/世代共线,"
      f"**识别不了就不能扣** —— 所以本轮回答的是「**非线性**年龄能不能解释掉它」,不是「年龄能不能」。")
json.dump(dict(n=len(W), raw=raw, adjusted=adjr, sham=shamr, delta=delta,
               delta_in_sd=abs(delta)/raw["sd"], true_age_share=true_share,
               band_width=BANDW, resolution_sweep=SWEEP,
               profile_range=float(prof.max()-prof.min()), nl_range=float(nl.max()-nl.min()),
               sham_nl_range=float(sprof.max()-sprof.min()),
               power=dict(plant=PLANT, got=got, ratio=got/PLANT, in_sd=abs(got)/raw["sd"]),
               admissible=adm, verdict=VERD, gate_ok=G.verdict()),
          open(OUT/"only_identified_half.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'only_identified_half.json'}")
