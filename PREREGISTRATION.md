# Preregistration

**Rounds r01–r52 were not preregistered.** They are exploratory, and every claim in
[`README.md`](README.md) carries that status. Retrofitting a preregistration onto completed work
is the thing this file exists to prevent, so nothing below refers to a round that has run.

Thresholds are fixed here before the data is touched. A threshold chosen after seeing a result is
a narrative.

---

## r53 · Prior art, and it runs before anything else

**Why first.** The CoVal programme spent 140 rounds discovering an algorithm documented in its own
dataset card. This release has a public explorer, a Zenodo record and a secondary literature, and
**no round here has checked whether any finding restates something already published.**

**Procedure.** For each of the eleven headline claims in `README.md`: search the release
documentation, the explorer's own analysis pages, and the literature. Code each claim
`NOVEL / VERIFICATION / SUPERSEDED`.

**Pre-registered consequence.** Any claim coded VERIFICATION is relabelled in `README.md` as a
verification that the object does what it says — **not** as a finding. No exceptions, including
for the three I most want to keep (the three role axes, the +0.815 coverage artifact, the
theta encapsulation).

## r54 · Multi-seed, because it is the cheapest thing absent

**Claim under test.** The eleven headline numbers are seed-stable.

**Procedure.** Re-run every round contributing a headline number at ≥5 seeds. Report
`seed_spread / |effect|` for each.

**Pre-registered kill.** Any headline whose `seed_spread / |effect| > 0.5` is downgraded to a
direction and loses its number in `README.md`.

## r55 · The norm manipulation, the one separator this release cannot provide

**Claim under test.** Retraction #14 killed `n_i` as a *trait*. It did not test the prediction the
parameter was invented for: that normalising a feature moves erotic value in **opposite
directions** for transgression-driven and feature-driven people.

**Why it is the highest-leverage new collection.** A manipulation where two subgroups move in
opposite signs is far more diagnostic than any mean shift, and it is unfakeable by demand
characteristics — a participant guessing the socially correct answer produces a uniform sign.

**Design.** Same physical scenario; the only manipulated variable is stated prevalence
("most people find this ordinary" vs "this is rare"). Outcome measured as arousal / wanting /
pleasure separately, never as one score. Between-subjects on the framing, within on the scenario.

**Pre-registered thresholds.** Support requires (a) a subgroup × framing interaction with the two
subgroup slopes of opposite sign, and (b) the interaction surviving the split defined on a
*held-out* half of each participant's scenarios. Fixed now: interaction p < .01 and both simple
slopes' 95% CIs excluding zero **in opposite directions**. A main effect with no sign split is
recorded as a refutation, not as partial support.

## r56 · Triple-blind reanalysis of the shared-grammar claim

**Claim under test.** The 0.200 demographic-adjusted cross-domain CCA.

**Procedure.** Two clean-context agents receive the data, the schema, the question in plain words,
this file, distinct seeds, and a list of files they may not open. They receive **no** estimand, no
statistic, no controls, no numbers. Designing the statistic is the task.

**Reading rule, fixed now.** All three agree → design-independent. Agree on sign, differ on size →
report the spread as the finding. **Disagree on sign → the framing is the finding**, and the
assumption they differ on becomes r57. The three will not be averaged.

---

## R359 · 有没有东西**放大**这份羞耻(2026-08-04,跑之前提交)

**为什么需要预注册**:缓冲与放大在**同一个交互系数**上是同一件事,所以这一轮
**新的不是方法,是候选集** —— 而候选集若从数据里挑就什么也不是。
`#309c` 里神经质那一格(19 格中最大的一格、无多重性校正)**不可用作本轮的理由**;
下面三个候选是**先验**的,理由写在跑之前。

**候选集(三个,新的)**
| 候选 | 先验理由(与本数据无关) |
|---|---|
| **无力感** `powerlessnessvariable` | 羞耻的核心是「我不该是这样」;一个觉得自己无力改变处境的人,更难把这种自我评价放下 |
| **神经质** `neuroticismvariable` | 负性情绪的素质性倾向,是几乎所有自我指向的负性情绪的一般性放大器 |
| **0–14 岁被打屁股** | 早年把身体当作惩罚对象的社会化,是「身体与羞耻」这条通路最直接的候选史 |

**已测过、本轮不重复计入新意**:成长期性开放度 · 关系风格 · 年龄 · 开放性(`#356`,联合零)。

**设计**:与 `#356` 同一套 —— 三个候选 × 两条路 = **6 个交互项**的联合检验,
统计量 = 去掉这 6 项后的 **R² 下降**,零 = `perm_finite` 打乱人(60 次)。
**同时报**七个变量(4 旧 + 3 新)全放的 **14 项联合检验**作为次要结果。

**预注册的判据**
- **KILL**:联合 R² 下降 > 置换零的 2×展布 -> 至少有一个在放大,下一轮定位(**本轮不定位**);
  否则 -> 报**逐变量 MDE**,并交给 guard 21 判这个零可不可发布。
- **正对照**:只在**一个**候选上植入放大,强度扫描 {0, 10%, 20%, 30%, 50%},过 guard 13;
  **逐变量各测一次 MDE**(`#311b` 的教训:不能用一个变量的 MDE 代表全体)。
- **负对照**:`perm_finite` 题内跨人打乱。
- **guard 21**:若结论是零,必须交出置换分位数 · 逐变量 MDE · 正对照灵敏度;
  「有意义的放大量」预先定为 **30%**(与 `#356` 同一口径)。
- **⚠ 本轮不报「哪一个」** —— 联合检验说不出,而事后挑最大的那个正是 `#309c` 的陷阱。

**污染声明**:我已经看过 `#309c` 里 19 个条件格的数值,其中包含神经质那两格
(RS 低 +0.0745 / 高 +0.1353)。**这构成对神经质这一个候选的污染**,
所以本轮**神经质的单变量结果不可作为独立证据**;它进入候选集的理由是先验的,
而它的**单独**读数按污染处理。另外两个候选(无力感 · 被打屁股)我没有看过与羞耻的任何数值。
