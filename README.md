# The Eroticization Operator

**"Sexual" — is it a content category the brain detects, or a value the brain assigns to ordinary
world-representations?** 189 self-attacking rounds against one public dataset, asking what can be
settled before anyone collects anything new. **144 numbered ledger entries** — [RETRACTIONS.md](RETRACTIONS.md)
carries the history; this page carries only what currently stands.

Private curiosity project. Not written for publication and not seeking any.

Object: the **Big Kink Survey** public subsample — [Zenodo 10.5281/zenodo.18625141](https://zenodo.org/records/18625141),
15,503 respondents × 365 columns, aggressively binned, demographically stripped, noise-injected.
Data is not committed; `E01…/A01…/R01_schema/run.py` re-downloads it.

## The three models being separated

Ivan's framing, kept verbatim because the whole repository is scoped by it:

| | |
|---|---|
| **A** dedicated sexual-content system | `Y_i(s) = f_i(h_sex(s))` — "sexual" is a content category like *face* in vision |
| **B** erotic valuation of ordinary representation | `v_i(s,c,t) = w_i(c,t)^T h(s)` — ordinary scene, a readout applied to it |
| **C** recursive mix | value re-enters and reshapes the representation it acted on |

---

---

# What this says about people

Each row is one sentence about people, at the strength it was actually measured. The ledger entry is
the receipt.

| | | strength |
|---|---|---|
| **性欲有一个发育顺序,先具体后关系,差 2–3 年** | 外观 14.0 · 身体部位 14.4 · 衣物 14.7 → 权力动态 16.8 · 束缚 16.9 · 精神改变 17.0。用它猜任意两个兴趣谁先谁后达 **66.852 ± 0.191**,而"全人群共享一个顺序"的理论上限是 66.5% —— 这个时间表几乎榨干了共享顺序能提供的全部信息 | `#75` 8 seeds · 35,438 对 · 6,230 人;全体对上界定 **[60.5%, 66.5%]** `#63` |
| **内容与个体化估值,对认可概率的影响一样大** | 选项基率 **±22.6 pp** · 人整体率 **±16.3 pp** · 人×选项交互 **±23.7 pp**。内容在每一个"可预测性"指标上碾压(3.5–5.6×),**纯粹因为它被估计得好 179 倍**(3,228 次观测 vs 18) | `#88` `#90` |
| **"这是性内容"在问卷里问不出来** | 二元认可矩阵里,一个选项的**基率**和它的**内容**是同一个数。题目主效应有闭式解,实测 R² = **0.994**、斜率 0.988,且对基率**形状**在 23/23 块中完全不敏感 | `#69` |
| **稀有偏好是一条所有人都在上面的连续维度,不是一群特殊的人** | 分布对称加宽:上尾挑更罕见的,下尾挑更常见的,**中位数纹丝不动**,两侧各 11–13× 自身自助 sd | `#99` |
| **它是一个可靠的人格维度,而且不是"你勾了多少"** | 跨不相交块集的分半信度 **+0.432**(地板 −0.022,种植天花板 +0.832);移除勾选数后存活 **67%**,而且它与勾选数的相关(+0.608)**低于**零的(+0.719) | `#100` `#104` |
| **一个人自己那些别人很少有的兴趣,排在他自己那些大家都有的兴趣前面** | 扣掉人群共享的时间表和这个人整体的早熟之后,人内「起始年龄残差 × 类别稀有度」= **−0.0328**,题内跨人置换零 +0.0018,**8.8×**。**离散版本(最早一格 vs 他自己的曲目库)独立复现:扣掉时间表后 +0.0767,14.5×**。`#114` 的回忆偏差实际只贡献 **4%**,且对主量符号相反 | `#128` `#130` |
| **版图真正的主分界线不是"先来后到",是"物件 vs 叙事"** | 6 个连贯分割 × 11 个非性变量的整格,最大统计量零给出全族阈值 0.0560。**PC4 点亮 6 个变量,最强一格 9.0×**:一侧是衣物·性别扮演·体液·角色·脏污·变形·玩具,另一侧是年龄·惊悚·非自愿·生物·乱伦·温柔·残暴·怀孕·神话。**五个锚在性别内部仍成立**(开放性 +0.083 · 外向性 +0.088 · 神经质 −0.064 · 无力感 −0.059,两性方向一致)。⚠ 命名是我对载荷的读法(D5);标度小(0.033–0.088)。⚠ **模型方向已撤回**:不依赖命名的结构检验(`#139`)显示 PC4 两侧**并不比别的分割更不同构**(0.8×),所以这条线与模型 A、B 都相容 | `#138` `#139` |
| **性别锚住的是"性"本身,不是其中某一族** | 早族 +0.0706 vs 晚族 +0.0816,差 **0.7×**。而开放性(+0.0593,4.1×)与无力感(−0.0441,2.8×)确实只锚住晚族 —— **但那个差不是这两族特有的**:同一矩阵的另外五个同样连贯的分割给出 0.0702–0.2106,**早/晚这条线是六个里最弱的一个**。两族在组织度与到达方式上不同,在**外部成因**上不是两样东西 | `#137` |
| **关系族同时到达,却装着好几个独立的维度** | 早/晚这条线的**结构不对称是六个连贯分割里最大的**(5.73,3.0×):晚族有效维度**多 2.42 个**、分半信度**低 0.35**、单因子占比**低 0.14**。**它在时间上是一个包,在偏好结构上不是一个东西** —— 而具体族相反:分散地到达,却更像单一的一个东西(拉平题目数后有效维度 7.64 vs 5.30)。**而它装的就是 `A02` 早就命名过的三条轴**(谁服从·谁被看·谁接受,留出相关 0.81/0.81/0.91,随机载荷地板 0.349)—— 关系族是三条近乎独立的轴的**交汇处** | `#139` `#140` |
| **早来的东西是散的,晚来的东西是一整套** | 人内置换零下的三个块:早×早 **−0.01272(27.9×)** · **晚×晚 +0.02103(23.1×)** · 晚×早 −0.00484(15.2×)。晚到的是关系族(权力动态·束缚·施虐受虐·精神改变·感官,26–31% 在 19 岁后),早到的是具体族(外观 8.6%·衣物 10.9%·身体部位 11.1%)。题目层:晚族彼此 **+0.0674**,早族彼此 −0.0075。**而且"一整套"是字面意义的**:关系族的获得年龄比同样多个随机挑的**靠拢 9.4%**(17.7×,在扣掉人群时间表之后),具体族反而**分散 9.0%**(12.1×)。**发育不是同一种东西按顺序到达,是两种不同组织度的东西以不同方式到达** —— 这是 `#75` 的时间表没有的维度。⚠ 连通度**没有拐点**(1.6×),所以"17 岁关上"只对类别数成立 | `#135` `#136` |
| **三分之一的人,版图在 17 岁就关上了;剩下的人,它是向外开的** | 29–32 岁档里 **33.3% 报告 17 岁后一个新兴趣都没有**。而在有晚期条目的人身上,晚获得的类别与他早期那些**更不相连**,配对层剥掉稀有度后保留 89%,`#114` 符号相反,**五个年龄档全部同号(5.3–10.9×)且无年龄梯度**(所以不是叙事整理)。**扩张,不是深化** —— 模型 C 不能被定域为"重排既有权重"。⚠ **只报方向,不报量级**:切点从 13.5 到 17.5 give −0.0165 → −0.0069(差 2.4 倍),且**人内中位数分割不可分辨(1.0×)—— 现象绑在绝对年龄上,不是"你序列里靠后的"** | `#133` `#134` |
| **性版图在青春期结束时就基本定型** | 29–32 岁的人,自己报告的性兴趣有 **68.4% 是 17 岁前获得的**,最晚的那一个平均在 **22.6 岁**(人内测量,不依赖横断面假设)。横断面上 15 年里类别数只从 12.4 长到 12.9(+4.5%,队列混淆,较弱)。**这给模型 C 划了时间边界:递归重塑若存在,作用窗口主要在青春期内** | `#132` |
| **人群层面,版图确实从大家都有的东西开始** | 一个人最早报告的那批兴趣,在他自己曲目库里按罕见度排落在**第 33 百分位**(49.4×);⚠ **机制 UNVERIFIED**:题目层 Spearman(稀有度, 起始年龄**中位数**)= +0.437,但中位数是一个**更差的时间表**(留出成对顺序 63.30% vs 均值 66.70%),而均值序与稀有度无关(+0.011)。左尾解释已死(−0.091)。**效应稳,机制未定** | `#130` `#131` |
| **越是最终口味罕见的人,这个提前量越大** | −0.0459;**按答题类别数卡钳 1:1 匹配后 −0.0417(3.1×,保留 91%)**,三个匹配规格 85–100% 一致,置换零 17%,种植正对照 +0.5878。`#114` 的回忆偏差实际贡献 39%。**换成离散统计量后复现:匹配后 +0.0532,5.1×** | `#128` `#129` `#130` |
| **它唯一挂得住的外部锚是性别** | +0.093(去衰减 +0.141)。人格五因素全部 \|r\| ≤ 0.056,**开放性只有 +0.023** —— 它不是一般性的求新 | `#101` `#102` |
| **人把最爱的性兴趣记得更早,约九个月** | −0.2000 年/评分标准差(**19.8× SE**),同类别内打乱零 = 效应的 **0%**,种植阶梯单调。**这是本数据集第一个被直接量出来的回忆偏差,而整个成熟时间表建立在它之上** | `#114` |
| **而且这个记忆畸变随时间加深,十五年里几乎翻倍** | 15 岁组 −0.0505 → 30 岁组 −0.0917,年龄趋势 **3.5×**;残余伪影**符号相反**,校正只会放大。**它不是在回答问卷那一刻生成的** | `#119` |
| **同类兴趣之间,先获得哪一个,预示你其余的整个偏好轮廓** | **+0.0159,6.1× 展布**,在减掉两项各自评分、它们的差、以及人均评分之后。**成对设计对"作答水平"结构上免疫** —— 它的兄弟命题("越早越在中心")正是死在这个混淆上 | `#107` `#110` `#116` `#117` |
| **而其余偏好被拉向先来的那一个** | 位移 **+0.0339,3.1×**,68 对里 46 对为正,生成式正对照单调开火。⚠ **加倍样本并没有加强它**(3.093 → 3.048),而源码里那句"加强成功"已撤回;强度的一位小数本身在自助噪声内(20 个种子:3.023–3.163) | `#118` `#142`,标 DESCRIPTION:方向不判别因果 |
| ~~**"上/下"一个词盖着三件近乎不相关的事**~~ **降级为:三条角色轴之间相关很强但不是同一个** | **`Entry 24` 早已撤回"三条近乎独立"与"有效维度 2.95/3"**。用验证过的信度阶梯(正对照:SUBSTANCE 对自身 r_true = **+1.018**;sham = +0.023)测得 **POWER–SUBSTANCE r_true = +0.605**,即**共享 37% 方差,不是 5%**。预注册判定落在中带(0.45 < 0.605 < 0.70)→ **UNVERIFIED**:既不是一个构念,也不是已确认可分。⚠ **`#141` 整条撤回** —— 它靠重跑 `A02/R10` 来"重新定价",却没读 `Entry 24`,于是把已撤回的「2.95/3」写回了这一行,并推荐引用**原始**相关(0.112),而原始相关正是**被衰减压平**的那个 | `Entry 24` · `#143` |
| **"极端"不是一条轴,是两个互斥的方向** | 卑贱污秽 ↔ 血/烧灼/武器 —— 喜欢一端的人系统性地不喜欢另一端 | A02 |
| **82.7% 说"色情给了我这个癖好",而这句话在他们的性癖结构和时间线上没有任何痕迹** | 无时序签名(1 年位移被排除,集中判别在 **8.8× 功率**下为零)· 无结构签名(错位 < 0.1 sd,极端度差 < 10%)· 它**追踪这个人整体勾了多少**(rho **+0.2922**)。⚠ 但"其中 85% 是作答风格"这一步 `#26` 已降级为 **UNVERIFIED** —— 全部题项都是情欲内容且无反向计分,"泛泛同意"与"泛泛认可情欲事物"在本 release 分不开 | A06 · `#26` |

---

# 被撤回的(留在这里,因为它们曾经被我报告过)

| 曾经的说法 | 现在 | 为什么 |
|---|---|---|
| "性是一个价值,不是一个类别"(epoch 标题) | **假,而且问题本身是两个** | 加载器第 1 行就删掉了对手;而块内基率≡内容,那个对比根本无法裁决 `#67` `#69` `#70` |
| "域一般胜过域特异,186/276" | **反转:0 一般 / 201 特异** | 那 186 是我给域特异侧的 −0.135 估计器坑 `#82` |
| "三个成分一样大(可预测性上),1.05×" | **撤回** | 减零不修正估计器,它把估计器的失败记在被估计量头上 `#86` `#87` `#88` |
| "少数人集中于稀有选项" | **措辞撤回,存在性保留** | 是对称加宽,不是子群 `#99` |
| "越早进入版图的东西越在中心" | **撤回,只剩 10%** | 那不是时间,是"什么都给高分的人" `#115` |
| "特质追踪更早的获取年龄" | **撤回** | −0.030 对零 −0.028,零是效应的 91% `#102` |
| "跨块累积的 √n 反映潜结构" | **反转** | 等预算下块越多越差;那是普通的 √N `#64` |

---

# 这套数据做不到的(每条都量过,不是猜的)

| 做不到 | 量出来的边界 |
|---|---|
| 看见 ≤30% 人携带的结构(方差解释法) | 5% 携带者在 **±50 pp** 下技能不离地板;需 ≈30% `#91` |
| 区分"一小群对某几样异常强烈"与"一条重尾连续谱" | 保边际的集中/弥散种植在 3000 交换/块下都在零抖动内;**出路是外部锚,不是更强的种植** `#122` |
| 区分记忆变模糊与故事被讲实 | 真实曲率 **0.3×**,推到 2× 需约 **48 倍**样本 `#120` |
| 直接测量交互幅度(不经模型) | **代数上不可能** —— 残差一阶矩是边际决定量,保边际零下恒等于零 `#105` |
| 数出维度 | 有诚实地板的估计器数不了;能数的地板是伪影 `#89` |
| 呈现顺序首因 | release **按字母序**导出多选答案(119 对一致性 1.0000),显示顺序已毁 `#69` |
| 区分广度与默许 | 需反向计分或强制选择题项,本 release 没有 |
| 因果方向 | 横断面,剥离过的 release |

**并且这是一条可迁移的方法学结论,不只对这个数据集**:任何在这个 release 上的分组比较,若没有按**块数**匹配,有一大半在测**问卷覆盖度**而不是性癖差异(`corr = +0.815`,9 次切分)。我没在任何已发表分析里见过有人提这一条。

---

## What each big round established (R = one iteration, one belief update, one commit)

**[R01 · the object](A01_object_and_structure)** — The item-level data is not the 68 category
ratings; it is 101 multiselect columns exploding to **1,332 options over 15,468 people**. Entry
to every block is gated on a parent rating: **P(enter | parent > 0) = 0.99**. This is
undocumented and it constrains every design downstream — naive cross-block holdout is
conditioned on liking the parent category. Person-profile split-half reliability **0.727**;
item base-rate reliability **0.999**. The noise injection did not destroy individual signal.

**[R02 · is there a shared grammar](A02_is_there_a_shared_grammar)** — Yes, and the first two
instruments could not see it. PC1-vs-PC1 gave 0.064 (retraction #2); held-out CCA gave 0.272 vs a
0.055 floor; leave-one-block-out gave **32/32 blocks positive**, median gain +0.0340 against a −0.0029 floor. **[BOUNDED — RETRACTIONS #25]** with all 32 coverage indicators in the baseline the gain is +0.0170 and 31/32 blocks stay positive; coverage is partly a mediator, so the honest statement is a bound of **[+0.017, +0.037]**. Factors learned from 31 domains predict a domain they never saw.

**[R03 · naming the coordinates](A03_naming_the_coordinates)** — Four coordinates survive a
**block** split-half (not a person split-half — the question is whether a coordinate is
recoverable from either half of the *domains*). Naming failed twice before succeeding. The
surviving three: *light restraint/toys ↔ insertive extremity and confinement* · *abjection and
filth ↔ blood, burning, weapons* (two **opposed** extremities, not one intensity axis) ·
*receiving a substance ↔ giving it*. None is the folk axis. **"Top/bottom" is one word over
three near-independent coordinates**: who submits, who is seen, who receives — disattenuated
**[DOWNGRADED — RETRACTIONS #24, re-affirmed #143]** "three near-independent axes / effective dimensionality 2.95 of 3" was withdrawn long ago. With a validated reliability ladder (positive control: SUBSTANCE against itself r_true = **+1.018**; sham +0.023), **POWER–SUBSTANCE r_true = +0.605** — they share **37% of variance, not 5%**. The pre-registered kill lands in the middle band (0.45 < 0.605 < 0.70): **UNVERIFIED**, neither one construct nor confirmed distinct. **RETRACTIONS #141 is withdrawn in full** — it re-priced the claim by re-running A02/R10 without reading the ledger entry that had already superseded it.

**[R04 · acquisition and time](A04_acquisition_and_time)** — Interests arrive on a
population-shared schedule (content-like early: appearance 14.0, body parts 14.4, clothing 14.7;
relational late: power dynamics 16.8, bondage 16.9, mental alteration 17.0). **[SURVIVES FRAMING SWAP — RETRACTIONS #50]** Within-person "acquired together" tracks "liked together" at RSA **+0.599**, and unlike the cross-domain grammar (#49) it carries predictable variance: person-level onset→preference R² = +0.0136 against a −0.0043 null, 31% of the same-domain ceiling after stripping intensity leakage
(−0.126 — **[VERIFICATION — see RETRACTIONS #16]**, published as `01-age-onset` Finding 1) and recall anchoring, and after near-synonymy is excluded (0.594 for
280 pairs sharing no content word vs 0.646 for 64 that do). **[CONFIRMED — RETRACTIONS #22]** onset carries structure preference does not: observed top residual eigenvalue 0.959 against a purpose-built rival world at 0.441 ± 0.040 (95% upper 0.532), positive control fires at injection 0.3 and stays silent at 0.0. The attached '80% of the sd' is descriptive, not tested — interests a person does *not* like together were acquired at the
same time. But the organising variable is arrival time, not coordinate membership (retraction #9).

**[R05 · group differences and the instrument](A05_group_differences_and_the_instrument)** —
**Any group comparison on this release must be block-count matched or it partly measures survey
coverage**: corr(congruence deficit, coverage gap) = **+0.815** across nine splits. This is the
transferable methodological result and it is not documented anywhere I have read about this
dataset. Drawn-vs-live-action consumers — whose content contains no real bodies — differ by **0.0204 ± 0.0265**, bounded below sex-sized (0.093). **[UNVERIFIED — RETRACTIONS #34]** the written-vs-visual deficit is unresolvable at 5 seeds (0 of 9 corpus cuts have effect > 2× seed spread, while sex is resolvable in 9/9); ~44 seeds per cell would be needed. That is a real, bounded constraint on A.

**[R06 · induction](A06_induction)** — 82.7% of the 13,530 respondents with fetishes say porn
induced ones that would not otherwise exist. That claim carries **no acquisition timing
signature** (a uniform 0.5-year shift of *all* their onsets; a 1-year shift is excluded; the
concentration discriminator is null at 8.8× the power needed for a single-interest induction)
and **no structural signature** (misfit < 0.1 sd, extremity < 10%). It tracks **breadth**
— **[CORRECTED — RETRACTIONS #26, caught by #144]** the previously stated "rho +0.2515, 85%
surviving response-style control" **does not appear anywhere in the ledger**; it came from a
pre-A09 draft of this page. Measured: **rho +0.2922 uncontrolled, +0.2523 with all indicators**,
and the "85% is response style" step is **UNVERIFIED** — every item is erotic content with no
reverse-keying, so general agreeableness and general erotic endorsement are not separable in this
release. Design consequence: retrospective self-report of induction is unusable as an outcome
measure. Phase 3 must be prospective.

**[R07 · breadth](A07_breadth)** — Quantity without shape. A person's set is **0.88%** more
concentrated in coordinate space than a size-matched base-rate set — real (t=−15.5) and
negligible. Sets are **24.2%** of the way from chance to perfect nesting. Breadth is moderately one trait (Spearman–Brown 0.557). **[UNVERIFIED — RETRACTIONS #26]** the '9–13% response style' figure comes from an index that is orthogonal to the POWER axis (+0.024) but correlates +0.385 with breadth itself; response style and erotic endorsement cannot be separated without balanced-keyed items, which this release lacks.

**[E01·R05 · is breadth the object](E01_sexual_as_a_value_not_a_category/A05_is_breadth_a_nuisance_or_the_object)** — Nearly domain-encapsulated. All 15 non-sexual
variables jointly: **R² = 0.012**. **[PARTLY RETRACTED — see RETRACTIONS #17]** That held only with survey progression controlled and a *gated* outcome. On the ungated `totalfetishcategory` with acquiescence controlled, childhood adversity → breadth is **r = +0.059, effect/floor 7.5** — real and small, matching the published d=0.151. Adult sexual assault, corporal punishment and upbringing remain near zero.

**[R09 · consumption](A09_consumption)** — Consumption touches both terms: it correlates with θ
(rho 0.17) *and* independently with the coordinates (**0.0439** after matching coverage, breadth
and sex — 3.6× the neuroticism reference). Direction is unavailable; this is exactly what a
prospective design would resolve.

**[R10 · additivity](A10_additivity)** — Ivan's model B requires it. The feature crosses every
substance boundary it meets (+0.24 to +0.58) **except source gender, where it is +0.017**.
Additivity is basis-dependent: in the folk basis (self/other) it fails; in a basis crossing role
with source gender it may hold. Small n on the decisive cells (3–4 pairs) — **PLAUSIBLE, not
CONFIRMED**.

---

## How the shared signal accumulates

**[RETRACTIONS #60]** `increment = 0.00723 × √(source domains)` — CV **6.4%** across n=1…31, against 11.2%
for log and 51.9% for linear, with a permuted-label null of exactly 0.0000 at every n. **[MECHANISM REVERSED — #61, #64]** blocks are *not* interchangeable (subset variance 58× seed variance), and at
**fixed total respondent-rows more blocks is worse**: n=8 gives 0.0042, n=16 gives 0.0028, n=31 gives 0.0025.
The accumulation tracks **total sample**, not block count — the √ shape is ordinary √N. Block boundaries are a
tax on the estimate, not a source of it. Price list for new collection:
0.06 needs 69 domains, 0.08 needs 122, 0.10 needs 191 — extrapolated beyond the measured range.

## Within-block vs cross-domain

**[RETRACTIONS #59]** Two largely independent contributions, 87% additive, both against two exactly-null
controls: **cross-block factors +0.0409** (fitted on 31 *other* blocks) and **within-block structure
+0.0290**, combined **+0.0606**. The aggregate cross-domain signal is the *larger* — pairwise block→block
is ~0 (#49), and it accumulates across blocks.

## The item margin

**[RETRACTIONS #57]** Fifty-seven rounds measured structure over PEOPLE. The ITEM side, measured once:
item-neighbour structure recovers **+0.0206** held-out R² over a marginals base, against the person
factors' **+0.0289** and a random-neighbour control's **+0.0006** — **71% of the person margin**, on
578,989 held-out cells. **[#58]** Fitted jointly, the two are **67% the same structure**: combined increment +0.0359 against the person margin's +0.0293, so the item side's unique contribution is +0.0066 — ratio 1.65, **below resolution**.

## Governance

| file | what it holds |
|---|---|
| [`RETRACTIONS.md`](RETRACTIONS.md) | **144 entries.** Every claim killed, scoped or corrected, with what killed it. Twelve are a later round of mine destroying an earlier round of mine |
| [`FROZEN.md`](FROZEN.md) | Lines where further computation cannot identify what it is measuring, each with its unfreeze condition |
| [`PREREGISTRATION.md`](PREREGISTRATION.md) | What the next rounds will test, with thresholds fixed in advance. Rounds r01–r52 were **not** preregistered and are labelled exploratory throughout |
| [`ADVERSARY_FORECAST.md`](ADVERSARY_FORECAST.md) | What I expect an outside challenger to overturn, written before one arrives |
| [`STANDARD_AUDIT.md`](STANDARD_AUDIT.md) | These 52 rounds scored against the campaign standard, including everything they fail |

## Layout

**E · R · r** — epoch (the ontology shifted) contains big rounds (a decision became safe) contains
sub-rounds (one belief update). One epoch, six big rounds, 54 sub-rounds. Every count is discovered,
never chosen; see `~/.claude/CLAUDE.md` §P16.

Previously described as ten campaigns, 52 rounds. Each round is a directory with `run.py`, `README.md` and `results/`.
`lib/rounds.py` maps a round name to its path, because several rounds reuse an earlier round's
loaders and that dependency is made explicit rather than hidden. Environment is a self-contained
`.venv` (system python 3.14 has no `ensurepip`; pip was bootstrapped).

## Scope, stated once

Cross-sectional · one instrument · one population (18–32, US/Canada/Europe) · aggressively
anonymised, with correlations attenuated roughly 25% by design · every measure self-report.
No causal claim in this repository is identified, and none is made.
