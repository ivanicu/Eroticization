# R686 · 她对家庭的七个判断,是一块东西吗

**类型:FRONTIER**。**判决:W2 —— 有一小块 + 三个散题。而它推翻了 `#649` 上一轮刚写下的 4.7 倍。**

## 硬规则①又抓到一次「变量名不是测量」

`chsuppor` **不是**「子女抚养」。题干是 **"Okay for unmarried woman to have and raise a child"**。

| 题 | 题干 |
|---|---|
| `okcohab` | A young couple should **not** live together unless married(**反向措辞**) |
| `chcohab` | OK for cohab couple to have and raise children |
| `chsuppor` | Okay for unmarried woman to have and raise a child |
| `prvntdiv` | Living together before marriage may prevent divorce |
| `staytog` | Divorce best solution when cannot work out marriage |
| `chunless` | People can't be really happy unless they have children |
| `marrfail` | Marriage has not worked out for most people R knows |

## 结果:最大连通子块 = 4 题,而它们是同一件事

**`okcohab` · `chcohab` · `chsuppor` · `prvntdiv`** —— 彼此相关 **0.32–0.40**。
**不是「家庭」,是「没结婚的人算不算一个家」。**

掉出来的三道各问各的,而 **`marrfail` 根本不是一个道德判断,是一句关于她身边的人的描述。**

**G4 规格曲线是平的**:阈 0.20 / 0.25 / 0.30 / 0.35 **全部给出 4**。
**正对照**:同一套算法作用在性三题上给出 **3** ✅ · **安慰剂**:打乱行后回到 **1** ✅
**方向对齐**:恰好翻了 `okcohab` 一道,**而它的题干本来就是反向措辞** —— 对齐步骤起作用了。

## 而这推翻了 `#649` 的 4.7 倍

| 对比 | 比值 | 判 |
|---|---|---|
| 性三题 vs **家庭七题** | **4.7 倍** | ⛔ **块选错了** —— 七题不是一个领域 |
| 性三题 vs 同居四题(全样本) | 1.10 倍 | ⚠ 在结果上选择 |
| **性三题 vs 同居四题(留出半样本)** | **1.10 倍** | **✅ 可报** |

**选择偏倚实测 −0.0019**:51 次「半样本选块、另半样本评」里选出的**都是同样这四题**(大小范围 [4,4]),
留出半样本上的块内中位 **+0.3750**(全样本 +0.3731)。

**归一后:性三题 +0.4679 · 同居四题 +0.4257 ⇒ 1.10 倍。**

> **「性道德比家庭道德更紧」这句话不成立。**
> **它看起来成立,只是因为「家庭道德」那一侧被当成了七道题的一个袋子,而袋子里有三道各走各的。**

## 范围(η)

**人群**:NSFG 2011–2013 女性回应者,n ≈ 5,590 · **仪器**:IH 系列,五级,同一批人同一份问卷 ·
**基线**:`#649` 的性三题 +0.4680 · **区间**:51 次半样本划分,3 个种子。

## 本轮结构上做不到的(不写 planned)

单波 ⇒ 无干预、非因果 · 七题是 NSFG 的设计不是我的,**「家庭道德」可能还有它没问的部分** ·
子块的**内容**是我读题干读出来的,不是数据说的 · `[unchallenged]`
