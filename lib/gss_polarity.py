"""GSS 容忍题组的**编码极性**,从码本查出来,放在一个所有轮次 `import` 的地方。

⚠⚠ **这个文件存在的理由,不是「方便」,是一次三度复发的缺陷:**

| 轮 | 发生了什么 |
|---|---|
| `#661` | 点名过「码 1 在不同题里意思不同」 |
| `#679` | 照样犯了,两条结论方向全反 |
| `#680` | 查码本、**把整张极性表写进账本**、并在代码里写下 `c != "colcom"` |
| `R123` · `R127` · `R151` · …(约 180 轮) | **都带着 `c != "colcom"`,一次没错** |
| **`#866`** | **从头重写了题目构造,那个例外静静地消失了** |
| `#868` | 我把它**当成新发现**重新查了一遍码本 |

**机制,而它比「colcom 是反的」这件事本身值钱得多:**
**那条守则活在一个被逐轮拷贝的 lambda 里。它靠拷贝续命,所以只要有一轮从头写,它就断了。**
**一条只靠拷贝存活的规则,离灭绝永远只差一次重写。**
⇒ 于是它搬到这里:**被 `import`,不被拷贝。**

**极性(从 `convert_categoricals=True` 的值标签直接读出,不是记忆):**
  `spk*`  1 = allowed to speak      · 2 = not allowed          ⇒ 宽容在**低**码
  `col*`  4 = allowed to teach      · 5 = not allowed          ⇒ 宽容在**低**码
  ⚠ `colcom`  **4 = yes, FIRED**    · **5 = not fired**        ⇒ 宽容在**高**码 —— **十五题里唯一一道**
  `lib*`  1 = remove the book       · 2 = not remove           ⇒ 宽容在**高**码

**为什么 `colcom` 不一样**:其余四靶问的是「**该不该允许他教书**」,
而共产主义者那一道问的是「假设他正在教书,**该不该把他开除**」——
**题干本身是反的,而取值编码(4/5)是平行的。**
⇒ **名字平行、码平行、意思相反。** 这正是硬规则①里最难看见的那一种。
"""

# 宽容(permissive)方向:+1 = 码越大越宽容 · −1 = 码越大越不宽容
PERMISSIVE_SIGN = {
    "spk": -1,   # 1=allowed → 低码宽容
    "col": -1,   # 4=allowed → 低码宽容
    "lib": +1,   # 1=remove  → 高码宽容
}
# 组内例外:键是**完整变量名**,值是它自己的方向,覆盖上面的前缀规则
ITEM_EXCEPTIONS = {
    "colcom": +1,   # 4=yes FIRED · 5=not fired ⇒ 高码宽容(`#680` 查的码本,`#868` 复查同一结论)
}
VALID_CODES = {"spk": (1, 2), "col": (4, 5), "lib": (1, 2)}


def permissive(series, name):
    """把一列原始 GSS 码变成 **0/1 的宽容指示**(1 = 宽容),缺失保持 NaN。

    `name` 必须是**完整变量名**(如 `colcom`),因为例外是逐题的,不是逐前缀的。
    """
    pre = name[:3]
    if pre not in PERMISSIVE_SIGN:
        raise ValueError(f"不认识的题目前缀:{name!r} —— 极性未知的题**不许**猜,请先查码本")
    lo, hi = VALID_CODES[pre]
    s = series.where(series.isin([lo, hi]))
    sign = ITEM_EXCEPTIONS.get(name, PERMISSIVE_SIGN[pre])
    return (s - lo) if sign > 0 else (hi - s)


def refusal(series, name):
    """1 = 不放行(拒绝)。`1 − permissive`。"""
    return 1 - permissive(series, name)


def audit(df):
    """自检:对一批已加载的原始列,把每一题的**宽容率**打出来。

    ⚠ 它不判断对错 —— 它只让极性错误**看得见**:一道极性反了的题,
    它的宽容率会和同组其余题**系统性地相反**(`colcom` 未修正时 `spk~col` 相关是 −0.509,
    修正后同组全部为正)。
    """
    out = {}
    for c in df.columns:
        if c[:3] in PERMISSIVE_SIGN:
            v = permissive(df[c], c)
            out[c] = (float(v.mean()), int(v.notna().sum()),
                      "⚠ 例外" if c in ITEM_EXCEPTIONS else "")
    return out
