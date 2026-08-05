"""页面上「广度」是否点名 —— 自带对照的扫描(`#422`,由 `#421` 的 NEXT 要求)。

`#420a`:两种广度**可以反号**(同一预测量 −0.076 vs +0.049)。
**⇒ 任何一句没点名的「广度」都可能被读成它的反面。** 这不是修辞问题。

⚠ **用前必须跑 `controls()`**(`P5★` + `#374b`):
一个从未开火过的扫描器,它的每一个「都点名了」都是沉默,不是无罪。
⚠ **覆盖率必须和结论一起报**(`#376b`):正则找不到的不是「没问题」,是「没看」。
"""
import re

# 「广度」的泛称(需要点名的)
GENERIC = re.compile(r'广度|breadth|勾选数|pick count|category count|类别数|计数')
# 具体的点名(五个量之一;来自 `#421a` 的列举,不是猜)
SPECIFIC = re.compile(
    r'性行为计数|sex[- ]act|Totalsexacts|'
    r'恋物类别|fetish[- ]categor|totalfetishcategory|'
    r'起始类别|onset[- ]categor|\bncat\b|'
    r'块覆盖|block coverage|COVB|'
    r'总勾选项|total picks|PICKS')


def scan(text, window=120):
    """返回每一处泛称的 (位置, 原文片段, 是否已点名)。

    ⚠ 原文**逐条返回**,因为正则是子串匹配的近亲(`#374b`):
    错配只有把原文摆出来才看得见。
    """
    out = []
    for m in GENERIC.finditer(text):
        a, b = max(0, m.start() - window), min(len(text), m.end() + window)
        ctx = text[a:b].replace('\n', ' ')
        out.append(dict(pos=m.start(), term=m.group(0), ctx=ctx,
                        named=bool(SPECIFIC.search(ctx))))
    return out


def coverage(text):
    """(泛称出现数, 页面总字数) —— 分母让「都点名了」无法冒充「整页干净」。"""
    return len(GENERIC.findall(text)), len(text)


def controls():
    """正/负对照。返回 dict,全 True 才可用。"""
    # 正对照:一句**已点名**的话必须被判为 named
    POS = "而认可的**性行为计数**与恋物类别计数只相关 +0.26,所以这里的广度是指后者。"
    # 正对照 2(英文)
    POS2 = "the endorsed-sex-act count is the outlier, so breadth here means fetish categories."
    # 负对照 1:一句**没点名**的话必须被判为 not named
    NEG = "控制了广度之后,这个系数仍然存在。"
    # 负对照 2:一段与广度无关的文字必须**一条都不匹配**
    NONE = "羞耻与疗愈这两维在这八个量上是相加的,没有哪个量专门属于那一格。"
    s1 = scan(POS); s2 = scan(POS2); s3 = scan(NEG); s4 = scan(NONE)
    return dict(
        pos_named=bool(s1) and all(x['named'] for x in s1),
        pos2_named=bool(s2) and all(x['named'] for x in s2),
        neg_unnamed=bool(s3) and not any(x['named'] for x in s3),
        none_no_match=(len(s4) == 0),
        _detail=dict(n_pos=len(s1), n_pos2=len(s2), n_neg=len(s3), n_none=len(s4)))
