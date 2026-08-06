"""tools/wall_audit.py — 我写在「结构上做不到」里的东西,有多少从来没有被量过?

造于 `E02·A243·R626`(`#581` 的 NEXT)。行动类型:**PRODUCTION**。

**由来(`#581c`):** 「三道性题只给三个配对 ⇒ 分辨率远低于家庭侧」这句话,
我在五份报告的「结构上做不到」里写过,**从来没有量过它**。
量了之后:**1.4 倍**,而效应是分辨率的 **8.3 倍** —— **它不是一堵墙,是一个系数。**
⚠ **一堵墙让停下来显得理所当然,所以它从不被审计。**

**这件工具做的事:** 从账本里机械抽出全部「结构上做不到 / IMPOSSIBLE」条款,按**可量性**分三类:
  `已量`         —— 该条款自己带着一个数
  `可量而未量`   —— 它谈的是**本地数据的一个可计算性质**(配对数 · n · 分辨率 · 覆盖 · 题数),却没有数
  `真结构不可量` —— 它要的东西**不在本地**(第二组编码者 · 干预 · 跨站点 · 系统发生树 · 联网 · 对抗 agent)

⚠ P6 代理账:
  PROPERTY   这条限制到底能不能量
  PROXY      条款文本里有没有数字,以及它提到的是「本地可算的性质」还是「不在本地的东西」
  IMPLICATION 只有一个方向可靠:**判「可量而未量」 -> 值得去量一次**(可靠);
             判「真结构不可量」**不证明**它真的不可量 —— 它只说这段文本里没有可量的抓手。
  SAFE SIDE  只把「可量而未量」当作行动项;从不用「真结构不可量」去关闭一个问题。
"""
import pathlib, re

# ⚠ 自检当场抓到第一版的弱点:『分辨率远低于家内的 **21 对**』带着**对方**的数,
#   于是被判『已量』—— 而那句话恰恰**没有量过两者的比**。
#   收紧:判『已量』须有**小数**,或**「倍 / %」这类比较量**。一个纯计数(21 对)不算量过。
NUM = re.compile(r'\d+\.\d+|\d+\s*(?:倍|%)')
LOCAL = re.compile(r'配对|对数|分辨率|MDE|精度|展布|n\s*[=<>]|题数|覆盖|基数|样本|格|池|k\s*[=<>]|区间')
ABSENT = re.compile(r'第二组|第二份|另一(?:批|组|份)|干预|因果|跨站点|跨国|系统发生树|Galton|联网|'
                    r'对抗\s*agent|unchallenged|复核|码本(?:不在|缺)|无码本|未发布|新采集|抽样框')


def extract(led_path="RETRACTIONS.md"):
    """返回 [(entry, clause)] —— 全部『结构上做不到 / IMPOSSIBLE』条款。

    ⚠ 抽取器改了三次,每次都是被对象改的(`#582b` 记了这个次数):
      ① 第一版取整行 -> 把 NEXT 里『读 #296 的 IMPOSSIBLE 栏』这类**引用**也算进来(`#546` 同型);
      ② 第二版要求该词在**行首** -> 只剩 2 条,因为真实形态是
         `**576e · 结构上做不到。** CLOSURE 不分离世界 · 仅女性 · …`,词在编号之后;
      ③ 第三版(本版):按**真实形态**取 —— 匹配标题,取标题**之后**的正文,按 `·` 切成条款。
    """
    led = pathlib.Path(led_path).read_text()
    marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
    HEAD = re.compile(r'(?:\*\*)?\d+[a-z]?\s*·\s*结构上做不到[^*]*(?:\*\*)?|^IMPOSSIBLE[::]')
    out = []
    for i, (n, s0) in enumerate(marks):
        body = led[s0:(marks[i + 1][1] if i + 1 < len(marks) else len(led))]
        for ln in body.splitlines():
            m = HEAD.search(ln)
            if not m: continue
            tail = ln[m.end():]
            for part in re.split(r'\s·\s', tail):
                part = part.strip(' *·。')
                if len(part) > 6 and not re.search(r'读\s*`?#|沿用|见\s*`?#', part):
                    out.append((n, part))
    return out


def classify(text):
    if NUM.search(text): return "已量"
    if ABSENT.search(text): return "真结构不可量"
    if LOCAL.search(text): return "可量而未量"
    return "判不了"


def self_test():
    """两个已知答案。#521 一个编码团队 -> 真结构不可量;#536e 三对 -> 可量而未量。"""
    a = "一个编码团队 ⇒ 无跨仪器复制,需要第二组编码者"
    b = "三道性题只给三对 ⇒ 性内的分辨率远低于家内的 21 对"
    ra, rb = classify(a), classify(b)
    # ⚠ b 里含「21 对」这个数 —— 它会被 NUM 命中。这正是判据的一个真实弱点:
    #   一句「A 远低于 B」即使带着 B 的数,也**没有量 A 与 B 的比**。
    #   所以 NUM 必须排除「只出现在被比较的另一方上的数」——用一个更严的判据:
    #   要判「已量」,文本里必须有**比较词 + 数**(倍/%)或一个**该条款自己的量**。
    print(f"  #521 型 -> {ra}(须『真结构不可量』){'✅' if ra=='真结构不可量' else '⛔'}")
    print(f"  #536e 型 -> {rb}(须『可量而未量』){'✅' if rb=='可量而未量' else '⛔'}")
    return ra == "真结构不可量" and rb == "可量而未量"


if __name__ == "__main__":
    import sys, collections, json
    print("=== 判据自检 ===")
    ok = self_test()
    if not ok:
        print("⛔ 判据未通过对照 —— 不输出任何一行(P5★)"); sys.exit(2)
    rows = extract()
    cnt = collections.Counter(classify(t) for _, t in rows)
    print(f"\n=== 全账本『结构上做不到』条款 {len(rows)} 条 ===")
    for k, v in cnt.most_common(): print(f"  {k:12s} {v:4d}")
    todo = [(e, t) for e, t in rows if classify(t) == "可量而未量"]
    print(f"\n=== 「可量而未量」逐条({len(todo)} 条,它们是下一批该跑的轮次)===")
    for e, t in todo[:25]: print(f"  #{e}  {t[:78]}")
    if len(todo) > 25: print(f"  …还有 {len(todo)-25} 条")
    out = pathlib.Path("E02_condemnation_is_not_rarity/A243_audit_the_walls/R626_measurability/results")
    out.mkdir(parents=True, exist_ok=True)
    json.dump(dict(total=len(rows), counts=dict(cnt), measurable_but_unmeasured=todo,
                   safe_side="只把『可量而未量』当作行动项;从不用『真结构不可量』去关闭一个问题",
                   unchallenged=True), open(out / "wall_audit.json", "w"), indent=1, ensure_ascii=False)
    print(f"\nwrote {out/'wall_audit.json'}")
