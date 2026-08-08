"""E02·A210·R570 — 一个错的锚,有多大概率照样通过?

`#525` 的 NEXT,**一轮打包多个操作**(Ivan 2026-08-05:「多个操作一起算一轮」)。
本轮同时做:① 识别性(短语锚的唯一性)② 主检验(碰撞率)③ 对照(篡改必须报错 + 未篡改必须通过)
④ 规格曲线(数字锚 × 短语锚 × 中英两页,逐格)⑤ 页面兑现(方法说明已写)。

行动类型:**FRONTIER** —— 两个世界在「这个仪器到底证明了什么」上不同,而不是参数不同。

G1 ESTIMAND(先于方法):
  **碰撞率 = P(一个随机选中的、错误的账本条目,仍然含有该锚所引的文字)**
  分母 = 账本全部条目 − 1(正确的那条);分子 = 其中含该文字的条目数。
  ⚠ 这是一个**关于仪器的量**,不是关于世界的量 —— 按 RULE 2 明写。

WORLDS(本体不同,不是参数不同):
  W-VERIFY  锚**接近验证**:碰撞率 ≈ 0 ⇒ 通过几乎就意味着「引对了条目」
  W-FALSIFY 锚**只能证伪**:碰撞率显著 > 0 ⇒ 通过什么也不意味着,只有失败有信息
  W-MIXED   两类锚不同:数字锚碰撞高、短语锚碰撞低 ⇒ **仪器的强度是可选择的**,
            那么页面上该写的不是「锚只能证伪」而是「这一类锚只能证伪」

PREDICTION MATRIX(粗数,形状是重点):
  数字碰撞率        短语碰撞率        ⇒
  ≈0               ≈0               W-VERIFY 存活,页面的方法说明**写错了**,要改
  高               高               W-FALSIFY 存活
  高               ≈0               W-MIXED 存活 ⇒ 应当把数字锚换成短语锚
⚠ **W-VERIFY 的正面结局是我不欢迎的**(它意味着我刚写上页面的那句话是错的)——
  这正是 frontier §3 BASIN RULE 要求的那种步:**设计一个我不希望它为真的阳性结局**。

CONTROLS(G2):
  正对照:**篡改任一锚,检查必须报错**(g=0 时必须不通过 —— 未篡改时必须全过)
  负对照/安慰剂:**一个不存在于任何条目的字符串**,碰撞率必须恰为 0
  ⚠ 「这个零该不该是零?」—— 该。所以是 negative_control,不是 offset_control。
KILL(条件式,预注册):
  if 篡改被抓 and 安慰剂碰撞率 == 0:
      max(碰撞率) < 0.01 -> W-VERIFY;>= 0.01 -> W-FALSIFY 或 W-MIXED(看两类是否分离)
  else: UNVERIFIED
IMPOSSIBLE(结构上做不到,不标「计划中」):
  · 碰撞率只说「错锚能否通过」,**不说锚指的条目是否真的在讲这件事** —— 语义正确性无法自动验,
    这正是 `#524d` 的第一个洞,本轮**不解决它,只给它一个数**
  · 单一账本 ⇒ 无跨站点复制 · 无干预 ⇒ 非因果 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
LED = (ROOT / "RETRACTIONS.md").read_text()
marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', LED, re.M)]
BODY = {n: LED[s:(marks[i+1][1] if i+1 < len(marks) else len(LED))] for i, (n, s) in enumerate(marks)}
print(f"账本条目数 N = {len(BODY)}  (#{min(BODY)}–#{max(BODY)})")

# ---- 规则①:锚是被声明的,不是被猜的。逐条打印它引的东西与所在条目
ANCH = [
    ("num",    498, "0 / 6",                    "预注册门槛 3,过 0"),
    ("num",    523, "0.20",                     "跨仪器一致性"),
    ("num",    523, "−0.019",                   "严格度斜率下端"),
    ("num",    523, "−0.027",                   "严格度斜率上端"),
    ("phrase", 498, "这条线按它自己写下的规则关闭", "做不到 ①"),
    ("phrase", 516, "结构上就是不可能的",          "做不到 ②"),
    ("phrase", 521, "零个第二性域",               "做不到 ③"),
    ("phrase", 522, "没有第二个人问过同一个问题",   "做不到 ④"),
]
print("\n=== 规则①:逐锚打印它引的文字、所在条目、以及在该条目内出现次数 ===")
for k, e, s, why in ANCH:
    inside = BODY.get(e, "").count(s)
    print(f"  {k:6s} #{e}  {s!r:34s} 条目内 {inside} 次  {why}")
    assert inside >= 1, f"锚 {s!r} 不在 #{e} 里 —— 锚是错的"

# ---- 主检验:碰撞率
print("\n=== 主检验:碰撞率 = 含有该文字的**错误**条目数 / (N−1) ===")
rows = []
for k, e, s, why in ANCH:
    others = [n for n in BODY if n != e and s in BODY[n]]
    rate = len(others) / (len(BODY) - 1)
    rows.append(dict(kind=k, entry=e, anchor=s, why=why, n_collide=len(others),
                     denom=len(BODY) - 1, rate=round(rate, 5),
                     collide_with=sorted(others)[:8],
                     inclusion=[f"账本 {len(BODY)} 条中除 #{e} 外的每一条", "纯子串匹配,无正则"]))
    print(f"  {k:6s} #{e:3d} {s[:22]:24s} 碰撞 {len(others):3d}/{len(BODY)-1}  "
          f"**{rate:.4f}** {'⛔' if rate > 0 else '✅'}  {('例:'+str(sorted(others)[:5])) if others else ''}")

by = {k: [r["rate"] for r in rows if r["kind"] == k] for k in ("num", "phrase")}
print(f"\n  数字锚 碰撞率 {[f'{x:.4f}' for x in by['num']]}   max={max(by['num']):.4f}")
print(f"  短语锚 碰撞率 {[f'{x:.4f}' for x in by['phrase']]}   max={max(by['phrase']):.4f}")

# ---- 对照
G = Gate("一个错的锚,有多大概率照样通过?(账本锚检查)")
# 正对照:篡改必须被抓;g=0(未篡改)必须不被抓
def check_all(led_text):
    bl = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led_text, re.M)]
    bd = {n: led_text[s:(bl[i+1][1] if i+1 < len(bl) else len(led_text))] for i, (n, s) in enumerate(bl)}
    return sum(1 for k, e, s, _ in ANCH if bd.get(e, "").count(s) == 0)   # 失败的锚数
fail_clean = check_all(LED)
tampered = LED.replace("零个第二性域", "零個第二性域", 1)          # 一个字改成繁体
fail_tamper = check_all(tampered)
print(f"\n=== 对照 ===\n  未篡改 -> 失败锚数 {fail_clean}(必须 0) · 篡改一个字 -> 失败锚数 {fail_tamper}(必须 ≥1)")
G.positive_control("篡改一个锚,检查必须报错", planted=float(fail_tamper), floor=0.5, spread=1e-9)
G.negative_control("g=0:未篡改时必须一个都不报",
                   null=float(fail_clean), effect=float(fail_tamper), null_spread=1e-9,
                   null_kind="不做任何篡改的同一份账本")
PLACEBO = "一个绝不会出现在任何条目里的字符串-zzq7"
pr = sum(1 for n in BODY if PLACEBO in BODY[n]) / len(BODY)
G.negative_control("安慰剂:不存在的字符串,碰撞率必须恰为 0",
                   null=pr, effect=max(by["num"]), null_spread=1e-9, null_kind="字面不存在的串")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['kind']}#{r['entry']}:{r['anchor'][:10]}":
                             dict(n=r["denom"], **r) for r in rows}, "8 个锚,每格 n = 分母")
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{r['kind']}#{r['entry']}:{r['anchor'][:10]}": r for r in rows})

print("\n" + "=" * 72)
ok = fail_tamper >= 1 and fail_clean == 0 and pr == 0.0
if ok:
    mx_n, mx_p = max(by["num"]), max(by["phrase"])
    if max(mx_n, mx_p) < 0.01:
        world, verdict = "W-VERIFY", "两类锚碰撞率都 <1% -> **锚接近验证**,页面上「只能证伪」那句要改"
    elif mx_n >= 0.01 and mx_p < 0.01:
        world, verdict = "W-MIXED", (f"**数字锚碰撞 {mx_n:.4f},短语锚碰撞 {mx_p:.4f} -> 两类分离**:"
                                     f"「只能证伪」只对**数字锚**成立;短语锚强得多")
    else:
        world, verdict = "W-FALSIFY", f"两类都碰撞(max {max(mx_n,mx_p):.4f}) -> **锚只能证伪**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:碰撞率只问「错锚能否通过」。"
          "一个**没有碰撞、却引错了条目**的锚 —— 引了一条恰好含这句话、但讲的是别的事的条目 —— "
          "它照样通过,而碰撞率对此完全沉默。**语义正确性仍然不可自动验(`#524d` 第一个洞)。**")
else:
    world, verdict = "UNVERIFIED", f"控制未齐:clean={fail_clean} tamper={fail_tamper} placebo={pr}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(n_entries=len(BODY), rows=rows, by_kind=by, world=world, verdict=verdict,
               controls=dict(fail_clean=fail_clean, fail_tampered=fail_tamper, placebo_rate=pr),
               estimand="P(错误条目仍含该锚文字)", instrument="纯子串匹配 over RETRACTIONS.md",
               impossible=["语义正确性不可自动验", "单一账本,无跨站点复制", "无干预,非因果"],
               unchallenged=True), open(OUT / "anchor_collision.json", "w"), indent=1)
print(f"\nwrote {OUT/'anchor_collision.json'}")
