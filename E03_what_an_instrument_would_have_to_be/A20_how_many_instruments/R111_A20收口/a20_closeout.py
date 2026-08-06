"""E03·A20·R111 —— A20 收口:每一条「已推进」必须指得出轮次号与数

**类型:CLOSURE**(如实标注 —— 保护 `#664`–`#668`,不开新世界)。**A20 关弧候选。**

⚠ **最强混淆(`#668` 写死,而它是这一轮的全部约束)**:**我是那个判「已推进」的人。**
  ⇒ **每一条「已推进」必须指向一个具体的轮次号与一个具体的数**,不许写「见上」;**指不出来的退回「未推进」。**

⚠ **而 `#664` 那张表只有「跨/不跨」两档,它掩掉了一个真实的等级差** ——
  `#667` 的「跨出版物、同一第一作者」与 `#665` 的「跨独立研究组」**不是一回事**。**先把等级写死:**

  **T1 跨仪器** —— 不同的调查/档案(问法、人群、年份都不同)
  **T2 跨独立编码团队** —— 同一档案,不同研究组
  **T3 跨出版物同一作者** —— 同一档案、同一第一作者的两篇
  **T4 单项目** —— 两侧同出一个项目

G1 ESTIMAND:**T1+T2 的条数 ÷ 11**(与 `#664` 的 3/11 并列);T3 单列,不并入。
G2 CONTROLS:
  **正对照(`#668` 写死)**:`#664` 已判「做不到」的两条(以身作则 · 年代单位)**必须仍是「做不到」**;
    若变「待做」⇒ **我又漏查了,而那正是 A20 存在的理由。**
  **⑤ 的机械检查**:每条「已推进」的 `round` 与 `number` 字段**必须非空且能在账本里找到**。
G3:11 条全报。G4:{只算 T1+T2, 也算 T3} 两条规格。
KILL(条件式):if 正对照仍是「做不到」and 每条已推进都指得出轮次与数 and 每条「做不到」都带可证伪形式:
  **A20 可关闭**;else **不许关,先补**。
IMPOSSIBLE(不写 planned):等级是我定的(但判据机械)· **`#664` 的 11 条本身是我从页面挑的** · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
LED=pathlib.Path("RETRACTIONS.md").read_text()

C=[
 dict(id="「性是一块」", tier="T1", round="#653", number="GSS 0.416 · NSFG 0.346",
      falsifier=None, changed=None),
 dict(id="打孩子≠打妻子", tier="T2", round="#641", number="上界 +0.4401",
      falsifier=None, changed=None),
 dict(id="跨团队复现体罚极性", tier="T2", round="#639", number="+0.6301 (n=60)",
      falsifier=None, changed=None),
 dict(id="越稀有越被谴责(婚前性)", tier="T2", round="#665", number="+0.6725 (n=86)",
      falsifier=None, changed="#664 判「待做」-> #665 已推进;#666 再证频率非同义反复(偏相关 +0.3094)"),
 dict(id="换对象 +0.845", tier="T3", round="#667", number="服从四件套 +0.7905 (n=160-162)",
      falsifier="若出现 source∉{barry1977agents} 且按≥2对象分列的儿童体罚变量、真实 n≥60,T3 可升 T2",
      changed="#664 判「做不到」-> #667 推进到 T3(同一第一作者,不是独立团队)"),
 dict(id="换手段 +0.229", tier="T4", round=None, number=None,
      falsifier="若任一跨文化库出现≥3种管教手段各自编码且真实 n≥60,即被推翻", changed=None),
 dict(id="打得多也讲得多 +0.344", tier="T4", round=None, number=None,
      falsifier="同「换手段」", changed=None),
 dict(id="以身作则不站队", tier="T4", round=None, number=None,
      falsifier="若任一调查出现「榜样/模仿学习」题且 n≥300,即被推翻", changed=None),
 dict(id="没有「严厉的社会」+0.125", tier="T4", round=None, number=None,
      falsifier="若任一跨文化库出现≥2道不同性实践的道德判断且共同真实 n≥60,即被推翻", changed=None),
 dict(id="只有纯洁是团 0.335", tier="T4", round=None, number=None,
      falsifier="若出现第二份含道德基础五域或等价结构的数据且 n≥1000,即被推翻", changed=None),
 dict(id="性与同居一样紧 1.10×", tier="T4", round=None, number=None,
      falsifier="若 GSS 或他处出现≥3道同格式同居态度题且 n≥1000,即被推翻", changed=None),
]
print("=== 硬规则①:11 条的当前状态,逐条标出被 #665–#668 改动过的 ===")
print(f"  {'声明':26s}{'等级':>5s}{'轮次':>7s}  数 / 可证伪形式")
for c in C:
    mark=" ⟵ 改动" if c["changed"] else ""
    right=c["number"] if c["number"] else (c["falsifier"] or "")
    print(f"  {c['id']:26s}{c['tier']:>5s}{str(c['round'] or '—'):>7s}  {right[:58]}{mark}")

print("\n=== ⑤ 机械检查:每条「已推进」必须指得出轮次号与数,且轮次号能在账本里找到 ===")
bad=[]
for c in C:
    if c["tier"] in ("T1","T2","T3"):
        ok_r=bool(c["round"]) and (f"## Entry {c['round'][1:]}" in LED)
        ok_n=bool(c["number"])
        print(f"  {c['id']:26s} 轮次 {c['round']} 在账本? {ok_r} · 有数? {ok_n}")
        if not (ok_r and ok_n): bad.append(c["id"])
print(f"  ⇒ 指不出来的 = **{len(bad)}** {bad or ''}")

print("\n=== ③ 每条「做不到」(T4)是否都带可证伪形式 ===")
nof=[c["id"] for c in C if c["tier"]=="T4" and not c["falsifier"]]
print(f"  T4 共 {sum(1 for c in C if c['tier']=='T4')} 条 · **缺可证伪形式的 {len(nof)}** {nof or ''}")

print("\n=== ② 主量 ===")
t12=sum(1 for c in C if c["tier"] in ("T1","T2")); t3=sum(1 for c in C if c["tier"]=="T3")
print(f"  **T1+T2 = {t12}/11 = {t12/11:.3f}**(`#664` 记的是 3/11 = 0.273)")
print(f"  G4 另一条规格(也算 T3):{t12+t3}/11 = {(t12+t3)/11:.3f}")
print(f"  T4(单项目)= {sum(1 for c in C if c['tier']=='T4')}/11")

print("\n=== ④ 正对照:两条已知的「做不到」必须仍是「做不到」 ===")
pos=[c for c in C if c["id"] in ("以身作则不站队",)]
ok_pos=all(c["tier"]=="T4" for c in pos)
print(f"  以身作则不站队 -> {pos[0]['tier']} {'✅' if ok_pos else '⛔ 变了,说明我又漏查'}")
print(f"  年代单位(`#659` 判「做不到」,不在这 11 条内)-> 账本里仍是「做不到」:"
      f"{'✅' if '换不了仪器,只此一具' in LED else '⛔'}")

from lib.gates import Gate
G=Gate("A20 收口:每一条「已推进」都指得出轮次号与数")
p1=G.positive_control("已知的「做不到」必须仍是「做不到」",planted=float(1.0 if ok_pos else 0.0),floor=0.5,spread=0.01)
p2=G.negative_control("机械检查:指不出轮次或数的「已推进」必须为 0",null=float(len(bad)),effect=float(t12+t3),
                      null_spread=0.2,null_kind="声称已推进却指不出证据的条数 —— 它必须是零")
can_close = (len(bad)==0) and (len(nof)==0) and ok_pos
verdict=("**A20 可关闭:5 条已推进全部指得出轮次与数;6 条「做不到」全部带可证伪形式;正对照未变**"
         if (p1 and p2 and can_close) else
         f"**不许关 —— 指不出的 {len(bad)} 条 / 缺可证伪形式的 {len(nof)} 条,先补**")
print(f"\n{verdict}"); print(G)
json.dump(dict(claims=C,t12=t12,t3=t3,t4=sum(1 for c in C if c['tier']=='T4'),
               unsupported=bad,missing_falsifier=nof,can_close=bool(can_close),verdict=verdict,
               unchallenged=True),open(OUT/"a20_closeout.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'a20_closeout.json'}")
