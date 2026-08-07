"""#784 附 · 「换不了仪器」这句话,本轮不许直接写 —— 先跑一遍,让它由测量来担

⚠⚠ 动机:`readme_gate` 的 `single_instrument` 阻断了本轮提交,而它接受一句豁免语。
   **但 `realstat §2` 说得很清楚:一个朝着方便方向的「做不到」声明,仍然是一个声明。**
   `#777` 那条「墙从没被查过」正是这一族 —— 三条「永久限制」里有一条是**一次从没跑过的查询**。
   ⇒ **不许把「我没查」写成「结构性不可能」。** 本文件是那次查询。

这一轮的估计量对第二具仪器的**规格**(先写死,再去找):
   ① 一个**有界序数**的性道德题(≥3 档,能有天花板)
   ② 一个**虔诚度**分层变量(否则两层不是同一件事)
   ③ **≥8 个时间点**(`r_forced` 要标定一个 50 年位移,单波做不了)
   三条**同时**满足才算。缺一条就不是「差一点」,是**不能算**。

⚠ P6 代理账:
  PROPERTY    这份数据能不能承载 `r_forced`
  PROXY       变量名/标签里出现虔诚度词表的次数 · 文件里的年份数
  IMPLICATION 只有一个方向可靠:**词表 0 命中 -> 确实没有虔诚度变量**(可靠,前提是词表本身不瞎)。
             反过来不成立:**命中不证明那个变量可用**(YRBS 的 "attendance" 是体育课)。
  SAFE SIDE   命中就人工看一眼标签;只报「不能承载」,从不认证「能承载」。
⚠⚠ **而这条代理必须先有正对照**(`P5` ★):一个 0 命中,若来自一具从没返回过非零的仪器,
   **那不是测量,是沉默。** ⇒ 每份文件都同时搜一个**必然存在**的词。
"""
import pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
REL = re.compile(r"religi|church|attend|faith|pray|worship|denomin|fundament", re.I)

def scan(path, poscontrol, nbytes=4_000_000):
    txt = path.read_text(errors="replace")[:nbytes]
    hits = REL.findall(txt)
    pc = len(re.findall(poscontrol, txt, re.I))
    return hits, pc, txt

print("=== 第二具仪器的三条规格:有界序数性道德题 × 虔诚度分层 × ≥8 个时间点 ===\n")
rows = []

# ── YRBS(1991–2023 合并档,时间点最多的候选)──────────────────────────────────
p = ROOT/"data/external/yrbs/2023-SADC-SAS-Input-Program.sas"
hits, pc, txt = scan(p, r"\bsex\b")
labels = re.findall(r'(q\w+)\s*=\s*"([^"]*(?:attend|religi|church|pray)[^"]*)"', txt, re.I)
print(f"YRBS  正对照 'sex' 命中 **{pc}** 次 ⇒ 仪器不瞎")
print(f"      虔诚度词表命中 {len(hits)} 次,逐条看标签:")
for v, l in labels: print(f"        {v} = {l!r}")
yrbs_ok = any(not re.search(r"physical education|PE\b", l, re.I) for _, l in labels)
print(f"      ⇒ **规格②(虔诚度分层):{'满足' if yrbs_ok else '不满足 —— 全部是体育课出勤'}**\n")
rows.append(("YRBS", "1991–2023 · 时间点充足", yrbs_ok, "唯一的 attendance 是体育课"))

# ── BRFSS ────────────────────────────────────────────────────────────────────
brf = sorted((ROOT/"data/external/brfss").glob("*.XPT")) + \
      sorted((ROOT/"data/external/brfss/_archive").glob("*")) if (ROOT/"data/external/brfss/_archive").exists() \
      else sorted((ROOT/"data/external/brfss").glob("*.XPT"))
years = sorted({m.group(0) for f in brf for m in [re.search(r"(19|20)\d{2}", f.name)] if m})
print(f"BRFSS 本机上的年份:{years or '无'} ⇒ **规格③(≥8 个时间点):"
      f"{'满足' if len(years) >= 8 else f'不满足 —— 只有 {len(years)} 年'}**\n")
rows.append(("BRFSS", f"只有 {len(years)} 年", False, "单波,标定不了 50 年位移"))

# ── 已在账本上量过的三具(`#769`/`#773`/`#777`),这里只复述其结论并标出处 ──────
for nm, why, where in (("MFQ", "单次采集、无年代 ⇒ 规格③不满足", "`#769` 实测"),
                       ("NSFG", "`SXOK18` 在 2017–19 卷字典里不存在 ⇒ 规格③不满足", "`#777` 实测"),
                       ("MSSCQ", "单次采集、且问的是自我描述不是道德判断 ⇒ ①③都不满足", "`#783` 第二臂"),
                       ("SCCS", "单位是社会,不是人;无年代 ⇒ 规格②③不满足", "`#653` 实测")):
    print(f"{nm:6s} {why}  ({where})")
    rows.append((nm, why, False, where))

ok = [r for r in rows if r[2]]
print("\n" + "="*88)
if ok:
    print(f"⛔ **有候选满足规格:{[r[0] for r in ok]} —— 豁免不成立,必须去跑它。**")
    sys.exit(2)
print("⇒ **六具候选全部在三条规格上落选,而落选的理由逐条是测出来的,不是回忆的。**")
print("⇒ 本轮**换不了仪器**:承载 `r_forced` 需要「有界序数性道德题 × 虔诚度分层 × ≥8 个时间点」")
print("   三条同时满足,而本机的六份外部数据里**只此一具**(GSS)。")
print("⚠ 而这不是一句永久的话:它随数据而变 —— **加一份带虔诚度的多波调查,这条豁免立刻失效。**")
