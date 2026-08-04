#!/usr/bin/env python3
"""
tools/guard_lint.py -- 哪些轮次缺了后来为它们那类失败而写的守卫。

#127 的 NEXT。六个守卫每一个都写于一次具体的失败之后,但**没有任何东西检查一个轮次
有没有用它该用的守卫**。这个 linter 做那件事。

⚠ P6 的代理账(必须先写,否则这个工具本身就是它要防的那种检查):

  PROPERTY   这一轮是否受到了它需要的那类保护
  PROXY      源码里有没有出现对应的调用,以及有没有出现触发它的模式
  IMPLICATION 只有一个方向可靠:**模式在而调用不在 -> 确实没保护**(可靠)。
              反过来"调用在 -> 受到了保护"**不可靠** —— 调用可能传了宽松的 tol,
              或者用在了错的量上(#119d、#125 都是这种)。
  WITNESS    A11R20 有 check_coverage 与 check_columns 却仍漏掉共享 item(#126c),
              所以"全绿"从来不等于"干净"。
  SAFE SIDE  只在**缺失**方向下结论;绿色一律报成 "未标记",不报成 "已保护"。
"""
import re,sys,pathlib,collections
ROOT=pathlib.Path(__file__).resolve().parents[1]
RULES=[
 # (名字, 触发模式, 应有的调用, 起因)
 ('columns',   r'\.groupby\(|\.agg\(|pivot_table',       r'check_columns',   '5 次 pandas 访问器撞名'),
 ('coverage',  r'\n\s*(if .*: *)?break\b|for .* in (IDENT|pairs|FC|SEL|CATS)\b', r'check_coverage', '6 次静默 cap 改变结论'),
 ('disjoint',  r'get_dummies|one[_-]?hot|np\.delete\(P|Xp=|H=pd',              r'check_disjoint_items', '#126c 设计时漏掉的共享 item'),
 ('degenerate',r'plant|种植|g=0\.0|\[0\.0,',                                    r'degenerate_matches_reference', '#124f 退化臂没复用种子'),
 ('resolvable',r'negative_control|artifact_cannot_explain',                     r'require_resolvable_first', '#120d 对未分辨的量问了形状'),
 ('nullspread',r'negative_control\(',                                           r'null_spread',      '#125 小效应下零按自身展布'),
]
rounds=sorted(ROOT.glob('E01*/A*/R*/run.py'))
rows=[]
for p in rounds:
    t=p.read_text(errors='ignore')
    miss=[n for n,trig,call,_ in RULES if re.search(trig,t) and not re.search(call,t)]
    rows.append((p.relative_to(ROOT), len(miss), miss, 'from lib.gates' in t))
n=len(rows)
print(f"{n} 个轮次脚本\n")
print(f"  用了 lib.gates 的:{sum(r[3] for r in rows)}   完全没用的:{n-sum(r[3] for r in rows)}")
cnt=collections.Counter(m for r in rows for m in r[2])
print(f"\n  各守卫的缺失数(仅在**触发了对应模式**的轮次里计):")
for name,trig,call,why in RULES:
    trg=sum(1 for p in rounds if re.search(trig,p.read_text(errors='ignore')))
    print(f"    {name:11s} 触发 {trg:3d} 轮  缺失 {cnt[name]:3d} 轮 ({100*cnt[name]/max(trg,1):3.0f}%)   起因:{why}")
print(f"\n  按缺失数排序的轮次(前 12):")
for rel,k,miss,_ in sorted(rows,key=lambda r:-r[1])[:12]:
    print(f"    {k}  {str(rel)[:72]:<72} {','.join(miss)}")
print(f"\n  缺 0 个的轮次:{sum(1 for r in rows if r[1]==0)}/{n}")
print("\n⚠ SAFE SIDE(#P6):以上只在**缺失**方向可读。绿色 = 未被标记,不等于已受保护 ——")
print("   A11R20 有 check_coverage 与 check_columns,仍漏掉了 #126c 的共享 item。")

def standing(paths):
    """只在**现存声明**的轮次上读:那 19 轮之外的缺失是"守卫还不存在",没有判别力。"""
    import re
    print("\n\n=== 只看现存声明背后的轮次(其余的缺失只说明守卫尚未存在) ===")
    for rel in paths:
        p=ROOT/rel
        if not p.exists(): print(f"  ?? {rel}"); continue
        t=p.read_text(errors='ignore')
        miss=[n for n,trig,call,_ in RULES if re.search(trig,t) and not re.search(call,t)]
        print(f"  {'OK ' if not miss else 'FLAG'} {rel.split('/')[1][:34]:<34} {rel.split('/')[2][:30]:<30} {','.join(miss) or '-'}")


# ---------------------------------------------------------------------------
# #142:写死的判定常数扫描(由 #141c 触发)
#
# A02/R10 的 `mx < 0.4 -> "still three axes"` 是一个**选定**的阈值,而 #141 证明
# 它的自助区间跨过它 30.2% 的时间 —— 也就是说,一个字面常数在源码里决定了一条
# 现存声明的判定,而没有人问过那个常数是量出来的还是选出来的。
#
# ⚠ P6 代理账:
#   PROPERTY   这个判定的阈值有没有对过一个**测量出来的**参照
#   PROXY      源码里"与字面浮点常数比较"的行,且该比较驱动一个结论字符串或 PASS/FAIL
#   IMPLICATION 只有一个方向可靠:**命中 -> 确实有一个字面常数在做判定**(可靠)。
#              反过来"没命中 -> 阈值是量出来的"**不可靠** —— 阈值可能来自一个变量,
#              而那个变量本身也可能是选的。
#   SAFE SIDE  只在**命中**方向下结论;未命中报"未标记",不报"已对过地板"。
#
# 白名单:与自身展布比较(`2*spread`、`boot`、`sd`、`null`)不是选定阈值,那正是正确做法。
VERDICT_WORDS = ('PASS','FAIL','still','collapse','UNVERIFIED','CONFIRMED','->','判定','结论')
MEASURED_HINTS = ('boot','spread','sd','std','null','thr','floor','perc','quantile','展布','地板','零')

def hardcoded_thresholds(paths=None):
    import re
    rounds = sorted(ROOT.glob('E01*/A*/R*/run.py')) if paths is None else [ROOT/p for p in paths]
    hits=[]
    for p in rounds:
        for i,line in enumerate(p.read_text(errors='ignore').splitlines(),1):
            s=line.strip()
            if s.startswith('#') or not re.search(r'[<>]=?\s*-?\d*\.\d+|[<>]=?\s*-?\d+\b', s): continue
            if not any(w in s for w in VERDICT_WORDS): continue
            if any(h in s.lower() for h in MEASURED_HINTS): continue      # 与自身展布比 = 正确做法
            hits.append((str(p.relative_to(ROOT)), i, s[:104]))
    return hits

def report_thresholds(paths=None, title=''):
    hits=hardcoded_thresholds(paths)
    print(f"\n=== 写死的判定常数{title}:{len(hits)} 处 ===")
    for rel,i,s in hits:
        parts=rel.split('/')
        print(f"  {parts[1][:26]:<26} {parts[2][:26]:<26} :{i:<4} {s}")
    print("\n⚠ SAFE SIDE(#P6):只在**命中**方向可读。未命中 = 未被标记,**不等于**阈值是量出来的。")
    return hits
