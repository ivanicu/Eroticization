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
# A02/R034 的 `mx < 0.4 -> "still three axes"` 是一个**选定**的阈值,而 #141 证明
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


# ============================================================================
# #168:哪些轮次的头条数字**没有误差棒**
# ----------------------------------------------------------------------------
# #167b 给出了可机械查的形状:如果一个数的"不确定度"只在零臂里出现,
# 而真实臂本身没有任何抖动来源,那它没有误差棒。
#
# P6 代理账:
#   PROPERTY    这一轮的头条数字有一个精度估计
#   PROXY       真实臂在持久化结果里出现 >=2 次且**取值不同**
#   IMPLICATION 真实臂只出现一次,或多次但逐行相同 => 该 artifact 里没有实现展布(**缺失方向可读**)
#   WITNESS     R147/grid.csv 带着 `r_half_sd` 列(有展布!)而它发表的是 `rel_resid`,
#               后者在 3 个 seed 上逐字节相同 —— **一个轮次可以带展布列却仍然没有误差棒**
#   SAFE SIDE   只在**缺失**方向判。带展布列 != 头条数字有误差棒 -> UNVERIFIED,永不判"有"。
# ============================================================================
import re as _re
_ARM  = _re.compile(r'real|真实|actual|observed|obs', _re.I)
_SEED = _re.compile(r'^(seed|sd_?|rep|draw|split|iter)\w*$', _re.I)
_SPRD = _re.compile(r'(_sd|_se|sd_|se_|std|spread|ci_|_err|_ci)', _re.I)

def error_bar_scan(root='.'):
    import pandas as pd, pathlib
    out=[]
    for csv in sorted(pathlib.Path(root).glob('E01_*/A*/R*/results/*.csv')):
        try: df=pd.read_csv(csv)
        except Exception: continue
        if df.empty: continue
        rnd=csv.parents[1].name
        armcol=next((c for c in df.columns if df[c].dtype==object
                     and df[c].astype(str).str.match(_ARM).any()), None)
        seedcol=next((c for c in df.columns if _SEED.match(str(c))), None)
        has_sp=[c for c in df.columns if _SPRD.search(str(c))]
        num=[c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])
             and c!=seedcol and not _SPRD.search(str(c))]
        sub=df[df[armcol].astype(str).str.match(_ARM)] if armcol else df
        if len(sub)==0: continue
        if len(sub)==1:
            verdict,why='NO_SPREAD','真实臂在 artifact 里只有一行 —— 没有实现展布'
        elif seedcol and sub[seedcol].nunique()>1 and num and \
             all(sub.groupby(seedcol)[c].first().nunique()==1 for c in num):
            verdict,why='SEED_IS_FAKE',(f'真实臂在 {sub[seedcol].nunique()} 个 {seedcol} 上'
                                        f'**逐行相同** —— 种子不驱动真实臂(#167a 的形状)')
        else:
            verdict,why='UNVERIFIED','有多行且取值不同 —— 但"有展布"不等于"头条数字有误差棒"'
        out.append(dict(round=rnd,artifact=csv.name,rows=len(sub),arm_col=armcol or '-',
                        seed_col=seedcol or '-',spread_cols=','.join(has_sp) or '-',
                        verdict=verdict,why=why))
    return pd.DataFrame(out)

def report_error_bars(root='.'):
    import pandas as pd
    T=error_bar_scan(root)
    if T.empty: print('没有可扫描的 results/*.csv'); return T
    print(f"\n扫描 {len(T)} 个 artifact,覆盖 {T['round'].nunique()} 个轮次")
    print(T.verdict.value_counts().to_string())
    for v in ['SEED_IS_FAKE','NO_SPREAD']:
        S=T[T.verdict==v]
        if len(S):
            print(f"\n── {v} ({len(S)}) ──")
            for _,r in S.iterrows():
                print(f"  {r['round'][:52]:<54}{r.artifact:<28}{r.why}")
    print("\n⚠ SAFE SIDE(#P6):只在**缺失**方向可读。UNVERIFIED **不等于**有误差棒 ——"
          "\n   R147 带着 r_half_sd 列却仍然没有 rel_resid 的误差棒(#167a)。")
    return T


# ============================================================================
# #178:打印出来的"事实"里,有多少是手写的
# ----------------------------------------------------------------------------
# `#177d`:`hardcoded_thresholds` 只扫**判定阈值**(驱动 verdict 字符串的字面量),
# 不扫 `print` 里那些**看起来像测量结果**的数字。
# 而 R222 的报告里 "只认反引号 -> 0" 就是硬编码的:我一补出处,那个 0 当场变成假话,
# **而没有任何检查会抓到它。一个印在报告里的常量,是一条不会被任何检查抓到的声明。**
#
# P6 代理账:
#   PROPERTY    这一行打印的数字是从数据算出来的
#   PROXY       它出现在 f-string 的**静态文本**里(不在 `{}` 表达式内),且形如一个量值
#               (带符号小数 ≥2 位 · 百分比 · 比值)
#   IMPLICATION 命中 => 这个数**没有**经过本次运行的计算(**命中方向可读**)
#   WITNESS     一个静态数字可以是**引用**(「注册阈值 20%」「#118a 的 +0.0339」)而非结果 ——
#               所以命中项必须人工分诊「它是在陈述本次结果,还是在引用别处」
#   SAFE SIDE   只在**命中**方向判。没命中 != 打印的都是算出来的:
#               一个从别处 import 的常量、一个上一轮抄下来的变量,这条规则都看不见。
# ============================================================================
import ast as _ast
_PRINTNUM = _re.compile(r'[+−-]\s?\d\.\d{2,4}|\b\d{1,3}(?:\.\d+)?\s?%|\b\d+(?:\.\d+)?\s?[×倍]')

def printed_literals(paths=None, root='.'):
    """扫 run.py 的 print(...),找 f-string **静态文本**里的量值字面量。"""
    import pathlib, pandas as pd
    if paths is None:
        paths=sorted(pathlib.Path(root).glob('E01_*/A*/R*/run.py'))
    out=[]
    for p in paths:
        try: tree=_ast.parse(pathlib.Path(p).read_text())
        except SyntaxError: 
            out.append(dict(path=str(p),line=0,literal='<SyntaxError>',text='')); continue
        for node in _ast.walk(tree):
            if not (isinstance(node,_ast.Call) and getattr(node.func,'id',None)=='print'): continue
            for a in node.args:
                pieces=[]
                if isinstance(a,_ast.JoinedStr):
                    pieces=[v.value for v in a.values if isinstance(v,_ast.Constant)
                            and isinstance(v.value,str)]
                elif isinstance(a,_ast.Constant) and isinstance(a.value,str):
                    pieces=[a.value]
                for s in pieces:
                    for m in _PRINTNUM.findall(s):
                        out.append(dict(path=str(p),line=getattr(node,'lineno',0),
                                        literal=m.strip(),text=s.strip()[:100]))
    return pd.DataFrame(out)
