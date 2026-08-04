#!/usr/bin/env python3
"""
tools/readme_ledger_audit.py -- README 上的数字,有没有被后来的账本条目撤回过?

由 #143 触发。#141 花了一整轮"重新定价"一条声明,而 **Entry 24 在 118 条之前就已经把它
杀掉了** —— 我用**重跑那一轮**的方式去审计它,却没读那条声明后来的账本条目。
一个轮次自己的输出**不是**那条声明的当前状态;账本才是。

⚠ P6 代理账:
  PROPERTY   README 上的这个数,是不是已经被后来的某条账本条目撤回或改写
  PROXY      该数字的字面串,出现在账本里一行**含撤回类词汇**的文本中
  IMPLICATION 只有一个方向可靠:**命中 -> 这个数确实出现在一条撤回语境里,必须人工读**。
             反过来"未命中 -> 这个数是当前的"**不可靠** —— 撤回可能换了写法或没写数字。
  WITNESS    大量命中是良性的:那个数字是撤回行里的**更正值**,不是被撤的值。必须人工分诊。
  SAFE SIDE  输出是**待读清单**,不是判决。
"""
import re,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
KILL=('撤回','withdrawn','WITHDRAWN','retracted','RETRACTED','killed','RETRACTION',
      'no longer follows','降级','INVALID','reduced')

def audit(readme='README.md',ledger='RETRACTIONS.md'):
    R=(ROOT/readme).read_text(); L=(ROOT/ledger).read_text().splitlines()
    nums=set(re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',R))
    nums|=set(re.findall(r'(?<![\w.])\d\.\d{2}\s*(?:of|/)\s*3',R))
    out=[]
    for n in sorted(nums):
        for i,l in enumerate(L,1):
            if n in l and any(k in l for k in KILL):
                out.append((n,i,l.strip())); break
    return len(nums),out

if __name__=='__main__':
    tot,hits=audit()
    print(f"README 里可辨识的数字 {tot} 个;{len(hits)} 个出现在账本的撤回语境里(**待读清单,非判决**)\n")
    for n,i,l in hits: print(f"  {n:<10} RETRACTIONS.md:{i}\n     {l[:170]}\n")
    print("⚠ SAFE SIDE(#P6):大量命中是良性的 —— 那个数字往往是撤回行里的**更正值**。")
    print("   未命中**不等于**这个数是当前的。输出必须人工分诊。")


# ---------------------------------------------------------------------------
# #144:README 内部一致性 —— 同一个事实被写了两遍,而修正只落在一遍上
#
# #143 找到的是"账本的撤回没走到前页"。更坏的一种是:**走到了前页的一半**。
# 本 README 有两套并行叙述(中文表 + 英文正文),同一条声明在两处各写一次:
#   色情诱导 82.7%    行 56(已修)  vs  行 131(仍写已被认定为"账本里根本没有"的 +0.2515)
#   覆盖度定律 0.815   行 87        vs  行 122
#   三条轴 2.95       行 54        vs  行 111
# 这是 §P16「一个事实一个家」的直接违反,而它的代价与 #143 相同:
# 读者拿到哪一份取决于他先读到哪一段。
#
# ⚠ P6 代理账:
#   PROPERTY   同一条声明的两处叙述,数字是否一致
#   PROXY      同一个**引用标记**(#NN / ANN / EntryNN)出现在两侧,而两侧的数字集合不同
#   IMPLICATION 只有一个方向可靠:**数字集合不同 -> 两处确实写了不同的数,必须人工读**。
#              反过来"相同 -> 一致"**不可靠** —— 两处可能都是旧的。
#   SAFE SIDE  输出是待读清单,不是判决。
import re as _re

def internal_consistency(readme='README.md'):
    lines=(ROOT/readme).read_text().splitlines()
    def nums(s): return set(_re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',s))
    def has_cjk(s): return bool(_re.search(r'[\u4e00-\u9fff]',s))
    cites={}
    for i,l in enumerate(lines,1):
        for c in set(_re.findall(r'#\d{1,3}\b|\bA\d{2}\b|Entry \d{1,3}',l)):
            cites.setdefault(c,[]).append((i,has_cjk(l),nums(l)))
    out=[]
    for c,rows in sorted(cites.items()):
        sides={True:set(),False:set()}
        for i,cj,ns in rows: sides[cj]|=ns
        if rows and len(set(cj for _,cj,_ in rows))==2 and sides[True]!=sides[False]:
            out.append((c,[(i,cj,sorted(ns)) for i,cj,ns in rows]))
    return out

def uncited_numbers(readme='README.md'):
    """#144b:带数字但**不带任何引用标记**的行 —— 它们对任何基于引用的审计都是不可见的。

    #143 的那个真错误(行 131 的 rho +0.2515 与 "85% surviving response-style control")
    正是这样一行:账本明说「账本里根本没有这个数」,而它活了下来,**因为没有任何东西
    把它连回账本**。一个不带出处的数字,是一个无法被撤回的数字。"""
    # ⚠ #145:第一版**按行**判,而引用往往在同一**段**的别处 —— 它把盲区高估了 14 倍
    #   (17 行 vs 实际 2 段,其中一段是 python 版本号)。检查的单位必须与"出处"的单位一致。
    txt=(ROOT/readme).read_text(); paras=_re.split(r'\n\s*\n',txt)
    off=1; out=[]
    for para in paras:
        n=para.count('\n')+2
        ns=_re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',para)
        if ns and not _re.search(r'#\d{1,3}\b|RETRACTIONS|Entry \d|\bA\d{2}\b|\bR\d{2}\b',para):
            out.append((off,ns,para.strip().replace('\n',' ')[:110]))
        off+=n
    return out

if __name__=='__main__' and '--internal' in __import__('sys').argv:
    hits=internal_consistency()
    print(f"\n=== README 内部一致性:{len(hits)} 个引用标记在两套叙述里带着不同的数字 ===")
    for c,rows in hits:
        print(f"  {c}")
        for i,cj,ns in rows:
            print(f"     行{i:>4} {'中文表' if cj else '英文正文'}  {ns}")
    print("\n⚠ SAFE SIDE(#P6):数字集合相同**不等于**一致 —— 两处可能都是旧的。待读清单。")
    un=uncited_numbers()
    print(f"\n=== 带数字但**不带出处**的行:{len(un)} 行(对任何基于引用的审计都不可见)===")
    for i,ns,s in un: print(f"  行{i:>4}  {ns}\n        {s}")
    print("\n⚠ 一个不带出处的数字,是一个**无法被撤回**的数字。#143 的那个真错误就在这类行里。")
