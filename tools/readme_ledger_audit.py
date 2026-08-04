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
