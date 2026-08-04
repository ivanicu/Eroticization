#!/usr/bin/env python3
"""
tools/repro_audit.py -- 账本条目声称的数字,那一轮的代码现在还跑得出来吗?

第四条边。前三条:`readme_ledger_audit`(README→账本)· `round_status`(轮次→账本)·
`guard_lint`(轮次→守卫)。**这一条是 账本→轮次。**

`#154a` 正是靠"跑一遍 `A11/R14` 本身"才分清了「我复现错了」和「README 错了」——
而那一步完全是手工的。`#141` 也在同一步上被绊住:它重跑了一个已作废的轮次,
却没有把输出与账本比对。

⚠ P6 代理账:
  PROPERTY   这条账本条目的数字,现在还是这一轮的输出吗
  PROXY      账本条目里出现的"有辨识度的数字",是否出现在该轮 stdout 里
  IMPLICATION 只有一个方向可靠:**账本里的数字在新 stdout 里找不到 -> 必须人工读**。
             反过来"找到了 -> 一致"**不可靠** —— 同一个数可能出自不同的量。
  SAFE SIDE  输出是**必读清单**,不是判决。轮次内部随机性、下采样、种子都会让数字
             轻微漂移,所以匹配用**容差**,而漂移本身也报出来。
"""
import re,sys,subprocess,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[1]
PY=str(ROOT/'.venv/bin/python')
# ⚠ #157:第一版要求 **>=2 位小数**,于是它看不见 `22.6` / `16.3` / `23.7` ——
#   而那正是 README 头条(「±22.6 / ±16.3 / ±23.7 pp」)用的格式。
#   **一个检查看不见它被造出来要检查的那种格式。** 放宽到 >=1 位。
NUM=re.compile(r'(?<![\w.])[-+]?\d+\.\d{1,4}(?![\w])')

def ledger_numbers(entry:int,ledger='RETRACTIONS.md'):
    L=(ROOT/ledger).read_text().splitlines()
    a=None
    for i,l in enumerate(L):
        m=re.match(r'^## Entry (\d+),',l)
        if m:
            if a is not None: return set(NUM.findall('\n'.join(L[a:i])))
            if int(m.group(1))==entry: a=i
    return set(NUM.findall('\n'.join(L[a:]))) if a is not None else set()

def run_round(path,timeout=1800):
    p=ROOT/path
    try:
        r=subprocess.run([PY,str(p)],cwd=str(ROOT),capture_output=True,text=True,timeout=timeout)
        return r.returncode,r.stdout+r.stderr
    except subprocess.TimeoutExpired:
        return -9,'TIMEOUT'

def compare(entry,path,tol=0.02,timeout=1800):
    want=ledger_numbers(entry)
    rc,out=run_round(path,timeout)
    got=set(NUM.findall(out))
    gotf=[float(x) for x in got]
    miss=[]
    for w in sorted(want):
        wf=float(w)
        if any(abs(wf-g)<=tol*max(abs(wf),1e-9) for g in gotf): continue
        miss.append(w)
    return dict(entry=entry,path=path,rc=rc,n_want=len(want),n_got=len(got),
                missing=miss,timeout=(out=='TIMEOUT'))
