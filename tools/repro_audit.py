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


# ---------------------------------------------------------------------------
# #158:反方向 —— **账本判为错的那个旧值,还在不在新输出里?**
#
# `#157` 的两个实例(`#117e` 崩溃 · `#90c` 静默)都是同一个病:
# 诊断写进账本、校正记进账本,**代码从没改**。而 v1 只查「账本的数在不在新输出里」,
# 于是它对 `#90` 的 ±30.8 完全无感 —— **那个数在新输出里,而它本该消失**。
#
# ⚠ P6 代理账:
#   PROPERTY   这一轮的代码,是否仍在产生一个账本已判为错的旧值
#   PROXY      账本条目里"校正类"句子中的**第一个数字**(旧值 A),是否出现在新 stdout 里
#   IMPLICATION 只有一个方向可靠:**A 出现在新输出里 -> 必须人工读**(可靠)。
#              反过来"A 不在 -> 代码已修"**不可靠** —— A 可能只是这一轮没打印。
#   SAFE SIDE  输出是必读清单,不是判决。
CORR=('WRONG','撤回','校正后','corrected','Corrected','falls from','改善到','降到',
      'overstat','withdrawn','WITHDRAWN','应作','而不是','实际是')
# ⚠ 第一版写成 `(\d+\.\d+)\s*(?:->|→|到|to)\s*(\d+\.\d+)` —— markdown 的 `**`、`±`、`%`
#   把匹配打断了,于是它连 `#90c` 自己那句 "falls from ±30.8 to **±23.7 pp**" 都抓不到。
#   允许两个数字之间与前后有少量非数字字符。
ARROW=re.compile(r'(\d+\.\d{1,4})[^0-9]{0,12}?(?:->|→|到|to)[^0-9]{0,12}?(\d+\.\d{1,4})')

def stale_pairs(entry:int,ledger='RETRACTIONS.md'):
    """返回该条目里 (旧值, 新值) 对 —— 只在含校正类词汇的行上抽。"""
    L=(ROOT/ledger).read_text().splitlines()
    a=None; body=[]
    for i,l in enumerate(L):
        m=re.match(r'^## Entry (\d+),',l)
        if m:
            if a is not None: break
            if int(m.group(1))==entry: a=i; continue
        if a is not None: body.append(l)
    # ⚠ #158:并非每个「A to B」都是校正。`Entry 69` 的 "the control at **g=0.15** moves ρ to
    #   −0.314" 是**参数 → 结果**,被 v1 当成了旧值 → 新值。排除第一个数字前带赋值/参数标记的。
    PARAM=re.compile(r'(?:g\s*=|scale|rank|at\s+g|=\s*|K\s*=|n\s*=|seed)\s*$')
    out=[]
    for l in body:
        if not any(w in l for w in CORR): continue
        for m in ARROW.finditer(l):
            old,new=m.group(1),m.group(2)
            if old==new: continue
            if PARAM.search(l[:m.start(1)][-16:]): continue      # 参数 -> 结果,不是校正
            out.append((old,new,l.strip()[:120]))
    return out


def stale_still_live(entry:int,round_paths,timeout=1500):
    """#158 的正确判据:旧值 A **只有在**「A 出现了 **且** 新值 B 在该条目的**任何一轮**里
    都不出现」时,才算"账本判为错的数还活着"。

    ⚠ 三种假阳性,都是实测踩出来的:
      ① 参数 → 结果(`Entry 69`:"at g=0.15 moves ρ to −0.314")—— 已在 `stale_pairs` 里排除
      ② 一条条目覆盖**多轮**,旧值在前一轮、新值在后一轮(`Entry 126`:`R20` 打 0.1822、
         `R21` 打 0.1449)—— 由本函数的"任何一轮"排除
      ③ 一轮**本来就该**同时报两个臂 —— 本函数不区分,留在必读清单里
    """
    outs=[]
    for p in round_paths:
        rc,o=run_round(p,timeout); outs.append(o)
    blob='\n'.join(outs)
    res=[]
    for old,new,ctx in stale_pairs(entry):
        a_live=old in blob; b_live=new in blob
        res.append(dict(old=old,new=new,old_in_output=a_live,new_in_output=b_live,
                        verdict=('⚠ 旧值还活着' if (a_live and not b_live)
                                 else ('✅ 新值已在输出里' if b_live else '— 两者都不在')),
                        ctx=ctx))
    return res
