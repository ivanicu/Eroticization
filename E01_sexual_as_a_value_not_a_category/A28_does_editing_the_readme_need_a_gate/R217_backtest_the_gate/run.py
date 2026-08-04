import os,sys,pathlib,subprocess,shutil
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A28 R217 -- 那道闸,用 git 历史回测它会拦下多少

`#171c`:一次修复必须**自带它自己的验收**。`tools/readme_gate.py` 把四条规则接成单一入口。
现在的问题是它值不值一道强制闸。

ESTIMAND        最近 N 个改动过 README 的提交里,**这次改动让任一规则的命中数上升**的比例。
IDENTIFICATION  判的是**增量**不是绝对值 —— 当前状态本来就有 2 个已确认的假阳性
                (`+0.023` 同串两指 `#170b`;`3.14` 是 python 版本号)。
                绝对值判会把每一个提交都拦下,那样拦截率 100% 而信息量 0。
KILL            **拦截率 > 20% -> "改 README" 一直缺一道闸,而不是偶尔出错。**
POSITIVE CTRL   `7b1951f`(`#170a` 的去重提交)**必须被拦下** —— `#171` 已确证它删掉了 `+0.2523`。
                拦不下它,整个回测不可读。
NEGATIVE CTRL   一个只改 RETRACTIONS.md、不碰 README 的提交必须**不被拦** —— 见下。
IMPOSSIBLE      规则只看数字。一句被删掉的**限定语**不带数字,这道闸看不见。
"""
import pandas as pd, hashlib
import readme_ledger_audit as A
import readme_gate as G
OUT=pathlib.Path(__file__).parent/'results'
WORK=OUT/'_bt'; WORK.mkdir(exist_ok=True)

def hits_at(rev, ref):
    """把 rev 的三份文件落到仓库内的临时目录,在那里跑四条规则。ref 供 numbers_that_left。"""
    # ⚠ 第一版:任一文件缺失就 return None -> 老提交里还没有 README_zh.md,
    #   于是 20 个样本静默变成 9 个,而脚本一声没吭。#118c 又一次。
    #   缺文件按**空文件**处理(那正是当时的事实:那份 README 还不存在),并记下来。
    missing=[]
    for f in ['README.md','README_zh.md','RETRACTIONS.md']:
        try: t=subprocess.run(['git','show',f'{rev}:{f}'],capture_output=True,text=True,check=True).stdout
        except subprocess.CalledProcessError:
            if f=='README.md': return None,['README.md']    # 这个缺了才是真的判不了
            t=''; missing.append(f)
        (WORK/f).write_text(t)
    cwd=os.getcwd(); os.chdir(WORK)
    try: _,h=G.run_gate(rev=ref,quiet=True)
    except Exception: h=None
    finally: os.chdir(cwd)
    return h,missing

dropped=[]
commits=subprocess.run(['git','log','--format=%h','-60','--','README.md'],
                       capture_output=True,text=True,check=True).stdout.split()
N_WANT=30
print(f"改动过 README 的提交(最近){len(commits)} 个,取前 {N_WANT} 个判")
rows=[]
for c in commits[:N_WANT]:
    par=subprocess.run(['git','rev-parse','--short',f'{c}^'],capture_output=True,text=True).stdout.strip()
    gpa=subprocess.run(['git','rev-parse','--short',f'{c}^^'],capture_output=True,text=True).stdout.strip()
    if not par: continue
    hc,mc=hits_at(c,par); hp,mp=hits_at(par,gpa or par)
    if hc is None or hp is None: dropped.append(c); continue
    up={k:hc.get(k,0)-hp.get(k,0) for k in hc}
    subj=subprocess.run(['git','log','-1','--format=%s',c],capture_output=True,text=True).stdout.strip()
    rows.append(dict(commit=c,blocked=any(v>0 for v in up.values()),
                     **{f'd_{k}':v for k,v in up.items()},subject=subj[:74]))
T=pd.DataFrame(rows); T.to_csv(OUT/'backtest.csv',index=False)
from lib.gates import check_coverage
print(f"覆盖:判了 {len(T)}/{min(N_WANT,len(commits))},丢弃 {len(dropped)} 个:{dropped}")
check_coverage(len(T),min(N_WANT,len(commits)),'R217 回测样本',tol=0.15)
print(f"\n{'提交':<10}{'拦':<5}{'Δ点名':>7}{'Δ离开':>7}{'Δ无出处':>8}{'Δ不一致':>8}  标题")
for _,r in T.iterrows():
    print(f"{r.commit:<10}{'YES' if r.blocked else '-':<5}{r.d_named_defects:>7}"
          f"{r.d_numbers_that_left:>7}{r.d_uncited_numbers:>8}{r.d_internal_consistency:>8}  {r.subject[:56]}")
rate=T.blocked.mean()
print(f"\n**拦截率 {T.blocked.sum()}/{len(T)} = {100*rate:.1f}%**   (注册阈值 20%)")

from lib.gates import Gate
g=Gate('改 README 要不要一道强制闸')
pc=T[T.commit.str.startswith('7b1951f')]
g.asserted('正对照:#170a 的去重提交 7b1951f 必须被拦下',
           bool(len(pc) and pc.blocked.all()),
           f"{'在样本里,blocked=' + str(bool(pc.blocked.all())) if len(pc) else '不在最近 20 个样本里'}")
g.asserted('可判前提:回测判的是增量,不是绝对值',True,
           '当前绝对命中 2 个,全为已确认假阳性(+0.023 同串两指;3.14 是 python 版本号)')
g.asserted('注册的 kill:拦截率 > 20%',rate>0.20,f"{100*rate:.1f}%")
g.threshold_outside_noise('拦截率离 20% 有多远',float(rate),0.20,
                          float((rate*(1-rate)/max(len(T),1))**0.5))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 被拦下的那几个,是真缺陷还是假阳性 -------------------------------------
# 一个拦截率没有精度就没有意义:5/30 若全是假阳性,这道闸只是噪声。
print("\n---- 被拦提交的逐条分诊 ----")
tri=[]
for _,r in T[T.blocked].iterrows():
    c=r.commit; par=subprocess.run(['git','rev-parse','--short',f'{c}^'],capture_output=True,text=True).stdout.strip()
    hc,_=hits_at(c,par)
    for f in ['README.md','README_zh.md','RETRACTIONS.md']:
        try: t=subprocess.run(['git','show',f'{c}:{f}'],capture_output=True,text=True,check=True).stdout
        except subprocess.CalledProcessError: t=''
        (WORK/f).write_text(t)
    cwd=os.getcwd(); os.chdir(WORK)
    try: left=A.numbers_that_left(rev=par)
    finally: os.chdir(cwd)
    toks=list(left.token) if len(left) else []
    tri.append(dict(commit=c,d_named=r.d_named_defects,d_left=r.d_numbers_that_left,
                    left_tokens=', '.join(toks[:6]),subject=r.subject[:60]))
    print(f"  {c}  Δ点名 {r.d_named_defects:+d}  Δ离开 {r.d_numbers_that_left:+d}"
          f"  离开的数:{', '.join(toks[:6]) if toks else '—'}")
    print(f"      {r.subject[:70]}")
TR=pd.DataFrame(tri); TR.to_csv(OUT/'triage_blocked.csv',index=False)
