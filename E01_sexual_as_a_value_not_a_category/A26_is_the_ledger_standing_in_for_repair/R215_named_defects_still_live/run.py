import os,sys,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A26 R215 -- 账本点名过的 README 缺陷,现在修好了没有

`#169c`:边跑了、输出了、没有人读。`#144d` 早就记过「两套并行叙述」,而 `+0.815` 的第二处断言
仍活到 149 条条目之后。**一条被记进账本的 README 缺陷,和一条被修好的,在账本里长得一模一样。**

ESTIMAND        账本正文点名过的 README 具体串,在当前 README 里**仍原样命中且邻近无撤回标记**的数量。
KILL            **≥3 条「记过但没修」-> 账本正在扮演修复的替身。**
POSITIVE CTRL   `+0.815`:`#169b` 刚把它三处全部改写并标记 —— 规则必须判它**已修**(marked=True)。
                规则若把它判成未修,说明标记检测失效,整轮不可读。
IMPOSSIBLE      按**字符串**匹配,不按**所指**:同一个数可以有两个所指(见 `#170b`),
                所以「同文件内重复」必须人工确认。这一步不是自动的。
"""
import pandas as pd, hashlib
import readme_ledger_audit as A
OUT=pathlib.Path(__file__).parent/'results'

def snapshot(rev):
    """把某个 revision 的两份 README 落到临时路径,对它跑同一条规则。"""
    d=OUT/f'_snap_{rev}'; d.mkdir(exist_ok=True)
    for f in ['README.md','README_zh.md']:
        try: t=subprocess.run(['git','show',f'{rev}:{f}'],capture_output=True,text=True,check=True).stdout
        except subprocess.CalledProcessError: return None
        (d/f).write_text(t)
    cwd=os.getcwd(); os.chdir(d)
    try: D=A.named_defects(ledger=str(ROOT/'RETRACTIONS.md'))
    finally: os.chdir(cwd)
    return D

def summarise(D,tag):
    if D is None or D.empty: print(f"{tag}: 无命中"); return 0,pd.DataFrame()
    live=D[~D.marked]
    dup=live.groupby(['entry','token','file']).size().reset_index(name='n')
    dup=dup[dup.n>1]
    print(f"\n{tag}: 命中 {len(D)}  无标记 {len(live)}  **同文件内重复 {len(dup)}**")
    for _,r in dup.iterrows():
        w=live[(live.entry==r.entry)&(live.token==r.token)&(live.file==r.file)]
        print(f"   #{r.entry} `{r.token}` -> {r.file} 的 {', '.join(str(x) for x in w.line)} 行各一遍")
    return len(dup),dup

before=snapshot('HEAD')          # 修之前(本轮的改动尚未提交)
after =A.named_defects()         # 修之后 = 工作区当前状态
nb,db=summarise(before,'修前(HEAD)')
na,da=summarise(after ,'修后(工作区)')
if before is not None and not before.empty: before.to_csv(OUT/'before.csv',index=False)
after.to_csv(OUT/'after.csv',index=False)

from lib.gates import Gate
g=Gate('账本是不是在扮演修复的替身')
pc=after[after.token.astype(str).str.contains('0.815')] if len(after) else pd.DataFrame()
g.asserted('正对照:`#169b` 刚改写的 +0.815 必须判为已修',
           len(pc)>0 and bool(pc.marked.all()),
           f"{len(pc)} 处命中,全部 marked={bool(pc.marked.all()) if len(pc) else '—'}")
g.asserted('可判前提:规则在修前确实读得到缺陷(否则修后的 0 无意义)',nb>0,f"修前 {nb} 处")
g.asserted('注册的 kill:修前 >=3 条「记过但没修」-> 账本在扮演修复的替身',nb>=3,
           f"修前 {nb} 处同文件内重复")
g.asserted('修完之后归零(剩余项须人工确认所指)',na<=1,f"修后 {na} 处")
print(g)
print(f"\nsha1 {hashlib.sha1(after.to_csv(index=False).encode()).hexdigest()[:12]}")
