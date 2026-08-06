import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A27 R216 -- 去重删掉的是重复,还是内容

`#170c` 排出三层(守卫没调用 / 输出没人读 / 记了没人修),**第四层没查:修得对不对。**
`#170a` 把 R06 正文段的数字删掉只留指针 —— 而没有任何检查能证明删掉的那些数
与声明表里留下的那些数是**同一批**。

ESTIMAND        修前 README 里出现过的**数量**,在修后整份文件里一次都不出现的集合。
IDENTIFICATION  参照提交**钉成 SHA**(`88b0993` = 去重提交 `7b1951f` 的父),
                否则每提交一次 `HEAD~1` 就漂一格,而这个检查就变成在查别的东西。
KILL            **存在任何一个离开这一页的数量 -> 去重删掉的是内容,不是重复。**
POSITIVE CTRL   规则必须在**修复之前**读到那个洞(否则修后的 0 无意义)——
                本轮直接对 `7b1951f` 的工作区状态跑一遍作对照。
IMPOSSIBLE      一个数可以被**改写**而不是删除。命中项必须逐条人工判。
"""
import pandas as pd, hashlib
import readme_ledger_audit as A
OUT=pathlib.Path(__file__).parent/'results'
PRE='88b0993'; DEDUP='7b1951f'

before=A.numbers_that_left(rev=PRE)          # 钉死的父提交 -> 工作区(已修)
print(f"参照 {PRE} (去重提交 {DEDUP} 的父) -> 工作区")
if before.empty: print("  没有数量离开这一页")
else:
    for _,r in before.iterrows(): print(f"  {r.file} `{r.token}`  {r.old_excerpt[:100]}")
before.to_csv(OUT/'after_fix.csv',index=False)

# 正对照:同一条规则,对**去重刚落地时**的那一版跑
import subprocess, tempfile
snap=OUT/'_ctrl'; snap.mkdir(exist_ok=True)
for f in ['README.md','README_zh.md']:
    (snap/f).write_text(subprocess.run(['git','show',f'{DEDUP}:{f}'],
                                       capture_output=True,text=True,check=True).stdout)
cwd=os.getcwd(); os.chdir(snap)
try: ctrl=A.numbers_that_left(rev=PRE)
finally: os.chdir(cwd)
print(f"\n正对照 —— 同一条规则对去重刚落地那一版({DEDUP})跑:")
for _,r in ctrl.iterrows(): print(f"  {r.file} `{r.token}`  {r.old_excerpt[:100]}")
ctrl.to_csv(OUT/'control_at_dedup.csv',index=False)

from lib.gates import Gate
g=Gate('去重删掉的是重复还是内容')
g.asserted('正对照:规则在修复之前确实读到了洞',len(ctrl)>0,
           f"{len(ctrl)} 个数量离开:{list(ctrl.token) if len(ctrl) else '—'}")
g.asserted('注册的 kill:去重刚落地时有数量离开这一页 -> 删的是内容',len(ctrl)>0,
           f"{list(ctrl.token) if len(ctrl) else '无'}")
g.asserted('修复之后归零',len(before)==0,f"修后 {len(before)} 个")
print(g)
print(f"\nsha1 {hashlib.sha1((before.to_csv(index=False)+ctrl.to_csv(index=False)).encode()).hexdigest()[:12]}")
