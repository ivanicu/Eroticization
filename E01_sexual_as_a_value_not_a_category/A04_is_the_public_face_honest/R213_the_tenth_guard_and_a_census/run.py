import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A25 R213 -- 第 10 条守卫,以及一次普查

`#167b`:如果一个数的"不确定度"只在零臂里出现,而真实臂本身没有抖动来源,那它没有误差棒。
`#167` 的 NEXT 注册:加守卫 `has_error_bar`,扫全部轮次,**若超过三分之一的现行声明
展布来源是零臂 -> 不是 `#100` 一处的问题。**

ESTIMAND        ① 有多少轮次的 artifact 里真实臂**没有实现展布**;② 有多少账本条目**不作任何精度陈述**。
IDENTIFICATION  ① 从 artifact 读(仪器:error_bar_scan);② 从账本文本读(仪器:正则普查)。
                **两台仪器盲区不同** —— 这是故意的:①看不见没落盘的自助,②看不见表格里的 sd 列。
KILL            条件式:先要**正对照开火**(R147 是已知案例,扫描器必须抓到它);
                再判注册阈值。**注册的量(展布来源=零臂)对 147 轮返回 UNVERIFIED -> 那个 kill 取不到数。**
                所以第二台仪器量的是**另一个量**,必须另行标注,不得默默顶替。
POSITIVE CTRL   R147/grid.csv —— `#167a` 已确证真实臂在 3 个 seed 上逐字节相同。
IMPOSSIBLE      判不了那 147 个 UNVERIFIED 轮次。**只在缺失方向可读**(P6)。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage, check_columns
import guard_lint as GL
OUT=pathlib.Path(__file__).parent/'results'

# ---- 仪器 ①:artifact 侧 ---------------------------------------------------
S=GL.error_bar_scan('.'); S.to_csv(OUT/'artifact_scan.csv',index=False)
fake=S[S.verdict=='SEED_IS_FAKE']; nosp=S[S.verdict=='NO_SPREAD']
print(f"仪器①  artifact {len(S)} 个,轮次 {S['round'].nunique()}")
print(f"  SEED_IS_FAKE {len(fake)}   NO_SPREAD {len(nosp)}   UNVERIFIED {int((S.verdict=='UNVERIFIED').sum())}")
print(f"  SEED_IS_FAKE 的轮次:{', '.join(sorted(r[:4] for r in fake['round']))}")

# ---- 仪器 ②:账本侧 --------------------------------------------------------
txt=pathlib.Path('RETRACTIONS.md').read_text()
# ⚠ 第一版只解析 `## Entry N` 散文条目 = 151 条,而账本最大编号是 167。
#   缺的 16 条是 `#1`-`#16`,住在账本顶部的**表格**里,一行一条 —— 也就是**最早的**那批,
#   最可能没有精度陈述。漏掉它们会让"无精度份额"偏低,**偏差方向对我有利**,所以必须补。
#   这正是 #118c 的静默截断;而第一版的 check_coverage 分子分母同源 = 一个不可能失败的检查。
def feats(b):
    return dict(has_pm=bool(re.search(r'±\s*\d*\.\d+|\+/-\s*\d*\.\d+',b)),
                has_ratio=bool(re.search(r'\d+(\.\d+)?×\s*(自身|its own|展布|spread)',b)),
                has_ci=bool(re.search(r'\[[−\-\+]?\d*\.\d+\s*,\s*[−\-\+]?\d*\.\d+\]',b)),
                has_floor=bool(re.search(r'地板|floor|零臂|置换零|null',b)))
rows=[]
for e in re.split(r'\n## Entry ',txt)[1:]:
    # ⚠ `(\d+),` 要求编号后紧跟逗号,漏掉 3 条头部写法不同的条目 —— 而独立分母用的是
    #   不带逗号的正则。**分子分母用了两把不同的尺**,这是覆盖检查最容易骗过自己的方式。
    m=re.match(r'(\d+)',e)
    if m: rows.append(dict(entry=int(m.group(1)),fmt='prose',**feats(e[:4000])))
# ⚠ 不限定区域时这个正则命中账本各处的其它表(69 行,还混进一个 `0`)——
#   把它锁在**第一条散文条目之前**的那张头表里。
head_tab=txt[:txt.index('\n## Entry ')]
for line in head_tab.split('\n'):
    m=re.match(r'\|\s*(\d+)\s*\|',line)
    if m and 1<=int(m.group(1))<=16: rows.append(dict(entry=int(m.group(1)),fmt='table',**feats(line)))
T=pd.DataFrame(rows).drop_duplicates('entry'); check_columns(T,'R213 普查')
NMAX=max(int(x) for x in re.findall(r'(?m)^## Entry (\d+)',txt))   # 分母独立于解析
check_coverage(len(T),NMAX,'账本条目普查',tol=0.0)
n_pro=int((T.fmt=='prose').sum()); n_tab=int((T.fmt=='table').sum())
print(f"\n仪器②  账本条目 {len(T)}/{NMAX}(散文 {n_pro} + 表格 {n_tab})")

T['any_precision']=T.has_pm|T.has_ratio|T.has_ci
T.to_csv(OUT/'ledger_census.csv',index=False)
NO=T[~T.any_precision]; share=len(NO)/len(T)
for c in ['has_pm','has_ratio','has_ci','any_precision','has_floor']:
    print(f"  {c:<16}{T[c].sum():>4}  ({100*T[c].mean():.1f}%)")
print(f"  **不作任何精度陈述:{len(NO)}/{len(T)} = {100*share:.1f}%**,其中 {int(NO.has_floor.sum())} 条同时引用了地板/零臂")
w=T.assign(bin=(T.entry-1)//20*20+1).groupby('bin').any_precision.mean()
print("\n  按条目区间(时间顺序)的精度陈述率:")
for k,v in w.items(): print(f"   #{k:>3}-{k+19:<4} {100*v:>5.1f}%  {'█'*int(v*30)}")

g=Gate('地板与误差棒的混用是不是系统性的')
g.asserted('正对照:扫描器抓到已知案例 R147',bool(fake['round'].str.contains('R147').any()),
           f"SEED_IS_FAKE 命中 {len(fake)} 轮,含 R147")
g.asserted('仪器①覆盖了全部持久化 artifact',len(S)>0 and S['round'].nunique()>=140,
           f"{S['round'].nunique()} 轮次")
g.asserted('注册的那个量取不到数,不得由第二台仪器默默顶替',
           int((S.verdict=='UNVERIFIED').sum())>0,
           f"仪器①对 {int((S.verdict=='UNVERIFIED').sum())} 个 artifact 返回 UNVERIFIED —— "
           f"「展布来源=零臂」无法逐轮判定")
g.asserted('第二台仪器量的是另一个量(无精度陈述),超过注册的 1/3',share>1/3,
           f"{100*share:.1f}% > 33.3%")
g.threshold_outside_noise('普查份额离 1/3 有多远',share,1/3,
                          float(np.sqrt(share*(1-share)/len(T))))
print(g)
print(f"\nsha1 {hashlib.sha1((S.to_csv(index=False)+T.to_csv(index=False)).encode()).hexdigest()[:12]}")
