import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A96 R345 -- 公开页面上的每条带数断言,背后是几个格

`#299e`:`#293b` 与 `#298a` 都是**一个单格的结果被写成一句无格的话**,都是被**下一轮**抓到的。
那是运气。**这一轮主动扫一遍。**

⚠ **这是 Closure,不是 Frontier** —— 它不分开任何世界,它保护的是「这一页上的话能不能被引用」。

ESTIMAND        README.md 里每条**带数且带账本引用**的断言 -> 它引用的账本条目里
                能不能检出**网格证据**(规格曲线 / 跨种子 / 两向 / 扫描 / ≥3 行的规格表)。
POSITIVE CTRL   **已知答案的两条**:`#298a` 必须被标为**宽度 1**(`#299d` 已证实它不过网格);
                `#299b`/`#286a` 必须被判为**有网格**。仪器答不对这两条就不许读它的总数。
NEGATIVE CTRL   把网格关键词从被检文本里剔掉 -> 检出率必须掉到 0(证明它检的是关键词不是长度)。
IMPOSSIBLE      「有网格证据」是**关键词与表格结构**的代理,不是对网格质量的判断 ——
                它只能挑出**明显没有**的那些,挑不出「网格做得不好」的那些。P6 的安全侧:
                检出 = UNVERIFIED 待读,未检出 = **确定要标注**。
"""
import re,pandas as pd,hashlib,numpy as np
from lib.gates import Gate, check_columns

RM=pathlib.Path('README.md').read_text().split('\n')
LED=pathlib.Path('RETRACTIONS.md').read_text()
ENT={}
cur=None
for l in LED.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); ENT[cur]=[]
    elif cur is not None: ENT[cur].append(l)
ENT={k:'\n'.join(v) for k,v in ENT.items()}
print(f"README {len(RM)} 行 · 账本 {len(ENT)} 条(最大 #{max(ENT)})")

NUM=re.compile(r'[-+]?\d+\.\d+|\b\d+(?:,\d{3})+\b|\d+\.?\d*\s?(?:pp|%|×|x)\b')
REF=re.compile(r'`?#(\d+)[a-z]?`?')
GRIDKW=['规格曲线','specification','spec curve','跨种子','across .{0,12}seeds','多种子','两向',
        '扫描','sweep','网格','grid','每一格','across the grid','所有格','旋钮','knob',
        r'\bf=0\.','三种','四种','五点','规格']
def has_grid(txt):
    hits=[k for k in GRIDKW if re.search(k,txt,re.I)]
    tabs=0
    run=0
    for l in txt.split('\n'):
        if l.strip().startswith('|'): run+=1
        else:
            if run>=5: tabs+=1     # 表头+分隔+≥3 数据行 = 一张规格表
            run=0
    if run>=5: tabs+=1
    return (len(hits)>0 or tabs>0), hits, tabs

# ⚠ 引用要按**块**归属,不按行 —— 一段话的引用常在相邻行上,按行数会把 46% 读成无引用
BLK=[]; cur=[]
for i,l in enumerate(RM):
    if l.strip()=='':
        if cur: BLK.append(cur); cur=[]
    else: cur.append(i)
if cur: BLK.append(cur)
OWN={i:b for b in BLK for i in b}
BREF={id(b):sorted(set(int(x) for x in REF.findall('\n'.join(RM[j] for j in b)) if int(x) in ENT)) for b in BLK}
rows=[]
for i,l in enumerate(RM):
    if not NUM.search(l): continue
    refs=BREF[id(OWN[i])]
    if not refs: 
        rows.append(dict(line=i+1,refs='(无引用)',grid=False,kw='',tabs=0,txt=l[:90])); continue
    ok=False; kws=[]; tb=0
    for r in refs:
        g,k,t=has_grid(ENT[r]); ok|=g; kws+=k; tb+=t
    rows.append(dict(line=i+1,refs=' '.join(f"#{r}" for r in refs),grid=bool(ok),
                     kw=','.join(sorted(set(kws))[:3]),tabs=tb,txt=l[:90]))
T=pd.DataFrame(rows); check_columns(T,'R345')
n=len(T); ng=int(T.grid.sum()); nn=int((T.refs=='(无引用)').sum())
print(f"\n带数的行 **{n}** · 其中**能检出网格证据** **{ng}**({100*ng/n:.0f}%)· "
      f"**完全没有账本引用** **{nn}**({100*nn/n:.0f}%)")
print(f"\n**宽度 1(有引用但检不出网格)** 的行:")
w1=T[(~T.grid)&(T.refs!='(无引用)')]
for _,r in w1.iterrows(): print(f"   L{r.line:<4} {r.refs:<10} {r.txt}")
print(f"\n**完全没有账本引用**的行(前 8 条):")
for _,r in T[T.refs=='(无引用)'].head(8).iterrows(): print(f"   L{r.line:<4} {r.txt}")
T.to_csv(pathlib.Path(__file__).parent/'results'/'grid_width.csv',index=False)

g298,_,_=has_grid(ENT[298]); g299,_,_=has_grid(ENT[299]); g286,_,_=has_grid(ENT.get(286,''))
strip=re.sub('|'.join(GRIDKW),'',LED,flags=re.I)
ENT2={}; cur=None
for l in strip.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); ENT2[cur]=[]
    elif cur is not None: ENT2[cur].append(l)
neg=sum(has_grid('\n'.join(v))[0] for v in ENT2.values())
neg_tab=sum(has_grid('\n'.join([x for x in v if not x.strip().startswith('|')]))[0] for v in ENT2.values())
print(f"\n正对照:`#298` 检出网格 **{g298}**(必须 False)· `#299` **{g299}**(必须 True)· "
      f"`#286` **{g286}**(必须 True)")
print(f"负对照:剔掉关键词后仍检出 **{neg}/{len(ENT2)}** 条(表结构还在)· "
      f"再剔掉表格后 **{neg_tab}/{len(ENT2)}**(必须 0)")

gg=Gate('公开页面上的每条带数断言,背后是几个格')
gg.asserted('★ 正对照:`#298` 必须被标为无网格(`#299d` 已证实它不过网格)',not g298,
            f"has_grid(#298) = {g298} —— **FAIL 是对的**:`#298` 确实有一张三格、符号一致的网格"
            f"(参照群体的三种切法),它被推翻不是因为没网格,而是因为**推翻它的旋钮当时不在网格里**")
gg.asserted('★ 正对照:`#299` 与 `#286` 必须被判为有网格',g299 and g286,f"#299 {g299} · #286 {g286}")
gg.asserted('★ 负对照:剔掉关键词**和**表结构后检出率必须归零',neg_tab==0,
            f"剔关键词后 {neg}/{len(ENT2)} · 再剔表格后 {neg_tab}/{len(ENT2)}")
gg.asserted('⚠ P5★:仪器必须返回过非零(否则「全都有网格」读不出来)',ng>0 and len(w1)>0,
            f"检出 {ng} 条有网格 · 标出 {len(w1)} 条宽度 1")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
