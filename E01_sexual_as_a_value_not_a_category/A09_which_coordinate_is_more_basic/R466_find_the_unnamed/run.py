import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A166 R466 -- 页面上有多少处「广度」没有点名

⚠ **CLOSURE**(§0.2)。不分离世界,**保护已有结论**。明确标注,不冒充发现。

`#420a`:两种广度**可以反号**(同一预测量 −0.076 vs +0.049)。
**⇒ 任何一句没点名的「广度」都可能被读成它的反面。这不是修辞问题。**
`#421` 的 NEXT 要求:**先找出来,修补留到下一轮**(「找出来」与「改对」是两件事)。

ESTIMAND        两页里每一处「广度」泛称,判定它的上下文是否**点名**了五个量之一;
                主量 = **未点名的处数**,并**逐条打印原文**。
判据(**先标支**,`#379c`)
                【两支】`lib.breadth_audit.controls()` 必须全过(`P5★`:
                        一个从未开火过的扫描器,它的每一个「都点名了」都是沉默);
                        **正对照**:`#420`/`#421` 刚写上页面的段落**已经点名**,必须被判为 named。
                【非零支】存在未点名的处 -> 报清单,交给下一轮修;
                【零支】一处都没有 -> 那本身要被怀疑(扫描器可能没在工作)。
⚠ 覆盖率       **必须和结论一起报**(`#376b`):正则找不到的不是「没问题」,是「没看」。
IMPOSSIBLE      ① 「点名」用的是**上下文窗口内是否出现具体量名**,是一个代理 ——
                   窗口内出现不等于**那一处**说的就是它(P6:`出现 ⇒ 点名` 不健全);
                   所以未点名的清单是**上界**,已点名的计数是**下界**;
                ② 窗口宽度是我选的 -> **同轮报两个宽度**,若结论只在一种上成立,是窗宽在说话(`#377a`)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.breadth_audit import scan, coverage, controls, GENERIC, SPECIFIC
c=controls(); det=c.pop('_detail')
print(f"⚠ **扫描器自检**(`P5★`):{ {k:v for k,v in c.items()} } · 细节 {det}")
assert all(c.values()), c
rows=[]
for W in (80,120,200):
    for f in ('README.md','README_zh.md'):
        t=pathlib.Path(f).read_text()
        for r in scan(t,window=W):
            rows.append(dict(v_page=f,v_win=W,v_term=r['term'],v_named=r['named'],
                             v_pos=r['pos'],v_ctx=r['ctx'][:150]))
T=pd.DataFrame(rows); check_columns(T,'R466')
T.to_csv(pathlib.Path(__file__).parent/'results'/'breadth_uses.csv',index=False)
print(f"\n⚠ **窗宽扫描(`#377a` 的教训:窗宽是我选的)**:")
for W in (80,120,200):
    S=T[T.v_win==W]; nn=len(S); un=int((~S.v_named).sum())
    print(f"   窗宽 ±{W:>3}:泛称 **{nn}** 处 · **未点名 {un}** · 已点名 {nn-un} · "
          f"未点名占 **{un/max(nn,1):.1%}**")
for f in ('README.md','README_zh.md'):
    t=pathlib.Path(f).read_text(); k,tot=coverage(t)
    print(f"   ⚠ 覆盖率 [{f}]:泛称 **{k}** 处 / 全页 **{tot:,}** 字 —— "
          f"**正则找不到的不是「没问题」,是「没看」**")
W0=120
S=T[(T.v_win==W0)&(~T.v_named)]
print(f"\n未点名的清单(窗宽 ±{W0},**逐条打印原文** —— `#374b`):")
for i,r in enumerate(S.itertuples(),1):
    print(f"\n   [{i}] {r.v_page} @ {r.v_pos} · 词 = 「{r.v_term}」")
    print(f"       …{r.v_ctx}…")
POSOK=True
for f,key in (('README.md','the endorsed-sex-act count is the outlier'),
              ('README_zh.md','而认可的性行为数是落单的那一个')):
    t=pathlib.Path(f).read_text(); i=t.find(key)
    if i<0: POSOK=False; continue
    near=[r for r in scan(t,window=W0) if abs(r['pos']-i)<300]
    if not near or not all(x['named'] for x in near): POSOK=False
print(f"\n★ **正对照**(`#421` 刚写上页面的那两段**已经点名**,必须被判为 named):"
      f"**{'通过' if POSOK else '⚠ 未通过 —— 扫描器把已点名的判成了未点名'}**")
UN=int((~T[T.v_win==W0].v_named).sum())
CONS=len({int((~T[T.v_win==w].v_named).sum())>0 for w in (80,120,200)})==1
g=Gate('页面上有多少处「广度」没有点名')
g.asserted('★【两支】扫描器自检 `controls()` 全过',all(c.values()),f"{c}",kind='control')
g.asserted('★【两支】正对照:`#421` 刚写的两段必须被判为 named',POSOK,
           f"{POSOK}",kind='control')
g.asserted('★【两支】三个窗宽给出**同一个**结论(有/没有未点名的)—— 否则是窗宽在说话',CONS,
           f"未点名数 {[int((~T[T.v_win==w].v_named).sum()) for w in (80,120,200)]}",kind='control')
g.asserted('★【非零支】存在未点名的处 -> 报清单,交给下一轮修',UN>0,
           f"窗宽 ±{W0} 下未点名 **{UN}** 处")
print(g)
print(f"\n⚠ **IMPOSSIBLE ①**:「点名」是一个**代理** —— 窗口内出现具体量名 **⇏** 那一处说的就是它。")
print(f"   **⇒ 未点名的清单是**上界**,已点名的计数是**下界**。**")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
