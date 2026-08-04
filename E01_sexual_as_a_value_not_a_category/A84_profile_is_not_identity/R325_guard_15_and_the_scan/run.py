import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A84 R325 -- 第 15 个守卫:剖面相似不是构念同一;并回头扫这个项目用过它几次

`#279b`:三次撞上同一个分离,第三次是**代数上保证**的
(`form_i` 与 D 分数层恰好正交,剖面相似 +0.7105)。

GATE            四端自检(独立 `Gate` 实例,`#249c`),只把「四端是否全对」汇报进主门:
                ① 剖面高 + 分数低 → **FAIL**(新)
                ② 剖面高 + 分数高 → 放行
                ③ 剖面低 → 放行(「不是同一个东西」这个方向不受限制)
                ④ **本轮这个代数保证的反例**(分数恰好 0,剖面 +0.7105)→ FAIL
外加             普查:账本里所有以**剖面相关**下过结论的条目。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

SELF=Gate('守卫 15 —— 四端自检(独立实例,#249c)')
CASES=[('①剖面高+分数低',0.7826,0.1589,False),
       ('②剖面高+分数高',0.7826,0.6500,True),
       ('③剖面低(方向不受限制)',0.1235,0.0100,True),
       ('④`#324` 的代数反例(分数恰好 0)',0.7105,0.0000,False)]
got=[SELF.profile_similarity_is_not_identity(n,pr,sr) for n,pr,sr,_ in CASES]
exp=[e for _,_,_,e in CASES]; allok=(got==exp)
print(SELF)
print(f"四端:期望 {exp} · 实测 {got} -> **{'全对' if allok else '不全对'}**")
LED=(ROOT/'RETRACTIONS.md').read_text()
ents=re.split(r'\n## Entry ',LED)[1:]
KEY=('剖面相关','剖面层','剖面相似','剖面的相关')
hits=[]
for e in ents:
    num=re.match(r'(\d+)',e)
    if not num: continue
    if any(k in e for k in KEY):
        has_score=('分数层' in e or '分数 ' in e)
        hits.append(dict(entry=int(num.group(1)),has_score_level=has_score,
                         snippet=re.sub(r'\s+',' ',e[:90])))
T=pd.DataFrame(hits); check_columns(T,'R325')
T.to_csv(pathlib.Path(__file__).parent/'results'/'profile_conclusions.csv',index=False)
bare=T[~T.has_score_level]
print(f"\n普查:账本 {len(ents)} 条,提到**剖面相关**的 **{len(T)}** 条;"
      f"其中**同时给了分数层的 {int(T.has_score_level.sum())} 条**,"
      f"**只给剖面、没给分数层的 {len(bare)} 条**")
if len(bare): print("  ⚠ 只给剖面的条目:"+' · '.join(f"#{int(x)}" for x in bare.entry))
print(f"  ✅ 已同时给分数层的:"+' · '.join(f"#{int(x)}" for x in T[T.has_score_level].entry))
g=Gate('守卫 15 + 普查')
g.asserted('⚠ 类型:CLOSURE —— 保护的是「剖面相似不能当同一构念的证据」',True,'§0 三类动作')
g.asserted('★ 四端全对(独立 Gate 实例;④ 是 `#324` 的代数保证反例,不是模拟)',
           allok, f"期望 {exp} · 实测 {got}")
g.asserted('★ 普查给出计数与名单',
           len(T)>0, f"{len(T)} 条提到剖面;只给剖面未给分数层的 {len(bare)} 条")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
