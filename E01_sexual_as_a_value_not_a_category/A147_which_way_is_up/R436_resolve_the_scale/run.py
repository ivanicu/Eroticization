import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A147 R436 -- 那个尺度哪一头是「更色」:`#388`–`#391` 全挂在这个符号上

`#391` 的 NEXT 让我去做无监督聚类。**做之前查了一件事,然后没走到聚类:**
`#388a` 的那五题**全部是 `RATING_NEG_FIB`,而全问卷只有 5 个这种题** ——
**我以为我从一堆里挑了五个,其实我拿的是那个类型的全部。**
而它们的取值只有 **{0, −1, −2, −3, −5, −8},从不为正**,众数是 **−8**(25–44%),**0 最罕见**(5–10%)。

⚠⚠ **那么「高 = 更色」这个我一直在用的读法,可能整个反了。**
若 **−8 = 最色**(0/1/2/3/5/8 的斐波那契强度,存成负数),则「值更高(更接近 0)」= **更不色**。
**`#388`/`#390`/`#391` 的每一句话都挂在这个符号上。**

两个活着的世界,**对同一个相关预测相反的符号**:
**A `0` = 最色**(值越高越色)-> `corr(题, 认可的性行为数)` 应当 **正**;
**B `−8` = 最色**(值越低越色)-> 应当 **负**。

ESTIMAND        用**方向已知**的锚变量定符号:
                `Totalsexacts`(认可的性行为**计数** —— 越大越色,方向由构造给出)·
                `Totalfetishcategory`(勾选的类别总数)。
                主量 = `corr(五题平均, 锚)` 的**符号**。
判据(**先标支**)
                【两支】两个锚必须**同号**(否则锚本身不可靠,不判);
                        `cunnilingus` 与 `"I find cunnilingus:"` 是同一问题的两次 -> 它们必须**高度正相关**
                        (这是一个**内部一致性**的正对照,方向已知)。
                【非零支】按符号判 A 或 B,并**当场重算** `#388a` 的那个 +0.1255 在正确符号下是什么。
⚠ 零的种类     `negative_control`:**这个零应该是零** —— 一列纯噪声与锚的相关。
                (不是 offset:锚与题之间不存在「本该有的非零基线」。)
IMPOSSIBLE      ① 计数型锚与「有多色」不是同一个构念,只用它定**符号**,不定大小;
                ② 若两个锚给出相反符号 -> 本轮不判,而那本身是关于**锚**的发现。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
print(f"`RATING_NEG_FIB` 全部 **{len(FIVE)}** 题(= `#388a` 的那五题,**一个不多一个不少**):")
for c in FIVE: print(f"   · {str(c)[:64]}")
print(f"\n取值 {sorted(set(V[np.isfinite(V)].tolist()))} · 众数 −8 · 0 最罕见\n")
def num(c): return pd.to_numeric(d[c],errors='coerce').values.astype(float)
ANCH={'Totalsexacts':num('Totalsexacts'),'Totalfetishcategory':num('totalfetishcategory')}
FIVEZ=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i])
                                  for i in range(V.shape[1])]),1)
print("① 用**方向已知**的锚定符号(锚:计数,越大越色 —— 方向由构造给出):")
res={}
for nm,a in ANCH.items():
    g=np.isfinite(FIVEZ)&np.isfinite(a)
    r=float(np.corrcoef(FIVEZ[g],a[g])[0,1]); res[nm]=r
    print(f"   corr(五题平均, {nm:<22}) = **{r:+.4f}**  (n={int(g.sum()):,})")
SAME=np.sign(list(res.values())[0])==np.sign(list(res.values())[1])
rbar=float(np.mean(list(res.values())))
print(f"   两锚同号 **{SAME}** · 平均 **{rbar:+.4f}** -> "
      f"**世界 {'A(0 = 最色,值越高越色)' if rbar>0 else 'B(−8 = 最色,值越**低**越色)'}**")
i1=[i for i,c in enumerate(FIVE) if str(c).strip()=='cunnilingus'][0]
i2=[i for i,c in enumerate(FIVE) if 'jn2b355' in str(c)][0]
g2=np.isfinite(V[:,i1])&np.isfinite(V[:,i2])
rdup=float(np.corrcoef(V[g2,i1],V[g2,i2])[0,1])
print(f"\n② 内部一致性正对照(方向已知):`cunnilingus` 与 `\"I find cunnilingus:\"` 是**同一问题的两次** -> "
      f"r = **{rdup:+.4f}**")
rgN=np.random.default_rng(7)
gA=np.isfinite(ANCH['Totalsexacts'])
rF=float(np.corrcoef(rgN.standard_normal(int(gA.sum())),ANCH['Totalsexacts'][gA])[0,1])
print(f"   负对照(纯噪声与锚):**{rF:+.4f}**(**这个零应该是零** -> `negative_control`,不是 offset)")
print(f"\n③ ⇒ 在**正确符号**下,`#388a`/`#390a` 的那个 +0.1255 该怎么读:")
SIGN=+1 if rbar>0 else -1
print(f"   我一直把「值更高」读成「更色」。锚说 **{'那是对的' if SIGN>0 else '那是**反的**'}**。")
if SIGN<0:
    print(f"   **⇒ 「常规的东西对他也很色」要改成「常规的东西对他**没那么**色」。**")
    print(f"   **⇒ 页面上 `#388`/`#390`/`#391` 的那句人话,方向要整个翻过来。**")
pd.DataFrame([dict(v_anchor=k,v_r=v) for k,v in res.items()]+
             [dict(v_anchor='dup_cunnilingus',v_r=rdup),dict(v_anchor='noise',v_r=rF)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'sign.csv',index=False)
g=Gate('那个尺度哪一头是「更色」')
g.asserted('★【两支】两个方向已知的锚必须同号(否则锚不可靠,不判)',SAME,
           ' · '.join(f"{k} {v:+.4f}" for k,v in res.items()),kind='control')
g.asserted('★【两支】内部一致性正对照:同一问题的两次必须高度正相关',rdup>0.5,
           f"r = {rdup:+.4f}",kind='control')
g.asserted('★【两支】负对照:纯噪声与锚 ≈ 0',abs(rF)<0.05,f"{rF:+.4f}",kind='control')
if SAME and rdup>0.5 and abs(rF)<0.05:
    g.asserted('★【非零支】我一直用的读法(值越高越色)是对的',rbar>0,
               f"锚平均 {rbar:+.4f} -> 世界 {'A' if rbar>0 else 'B'}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
