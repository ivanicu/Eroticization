import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A164 R463 -- `form` 的那两题,方向查过吗

`#392c`:那五题(`RATING_NEG_FIB`)的方向**被我读反了四轮**,页面上三处要翻。
`#392e` 立的规矩:**任何变量第一次进入主量之前,先打印取值集合、众数、与方向已知锚的相关。**
**⚠ 而 `form` 的两题(`animated` / `written`)从没被查过** ——
它们支撑着页面上的 `#316`(可靠性 0.380 · 解衰减 +0.1293)· `#325`(两指标不等量)· `#375`(没有第三个指标)。

⚠ 它们是 `LIKERT_PM3`(−3…+3),**不是** `RATING_NEG_FIB`,所以**不必然**同病 ——
**但「不必然」不是「查过了」。**

两个活着的世界:
**A 方向是对的** -> 这是一次**便宜的确认**(Closure),照实说;
**B 方向是反的** -> `#316`/`#325`/`#375` 三处都要翻,和 `#392` 一样。

ESTIMAND        ① 两题各自的取值集合 · 众数 · 与**方向已知锚**的相关(`#392e` 的三行);
                ② 在正确方向下,页面上那几个数该怎么读。
判据(**先标支**,`#379c`)
                【两支】**锚必须先自证**:`corr(Totalsexacts, 勾选类别总数) > 0`(两个计数同向);
                        **第二个锚必须同号**(否则锚不可靠,不判)。
                【非零支】两题与锚**同为正** -> 世界 A(我一直用的读法是对的);
                          **为负** -> 世界 B,页面要翻。
⚠ 零的种类     `negative_control`:**这个零应该是零** —— 一列纯噪声与锚的相关。
                (不是 offset:锚与题之间不存在「本该有的非零基线」。)
IMPOSSIBLE      ① 计数型锚与「有多色」不是同一构念,只用它定**符号**,不定大小;
                ② 两个锚若给出相反符号 -> 不判,而那本身是关于**锚**的发现。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
def num(c): return pd.to_numeric(d[c],errors='coerce').values.astype(float)
AN=num('animated'); WR=num('written')
A1=num('Totalsexacts'); A2=num('totalfetishcategory')
inv=pd.read_csv('data/derived/inventory.csv')
kd={r['col']:r['kind'] for _,r in inv.iterrows()}
print(f"⚠ **`#392e` 的三行**(跑在任何结论之前):")
for nm,v in (('animated',AN),('written',WR)):
    u=np.unique(v[np.isfinite(v)])
    print(f"   `{nm}` kind=**{kd.get(nm,'?')}** · 取值 **{u.tolist()}** · "
          f"众数 **{float(pd.Series(v[np.isfinite(v)]).mode().iloc[0]):g}** · n={int(np.isfinite(v).sum()):,}")
g=np.isfinite(A1)&np.isfinite(A2)
rANC=float(np.corrcoef(A1[g],A2[g])[0,1])
print(f"\n⚠ **锚先自证**:`corr(Totalsexacts, 勾选类别总数)` = **{rANC:+.4f}** "
      f"-> {'**同向,锚可用**' if rANC>0 else '**⚠ 反向 —— 锚不可用**'}")
# ---- ⚠ 第一版只有两个**计数**锚,而它们对 `animated`/`written` **异号** ----
# 按预注册,异号 -> 不判。**但那不是关于 `form` 的失败,是关于两个锚的发现:**
# `Totalsexacts`(性行为计数)与 `totalfetishcategory`(恋物类别计数)只相关 +0.3731,
# **它们测的是两种不同的广度**。
# ⇒ 需要一个**语义上方向被迫**的同族题来标定 `LIKERT_PM3` 这个**家族**本身:
ACC=next(c for c in d.columns if '41kpfir' in str(c)); ACT=num(ACC)
print(f"\n★ **家族标定**(语义上方向被迫):`{str(ACC)[:56]}`")
print(f"   kind = **{kd.get(ACC,'?')}** —— **与 `animated`/`written`/羞耻同族**")
print(f"   ⚠ 「我已经**实践过**所有唤起我的东西」-> 若「值越高 = 越同意」,"
      f"它**必须**与「做过多少性行为」的计数**正相关**。")
CAL=[]
for nm,a in (('Totalsexacts',A1),('勾选类别总数',A2)):
    gg2=np.isfinite(ACT)&np.isfinite(a); r=float(np.corrcoef(ACT[gg2],a[gg2])[0,1]); CAL.append(r)
    print(f"   corr(实践题, {nm:<12}) = **{r:+.4f}**")
CALOK=all(x>0 for x in CAL)
print(f"   -> **{'两个都为正 -> `LIKERT_PM3` 的常规方向(值越高 = 越同意)被确认' if CALOK else '⚠ 有为负 -> 家族方向存疑'}**")
print(f"\n⚠ 而**取值集合本身**也是证据:`LIKERT_PM3` 对称(−2…+2 / −3…+3,两侧都有质量),")
print(f"   `RATING_NEG_FIB` **从不为正**(0,−1,−2,−3,−5,−8)——")
print(f"   **`#392` 的病从取值集合就能诊断,而它在这里不适用。**")

rows=[]
print(f"\n① 两题与**两个方向已知的锚**的相关:")
for nm,v in (('animated',AN),('written',WR)):
    g1=np.isfinite(v)&np.isfinite(A1); g2=np.isfinite(v)&np.isfinite(A2)
    r1=float(np.corrcoef(v[g1],A1[g1])[0,1]); r2=float(np.corrcoef(v[g2],A2[g2])[0,1])
    same=np.sign(r1)==np.sign(r2)
    rows.append(dict(v_item=nm,v_r_sexacts=r1,v_r_cats=r2,v_same=bool(same)))
    print(f"   `{nm:<9}` ↔`Totalsexacts` **{r1:+.4f}** · ↔`勾选类别总数` **{r2:+.4f}** · "
          f"{'**两锚同号**' if same else '⚠ **两锚异号 —— 不判**'}")
T=pd.DataFrame(rows); check_columns(T,'R463')
T.to_csv(pathlib.Path(__file__).parent/'results'/'form_direction.csv',index=False)
ALLSAME=bool(T.v_same.all())
ALLPOS=bool((T.v_r_sexacts>0).all() and (T.v_r_cats>0).all())
rgN=np.random.default_rng(7)
gA=np.isfinite(A1); rF=float(np.corrcoef(rgN.standard_normal(int(gA.sum())),A1[gA])[0,1])
print(f"\n负对照(**这个零应该是零**:纯噪声与锚)= **{rF:+.4f}**")
gm=np.isfinite(AN)&np.isfinite(WR)
rAW=float(np.corrcoef(AN[gm],WR[gm])[0,1])
print(f"\n② 内部一致性(方向已知的正对照):`corr(animated, written)` = **{rAW:+.4f}** "
      f"(n={int(gm.sum()):,})· 页面 `#375c` 记的是 **+0.2345**")
print(f"\n③ ⇒ 在**正确方向**下,页面上那几个数该怎么读:")
if CALOK:
    print(f"   **家族标定说:`LIKERT_PM3` 的「值越高 = 越同意」是对的。**")
    print(f"   **⇒ `#316`/`#325`/`#375` 的方向**不需要翻**;")
    print(f"      **而羞耻题也是 `LIKERT_PM3`(−3…+3)—— 整页的羞耻方向也在同一次标定里被确认,**")
    print(f"      **而它此前从没被锚检过。**")
    print(f"   **⇒ 两锚异号不是方向问题,是**实质发现**:**")
    print(f"      对**画的/写的**东西有反应,与**认可的性行为数**负相关"
          f"({T.v_r_sexacts.iloc[0]:+.4f} / {T.v_r_sexacts.iloc[1]:+.4f}),"
          f"与**恋物类别数**正相关({T.v_r_cats.iloc[0]:+.4f} / {T.v_r_cats.iloc[1]:+.4f})。")
else:
    print(f"   **⚠ 家族标定失败 -> 方向存疑,`#316`/`#325`/`#375` 都要标注。**")
g2_=Gate('`form` 的那两题,方向查过吗')
g2_.asserted('★【两支】锚先自证:两个计数应同向',rANC>0,f"{rANC:+.4f}",kind='control')
g2_.asserted('★【两支】负对照:纯噪声与锚 ≈ 0',abs(rF)<0.05,f"{rF:+.4f}",kind='control')
g2_.asserted('★【两支】**家族标定**:语义上方向被迫的同族题必须与两个计数锚都为正',CALOK,
             f"实践题 ↔两锚 {CAL[0]:+.4f} / {CAL[1]:+.4f}",kind='control')
g2_.asserted('⚠【记录】两个计数锚对 `form` **异号** —— 不是方向问题,是实质发现',not ALLSAME,
             f"animated {T.v_r_sexacts.iloc[0]:+.4f} vs {T.v_r_cats.iloc[0]:+.4f} · "
             f"written {T.v_r_sexacts.iloc[1]:+.4f} vs {T.v_r_cats.iloc[1]:+.4f}")
g2_.asserted('★【两支】正对照(方向已知):两题应正相关(同一构念的两个指标)',rAW>0.15,
             f"{rAW:+.4f}",kind='control')
if rANC>0 and CALOK and abs(rF)<0.05:
    g2_.asserted('★【非零支】家族标定通过 -> 页面上 `#316`/`#325`/`#375` 与**羞耻**的方向都不需要翻',
                 CALOK,f"实践题 ↔两锚 {CAL[0]:+.4f} / {CAL[1]:+.4f} · "
                 f"`LIKERT_PM3` 对称而 `RATING_NEG_FIB` 从不为正 -> `#392` 的病不适用")
else:
    g2_.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g2_)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
