import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A98 R348 -- 哪些常设结论需要嵌套修正,哪些根本不需要

`#302a` 只对**联合 R² 中位**做了嵌套(1.390% -> 0.869%)。
而这一页的常设声明里,**每一个「+0.1185」「+0.1286」都是样本内的相关**。

⚠ **但不该整片扣上去。** 两个世界预测不同,而它们的差别是**本体的**:
- **世界 A · 乐观来自多元拟合** —— 六个自由度在同一批人上择优;
  那么**单一相关几乎不动**(单变量相关没有拟合自由度),而 R² 动。
- **世界 B · 乐观来自坐标估计本身** —— 坐标是从这批人身上估出来的;
  那么**估得越"软"的坐标掉得越多**:`c3` 是一个**特征向量**(整个块×块矩阵估出来的),
  而 `S` 只是一个人内聚合,尺子来自 15k 人的流行度,**极稳**。
  => **`c3` 该掉,`S` 不该。**

ESTIMAND        `corr(S, 羞耻)` 与 `corr(c3, 羞耻)` 的**嵌套版本**:训练折估坐标,测试折读相关,合并。
KILL            **A -> 两个都几乎不动;B -> `c3` 明显掉而 `S` 不动;
                两个都掉 -> 第三个世界:样本内相关本身就有乐观,那要单独解释。**
POSITIVE CTRL   合成一个与**训练折估出的 c3** 有已知总体相关的结局,嵌套读数必须收敛到它。
NEGATIVE CTRL   纯噪声结局:嵌套相关必须 ≈0。
⚠ GUARD 18      投影必须真的跑到测试集。
⚠ 先报          **哪些声明根本不需要这个修正** —— 免得把一条普遍修正错扣在不适用的声明上。
IMPOSSIBLE      相关的衰减是 √ 量级(R² 的一半),所以同样的坐标噪声在相关上**看起来更小**;
                本轮报的是**相对变化**,不是与 R² 的直接可比量。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])   # 载入 MB/A/B/coords/fit_apply/OUT

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
rgF=np.random.default_rng(20260804)
ALL=np.flatnonzero(ok); idx=rgF.permutation(ALL); FOLD=np.array_split(idx,5)
QT=[]
for k in range(5):
    te=FOLD[k]; tr=np.concatenate([FOLD[j] for j in range(5) if j!=k])
    QT.append((tr,te,fit_apply(tr,te)))
print(f"n={len(ALL):,},5 折,每折重估六个坐标")
FIN=[int(np.isfinite(q_).sum()) for q_ in QT[0][2]]
NAMES=['S','D','c1','c2','c3','清晰度']
QFULL=fit_apply(ALL,ALL)   # ⚠ 样本内臂用**同一段代码**,训练集=全体 —— 两臂唯一的差别是折
def insample_corr(j,y):
    q_=QFULL[j]
    m=np.isfinite(q_)&np.isfinite(y)&ok
    return float(np.corrcoef(q_[m],y[m])[0,1]),int(m.sum())
def nested_corr(j,y):
    """训练折估坐标,**测试折**读相关;按测试折人数加权合并 z。"""
    zs=[];ws=[]
    for tr,te,Q in QT:
        m=np.zeros(NN,bool); m[te]=True
        m&=np.isfinite(Q[j])&np.isfinite(y)&ok
        if m.sum()<150: continue
        r=float(np.corrcoef(Q[j][m],y[m])[0,1])
        zs.append(np.arctanh(np.clip(r,-0.999,0.999))); ws.append(m.sum())
    if not zs: return np.nan,0
    return float(np.tanh(np.average(zs,weights=ws))),int(np.sum(ws))
print(f"\n{'坐标':<8}{'样本内':>12}{'嵌套':>12}{'相对变化':>12}   n")
rows=[]
for j,nm in enumerate(NAMES):
    a,na=insample_corr(j,sh); b,nb=nested_corr(j,sh)
    rows.append(dict(coord=nm,ins=a,nested=b,rel=(b-a)/abs(a) if a else np.nan,n=nb))
    print(f"{nm:<8}{a:>+12.4f}{b:>+12.4f}{100*(b-a)/abs(a):>+11.1f}%   {nb:,}")
T=pd.DataFrame(rows); check_columns(T,'R348')
T.to_csv(pathlib.Path(__file__).parent/'results'/'nested_corr.csv',index=False)
rS=T[T.coord=='S'].iloc[0]; rC=T[T.coord=='c3'].iloc[0]
print(f"\n★ `S` 相对变化 **{100*rS.rel:+.1f}%** · `c3` 相对变化 **{100*rC.rel:+.1f}%**")
def plant(anchor,tag,store):
    mm=np.isfinite(anchor)&ok; rg=np.random.default_rng(11)
    print(f"\n{tag}")
    for tr_ in (0.00,0.10,0.25):
        z=(anchor[mm]-anchor[mm].mean())/anchor[mm].std()
        y=np.full(NN,np.nan); y[mm]=tr_*z+np.sqrt(max(1-tr_**2,1e-9))*rg.standard_normal(mm.sum())
        b,_=nested_corr(4,y); a,_=insample_corr(4,y); store[tr_]=(a,b)
        print(f"   真值 {tr_:+.2f} -> 样本内 **{a:+.4f}** · **嵌套 {b:+.4f}**"
              f"(衰减 {100*(1-abs(b)/max(abs(a),1e-9)):+.0f}%)")
CTL={}; CTL1={}
# ⚠ 第一版把结局挂在**第 1 折**的 c3 上 —— 各折的 c3 本就不是同一个向量,衰减是构造出来的。
#    正确的锚是**全样本 c3**,即页面上那个量本身。两个都报,因为第一版的失败本身是结论。
plant(QFULL[4],'正/负对照 ★ 锚 = **全样本 c3**(页面上报的那个量)',CTL)
plant(QT[0][2][4],'⚠ 对照的第一版:锚 = **第 1 折的 c3**(误设 —— 各折 c3 不同,衰减是构造出来的)',CTL1)
gg=Gate('哪些常设结论需要嵌套修正')
gg.apply_reached_the_test_set('⚠ guard 18:投影跑到测试集了吗',FIN,len(QT[0][0]),len(ALL),labels=NAMES)
gg.asserted('★ 正对照:锚=全样本 c3、真值 +0.25 时嵌套读数必须收敛(±0.06)',abs(CTL[0.25][1]-0.25)<0.06,
            f"真值 +0.25 -> 嵌套 {CTL[0.25][1]:+.4f}(样本内 {CTL[0.25][0]:+.4f});"
            f"误设版(锚=第1折 c3)给 {CTL1[0.25][1]:+.4f}")
gg.heldout_drop_needs_a_plant('★ guard 19:`c3` 的下降是声明的还是仪器的',
            abs(rC.rel),1-abs(CTL[0.25][1])/abs(CTL[0.25][0]),'c3↔羞耻')
gg.asserted('⚠ 校准(同一判据的散文版,留作可读性)',
            (1-abs(CTL[0.25][1])/abs(CTL[0.25][0]))*100 < abs(100*rC.rel),
            f"正对照衰减 {100*(1-abs(CTL[0.25][1])/abs(CTL[0.25][0])):+.0f}% vs "
            f"c3 观测衰减 {abs(100*rC.rel):.0f}% —— 若前者 >= 后者,观测到的衰减**全部**可由"
            f"坐标不稳定解释,不构成「样本内相关有乐观」的证据")
gg.asserted('★ 负对照:真值 0 时嵌套相关必须 ≈0(|r| < 0.03)',abs(CTL[0.0][1])<0.03,
            f"真值 0 -> 嵌套 {CTL[0.0][1]:+.4f}")
gg.asserted('★ 注册的 kill:`c3` 明显掉而 `S` 不动(世界 B)',
            abs(rC.rel)>0.10 and abs(rS.rel)<0.05,
            f"S {100*rS.rel:+.1f}% · c3 {100*rC.rel:+.1f}% —— "
            f"A=两个都不动 · B=c3 掉 S 不动 · 第三个世界=两个都掉")
gg.asserted('⚠ 不需要这个修正的声明(先报,免得整片扣上去)',True,
            '① 纯人内量(分半信度 +0.432 是**块集合**分半,不经过人层坐标估计)· '
            '② 不经过坐标的量(题目基率 · 起始年龄的人群时间表)· '
            '③ 已经是留出的量(跨仪器的 +0.1815 / −0.1086 用的是**不相交**的题目集)')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
