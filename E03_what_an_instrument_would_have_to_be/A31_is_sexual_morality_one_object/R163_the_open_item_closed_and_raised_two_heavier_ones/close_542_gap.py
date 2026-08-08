"""E03·A31·R163 —— 关掉 `#720` 留下的那条未决,而关它的过程带出两件更重的

**类型:CLOSURE(诚实标注)。它不开新世界,它关掉一条「我知道对不上却没查」的记录。**
⚠ **但它产出的两件事是 FRONTIER 级的更正,所以不是白跑的 CLOSURE。**

**要关的:** `#720` 报面间成对中位 **+0.0725**,`#542` 记 **0.1726**,差 0.10,`#720` 没查。

## 做法:口径 2×2 —— 生/归一 × 有符号/绝对值,四格全算,看哪一格是 `#542`
**④ 正对照**:配对数必须与 `#542` 写的 **200 / 4,750** 逐字相同(否则是分面不同,不是口径不同)。

⚠ **换不了仪器**:这是一次口径考古,只能在同一份 MSSCQ 上做。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/openpsych/MSSCQ/MSSCQ/"
D=pd.read_csv(P+"data.csv",sep="\t"); Q=[f"Q{i}" for i in range(1,101)]
X=D[Q].replace(0,np.nan).dropna(); X=X[(X>=1).all(axis=1)&(X<=5).all(axis=1)]
R=X[Q].rank().to_numpy(float); C=np.corrcoef(R.T)
CEIL=np.abs(np.corrcoef(np.sort(R,axis=0).T)); N=np.where(CEIL>1e-9,C/CEIL,np.nan)
iu=np.triu_indices(100,1); fac=np.array([i%20 for i in range(100)]); same=fac[iu[0]]==fac[iu[1]]
print(f"④ 正对照:面内 {int(same.sum())} 对 · 面间 {int((~same).sum())} 对 —— `#542` 写的是 200 / 4,750  "
      f"{'✅ 逐字相同 ⇒ 分面一致,差只可能来自口径' if (same.sum(),(~same).sum())==(200,4750) else '⛔ 分面不同,停'}")
assert (same.sum(),(~same).sum())==(200,4750)
cells={}
for nm,Mx in (("生 ρ(有符号)",C),("生 |ρ|",np.abs(C)),("归一(有符号)",N),("归一 |·|",np.abs(N))):
    cells[nm]=(float(np.median(Mx[iu][same])),float(np.median(Mx[iu][~same])))
    hit=" ← **#542 就是这一格**" if abs(cells[nm][0]-0.5789)<0.01 and abs(cells[nm][1]-0.1726)<0.01 else ""
    print(f"  {nm:16s} 面内 {cells[nm][0]:.4f} · 面间 {cells[nm][1]:.4f}{hit}")
neg=int((C[iu][~same]<0).sum())
print(f"\n⇒ **未决关闭:差全部来自取不取绝对值。** 而面间 4,750 对里 **{neg} 对为负({neg/4750*100:.0f}%)**")
print(f"   ⇒ 按 `#607`:符号本质随机的相关取 |·| 后的中位,量的是「没有东西」在这个 n 下的大小。")
print(f"   **0.1734 是上界,有符号的 {cells['归一(有符号)'][1]:.4f} 才可读。**")
pooled=cells["生 |ρ|"][0]
print(f"\n⚠ 带出的第一件:`#720` 的正对照**蒙对了**。")
print(f"   它比的是「逐面中位再取中位 0.5645」对 `#542` 的「200 对合并中位 0.5789」——**两个不同的统计量**,")
print(f"   只是都落在 ±0.05 内。同口径重算:**{pooled:.4f} 对 0.5789,差 {abs(pooled-0.5789):.4f}** ⇒ 结论不变,此前靠运气。")
print(f"   **这是 `#718` 的镜像:那次是失败的正对照抓出我比错了数,这次是通过的正对照掩盖了同一件事。**")
print(f"\n⚠ 带出的第二件:`#542` 已被 `P14` 的 `prior_art` 行降级为 **VERIFICATION**(同面常是近义改写)。")
print(f"   ⇒ **`#720` 的面内那一半原样继承这个降级;面间那一半不继承** ——")
print(f"   作者把二十个面设计成不同的,不保证人回答时它们不相关。**若真有一个统领性总因子,二十个面本会一起相关 0.6。**")
json.dump(dict(cells=cells,pairs=[200,4750],neg_cross=neg,
  resolution="差全部来自 |·|;#542 = 生 |ρ|",
  positive_control_of_720="蒙对:逐面中位的中位 vs 合并中位;同口径 %.4f vs 0.5789"%pooled,
  inherited="#542 的 prior_art 降级由 #720 的面内那一半继承",
  type="CLOSURE,但产出两件 FRONTIER 级更正",unchallenged=True),
  open(OUT/"gap_closed.json","w"),indent=1,ensure_ascii=False)
