"""E03·A25·R132 —— 第三具仪器,而这是 `#111c` 允许的最后一次机会

**类型:FRONTIER。** `#689` 的 NEXT,判据由 `#689` 在跑之前写死,**不得改**。
**心理学的那一句:一个人答得更平,到底是「更少极端」还是「更一致」——
这两者在总变化量上一模一样,而在情欲化算子里是相反的。**

## 硬规则①(已跑,全部从对象读出,不靠码本也不靠记忆)
`data/external/` 四库逐个查的结果:**第三具仪器存在 —— `openpsych/RWAS`**。
- **22 道同刻度态度题 `Q1–Q22`,九档(1–9)** —— ⚠ **码本没写档数,是从数据读出来的**;
- `education` **四档**(1=不足高中 2=高中 3=大学 4=研究生);n 原始 **9,881**;
- ⚠ **自带假词效度检查**(`VCL6/9/12` 是假词):**24.3% 的人至少勾了一个假词**;
- **样本定义(在跑主分析之前写死)**:22 题全答 + `education∈[1,4]` + **假词一个都没勾**
  ⇒ **n = 7,258**(教育分布 754 / 2,510 / 2,337 / 1,657);
- **反向计分题由题—总分相关实测**(`#680` 的教训:不许凭记忆定极性):
  **Q4 Q6 Q8 Q9 Q11 Q13 Q15 Q18 Q20 Q21 共 10 题为负,已反向**;`|r|<0.15` 的填充题:**无**。
- ⚠ **⑤ 预注册的混淆已量:这是网络自选样本,非概率样本**
  ⇒ **结论只能写成「在自选样本上也成立」,不许升格成人群声明。**

## G1 ESTIMAND(沿用 `#689`,不改)
`ext` = 答在刻度两端(1 或 9)的题占比;`totvar` = 22 题标准化答案的人内标准差;
`disp` = `totvar` 秩对 `ext` 秩回归后的残差。**主量:`ρ(education, ext)` 与 `ρ(education, disp)`。**
## ⑧ 判据(`#689` 在跑之前写死)
**`|ρ|` 超过它自己的打乱零的 95% 分位即算「动了」**,两个分量各自判定;
**同号且都动 ⇒ 整体压缩;只有 `ext` 动 ⇒ 更少极端;只有 `disp` 动 ⇒ 更一致。**
## G2 CONTROLS
**正对照(`#689` 指定)**:`ρ(education, totvar) < 0` 且超零 ——
**⚠ 而这同时是第三具仪器上的复制检验:它失败就意味着 `#686` 有麻烦。设计允许它失败。**
**零**:打乱 `education`。
**⑤ 合成基线**:每题在人之间独立打乱 ⇒ `corr(ext, totvar)` 的必然部分,**在结论里减掉**。
## G3/G4 规格:{不过滤假词 · 过滤假词(主) · 假词+`surveyaccurate`} × 三个量 = 9 格,**全部照登**。
## KILL(条件式)
if 正对照复现: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**RWAS 是单次横断面、自选样本、教育只有四档** ⇒ 拿不到世代分层,`#683` 无法在此复制;
**RWAS 题目内容是威权主义,与 GSS/MFQ 都不同** ⇒ 本轮是「同一问题在第三具仪器上」,不是同一测量的复制;
因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_csv("data/external/openpsych/RWAS/RWAS/data.csv",sep=",",low_memory=False)
assert d.shape[1]==90, d.shape
Q=[f"Q{i}" for i in range(1,23)]; FAKE=["VCL6","VCL9","VCL12"]
REV=["Q4","Q6","Q8","Q9","Q11","Q13","Q15","Q18","Q20","Q21"]   # 由题—总分相关实测
rc=lambda a,b:float(np.corrcoef(pd.Series(np.asarray(a)).rank(),pd.Series(np.asarray(b)).rank())[0,1])
def build(sub):
    V=sub[Q].copy()
    for q in REV: V[q]=10-V[q]                     # 九档反向
    Z=pd.DataFrame({c:(V[c]-V[c].mean())/V[c].std() for c in Q})
    ext=((sub[Q]==1)|(sub[Q]==9)).mean(axis=1).to_numpy(float)
    tot=Z.std(axis=1).to_numpy(float)
    A=np.column_stack([pd.Series(ext).rank().to_numpy(),np.ones(len(ext))])
    tr=pd.Series(tot).rank().to_numpy()
    disp=tr-A@np.linalg.lstsq(A,tr,rcond=None)[0]
    return ext,tot,disp,sub.education.to_numpy(float),V
base=(d[Q]!=0).all(axis=1)&d.education.between(1,4)
SPECS={"不过滤假词":base,
       "过滤假词(主)":base&(d[FAKE].sum(axis=1)==0),
       "假词+surveyaccurate":base&(d[FAKE].sum(axis=1)==0)&(d.surveyaccurate==1)}
rng=np.random.default_rng(20260806); grid={}
print(f"{'规格':22s} {'n':>6s} {'ρ(edu,totvar)':>14s} {'ρ(edu,ext)':>11s} {'ρ(edu,disp)':>12s}")
for name,mask in SPECS.items():
    sub=d[mask]
    if len(sub)<500: print(f"{name:22s} n={len(sub)} 太小,跳过并记下"); continue
    ext,tot,disp,E,_=build(sub)
    rt,re_,rd=rc(E,tot),rc(E,ext),rc(E,disp)
    qt=float(np.quantile([abs(rc(rng.permutation(E),tot)) for _ in range(300)],.95))
    qe=float(np.quantile([abs(rc(rng.permutation(E),ext)) for _ in range(300)],.95))
    qd=float(np.quantile([abs(rc(rng.permutation(E),disp)) for _ in range(300)],.95))
    grid[name]=dict(n=int(len(sub)),r_tot=rt,r_ext=re_,r_disp=rd,q_tot=qt,q_ext=qe,q_disp=qd)
    print(f"{name:22s} {len(sub):>6,} {rt:>+9.4f}{'✅' if abs(rt)>qt else '⛔':>5s} "
          f"{re_:>+7.4f}{'✅' if abs(re_)>qe else '⛔':>4s} {rd:>+8.4f}{'✅' if abs(rd)>qd else '⛔':>4s}")
    print(f"{'':22s}        (零 95%: {qt:.4f} / {qe:.4f} / {qd:.4f})")
main=grid["过滤假词(主)"]
sub=d[SPECS["过滤假词(主)"]]; ext,tot,_,E,V=build(sub)
obs_et=rc(ext,tot); sims=[]
for _ in range(50):
    Vp=pd.DataFrame({c:rng.permutation(V[c].to_numpy()) for c in Q})
    Zp=pd.DataFrame({c:(Vp[c]-Vp[c].mean())/Vp[c].std() for c in Q})
    raw=pd.DataFrame({c:rng.permutation(sub[c].to_numpy()) for c in Q})
    sims.append(rc(((raw==1)|(raw==9)).mean(axis=1).to_numpy(),Zp.std(axis=1).to_numpy()))
forced=float(np.median(sims))
print(f"\n⑤ 合成基线:每题独立打乱 ⇒ corr(ext, totvar) 必然部分 **{forced:+.4f}** · "
      f"实测 **{obs_et:+.4f}** ⇒ 必然占 **{abs(forced)/abs(obs_et)*100:.0f}%**")
G=Gate("第三具仪器上的分解")
p1=G.positive_control("RWAS 上必须复现 ρ(education, totvar) < 0 且超零(这同时是第三具仪器的复制)",
                      planted=1.0 if (main["r_tot"]<0 and abs(main["r_tot"])>main["q_tot"]) else 0.0,
                      floor=0.5,spread=0.01)
me,md=abs(main["r_ext"])>main["q_ext"], abs(main["r_disp"])>main["q_disp"]
if p1:
    v=("**整体压缩:极端度与离散度都动且同号**" if (me and md and main["r_ext"]*main["r_disp"]>0) else
       f"**更少极端:极端度 {main['r_ext']:+.4f} 动、离散度 {main['r_disp']:+.4f} 不动 ⇒ "
       f"教育让人更少走极端,而不是更一致**" if (me and not md) else
       f"**更一致:离散度 {main['r_disp']:+.4f} 动、极端度 {main['r_ext']:+.4f} 不动**" if (md and not me) else
       f"**两者都不动:极端度 {main['r_ext']:+.4f} / 离散度 {main['r_disp']:+.4f} 都在零里 ⇒ "
       f"分解在本仪器上没有分辨力**" if not (me or md) else
       f"**两者都动但反号,照登:{main['r_ext']:+.4f} / {main['r_disp']:+.4f}**")
else: v="UNVERIFIED —— 正对照失败,而它同时意味着 `#686` 在第三具仪器上不复制"
print(f"\n{v}"); print(G)
json.dump(dict(grid=grid,forced=forced,obs_ext_tot=obs_et,rev_items=REV,verdict=v,unchallenged=True),
          open(OUT/"third_instrument.json","w"),indent=1,ensure_ascii=False)
