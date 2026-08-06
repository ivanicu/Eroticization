"""E03·A25·R133 —— `#687` 那个 2.4 倍,到底是「更极端」还是「更一致」的?

**类型:FRONTIER。** `#690` 的 NEXT。
`#687` 量到「压缩在宽容侧强 2.4 倍」,而 `#690` 证明「压缩」是**两件反向的事**压在一起
⇒ **那 2.4 倍必须归到其中一件头上,否则它是一句没有指称的话。**

## 硬规则①(已跑),而它先封死了一条路
`#687` ⑤ 写死:`|corr(stance, totvar)| > 0.3` ⇒ 分侧受污染,须改用「仅由外部题定义的立场」。
逐具仪器实测:

| 仪器 | `corr(stance, totvar)` | 能否分侧 | `ext` 可算? |
|---|---|---|---|
| GSS | **−0.2189** ✅ | 可 | ⛔ **27/33 二值**(`#689` 已登记) |
| **MFQ** | **−0.1650** ✅ | **可** | ✅ 34 题六档 |
| RWAS | **+0.7656** ⛔ **反号且四倍** | **不可** | ✅ 22 题九档 |

**⇒ RWAS 上这条路结构性封死**:`stance` 与 `totvar` 相关 0.77,而 **RWAS 的 22 题就是全部仪器,
没有外部题可用来另定义 `stance`** —— 预注册的改道方案在这具仪器上不存在。
**⇒ 只有 MFQ 同时满足两个条件。** n = 3,820,两侧各 **1,910 ≥ 预注册的 1,000 底线**。
⚠ **而那张表本身是一条发现:同一个量在三具仪器上是 −0.22 / −0.17 / +0.77,
RWAS 反号且大四倍 —— 它是一条带反向计分的单一意识形态量表,反向计分本身制造了这个耦合。**

## G1 ESTIMAND
在 MFQ 的**高立场侧**与**低立场侧**各自算 `ρ(educ_num, ext)` 与 `ρ(educ_num, disp)`,n = 人数。
## ⑧ 判据(`#690` 在跑之前写死,不得改)
**两侧的 `ext` 与 `disp` 若都同号且量级比在 [0.5, 2] 内 ⇒ 与方向无关;
若某一分量的比值落在 [0.5,2] 之外 ⇒ 那一分量是有方向的,而另一分量不是 —— 这两句必须分开写。**
## G2 CONTROLS
**正对照**:`ρ(educ_num, totvar) < 0` 且超零(`#689` 已在全样本上建立 −0.0542)。
**安慰剂(沿用 `#687` 的两个)**:纯随机分侧;以及**与 `totvar` 相关度匹配的 sham 分侧**
—— 后者复制选择结构而不携带立场内容 ⇒ `offset_control`,**零的种类 = 相关度匹配的随机分侧下的同一比值**。
**⑤ 合成基线**:**同一次打乱**同时算 `stance` 与 `ext`,两侧各自减掉必然部分。
## KILL(条件式)
if 正对照复现 and 纯随机分侧的比值≈1: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**每侧 n = 1,910 ⇒ 功率只有全样本的一半**,零的 95% 分位约 0.045,
**而效应量在 0.05–0.10 量级 ⇒ 本设计只能分辨大的两侧差异,分辨不了小的**;
**RWAS 上此分侧结构性不可能**(见上);**GSS 上 `ext` 结构性不可算**(`#689`);
MFQ 是网络自选样本 ⇒ 结论只在自选样本上成立;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
d,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
ITEMS=[c for c in ["harm","emotionally","weak","cruel","compassion","animal","kill",
  "fairly","unfairly","treated","justice","rights","rich","loyalty","betray","yourgroup",
  "lovecountry","family","team","history","duties","traditions","respect","chaos","kidrespect",
  "soldier","shutup","disgusting","decency","desires","god","harmlessdg","unnatural","chastity"]
  if c in d.columns and d[c].nunique()>=5]
j=d.dropna(subset=ITEMS+["educ_num"]); j=j[j.educ_num<=7].reset_index(drop=True)
Z=pd.DataFrame({c:(j[c]-j[c].mean())/j[c].std() for c in ITEMS})
stance=Z.mean(axis=1).to_numpy(float); tot=Z.std(axis=1).to_numpy(float)
lo_,hi_=j[ITEMS].min().min(),j[ITEMS].max().max()
ext=((j[ITEMS]==lo_)|(j[ITEMS]==hi_)).mean(axis=1).to_numpy(float)
A=np.column_stack([pd.Series(ext).rank().to_numpy(),np.ones(len(ext))])
tr=pd.Series(tot).rank().to_numpy(); disp=tr-A@np.linalg.lstsq(A,tr,rcond=None)[0]
E=j.educ_num.to_numpy(float)
rc=lambda a,b:float(np.corrcoef(pd.Series(np.asarray(a)).rank(),pd.Series(np.asarray(b)).rank())[0,1])
rng=np.random.default_rng(20260806)
qt=float(np.quantile([abs(rc(rng.permutation(E),tot)) for _ in range(300)],.95))
print(f"n = {len(j):,} · 正对照 ρ(educ, totvar) = **{rc(E,tot):+.4f}** (零 95% {qt:.4f}) "
      f"{'✅' if (rc(E,tot)<0 and abs(rc(E,tot))>qt) else '⛔ 当场停'}")
def sides(sel):
    out={}
    for nm,m in [("高立场侧",sel),("低立场侧",~sel)]:
        re_,rd=rc(E[m],ext[m]),rc(E[m],disp[m])
        qe=float(np.quantile([abs(rc(rng.permutation(E[m]),ext[m])) for _ in range(200)],.95))
        qd=float(np.quantile([abs(rc(rng.permutation(E[m]),disp[m])) for _ in range(200)],.95))
        out[nm]=dict(n=int(m.sum()),ext=re_,disp=rd,q_ext=qe,q_disp=qd)
    return out
m=np.median(stance); obs=sides(stance>m)
print(f"\n{'侧':10s} {'n':>6s} {'ρ(edu,ext)':>12s} {'ρ(edu,disp)':>13s}")
for k,v in obs.items():
    print(f"{k:10s} {v['n']:>6,} {v['ext']:>+8.4f}{'✅' if abs(v['ext'])>v['q_ext'] else '⛔':>4s} "
          f"{v['disp']:>+9.4f}{'✅' if abs(v['disp'])>v['q_disp'] else '⛔':>4s}   (零 {v['q_ext']:.4f} / {v['q_disp']:.4f})")
re_ratio=abs(obs["高立场侧"]["ext"])/abs(obs["低立场侧"]["ext"])
rd_ratio=abs(obs["高立场侧"]["disp"])/abs(obs["低立场侧"]["disp"])
print(f"\n量级比:**ext {re_ratio:.3f}** {'✅ 在 [0.5,2]' if 0.5<=re_ratio<=2 else '⛔ 出界 ⇒ 这一分量是有方向的'} · "
      f"**disp {rd_ratio:.3f}** {'✅ 在 [0.5,2]' if 0.5<=rd_ratio<=2 else '⛔ 出界 ⇒ 这一分量是有方向的'}")
ra=[];rb=[]
target=rc(stance,tot)
for _ in range(200):
    s=rng.permutation(stance>m); o=sides(s)
    ra.append((abs(o["高立场侧"]["ext"])/abs(o["低立场侧"]["ext"]),
               abs(o["高立场侧"]["disp"])/abs(o["低立场侧"]["disp"])))
    g=rng.standard_normal(len(tot)); trn=(tr-tr.mean())/tr.std()
    sham=target*trn+np.sqrt(max(0.,1-target**2))*g; o2=sides(sham>np.median(sham))
    rb.append((abs(o2["高立场侧"]["ext"])/abs(o2["低立场侧"]["ext"]),
               abs(o2["高立场侧"]["disp"])/abs(o2["低立场侧"]["disp"])))
ra=np.array(ra); rb=np.array(rb)
print(f"安慰剂 A(纯随机分侧)比值中位 ext **{np.nanmedian(ra[:,0]):.3f}** · disp **{np.nanmedian(ra[:,1]):.3f}**")
print(f"安慰剂 B(相关度匹配 sham)比值中位 ext **{np.nanmedian(rb[:,0]):.3f}** · disp **{np.nanmedian(rb[:,1]):.3f}**")
G=Gate("那 2.4 倍是哪一件的")
p1=G.positive_control("MFQ 上 ρ(educ, totvar) < 0 且超零",planted=1.0 if (rc(E,tot)<0 and abs(rc(E,tot))>qt) else 0.0,floor=0.5,spread=0.01)
p2=G.negative_control("纯随机分侧下两侧比值应≈1",
                      null=max(abs(np.nanmedian(ra[:,0])-1),abs(np.nanmedian(ra[:,1])-1)),
                      effect=max(abs(re_ratio-1),abs(rd_ratio-1)),null_spread=0.05,
                      null_kind="随机分侧下的同一个两侧比值 —— 若分侧本身不制造差异,它应当等于 1")
if p1 and p2:
    de=not (0.5<=re_ratio<=2); dd=not (0.5<=rd_ratio<=2)
    v=("**与方向无关:两个分量的两侧比值都在 [0.5,2] 内**" if (not de and not dd) else
       f"**是「更极端」那一件有方向:ext 比值 {re_ratio:.3f} 出界,而 disp {rd_ratio:.3f} 在界内**" if (de and not dd) else
       f"**是「更一致」那一件有方向:disp 比值 {rd_ratio:.3f} 出界,而 ext {re_ratio:.3f} 在界内**" if (dd and not de) else
       f"**两件都有方向:ext {re_ratio:.3f} / disp {rd_ratio:.3f} 都出界 —— 分侧本身可能不是一个好切法**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(j)),r_tot=rc(E,tot),q_tot=qt,sides=obs,ratio_ext=re_ratio,ratio_disp=rd_ratio,
               placeboA=[float(np.nanmedian(ra[:,0])),float(np.nanmedian(ra[:,1]))],
               placeboB=[float(np.nanmedian(rb[:,0])),float(np.nanmedian(rb[:,1]))],
               corr_stance_tot={"GSS":-0.2189,"MFQ":rc(stance,tot),"RWAS":0.7656},
               verdict=v,unchallenged=True),open(OUT/"which_component.json","w"),indent=1,ensure_ascii=False)
