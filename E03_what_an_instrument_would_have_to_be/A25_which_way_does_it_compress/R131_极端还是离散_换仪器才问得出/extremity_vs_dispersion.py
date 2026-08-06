"""E03·A25·R131 —— 压缩的是「更少极端」还是「更一致」——而这个问题只有换仪器才问得出

**类型:FRONTIER。** `#688` 的 NEXT。
**心理学的那一句(本轮要分开的两件):一个人答得更平,可以是立场更笃定(极端且一致),
也可以是更少极端(全部往中间收)。两者在 `totvar` 上一模一样,而在情欲化算子里是相反的。**

## ⚠ 硬规则①当场杀掉预注册的分解,而这是一条结构性登记,不是失败
GSS 的 33 题里 **27 题是二值** ⇒ **极端度在二值题上无定义**(每个答案都在端点);
只有 **6 题四档**(`fefam · fepresch · homosex · premarsx · teensex · xmarsex`)。
**⇒ 「极端度 vs 离散度」的分解在 GSS 上结构性做不到。** 按 realstat,**不许写「planned」。**
(顺带量到:在那 6 题上 `corr(极端度, totvar) = +0.2713`,|r| ≤ 0.7,两者确实不是同一个量 ——
**所以问题是真的,只是这具仪器答不了。**)

## ⇒ 按硬规则④改仪器(跨仪器复制优于同一仪器再来一轮)
**MFQ(Graham/Haidt/Nosek 2009 Study 3)的 30 道道德基础题全部是六档**,`educ_num` 七档 ——
**极端度在那里有定义。这不是绕路,是这份数据唯一允许的走法。**

## G1 ESTIMAND
`ext` = 该人答在刻度两端(1 或 6)的题占比;`totvar` = 该人 30 题标准化答案的标准差;
`disp` = **给定极端度下的离散度** = `totvar` 对 `ext` 回归后的残差。
**主量:`ρ(educ_num, ext)` 与 `ρ(educ_num, disp)`,n = 人数。**
## ⑧ 判据(`#688` 在跑之前写死,不得改)
**两者同号 ⇒ 压缩是整体的;只有 `ext` 动 ⇒ 教育让人更少极端(完全不同、且更接近本项目的问题);
只有 `disp` 动 ⇒ 是一致性。**
## G2 CONTROLS
**正对照**:MFQ 上必须复现 `ρ(educ_num, totvar) < 0`(即 `#686` 的方向在第二具仪器上也在),
**否则本轮不是在测同一件事,当场停。**
**⑤ 最强混淆(`#688` 预注册)**:极端度与离散度在有界离散刻度上也必然相关 ⇒
**每题在人之间独立打乱**(边际保留、跨题结构全毁)得到 `corr(ext, totvar)` 的**必然部分**,
**在结论里减掉它,不许算成教育的作用。**
**零**:打乱 `educ_num`。
## KILL(条件式)
if 正对照复现 and 打乱 educ 的零确实为零: evaluate(判据⑧) else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**GSS 上这个分解结构性做不到**(27/33 二值)—— 已登记;
MFQ 是**单次横断面、无 cohort、教育只有七档** ⇒ **拿不到世代分层,也无法复制 `#683`**;
MFQ 与 GSS 题目内容不同 ⇒ **本轮是「同一问题在第二具仪器上」,不是「同一测量的复制」**;
因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
d,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
ITEMS=["harm","emotionally","weak","cruel","compassion","animal","kill",
       "fairly","unfairly","treated","justice","rights","rich",
       "loyalty","betray","yourgroup","lovecountry","family","team","history",
       "duties","traditions","respect","chaos","kidrespect","soldier","shutup",
       "disgusting","decency","desires","god","harmlessdg","unnatural","chastity"]
ITEMS=[c for c in ITEMS if c in d.columns and d[c].nunique()>=5]
print(f"MFQ 可用六档题 **{len(ITEMS)}** 道 · 各题取值数 {sorted(set(int(d[c].nunique()) for c in ITEMS))}")
j=d.dropna(subset=ITEMS+["educ_num"]); j=j[j.educ_num<=7].reset_index(drop=True)
print(f"皆答 + educ_num ⇒ **n = {len(j):,}** · educ_num 取值 {sorted(j.educ_num.unique().astype(int))}")
V=j[ITEMS]; lo=V.min().min(); hi=V.max().max()
ext=((V==lo)|(V==hi)).mean(1).to_numpy(float)
Z=pd.DataFrame({c:(j[c]-j[c].mean())/j[c].std() for c in ITEMS})
tot=Z.std(1).to_numpy(float); E=j.educ_num.to_numpy(float)
rc=lambda a,b:float(np.corrcoef(pd.Series(a).rank(),pd.Series(b).rank())[0,1])
r_et=rc(ext,tot)
rng=np.random.default_rng(20260806); sims=[]
for _ in range(50):
    Vp=pd.DataFrame({c:rng.permutation(V[c].to_numpy()) for c in ITEMS})
    Zp=pd.DataFrame({c:(Vp[c]-Vp[c].mean())/Vp[c].std() for c in ITEMS})
    sims.append(rc(((Vp==lo)|(Vp==hi)).mean(1).to_numpy(),Zp.std(1).to_numpy()))
forced=float(np.median(sims))
print(f"\n⑤ 合成基线:每题独立打乱 ⇒ **corr(ext, totvar) 必然部分 = {forced:+.4f}** · "
      f"实测 **{r_et:+.4f}** ⇒ 必然占 **{abs(forced)/abs(r_et)*100:.0f}%**")
A=np.column_stack([pd.Series(ext).rank().to_numpy(),np.ones(len(ext))])
tr=pd.Series(tot).rank().to_numpy()
disp=tr-A@np.linalg.lstsq(A,tr,rcond=None)[0]
r_tot=rc(E,tot); r_ext=rc(E,ext); r_disp=rc(E,disp)
nul_t=np.array([abs(rc(rng.permutation(E),tot)) for _ in range(300)]); q_t=float(np.quantile(nul_t,.95))
nul_e=np.array([abs(rc(rng.permutation(E),ext)) for _ in range(300)]); q_e=float(np.quantile(nul_e,.95))
nul_d=np.array([abs(rc(rng.permutation(E),disp)) for _ in range(300)]); q_d=float(np.quantile(nul_d,.95))
print(f"\n正对照:MFQ 上 ρ(educ_num, totvar) = **{r_tot:+.4f}** (零 95% {q_t:.4f}) "
      f"{'✅ 与 #686 同向且超零' if (r_tot<0 and abs(r_tot)>q_t) else '⛔ 不复现 —— 当场停'}")
print(f"\n=== 分解 ===")
print(f"  ρ(educ_num, **极端度**)   = **{r_ext:+.4f}** · 零 95% {q_e:.4f} {'✅' if abs(r_ext)>q_e else '⛔ 在零里'}")
print(f"  ρ(educ_num, **给定极端度下的离散度**) = **{r_disp:+.4f}** · 零 95% {q_d:.4f} {'✅' if abs(r_disp)>q_d else '⛔ 在零里'}")
G=Gate("压缩的是更少极端还是更一致")
p1=G.positive_control("MFQ 上必须复现 #686 的方向(ρ(educ,totvar)<0 且超零)",
                      planted=1.0 if (r_tot<0 and abs(r_tot)>q_t) else 0.0,floor=0.5,spread=0.01)
p2=G.negative_control("打乱 educ_num 后两条关系都应消失",null=max(q_e,q_d),
                      effect=max(abs(r_ext),abs(r_disp)),null_spread=0.005,
                      null_kind="人层打乱 educ_num —— 若教育与两个分量都无关,打乱后应无差别")
se,sd_=abs(r_ext)>q_e, abs(r_disp)>q_d
if p1 and p2:
    v=("**压缩是整体的:极端度与离散度都随教育动且同号**" if (se and sd_ and r_ext*r_disp>0) else
       f"**教育让人更少极端,而不是更一致:极端度 {r_ext:+.4f} 超零、离散度 {r_disp:+.4f} 在零里 ⇒ "
       f"这是完全不同、且更接近本项目问题的一句话**" if (se and not sd_) else
       f"**是一致性:离散度 {r_disp:+.4f} 超零而极端度 {r_ext:+.4f} 在零里**" if (sd_ and not se) else
       f"**两者反号,照登:极端度 {r_ext:+.4f} / 离散度 {r_disp:+.4f} —— 分解本身需要重想**")
else: v="UNVERIFIED"
print(f"\n{v}"); print(G)
json.dump(dict(n=int(len(j)),n_items=len(ITEMS),r_tot=r_tot,r_ext=r_ext,r_disp=r_disp,
               q_tot=q_t,q_ext=q_e,q_disp=q_d,corr_ext_tot=r_et,forced=forced,
               verdict=v,unchallenged=True),open(OUT/"extremity_vs_dispersion.json","w"),indent=1,ensure_ascii=False)
