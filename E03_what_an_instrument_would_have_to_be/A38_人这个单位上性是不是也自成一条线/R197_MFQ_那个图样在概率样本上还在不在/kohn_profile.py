"""#756 · E03·A38·R197 —— MFQ 那个图样,在概率样本上还在不在?

`#752` 的人层结论(性独立于伤害-公平,与权威-内群体强相关)是在 **MFQ 这个自选样本**
(YourMorals.org)上拿到的 —— 那是它写在账本里的最大局限。
GSS 的 Kohn 育儿价值电池是**概率样本**,与四题性两两重叠 **n≈16,000**,
而且给出一个粗糙但对得上的基础类比:
  `obey` 服从 ≈ 权威 · `helpoth` 助人 ≈ 关怀/伤害 · `workhard` 勤奋 ≈ 成就 · `popular` 受欢迎 ≈ 内群体 · `thnkself` 独立思考 ≈ 自主
⇒ **从 MFQ 迁移过来的、可证伪的预测:`obey` 与性的严厉最强,`helpoth` 近零或反向。**

⚠⚠ **算术陷阱,写在最前面,因为它决定了什么能读什么不能读。**
这个电池是**严格 ipsative**:每人把五项排 1–5,**每人五题之和恒定(实测 std = 0.0000)**。
⇒ `Σᵢ cov(vᵢ, x) = cov(Σvᵢ, x) = cov(const, x) = 0` 对**任何**外部 x 成立。
**所以「一个正、一个负」这件事本身是被约束逼出来的,不是证据。**
可读的是**排名**(约束不决定谁最高);**不可读的是符号本身**。
⇒ 本轮把「和为零」当作**正控**:若实测之和不≈0,说明我对这具仪器的理解是错的。
⚠⚠ **第一版正控没过(Σ = −0.083…−0.113),而查因发现是我把它指向了错的统计量。**
   恒等式管的是**原始值协方差**;我算的是 Spearman —— **逐题秩化会摧毁 ipsative**
   (实测:秩化后五题之和的 std 从 **0.0000000000** 变成 **2620.48**),再除以各自 sd 又加了不同权重。
   实测 Σcov(原始值) = **−1.96e-15**,精确为零。
   ⇒ **控制是有效的,是我指错了对象** —— 与 `#750` 同一族,而这次是控制自己。
   ⇒ 正控改到协方差上;三种统计量(cov · Pearson · Spearman)一起报,当规格曲线读(G4)。

⚠ 编码方向(值标签实读,符号族第五次,这次写在跑之前):
  性题   `1=always wrong → 4=not wrong at all`  ⇒ 高 = **宽容**
  育儿题 `1=most important → 5=least important` ⇒ 高 = **不看重**
  ⇒ 「越不看重服从的人越宽容」在**原始编码下是正相关**。两边都是「越高越松」,方向一致。

预注册判词(阈值写在跑之前,且**分支覆盖整个区间**,`#754` 那一族):
  W1 图样复现:四题性中**至少三题**上 `obey` 是五个里绝对值最大的,且 `helpoth` 的绝对值排在后两位
  W2 图样不复现:`obey` 在**至多一题**上最大,或 `helpoth` 在**至少三题**上进前两位
  W3 其余一切情况 ⇒ 判不了,并**打印落在哪一段**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp
RNG=np.random.default_rng(197)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
K=["obey","thnkself","workhard","helpoth","popular"]; SEX=["premarsx","xmarsex","homosex","teensex"]
d=pd.read_stata(ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta",
                columns=["year"]+K+SEX,convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(lambda x:x>0) for c in K+SEX})
M["homosex"]=M["homosex"].where(M["homosex"]<=4); M["year"]=d.year
print("=== 硬规则①:变量名不是测量 ===")
for c in K+SEX:
    yrs=sorted(M.loc[M[c].notna(),"year"].unique().astype(int))
    print(f"  {c:9s} n={int(M[c].notna().sum()):6d} 年份 {len(yrs)} 个 {yrs[0]}–{yrs[-1]}")

def prof(sub,s,kind):
    if kind=="cov":  return {k:float(np.cov(sub[k],sub[s])[0,1]) for k in K}
    if kind=="pear": return {k:float(np.corrcoef(sub[k],sub[s])[0,1]) for k in K}
    return {k:sp(sub[k],sub[s]) for k in K}
PROF={}; rows={}; sums={}
for s in SEX:
    sub=M[K+[s,"year"]].dropna()
    PROF[s]={kind:prof(sub,s,kind) for kind in ("cov","pear","spear")}
    rows[s]=dict(PROF[s]["cov"]); rows[s]["_n"]=len(sub)
    sums[s]=sum(PROF[s]["cov"][k] for k in K)
print(f"\n=== 正控:五个**原始值协方差**之和必须 ≈ 0(ipsative 的算术后果)===")
for s in SEX: print(f"  {s:9s} Σcov = {sums[s]:+.3e}   {'通过' if abs(sums[s])<1e-9 else '⚠ 不通过 —— 我对这具仪器的理解是错的'}")
print(f"  ⚠ 对照:Σ Pearson 与 Σ Spearman **不**继承这个恒等式 ——")
for s in SEX:
    print(f"    {s:9s} ΣPearson {sum(PROF[s]['pear'][k] for k in K):+.4f} · ΣSpearman {sum(PROF[s]['spear'][k] for k in K):+.4f}")
if max(abs(v) for v in sums.values())>=1e-9:
    print("\n**UNVERIFIED:正控没过,后面的排名不许读**"); sys.exit(0)

for kind,nm in (("cov","原始值协方差(恒等式管的就是它)"),("pear","Pearson"),("spear","Spearman")):
    print(f"\n=== 五个育儿价值 × 四题性 · {nm}(高=松,两边一致)===")
    print(f"  {'':10s}"+"".join(f"{s[:9]:>11s}" for s in SEX))
    for k in K: print(f"  {k:10s}"+"".join(f"{PROF[s][kind][k]:+11.4f}" for s in SEX))
print(f"  {'n':10s}"+"".join(f"{rows[s]['_n']:11d}" for s in SEX))

# 每题性上,五个价值按绝对值排名 —— ⚠ G4:三种统计量各排一次,不一致就是发现
rank={s:sorted(K,key=lambda k:-abs(PROF[s]["cov"][k])) for s in SEX}
agree=all(sorted(K,key=lambda k:-abs(PROF[s][kd][k]))[0]==rank[s][0] for s in SEX for kd in ("pear","spear"))
print(f"\n  ⚠ G4 规格曲线:三种统计量在四题上是否都给出同一个第一名 -> {'是' if agree else '**否 —— 排名依赖统计量的选择,这本身是发现**'}")
print(f"\n  按|相关|排名:")
for s in SEX: print(f"    {s:10s} {' > '.join(rank[s])}")
top_obey=sum(1 for s in SEX if rank[s][0]=="obey")
help_top2=sum(1 for s in SEX if "helpoth" in rank[s][:2])
help_bot2=sum(1 for s in SEX if "helpoth" in rank[s][-2:])
print(f"\n  obey 排第一的题数 {top_obey}/4 · helpoth 进前二的题数 {help_top2}/4 · helpoth 落后二的题数 {help_bot2}/4")

# 零:年内打乱性题(保住 ipsative 结构与年代构成)
NP=400
nul={k:[] for k in K}
for _ in range(NP):
    sub=M[K+["premarsx","year"]].dropna().copy()
    sub["p"]=sub.groupby("year")["premarsx"].transform(lambda v: RNG.permutation(v.to_numpy()))
    for k in K: nul[k].append(abs(sp(sub[k],sub.p)))
print(f"\n=== 零(年内打乱 premarsx,{NP} 次,保住 ipsative 结构与年代构成)===")
for k in K: print(f"  {k:10s} 零的 95% 分位 {np.quantile(nul[k],.95):.4f} · 观测 |{abs(rows['premarsx'][k]):.4f}| ⇒ {abs(rows['premarsx'][k])/np.quantile(nul[k],.95):.1f}×")

print("\n"+"="*66)
if top_obey>=3 and help_bot2>=3:
    v=(f"**W1:`obey` 在 {top_obey}/4 题上是五个里最强的,`helpoth` 在 {help_bot2}/4 题上落在后两位 "
       f"⇒ MFQ 那个图样在**概率样本**上复现了**")
elif top_obey<=1 or help_top2>=3:
    v=(f"**W2:`obey` 只在 {top_obey}/4 题上最强,`helpoth` 在 {help_top2}/4 题上进前二 "
       f"⇒ MFQ 的图样很可能是自选样本的产物**")
else:
    v=(f"**W3:obey 第一 {top_obey}/4 · helpoth 后二 {help_bot2}/4 · helpoth 前二 {help_top2}/4 "
       f"—— 两条预注册都没到,判不了;落在 W1({'≥3 且 ≥3'}) 与 W2('≤1 或 ≥3')之间的那一段**")
print(v)
print("\n⚠ 无论判词是哪个:**符号本身不可读**(ipsative 强制 Σ=0),可读的只有排名。")
json.dump(dict(rows=rows,sums=sums,rank=rank,top_obey=top_obey,help_top2=help_top2,
               help_bot2=help_bot2,null95={k:float(np.quantile(nul[k],.95)) for k in K},verdict=v),
          open(OUT/"kohn_profile.json","w"),ensure_ascii=False,indent=1)
