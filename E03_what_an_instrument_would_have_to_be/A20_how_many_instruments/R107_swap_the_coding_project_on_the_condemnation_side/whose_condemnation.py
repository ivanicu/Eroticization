"""E03·A20·R107 —— 那条 +0.80,换掉谴责那一侧的编码项目还在不在

**类型:FRONTIER**。`#664` 把「越稀有越被谴责(婚前性)」从「做不到」改判「待做」。
兑现它时,硬规则①把这一轮变成了一个更锋利的东西。

## 硬规则①做完之后,这一轮不是原来那一轮

**`EA078` 已经在 SCCS 里** —— 它就是 `SCCS282`(来源栏写着 Murdock (1962-1971); Gray (1999),
标题里直接注明 *identical to EA078*)。**所以「第二具仪器」不是新的库,而是同一批社会上的另一个编码项目。**
而 SCCS 里**同一个构念被四个项目各编了一次**:
`SCCS165` Broude 1976 · `SCCS282` Murdock/Gray · `SCCS961` Frayser 1985 · `SCCS596` Whyte 1978。

**⚠ 而页面那条 +0.80 是 `SCCS165 × SCCS167` —— 谴责与频率都出自 Broude 同一个项目。**
`#528` 已经量过:**同队变量对的相关强度是跨队的 2.14 倍。**
⇒ **真正该问的是:把谴责那一侧换成另一个项目,+0.80 还在不在。**

⚠ **BASIN**:+0.80 是这一页最老、最强的声明之一 ⇒ **下注 W2(它是同项目伪影)。**
W1 换项目后仍强 ⇒ 构念稳健。 W2 换项目后塌 ⇒ **+0.80 是同项目伪影,页面必须改写。**
**W3 = meta-separator:若两个项目的编码几乎相同(ρ≈1),它们是抄本而不是两次测量,
  那么「跨编码项目复制」这个类别本身在民族志数据上就不成立,而 A20 整条线要重估。**

## 硬规则①救了这一轮两次(码,不是变量名)
- **`SCCS165`** n=130:`1` Expected → `6` Strongly disapproved,**越大越谴责**。
- **`SCCS167`** n=109:`1` Universal → `4` Uncommon,**越大越稀有**。
- **`SCCS282`** n=146:`2` 禁止严惩 · `3` 禁止轻罚 · `4` 允许但怀孕则罚 · `6` 允许无制裁 ——
  **越大越宽容,极性与 `SCCS165` 相反**;而 `1`(早婚排除)与 `5`(试婚)**不在这条尺上,剔除**
  (与 `#529` 丢 `SCCS172` 同一条规矩)。
- **⚠ `SCCS961`** n=61:码 `7` = *strongly disapproved **and rare*** ——
  **码本身把谴责与稀有捆在一起**,与稀有度相关是**定义上的必然**(realstat 的「算术陷阱」)
  ⇒ **主检验剔除,单列并标注。**

## G1 ESTIMAND(先于方法)
稀有度固定为 **`SCCS167`(Broude)**;谴责侧取 **`SCCS165`(Broude,同项目)** 与
**`SCCS282`(Murdock,跨项目,已定向)**。**主量 = 跨项目那一个的 ρ。**
## G2 CONTROLS
**正对照**:同项目 `SCCS167 × SCCS165` 必须复现 `#528` 的 **+0.80**(容差 0.10)。
**安慰剂**:`SCCS282`(定向后)× 与性无关的 `SCCS31`(社区规模类)必须 ≈0。
  **这个零该不该是零?** 该 ⇒ `negative_control`。
**W3 检查**:`ρ(SCCS165, SCCS282定向)` —— **≈1 ⇒ 抄本,判不了**。
## G3/G4:三个谴责编码 × {原始, 定向} 全报,含被剔除的 Frayser(标注定义污染)。
## KILL(条件式)
if 正对照复现 and 安慰剂≈0 and W3 检查 ρ(165,282) < 0.90:
  跨项目 |ρ| **≥ 0.40** -> W1 构念稳健 · **< 0.40 且同项目 ≥0.70** -> **W2:+0.80 是同项目伪影**
  两者之间 -> 判不了
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**四个项目读的是同一批民族志** ⇒ 独立性只是编码独立,不是观察独立 · 无干预 ·
**跨库:EA078 已在 SCCS 内,故「换库」在这个构念上等价于「换项目」** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
B="data/external/dplace/repo/datasets/SCCS/"
D=pd.read_csv(B+"data.csv"); SOC=pd.read_csv(B+"societies.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BLK={r.id:(int(np.floor(r.Lat/10)),int(np.floor(r.Long/10))) for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}

RAR="SCCS167"
COND={"Broude 1976(同项目)":("SCCS165",+1,None),
      "Murdock/Gray(跨项目,=EA078)":("SCCS282",-1,{1.0,5.0}),
      "⚠Frayser 1985(码 7 含 rare,定义污染)":("SCCS961",+1,None)}
def series(v,sign,drop):
    s=W[v].dropna()
    if drop: s=s[~s.isin(drop)]
    return sign*s
print("=== 硬规则①:各列实际可用 n 与联合 n ===")
print(f"  {RAR} 稀有度 n={int(W[RAR].notna().sum())}")
rows=[]
for name,(v,sign,drop) in COND.items():
    s=series(v,sign,drop); m=pd.concat([W[RAR],s.rename("c")],axis=1).dropna()
    print(f"  {name:34s} {v} 可用 n={len(s):4d} · 与稀有度联合 n=**{len(m)}**")
for name,(v,sign,drop) in COND.items():
    s=series(v,sign,drop); m=pd.concat([W[RAR].rename("r"),s.rename("c")],axis=1).dropna()
    if len(m)<30: rows.append(dict(name=name,var=v,n=len(m),rho=np.nan,note="判不了(n<30)")); continue
    r=float(spearmanr(m.r,m.c).statistic)
    def f(d): return float(spearmanr(d.r,d.c).statistic)
    bl=sorted({BLK.get(x,("na","na")) for x in m.index}); by={b:[x for x in m.index if BLK.get(x,("na","na"))==b] for b in bl}
    bs=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(200):
            socs=[x for i in rng.integers(0,len(bl),len(bl)) for x in by[bl[i]]]
            v_=f(m.loc[socs])
            if np.isfinite(v_): bs.append(v_)
    lo,hi=np.quantile(bs,[.025,.975])
    rows.append(dict(name=name,var=v,n=int(len(m)),rho=r,lo=float(lo),hi=float(hi)))
print("\n=== G3:稀有度(Broude SCCS167)× 三个谴责编码 ===")
for x in rows:
    if np.isnan(x["rho"]): print(f"  {x['name']:34s} n={x['n']:3d}  {x['note']}"); continue
    print(f"  {x['name']:34s} n={x['n']:3d}  ρ = **{x['rho']:+.4f}**  95%CI [{x['lo']:+.4f}, {x['hi']:+.4f}]")
same=[x for x in rows if "同项目" in x["name"]][0]
cross=[x for x in rows if "跨项目" in x["name"]][0]
print(f"\n  **同项目 {same['rho']:+.4f} -> 跨项目 {cross['rho']:+.4f} · 掉了 {abs(same['rho'])-abs(cross['rho']):+.4f}**")

# W3:两个谴责编码彼此有多像
m2=pd.concat([W["SCCS165"].rename("a"),series("SCCS282",-1,{1.0,5.0}).rename("b")],axis=1).dropna()
w3=float(spearmanr(m2.a,m2.b).statistic)
print(f"\n=== W3(meta-separator):两个项目的谴责编码彼此 ρ = **{w3:+.4f}**(n={len(m2)})")
print(f"  {'⚠ ≈1 ⇒ 抄本,判不了' if abs(w3)>=0.90 else '**< 0.90 ⇒ 是两次真正不同的编码,不是抄本**'}")
# 安慰剂
alt=[c for c in ["SCCS31","SCCS63","SCCS64","SCCS1"] if c in W.columns]
pl=np.nan
for c in alt:
    mm=pd.concat([series("SCCS282",-1,{1.0,5.0}).rename("a"),W[c].rename("b")],axis=1).dropna()
    if len(mm)>=40: pl=abs(float(spearmanr(mm.a,mm.b).statistic)); print(f"  安慰剂 SCCS282定向 × {c} (n={len(mm)}) = {pl:.4f}"); break

G=Gate("那条 +0.80,换掉谴责那一侧的编码项目还在不在")
p1=G.positive_control("同项目 SCCS167×SCCS165 必须复现 #528 的 +0.80(容差 0.10)",
                      planted=float(0.10-abs(abs(same["rho"])-0.80)),floor=0.0,spread=0.005)
p2=G.negative_control("安慰剂:跨项目谴责 × 与性无关的社会变量",null=float(pl) if np.isfinite(pl) else 1.0,
                      effect=abs(cross["rho"]),null_spread=0.05,null_kind="与性道德无关的社会结构变量")
if p1 and p2 and abs(w3)<0.90:
    if abs(cross["rho"])>=0.40: verdict=f"**W1 —— 构念稳健:跨项目仍 {cross['rho']:+.4f}**"
    elif abs(same["rho"])>=0.70: verdict=(f"**W2 —— +0.80 是同项目伪影:换一个编码项目,{same['rho']:+.4f} -> {cross['rho']:+.4f}**")
    else: verdict="**判不了 —— 落在两条判据之间**"
elif abs(w3)>=0.90: verdict=f"**判不了 —— 两个编码 ρ={w3:+.4f} ≈ 抄本(W3)**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(rows=rows,w3=w3,placebo=float(pl) if np.isfinite(pl) else None,verdict=verdict,
               note="SCCS282 = EA078;SCCS961 码7 含 rare,定义污染,单列不入主检验",unchallenged=True),
          open(OUT/"whose_condemnation.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'whose_condemnation.json'}")
