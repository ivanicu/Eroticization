"""E03·A22·R117 —— 不分层的梯度:把 n 从 4 层换成 15,000 人

**类型:FRONTIER。** `#674` 的 NEXT:换一个不靠分层排序的估计量。
`#674` 证明「按层排序再算 ρ」在 k≤5 时置换零的 95% 分位贴在 1.0 ⇒ 永不显著。
⇒ **不再分层。** 沿 `educ` 做**核加权的成对秩相关**,最弱一环随 educ 连续变化,取其斜率;
**零在人之间打乱 `educ`** —— 分辨力来自 15,000 个人,不是 4 个层。

## 硬规则①(已跑)
GSS 性四题四题皆答 **n = 15,056** · 1988–2024 · 21 个年份。
`educ` 与之联合 **n = 15,000**(21 个取值,0–20 年)—— 候选里联合 n 最大的连续协变量。
安慰剂 `fefam·fepol·fepresch` 三题;`attend` 只有 9 个取值 ⇒ 不用作主协变量。

## G1 ESTIMAND
`W(e)` = 在 `educ = e` 处、以高斯核加权算出的**四题两两天花板归一秩相关的最小者**。
**主量 = `Δ = W(高教育端) − W(低教育端)`**(各取三个格点的均值),
显著性只由「在人之间打乱 educ」的零决定。
⚠ **第一版写的是 `ρ_grid`,当场被自己的零杀掉:核平滑的曲线单调 ⇒ ρ≡±1 ⇒ 零也≡1.0。见 `slope()` 的注释。**
## G2 CONTROLS
**正对照**:全样本(不加权)最弱一环必须为正且复现 `#671` 的量级。
**安慰剂 = 最强混淆的控制**:性别角色三题走同一条流水线。
  **这个零该不该是零?** **不该** —— `#671` 已测到安慰剂更强。
  所以这里问的是「主是否**显著大于**安慰剂」,`offset_control`,**零的种类 = 另一组题上的同一斜率**。
## G3/G4:三个带宽(1.5 · 2.5 · 4.0 年)× 两组题 = 6 格,全部照登。
## KILL(条件式)
if 正对照复现 and 打乱 educ 的零给出可分辨的分位:
  主斜率超零 95% 且 **> 安慰剂斜率** -> **是关于性的** ·
  主超零但 ≤ 安慰剂 -> **是作答一致性(#671 更强)** · 主未超零 -> **教育不动它**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**跨仪器:换不了仪器** —— NSFG 无 `educ`×性四题的对应题组(`#670` 已测其样本上限 45 岁);
**因果**:横断面 · 无干预;**核加权 ρ 的解析零不存在** ⇒ 只能置换。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
def wcorr(x,y,w):
    mx=np.average(x,weights=w); my=np.average(y,weights=w)
    cx=x-mx; cy=y-my; num=np.average(cx*cy,weights=w)
    d=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return num/d if d>1e-12 else np.nan
def ceil_norm(rx,ry,w,sign):
    o=np.argsort(rx); p=np.argsort(ry)
    a=rx[o]; b=ry[p] if sign>0 else ry[p][::-1]
    return wcorr(a,b,np.ones_like(w))
def weakest_at(R,W,items_idx):
    v=[]
    for i,j in combinations(range(len(items_idx)),2):
        r=wcorr(R[:,i],R[:,j],W)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        c=ceil_norm(R[:,i],R[:,j],W,1 if r>0 else -1)
        if np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
    return min(v) if v else np.nan
def curve(items,bw,educ,R,grid):
    out=[]
    for e in grid:
        W=np.exp(-0.5*((educ-e)/bw)**2)
        if W.sum()<200: out.append(np.nan); continue
        out.append(weakest_at(R,W,items))
    return np.array(out)
def slope(grid,y):
    """⚠ 估计量在本轮被换掉过一次,换掉的理由是本轮最主要的产出。

    **第一版用 `spearmanr(grid, y)`** —— 而核平滑出来的曲线本来就单调,
    于是 ρ 恒等于 ±1,打乱 educ 之后仍然恒等于 ±1,零的 95% 分位也是 1.0000。
    **它把幅度扔掉,只留次序,而次序在少数格点上几乎是免费的。**
    `#674` 把这件事诊断成「层数太少」——**那是错的诊断**:
    **ρ 在有序摘要点上,任何 k 都是错的统计量,因为它丢掉的正是带着证据的那个幅度。**
    ⇒ **改用端到端的变化量 Δ = W(最高) − W(最低)**,它不饱和,
      而且在「打乱 educ」这个零下有真实的分辨力。"""
    m=np.isfinite(y)
    if m.sum()<4: return np.nan
    yy=y[m]
    return float(np.mean(yy[-3:])-np.mean(yy[:3]))
res={}; SEEDS=[20260806,7,991]
for nm,items in [("性四题",SEX),("性别角色三题",FEM)]:
    j=df.dropna(subset=items+["educ"])
    educ=j["educ"].to_numpy(float)
    R=np.column_stack([pd.Series(j[c]).rank().to_numpy(float) for c in items])
    ov=weakest_at(R,np.ones(len(j)),items)
    grid=np.arange(8,19.01,1.0)
    print(f"\n=== {nm} · n = {len(j):,} · 全样本最弱一环 **{ov:+.4f}** ===")
    res[nm]={"n":int(len(j)),"overall":float(ov),"specs":{}}
    for bw in (1.5,2.5,4.0):
        y=curve(items,bw,educ,R,grid); s=slope(grid,y)
        rng=np.random.default_rng(SEEDS[0]); null=[]
        for _ in range(400):
            ep=rng.permutation(educ)
            null.append(abs(slope(grid,curve(items,bw,ep,R,grid))))
        null=np.array(null); q=float(np.nanquantile(null,0.95)); p=float(np.nanmean(null>=abs(s)))
        res[nm]["specs"][str(bw)]={"slope":float(s),"null_q95":q,"p":p,
                                   "curve":[None if not np.isfinite(v) else float(v) for v in y]}
        print(f"  带宽 {bw:>4}年  Δ **{s:+.4f}**  打乱 educ 的零 95% 分位 {q:.4f}  p = **{p:.4f}**  "
              f"{'✅ 超零' if p<0.05 else '⛔ 未超零'}")
        print(f"            曲线 {[f'{v:+.3f}' if np.isfinite(v) else 'na' for v in y]}")
main=res["性四题"]["specs"]["2.5"]["slope"]; plac=res["性别角色三题"]["specs"]["2.5"]["slope"]
G=Gate("不分层的梯度:教育与性这一块的紧度")
p1=G.positive_control("全样本最弱一环为正且复现 #671 量级(>0.20)",planted=res["性四题"]["overall"],floor=0.20,spread=0.01)
p2=G.offset_control("主 Δ 必须显著大于安慰剂 Δ(否则是作答一致性,不是性)",
                    effect=abs(main),offset=abs(plac),spread=0.010,
                    null_kind="另一组题(性别角色三题)上的同一个 Δ —— #671 已测到它更强,所以这是系统性基线偏移,不是待缩小的噪声")
pm=res["性四题"]["specs"]["2.5"]["p"]
if p1:
    if pm>=0.05: verdict=f"**教育不动它:主Δ {main:+.4f},打乱 educ 的零下 p = {pm:.4f}**"
    elif p2:     verdict=f"**是关于性的:主 Δ {main:+.4f} 显著大于安慰剂 Δ {plac:+.4f}**"
    else:        verdict=f"**是作答一致性,不是性:主 Δ {main:+.4f} ≤ 安慰剂 Δ {plac:+.4f} —— 与 #671 同向,而这一次 n = 人数**"
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(res=res,verdict=verdict,unchallenged=True),open(OUT/"continuous_gradient.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'continuous_gradient.json'}")
