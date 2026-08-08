"""#768 第三臂 —— 「去年 vs 一生」这个对比,在第二具仪器上还成不成立?

GSS 臂量到:本人「去年」的性经历只解释 **1.9%**,而「一生伴侣数」解释 **9.1%** —— **差 4.8 倍**。
⚠ 那是**一个仪器上的一次对比**。而本轮真正的方法论断言是
**「同一个名字下的两个时间窗,是两个构念」** —— 这句话必须换仪器。

NSFG 2011–13 女性卷:不同问卷、不同编码团队、不同人群(15–44 岁女性),
**同时有 `PARTS1YR`(去年伴侣数)与 `LIFPRTNR`(一生伴侣数)**,
以及**对别人性行为的态度**:`SXOK18`/`SXOK16`(未婚 18/16 岁发生性行为可不可以)· `GAYADOPT` · `OKCOHAB`。

⚠ **它没有 `obey`** ⇒ **不能复现「偏掉 obey 那条线」的份额**;
   能复现的是本轮的**核心测量事实**:同一个态度题,与「一生」的关联比与「去年」的强多少。
   **这一格如实标注为「部分复现」,不是完全复现。**

G1 估计量:`|ρ(态度, LIFPRTNR)| − |ρ(态度, PARTS1YR)|`,逐题算,报整张网格(G3)。
预注册:**若每一题上「一生」都比「去年」强 ⇒ 方法论断言在第二具仪器上重现;
若符号不一 ⇒ 那是 GSS 那一次的特殊性,第一/二臂的推广要收窄。**
"""
import pandas as pd, numpy as np, re, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
RNG=np.random.default_rng(2207)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
NS=ROOT/"data/external/nsfg"
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
ATT=["sxok18","sxok16","gayadopt","okcohab"]; BIO=["lifprtnr","parts1yr"]; EX=["ager"]
need=[c for c in ATT+BIO+EX if c in LAY]
print("=== 硬规则①:变量名不是测量 ===")
for c in need: print(f"  {c:10s} 列位 {LAY[c][0]+1} 宽 {LAY[c][1]}  {LAY[c][2][:60]}")
buf={n:[] for n in need}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n in need:
        s,w,_=LAY[n]; v=line[s:s+w].strip()
        buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame(buf)
print(f"\n原始行数 {len(X)}")
for c in ATT:
    print(f"  {c:10s} 取值 {sorted(pd.Series(X[c].dropna().unique()).astype(int))[:12]}")
for c in BIO:
    v=X[c]; print(f"  {c:10s} 非缺失 {int(v.notna().sum())} · 分位 "
                  f"{[float(np.nanquantile(v,q)) for q in (.5,.9,.99)]} · 最大 {float(np.nanmax(v))}")
# ⚠ NSFG 的态度题 1–5 有效,其余是 DK/RF 的大码;伴侣数同理有大码哨兵
for c in ATT: X[c]=X[c].where(X[c].between(1,5))
for c in BIO: X[c]=X[c].where(X[c].between(0,200))
print(f"  ⚠ 态度题限 1–5、伴侣数限 0–200(大码是 DK/RF 哨兵)")
sp=lambda a,b: float(pd.Series(a).corr(pd.Series(b),method="spearman"))
print(f"\n=== G3 全网格:同一态度题,与「一生」比与「去年」强多少 ===")
print(f"  {'态度题':12s}{'n':>7s}{'ρ 一生':>10s}{'ρ 去年':>10s}{'|Δ|':>9s}{'倍数':>8s}")
rows={}
for a in ATT:
    sub=X[[a]+BIO].dropna()
    rl,rp=sp(sub[a],sub.lifprtnr),sp(sub[a],sub.parts1yr)
    rows[a]=dict(n=len(sub),rho_life=rl,rho_year=rp,d=abs(rl)-abs(rp),
                 ratio=(abs(rl)/abs(rp) if abs(rp)>1e-6 else None))
    print(f"  {a:12s}{len(sub):7d}{rl:+10.4f}{rp:+10.4f}{abs(rl)-abs(rp):+9.4f}"
          +(f"{abs(rl)/abs(rp):8.2f}×" if abs(rp)>1e-6 else f"{'-':>8s}"))
up=sum(1 for r in rows.values() if r["d"]>0)
# ---- 闸 ----
sub=X[["sxok18"]+BIO].dropna()
G=Gate("#768 第三臂 · NSFG 去年 vs 一生")
G.identity_control("① 同一列与自己的相关必须是 1(仪器活着吗)",
                   observed=sp(sub.lifprtnr,sub.lifprtnr), expected=1.0, tol=1e-9, what="秩相关的恒等检查")
nul=[abs(sp(sub.sxok18,RNG.permutation(sub.lifprtnr.to_numpy()))) for _ in range(300)]
G.negative_control("② 打乱一生伴侣数与人的配对后,关联须回到零",
                   null=float(np.median(nul)), effect=abs(sp(sub.sxok18,sub.lifprtnr)),
                   null_spread=float(np.std(nul)),
                   null_kind="在人之间打乱 `LIFPRTNR`,保住它的边际,只毁掉「谁的性史配谁的态度」")
print(); print(G)
print("\n"+"="*70)
if not all(r[2] for r in G.rows): v="**UNVERIFIED:闸没全过**"
elif up==len(rows):
    v=(f"**方法论断言重现:NSFG 上 {up}/{len(rows)} 题「一生」都强过「去年」"
       f"(倍数 {min(r['ratio'] for r in rows.values() if r['ratio']):.2f}×–{max(r['ratio'] for r in rows.values() if r['ratio']):.2f}×)"
       f" ⇒ 「同一个名字下的两个时间窗是两个构念」不是 GSS 那一次的特殊性**")
elif up==0: v=f"**反向:0/{len(rows)} 题 ⇒ 推广收窄到 GSS**"
else: v=f"**符号不一:{up}/{len(rows)} 题「一生」更强 ⇒ 不是普遍规律,整张网格已全列**"
print(v)
print("⚠ 本臂**不能**复现「偏掉 obey 那条线」的份额 —— NSFG 没有 `obey`;这是**部分复现**。")
json.dump(dict(rows=rows,n_up=up,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"nsfg_arm.json","w"),ensure_ascii=False,indent=1)
