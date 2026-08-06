"""E03·A34·R178 —— 对着信度读那个 0.47

**类型:FRONTIER。新弧 `A34`。**

**心理学的那一句(本轮要判的):一个人对同一道性道德题换个问法只答到 0.84 的一致。
那么她对「婚前」和「同性」的判断相关 0.47,到底是「差得远」还是「已经顶到测量的天花板」?**

## 缺口(`#731` 打开的那个)
`#731` 证明**比值是相对池的**(同一组题换陪衬域,0.93× ↔ 4.56×),`#734` 又证明在 NSFG 那种
题题相关的池里比值整个失效。⇒ **「一个人有一套性道德」不能再靠比值来说。**
⚠ **而 `#732` 已经把一把与池无关的尺子交到手上:重测信度。**
一个人对**同一个构念、换一套问法**只答到 **0.73–0.91** 的一致 ⇒
**跨题相关的上限不是 1.0,是 √(r_xx·r_yy)。**

## G1 ESTIMAND
**衰减校正后的跨题相关** `r_xy / √(r_xx · r_yy)`,以及**校正后的最弱一环**。
**这是与池无关的量** —— 分母来自这具仪器自己的重测,不来自我挑的陪衬题。
## IDENTIFICATION
`r_xx` 只在 **1994 年、368 个两版都答过的人** 身上可估 ⇒ **本轮的 scope 就到这里,不外推。**

## W1 / W2 / W3
| 世界 | 校正后最弱一环 | 读法 |
|---|---|---|
| **W1 就是一件事** | **≈ 1** | 四道题在测量误差之内是同一道题 —— 「一套性道德」是字面真,而「块」这个说法变得平凡 |
| **W2 相关但可分** | **0.5–0.8** | 是一套,但四个判断仍是四件事 |
| **W3 只是弱相关** | **< 0.5** | **即使测量完美,四题也只是松散相关 ——「一套性道德」是比喻,不是测量结论** |

⚠ **W3 是我不高兴的那个(削页上人层招牌);而 W1 也不是好消息** ——
它会把「一块」变成一句同义反复。**两端都不舒服,中间才是有内容的,而这正是设计成能三分的理由。**

## G2 CONTROLS
**④ 正对照**:重测系数必须复现 `#732` 的 **0.8372 / 0.8428 / 0.9106 / 0.7341**(容差 0.005)。
⚠ **且必须在 g=0 时失败**:把两版的配对**打乱**(A 的第 i 个人配 B 的随机一个人)后,
重测系数必须回到 0 —— **否则它测的不是同一个人的两次作答。**
**零** = `negative_control`,**零的种类 = 打乱两版之间的个体配对 ——
保住两版各自的边际与每个人的作答,只毁掉「哪两条记录是同一个人」。**
**SHAM**:对**跨构念**的一对做同样的校正(如 `premarsx` × `homosex1`)——
**若校正把任何一对都推到 1,那是校正的性质,不是这四题的性质。**
## G3:六对全报。G4:用 A 版 / B 版 / 两版平均的 `r_xy` 各算一次。
## ⑤ 停止条件(跑之前写死)
- **重测系数复现不到 0.005,或打乱配对后不回到零 ⇒ UNVERIFIED 并停。**
- **校正后最弱一环 ≥ 0.90 ⇒ 判 W1**(并写明「块」在此变成平凡);
  **0.50–0.90 ⇒ W2**;**< 0.50 ⇒ W3,页上人层那一行要加限定。**
- **任何一对校正后 > 1.05 ⇒ 该对记「判不了」**(衰减校正在 r_xx 被低估时会越界)。
## IMPOSSIBLE(不写 planned)
⚠ **重测是跨「问法」的,不是同一套题的重复施测** ⇒ 它把**信度**与**措辞差异**混在一起,
**因此它低估信度,而校正后的值是一个上界。这一条不许省略。**
`teensex1` 只在 1994 问过一次;368 人是**两版都答过**的子样本,不是随机子样;
**换不了仪器**:GSS 内部只有这一次重复施测(`#732` 已枚举)。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A=["premarsx","xmarsex","homosex","teensex"]; Bv=["premars1","xmarsex1","homosex1","teensex1"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+A+Bv,encoding="latin1")
ov=g.dropna(subset=A+Bv); g94=g[g.year==1994]
print(f"=== 硬规则①:1994 年 · 两版都答过的 **n = {len(ov)}** · 只用 1994 的原版 n = {len(g94.dropna(subset=A)):,} ===")
def sp(a,b):
    # ⚠ 这版 numpy 上 spearmanr 直接吃 pandas Series 会崩,必须显式转 float 数组
    return float(spearmanr(np.asarray(a,dtype=float),np.asarray(b,dtype=float)).statistic)
LED={"premarsx":0.8372,"xmarsex":0.8428,"homosex":0.9106,"teensex":0.7341}
rxx={}
print("\n=== ④ 正对照:重测系数必须复现 `#732`(容差 0.005)===")
for a,b in zip(A,Bv):
    m=ov[[a,b]].dropna(); r=sp(m[a],m[b]); rxx[a]=r
    print(f"  {a:10s} r_xx = **{r:+.4f}** · 账本 {LED[a]:+.4f} · 差 {abs(r-LED[a]):.4f} {'✅' if abs(r-LED[a])<=0.005 else '⛔'}")
maxd=max(abs(rxx[a]-LED[a]) for a in A)
rng=np.random.default_rng(20260806)
nul=[]
for _ in range(2000):
    p=rng.permutation(len(ov))
    nul.append(np.mean([abs(sp(ov[a].to_numpy(),ov[b].to_numpy()[p])) for a,b in zip(A,Bv)]))
q95=float(np.quantile(nul,0.95))
print(f"\n  零(打乱两版之间的个体配对,B=2000):**零的 95% 分位 {q95:.4f}** · 中位 {np.median(nul):.4f}")
print(f"  ⇒ 实测均值 {np.mean(list(rxx.values())):.4f} vs 零 {q95:.4f} · "
      f"{'✅ 配对是真的' if np.mean(list(rxx.values()))>q95 else '⛔ ⑤ 触发'}")
if maxd>0.005 or np.mean(list(rxx.values()))<=q95:
    print("\n⛔ 停"); json.dump(dict(stop="正对照或零不合格",rxx=rxx),open(OUT/"dis.json","w"),indent=1,ensure_ascii=False); sys.exit(0)
print(f"\n=== G3 六对全报:原始 → 衰减校正(1994 原版 n={len(g94.dropna(subset=A)):,})===")
sub=g94.dropna(subset=A)
raw={}; cor={}
for a,b in itertools.combinations(A,2):
    r=sp(sub[a],sub[b]); c=r/np.sqrt(rxx[a]*rxx[b])
    raw[(a,b)]=r; cor[(a,b)]=c
    flag=" ⚠ >1.05,判不了" if c>1.05 else ""
    print(f"  {a:10s} × {b:10s}  原始 {r:+.4f} → **校正 {c:+.4f}**{flag}")
ok={k:v for k,v in cor.items() if v<=1.05}
wl_raw=min(raw.values()); wl_cor=min(ok.values()) if ok else np.nan
print(f"\n  最弱一环:原始 **{wl_raw:+.4f}** → **校正 {wl_cor:+.4f}**(可用对 {len(ok)}/6)")
print("\n=== SHAM:跨版、跨构念的一对做同样校正(校正会不会把任何一对都推到 1)===")
for a,b in (("premarsx","homosex1"),("xmarsex","teensex1")):
    m=ov[[a,b]].dropna(); r=sp(m[a],m[b]); c=r/np.sqrt(rxx[a]*rxx[b.replace('1','')])
    print(f"  {a:10s} × {b:10s} 原始 {r:+.4f} → 校正 **{c:+.4f}**")
G=Gate("对着信度读那个 0.47")
p1=G.positive_control("重测系数必须复现 #732,且打乱配对后回到零",
    planted=float(np.mean(list(rxx.values()))-q95),floor=0.0,spread=0.01)
p2=G.negative_control("打乱两版之间的个体配对后重测系数应回到零",null=q95,
    effect=float(np.mean(list(rxx.values()))),null_spread=0.005,
    null_kind="打乱两版之间的个体配对 —— 保住两版各自的边际与每个人的作答,只毁掉「哪两条记录是同一个人」")
if not p1: v="**UNVERIFIED:正对照没过**"
elif wl_cor>=0.90: v=f"**W1:校正后最弱一环 {wl_cor:.4f} ≥ 0.90 —— 在测量误差之内,四道题几乎是同一道题**"
elif wl_cor>=0.50: v=f"**W2:校正后最弱一环 {wl_cor:.4f} —— 是一套,而四个判断仍是四件事**"
else: v=f"**W3:校正后最弱一环 {wl_cor:.4f} < 0.50 —— 即使测量完美,四题也只是松散相关,页上人层那一行要加限定**"
print(f"\n{v}")
print("⚠ 必报的上界声明:**重测是跨「问法」的,把信度与措辞差异混在一起 ⇒ 它低估信度 ⇒ 校正值是上界。**")
print(G)
json.dump(dict(n_retest=int(len(ov)),n_1994=int(len(sub)),rxx=rxx,null_q95=q95,
  raw={f"{a}×{b}":raw[(a,b)] for a,b in raw},corrected={f"{a}×{b}":cor[(a,b)] for a,b in cor},
  weakest_raw=wl_raw,weakest_corrected=float(wl_cor),verdict=v,
  upper_bound_note="重测跨问法 ⇒ 低估信度 ⇒ 校正值是上界",unchallenged=True),
  open(OUT/"dis.json","w"),indent=1,ensure_ascii=False)
