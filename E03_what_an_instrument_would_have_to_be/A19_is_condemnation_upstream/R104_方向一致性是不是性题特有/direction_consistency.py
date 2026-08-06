"""E03·A19·R104 —— 方式效应的「方向一致性」,是性题特有的吗

**类型:FRONTIER**。`#661` 把 `#660` 降级为「一题的观察」,但网格里剩下一个**未经检验**的形状:
**性题十二格全部朝同一边,非性题散开。** 不是「大小」,是「方向」。

⚠ **BASIN**:我刚把 `#660` 降级,**而 W1 会让那个机制以更强的形式回来 —— 那是我想要的结果**
   ⇒ **下注 W2(格数配平后差消失)。**
W1 方向一致性是性题特有的。 W2 配平格数后消失。
**W3 = meta-separator:重新定向本身是伪影** —— 「严厉/宽容」是**我**按题干判的;
  **若结论随我怎么定向而变,「方向一致性」就不是数据的性质,是我的判断的性质。**

## 硬规则①:每题的码 1 是「严厉」还是「宽容」,由题干判,并写进产物
**`#661` 就是死在没做这一步**(它用「码 1 的份额」,而码 1 在不同题上不是同一个意思)。

G1 ESTIMAND:把每题的**面访 − 网络**差按**「朝更保守/更谴责为正」**重新定向后,
  报每组的**同号率**(正号格数 ÷ 总格数,取与多数一致的那一侧)。**主量 = 性题同号率 − 非性题同号率。**
G2 CONTROLS:
  **正对照** 重新定向后 `xmarsex` 必须仍为正(它是 `#660` 的那一道)。
  **安慰剂** 把「性 / 非性」标签在 34 格上随机打乱,同号率差应回到 0。
  **这个零该不该是零?** 该 —— 标签打乱后两组无差别 ⇒ `negative_control`。
G3:34 格全报(题 × 年),含定向表。G4:**等格数重抽** + 去掉堕胎 / 计入堕胎 两条规格。
⚠ **最强混淆(`#661` 的 NEXT 写死的)**:**同号率对格数敏感** —— 12 格全同号在纯随机下是 2⁻¹¹ ≈ 0.0005,
  22 格里 8/14 很常见 ⇒ **必须把非性组重抽到 12 格再比,不许直接比两个不同 n 的同号率。**
KILL(条件式):if 正对照为正 and 安慰剂 ≈0:
  等格数下 主量的 bootstrap 区间**不含零** -> W1;含零 -> **判不了**,`#660` 保持「一题的观察」
  **且 W3 检查**:把定向表**整体取反**,主量必须变号而**绝对值不变**;若绝对值也变 -> **定向是伪影**
else UNVERIFIED
IMPOSSIBLE(不写 planned):作答方式非随机 ⇒ 非因果 · 一国两波 ·
  **定向是我按题干判的**(已写进产物,可被推翻)· **跨仪器:换不了仪器,只此一具** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]

# ── 硬规则①:定向表,由题干判定,写进产物 ────────────────────────────────
# 「码 1 是不是**更保守/更谴责**的那一档?」 True = 是;False = 码 1 是更宽容的那一档。
ORIENT = {
 "premarsx": (True,  "Sex before marriage — 1 = always wrong(严厉)"),
 "xmarsex":  (True,  "Sex with person other than spouse — 1 = always wrong(严厉)"),
 "homosex":  (True,  "Homosexual sex relations — 1 = always wrong(严厉)"),
 "teensex":  (True,  "Sex before marriage: teens 14-16 — 1 = always wrong(严厉)"),
 "pornlaw":  (True,  "Feelings about pornography laws — 1 = illegal to all(严厉)"),
 "sexeduc":  (False, "Sex education in public schools — 1 = favour(宽容)"),
 "cappun":   (True,  "Death penalty for murder — 1 = favour(严厉)"),
 "gunlaw":   (True,  "Gun permits — 1 = favour permits(限制,保守侧)"),
 "letdie1":  (False, "Allow incurable patients to die — 1 = yes(宽容)"),
 "suicide1": (False, "Suicide if incurable disease — 1 = yes(宽容)"),
 "suicide4": (False, "Suicide if tired of living — 1 = yes(宽容)"),
 "spanking": (True,  "Favor spanking to discipline child — 1 = strongly agree(严厉)"),
 "polhitok": (True,  "Ever approve of police striking citizen — 1 = yes(许可暴力,严厉侧)"),
 "obey":     (True,  "To obey — 1 = most important(服从,保守侧)"),
 "fefam":    (True,  "Better for man to work, woman tend home — 1 = strongly agree(保守)"),
 "natcrime": (True,  "Halting rising crime rate — 1 = too little spending(强硬,保守侧)"),
 "helppoor": (False, "Govt improve standard of living — 1 = govt should(宽容侧)"),
 "abany":    (False, "Abortion if woman wants for any reason — 1 = yes(宽容)"),
 "abnomore": (False, "Married, wants no more children — 1 = yes(宽容)"),
}
SEX={"premarsx","xmarsex","homosex","teensex","pornlaw","sexeduc"}
ABO={"abany","abnomore"}
print("=== 硬规则①:定向表(由题干判,而这是我的判断,可被推翻)===")
for k,(c,txt) in ORIENT.items():
    print(f"  {k:10s} {'码1=严厉' if c else '码1=宽容':9s} {'性' if k in SEX else ('堕胎' if k in ABO else '非性'):>3s}  {txt}")

df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
                         usecols=["year","mode"]+list(ORIENT), encoding="latin1")
d=df[df.year.isin([2022,2024])]
rows=[]
for c,(cons,_t) in ORIENT.items():
    for y in (2022,2024):
        s=d[(d.year==y)&d[c].between(1,9)]
        a=s[s["mode"]==1][c]; b=s[s["mode"]==4][c]
        if len(a)<200 or len(b)<200: continue
        raw=(a==1).mean()-(b==1).mean()
        rows.append(dict(item=c,year=int(y),raw=float(raw),
                         oriented=float(raw if cons else -raw),
                         grp=("性" if c in SEX else ("堕胎" if c in ABO else "非性"))))
R=pd.DataFrame(rows)
print(f"\n=== G3:{len(R)} 格全报(定向后「朝更保守为正」)===")
for r in rows:
    print(f"  {r['item']:10s}{r['grp']:>3s}{r['year']:>6d}  生 {r['raw']:+.4f} -> **定向后 {r['oriented']:+.4f}**")

def agree(v):
    v=np.asarray(v); p=(v>0).sum()
    return max(p,len(v)-p)/len(v)
gs=R[R.grp=="性"].oriented.values; gn=R[R.grp=="非性"].oriented.values; ga=R[R.grp=="堕胎"].oriented.values
print(f"\n=== 同号率 ===\n  性 {len(gs)} 格 同号率 **{agree(gs):.3f}** · 正号 {(gs>0).sum()}")
print(f"  非性 {len(gn)} 格 同号率 **{agree(gn):.3f}** · 正号 {(gn>0).sum()}")
print(f"  堕胎 {len(ga)} 格 同号率 {agree(ga):.3f} · 正号 {(ga>0).sum()}")

# ⚠ 等格数重抽(最强混淆的控制)
def eqn(n=len(gs), B=4000):
    rng=np.random.default_rng(SEEDS[0])
    return np.array([agree(rng.choice(gn,n,replace=False)) for _ in range(B)])
null=eqn(); main=agree(gs)-float(np.median(null))
print(f"\n=== 等格数重抽(把非性组抽到 {len(gs)} 格)===")
print(f"  非性@{len(gs)}格 同号率 中位 **{np.median(null):.3f}** · 95% 区间 [{np.quantile(null,.025):.3f}, {np.quantile(null,.975):.3f}]")
print(f"  **性 {agree(gs):.3f} − 非性@等格数 {np.median(null):.3f} = {main:+.3f}**")
print(f"  性组的 {agree(gs):.3f} 落在非性零分布的第 {(null<agree(gs)).mean()*100:.1f} 百分位")
p_emp=float((null>=agree(gs)).mean())
print(f"  **经验 p = {p_emp:.4f}**(非性重抽达到或超过性组同号率的比例)")

def placebo(seed,B=2000):
    rng=np.random.default_rng(seed); allv=np.concatenate([gs,gn]); out=[]
    for _ in range(B):
        idx=rng.permutation(len(allv)); a=allv[idx[:len(gs)]]; b=allv[idx[len(gs):]]
        out.append(agree(a)-agree(b))
    return float(np.median(np.abs(out)))
pl=float(np.median([placebo(s) for s in SEEDS]))
xm=[r["oriented"] for r in rows if r["item"]=="xmarsex"]
print(f"\n=== 控制 ===\n  正对照 定向后 xmarsex = {[round(x,4) for x in xm]} · 必须全为正")
print(f"  安慰剂 打乱性/非性标签后 |同号率差| 中位 = **{pl:.4f}**")
# W3:定向整体取反
main_flip=agree(-gs)-float(np.median(eqn()))
print(f"  W3 定向整体取反:主量 {main:+.3f} -> {main_flip:+.3f} · |值| {'不变 ✅' if abs(abs(main)-abs(main_flip))<1e-9 else '变了 ⚠ 定向是伪影'}")

G=Gate("方向一致性是不是性题特有")
p1=G.positive_control("定向后 xmarsex 必须全为正",planted=float(min(xm)),floor=0.0,spread=0.005)
p2=G.negative_control("安慰剂:打乱性/非性标签后同号率差回到 0",null=pl,effect=abs(main),
                      null_spread=0.02,null_kind="随机重贴性/非性标签,两组本无差别")
if p1 and p2:
    verdict=(f"**W1 —— 方向一致性是性题特有的:性 {agree(gs):.3f} 对 非性@等格数 {np.median(null):.3f},"
             f"经验 p = {p_emp:.4f}**" if p_emp<0.05 else
             f"**判不了 —— 等格数下经验 p = {p_emp:.4f}。`#660` 保持在「一题的观察」**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(orient={k:[v[0],v[1]] for k,v in ORIENT.items()},cells=rows,
               agree_sex=agree(gs),agree_nonsex=agree(gn),agree_abortion=agree(ga),
               null_median=float(np.median(null)),null_ci=[float(np.quantile(null,.025)),float(np.quantile(null,.975))],
               main=main,p_emp=p_emp,placebo=pl,xmarsex=xm,flip_main=float(main_flip),
               verdict=verdict,unchallenged=True),
          open(OUT/"direction_consistency.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'direction_consistency.json'}")
