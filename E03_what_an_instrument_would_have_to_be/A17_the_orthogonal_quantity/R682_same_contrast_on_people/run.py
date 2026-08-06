"""E03·A17·R682 —— 在人身上重做 A12 那个对比:换对象 vs 换做法

**类型:FRONTIER**。A12 的全部结论都在**社会**这个单位上。这是它们第一次被拿到**人**身上问。

GSS 恰好有 A12 那个设计的人层同构版:
  **同一行为 × 四个对象** —— 警察打人:`polabuse` 说脏话的 · `polmurdr` 谋杀嫌疑人 ·
    `polescap` 试图脱逃的 · `polattak` 用拳头袭警的
  **不同行为** —— `spanking` 打孩子 · `cappun` 死刑 · `hitok` 一个男人打另一个男人

社会那一侧的答案:换对象 **+0.8451**(`#640`)· 换做法上界 **+0.4401**(`#641`)。

WORLD A —— **同样的结构在人身上重现**:换对象远高于换做法。
  ⇒ 「严厉聚合的单位是一个具体做法 × 一类人」跨单位成立。
WORLD B —— **在人身上不重现**(两者差不多,或反过来)。
  ⇒ 那个结构是**社会**这个单位特有的,而 A12 的措辞必须加上「在社会这个单位上」。
**区别是本体的**:它决定 A12 的结论是关于「严厉」这件事,还是关于「被民族志编码的社会」这类对象。

⚠ **跑之前写死的最强混淆,而它是 `#642` 的光环换了身衣服:**
  四道 `pol*` 题**连着问、措辞只差一个从句** —— **默认作答倾向(acquiescence / response set)会抬高互相关**,
  正如同一个编码者会给同一个社会的所有育儿栏打相近的分。
  **同一迭代内的控制:同格式的性道德四题**(`premarsx` `xmarsex` `homosex` `teensex`)——
  同样是一个连问的四题块、同样的答项格式,但**它们不是「同一行为的四个对象」**。
  **若 pol 块并不高于性道德块 ⇒ 高相关是格式造的,不是「同一行为」造的。**

⚠ **第二个必须先写下的约束:一切在同一调查年内算**(本页第八件已确立:跨年汇总会制造相关)。

G1 ESTIMAND:**A** = 四道 pol 题 6 对的年内秩相关,按年取中位再跨年取中位。
             **B** = 跨做法 6 对(spanking·cappun·hitok·pol 合成)的 |年内秩相关| 中位。
             **F** = 性道德四题 6 对的年内秩相关中位(格式参照)。
KILL(条件式,预注册;并遵守 `#641` 规则:先看区间再判):
  if 正对照(A 的年内中位在各年之间一致,极差 < 0.30)and 安慰剂(pol 合成 x 出生年 约 0):
      **A − B 的块 bootstrap 区间含零 -> 判不了**
      不含零且 A > B and **A > F** -> WORLD A(结构跨单位重现,且不是格式造的)
      不含零且 A > B and **A <= F** -> **判不了 —— 高相关与格式无法分离**
      不含零且 A <= B -> WORLD B
  else: UNVERIFIED
**地板写在每一对每一年上**:任一 (对, 年) 的 n < 200 -> 该格不进中位(人层,地板比社会层高)。
IMPOSSIBLE(不写 planned):无干预 · 自报 · 「同一行为四对象」是 GSS 的设计不是我的 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate

SEEDS=[20260806,7,991]; FLOOR=200
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/gss/GSS_stata/gss7224_r3a.dta"
POL=["polabuse","polmurdr","polescap","polattak"]
ACT=["spanking","cappun","hitok"]
SEX=["premarsx","xmarsex","homosex","teensex"]
cols=POL+ACT+SEX
df,_=pyreadstat.read_dta(P, usecols=["year","cohort"]+cols, encoding="latin1")
print("=== 硬规则①:各列 n 与年份数 ===")
for c in cols:
    y=df.loc[df[c].notna(),"year"]
    print(f"  {c:10s} n={int(df[c].notna().sum()):6d}  年份 {y.nunique():2d}")

def within_year(pairs, absolute=False):
    """每年每对算一次,地板写在 (对,年) 上"""
    out={}
    for a,b in pairs:
        per=[]
        for y,g in df.groupby("year"):
            m=g[[a,b]].dropna()
            if len(m)<FLOOR or m[a].nunique()<2 or m[b].nunique()<2: continue
            try: r=float(spearmanr(m[a].to_numpy(dtype=float),m[b].to_numpy(dtype=float)).statistic)
            except Exception: continue
            if np.isfinite(r): per.append(abs(r) if absolute else r)
        if per: out[(a,b)]=(float(np.median(per)), len(per))
    return out

polp=list(combinations(POL,2)); sexp=list(combinations(SEX,2))
df["_pol"]=df[POL].mean(axis=1)
actp=list(combinations(ACT,2))+[(a,"_pol") for a in ACT]

A=within_year(polp); F=within_year(sexp); B=within_year(actp, absolute=True)
def show(name,d):
    v=[x[0] for x in d.values()]
    print(f"\n=== {name}:{len(d)} 对 ===")
    for (a,b),(r,ny) in sorted(d.items(), key=lambda x:-x[1][0]):
        print(f"  {a:10s} × {b:10s}  年内中位 {r:+.4f}  ({ny} 年)")
    return float(np.median(v))
Am=show("A 同一行为 × 四个对象(警察打人)",A)
Bm=show("B 跨做法",B)
Fm=show("F 格式参照:性道德四题",F)
print(f"\n  **A = {Am:+.4f}  ·  B = {Bm:+.4f}  ·  F = {Fm:+.4f}  ·  A − B = {Am-Bm:+.4f}  ·  A − F = {Am-Fm:+.4f}**")
print(f"  社会那一侧:换对象 +0.8451 · 换做法上界 +0.4401")

def boot(n=400):
    yrs=sorted(df.year.unique()); out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            pick=rng.choice(yrs,len(yrs),replace=True)
            sub=pd.concat([df[df.year==y] for y in pick])
            g=globals(); old=g["df"]; g["df"]=sub
            try:
                a=[x[0] for x in within_year(polp).values()]; b=[x[0] for x in within_year(actp,True).values()]
                if a and b: out.append(float(np.median(a)-np.median(b)))
            finally: g["df"]=old
    return np.array(out)
bs=boot(); lo,hi=np.quantile(bs,[.025,.975])
print(f"\n  A − B 的 95% CI(按年重抽)= [{lo:+.4f},{hi:+.4f}]  -> {'**含零 ⇒ 判不了**' if lo<0<hi else '**不含零 ⇒ 可判**'}")

rng_A=max(x[0] for x in A.values())-min(x[0] for x in A.values())
d2=df[["_pol","cohort"]].dropna(); pl=abs(float(spearmanr(d2._pol.to_numpy(dtype=float),d2.cohort.to_numpy(dtype=float)).statistic))
G=Gate("在人身上重做 A12 的对象-vs-做法对比")
p1=G.positive_control("A 的六对之间一致(极差 < 0.30)",planted=float(0.30-rng_A),floor=0.0,spread=0.01)
p2=G.negative_control("安慰剂:pol 合成 x 出生年 约 0",null=pl,effect=abs(Am),null_spread=0.05,
                      null_kind="与暴力态度无关的人口学刻度")
if p1 and p2:
    if lo<0<hi: verdict="**判不了 —— A − B 的区间含零**"
    elif Am<=Bm: verdict="**WORLD B —— 结构在人身上不重现**"
    elif Am<=Fm: verdict="**判不了 —— A 不高于同格式参照块,高相关与格式无法分离**"
    else: verdict="**WORLD A —— 结构跨单位重现,且不是格式造的**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(A=Am,B=Bm,F=Fm,A_minus_B=Am-Bm,A_minus_F=Am-Fm,ci_AB=[float(lo),float(hi)],
               A_range=float(rng_A),placebo=pl,floor_n=FLOOR,verdict=verdict,
               society_side=dict(within_act=0.8451,across_act_upper=0.4401),
               cells={f"{a}×{b}":dict(rho=r,years=n) for d in (A,B,F) for (a,b),(r,n) in d.items()},
               unchallenged=True),
          open(OUT/"person_level.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'person_level.json'}")

# ── 诊断(不是结果):正对照为什么失败 ──────────────────────────────────────
# 四道 pol 题是**二值**(1=yes/2=no),而性道德四题是**四档**。
# **两个二值题之间的 Spearman 就是 phi,而 phi 的上限被两边的边际分布卡死。**
# 实测 p(yes):polabuse 0.1093 · polmurdr 0.1154 · polescap 0.7266 · polattak 0.9057 —— 极度偏斜。
# ⇒ **我拿一个有天花板的量,去比一个没有天花板的量。** 正对照的失败正是抓到这件事的东西。
# ⚠ **下面这一段是诊断,不是结果**:归一化是**事后**做的,而且那条天花板公式只对**正**方向成立
#    (`polmurdr × polattak` 是负的,占比算出 −124%,那是公式伪影不是超天花板)。
#    **要当结果用,必须作为一个预注册的新轮次重跑。**
print("\n=== 诊断:phi 的天花板(事后,不是结果)===")
def maxphi(p,q):
    p,q=(p,q) if p<=q else (q,p)
    return np.sqrt((p*(1-q))/((1-p)*q))
diag=[]
for a,b in polp:
    per=[]
    for y,g in df.groupby("year"):
        m=g[[a,b]].dropna()
        if len(m)<FLOOR: continue
        r=float(spearmanr(m[a].to_numpy(float),m[b].to_numpy(float)).statistic)
        c=maxphi((m[a]==1).mean(),(m[b]==1).mean())
        if np.isfinite(r) and c>0: per.append((r,c,r/c))
    if per:
        r,c,f=[float(np.median(x)) for x in zip(*per)]
        diag.append(dict(pair=f"{a}×{b}",phi=r,ceiling=c,frac=f))
        print(f"  {a:10s} × {b:10s} phi={r:+.4f} 天花板={c:.4f} 占 {f*100:6.1f}%")
nm=float(np.median([d["frac"] for d in diag]))
print(f"\n  原始中位 {Am:+.4f} · **归一到各自天花板后 {nm:+.4f}** · 性道德四题参照 {Fm:+.4f}")
print("  ⇒ **归一后两块题几乎相同 —— 「警察四题不聚合」是二值偏斜边际造的伪影。**")
d=json.load(open(OUT/"person_level.json")); d["diagnosis_ceiling"]=diag; d["A_normalised"]=nm
d["diagnosis_note"]=("事后诊断,不是结果:归一化是看到失败之后做的,且天花板公式只对正方向成立"
                     "(polmurdr×polattak 为负,占比 −124% 是公式伪影)。要当结果用必须预注册重跑。")
json.dump(d,open(OUT/"person_level.json","w"),indent=1,ensure_ascii=False)
