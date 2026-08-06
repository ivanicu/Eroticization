"""E03·A32·R167 —— 一个社会教孩子的十件事里,它最按性别区别对待哪一件

**类型:FRONTIER。**

**心理学的那一句:一个社会把男孩女孩分开来教的,首先不是勇气也不是服从 ——
要看数据才知道是不是性。而「双重标准」若只是十件事里普通的一件,页上那句新话就要收窄。**

## 缺口:`#724`① 的地板问题被**解散**,不是绕过
`#724` 用 Whyte 1978 的跨队臂量出「有双标的社会 |女−男| 中位 0.750」,而**某一臂只有 22–28 个社会**,
按预注册只作观察。`#724`① 要求「要么找第三个编码,要么记为测不到」。
⚠ **第三条路:`SCCS330–333` 自己就是男女差的直接测量**,联合 **n=147**,**根本不需要 Whyte**。
而且码本给了更好的东西:**十个品格共用同一条 0–10「灌输强度」量表** ⇒
**可以问「十件事里社会最按性别区别对待的是哪一件」,而这是一个比 0.750 大得多的问题。**

## G1 ESTIMAND(先命名,后选统计量)
① **有符号男女差** `g = 女 − 男`,逐社会、逐品格、逐年龄段(早/晚)。
② **普遍还是两极**:`sign(g)` 的分布 —— 正/零/负各多少社会。
③ **十个品格的性别化名次**:按 `median(g)` 与 `median(|g|)` 各排一次(**两个都排,不许只报一个**)。
④ **差与水平的关系**:`g` 对该社会该品格的**平均水平**回归的斜率。

## W1–W4(预测矩阵;`#724` 之后的世界)
| 世界 | ④ 斜率 | ③ 性克制的名次 | 读法 |
|---|---|---|---|
| **W1 强化** | **正**,且大于 sham | 高 | 越严的社会,额外加在女孩身上的越多 |
| **W2 独立轴** | ≈0(相对 sham) | 高 | 双标是一条与严厉度正交的轴 |
| **W3 天花板** | **负**,但 sham 也负 | 高 | 斜率是量表的,不是社会的 |
| **W4 元分离器** | 任意 | **不在前三** | **性不是被性别化最狠的那一件 ⇒ 页上那句新话要收窄** |

⚠ **W4 的正结果我不高兴** —— 它直接削 `#724` 刚上页的那一行。
⚠ **而这就是元分离器本身**:前三个世界都预设「双标是一件关于性的事」;W4 说这个分解方式本身可能是错的。

## G2 CONTROLS
- **④ 正对照(来自文献,不是我编的)**:**攻击性(Aggression)必须是男 > 女** ——
  这是民族志里最被记录的性别差之一。**若这具仪器连它都测不出,它测不了任何性别化。**
  ⚠ **且必须在 g=0 时失败**:把每个社会内的男/女标签**打乱**后,攻击性的差必须回到零。
- **零** = `negative_control`,**零的种类 = 在每个社会内部打乱「男/女」这个标签,
  保住该社会该品格的两个取值与它的总体水平,只毁掉「哪个值属于哪个性别」。**
- **SHAM(最强混淆的控制,同一迭代内)**:**九个非性品格**走完全相同的流程 ——
  同一批编码者、同一条 0–10 量表、同样的有界差压缩。**「两个有界分数之差在量表两端被机械压缩」
  这条代数因此在 sham 上同样成立,所以任何「性克制特殊」的说法必须相对 sham 成立。**
- **PLACEBO**:`早男 − 晚男` 是一个**年龄**差,不是性别差 ⇒ 它在性别零下**不该**回到零,
  用来证明零只毁掉了性别标签而没有毁掉别的。
## G3 全格:10 品格 × 2 年龄 = **20 格全报**,含不支持结论的。
## G4 规格曲线:`median(g)` / `median(|g|)` / **按剩余空间归一的 g**(应对天花板)× 早/晚 × 有符号/绝对值。

## ⑤ 停止条件(跑之前写死,不许跑完再找理由)
- **攻击性的男>女差在真数据上不成立,或在打乱标签后仍成立 ⇒ 记 UNVERIFIED 并停。**
- **性克制若不在 `median(|g|)` 的前三 ⇒ 判 W4**,如实报名次,**页上那一行要收窄**。
- **④ 的斜率:性克制的斜率若落在九个 sham 品格斜率的展布之内 ⇒ 「双标随严厉度变化」判不了。**

## IMPOSSIBLE(不写 planned)
**换不了仪器**:只有 `barry1976traits` 编过「品格 × 四对象」的 0–10 网格(`#700` 已枚举);
Whyte 1978 只有二值/三值,**不能替代**。横断面无因果。**「灌输强度」是编码者读民族志给的分,
不是被观察的行为** —— 这是仪器,必须点名(硬规则②)。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
B="data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B+"variables.csv",low_memory=False); Dd=pd.read_csv(B+"data.csv")
W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),{})[f"{m.group(2)}_{m.group(3)}"]=r.id
FAM={k:v for k,v in fam.items() if len(v)==4}
print(f"=== 硬规则①:`barry1976traits` · 十个品格 × 四对象 · 共用 0–10「灌输强度」量表 ===")
print(f"{'品格':18s}{'早男/早女 n':>12s}{'晚男/晚女 n':>12s}{'早配对':>8s}{'晚配对':>8s}")
PAIR={}
for k,d in FAM.items():
    e=W[[d["Early_Boy"],d["Early_Girl"]]].dropna(); l=W[[d["Late_Boy"],d["Late_Girl"]]].dropna()
    PAIR[k]=dict(early=e.rename(columns={d["Early_Boy"]:"boy",d["Early_Girl"]:"girl"}),
                 late=l.rename(columns={d["Late_Boy"]:"boy",d["Late_Girl"]:"girl"}))
    print(f"{k:18s}{int(W[d['Early_Boy']].notna().sum()):5d}/{int(W[d['Early_Girl']].notna().sum()):<6d}"
          f"{int(W[d['Late_Boy']].notna().sum()):5d}/{int(W[d['Late_Girl']].notna().sum()):<6d}{len(e):>8d}{len(l):>8d}")
MAXC=10.0
def cells(fr):
    g=(fr.girl-fr.boy).to_numpy(float); lev=((fr.girl+fr.boy)/2).to_numpy(float)
    head=MAXC-lev                                   # 剩余空间(天花板对照用)
    return g,lev,head
res={}
for k in FAM:
    for age in ("early","late"):
        g,lev,head=cells(PAIR[k][age])
        sl=np.polyfit(lev,g,1)[0] if len(g)>10 and np.std(lev)>0 else np.nan
        res[(k,age)]=dict(n=len(g),med=float(np.median(g)),med_abs=float(np.median(np.abs(g))),
            mean=float(g.mean()),pos=int((g>0).sum()),zero=int((g==0).sum()),neg=int((g<0).sum()),
            slope=float(sl),med_norm=float(np.median(g/np.maximum(head,1e-9))))
print(f"\n=== G3 全格 20 格(10 品格 × 2 年龄),含不支持结论的 ===")
print(f"{'品格':18s}{'年龄':>5s}{'n':>5s}{'中位 g':>9s}{'中位|g|':>9s}{'女>男':>7s}{'=':>5s}{'男>女':>7s}{'斜率':>9s}")
for k in FAM:
    for age in ("early","late"):
        r=res[(k,age)]
        star=" ★" if k=="Sexual Restraint" else (" ⊕正对照" if k=="Aggression" else "")
        print(f"{k:18s}{age:>5s}{r['n']:>5d}{r['med']:>+9.3f}{r['med_abs']:>9.3f}"
              f"{r['pos']:>7d}{r['zero']:>5d}{r['neg']:>7d}{r['slope']:>+9.4f}{star}")
# ④ 正对照 + 零
rng=np.random.default_rng(20260806)
def shuffled_stat(fr,B=4000):
    a=fr.boy.to_numpy(float); b=fr.girl.to_numpy(float)
    out=[]
    for _ in range(B):
        f=rng.random(len(a))<0.5
        gb=np.where(f,b,a); gg=np.where(f,a,b)
        out.append(float(np.median(gg-gb)))
    return np.array(out)
print(f"\n=== ④ 正对照(文献先验):攻击性必须是 **男 > 女**,且打乱性别标签后必须回到零 ===")
agg_ok=True
for age in ("early","late"):
    fr=PAIR["Aggression"][age]; obs=res[("Aggression",age)]["med"]
    nul=shuffled_stat(fr); q=float(np.quantile(np.abs(nul),0.95))
    ok=(obs<0) and (abs(obs)>q)
    agg_ok&=ok
    print(f"  Aggression {age:5s}: 中位 g = **{obs:+.3f}**(负 = 男>女)· 打乱标签后零的 95% 分位 {q:.3f} · "
          f"中位 {np.median(nul):+.4f} ⇒ {'✅' if ok else '⛔'}")
print(f"\n=== PLACEBO:早男 − 晚男 是**年龄**差,性别零毁不掉它 ===")
for k in ("Aggression","Sexual Restraint"):
    d=FAM[k]; J=W[[d["Early_Boy"],d["Late_Boy"]]].dropna()
    ag=(J[d["Late_Boy"]]-J[d["Early_Boy"]]).median()
    print(f"  {k:18s} 晚男 − 早男 中位 = **{ag:+.3f}**(n={len(J)})—— 性别零不触及这一维")
print(f"\n=== 性克制在十个品格里的性别化名次(G4:两种口径各排一次)===")
for age in ("early","late"):
    for key,nm in (("med_abs","中位|g|"),("med","中位 g(有符号,负=男>女)")):
        od=sorted(FAM,key=lambda k:-abs(res[(k,age)][key]) if key=="med" else -res[(k,age)][key])
        rk=od.index("Sexual Restraint")+1
        print(f"  {age:5s} 按 {nm:22s} ⇒ 性克制第 **{rk}/10** · 前三:{', '.join(od[:3])}")
sr_ranks=[]
for age in ("early","late"):
    od=sorted(FAM,key=lambda k:-res[(k,age)]["med_abs"]); sr_ranks.append(od.index("Sexual Restraint")+1)
# 斜率 sham 对照
sham=[res[(k,a)]["slope"] for k in FAM if k!="Sexual Restraint" for a in ("early","late")]
sr_sl=[res[("Sexual Restraint",a)]["slope"] for a in ("early","late")]
print(f"\n=== ④ 斜率 vs SHAM(九个非性品格,同量表同编码者)===")
print(f"  sham 斜率 18 格:中位 {np.median(sham):+.4f} · 范围 [{min(sham):+.4f}, {max(sham):+.4f}]")
print(f"  性克制:早 {sr_sl[0]:+.4f} · 晚 {sr_sl[1]:+.4f} ⇒ "
      f"{'**落在 sham 展布之内 ⇒ 「双标随严厉度变化」判不了**' if all(min(sham)<=s<=max(sham) for s in sr_sl) else '**在 sham 之外**'}")
G=Gate("十件事里社会最按性别区别对待哪一件")
p1=G.positive_control("攻击性必须是男>女,且打乱性别标签后回到零(文献先验)",
    planted=1.0 if agg_ok else 0.0,floor=0.0,spread=0.1)
fr=PAIR["Sexual Restraint"]["late"]; nul=shuffled_stat(fr)
p2=G.negative_control("打乱每个社会内的男/女标签后,性克制的差应回到零",
    null=float(np.quantile(np.abs(nul),0.95)),effect=abs(res[("Sexual Restraint","late")]["med"]),
    null_spread=0.02,null_kind="在每个社会内部打乱「男/女」标签 —— 保住该社会的两个取值与总体水平,只毁掉「哪个值属于哪个性别」")
if not p1: v="**UNVERIFIED:正对照没过 —— 这具仪器测不出已知的性别差,任何名次都不可采信**"
elif min(sr_ranks)<=3: v=f"**性克制在 median(\\|g\\|) 上排第 {sr_ranks[0]}(早)/ {sr_ranks[1]}(晚)—— 它在前三**"
else: v=(f"**W4:性克制排第 {sr_ranks[0]}(早)/ {sr_ranks[1]}(晚),**不在前三** ⇒ "
         f"性不是被性别化最狠的那一件,页上那一行要收窄**")
print(f"\n{v}"); print(G)
json.dump(dict(cells={f"{k}|{a}":res[(k,a)] for k,a in res},sr_ranks=sr_ranks,
   sham_slopes=sham,sr_slopes=sr_sl,agg_positive_control=bool(agg_ok),verdict=v,unchallenged=True),
   open(OUT/"gendered.json","w"),indent=1,ensure_ascii=False)
