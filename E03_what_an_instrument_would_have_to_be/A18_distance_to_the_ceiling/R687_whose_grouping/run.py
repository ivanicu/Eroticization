"""E03·A18·R687 —— 这一页的每一个「领域」,是谁分的组,而它在数据里成不成立

**类型:FRONTIER**。`#650` 点出一个跨三轮复发的错误类(`#528` 笔迹 · `#642` 光环 · `#650` 袋子)——
**把仪器的组织方式当成世界的组织方式**。本轮把它做成一个可测的量。

⚠ **BASIN 检查,写在最前面**:最近几轮**都在确认同一个故事**。这就是盆地。
  **所以本轮下注 W1** —— 即下注我上一轮那句「三次同一个错误类」是**一个实例套上了一个模式**。

W1 **孤例** —— 除 NSFG 家庭七题外,所有块给 1.00 ⇒ **「错误类」只有一个测过的实例,另外两个是类比。**
W2 **系统性** —— 多块名义大于真实 ⇒ 多条声明要重看。
W3 **算法是伪影** —— 比值随阈乱跳(由规格曲线 + 两个控制臂抓)。

G1 ESTIMAND(先于方法):对每一个**继承自仪器分组**的块,
  **最大连通子块题数 ÷ 名义块题数**(`#686` 的算法:按与块均值的相关符号对齐 -> 阈上建图 -> 最大连通分量)。
G2 CONTROLS:
  **臂 A** MFQ 五域必须给 **1.00**(`#608`/`#648` 已证成块)。
  **臂 B** NSFG 家庭七题必须给 **4/7 = 0.571**(`#686` 实测)。**分不出这两块 ⇒ 算法坏 ⇒ UNVERIFIED。**
  **安慰剂** 任一块打乱行后最大子块必须回到 **1**。
G3:14 个块全报,包括给 1.00 的。G4:阈 {0.20, 0.25, 0.30, 0.35} × {生相关, 归一相关} 八条规格。
KILL(条件式):if 臂A全=1.00 and 臂B=0.571 and 安慰剂=1:
  **除臂 B 外还有 >= 1 块 < 1.00 -> W2** · **全部 = 1.00 -> W1(我上一轮过度声称了)**
  else: UNVERIFIED
⚠ **两个必须先写下的混淆:**
  ① **阈敏感** —— `#686` 里 0.20–0.35 恰好都给 4,**那是那一块的性质,不是算法的性质** ⇒ 每块四阈全报。
  ② **阈是绝对的,而各块天花板不同** —— GSS 二值块天花板 **0.3552**(`#685`),
     **生相关永远够不到 0.30,会被系统性误判成「不成块」**。
     ⇒ **归一阈与生阈两条都跑**;控制臂在**生阈**上评(才对得上 `#686` 的已知答案),
     **而结论以归一阈为准,并把两者的分歧公布。**
IMPOSSIBLE(不写 planned):三个仪器都是单波或汇总,**无跨时间复核** · MFQ 非概率 ·
  SCCS 是同一编码团队(`#642` 已量) · 「谁分的组」是我读元数据读出来的,不是数据说的 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; THR=(0.20,0.25,0.30,0.35)

def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,sign=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if sign<0: y=y[::-1]
    return sp(x,y)
def align(F,items):
    Z=F[items].dropna(); m=Z.mean(axis=1)
    flip=[i for i in items if sp(Z[i],m)<0]
    A=Z.copy()
    for i in flip: A[i]=-A[i]
    return A,flip
def pairmap(F,items,year=None,floor=150):
    """返回 {(a,b):(raw, norm)};year 非空则逐年算再取中位。"""
    out={}
    groups=[(None,F)] if year is None else list(F.groupby(year))
    for a,b in combinations(items,2):
        per=[]
        for _,g in groups:
            m=g[[a,b]].dropna()
            if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append((r,r/abs(c)))
        if per: out[(a,b)]=(float(np.median([p[0] for p in per])),float(np.median([p[1] for p in per])))
    return out
def biggest(pm,items,thr,idx):
    adj={i:set() for i in items}
    for (a,b),v in pm.items():
        if v[idx]>=thr: adj[a].add(b); adj[b].add(a)
    seen=set(); best=1
    for i in items:
        if i in seen: continue
        st=[i]; c=0
        while st:
            u=st.pop()
            if u in seen: continue
            seen.add(u); c+=1; st+=[v for v in adj[u] if v not in seen]
        best=max(best,c)
    return best

BLK={}   # name -> (frame, items, 谁分的组, year_col)
# ── MFQ(问卷作者 Graham/Haidt/Nosek)
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
Xm=mfq[list(ITEM)].dropna()
for dom in sorted(set(ITEM.values())):
    BLK[f"MFQ·{dom}"]=(Xm,[k for k,v in ITEM.items() if v==dom],"问卷作者 GHN 2009",None)
# ── GSS(问卷作者 NORC)
POL=["polabuse","polmurdr","polescap","polattak"]; SEXG=["premarsx","xmarsex","homosex","teensex"]
gss,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+POL+SEXG,encoding="latin1")
BLK["GSS·警察四题"]=(gss,POL,"问卷作者 NORC","year"); BLK["GSS·性道德四题"]=(gss,SEXG,"问卷作者 NORC","year")
# ── NSFG(问卷作者 NCHS)
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEXN+FAMN if n in LAY}; buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
Xn=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
BLK["NSFG·性三题"]=(Xn,SEXN,"问卷作者 NCHS",None); BLK["NSFG·家庭七题"]=(Xn,FAMN,"问卷作者 NCHS",None)
# ── SCCS(编码团队 Barry 1977):每种手段的四个对象
S=pathlib.Path("data/external/dplace/repo/datasets/SCCS")
D=pd.read_csv(S/"data.csv"); W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
for nm,b in [("体罚",453),("以身作则",429),("讲课",437),("放任",465),("疼爱",469)]:
    BLK[f"SCCS·{nm}四对象"]=(W,[f"SCCS{b+i}" for i in range(4)],"编码团队 Barry 1977",None)

print("=== 硬规则①:每个块是谁分的组,名义大小,以及它的对数 ===")
ROWS={}
for name,(F,items,who,yc) in BLK.items():
    A,flip=align(F,items)
    Fa=A.join(F[[yc]]) if yc else A
    pm=pairmap(Fa,items,year=yc)
    row=dict(who=who,nominal=len(items),npairs=len(pm),flipped=flip,
             spec={}, raw_med=float(np.median([v[0] for v in pm.values()])) if pm else np.nan,
             norm_med=float(np.median([v[1] for v in pm.values()])) if pm else np.nan)
    for idx,tag in ((0,"raw"),(1,"norm")):
        for t in THR: row["spec"][f"{tag}@{t:.2f}"]=biggest(pm,items,t,idx)
    row["ratio_raw"]=row["spec"]["raw@0.30"]/len(items)
    row["ratio_norm"]=row["spec"]["norm@0.30"]/len(items)
    ROWS[name]=row
    print(f"  {name:18s} {who:16s} 名义 {len(items)}  对 {len(pm):2d}  翻向 {flip or '无'}")

print(f"\n=== G3/G4:14 块 × 4 阈 × {{生, 归一}} ===")
print(f"  {'块':18s}{'名义':>4s}{'生中位':>8s}{'归一中位':>9s}   生@.20/.25/.30/.35   归一@.20/.25/.30/.35   {'比值(归一@.30)':>14s}")
for name,r in ROWS.items():
    sr="/".join(str(r["spec"][f"raw@{t:.2f}"]) for t in THR)
    sn="/".join(str(r["spec"][f"norm@{t:.2f}"]) for t in THR)
    print(f"  {name:18s}{r['nominal']:>4d}{r['raw_med']:>8.3f}{r['norm_med']:>9.3f}        {sr:11s}          {sn:11s}   **{r['ratio_norm']:.3f}**")

print("\n=== 控制 ===")
armA=[k for k in ROWS if k.startswith("MFQ·")]
okA=all(ROWS[k]["ratio_raw"]==1.0 for k in armA)
armB=ROWS["NSFG·家庭七题"]["ratio_raw"]
print(f"  臂 A · MFQ 五域生阈比值必须全 = 1.00:{[ROWS[k]['ratio_raw'] for k in armA]} -> **{okA}**")
print(f"  臂 B · NSFG 家庭七题必须 = 4/7 = 0.571:**{armB:.3f}** -> **{abs(armB-4/7)<0.01}**")
def placebo(seed):
    rng=np.random.default_rng(seed); F,items,_,_=BLK["MFQ·PURITY"]
    Z=F[items].dropna().copy()
    for c in items: Z[c]=rng.permutation(Z[c].to_numpy())
    A,_=align(Z,items)
    return biggest(pairmap(A,items),items,0.30,0)
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"  安慰剂 · 打乱行后最大子块必须 = 1:**{pl:.1f}**")

G=Gate("这一页的每一个领域,是谁分的组")
p1=G.positive_control("臂 A:MFQ 五域生阈比值全 = 1.00",planted=float(min(ROWS[k]['ratio_raw'] for k in armA)),
                      floor=0.99,spread=0.002)
p2=G.negative_control("安慰剂:打乱行后最大子块回到 1",null=pl,effect=float(ROWS["MFQ·PURITY"]["spec"]["raw@0.30"]),
                      null_spread=0.2,null_kind="行内打乱,保留每题边际")
broken=[k for k,r in ROWS.items() if r["ratio_norm"]<0.999]
armB_ok=abs(armB-4/7)<0.01
if p1 and p2 and armB_ok:
    others=[k for k in broken if k!="NSFG·家庭七题"]
    verdict=(f"**W2 —— 系统性:除臂 B 外还有 {len(others)} 块名义大于真实** {others}" if others
             else "**W1 —— 孤例。我上一轮那句「三次同一个错误类」是一个实例套上了一个模式。**")
else:
    verdict=f"UNVERIFIED —— 控制未齐(臂A {p1} · 安慰剂 {p2} · 臂B {armB_ok})"
print(f"\n{verdict}"); print(G)
json.dump(dict(blocks=ROWS,armA_ok=bool(okA),armB=armB,placebo=pl,broken=broken,verdict=verdict,
               unchallenged=True),open(OUT/"whose_grouping.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'whose_grouping.json'}")
