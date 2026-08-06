"""E03·A18·R688 —— 修完仪器再问一次:这一页的每一个「领域」,是谁分的组

**类型:FRONTIER**。`#651` 是 UNVERIFIED,两处坏都是我的。本轮**两处修法先写死**,重跑。

⚠ **BASIN 仍然下注 W1**(孤例),即仍然下注 `#650` 那句「三次同一个错误类」是过度声称。
  **而 MFQ 五域从此不再是正对照,它们变成被测对象** —— 于是「FAIRNESS/INGROUP 只有 5/6」
  第一次**有资格成为一个结果**,而它若成立,**W1 就输了**。**我下注的仍是我不希望赢的那一面的反面。**

## 修法①:地板按站点,而且不可算就整块判不了

`#651` 的地板写在 150,而 SCCS 最多只有 186 个社会 -> 两块**一对都没算出来**,
**而它没有报错,它报了 `1/1/1/1`,看起来像一个结论。**
⇒ **人层 150 · 社会层 30**(A12 全程用的就是 30);
⇒ **任何块若有任一对不可算,整块记「判不了」,不许输出一个数。**

## 修法②:正对照换成有独立已知答案的块

| 臂 | 块 | 独立已知答案 | 要求 |
|---|---|---|---|
| **A'** | `SCCS·以身作则四对象` | `#642` 实测六对全在 **+0.85** 以上 | **必须 4/4** |
| **A''** | `SCCS·体罚四对象` | `#640` 实测六对全在 **+0.79** 以上 | **必须 4/4**(A' 的复核) |
| **B'** | `NSFG·家庭七题` | `#686` 实测 **4/7** | **必须 4/7** |

**MFQ 五域不再当正对照** —— 它们的**完全连通性从来没有被独立确立过**(`#651` 学到的)。

## 修法③:判据以归一阈为准

`#651` 第四节已证生阈会把天花板低的块系统性误判(GSS 警察块 生 2/4 · 归一 4/4)。
**而生阈与归一阈分歧的块必须单独列出。**

G1 ESTIMAND:最大连通子块 ÷ 名义块大小,**归一阈 0.30**。
KILL(条件式):if A'=1.00 and A''=1.00 and B'=4/7 and 安慰剂=1:
  **除 B' 外还有 >= 1 块 < 1.00 -> W2(系统性)** · **全部 = 1.00 -> W1(孤例,我过度声称了)**
  else: UNVERIFIED
G3:14 块全报。G4:阈 {0.20,0.25,0.30,0.35} × {生, 归一} 八条规格,每块全报。
IMPOSSIBLE(不写 planned):单波/汇总,无跨时间复核 · MFQ 非概率 · SCCS 同一编码团队 ·
  「谁分的组」是我读元数据读出来的 · `[unchallenged]`
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
    # 修法①:floor 由调用方按站点给;缺一对即整块判不了(见下)
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
FLOOR={k:(30 if k.startswith("SCCS·") else 150) for k in BLK}

print("=== 硬规则①:每个块是谁分的组,名义大小,以及它的对数 ===")
ROWS={}
for name,(F,items,who,yc) in BLK.items():
    A,flip=align(F,items)
    Fa=A.join(F[[yc]]) if yc else A
    pm=pairmap(Fa,items,year=yc,floor=FLOOR[name])
    need=len(items)*(len(items)-1)//2
    if len(pm)<need:                      # 修法①:缺一对即整块判不了,不许输出一个数
        ROWS[name]=dict(who=who,nominal=len(items),npairs=len(pm),undecidable=True,
                        ratio_raw=None,ratio_norm=None,spec={},flipped=flip,
                        raw_med=None,norm_med=None)
        print(f"  {name:18s} {who:16s} 名义 {len(items)}  对 {len(pm)}/{need}  **判不了(缺对)**"); continue
    row=dict(who=who,nominal=len(items),npairs=len(pm),flipped=flip,
             spec={}, raw_med=float(np.median([v[0] for v in pm.values()])) if pm else np.nan,
             norm_med=float(np.median([v[1] for v in pm.values()])) if pm else np.nan)
    for idx,tag in ((0,"raw"),(1,"norm")):
        for t in THR: row["spec"][f"{tag}@{t:.2f}"]=biggest(pm,items,t,idx)
    row["ratio_raw"]=row["spec"]["raw@0.30"]/len(items)
    row["ratio_norm"]=row["spec"]["norm@0.30"]/len(items)
    row["undecidable"]=False
    ROWS[name]=row
    print(f"  {name:18s} {who:16s} 名义 {len(items)}  对 {len(pm):2d}  翻向 {flip or '无'}")

OK={k:r for k,r in ROWS.items() if not r["undecidable"]}
print(f"\n=== G3/G4:{len(OK)} 块可算 × 4 阈 × {{生, 归一}} ===")
print(f"  {'块':18s}{'名义':>4s}{'生中位':>8s}{'归一中位':>9s}   生@.20/.25/.30/.35   归一@.20/.25/.30/.35   {'比值(归一)':>12s}")
for name,r in OK.items():
    sr="/".join(str(r["spec"][f"raw@{t:.2f}"]) for t in THR)
    sn="/".join(str(r["spec"][f"norm@{t:.2f}"]) for t in THR)
    print(f"  {name:18s}{r['nominal']:>4d}{r['raw_med']:>8.3f}{r['norm_med']:>9.3f}        {sr:11s}          {sn:11s}   **{r['ratio_norm']:.3f}**")
und=[k for k,r in ROWS.items() if r["undecidable"]]
print(f"  判不了的块:{und or '无'}")

print("\n=== 控制:三块独立已知答案 ===")
def get(k,f="ratio_norm"): return OK[k][f] if k in OK else None
A1=get("SCCS·以身作则四对象"); A2=get("SCCS·体罚四对象"); B1=get("NSFG·家庭七题")
print(f"  臂 A'  SCCS·以身作则四对象(#642 六对全 >= +0.85)必须 1.00 -> **{A1}**")
print(f"  臂 A'' SCCS·体罚四对象  (#640 六对全 >= +0.79)必须 1.00 -> **{A2}**")
print(f"  臂 B'  NSFG·家庭七题    (#686 实测 4/7)必须 0.571 -> **{None if B1 is None else round(B1,3)}**")
def placebo(seed):
    rng=np.random.default_rng(seed); F,items,_,_=BLK["MFQ·PURITY"]
    Z=F[items].dropna().copy()
    for c in items: Z[c]=rng.permutation(Z[c].to_numpy())
    A,_=align(Z,items)
    return biggest(pairmap(A,items,floor=150),items,0.30,0)
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"  安慰剂 打乱行后最大子块必须 = 1 -> **{pl:.1f}**")

G=Gate("修完仪器再问一次:这一页的每一个领域,是谁分的组")
armsok = (A1==1.0 and A2==1.0 and B1 is not None and abs(B1-4/7)<0.01)
p1=G.positive_control("三块独立已知答案全部复现",planted=float(1.0 if armsok else 0.0),floor=0.5,spread=0.01)
p2=G.negative_control("安慰剂:打乱行后最大子块回到 1",null=pl,effect=6.0,null_spread=0.2,
                      null_kind="行内打乱,保留每题边际")
partial=[k for k,r in OK.items() if r["ratio_norm"]<0.999]
if p1 and p2:
    others=[k for k in partial if k!="NSFG·家庭七题"]
    verdict=(f"**W2 —— 系统性:除臂 B' 外还有 {len(others)} 块名义大于真实** -> {others}" if others
             else "**W1 —— 孤例。`#650` 那句「三次同一个错误类」是一个实例套上了一个模式,我过度声称了。**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)

print("\n=== 生阈与归一阈分歧的块(必须单独列出)===")
dis=[(k,r["spec"]["raw@0.30"],r["spec"]["norm@0.30"]) for k,r in OK.items()
     if r["spec"]["raw@0.30"]!=r["spec"]["norm@0.30"]]
for k,a,b in dis: print(f"  {k:18s} 生 {a} -> 归一 {b}")
print(f"  分歧块数 = **{len(dis)}**" + ("(而生阈会系统性低估,见 `#651`)" if dis else ""))
json.dump(dict(blocks=ROWS,armA1=A1,armA2=A2,armB=B1,placebo=pl,partial=partial,
               undecidable=und,disagree=[[k,a,b] for k,a,b in dis],verdict=verdict,unchallenged=True),
          open(OUT/"whose_grouping_fixed.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'whose_grouping_fixed.json'}")

# ── 哪一道题掉出来了,而硬规则①要求先读题干 ──────────────────────────────────
print("\n=== 掉出来的是哪一道题(硬规则①:先读题干)===")
_,mm=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",metadataonly=True)
LB=mm.column_names_to_labels
dropped={}
for name in ["MFQ·FAIRNESS","MFQ·INGROUP"]:
    F,items,_,_=BLK[name]; A,_=align(F,items); pm=pairmap(A,items,floor=150)
    adj={i:0 for i in items}
    for (a,b),v in pm.items():
        if v[1]>=0.30: adj[a]+=1; adj[b]+=1
    out=[i for i in items if adj[i]==0]
    dropped[name]=out
    print(f"  {name}:")
    for i in items:
        mark="  **<- 掉出来**" if adj[i]==0 else ""
        print(f"    {i:12s} 连到 {adj[i]} 题  「{str(LB.get(i,''))[:66]}」{mark}")
d=json.load(open(OUT/"whose_grouping_fixed.json"))
d["dropped_items"]={k:{i:str(LB.get(i,'')) for i in v} for k,v in dropped.items()}
json.dump(d,open(OUT/"whose_grouping_fixed.json","w"),indent=1,ensure_ascii=False)
