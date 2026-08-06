"""E03·A18·R685 —— 分诊,而分诊器自己要先能失败

**类型:FRONTIER**。`#648`:天花板修正只在边际偏斜的地方要紧 ⇒ **按天花板分诊,不是一律重算。**
**而一句「分诊」如果不能失败,它就只是一次重新措辞。** 所以本轮先给分诊器造一个能失败的正对照。

分诊判据(先写死,只留能分辨的那一刀 —— `#641` 规则,本会话第六次):
  天花板中位 >= 0.85 -> **无需重算**
  天花板中位 <  0.60 -> **必须重算,可能改写**
  0.60 – 0.85        -> **待定,不硬判**

**正对照(两块已知答案的样本)**:
  `#648` 实测 MFQ 五域 **0.913–0.950** -> 必须全部落进 `>=0.85`;
  `#647` 实测 GSS 四道二值警察题 **0.095–0.769**(中位 0.185)-> 必须落进 `<0.60`。
  **分不出这两块,分诊器就是坏的 ⇒ 整轮 UNVERIFIED。**

被分诊的对象(G3:全部报,包括落进「无需重算」的):
  ① MFQ 五域(控制臂 A)· ② GSS 警察四题(控制臂 B)
  ③ **GSS 性道德四题** —— `#647` 实测天花板 0.4957–0.8373,**它跨在两档之间**
  ④ **NSFG 性三题 vs 家庭七题** —— 这一页「第九件」的来源,**从未被分诊过**,
     而两块题的边际差得远(`samesex` 是极化的,`chsuppor` 近乎一致)⇒ **排序完全可能翻**。

⚠ **最强混淆,跑之前写死**:天花板由**观测边际**估出,而有些数是**跨年汇总**的 ——
  `#647` 实测 `polescap` 的 p(yes) 在 **0.211–0.775** 之间移动。
  **在合并边际上算一次天花板会算错。** ⇒ **凡有时间维的逐年算再取中位**(GSS 两块);
  NSFG 与 MFQ 是单波,**如实标注这一条限制而不是假装它不在**。
KILL(条件式):if 控制臂 A 全进 `>=0.85` and 控制臂 B 进 `<0.60`: 按判据分诊 ④ 与 ③,
  **并对落进「必须重算」的块实际重算,报排序是否改变**;else: UNVERIFIED。
IMPOSSIBLE(不写 planned):NSFG/MFQ 单波 ⇒ 无法逐年算天花板 · 非概率(MFQ)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr, rankdata
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)

def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,sign=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if sign<0: y=y[::-1]
    return sp(x,y)
def block(frame, items, year=None, floor=200):
    """返回 (raw 中位, ceiling 中位, norm 中位, 对数)。year 非空则逐年算再取中位。"""
    raw=[];cei=[];nor=[]
    groups=[(None,frame)] if year is None else list(frame.groupby(year))
    for a,b in combinations(items,2):
        per=[]
        for _,g in groups:
            m=g[[a,b]].dropna()
            if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append((r,abs(c),r/abs(c)))
        if per:
            raw.append(np.median([p[0] for p in per])); cei.append(np.median([p[1] for p in per]))
            nor.append(np.median([p[2] for p in per]))
    if not raw: return None
    return dict(raw=float(np.median(raw)),ceil=float(np.median(cei)),
                norm=float(np.median(nor)),npairs=len(raw),ceil_min=float(min(cei)),ceil_max=float(max(cei)))

BLOCKS={}
# ① MFQ 五域(控制臂 A)
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
Xm=mfq[list(ITEM)].dropna()
for dom in sorted(set(ITEM.values())):
    BLOCKS[f"①MFQ·{dom}"]=block(Xm,[k for k,v in ITEM.items() if v==dom])
# ②③ GSS(逐年)
POL=["polabuse","polmurdr","polescap","polattak"]; SEXG=["premarsx","xmarsex","homosex","teensex"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta", usecols=["year"]+POL+SEXG, encoding="latin1")
BLOCKS["②GSS·警察四题(二值)"]=block(g,POL,year="year")
BLOCKS["③GSS·性道德四题(四档)"]=block(g,SEXG,year="year")
# ④ NSFG 性 vs 家庭(单波)
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct", errors="replace"):
    mm=pat.search(line)
    if mm: LAY[mm.group(2).lower()]=(int(mm.group(1))-1,int(mm.group(3)),mm.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEXN+FAMN if n in LAY}
buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat", errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
Xn=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
BLOCKS["④NSFG·性三题"]=block(Xn,SEXN); BLOCKS["④NSFG·家庭七题"]=block(Xn,FAMN)

print("=== G3:全部报,包括落进「无需重算」的 ===")
print(f"  {'块':24s}{'raw':>9s}{'天花板':>9s}{'norm':>9s}{'天花板范围':>20s}{'对':>4s}  分诊")
def triage(c): return "无需重算" if c>=0.85 else ("必须重算" if c<0.60 else "待定")
rows={}
for k,v in BLOCKS.items():
    if v is None:
        print(f"  {k:24s} 判不了"); continue
    t=triage(v["ceil"]); rows[k]=dict(v,triage=t)
    rng=f"[{v['ceil_min']:.3f}, {v['ceil_max']:.3f}]"
    print(f"  {k:24s}{v['raw']:>9.4f}{v['ceil']:>9.4f}{v['norm']:>9.4f}{rng:>20s}{v['npairs']:>4d}  **{t}**")

print("\n=== 正对照:分诊器分得出两块已知答案的样本吗 ===")
armA=[k for k in rows if k.startswith("①")]
armB="②GSS·警察四题(二值)"
okA=all(rows[k]["ceil"]>=0.85 for k in armA)
okB=rows[armB]["ceil"]<0.60
print(f"  臂 A · MFQ 五域必须全 >=0.85:{[round(rows[k]['ceil'],4) for k in armA]} -> **{okA}**")
print(f"  臂 B · GSS 警察四题必须 <0.60:{rows[armB]['ceil']:.4f} -> **{okB}**")

G=Gate("分诊器自己能不能失败")
p1=G.positive_control("臂 A:MFQ 五域全部落进 >=0.85",planted=float(min(rows[k]["ceil"] for k in armA)),
                      floor=0.85,spread=0.005)
p2=G.negative_control("臂 B:GSS 二值块必须落进 <0.60(它不该被判成「无需重算」)",
                      null=float(rows[armB]["ceil"]),effect=float(min(rows[k]["ceil"] for k in armA)),
                      null_spread=0.02,null_kind="已知偏斜的二值块,天花板必然低")

print("\n=== ④ 这一页「第九件」的两块题:排序会不会翻 ===")
s3,f7=rows["④NSFG·性三题"],rows["④NSFG·家庭七题"]
print(f"  未归一:性 {s3['raw']:+.4f} vs 家庭 {f7['raw']:+.4f}  ->  性{'更紧' if s3['raw']>f7['raw'] else '更松'}")
print(f"  归一后:性 {s3['norm']:+.4f} vs 家庭 {f7['norm']:+.4f}  ->  性{'更紧' if s3['norm']>f7['norm'] else '更松'}")
flip = (s3["raw"]>f7["raw"]) != (s3["norm"]>f7["norm"])
print(f"  天花板:性 {s3['ceil']:.4f} · 家庭 {f7['ceil']:.4f}  ->  **排序{'翻了 ⇒ 第九件必须改写' if flip else '没翻'}**")

if p1 and p2:
    verdict=("**分诊器可用。第九件的排序" + ("**翻了 ⇒ 必须改写**" if flip else "未翻 ⇒ 结论不变**") )
else:
    verdict="UNVERIFIED —— 分诊器分不出已知答案的两块,它是坏的"
print(f"\n{verdict}"); print(G)
json.dump(dict(blocks=rows,armA_ok=bool(okA),armB_ok=bool(okB),nsfg_flip=bool(flip),
               verdict=verdict,unchallenged=True,
               limit="NSFG/MFQ 单波 ⇒ 天花板无法逐年算;GSS 两块已逐年算再取中位"),
          open(OUT/"triage.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'triage.json'}")
