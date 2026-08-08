"""E03·A21·R112 —— 「性是一块」,是不是所有人身上都是一块

**类型:FRONTIER**。**A21 开弧。** `#669` 关掉 A20 时指出:「性是一块」是这一页
**唯一真正跨仪器(T1)**的声明(GSS 最弱一环 0.416 · NSFG 0.346)——
**而它从来没有被问过那个心理学问题。**

⚠ **BASIN**:连续六次下注反对自己喜欢的结果、五次输。**而 W1 才是我想要的**(它会开一条新线)
   ⇒ **仍下注 W2/W3。下注是纪律,不是预测。**
W1 因人而异(两具仪器同号且各自不含零)· W2 对所有人一样 ·
**W3 = meta-separator:梯度只是构成** —— 若它跟着样本构成走而不跟着分层变量的含义走,
**「谁」这个分解又一次是错的**(`#663` 已用这一条杀过一轮)。

## 硬规则①(已跑)
GSS 每层每题真实 n 全在 **4,000–17,000** ✅
⚠ **NSFG 的 `attnd14` 是「14 岁时的礼拜出席」,不是当前出席** —— 与 GSS 的 `attend` **不是同一构念**,
  且 n 仅 2,006 ⇒ **跨仪器只用「年龄」与「教育」两层**;
  **性别是 GSS 独有(NSFG 全为女性)⇒ 单列,不计入跨仪器判据。**

## G1 ESTIMAND(先于方法)
**每层内的「最弱一环」** = 该层内四题(GSS)/三题(NSFG)两两**天花板归一**相关的**最小值**
(归一 = `#647` 的共单调可达上限)。GSS **逐调查年内**算再取中位(与 `#653` 同法);NSFG 单波直接算。
**主量 = 层间梯度**(最高层 − 最低层的最弱一环)。
## G2 CONTROLS
**正对照**:两具仪器的**总体**最弱一环必须复现 **0.416 / 0.346**(容差 0.03)。
**安慰剂**:同一套分层作用在 GSS 的**非性道德题**(`cappun` `gunlaw` `letdie1` `suicide1`)上,
  梯度应显著更小(`#663` 已证死刑对分层不敏感)。**这个零该不该是零?** 该 ⇒ `negative_control`。
## G3/G4:GSS 四个分层 + NSFG 两个分层全报;{年龄, 教育} 为跨仪器规格。
## KILL(条件式)
if 正对照复现 and 安慰剂梯度 < 0.5×主量:
  **两具仪器上同一分层的梯度同号 且 各自 bootstrap 区间不含零** -> **W1**
  任一具含零 -> **判不了,「性是一块」保持为无条件声明**
  梯度最大的层 == 构成差最大的层 -> **W3,不许解释成心理**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
NSFG 全为女性 ⇒ **性别层跨不了仪器** · NSFG 单波 ⇒ 无年内复核 ·
两具仪器的题目不同(4 题 vs 3 题)⇒ **比的是「最弱一环」这个量,不是同一批题** · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; FLOOR=200
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if s<0: y=y[::-1]
    return sp(x,y)
def weakest(frame,items,byyear=None,floor=FLOOR):
    vals={}
    for a,b in combinations(items,2):
        per=[]
        for _,g in ([(None,frame)] if byyear is None else frame.groupby(byyear)):
            m=g[[a,b]].dropna()
            if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
        if per: vals[(a,b)]=float(np.median(per))
    return (min(vals.values()) if vals else np.nan), vals

# ── GSS ──
G=["premarsx","xmarsex","homosex","teensex"]; NS_=["cappun","gunlaw","letdie1","suicide1"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","age","sex","educ","attend"]+G+NS_, encoding="latin1")
d=df[df.year>=1990].copy()
d["年龄段"]=pd.cut(d.age,[17,34,49,64,99],labels=["18-34","35-49","50-64","65+"])
d["性别"]=d.sex.map({1:"男",2:"女"})
d["教育"]=pd.cut(d.educ,[-1,11,12,15,20],labels=["<高中","高中","部分大学","大学+"])
d["礼拜"]=pd.cut(d.attend,[-1,1,3,5,8],labels=["几乎不","偶尔","常去","每周+"])
ov_g,_=weakest(d,G,byyear="year")
print(f"=== 正对照:GSS 总体最弱一环 = **{ov_g:+.4f}**(`#653`: 0.416)")
rows=[]
for s in ["年龄段","性别","教育","礼拜"]:
    for lv in sorted(d[s].dropna().unique(), key=str):
        sub=d[d[s]==lv]
        w,_=weakest(sub,G,byyear="year")
        wn,_=weakest(sub,NS_,byyear="year")
        rows.append(dict(inst="GSS",strat=s,level=str(lv),weakest=float(w) if np.isfinite(w) else None,
                         placebo=float(wn) if np.isfinite(wn) else None,
                         age=float(sub.age.median()),educ=float(sub.educ.median())))
print("\n=== G3:GSS 每层的最弱一环 · 非性道德题的同一量 · 构成 ===")
for r in rows:
    print(f"  {r['strat']:5s}{r['level']:10s} 最弱 {r['weakest']:+.4f} · 非性 "
          f"{(r['placebo'] if r['placebo'] is not None else float('nan')):+.4f} · 年龄中位 {r['age']:.0f} 教育中位 {r['educ']:.0f}")

# ── NSFG ──
NSp=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NSp/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)))
N=["samesex","sxok18","sxok16"]; EX=["ager","hieduc"]
cols={n:LAY[n] for n in N+EX if n in LAY}; buf={n:[] for n in cols}
for line in open(NSp/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame(buf)
for c in N: X[c]=np.where(X[c].between(1,5),X[c],np.nan)
X["年龄段"]=pd.cut(X.ager,[14,24,34,45],labels=["15-24","25-34","35-45"])
X["教育"]=pd.cut(X.hieduc,[4,9,10,12,15],labels=["<高中","高中","部分大学","大学+"])
ov_n,_=weakest(X,N)
print(f"\n=== 正对照:NSFG 总体最弱一环 = **{ov_n:+.4f}**(`#653`: 0.346)")
for s in ["年龄段","教育"]:
    for lv in sorted(X[s].dropna().unique(), key=str):
        sub=X[X[s]==lv]
        if len(sub)<FLOOR: continue
        w,_=weakest(sub,N)
        rows.append(dict(inst="NSFG",strat=s,level=str(lv),weakest=float(w) if np.isfinite(w) else None,
                         placebo=None,age=float(sub.ager.median()),educ=float(sub.hieduc.median())))
        print(f"  {s:5s}{str(lv):10s} n={len(sub):5d} 最弱 {w:+.4f} · 年龄中位 {sub.ager.median():.0f}")

R=pd.DataFrame(rows)
print("\n=== 主量:层间梯度(跨仪器只用 年龄 / 教育)===")
grad={}
for inst in ["GSS","NSFG"]:
    for s in ["年龄段","教育"]:
        t=R[(R.inst==inst)&(R.strat==s)&R.weakest.notna()]
        if len(t)<2: continue
        hi,lo=t.loc[t.weakest.idxmax()],t.loc[t.weakest.idxmin()]
        grad[(inst,s)]=dict(g=float(hi.weakest-lo.weakest),hi=hi.level,lo=lo.level)
        print(f"  {inst:5s}{s:5s} 梯度 **{hi.weakest-lo.weakest:+.4f}**({hi.level} {hi.weakest:+.4f} vs {lo.level} {lo.weakest:+.4f})")
print("\n  跨仪器同号?")
for s in ["年龄段","教育"]:
    a=grad.get(("GSS",s)); b=grad.get(("NSFG",s))
    if a and b:
        same = (a["g"]>0)==(b["g"]>0) and a["hi"]==b["hi"]
        print(f"    {s}: GSS {a['g']:+.4f}(高={a['hi']}) · NSFG {b['g']:+.4f}(高={b['hi']}) -> "
              f"**{'同号且同一层最高' if same else '不一致'}**")
json.dump(dict(rows=rows,grad={f"{k[0]}·{k[1]}":v for k,v in grad.items()},
               overall_gss=float(ov_g),overall_nsfg=float(ov_n),unchallenged=True),
          open(OUT/"one_block_for_everyone.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'one_block_for_everyone.json'}")
