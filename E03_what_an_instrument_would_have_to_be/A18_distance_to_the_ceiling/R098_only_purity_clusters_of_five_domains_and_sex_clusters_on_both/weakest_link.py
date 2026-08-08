"""E03·A18·R689 —— 那十一个 1.00,是团还是链

**类型:FRONTIER**。`#652` 自己指出:连通子块只要求**每题至少连到一题**,这是一个**非常弱的条件** ——
一条链 `a—b—c—d` 也给 1.00,而它和一个真正的**团**(每对都达阈)是两回事。

⚠ **判据可以化简成一个不含阈的量**:「边密度 = 1.00 所能容忍的最大阈」**就是块内最小的那一对相关**。
  ⇒ **估计量 = 最弱一环(min pairwise normalised r)**,而 `#652` 要求的「每块的最大阈」正是它。
  **这不是偏离,是把预注册的量写成它的闭式。**

W1 **十一块都是团** —— 最弱一环 >= 0.30 ⇒ 「一个维度」这个说法在它们身上成立。
W2 **有些只是链** —— 最弱一环 < 0.30 ⇒ **它只是连通,不是一块**,凡把它当「一个维度」用过的声明要重看。
⚠ **BASIN**:上两轮我连输两次(下注孤例、下注天花板伪影),**而连输本身也是一种盆地** ——
  所以本轮**不下注**,只把两个世界的判据都写死,并且**先写下哪一块我最不希望它是链**:
  **`MFQ·PURITY`** —— 这一页有五条声明立在「纯洁是一个域」上。**若它是链,那五条都要重看。**

G1 ESTIMAND:每块的 **min 归一对相关**(最弱一环)+ **0.30 处的边密度**(达阈对数 ÷ 全部对数)。
G2 CONTROLS:**正对照** `SCCS·以身作则四对象`(`#642` 六对全 >= +0.85)最弱一环必须 >= 0.85(生);
  **安慰剂** 打乱行后最弱一环 -> 约 0 且边密度 -> 0。
G3:14 块全报。G4:边密度在 4 阈上各报一次。
KILL(条件式):if 正对照 >= 0.85 and 安慰剂边密度 ≈0:
  逐块按「最弱一环 >= 0.30」判团/链;else UNVERIFIED
IMPOSSIBLE(不写 planned):同 `#652` · `[unchallenged]`
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

def linkstats(name):
    F,items,who,yc=BLK[name]
    A,_=align(F,items); Fa=A.join(F[[yc]]) if yc else A
    pm=pairmap(Fa,items,year=yc,floor=FLOOR[name])
    raws=[v[0] for v in pm.values()]; norms=[v[1] for v in pm.values()]
    dens={t: sum(1 for v in pm.values() if v[1]>=t)/len(pm) for t in THR}
    return dict(min_raw=float(min(raws)),min_norm=float(min(norms)),
                max_norm=float(max(norms)),dens={f"{t:.2f}":dens[t] for t in THR},
                npairs=len(pm),weakest=min(pm,key=lambda k:pm[k][1]))

print("\n=== G3/G4:14 块 —— 最弱一环 与 边密度 ===")
print(f"  {'块':18s}{'名义':>4s}{'最弱(生)':>10s}{'最弱(归一)':>11s}{'最强(归一)':>11s}   边密度@.20/.25/.30/.35   最弱的那一对")
L={}
for name in ROWS:
    if ROWS[name]["undecidable"]: print(f"  {name:18s} 判不了"); continue
    s=linkstats(name); L[name]=s
    dd="/".join(f"{s['dens'][f'{t:.2f}']:.2f}" for t in THR)
    print(f"  {name:18s}{ROWS[name]['nominal']:>4d}{s['min_raw']:>10.3f}{s['min_norm']:>11.3f}{s['max_norm']:>11.3f}        {dd:19s}  {s['weakest'][0]}×{s['weakest'][1]}")

print("\n=== 控制 ===")
pc=L["SCCS·以身作则四对象"]["min_raw"]
print(f"  正对照 SCCS·以身作则最弱一环(生)必须 >= 0.85 -> **{pc:.4f}**")
def placebo(seed):
    rng=np.random.default_rng(seed); F,items,_,_=BLK["MFQ·PURITY"]
    Z=F[items].dropna().copy()
    for c in items: Z[c]=rng.permutation(Z[c].to_numpy())
    A,_=align(Z,items); pm=pairmap(A,items,floor=150)
    return sum(1 for v in pm.values() if v[1]>=0.30)/len(pm)
pd_=float(np.median([placebo(s) for s in SEEDS]))
print(f"  安慰剂 打乱行后 0.30 处边密度必须 -> 0 -> **{pd_:.3f}**")

G=Gate("那十一个 1.00,是团还是链")
p1=G.positive_control("SCCS·以身作则最弱一环(生) >= 0.85",planted=pc,floor=0.85,spread=0.005)
p2=G.negative_control("安慰剂:打乱行后边密度 -> 0",null=pd_,effect=float(np.median([s["dens"]["0.30"] for s in L.values()])),
                      null_spread=0.02,null_kind="行内打乱,保留每题边际")
if p1 and p2:
    clique=[k for k,s in L.items() if s["min_norm"]>=0.30]
    chain=[k for k,s in L.items() if s["min_norm"]<0.30]
    verdict=(f"**{len(clique)}/{len(L)} 块是团(最弱一环 >= 0.30)· {len(chain)} 块只是连通** -> 链:{chain}")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
pur=L["MFQ·PURITY"]
print(f"\n=== 事先点名的那一块:MFQ·PURITY 最弱一环(归一)= **{pur['min_norm']:.4f}** "
      f"({pur['weakest'][0]}×{pur['weakest'][1]}) -> **{'团' if pur['min_norm']>=0.30 else '链 ⇒ 五条声明要重看'}**")
json.dump(dict(links=L,positive=pc,placebo=pd_,verdict=verdict,
               purity_weakest=[pur["weakest"][0],pur["weakest"][1],pur["min_norm"]],unchallenged=True),
          open(OUT/"weakest_link.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'weakest_link.json'}")
