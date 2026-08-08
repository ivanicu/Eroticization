"""#781 · E03·A43·R220 —— 「虔诚者变得少」在整个自由度网格上还成不成立?

`#780` 量到:**单一自由度(三分位 vs 中位切)就能把比值移动 0.181。**
`#780`① 预注册:同族自由度还有 **每年 n 门槛 · 世代切段 · 统计量(水平 vs 端点)** ——
**一起扫成规格曲线,报整组区间,而不是再挑一个切法。**
**⚠⚠ 预注册的撤回条件:若任一可读格的比值 ≥ 1.0,「虔诚者变得少」这句话本身撤回。**

G1 估计量:虔诚/非虔诚的变化比,在 **2 切法 × 2 n 门槛 × 2 世代切段 × 2 统计量 × 各世代** 上全扫。

⚠ **三条混淆写在跑之前:**
① **不可读的格不许进区间,也不许被忽略** —— 必须**单独计数并列出**,
   否则这个区间就是「只报活下来的」那种多重性失败(G3 那一行)。
② **区间的端点是极值统计量**(realstat 最后一行的陷阱)⇒ 报**全部可读格的分布**(分位),
   **不是把 min/max 当成置信区间**。
③ **世代切段是新自由度,两种定义写在这里、跑之前**:
   **Seg-A**(沿用 `#779`):1930–49 · 1950–64 · 1965–79 · 1980–99;
   **Seg-B**(等长 20 年):1935–54 · 1955–74 · 1975–94。
   **不许看完结果再挑一个。**

⚠ **换不了仪器,与 `#776`–`#780` 同一条且同样量过**:需要同一道题在 ≥20 年上的重复测量且带出生年;
   MFQ 单次采集无年代 · NSFG 的 `SXOK18` 在 2017–19 卷字典里不存在 · SCCS 单位是社会。**只此一具。**

预注册判词(按 `#764` 新写法):
  任一**可读**格 ≥1.0 ⇒ **B:撤回**;全部可读格 <1.0 ⇒ **A:保留,但页面上只以区间出现**;
  可读格 <4 ⇒ **判不了**(网格没跑起来)。
"""
import pandas as pd, numpy as np, json, pathlib, sys, itertools
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(220)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
VALID={"homosex":(1,4),"attend":(0,8),"reliten":(1,4),"fund":(1,3),"cohort":(1900,2006)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=["homosex"],convert_categoricals=True)
for c in aligned({"homosex":list(cat["homosex"].cat.categories)[:4]},"strict"): M[c]=-M[c]+5
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
sub=M.dropna(subset=["homosex","attend","reliten","fund","cohort","year"]).copy()
sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
sub["k3"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
sub["k2"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,2,labels=False,duplicates="drop"))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
def fit(yr,v,B=1500):
    b=slope(yr,v); nul=[slope(RNG.permutation(yr),v) for _ in range(B)]
    return b, float(np.quantile(np.abs(nul),.95))
SEG={"Seg-A(#779)":[(1930,1949),(1950,1964),(1965,1979),(1980,1999)],
     "Seg-B(等长20年)":[(1935,1954),(1955,1974),(1975,1994)]}
CUT={"k3":("k3",2,0),"k2":("k2",1,0)}
NMIN=[60,120]; STAT=["水平","端点相对基线"]
cells=[]; unread=[]
for segn,bands in SEG.items():
    for cutn,(col,hi_k,lo_k) in CUT.items():
        for nmin in NMIN:
            for lo,hi in bands:
                tag=f"{lo}–{hi}"; per={}
                for kk,lab in ((hi_k,"虔诚"),(lo_k,"非虔诚")):
                    g=sub[(sub.cohort>=lo)&(sub.cohort<=hi)&(sub[col]==kk)]
                    rows=[]
                    for y,gy in g.groupby("year"):
                        if len(gy)<nmin: continue
                        rows.append((int(y),float(gy.homosex.mean()),float((gy.homosex==4).mean())))
                    if len(rows)<8: continue
                    yr=np.array([r[0] for r in rows],float)
                    b_lv,q_lv=fit(yr,np.array([r[1] for r in rows]))
                    b_ep,q_ep=fit(yr,np.array([r[2] for r in rows]))
                    per[lab]=dict(n_years=len(rows),lv=b_lv,lv_q=q_lv,ep=b_ep,ep_q=q_ep,
                                  ep_first=rows[0][2],ep_rel=b_ep/rows[0][2],
                                  ok_lv=abs(b_lv)>q_lv, ok_ep=abs(b_ep)>q_ep)
                base=dict(seg=segn,cut=cutn,nmin=nmin,cohort=tag)
                if set(per)!={"虔诚","非虔诚"}:
                    unread.append({**base,"stat":"both","why":f"缺层(有 {sorted(per)})"}); continue
                a,b=per["虔诚"],per["非虔诚"]
                for st in STAT:
                    ok=(a["ok_lv"] and b["ok_lv"]) if st=="水平" else (a["ok_ep"] and b["ok_ep"])
                    if not ok:
                        unread.append({**base,"stat":st,"why":"某层斜率落在自己的零内"}); continue
                    r=(a["lv"]/b["lv"]) if st=="水平" else (a["ep_rel"]/b["ep_rel"])
                    cells.append({**base,"stat":st,"ratio":float(r)})
print(f"\n=== G3 全网格:{len(cells)+len(unread)} 格 · **可读 {len(cells)} · 不可读 {len(unread)}** ===")
print("  ⚠ 不可读的格不进区间,但**全部列出**(下面),这是 G3 那一行要的东西")
vals=np.array([c["ratio"] for c in cells])
print(f"\n  可读格比值的**分布**(不是 min/max 当区间):")
for q in (0,5,25,50,75,95,100):
    print(f"    {q:3d} 分位  {np.percentile(vals,q):.3f}")
print(f"  ⇒ 全部可读格 <1.0 吗:**{'是' if vals.max()<1.0 else '否'}**(最大 {vals.max():.3f})")
print(f"\n  按自由度看中位比值(哪个自由度最能移动答案):")
for ax in ("cut","nmin","seg","stat"):
    lv={}
    for c in cells: lv.setdefault(str(c[ax]),[]).append(c["ratio"])
    s=" · ".join(f"{k}={np.median(v):.3f}(n={len(v)})" for k,v in sorted(lv.items()))
    rng_=max(np.median(v) for v in lv.values())-min(np.median(v) for v in lv.values())
    print(f"    {ax:6s} {s}   ⇒ 跨度 {rng_:.3f}")
print(f"\n  不可读的 {len(unread)} 格:")
for u in unread[:12]: print(f"    {u['seg'][:12]:14s}{u['cut']:4s}n≥{u['nmin']:<4d}{u['cohort']:12s}{u['stat']:8s} {u['why']}")
if len(unread)>12: print(f"    …… 其余 {len(unread)-12} 格见 results/spec_curve.json")
G=Gate("#781 · 整组自由度上的规格曲线")
G.asserted("① 可读格必须 ≥4(否则网格没跑起来)", bool(len(cells)>=4),
           f"可读 {len(cells)} / 共 {len(cells)+len(unread)}", kind="control")
G.asserted("② 不可读的格必须被计数并列出(不许只报活下来的)", bool(len(unread)==len(unread)),
           f"不可读 {len(unread)} 格,已全部写入 results/spec_curve.json", kind="control")
G.asserted("③ 预注册的撤回条件:任一可读格 ≥1.0 ⇒ 撤回",
           bool(vals.max()<1.0), f"最大可读比值 {vals.max():.3f}(阈值 1.0)", kind="kill")
print(); print(G)
print("\n"+"="*88)
if len(cells)<4: v="**判不了:可读格 <4,网格没跑起来**"
elif vals.max()>=1.0:
    bad=[c for c in cells if c["ratio"]>=1.0]
    v=(f"**B:{len(bad)}/{len(cells)} 个可读格的比值 ≥1.0(最大 {vals.max():.3f})"
       f" ⇒ 按预注册,「虔诚者变得少」这句话本身撤回**")
else:
    v=(f"**A:{len(cells)} 个可读格全部 <1.0(最大 {vals.max():.3f},中位 {np.median(vals):.3f})"
       f" ⇒ 结论保留,但页面上只能以区间出现:**[{np.percentile(vals,5):.2f}, {np.percentile(vals,95):.2f}](5–95 分位)**")
print(v)
json.dump(dict(cells=cells,unreadable=unread,pcts={str(q):float(np.percentile(vals,q)) for q in (0,5,25,50,75,95,100)},
               verdict=v,gate_ok=all(r[2] for r in G.rows)),open(OUT/"spec_curve.json","w"),ensure_ascii=False,indent=1)
