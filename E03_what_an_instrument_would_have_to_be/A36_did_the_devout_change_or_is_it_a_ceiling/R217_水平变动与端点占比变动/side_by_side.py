"""#778 · E03·A43·R217 —— 「掉了 18.5 个百分点」相对谁?两条已发表的话从没并排读过

`#777` 报:最虔诚层里说同性恋「总是错」的从 **85.9% 掉到 67.4%**(2.38× 自己的零)。
`#747` 报:**常去教堂的人在同性恋态度上的水平变动只有别人的一半**(+0.640 vs +1.187,比值 **0.539**)。
⚠ **这两条是同一个对象的两种切法,而我从没把它们放在一起读** ——
   `#777` 的「真的改了主意」**没有配对照**:掉 18.5 点,相对谁?

G1 估计量:同一批年份、同一分层下,**两个统计量各自的「虔诚/非虔诚」比**:
  ① **水平**:该层 `homosex` 均值对年份的斜率(`#747` 用的量)
  ② **端点占比**:该层说「总是错」的占比对年份的斜率(`#777` 用的量)
  ③ ⚠ **相对基线的端点跌幅**:斜率 ÷ 首年占比 —— 因为**基线不同,绝对点数不可比**

⚠ **两条混淆写在跑之前:**
① **基线不同 ⇒ 绝对点数不可比。** 非虔诚层的端点占比起点低得多,**机械上就掉不了那么多点** ——
   所以 ③ 必须与 ② 同报,而不是留给读者自己换算。
② **均值在有界量表上同样被端点挤压**(`#776` 那条算术)⇒ **两个统计量都带着有界性**,
   **这不是「一个干净一个脏」**,而是两种不同的挤压方式。

⚠ 分层两种都跑(G4):**(a) `attend` 三档**(`#747` 的原分层)· **(b) 三题虔诚度合成三分位**(`#776`/`#777` 的)。
⚠ **换不了仪器,与 `#776`/`#777` 同一条且同样量过**:需要同一道题在 ≥20 年上的重复测量;
   MFQ 单次采集无年代 · NSFG 的 `SXOK18` 在 2017–19 卷字典里不存在 · SCCS 的单位是社会。**只此一具。**

预注册判词(按 `#764` 新写法:只比已测量的量,各带自己的零):
  **两个统计量的「虔诚/非虔诚」比相差 ≤2× ⇒ A(一致),而 `#777` 的措辞必须挂上比较;
  相差 >2× ⇒ B,统计量的选择本身是发现,不选边、报整张网格。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(217)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
VALID={"homosex":(1,4),"attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=["homosex"],convert_categoricals=True)
cats={"homosex":list(cat["homosex"].cat.categories)[:4]}
# ⚠ 统一成「高=严」并留在 1..4,端点 4 = 「总是错」
for c in aligned(cats,"strict"): M[c]=-M[c]+5
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
def fit(yr,v):
    b=slope(yr,v); nul=[slope(RNG.permutation(yr),v) for _ in range(3000)]
    return b, float(np.quantile(np.abs(nul),.95))

CUTS={}
sub=M.dropna(subset=["homosex","attend","reliten","fund","year"]).copy()
sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
# (a) #747 的原分层:attend 三档
sub["a"]=pd.cut(sub.attend,[-1,1,5,8],labels=[0,1,2]).astype(float)
# (b) #776/#777 的:三题合成的逐年三分位
sub["b"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
CUTS["(a) attend 三档(#747 原分层)"]="a"; CUTS["(b) 三题合成三分位(#776/#777)"]="b"

print("\n=== 硬规则①:两种分层各层的年数与 n ===")
res={}
for name,col in CUTS.items():
    res[name]={}
    for k,lab in ((2,"虔诚"),(0,"非虔诚")):
        g=sub[sub[col]==k]; rows=[]
        for y,gy in g.groupby("year"):
            if len(gy)<120: continue
            rows.append((int(y),float(gy.homosex.mean()),float((gy.homosex==4).mean()),len(gy)))
        if len(rows)<8: continue
        yr=np.array([r[0] for r in rows],float)
        b_lv,q_lv=fit(yr,np.array([r[1] for r in rows]))
        b_ep,q_ep=fit(yr,np.array([r[2] for r in rows]))
        res[name][lab]=dict(n_years=len(rows),n_med=int(np.median([r[3] for r in rows])),
                            lv_first=rows[0][1],lv_last=rows[-1][1],lv_slope=b_lv,lv_null=q_lv,
                            ep_first=rows[0][2],ep_last=rows[-1][2],ep_slope=b_ep,ep_null=q_ep,
                            ep_rel=b_ep/rows[0][2])
        r=res[name][lab]
        print(f"  {name[:22]:24s}{lab:5s} {len(rows):2d} 年 · 每年中位 n={r['n_med']:4d} · "
              f"水平 {r['lv_first']:.3f}→{r['lv_last']:.3f} · 端点 {r['ep_first']*100:.1f}%→{r['ep_last']*100:.1f}%")

print(f"\n=== G3 全网格:两个统计量各自的斜率、零、以及「虔诚/非虔诚」比 ===")
print(f"  {'分层':26s}{'层':6s}{'水平斜率/年':>13s}{'零95%':>9s}{'倍':>6s}"
      f"{'端点斜率/年':>13s}{'零95%':>9s}{'倍':>6s}{'相对基线':>10s}")
ratios={}
for name in res:
    for lab in ("虔诚","非虔诚"):
        if lab not in res[name]: continue
        r=res[name][lab]
        print(f"  {name[:24]:26s}{lab:6s}{r['lv_slope']:+13.6f}{r['lv_null']:9.6f}{abs(r['lv_slope'])/r['lv_null']:5.2f}×"
              f"{r['ep_slope']:+13.6f}{r['ep_null']:9.6f}{abs(r['ep_slope'])/r['ep_null']:5.2f}×{r['ep_rel']:+10.5f}")
    if "虔诚" in res[name] and "非虔诚" in res[name]:
        a,b=res[name]["虔诚"],res[name]["非虔诚"]
        ratios[name]=dict(level=a["lv_slope"]/b["lv_slope"], endpoint=a["ep_slope"]/b["ep_slope"],
                          endpoint_rel=a["ep_rel"]/b["ep_rel"])
        print(f"    -> 虔诚/非虔诚 比:**水平 {ratios[name]['level']:.3f}** · "
              f"端点(绝对){ratios[name]['endpoint']:.3f} · **端点(相对基线){ratios[name]['endpoint_rel']:.3f}**")

G=Gate("#778 · 两个统计量的比")
nm=list(ratios)[0]
G.asserted("① 两层的年数必须相同(否则两个比不可比)",
           bool(res[nm]["虔诚"]["n_years"]==res[nm]["非虔诚"]["n_years"]),
           f"虔诚 {res[nm]['虔诚']['n_years']} 年 vs 非虔诚 {res[nm]['非虔诚']['n_years']} 年", kind="control")
G.offset_control("② 虔诚层的端点斜率须超出它自己的零(否则这个比的分子是噪声)",
                 effect=abs(res[nm]["虔诚"]["ep_slope"]), offset=res[nm]["虔诚"]["ep_null"],
                 spread=res[nm]["虔诚"]["ep_null"]*0.1,
                 null_kind="打乱年份标签与该层端点占比的配对,保住每年的 n,只毁掉「哪一年配哪个占比」")
G.offset_control("③ 非虔诚层同样须超零(否则这个比的分母是噪声)",
                 effect=abs(res[nm]["非虔诚"]["ep_slope"]), offset=res[nm]["非虔诚"]["ep_null"],
                 spread=res[nm]["非虔诚"]["ep_null"]*0.1,
                 null_kind="同上,换成非虔诚层")
print(); print(G)
print("\n"+"="*94)
if not all(r[2] for r in G.rows): v="**UNVERIFIED:闸没全过**"
else:
    lv=np.median([r["level"] for r in ratios.values()]); ep=np.median([r["endpoint_rel"] for r in ratios.values()])
    fac=max(lv,ep)/min(lv,ep) if min(lv,ep)>0 else float("inf")
    if fac<=2.0:
        v=(f"**A 一致:水平比 {lv:.3f} vs 端点(相对基线)比 {ep:.3f},相差 {fac:.2f}× ≤2 "
           f"⇒ 两条已发表的话说的是同一件事,而 `#777` 的「真的改了主意」必须挂上比较:"
           f"**他们改了,但只有别人的 {lv:.0%}–{ep:.0%}**")
    else:
        v=(f"**B 不一致:水平比 {lv:.3f} vs 端点(相对基线)比 {ep:.3f},相差 {fac:.2f}× >2 "
           f"⇒ 统计量的选择本身是发现,不选边,整张网格已全列**")
print(v)
json.dump(dict(res=res,ratios=ratios,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"side_by_side.json","w"),ensure_ascii=False,indent=1)
