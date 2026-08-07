"""#780 · E03·A43·R219 —— 换两层切法能让年轻世代可读吗?而切法本身会不会在决定答案?

`#779` 只有 **2/4** 世代可读(1930–49 · 1950–64),刚好卡在预注册下限;
1965–79 的虔诚层**每年中位 n 仅 97、年数 15**,斜率落在自己的零内(0.54–0.64×)。
`#779`① 预注册:**两层(中位切)会让 n 翻倍 —— 够不够让它进可读集合?**

⚠⚠ **而这一步天然长着「换个更宽松的切法去凑显著」的样子。护栏是预注册的顺序,不是我的自觉:**
**先看老世代的比值移动量,再看年轻世代能不能读。若老世代移动 >0.1,不许读年轻世代的结果。**

⚠ **第二条混淆,同样写在跑之前**:两层把**中间三分之一并进两侧** ⇒ **对比度天然下降**
⇒ **比值天然向 1 靠**。**这是算术不是发现** —— 所以必须在**老世代**上先把这个基线量出来,
再拿它当尺子读年轻世代。

G1 估计量:同 `#779`,但分层是 **k=2(中位切)** 与 **k=3(三分位)** 两种;
每格报斜率、它自己的零 95% 分位、以及可读性;老世代额外报**两种切法下比值的移动量**。

⚠ **换不了仪器,与 `#776`–`#779` 同一条且同样量过**:需要同一道题在 ≥20 年上的重复测量且带出生年;
   MFQ 单次采集无年代 · NSFG 的 `SXOK18` 在 2017–19 卷字典里不存在 · SCCS 单位是社会。**只此一具。**

预注册判词(按 `#764` 新写法):
  ① 老世代(1930–49 · 1950–64)在两种切法下的比值移动 **>0.1 ⇒ B**:切法在决定答案,
     整组比值改报成「按切法的区间」,**且不许读年轻世代**;
  ② 移动 **≤0.1** 且 1965–79 的**两层都超各自的零 ⇒ A**;
  ③ 移动 ≤0.1 但 1965–79 仍不超零 ⇒ **C:能力边界**(两层也救不回)。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(219)
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
def fit(yr,v,B=3000):
    b=slope(yr,v); nul=[slope(RNG.permutation(yr),v) for _ in range(B)]
    return b, float(np.quantile(np.abs(nul),.95))
COH=[(1930,1949),(1950,1964),(1965,1979),(1980,1999)]
CUTS={"k3 三分位(#779)":("k3",2,0),"k2 中位切(本轮)":("k2",1,0)}
res={}
print(f"\n=== 逐切法 · 逐世代 · 逐层(每年 n≥60,与 `#779` 同门槛)===")
print(f"  {'切法':18s}{'世代':12s}{'层':6s}{'年数':>5s}{'中位n':>6s}{'端点斜率':>11s}{'零95%':>9s}{'倍':>6s}")
for cname,(col,hi_k,lo_k) in CUTS.items():
    res[cname]={}
    for lo,hi in COH:
        tag=f"{lo}–{hi}"; res[cname][tag]={}
        for kk,lab in ((hi_k,"虔诚"),(lo_k,"非虔诚")):
            g=sub[(sub.cohort>=lo)&(sub.cohort<=hi)&(sub[col]==kk)]
            rows=[]
            for y,gy in g.groupby("year"):
                if len(gy)<60: continue
                rows.append((int(y),float(gy.homosex.mean()),float((gy.homosex==4).mean()),len(gy)))
            if len(rows)<8: continue
            yr=np.array([r[0] for r in rows],float)
            b_lv,q_lv=fit(yr,np.array([r[1] for r in rows]))
            b_ep,q_ep=fit(yr,np.array([r[2] for r in rows]))
            res[cname][tag][lab]=dict(n_years=len(rows),n_med=int(np.median([r[3] for r in rows])),
                lv_slope=b_lv,lv_null=q_lv,ep_slope=b_ep,ep_null=q_ep,ep_first=rows[0][2],ep_rel=b_ep/rows[0][2],
                readable=bool(abs(b_ep)>q_ep and abs(b_lv)>q_lv))
            r=res[cname][tag][lab]
            print(f"  {cname[:16] if (tag==f'{COH[0][0]}–{COH[0][1]}' and lab=='虔诚') else '':18s}"
                  f"{tag if lab=='虔诚' else '':12s}{lab:6s}{r['n_years']:5d}{r['n_med']:6d}"
                  f"{r['ep_slope']:+11.6f}{r['ep_null']:9.6f}{abs(r['ep_slope'])/r['ep_null']:5.2f}×")
def ratio(cname,tag):
    c=res[cname].get(tag,{})
    if set(c)!={"虔诚","非虔诚"} or not (c["虔诚"]["readable"] and c["非虔诚"]["readable"]): return None
    return dict(level=c["虔诚"]["lv_slope"]/c["非虔诚"]["lv_slope"],
                endpoint_rel=c["虔诚"]["ep_rel"]/c["非虔诚"]["ep_rel"])
print(f"\n=== ① 先看老世代:两种切法下比值移动了多少(护栏在此,不是在年轻世代)===")
moves=[]
for tag in ("1930–1949","1950–1964"):
    r3,r2=ratio("k3 三分位(#779)",tag),ratio("k2 中位切(本轮)",tag)
    if r3 and r2:
        dl=abs(r2["level"]-r3["level"]); de=abs(r2["endpoint_rel"]-r3["endpoint_rel"])
        moves += [dl,de]
        print(f"  {tag:12s} 水平 {r3['level']:.3f} → {r2['level']:.3f}(移动 {dl:.3f})· "
              f"端点 {r3['endpoint_rel']:.3f} → {r2['endpoint_rel']:.3f}(移动 {de:.3f})")
    else: print(f"  {tag:12s} 某切法下不可读")
mx=max(moves) if moves else float("inf")
print(f"  ⇒ **最大移动 {mx:.3f}(预注册阈值 0.1)**")
print(f"\n=== ② 年轻世代在两层下可读了吗 ===")
young={}
for tag in ("1965–1979","1980–1999"):
    for cname in CUTS:
        c=res[cname].get(tag,{})
        st=[f"{lab}:{'可读' if v['readable'] else f'不可读({abs(v[chr(101)+chr(112)+chr(95)+chr(115)+chr(108)+chr(111)+chr(112)+chr(101)])/v[chr(101)+chr(112)+chr(95)+chr(110)+chr(117)+chr(108)+chr(108)]:.2f}×)'}"
              for lab,v in c.items()] or ["缺层"]
        young[(tag,cname)]=all(v["readable"] for v in c.values()) and len(c)==2
        print(f"  {tag:12s}{cname:18s} {' · '.join(st)}")
G=Gate("#780 · 切法会不会在决定答案")
G.asserted("① 护栏:老世代在两种切法下的比值移动必须 ≤0.1(否则不许读年轻世代)",
           bool(mx<=0.10), f"最大移动 {mx:.3f}(阈值 0.10)", kind="kill")
t0="1930–1949"
if ratio("k2 中位切(本轮)",t0):
    G.offset_control("② k2 下该世代虔诚层端点斜率仍须超它自己的零",
                     effect=abs(res["k2 中位切(本轮)"][t0]["虔诚"]["ep_slope"]),
                     offset=res["k2 中位切(本轮)"][t0]["虔诚"]["ep_null"],
                     spread=res["k2 中位切(本轮)"][t0]["虔诚"]["ep_null"]*0.1,
                     null_kind="在该世代该层内打乱年份标签与端点占比的配对,保住每年的 n,只毁掉「哪一年配哪个占比」")
print(); print(G)
print("\n"+"="*96)
if mx>0.10:
    v=(f"**B:老世代比值在两种切法下最大移动 {mx:.3f} > 0.10 ⇒ **切法在决定答案** —— "
       f"整组比值改报成「按切法的区间」,而年轻世代的结果**按预注册不读**")
elif young[("1965–1979","k2 中位切(本轮)")]:
    r=ratio("k2 中位切(本轮)","1965–1979")
    v=(f"**A:老世代移动 {mx:.3f} ≤0.10,而 1965–1979 在两层下可读 ⇒ 水平比 {r['level']:.3f} · "
       f"端点(相对基线)比 {r['endpoint_rel']:.3f}**")
else:
    v=(f"**C:老世代移动 {mx:.3f} ≤0.10,但 1965–1979 在两层下仍不可读 ⇒ 能力边界 —— "
       f"两层也救不回,而这不是切法的问题,是那个世代的年份跨度与 n**")
print(v)
json.dump(dict(res=res,max_move=mx,young=({f"{k[0]}|{k[1]}":vv for k,vv in young.items()}),
               verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"cut_or_answer.json","w"),ensure_ascii=False,indent=1)
