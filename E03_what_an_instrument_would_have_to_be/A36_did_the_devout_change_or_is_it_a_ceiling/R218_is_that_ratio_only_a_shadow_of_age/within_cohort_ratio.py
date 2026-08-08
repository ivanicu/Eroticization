"""#779 · E03·A43·R218 —— 「虔诚者变得少」会不会只是「老人变得少」?

`#747`/`#777`/`#778` 三条都说虔诚者在同性恋态度上变得少(虔诚/非虔诚比 0.33–0.43)。
**⚠ 而虔诚者更老,老人本来就变得少 —— 这条弧从没做过世代控制。**
`#776`·`#777`·`#778` 连着三轮在同一方向上 ⇒ **basin**。本轮是那个**正面结果我不欢迎**的步:
  **A 宗教是真的轴**:固定出生世代内,比值仍 ≈0.33–0.43;
  **B 年龄/世代的影子**:固定世代后比值**向 1 靠** ⇒ **三条已发表的行都要收窄成「老人变得少」。**

G1 估计量:在**固定出生世代**内,虔诚层与非虔诚层各自的
  ① `homosex` 均值对年份的斜率 · ② 说「总是错」的占比对年份的斜率,
以及两者的**虔诚/非虔诚比**;与 `#778` 的合并比(水平 0.409–0.431 · 端点相对基线 0.328–0.341)并排。

⚠ **三条混淆写在跑之前:**
① ⚠⚠ **APC 不可识别**:世代内 `age = year − cohort` 是**恒等式** ⇒ **世代内的斜率是「时期+年龄」的合成**,
   **不能声称分离了年龄**。⇒ **判词只能说「固定世代后那个比动不动」**,不能说「所以不是年龄」。
② **世代切细 ⇒ 每格 n 掉、零变宽。先算 MDE** —— 每格报它自己的零 95% 分位与年数;
   若两层的斜率都落在各自零内,**该世代不可读**,如实标注,不许平均进去。
③ **虔诚度本身随年龄变**(老人更虔诚)⇒ 世代内分层仍是**同期分位**,
   **不是「同一个人一生的虔诚度」**,如实记。

⚠ **换不了仪器,与 `#776`/`#777`/`#778` 同一条且同样量过**:需要同一道题在 ≥20 年上的重复测量
   **且带出生年**;MFQ 单次采集无年代 · NSFG 的 `SXOK18` 在 2017–19 卷字典里不存在 · SCCS 单位是社会。**只此一具。**
⚠ **打印小瑕疵(留给读者)**:世代标签只在该世代第一行打印,所以 1980–1999 的「非虔诚」行看起来像挂在
   1965–1979 下面。**数据不受影响**(JSON 按世代键存),但表读起来会误导 —— 记在这里而不是悄悄改。

预注册判词(按 `#764` 新写法:只比已测量的量,各带自己的零):
  **可读世代的比值中位 ≤0.7 ⇒ A**(宗教这条轴在固定世代内仍在);
  **>0.7 ⇒ B**(向 1 靠,三条行要收窄);
  **可读世代少于 2 个 ⇒ 判不了**,不许拿一个世代下结论。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(218)
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
sub["k"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
def fit(yr,v,B=3000):
    b=slope(yr,v); nul=[slope(RNG.permutation(yr),v) for _ in range(B)]
    return b, float(np.quantile(np.abs(nul),.95))

COH=[(1930,1949),(1950,1964),(1965,1979),(1980,1999)]
print(f"\n=== 逐世代 · 逐层(每年 n≥60;⚠ 比 `#778` 的 120 低,因为世代切细,如实记)===")
print(f"  {'世代':12s}{'层':6s}{'年数':>5s}{'中位n':>6s}{'水平斜率':>11s}{'零95%':>9s}{'倍':>6s}"
      f"{'端点斜率':>11s}{'零95%':>9s}{'倍':>6s}")
res={}
for lo,hi in COH:
    tag=f"{lo}–{hi}"; res[tag]={}
    for kk,lab in ((2,"虔诚"),(0,"非虔诚")):
        g=sub[(sub.cohort>=lo)&(sub.cohort<=hi)&(sub.k==kk)]
        rows=[]
        for y,gy in g.groupby("year"):
            if len(gy)<60: continue
            rows.append((int(y),float(gy.homosex.mean()),float((gy.homosex==4).mean()),len(gy)))
        if len(rows)<8: continue
        yr=np.array([r[0] for r in rows],float)
        b_lv,q_lv=fit(yr,np.array([r[1] for r in rows]))
        b_ep,q_ep=fit(yr,np.array([r[2] for r in rows]))
        res[tag][lab]=dict(n_years=len(rows),n_med=int(np.median([r[3] for r in rows])),
                           lv_slope=b_lv,lv_null=q_lv,ep_slope=b_ep,ep_null=q_ep,
                           ep_first=rows[0][2],ep_rel=b_ep/rows[0][2])
        r=res[tag][lab]
        print(f"  {tag if lab=='虔诚' else '':12s}{lab:6s}{r['n_years']:5d}{r['n_med']:6d}"
              f"{r['lv_slope']:+11.6f}{r['lv_null']:9.6f}{abs(r['lv_slope'])/r['lv_null']:5.2f}×"
              f"{r['ep_slope']:+11.6f}{r['ep_null']:9.6f}{abs(r['ep_slope'])/r['ep_null']:5.2f}×")

print(f"\n=== 比值(仅在**两层都超各自零**的世代上计算;其余标为不可读)===")
ratios={}; unread=[]
for tag,cells in res.items():
    if set(cells)!={"虔诚","非虔诚"}: unread.append((tag,"缺层")); continue
    a,b=cells["虔诚"],cells["非虔诚"]
    ok=(abs(a["ep_slope"])>a["ep_null"]) and (abs(b["ep_slope"])>b["ep_null"]) \
       and (abs(a["lv_slope"])>a["lv_null"]) and (abs(b["lv_slope"])>b["lv_null"])
    if not ok: unread.append((tag,"某层斜率落在自己的零内")); continue
    ratios[tag]=dict(level=a["lv_slope"]/b["lv_slope"], endpoint_rel=a["ep_rel"]/b["ep_rel"])
    print(f"  {tag:12s} 水平比 {ratios[tag]['level']:.3f} · 端点(相对基线)比 {ratios[tag]['endpoint_rel']:.3f}")
for tag,why in unread: print(f"  {tag:12s} **不可读** —— {why}")
print(f"  ⇒ 可读世代 {len(ratios)}/{len(COH)};`#778` 的合并比:水平 0.409–0.431 · 端点(相对基线)0.328–0.341")

G=Gate("#779 · 固定世代后那个比还在不在")
G.asserted("① 可读世代必须 ≥2(否则不许拿一个世代下结论)", bool(len(ratios)>=2),
           f"可读 {len(ratios)}/{len(COH)}:{list(ratios)}", kind="control")
if ratios:
    t0=list(ratios)[0]
    G.offset_control("② 该世代虔诚层的端点斜率须超它自己的零(比值的分子不能是噪声)",
                     effect=abs(res[t0]["虔诚"]["ep_slope"]), offset=res[t0]["虔诚"]["ep_null"],
                     spread=res[t0]["虔诚"]["ep_null"]*0.1,
                     null_kind="在该世代该层内打乱年份标签与端点占比的配对,保住每年的 n,只毁掉「哪一年配哪个占比」")
    G.offset_control("③ 同世代非虔诚层同样须超零(分母不能是噪声)",
                     effect=abs(res[t0]["非虔诚"]["ep_slope"]), offset=res[t0]["非虔诚"]["ep_null"],
                     spread=res[t0]["非虔诚"]["ep_null"]*0.1, null_kind="同上,换成非虔诚层")
print(); print(G)
print("\n"+"="*92)
if not all(r[2] for r in G.rows): v="**UNVERIFIED:闸没全过(可读世代不足或某层是噪声)**"
else:
    med=float(np.median([r["level"] for r in ratios.values()]))
    med_ep=float(np.median([r["endpoint_rel"] for r in ratios.values()]))
    if max(med,med_ep)<=0.7:
        v=(f"**A:固定出生世代后,比值仍是 水平 {med:.3f} · 端点(相对基线) {med_ep:.3f} —— "
           f"与 `#778` 的合并比(0.409–0.431 / 0.328–0.341)同量级 ⇒ 那个比不是年龄的影子**")
    else:
        v=(f"**B:固定世代后比值升到 水平 {med:.3f} · 端点 {med_ep:.3f} > 0.7 ⇒ 向 1 靠,"
           f"三条已发表的行要收窄成「老人变得少」**")
print(v)
print("⚠ 无论哪一支:**世代内 `age = year − cohort` 是恒等式 ⇒ 这不是「分离了年龄」,只是「固定世代后比值动不动」。**")
json.dump(dict(res=res,ratios=ratios,unreadable=unread,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"within_cohort.json","w"),ensure_ascii=False,indent=1)
