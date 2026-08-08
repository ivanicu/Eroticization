"""#777 · E03·A42·R216 —— 婚外性纹丝不动,是真稳定,还是早就到顶所以量不动?

`#776` 量到:扣掉有界量表的算术后,最虔诚层 1988–2024 的 `var_obs/var_max` 在
`homosex` +2.20×、`premarsx` +1.85× 自己的零上**上升**,而 **`xmarsex` 0.13×,纹丝不动**。
`#776`① 预注册:**那两件事可能是同一件事,也可能不是** ——
  **A 真稳定**:婚外性上的共识几十年不动,不是量表的事;
  **B 早熟的天花板**:它只是**更早**到顶,到顶之前 ratio 也在升(与其余三题同族)。

⚠⚠ **我不欢迎的是 B** —— 它会把我上一轮刚写上页面的「唯一例外」收窄成「更早到顶」。

G1 估计量:在 **`xmarsex` 端点占比 < C 的那些年**里,虔诚层 ratio 对年份的斜率;
**C 扫成规格曲线**(窗口是我选的,不能只报一格)。

⚠ **三条混淆写在跑之前:**
① **年份少一半 ⇒ 零变宽。先算 MDE** —— 每个窗口报它自己的零 95% 分位与年数,撑不起就判不了。
② **早期年份 n 更小 ⇒ ratio 本身噪声更大 ⇒ 零更宽**。**这个方向对我不利**(更难判 B),如实记。
③ **最锋利的控制:同一窗口内比其余三题。** 若 `xmarsex` 平,而 `homosex` 在**同样的年数**下
   显著为正,**那「平」就不是功效问题** —— 这一条把窗口长度这个混淆整个消掉。
⚠ **正控**:往该窗口的 ratio 序列里种一个已知斜率,必须取回;**且在 g=0 时必须不过。**

⚠ **换不了仪器,与 `#776` 同一条且同样量过**:本轮需要**同一道题在 ≥20 年上的重复测量**。
   MFQ 单次采集无年代;NSFG 的 `SXOK18` 在 2017–19 卷字典里不存在;SCCS 的单位是社会。
   ⇒ **GSS 只此一具。**

预注册判词(按 `#764` 新写法:只比已测量的量,各带自己的零):
  在**至少一个**端点占比 <C 的窗口里,`xmarsex` 斜率为正且超它自己的零 ⇒ **B**;
  在**所有**这样的窗口里都不超零,**而同窗口的 `homosex` 超零** ⇒ **A**(平不是功效问题);
  两者都不超零 ⇒ **判不了**(窗口撑不起),不许解释。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(216)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["xmarsex","homosex","premarsx","teensex"]
VALID={"xmarsex":(1,4),"homosex":(1,4),"premarsx":(1,4),"teensex":(1,4),
       "attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict"): M[c]=-M[c]+5   # 高=严,留在 1..4
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
vmax=lambda mu:(9*((mu-1)/3.0)*(1-(mu-1)/3.0))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))

# 逐题 · 逐年:虔诚层的 ratio 与端点占比
SER={}
for s in SEX:
    sub=M[[s,"attend","reliten","fund","year"]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    sub["k"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
    g=sub[sub.k==2]; rows=[]
    for y,gy in g.groupby("year"):
        if len(gy)<120: continue
        mu=float(gy[s].mean()); vm=vmax(mu)
        if vm<=1e-6: continue
        rows.append((int(y),float(gy[s].var(ddof=1))/vm,float((gy[s]==4).mean()),len(gy)))
    SER[s]=rows
print(f"\n=== 各题的年数(虔诚层,每年 n≥120)===")
for s in SEX: print(f"  {s:10s} {len(SER[s])} 年 · 端点占比 {SER[s][0][2]*100:.1f}% → {SER[s][-1][2]*100:.1f}%")

def fit(rows):
    yr=np.array([r[0] for r in rows],float); rt=np.array([r[1] for r in rows])
    b=slope(yr,rt); nul=[slope(RNG.permutation(yr),rt) for _ in range(4000)]
    return b, float(np.quantile(np.abs(nul),.95))

CUTS=[0.45,0.50,0.55,0.60,0.65,0.70]
print(f"\n=== G4 窗口扫描:只留 `xmarsex` 端点占比 < C 的那些年,四题在**同一批年份**上各自估斜率 ===")
print(f"  {'C':>5s}{'年数':>5s}{'年份范围':>14s}"+"".join(f"{s[:8]:>20s}" for s in SEX))
grid={}
xm_years={r[0]:r[2] for r in SER["xmarsex"]}
for C in CUTS:
    keep=sorted([y for y,t in xm_years.items() if t<C])
    if len(keep)<6: print(f"  {C:5.2f}{len(keep):5d}   年数不足,跳过"); continue
    cell={}
    for s in SEX:
        rows=[r for r in SER[s] if r[0] in keep]
        if len(rows)<6: continue
        b,q=fit(rows); cell[s]=dict(n_years=len(rows),slope=b,null95=q,ratio=abs(b)/q,sign=np.sign(b))
    grid[C]=dict(years=keep,cell=cell)
    print(f"  {C:5.2f}{len(keep):5d}{f'{keep[0]}–{keep[-1]}':>14s}"
          +"".join(f"{cell[s]['slope']:+11.6f}{cell[s]['ratio']:8.2f}×" if s in cell else f"{'-':>20s}" for s in SEX))

# ---- 闸 ----
C0=0.60; base=[r for r in SER["xmarsex"] if r[0] in grid[C0]["years"]] if C0 in grid else SER["xmarsex"]
yr=np.array([r[0] for r in base],float); rt=np.array([r[1] for r in base])
G=Gate("#777 · 婚外性:真稳定还是早熟的天花板")
_,q0=fit(base)
plant=q0*3.0
b_plant=slope(yr, rt+plant*(yr-yr.mean()))
G.identity_control("① 正控:种一个 3× 零的斜率必须原样取回",
                   observed=b_plant-slope(yr,rt), expected=plant, tol=abs(plant)*0.02,
                   what=f"在 C={C0} 的窗口({len(base)} 年)里种 {plant:+.6f}/年")
G.asserted("② 正控必须在 g=0 时不过(能失败)",
           bool(abs(slope(yr,rt+0.0*(yr-yr.mean()))-slope(yr,rt))<abs(plant)*0.02 and plant>q0),
           f"g=0 时取回 0(≤{abs(plant)*0.02:.2e});而 g=3×零 时取回 {plant:+.6f} > 零 {q0:.6f}", kind="control")
# ⚠⚠ 第一版这一条在窗口为空时用 `effect=0.0, offset=1.0` 造了一个**假的数值比较**,
#    于是 FAIL 读起来像「homosex 没功效」,而真相是**根本没有窗口**。
#    正是 realstat「控制因自己的理由而失败」那一行 —— 改成如实报「窗口是空的」。
ho=grid.get(C0,{}).get("cell",{}).get("homosex")
G.asserted("③ 前提检查:`xmarsex` 端点占比 <C 的窗口必须非空(否则整个设计问不了这个问题)",
           bool(ho is not None),
           f"C={C0} 时可用年数 {len(grid.get(C0,{}).get('years',[]))};"
           f"而 `xmarsex` 的端点占比全期最低 {min(r[2] for r in SER['xmarsex'])*100:.1f}% —— "
           f"**1988 年就已经到顶,「到顶之前」在 GSS 里不存在**", kind="control")
print(); print(G)

# ============================ 副产品:端点占比自己的走势 ============================
# ⚠ 预注册的问题问不了,而**数据在回答一个我没问的问题**:虔诚层的端点占比各自怎么走。
#   它必须自带零才能报 —— 否则就是「看图说话」。
print("\n=== 副产品 · 虔诚层端点占比的年趋势(各带自己的零)===")
print(f"  {'题':10s}{'年数':>5s}{'首年':>9s}{'末年':>9s}{'斜率/年':>11s}{'零95%':>10s}{'倍数':>8s}")
END={}
for s in SEX:
    rows=SER[s]; yr=np.array([r[0] for r in rows],float); tp=np.array([r[2] for r in rows])
    b=slope(yr,tp); nul=[slope(RNG.permutation(yr),tp) for _ in range(4000)]
    q=float(np.quantile(np.abs(nul),.95))
    END[s]=dict(n_years=len(rows),first=tp[0],last=tp[-1],slope=b,null95=q,ratio=abs(b)/q)
    print(f"  {s:10s}{len(rows):5d}{tp[0]*100:8.1f}%{tp[-1]*100:8.1f}%{b:+11.6f}{q:10.6f}{abs(b)/q:7.2f}×")
G2=Gate("#777 副产品 · 端点占比的年趋势")
G2.offset_control("④ `homosex` 端点占比的斜率须超出它自己的零",
                  effect=abs(END["homosex"]["slope"]), offset=END["homosex"]["null95"],
                  spread=END["homosex"]["null95"]*0.1,
                  null_kind="打乱年份标签与端点占比的配对,保住每年的 n,只毁掉「哪一年配哪个占比」")
G2.asserted("⑤ 而 `xmarsex` 的同一统计量必须**不**超零,否则「只有它没动」这句话不成立",
            bool(END["xmarsex"]["ratio"]<=1.0),
            f"`xmarsex` {END['xmarsex']['ratio']:.2f}× vs `homosex` {END['homosex']['ratio']:.2f}×", kind="control")
print(); print(G2)

pos=[(C,g["cell"]["xmarsex"]) for C,g in grid.items()
     if "xmarsex" in g["cell"] and g["cell"]["xmarsex"]["slope"]>0 and g["cell"]["xmarsex"]["ratio"]>1.0]
ho_ok=[(C,g["cell"]["homosex"]) for C,g in grid.items()
       if "homosex" in g["cell"] and g["cell"]["homosex"]["ratio"]>1.0]
print("\n"+"="*80)
if not all(r[2] for r in G.rows):
    v=("**UNVERIFIED —— 而原因不是仪器坏了,是设计的前提被数据推翻**:\n"
       f"  `xmarsex` 的端点占比全期最低 {min(r[2] for r in SER['xmarsex'])*100:.1f}%(1988 年),最高 "
       f"{max(r[2] for r in SER['xmarsex'])*100:.1f}% ⇒ **「到顶之前」那个窗口在 GSS 里不存在。**\n"
       "  ⇒ **A(真稳定)与 B(早熟的天花板)在这套数据上不可分** —— 记为能力边界,不是判词。")
elif pos:
    v=(f"**B 早熟的天花板:`xmarsex` 在 {len(pos)}/{len(grid)} 个窗口里斜率为正且超自己的零 "
       f"(C={[c for c,_ in pos]})⇒ `#776` 的「唯一例外」要收窄为「更早到顶」**")
elif ho_ok:
    v=(f"**A 真稳定:`xmarsex` 在 {len(grid)} 个窗口里全部不超零,而同窗口的 `homosex` 在 "
       f"{len(ho_ok)}/{len(grid)} 个窗口里超零 ⇒ 那个「平」不是功效问题**")
else: v="**判不了:窗口里连 `homosex` 都不超零 ⇒ 这些窗口撑不起这个问题**"
print(v)
json.dump(dict(series={s:SER[s] for s in SEX},grid={str(k):v for k,v in grid.items()},
               endpoint_trend=END,verdict=v,gate_ok=all(r[2] for r in G.rows),
               gate2_ok=all(r[2] for r in G2.rows)),
          open(OUT/"early_window.json","w"),ensure_ascii=False,indent=1)
