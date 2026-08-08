"""#774 · E03·A42·R213 —— 页面上每一处「解释份额」背后的那三条相关

⚠ **Production(如实标注)**:不产生新的关于人的信念,**把页面从错的读法改成对的读法**。

`#770` 证明:单控制量下 `share = 1 − (1 − r_XZ·r_YZ/r_XY) / √((1−r_XZ²)(1−r_YZ²))`,
是三条双变量相关的**确定函数**(解析 vs 数值差 1.9e−16),而只匹配那两条相关的高斯给出同一份额。
`#770`① 因此要求:**页面上每处份额都要写出它等价于哪两个相关。**

⚠ **而 `#760` 的产物里没有存政治那三条相关** —— 按 `#767`①(账本/页面上的每个数必须能在
`results/*.json` 里找到),**不许凭空写**。⇒ 本轮把三元组全部算出来并持久化,页面才有得引。

G1 估计量:对每个(控制构念 × 性态度题)格子,报 `(r_XY, r_XZ, r_YZ)` 与由它们**解析算出**的份额,
并与**数值管线**的份额对齐 —— 两者必须一致到 1e−9,否则页面引的数与脚本算的不是同一个东西。
⚠ 控制构念:政治(`polviews`)· 宗教(`attend`·`reliten`·`fund` 的 z 均值)。
⚠ **换不了仪器**:本轮的对象是**页面上已发表的那些数字自己**,不是任何外部数据 —— 只此一具。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "polviews":(1,7),"attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
def sp(a,b): return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)),method="spearman"))
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho_num(a,b,c):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    return float(np.corrcoef(resid(r(a),r(c).reshape(-1,1)),resid(r(b),r(c).reshape(-1,1)))[0,1])
def prho_ana(rxy,rxz,ryz): return (rxy-rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))

CTRL={"政治 polviews":["polviews"],"宗教 attend+reliten+fund":["attend","reliten","fund"]}
out={}; maxd=0.0
for cname,cols in CTRL.items():
    out[cname]={}
    print(f"\n=== 控制构念:{cname} ===")
    print(f"  {'题':10s}{'n':>7s}{'r_XY':>9s}{'r_XZ':>9s}{'r_YZ':>9s}{'份额(解析)':>12s}{'份额(数值)':>12s}")
    for s in SEX:
        sub=M[["obey"]+cols+[s]].dropna().copy()
        Zc=z(sub[cols]).mean(axis=1) if len(cols)>1 else sub[cols[0]]
        rxy,rxz,ryz=sp(sub.obey,sub[s]),sp(sub.obey,Zc),sp(sub[s],Zc)
        ana=prho_ana(rxy,rxz,ryz); num=prho_num(sub.obey.to_numpy(),sub[s].to_numpy(),Zc.to_numpy())
        sa,sn=1-ana/rxy,1-num/rxy; maxd=max(maxd,abs(sa-sn))
        out[cname][s]=dict(n=len(sub),rxy=rxy,rxz=rxz,ryz=ryz,share_analytic=sa,share_numeric=sn)
        print(f"  {s:10s}{len(sub):7d}{rxy:+9.4f}{rxz:+9.4f}{ryz:+9.4f}{sa*100:11.1f}%{sn*100:11.1f}%")
    med=float(np.median([out[cname][s]["share_analytic"] for s in SEX]))
    out[cname]["_median_share"]=med
    print(f"  中位份额 {med*100:.1f}%  ·  ⇒ 这一行页面上写的就是它,而它等价于上面那两列 r_XZ / r_YZ")

G=Gate("#774 · 页面份额的三元组")
# ⚠⚠ 第一版又把「差」与 0 比 —— **第四次**,而这条规则正是 `#770` 自己写下的:
#    「等式检查要比两个值,不要比它们的差与零」。**我写完那条规则,下一个脚本就破了它。**
#    ⇒ 逐格比两个份额本身(各约 0.1–0.4),非退化。
for _c in CTRL:
    for _s in SEX:
        G.identity_control(f"① 解析 vs 数值份额 · {_c[:2]}·{_s}",
                           observed=out[_c][_s]["share_numeric"], expected=out[_c][_s]["share_analytic"],
                           tol=1e-9, what="页面引的数与脚本算的必须是同一个量")
_bad=prho_ana(out["政治 polviews"]["premarsx"]["rxy"],out["政治 polviews"]["premarsx"]["rxz"],
              out["政治 polviews"]["premarsx"]["ryz"])
G.asserted("② 故意写错的式子(漏分母)必须给出不同的数 —— 证明 ① 能失败",
           bool(abs((out["政治 polviews"]["premarsx"]["rxy"]
                     -out["政治 polviews"]["premarsx"]["rxz"]*out["政治 polviews"]["premarsx"]["ryz"])-_bad)>1e-6),
           "漏掉分母后与正确式子的差必须 >1e−6", kind="control")
print(); print(G)
json.dump(out,open(OUT/"triples.json","w"),ensure_ascii=False,indent=1)
print(f"\n持久化 -> {OUT/'triples.json'}")
