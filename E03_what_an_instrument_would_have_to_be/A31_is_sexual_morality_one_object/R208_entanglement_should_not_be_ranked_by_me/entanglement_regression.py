"""#769 · E03·A40·R208 —— 纠缠度不该由我来排

`#768` 报的是 **[5.8%, 12.4%]**,而那个区间的下界是**我事后挑的**:
我看完结果,认定 `homosex` 与「本人异性伴侣数」语义纠缠最小,于是取它作下界。
⚠ **那是一次拍板,不是一次测量** —— 与 `#764` 那条(阈值不该由我拍)同一族。

G1 估计量:**零纠缠处的解释份额** —— 把每题的份额对**它自己的纠缠度**回归,取截距。
纠缠度的可测定义:**该态度题与该行为变量的双变量 |ρ|**
(伴侣越多越可能有过婚前性行为 ⇒ `premarsx` 纠缠高;同性恋态度与本人异性伴侣数无关 ⇒ 纠缠低)。

⚠ **`xmarsex` 在本轮是可用的**:`#768` 把它排除是因为 `evstray`(出过轨)与它是同一件事的行为面/态度面;
   而本轮的行为变量是**一生伴侣数**,不是出轨 ⇒ **四题全用**,这一条写在跑之前。

⚠ **本轮结构上只此一具仪器,而这句话是量出来的**:估计量同时需要「被解释的那条线」(`obey` 一类)
   与本人性史。MFQ 无任何伴侣数变量(命中 0);NSFG 无 obey/authority 类(3,897 个变量名上命中 0);
   MSSCQ 是性自我概念,无权威类。**换不了仪器** —— 而能换的那一半(去年 vs 一生)`#768` 第三臂已在 NSFG 上换过。

⚠ **识别先于功效(G1)**:四个点、两个参数 ⇒ **只剩 2 个自由度**,截距的区间必然很宽。
   **先算出来再看要不要读它** —— 若截距的 95% 区间宽到跨越 `#768` 的整个 [5.8%, 12.4%],
   那本轮只是把「我拍的」换成「测不动的」,**必须如实说,不许当成收窄。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(208)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
ITEMS=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "numwomen":(0,988),"nummen":(0,988)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c:9s} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+ITEMS,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in ITEMS},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["life"]=M.numwomen+M.nummen
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
def cell(s, idx=None):
    sub=M[["obey","life",s]].dropna()
    if idx is not None: sub=sub.iloc[idx]
    raw=prho(sub.obey.to_numpy(),sub[s].to_numpy())
    par=prho(sub.obey.to_numpy(),sub[s].to_numpy(),sub[["life"]].to_numpy())
    ent=abs(prho(sub[s].to_numpy(),sub.life.to_numpy()))
    return dict(n=len(sub),share=1-par/raw,ent=ent)
print(f"\n=== 每题:纠缠度(该题与一生伴侣数的 |ρ|)与解释份额 ===")
print(f"  {'题':10s}{'n':>8s}{'纠缠度':>9s}{'份额':>9s}")
base={s:cell(s) for s in ITEMS}
for s in ITEMS: print(f"  {s:10s}{base[s]['n']:8d}{base[s]['ent']:9.4f}{base[s]['share']*100:8.1f}%")
x=np.array([base[s]["ent"] for s in ITEMS]); y=np.array([base[s]["share"] for s in ITEMS])
A=np.c_[np.ones(4),x]; b=np.linalg.lstsq(A,y,rcond=None)[0]
print(f"\n=== 份额 = a + b·纠缠度 ===")
print(f"  斜率 {b[1]:+.4f} · **截距(零纠缠处的份额) = {b[0]*100:.1f}%**")
print(f"  ⚠ 四个点、两个参数 ⇒ 只剩 2 个自由度")
B=800; ints=[]
for _ in range(B):
    idx=RNG.integers(0,10**9,1)  # 占位,下面按题各自 bootstrap 行
    yy=[]
    for s in ITEMS:
        n=base[s]["n"]; ii=RNG.integers(0,n,n); c=cell(s,ii); yy.append((c["ent"],c["share"]))
    xs=np.array([t[0] for t in yy]); ys=np.array([t[1] for t in yy])
    ints.append(np.linalg.lstsq(np.c_[np.ones(4),xs],ys,rcond=None)[0][0])
lo,hi=np.quantile(ints,[.025,.975])
print(f"  截距的 95% 自助区间 **[{lo*100:.1f}%, {hi*100:.1f}%]**(B={B})")
G=Gate("#769 · 零纠缠处的份额")
# ⚠⚠ 第一版这一条是**空洞的**:我比的是「残差和」与 0,而最小二乘的残差和**恒等于 0**,
#    所以它无论设计矩阵对不对都 PASS。`lib/gates.py` 的 `_degenerate` 当场把它判为
#    DEGENERATE(0 比 0)—— **「不会失败的检查」那一族,这次是库替我抓的,不是我。**
#    改成真检查:**已知斜率的合成数据上,回归必须把那个斜率取回来**;取不回来就是设计矩阵错了。
_xs=np.array([0.10,0.20,0.30,0.40]); _true_a,_true_b=0.05,0.30
_ys=_true_a+_true_b*_xs
_fit=np.linalg.lstsq(np.c_[np.ones(4),_xs],_ys,rcond=None)[0]
G.identity_control("① 合成数据上须取回已知斜率(仪器检查,能失败)",
                   observed=float(_fit[1]), expected=_true_b, tol=1e-9,
                   what="种一个斜率 0.30 的世界;取不回来就是设计矩阵错了")
G.identity_control("② 合成数据上须取回已知截距", observed=float(_fit[0]), expected=_true_a, tol=1e-9,
                   what="同上,截距 0.05")
print(); print(G)
print("\n"+"="*70)
w_old=0.124-0.058; w_new=hi-lo
if w_new>=w_old:
    v=(f"**测不动:截距区间宽 {w_new*100:.1f}pp,而 `#768` 的手排区间宽 {w_old*100:.1f}pp "
       f"⇒ 本轮把「我拍的」换成了「测不动的」,**不是收窄** —— `#768` 的区间维持**")
else:
    v=(f"**收窄:零纠缠处的份额 {b[0]*100:.1f}% [{lo*100:.1f}%, {hi*100:.1f}%],"
       f"区间从 {w_old*100:.1f}pp 收到 {w_new*100:.1f}pp,而且我的手不在排序里**")
print(v)
json.dump(dict(per_item=base,slope=float(b[1]),intercept=float(b[0]),ci=[float(lo),float(hi)],
               width_old=w_old,width_new=float(w_new),verdict=v),
          open(OUT/"entanglement.json","w"),ensure_ascii=False,indent=1)
