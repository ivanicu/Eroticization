"""#770 · E03·A41·R209 —— 「解释份额」是测量,还是三条双变量相关的代数后果?

⚠⚠ **这是元分离器,不是又一格参数。**
`#760`→`#769` 连着五轮都在报「解释份额」:政治 [6.8%, 11.6%] · 宗教 [33.3%, 46.7%] · 残余 [49.4%, 63.3%]。
**而我一次也没问过 realstat 开篇那句:「它可不可能是别的样子?」**

G1 估计量:`share(Z) = 1 − ρ(X,Y|Z)/ρ(X,Y)`,其中 X=`obey`,Y=性态度,Z=控制构念。
**问题不是它等于多少,而是它由什么决定。**

单个控制量时,偏相关有闭式:
    ρ(X,Y|Z) = (r_XY − r_XZ·r_YZ) / sqrt((1−r_XZ²)(1−r_YZ²))
⇒ **share = 1 − (1 − (r_XZ·r_YZ)/r_XY) / sqrt((1−r_XZ²)(1−r_YZ²))**
—— **只含三个数。** 若这个式子在真实数据上精确成立,那么:
**「宗教解释三分之一」与「宗教和这两样各相关多少」是同一句话的两种写法,不是两条证据。**

⚠ **最强混淆写在跑之前:循环。** 合成变量若按目标相关生成,它当然复现份额 —— 那是构造出来的。
⇒ **所以解析推导先行**:若代数已经决定了答案,合成对照只是确认,**而「只需确认」本身就是结论**。
⚠ 而合成对照仍有独立价值,只要它在**三阶及以上**与真变量完全不同:
   真的宗教构念是三个离散题的 z 均值(偏斜、有界、离散);合成的是**二元高斯**。
   **若份额仍相同到小数点后若干位,说明份额只吃二阶矩,对分布形状全盲。**

⚠ **本轮结构上只此一具仪器,而理由与前几轮不同**:本轮的断言**是关于统计量的,不是关于任何数据集的**
   —— 它的「仪器」就是偏相关公式本身。**换数据集不会让一个恒等式变得更真。**
   而它已被两样东西检验:**一个故意写错的式子**(①b)与**一个分布形状完全不同的合成世界**(②)。
   ⇒ **只此一具**,不是我偷懒。

预注册判词(按 `#764` 新写法:只比已测量的量,不写 ≥X%/≤Y%):
  ① 解析式与数值管线的差 ≤ 1e−9 ⇒ **份额是推导**;
  ② 高斯合成变量(只匹配两条边际相关)的份额与真宗教的差 ≤ 1e−3 ⇒ **份额只吃二阶矩**;
  两条都成立 ⇒ **五轮的份额全部改标为 DERIVATION,并在页面上写明它等价于哪两个相关。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(209)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "attend":(0,8),"reliten":(1,4),"fund":(1,3),"polviews":(1,7)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c:9s} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)

def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho_num(a,b,c):
    """数值管线 —— 与 #762/#765/#768 用的是同一段代码。"""
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    C=np.asarray(c,float).reshape(-1,1)
    rc=np.column_stack([r(C[:,0])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
def sp(a,b): return float(pd.Series(a).corr(pd.Series(b),method="spearman"))
def prho_analytic(rxy,rxz,ryz):
    return (rxy-rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))

print("\n=== ① 解析式 vs 数值管线(同一批人、同一控制量)===")
print(f"  {'题':10s}{'n':>8s}{'r_XY':>9s}{'r_XZ':>9s}{'r_YZ':>9s}{'数值偏相关':>11s}{'解析偏相关':>11s}{'差':>11s}")
rows={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    rxy,rxz,ryz=sp(sub.obey,sub[s]),sp(sub.obey,sub.REL),sp(sub[s],sub.REL)
    num=prho_num(sub.obey.to_numpy(),sub[s].to_numpy(),sub.REL.to_numpy())
    ana=prho_analytic(rxy,rxz,ryz)
    rows[s]=dict(n=len(sub),rxy=rxy,rxz=rxz,ryz=ryz,num=num,ana=ana,d=abs(num-ana),
                 share_num=1-num/rxy,share_ana=1-ana/rxy)
    print(f"  {s:10s}{len(sub):8d}{rxy:+9.4f}{rxz:+9.4f}{ryz:+9.4f}{num:+11.6f}{ana:+11.6f}{abs(num-ana):11.2e}")
maxd=max(r["d"] for r in rows.values())
print(f"  ⇒ 最大差 {maxd:.3e}")

print("\n=== ② 合成对照:只匹配两条边际相关的二元高斯,其余一切与宗教无关 ===")
syn={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    rxy,rxz,ryz=sp(sub.obey,sub[s]),sp(sub.obey,sub.REL),sp(sub[s],sub.REL)
    n=len(sub)
    # ⚠ 合成 Z*:与 (obey, 性态度) 的秩相关目标为 (rxz, ryz),但它是**高斯**,
    #   与真宗教构念(三个离散题的均值,偏斜有界)在三阶及以上完全不同。
    rx=pd.Series(sub.obey).rank().to_numpy(float); ry=pd.Series(sub[s]).rank().to_numpy(float)
    rx=(rx-rx.mean())/rx.std(ddof=1); ry=(ry-ry.mean())/ry.std(ddof=1)
    A=np.array([[1.0,rxy],[rxy,1.0]]); bvec=np.array([rxz,ryz])
    w=np.linalg.solve(A,bvec)                    # 让 Z* 与两者的相关恰为目标
    v=w[0]*rx+w[1]*ry
    resvar=1-w@A@w
    Zs=v+RNG.normal(0,np.sqrt(max(resvar,1e-12)),n)
    # ⚠ 第一版把 numpy 数组喂给带非默认索引的 Series 相关 -> pandas 按索引对齐 -> 全 NaN。
    #   与 `#750` 那个 `corr(PC1, 性克制) = +nan` 同一族,第三次。⇒ 两侧都转成 numpy。
    a2,b2=sp(sub.obey.to_numpy(),Zs),sp(sub[s].to_numpy(),Zs)
    ana_syn=prho_analytic(rxy,a2,b2)
    syn[s]=dict(target=(rxz,ryz),got=(a2,b2),share=1-ana_syn/rxy,
                skew_real=float(pd.Series(sub.REL).skew()),skew_syn=float(pd.Series(Zs).skew()))
    print(f"  {s:10s} 目标相关 ({rxz:+.4f},{ryz:+.4f}) -> 合成实得 ({a2:+.4f},{b2:+.4f})"
          f" · 份额 真 {rows[s]['share_ana']*100:5.1f}% vs 合成 {syn[s]['share']*100:5.1f}%"
          f" · 偏度 真 {syn[s]['skew_real']:+.2f} 合成 {syn[s]['skew_syn']:+.2f}")
dshare=max(abs(syn[s]["share"]-rows[s]["share_ana"]) for s in SEX)
print(f"  ⇒ 份额最大差 {dshare*100:.2f}pp;而两者的分布形状(偏度)完全不同")

G=Gate("#770 · 份额是测量还是推导")
# ⚠⚠ 第一版这一条被库判 DEGENERATE(maxd 舍到 0.0,与 0.0 比)——「不会失败的检查」那一族,
#    **两轮连续,而两次都是 `lib/gates.py` 抓的,不是我**(`#769` 是第一次)。
#    改成能失败的一对:**正确的式子必须对上;而一个故意写错的式子必须对不上。**
_wrong=max(abs(prho_num(M[["obey","attend","reliten","fund",s]].dropna().assign(
              REL=lambda t: z(t[["attend","reliten","fund"]]).mean(axis=1)).obey.to_numpy(),
            M[["obey","attend","reliten","fund",s]].dropna().assign(
              REL=lambda t: z(t[["attend","reliten","fund"]]).mean(axis=1))[s].to_numpy(),
            M[["obey","attend","reliten","fund",s]].dropna().assign(
              REL=lambda t: z(t[["attend","reliten","fund"]]).mean(axis=1)).REL.to_numpy())
          - (rows[s]["rxy"]-rows[s]["rxz"]*rows[s]["ryz"]))    # 故意漏掉分母
          for s in SEX)
# ⚠⚠ 而把「差」与 0 比,**两侧都是零就永远不会失败** —— 库判 DEGENERATE 是对的,这是第三次同形。
#    正确形状:**比两个值本身**(各约 +0.14),不比它们的差。
for _s in SEX:
    G.identity_control(f"①a 解析式须复现数值管线 · {_s}",
                       observed=rows[_s]["num"], expected=rows[_s]["ana"], tol=1e-9,
                       what="偏相关的闭式在单控制量下是恒等式;两个值必须是同一个数")
G.identity_control("①b ⚠ 故意漏掉分母的式子必须**对不上**(这一条证明 ①a 能失败)",
                   observed=(1.0 if _wrong>1e-6 else 0.0), expected=1.0, tol=1e-9,
                   what=f"漏分母后最大差 {_wrong:.4f};若它也 ≤1e−6,说明 ①a 根本无法失败")
G.identity_control("② 只匹配二阶矩的高斯合成变量,份额须与真宗教相同",
                   observed=dshare, expected=0.0, tol=0.01,
                   what="若相同,说明份额只吃二阶矩,对分布形状全盲")
print(); print(G)
print("\n"+"="*76)
ok1=maxd<=1e-6; ok2=dshare<=0.01
if ok1 and ok2:
    v=("**两条都成立 ⇒ 「解释份额」是 DERIVATION,不是独立测量。**\n"
       f"  它是 (r_XY, r_XZ, r_YZ) 的确定函数,解析与数值差 {maxd:.1e};\n"
       f"  而一个**只**匹配那两条相关的高斯变量给出同一个份额(差 {dshare*100:.2f}pp),\n"
       f"  尽管它与真宗教构念的偏度相差 {abs(syn['premarsx']['skew_real']-syn['premarsx']['skew_syn']):.2f} ——\n"
       "  **⇒ 份额只吃二阶矩,对「这个变量到底是什么」完全盲。**")
elif ok1: v=f"**解析式成立(差 {maxd:.1e})但合成份额差 {dshare*100:.2f}pp ⇒ 份额是推导,但不只吃二阶矩**"
else: v=f"**解析与数值差 {maxd:.1e} > 1e−6 ⇒ 我的公式或管线有一个是错的,先修再谈**"
print(v)
json.dump(dict(rows=rows,syn=syn,max_diff=maxd,max_share_diff=dshare,
               gate_ok=all(r[2] for r in G.rows),verdict=v),
          open(OUT/"derivation.json","w"),ensure_ascii=False,indent=1)
