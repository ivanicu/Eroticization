"""#772 · E03·A41·R211 —— 换一个不吃方差的统计量,方向落在哪一侧?

`#771` 两臂方向相反:GSS 说越虔诚 ρ(obey, 性态度) 越低,MFQ 说越虔诚 ρ(chastity, 权威) 越高。
**而分歧的假设就在每层 sd:GSS 的下降伴着 37% 的方差塌缩,MFQ 的上升伴着平坦的方差。**
⇒ 本轮换**不吃 sd(Y)** 的统计量,看方向落在哪一侧。

⚠⚠ **而我第一反应是「分层内秩回归斜率」—— 那是错的,写在这里免得下次再犯:**
层内秩化会把 `sd(rank Y)` 固定成 n 的函数,于是**斜率退化回相关本身**,等于什么都没修。
**不吃 sd(Y) 的是原始值上的 OLS 斜率** `cov(X,Y)/var(X)` —— 它是「每多一分服从,态度严多少分」。

G1 估计量,三个,同轮报(G4 规格曲线):
  ① **原始值 OLS 斜率**(层内),单位 = 态度分/服从分 —— 不除以 sd(Y)
  ② **边际匹配**:把各层 Y 的类别分布重抽成与**合并分布**相同,再算层内 ρ —— 直接消掉边际差异
  ③ **每层的量表端点占比** —— 世界 C(删失)的直接证据

⚠ **最强混淆写在跑之前:斜率修方差,修不了删失。**
若最虔诚层挤在「总是错」这个端点,**斜率同样被压平**,而那不是关系变弱,是量表到顶了。
⇒ **③ 必须与 ①② 同轮报**;若端点占比在高虔诚层 **>60%**,**判世界 C:①② 都不可读**,记为能力边界。

预注册判词(按 `#764` 新写法:只比已测量的量,各带自己的零):
  先看 ③。端点占比 >60% ⇒ **C:能力边界,不许解释方向**;
  否则比 ① 与 ② 的方向:两者都随虔诚下降 ⇒ A(GSS 方向真);两者都平坦 ⇒ B(方差假象);
  不一致 ⇒ 报出来,不选边。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(211)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
flip=aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important")
for c in flip: M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
def sp(a,b): return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)),method="spearman"))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))

print("\n=== ③ 先看删失:每层 Y 在量表端点的占比(世界 C 的直接证据)===")
print(f"  {'题':10s}{'层':>4s}{'n':>7s}{'最严端点%':>10s}{'最松端点%':>10s}{'sd(Y)':>8s}{'sd(X)':>8s}")
CEIL={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    sub["k"]=pd.qcut(sub.REL,3,labels=False,duplicates="drop")
    hi_end=sub[s].max(); lo_end=sub[s].min()
    CEIL[s]=[]
    for k in sorted(sub.k.unique()):
        g=sub[sub.k==k]
        CEIL[s].append(dict(n=len(g),top=float((g[s]==hi_end).mean()),bot=float((g[s]==lo_end).mean()),
                            sdy=float(g[s].std(ddof=1)),sdx=float(g.obey.std(ddof=1))))
        print(f"  {s if k==0 else '':10s}{int(k):4d}{len(g):7d}{CEIL[s][-1]['top']*100:9.1f}%{CEIL[s][-1]['bot']*100:9.1f}%"
              f"{CEIL[s][-1]['sdy']:8.3f}{CEIL[s][-1]['sdx']:8.3f}")
maxtop=max(c["top"] for s in SEX for c in CEIL[s])
print(f"  ⇒ 端点占比最高的一格 {maxtop*100:.1f}%(预注册阈值 60%)")

print("\n=== ① 原始值 OLS 斜率(层内;单位 = 态度分/服从分,不除以 sd(Y))===")
print(f"  {'题':10s}{'低':>10s}{'中':>10s}{'高':>10s}{'极差':>9s}")
SL={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    sub["k"]=pd.qcut(sub.REL,3,labels=False,duplicates="drop")
    b=[slope(sub[sub.k==k].obey,sub[sub.k==k][s]) for k in sorted(sub.k.unique())]
    SL[s]=b; print(f"  {s:10s}"+"".join(f"{v:+10.4f}" for v in b)+f"{max(b)-min(b):9.4f}")

print("\n=== ② 边际匹配:把各层 Y 的类别分布重抽成与合并分布相同,再算层内 ρ ===")
print(f"  {'题':10s}{'低':>10s}{'中':>10s}{'高':>10s}{'极差':>9s}   (原始 ρ 极差)")
MM={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    sub["k"]=pd.qcut(sub.REL,3,labels=False,duplicates="drop")
    pooled=sub[s].value_counts(normalize=True)
    raw=[sp(sub[sub.k==k].obey,sub[sub.k==k][s]) for k in sorted(sub.k.unique())]
    reps=[]
    for _ in range(200):
        rs=[]
        for k in sorted(sub.k.unique()):
            g=sub[sub.k==k]; take=[]
            for cval,share in pooled.items():
                pool=g[g[s]==cval]
                need=int(round(share*len(g)))
                if len(pool)==0 or need==0: continue
                take.append(pool.iloc[RNG.integers(0,len(pool),need)])
            gg=pd.concat(take)
            rs.append(sp(gg.obey,gg[s]))
        reps.append(rs)
    mm=np.median(np.array(reps),axis=0)
    MM[s]=dict(matched=mm.tolist(),raw=raw)
    print(f"  {s:10s}"+"".join(f"{v:+10.4f}" for v in mm)+f"{mm.max()-mm.min():9.4f}   ({max(raw)-min(raw):.4f})")

# ⚠⚠ 第一版这两条都是**把常数与自己比**(1.0 vs 1.0)——空洞检查,而 `lib/gates.py` 的
#    `_degenerate` **没有抓到**,因为它只认「全零」。⇒ 库有一个缺口,已在 `#773` 补。
#    改成用**实测量**:匹配前后的层间边际最大差。
G=Gate("#772 · 不吃方差的统计量")
sub=M[["obey","attend","reliten","fund","xmarsex"]].dropna().copy()
sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1); sub["k"]=pd.qcut(sub.REL,3,labels=False,duplicates="drop")
pooled=sub["xmarsex"].value_counts(normalize=True).sort_index()
def marg_gap(frame):
    return max(abs(frame[frame.k==k]["xmarsex"].value_counts(normalize=True).reindex(pooled.index).fillna(0)-pooled).max()
               for k in sorted(frame.k.unique()))
before=marg_gap(sub)
take=[]
for k in sorted(sub.k.unique()):
    g=sub[sub.k==k]
    for cval,share in pooled.items():
        pool=g[g["xmarsex"]==cval]; need=int(round(share*len(g)))
        if len(pool) and need: take.append(pool.iloc[RNG.integers(0,len(pool),need)])
after=marg_gap(pd.concat(take))
G.identity_control("① 匹配**后**层间边际须一致(操作成功了吗)", observed=after, expected=0.0, tol=0.01,
                   what=f"匹配后最大边际差 {after:.4f}")
# ⚠⚠ 第二版把这一条写成「必须 FAIL 的 identity_control」—— 于是汇总又变成 UNVERIFIED,
#    而 `#767` 正是记下这一条的那一轮:**故意反向的检查不能放进同一个 Gate**。我又犯了一次。
#    ⇒ 改成 PASS 形状的断言:**「匹配前的边际差必须大于容差」本身就是一个可失败的命题。**
G.asserted("② 匹配**前**层间边际必须大于容差(证明 ① 不是空操作)",
           bool(before>0.01),
           f"匹配前最大边际差 {before:.4f} > 0.01;若不大于,①(匹配后≈0)就无从失败", kind="control")
print(); print(G)
print("\n"+"="*72)
if maxtop>0.60:
    v=(f"**世界 C:端点占比最高 {maxtop*100:.1f}% > 60% ⇒ 量表到顶,①② 都不可读 —— 记为能力边界**")
else:
    dn_sl=sum(1 for s in SEX if SL[s][0]>SL[s][-1]); dn_mm=sum(1 for s in SEX if MM[s]["matched"][0]>MM[s]["matched"][-1])
    if dn_sl>=3 and dn_mm>=3: v=f"**A:斜率 {dn_sl}/4、边际匹配后 {dn_mm}/4 仍随虔诚下降 ⇒ GSS 的方向不是方差假象**"
    elif dn_sl<=1 and dn_mm<=1: v=f"**B:斜率 {dn_sl}/4、匹配后 {dn_mm}/4 下降 ⇒ GSS 那个方向是方差假象,撤销**"
    else: v=f"**不一致:斜率 {dn_sl}/4 下降、匹配后 {dn_mm}/4 下降 ⇒ 报出来,不选边**"
print(v)
json.dump(dict(ceiling=CEIL,slopes=SL,matched=MM,maxtop=maxtop,verdict=v,
               gate_ok=all(r[2] for r in G.rows)),open(OUT/"variance_robust.json","w"),ensure_ascii=False,indent=1)
