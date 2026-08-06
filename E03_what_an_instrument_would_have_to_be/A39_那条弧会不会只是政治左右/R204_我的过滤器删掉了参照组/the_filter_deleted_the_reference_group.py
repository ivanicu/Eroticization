"""#765 · E03·A39·R204 —— 我的过滤器删掉了参照组,而残余正是被它抬高的那个数

⚠⚠ **本轮起于一个已发表结果里的真缺陷,不是一个新问题。**
`#762`/`#763` 用 `.where(x>0)` 做笼统缺失过滤。而 GSS 的 `.dta` **本来就把 DK/NA 设成缺失了**
(`polviews` 只有 1–7,没有 8/9)—— **是我又「清洗」了一遍,清掉的是真数据**:
  `attend` 码 **0 = "never",n=14,883 —— 最大的一档**,被整档删除;
  `partyid` 码 0 = "strong democrat",n=12,327(本轮才引入,尚未进过任何已发表结果)。

**跑之前写下的方向**:删掉从不去教堂的人 = **对正在被偏掉的那个变量做全距截断**
⇒ 一个被截断的控制量能解释的更少 ⇒ **宗教的份额被低估、残余被高估。**
而残余(62.4%)正是 `#764` 刚提为主角的那个数 —— **所以这个缺陷的方向指向我的头号声明。**

**暴露面已量,并且锁死在一处**(不是猜的):全仓 13 个脚本用 `x>0`,
其中同时碰 `attend`/`partyid` 的**只有 `#762`**;`R112`/`R113` 用的是 `pd.cut(attend,[-1,1,3,5,8])`,
**下界 −1,码 0 被包含,干净**;`#763` 的 MFQ 臂用 `religatt_num`(有效值含 0,n=4,142)且**直接 dropna,干净**。

G1 估计量:与 `#762` **完全相同** —— `ρ(obey, 性条目 | 控制量)` 的保留率网格。
**唯一的差别是缺失码的处理。** ⇒ 这是一次 A/B:**同一段代码,两种过滤。**

⚠ 判词**按 `#764` 新写法**:只比已测量的量,各带自己的零;**不写 ≥X%/≤Y%。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned
from lib.gates import Gate
RNG=np.random.default_rng(204)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]; REL=["attend","reliten","fund"]

# ⚠ 逐变量的有效码,从值标签读出来的,不是笼统 >0(硬规则①)
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"teensex":(1,4),"homosex":(1,4),
       "polviews":(1,7),"reliten":(1,4),"fund":(1,3),
       "attend":(0,8),          # ⚠ 0 = never,有效
       "partyid":(0,6)}         # ⚠ 0 = strong democrat 有效;7 = other party **不在连续统上**,剔除
cols=list(VALID)
d=pd.read_stata(gp,columns=["year"]+cols,convert_categoricals=False)
def clean(c, blanket):
    v=pd.to_numeric(d[c],errors="coerce")
    if blanket: return v.where(v>0).where(lambda x: x<=VALID[c][1])
    lo,hi=VALID[c]; return v.where((v>=lo)&(v<=hi))
for tag,blanket in (("旧(笼统 x>0)",True),("新(逐变量有效码)",False)):
    M=pd.DataFrame({c:clean(c,blanket) for c in cols}); M["year"]=d.year
    n_attend=int(M.attend.notna().sum()); n0=int((M.attend==0).sum())
    print(f"{tag:18s} attend n={n_attend:6d} 其中「从不」={n0:6d}")

cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
flip=aligned({c:cats[c] for c in SEX},"strict") | aligned({c:cats[c] for c in ["obey"]},"important")
print(f"\n方向由 aligned() 定 -> 要翻 {sorted(flip)}")

def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])

SPECS={"政治 polviews":["polviews"],"出席 attend":["attend"],"认同强度 reliten":["reliten"],
       "原教旨 fund":["fund"],"三个宗教量一起":REL,"宗教 + 政治全放":REL+["polviews"]}
ALL={}
for tag,blanket in (("旧(笼统 x>0)",True),("新(逐变量有效码)",False)):
    M=pd.DataFrame({c:clean(c,blanket) for c in cols}); M["year"]=d.year
    for c in flip: M[c]=-M[c]
    keeps={}; ns={}
    for nm,cc in SPECS.items():
        ks=[]
        for s in SEX:
            sub=M[["obey"]+REL+["polviews",s]].dropna()
            raw=prho(sub.obey.to_numpy(),sub[s].to_numpy())
            v=prho(sub.obey.to_numpy(),sub[s].to_numpy(),sub[cc].to_numpy())
            floor=3*1.65/np.sqrt(len(sub))
            if abs(raw)>=floor: ks.append(v/raw)
            ns[s]=len(sub)
        keeps[nm]=float(np.median(ks)) if ks else None
    ALL[tag]=dict(keeps=keeps,n=ns)
    print(f"\n=== {tag} · n={ns} ===")
    for nm in SPECS: print(f"  {nm:20s} 保留 {keeps[nm]*100:6.1f}%  解释 {(1-keeps[nm])*100:5.1f}%")

old,new=ALL["旧(笼统 x>0)"]["keeps"],ALL["新(逐变量有效码)"]["keeps"]
print("\n=== A/B:同一段代码,两种过滤 ===")
print(f"  {'':20s}{'旧 解释':>9s}{'新 解释':>9s}{'移动':>9s}")
for nm in SPECS:
    print(f"  {nm:20s}{(1-old[nm])*100:8.1f}%{(1-new[nm])*100:8.1f}%{((old[nm]-new[nm]))*100:+8.1f}pp")
res_old,res_new=old["宗教 + 政治全放"],new["宗教 + 政治全放"]
rel_new,pol_new=1-new["三个宗教量一起"],1-new["政治 polviews"]

# ---- 闸 ----
M=pd.DataFrame({c:clean(c,False) for c in cols}); M["year"]=d.year
for c in flip: M[c]=-M[c]
sub=M[["obey"]+REL+["polviews","premarsx"]].dropna()
raw=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy())
pc=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),RNG.normal(0,1e-9,(len(sub),1)))
nul=[prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[REL].to_numpy()[RNG.permutation(len(sub))]) for _ in range(300)]
G=Gate("#765 · 修好过滤器之后")
G.identity_control("① 常数控制须回到偏前", observed=pc, expected=raw, tol=0.005, what="仪器活着吗")
G.identity_control("② 打乱三个宗教量须回到偏前", observed=float(np.median(nul)), expected=raw, tol=0.005,
                   what="打乱只毁掉「谁的宗教配谁的态度」")
G.offset_control("③ 偏掉真宗教后须显著低于偏前",
                 effect=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[REL].to_numpy()),
                 offset=raw, spread=float(np.std(nul)),
                 null_kind="同一批人、同一对题,唯一差别是控制量是真的三个宗教量还是打乱后的三个宗教量")
print(); print(G)

print("\n"+"="*72)
if not all(r[2] for r in G.rows): print("**UNVERIFIED:闸没全过**")
else:
    # ⚠ #764 新写法:只比已测量的量,各带自己的零 —— 不写 ≥X%/≤Y%
    print(f"**修好过滤器后(n 从 {ALL['旧(笼统 x>0)']['n']['premarsx']} 涨到 {ALL['新(逐变量有效码)']['n']['premarsx']}):**")
    print(f"  政治解释 {pol_new*100:.1f}% · 宗教解释 {rel_new*100:.1f}% · 残余 {res_new*100:.1f}%")
    print(f"  残余是宗教的 {res_new/rel_new:.2f}× · 是政治的 {res_new/pol_new:.2f}×")
    print(f"  ⚠ 而 `#764` 报的是「残余是宗教的 1.89× · 是政治的 6.61×」 ——")
    print(f"     移动了 {res_new/rel_new-1.89:+.2f} 与 {res_new/pol_new-6.61:+.2f}")
json.dump(dict(old=ALL["旧(笼统 x>0)"],new=ALL["新(逐变量有效码)"],
               res_old=res_old,res_new=res_new,rel_new=rel_new,pol_new=pol_new,
               gate_ok=all(r[2] for r in G.rows)),open(OUT/"filter_ab.json","w"),ensure_ascii=False,indent=1)
