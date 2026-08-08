"""#762 · E03·A39·R202 —— 宗教能不能解释掉「服从」?

`#760` 量到:偏掉政治左右后,`obey` 与管性的关系还剩 **87.7%**,政治只解释 12.1%。
`#760`① 预注册的下一问:**那剩下的九成是什么?最强候选是宗教** ——
而 `#747` 已量到常去教堂的人在同性恋态度上的变动只有别人的一半。

⚠⚠ **我不欢迎的正面结果是「宗教解释掉它」** —— 那会把发现压成「性绑在虔诚上」,同样是旧话。
**预注册的强结论则是它的反面:若宗教也只解释一成,那「服从」就不是任何一个已知社会分类的影子。**

G1 估计量:`ρ(obey, 性条目 | 控制量)`,控制量依次为
  空 · 政治 · 出席频率 · 宗教认同强度 · 原教旨-自由 · 三个宗教量一起 · 宗教+政治全放。
**全部在同一批行上算**(n=9,580–10,467),否则比的是不同的人。

⚠ **共线性预检(跑之前,`#750`① 的规矩,|ρ|>0.8 ⇒ 判不了)**:
  ρ(obey, attend) = −0.1980 · ρ(obey, reliten) = +0.1629 · ρ(obey, fund) = +0.2328 · ρ(obey, polviews) = −0.1191
  **最大 0.2328,远低于 0.8 ⇒ 偏掉是有意义的操作,不是把自己减自己。**

预注册判词(阈值写在跑之前;**分支覆盖整个区间**):
  W-R 宗教解释掉它:三个宗教量一起偏掉后,`obey` 保留 **≤50%** ⇒ 「服从」是虔诚的影子
  W-O 服从是它自己:保留 **≥75%** ⇒ **它不是任何一个已知社会分类的影子**(强结论)
  W-M 之间(50%–75%)⇒ 判不了,打印落在哪一段,**不许事后改阈值**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, label_pole
from lib.gates import Gate
RNG=np.random.default_rng(202)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
K=["obey"]; SEX=["premarsx","xmarsex","homosex","teensex"]; REL=["attend","reliten","fund"]

cat=pd.read_stata(gp,columns=K+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
flip=aligned({c:cats[c] for c in SEX},"strict") | aligned({c:cats[c] for c in K},"important")
print(f"=== 方向由值标签决定(`aligned()`,#759①)-> 要翻 {sorted(flip)} ===")

d=pd.read_stata(gp,columns=["year","polviews"]+REL+K+SEX,convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(lambda x:x>0) for c in ["polviews"]+REL+K+SEX})
M["homosex"]=M["homosex"].where(M["homosex"]<=4); M["year"]=d.year
for c in flip: M[c]=-M[c]

def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])

SPECS={"(空,偏前)":[], "政治 polviews":["polviews"], "出席 attend":["attend"],
       "认同强度 reliten":["reliten"], "原教旨 fund":["fund"],
       "三个宗教量一起":REL, "宗教 + 政治全放":REL+["polviews"]}
res={}; keeps={}
print(f"\n=== G3 全网格(同一批行;偏前在第一列)===")
print(f"  {'控制量':20s}"+"".join(f"{s[:9]:>12s}" for s in SEX)+f"{'保留中位':>10s}")
for nm,cols in SPECS.items():
    row={}; ks=[]
    for s in SEX:
        sub=M[K+REL+["polviews",s]].dropna()
        raw=prho(sub.obey.to_numpy(),sub[s].to_numpy())
        v=raw if not cols else prho(sub.obey.to_numpy(),sub[s].to_numpy(),sub[cols].to_numpy())
        floor=3*1.65/np.sqrt(len(sub))      # #745:近零分母的比值不可读
        row[s]=dict(val=v,raw=raw,n=len(sub),keep=(v/raw if abs(raw)>=floor else None))
        if row[s]["keep"] is not None: ks.append(row[s]["keep"])
    res[nm]=row; keeps[nm]=float(np.median(ks)) if ks else None
    print(f"  {nm:20s}"+"".join(f"{row[s]['val']:+12.4f}" for s in SEX)
          + (f"{keeps[nm]*100:9.1f}%" if keeps[nm] is not None else f"{'不可读':>10s}"))
print(f"  {'n':20s}"+"".join(f"{res['(空,偏前)'][s]['n']:12d}" for s in SEX))

# ---- 闸:两条 identity_control + 一条 offset_control ----
sub=M[K+REL+["polviews","premarsx"]].dropna()
raw=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy())
pc=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),RNG.normal(0,1e-9,(len(sub),1)))
nul=[prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),
          sub[REL].to_numpy()[RNG.permutation(len(sub))]) for _ in range(300)]
rel_only=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[REL].to_numpy())
G=Gate("#762 · 宗教能不能解释掉服从")
G.identity_control("① 常数控制须回到偏前(仪器活着吗)", observed=pc, expected=raw, tol=0.005,
                   what="控制量取常数时偏出来的东西是零,偏相关必须原样等于偏前")
G.identity_control("② 打乱三个宗教量须回到偏前(毁的是配对)", observed=float(np.median(nul)), expected=raw, tol=0.005,
                   what="打乱只毁掉「谁的宗教配谁的态度」;若没回到偏前,说明我毁掉的不止是配对")
# ⚠ 这个零该不该是零?不该 —— 偏掉真宗教之后应当**低于偏前**,而不是回到 0。
G.offset_control("③ 偏掉真宗教后须显著低于偏前,才算宗教解释了东西",
                 effect=rel_only, offset=raw, spread=float(np.std(nul)),
                 null_kind="同一批人、同一对题,唯一差别是控制量是真的三个宗教量还是打乱后的三个宗教量")
print(); print(G)

k=keeps["三个宗教量一起"]
print("\n"+"="*70)
gate_ok=all(r[2] for r in G.rows)
if not gate_ok: v="**UNVERIFIED:闸没全过**"
elif k is None: v="**UNVERIFIED:偏前值近零,保留率不可读**"
elif k<=0.50: v=f"**W-R:三个宗教量一起偏掉后 `obey` 只保留 {k*100:.1f}% ⇒ 「服从」是虔诚的影子**"
elif k>=0.75: v=(f"**W-O:三个宗教量一起偏掉后 `obey` 仍保留 {k*100:.1f}%,再加政治也只到 "
                 f"{keeps['宗教 + 政治全放']*100:.1f}% ⇒ 「服从」不是任何一个已知社会分类的影子**")
else: v=f"**W-M:保留 {k*100:.1f}%,落在预注册的 50%–75% 之间 —— 判不了,不改阈值**"
print(v)
json.dump(dict(specs={n:{s:res[n][s] for s in SEX} for n in SPECS},keeps=keeps,
               raw=raw,pos_control=pc,perm_median=float(np.median(nul)),rel_only=rel_only,
               gate_ok=bool(gate_ok),verdict=v),open(OUT/"is_it_religion.json","w"),ensure_ascii=False,indent=1)
