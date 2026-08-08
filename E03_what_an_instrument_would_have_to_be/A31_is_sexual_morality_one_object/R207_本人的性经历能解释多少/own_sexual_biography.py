"""#768 · E03·A40·R207 —— 管别人的性,是不是在管自己做过或没做过的事?

残余 **[49.4%, 63.3%]** 是本项目当前最大的空缺。政治只解释 [6.8%, 11.6%],宗教 [33.3%, 46.7%]。
⚠ **硬规则③禁止再从社会分类里找第四个候选** —— 那是句法,答案是语义的。
⇒ 换**一类不同的变量**:不是这个人属于哪个群体,而是**这个人自己做过什么**。

G1 估计量:`ρ(obey, 性态度 | 本人性经历)` 的保留率 —— 与政治、宗教**同一把尺子**。
控制量:`partners`(去年伴侣数 0–6)· `sexfreq`(去年性频率 0–6)· `evstray`(婚内出轨 是/否)。

⚠⚠ **两个最强混淆,都写在跑之前:**
① **对结局取条件(Oldham)**:`evstray`(出过轨)与 `xmarsex`(婚外性错不错)是**同一件事的行为面与态度面**。
   偏掉它等于偏掉结局本身。⇒ **`xmarsex` 那一格单列并标注「不可读」,主判据只用 `premarsx`/`homosex`/`teensex`。**
② **年龄**:伴侣数与频率**强烈随年龄**,而年龄也预测性态度 ⇒ 「性经历解释 X%」可能是年龄。
   ⇒ **网格里必须有一格是「性经历 + 年龄」**,并与「只年龄」对照。

⚠ **码位陷阱,由 `#766` 的 `check_kept_codes()` 在跑前摊开**(它的第三、第四次拦截):
  `partners` 码 **0 = no partners,n=9,180 有效**;码 **9 = "1 or more, (unspecified)",n=219,不在序上 ⇒ 剔除**
  `sexfreq`  码 **0 = not at all,n=8,593 有效**
  `evstray`  码 **3 = never married,n=9,509,不在 yes/no 连续统上 ⇒ 剔除**(并在 G4 里单列一格作对照)

预注册判词(按 `#764` 新写法:只比已测量的量,不写 ≥X%/≤Y%):
  与**政治的上界 11.6%** 和**宗教的下界 33.3%** 并排比,报出性经历落在哪一段,并各带自己的零。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(207)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","homosex","teensex"]; FLAG="xmarsex"
VALID={"obey":(1,5),"premarsx":(1,4),"homosex":(1,4),"teensex":(1,4),"xmarsex":(1,4),
       "partners":(0,6),"sexfreq":(0,6),"evstray":(1,2),"age":(18,89)}
print("=== #766 前瞻使用:我保留的范围排除了哪些带标签的档 ===")
for c,rng in VALID.items():
    dr,tot=check_kept_codes(gp,c,rng)
    if dr: print(f"  {c:9s} keep={rng} -> "+" · ".join(f"删 码{int(a)} {b!r} {n}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
print("  (未列出的列 = 没有排除任何带标签的档)")

d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=["obey"]+SEX+[FLAG],convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}
cats["homosex"]=cats["homosex"][:4]
flip=aligned({c:cats[c] for c in SEX+[FLAG]},"strict")|aligned({c:cats[c] for c in ["obey"]},"important")
print(f"\n方向由 `aligned()` 从值标签定 -> 要翻 {sorted(flip)}")
for c in flip: M[c]=-M[c]

def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])

BIO=["partners","sexfreq","evstray"]
SPECS={"只伴侣数":["partners"],"只性频率":["sexfreq"],"只出轨":["evstray"],
       "三个性经历一起":BIO,"只年龄":["age"],"性经历 + 年龄":BIO+["age"]}
print(f"\n=== G3 全网格 · 保留率(偏后/偏前);每格同一批行 ===")
print(f"  {'控制量':16s}"+"".join(f"{s[:9]:>11s}" for s in SEX)+f"{'⚠xmarsex':>11s}{'中位(主)':>10s}")
res={}; ns={}
for nm,cols in SPECS.items():
    row={}
    for s in SEX+[FLAG]:
        sub=M[["obey"]+BIO+["age",s]].dropna()
        raw=prho(sub.obey.to_numpy(),sub[s].to_numpy())
        v=prho(sub.obey.to_numpy(),sub[s].to_numpy(),sub[cols].to_numpy())
        floor=3*1.65/np.sqrt(len(sub))
        row[s]=dict(val=v,raw=raw,n=len(sub),keep=(v/raw if abs(raw)>=floor else None))
        ns[s]=len(sub)
    res[nm]=row
    ks=[row[s]["keep"] for s in SEX if row[s]["keep"] is not None]
    med=float(np.median(ks)) if ks else None
    res[nm]["_med"]=med
    print(f"  {nm:16s}"+"".join(f"{row[s]['keep']*100:10.1f}%" if row[s]["keep"] is not None else f"{'不可读':>11s}" for s in SEX)
          +(f"{row[FLAG]['keep']*100:10.1f}%" if row[FLAG]["keep"] is not None else f"{'不可读':>11s}")
          +(f"{med*100:9.1f}%" if med else f"{'-':>10s}"))
print(f"  {'n':16s}"+"".join(f"{ns[s]:11d}" for s in SEX)+f"{ns[FLAG]:11d}")
print(f"  ⚠ `xmarsex` 一列**不参与判词**:`evstray` 与它是同一件事的行为面与态度面 ⇒ 偏掉它是对结局取条件。")

# ---- 闸 ----
sub=M[["obey"]+BIO+["age","premarsx"]].dropna()
raw=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy())
pc=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),RNG.normal(0,1e-9,(len(sub),1)))
nul=[prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[BIO].to_numpy()[RNG.permutation(len(sub))]) for _ in range(300)]
G=Gate("#768 · 本人的性经历能解释多少")
G.identity_control("① 常数控制须回到偏前", observed=pc, expected=raw, tol=0.005, what="仪器活着吗")
G.identity_control("② 打乱三个性经历须回到偏前", observed=float(np.median(nul)), expected=raw, tol=0.005,
                   what="打乱只毁掉「谁的性经历配谁的态度」;没回到偏前就说明我毁掉的不止是配对")
# ⚠ 这个零该不该是零?不该 —— 偏掉真性经历后应当**低于偏前**,不是回到 0
G.offset_control("③ 偏掉真性经历后须显著低于偏前,才算它解释了东西",
                 effect=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[BIO].to_numpy()),
                 offset=raw, spread=float(np.std(nul)),
                 null_kind="同一批人、同一对题,唯一差别是控制量是真的三个性经历变量还是打乱后的三个")
print(); print(G)

bio=1-res["三个性经历一起"]["_med"]; age=1-res["只年龄"]["_med"]; both=1-res["性经历 + 年龄"]["_med"]
POL=(0.068,0.116); REL=(0.333,0.467)
print("\n"+"="*72)
if not all(r[2] for r in G.rows): print("**UNVERIFIED:闸没全过**")
else:
    print(f"**按 `#764` 新写法 —— 与已测量的量并排比:**")
    print(f"  政治     [{POL[0]*100:.1f}%, {POL[1]*100:.1f}%]")
    print(f"  宗教     [{REL[0]*100:.1f}%, {REL[1]*100:.1f}%]")
    print(f"  **本人性经历  {bio*100:.1f}%**   · 只年龄 {age*100:.1f}% · 性经历+年龄 {both*100:.1f}%")
    where=("**低于政治的上界**" if bio<POL[1] else "**落在政治与宗教之间**" if bio<REL[0] else "**达到或超过宗教的下界**")
    print(f"  ⇒ 本人性经历 {where}")
    print(f"  ⚠ 而它与年龄高度重叠:单独年龄已解释 {age*100:.1f}%,两者一起只到 {both*100:.1f}%"
          f" ⇒ **性经历超出年龄的增量只有 {(both-age)*100:+.1f}pp**")
json.dump(dict(n=ns,specs={k:{s:res[k][s] for s in SEX+[FLAG]} for k in SPECS},
               med={k:res[k]["_med"] for k in SPECS},bio=bio,age=age,both=both,
               pol_ref=POL,rel_ref=REL,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"own_biography.json","w"),ensure_ascii=False,indent=1)
