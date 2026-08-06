"""#760 · E03·A39·R201 —— 把政治意识形态偏掉之后,「性绑在权威上」还剩多少?

⚠⚠ **这一步是故意挑的、正面结果我不欢迎的一步(basin 规则)。**
`#750`→`#759` 一路下来,三具仪器同向支持「性绑在权威/服从上,独立于伤害-公平」。
**而这条弧的最强反对解释是:一条左右政治轴同时造出了 +0.5001(权威)与 −0.1685(公平)。**
若是这样,「性绑在权威上」就压缩成「性绑在保守主义上」—— **一句已知的旧话(prior_art)**。
**B 世界成立会毁掉我这条弧的头号发现,所以它必须先被测。**

G1 估计量(先于方法):
  `ρ(性条目, 道德维度 | 政治意识形态)` —— 偏掉政治之后,每个道德维度还剩多少。
  MFQ 臂:chastity × {伤害 · 公平 · 内群体 · 权威},控制 `politics`(1–7)
  GSS 臂:四题性 × 五个育儿价值,控制 `polviews`(1–7)

⚠ **跑之前写下的最强混淆,并且量过了**:
  MFQ 是 YourMorals.org 自选样本,`politics` 均值 **2.884**、最保守一格只有 139 人 ——
  **政治的全距受限 ⇒ 偏掉它会校正不足 ⇒ 方向偏向「A 存活」,即偏向我想要的结论。**
  控制:**GSS `polviews` 是概率样本**(均值≈4,七档 2221/8044/8268/25140/9977/9877/2351,分布饱满),
  **同一操作在它上面再做一次** —— 两具仪器不一致时,以没有全距受限的那一具为准。

预注册判词(阈值写在跑之前;**分支必须覆盖整个区间**,`#754` 那一族):
  W-B **政治排序**:GSS 臂上 `obey` 偏后保留 **≤50%**  ⇒ 头号声明降级为「性绑在保守主义上」
  W-A **道德结构**:GSS 臂上 `obey` 偏后保留 **≥75%**  ⇒ 权威不是政治的影子
  W-M 之间(50%–75%)⇒ **判不了,并打印落在哪一段** —— 不许事后改阈值
  ⚠ 另记(不参与判词,因为它是本轮的第二问):**公平/助人那一格偏后是零还是负。**

⚠ `#759`① 说 `aligned()` 还没有任何调用者,而没人用的工具与没造一样 —— **本轮真的用它**。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp, aligned, label_pole
RNG=np.random.default_rng(201)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)

def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])

# ============================ GSS 臂(概率样本,判词判在这里)============================
K=["obey","thnkself","workhard","helpoth","popular"]; SEX=["premarsx","xmarsex","homosex","teensex"]
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
cat=pd.read_stata(gp,columns=K+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}
print("=== #759① · 真的调用 aligned():方向从值标签读,不从记忆 ===")
try:
    aligned({c:cats[c] for c in SEX},"strict")
except ValueError as e:
    print(f"  ⚠ 未剔除离表码时 aligned 直接 raise:{str(e)[:66]}")
cats["homosex"]=cats["homosex"][:4]           # 剔除 'other'(码 5)
flip_sex=aligned({c:cats[c] for c in SEX},"strict")
flip_val=aligned({c:cats[c] for c in K},"important")
print(f"  四题性统一成「高=严」 -> 要翻 {sorted(flip_sex)}")
print(f"  五育儿价值统一成「高=看重」-> 要翻 {sorted(flip_val)}")
for c in K: print(f"    {c:9s} 高值端 = {label_pole(cats[c])}")

g=pd.read_stata(gp,columns=["year","polviews"]+K+SEX,convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(g[c],errors="coerce").where(lambda x:x>0) for c in ["polviews"]+K+SEX})
M["homosex"]=M["homosex"].where(M["homosex"]<=4); M["year"]=g.year
for c in (flip_sex|flip_val): M[c]=-M[c]        # ⚠ 翻的动作写在调用者这里,库不改数据
print(f"\n  ⚠ polviews 高值端(GSS 编码 1=extremely liberal … 7=extremely conservative)= 保守")

print("\n=== GSS 臂:偏掉 polviews 之前 / 之后 ===")
G_rows={}
for s in SEX:
    sub=M[K+[s,"polviews","year"]].dropna()
    row={}
    for k in K:
        b=prho(sub[k].to_numpy(),sub[s].to_numpy())
        a=prho(sub[k].to_numpy(),sub[s].to_numpy(),sub[["polviews"]].to_numpy())
        # ⚠ #745 那一族：分母近零时的比值不可读。判据不是 1e-9，是**偏前值自己的噪声**：
        #   n≈16,000 时秩相关的零 95% 分位约 1.65/sqrt(n) ≈ 0.013；偏前 |b| ≥ 3× 它才算有分母。
        floor=3*1.65/np.sqrt(len(sub))
        row[k]=dict(before=b,after=a,keep=(a/b if abs(b)>=floor else None),
                    keep_inadmissible=bool(abs(b)<floor), floor=float(floor))
    row["_n"]=len(sub); G_rows[s]=row
hdr=f"  {'':10s}"+"".join(f"{s[:9]:>22s}" for s in SEX); print(hdr)
print(f"  {'':10s}"+"".join(f"{'偏前':>7s}{'偏后':>8s}{'保留':>7s}" for _ in SEX))
for k in K:
    print(f"  {k:10s}"+"".join(
        f"{G_rows[s][k]['before']:+7.3f}{G_rows[s][k]['after']:+8.3f}"
        + (f"{G_rows[s][k]['keep']*100:6.0f}%" if G_rows[s][k]['keep'] is not None else f"{chr(19981)+chr(21487)+chr(35835):>7s}")
        for s in SEX))
print(f"  ⚠ 「不可读」= 偏前值本身小于 3× 它自己的零（约 {3*1.65/np.sqrt(15571):.4f}），**比值的分母近零 ⇒ 不许报**（`#745`）")
print(f"  {'n':10s}"+"".join(f"{G_rows[s]['_n']:22d}" for s in SEX))
obey_keep=float(np.median([G_rows[s]["obey"]["keep"] for s in SEX]))
help_after=[G_rows[s]["helpoth"]["after"] for s in SEX]
print(f"\n  obey 保留率中位数 = {obey_keep*100:.1f}%")
print(f"  helpoth 偏后 = {['%+.4f'%x for x in help_after]}")

# --- G4 规格曲线:政治当 7 个哑变量(非线性)· 年内 ---
print("\n=== G4 规格曲线(GSS 臂,obey × premarsx)===")
sub=M[K+["premarsx","polviews","year"]].dropna()
specs={}
specs["线性 polviews"]=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),sub[["polviews"]].to_numpy())
D7=pd.get_dummies(sub.polviews.astype(int),drop_first=True).to_numpy(float)
specs["polviews 7 档哑变量"]=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),D7)
DY=pd.get_dummies(sub.year.astype(int),drop_first=True).to_numpy(float)
specs["polviews + 年份哑变量"]=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),np.c_[D7,DY])
raw=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy())
for nm,v in specs.items(): print(f"  {nm:22s} {v:+.4f}  保留 {v/raw*100:5.1f}%   (偏前 {raw:+.4f})")

# --- G2 控制 ---
print("\n=== G2 控制(GSS 臂)===")
z=RNG.normal(0,1e-9,(len(sub),1))
pc=prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),z)
print(f"  正控:控制量为常数须退回偏前 {raw:+.4f},实测 {pc:+.4f} ⇒ {'通过' if abs(pc-raw)<0.005 else '**不通过**'}")
nul=[prho(sub.obey.to_numpy(),sub.premarsx.to_numpy(),RNG.permutation(sub.polviews.to_numpy()).reshape(-1,1)) for _ in range(300)]
print(f"  阴性:打乱 polviews 与人的配对后,偏相关应回到偏前 —— 实测中位 {np.median(nul):+.4f}(偏前 {raw:+.4f})")
print(f"     ⇒ **这个零该不该是零?不该** —— 打乱后偏相关**必须回到偏前**,所以用的是 offset_control 的逻辑:")
print(f"       零的种类 = 在人之间打乱 `polviews`,保住它的边际,只毁掉「谁的政治配谁的态度」")

# ============================ MFQ 臂(自选样本,全距受限)============================
print("\n=== MFQ 臂(⚠ 政治全距受限,只作参照)===")
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
F=["HARM_AVG","FAIRNESS_AVG","INGROUP_AVG","AUTHORITY_AVG"]
Md=d[["chastity","politics"]+F].apply(pd.to_numeric,errors="coerce").dropna()
print(f"  n={len(Md)} · politics 均值 {Md.politics.mean():.3f} sd {Md.politics.std():.3f}"
      f"  ⚠ 对照 GSS polviews sd {M.polviews.std():.3f} ⇒ **全距比 {Md.politics.std()/M.polviews.std():.2f}×**")
Mrows={}
for f in F:
    b=prho(Md.chastity.to_numpy(),Md[f].to_numpy())
    a=prho(Md.chastity.to_numpy(),Md[f].to_numpy(),Md[["politics"]].to_numpy())
    fl=3*1.65/np.sqrt(len(Md)); kp=(a/b) if abs(b)>=fl else None
    Mrows[f]=dict(before=b,after=a,keep=kp,keep_inadmissible=kp is None,floor=float(fl))
    print(f"  {f:14s} 偏前 {b:+.4f} -> 偏后 {a:+.4f}  保留 "
          + (f"{kp*100:6.1f}%" if kp is not None else "不可读（分母近零）"))

from lib.gates import Gate
GG=Gate("#760 · 把政治偏掉之后权威还剩多少")
# ⚠ 先问：**这个零该不该是零？** 打乱 polviews 之后偏相关**必须回到偏前**，不是回到 0
#   ⇒ 所以是 offset_control（系统性基线偏移），不是 negative_control。
# ⚠⚠ 第一版把这一条写成了 offset_control —— 而 offset_control 测的是「差要够大」，
#    我要断言的却是「打乱后必须**回到**偏前」，即**相等**。**用差值检测器去断言等式**，
#    是这几轮反复犯的同一族（`#728`·`#748`·`#750`·`#758`），这次犯在闸的选型上。
#    ⇒ 这一条是**仪器检查**，不是零；它属于 positive_control 的形状。
# ⚠⚠ 第一版把这两条写成了 offset_control —— 而 offset_control 测的是「差要够大」，
#    我要断言的却是「换掉控制量后偏相关必须**回到**偏前」，即**相等**。
#    **用差值检测器去断言等式**，判词当场 FAIL 在一个完全正常的仪器上 ——
#    这是 `#728`·`#748`·`#750`·`#758` 那一族的第七次，这次犯在**闸的选型**上。
#    ⇒ 库里本来就没有「这两个量必须相等」的形状。`#761` 补上了 `identity_control`，本轮是它的第一个调用者。
GG.identity_control("① 常数控制须回到偏前（仪器活着吗）",
                    observed=pc, expected=raw, tol=0.005,
                    what="控制量取常数时偏出来的东西是零，偏相关必须原样等于偏前")
GG.identity_control("② 打乱 polviews 须回到偏前（毁的是配对，不是 polviews 本身）",
                    observed=float(np.median(nul)), expected=raw, tol=0.005,
                    what="打乱只毁掉「谁的政治配谁的态度」；若偏相关没回到偏前，说明我毁掉的不止是配对")
# ⚠ 而「偏掉真政治之后还剩多少」这一格，零**应该**是零吗？不 —— 它应该是**偏前值**。
#   若拿它去比 0，任何非零都会被判成「效应还在」，那是不会失败的检查。
GG.offset_control("偏掉真 polviews 后的 obey×premarsx，须显著低于偏前才算政治解释了东西",
                  effect=specs["线性 polviews"], offset=raw, spread=float(np.std(nul)),
                  null_kind="同一批人、同一对题，唯一差别是控制量是真 `polviews` 还是打乱后的 `polviews`")
print(GG)

print("\n"+"="*70)
if abs(pc-raw)>=0.005: v="**UNVERIFIED:正控没过,偏相关这具仪器本身不可信**"
elif obey_keep<=0.50:
    v=(f"**W-B 政治排序:GSS 概率样本上 `obey` 偏掉政治后只保留 {obey_keep*100:.1f}% "
       f"⇒ 「性绑在权威上」降级为「性绑在保守主义上」,而那是已知的**")
elif obey_keep>=0.75:
    v=(f"**W-A 道德结构:GSS 概率样本上 `obey` 偏掉政治后仍保留 {obey_keep*100:.1f}% "
       f"⇒ 权威不是政治的影子**")
else:
    v=(f"**W-M:`obey` 保留 {obey_keep*100:.1f}%,落在预注册的 50%–75% 之间 —— 判不了,不改阈值**")
print(v)
json.dump(dict(gss=G_rows,gss_specs=specs,gss_raw=raw,obey_keep=obey_keep,
               mfq={k:vv for k,vv in Mrows.items()},mfq_n=len(Md),
               sd_ratio=float(Md.politics.std()/M.polviews.std()),
               pos_control=pc,perm_median=float(np.median(nul)),verdict=v),
          open(OUT/"is_it_just_politics.json","w"),ensure_ascii=False,indent=1)
