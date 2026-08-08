"""#751 · E03·A38·R194 —— 一个人对性的严,是不是他一般有多严的一个侧面?

⚠ 这一步是**故意选的、正面结果我不欢迎的一步**(basin 规则)。
   `#749` 与 `#750` 连着两轮确认「性在社会这个单位上自成一条线」。两轮同向 = basin。
   所以本轮问它在**人**这个单位上成不成立,而**我不欢迎的答案(W2:人层上性只是一般道德严厉的一个侧面)
   会把前两轮的话压缩成「只在社会层成立」** —— 那是本体论上的变化,不是参数上的。

G1 估计量(先于方法):
   S = 四题性的 z 均值(premarsx · xmarsex · homosex · teensex)
   G = 两题非性道德的 z 均值(taxcheat 逃税 · govcheat 骗福利)—— 同一份问卷、同一套「错不错」尺度
   ① 收敛那一半:ρ(S, G) —— 性的严厉里有多少就是一般道德的严厉
   ② 判别那一半:**偏掉 G 之后,四题性彼此之间还剩多少一致性**
      —— 若性只是一般严厉的一个侧面,偏掉 G 后它们应当散开

⚠⚠ 本轮真正的控制,针对的是**会替我造出想要答案的那个机制**:
   G 只有两题,信度低;**偏掉一个不可靠的控制量必然校正不足**,而校正不足**偏向我想要的结论**。
   所以除了原始偏相关,还算**信度校正后的上界**(Spearman-Brown + 去衰减):
   若 G 完美可靠,偏掉它之后性还剩多少。**W2 用的是这个上界那一侧,不是原始值。**

硬规则①(变量名不是测量):所有列的 n 与真正问过的年份都在下面打印,不引用记忆。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import spearman as sp

RNG=np.random.default_rng(194)
ROOT=pathlib.Path(__file__).resolve().parents[3]
P=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]; NON=["taxcheat","govcheat"]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)

d=pd.read_stata(P,columns=["year"]+SEX+NON,convert_categoricals=False)
print("=== 硬规则①:变量名不是测量 ===")
raw={}
for c in SEX+NON:
    v=pd.to_numeric(d[c],errors="coerce").where(lambda x:x>0)
    # ⚠ homosex 有第 5 个码。GSS 的 1..4 是「总是错→完全不错」,5 不在这条量表上 —— 必须剔除。
    if c=="homosex":
        n5=int((v==5).sum()); v=v.where(v<=4); print(f"  ⚠ homosex 码 5 共 {n5} 例,不在 1..4 的错度量表上,剔除")
    yrs=sorted(d.loc[v.notna(),"year"].unique().astype(int))
    raw[c]=v; print(f"  {c:9s} n={int(v.notna().sum()):6d} 年份 {len(yrs)} 个:{yrs[0]}–{yrs[-1]}")
X=pd.DataFrame(raw); X["year"]=d.year
D=X.dropna().copy()
print(f"\n完整个案 n={len(D)} · 年份 {sorted(D.year.unique().astype(int))}")
print("  ⚠ 两条非性条目只在 1991 与 1998 问过 —— **这不是 GSS 的 4.6 万,是 734。**")

z=lambda s:(s-s.mean())/s.std(ddof=1)
# ⚠⚠ 编码方向:从 .dta 的值标签读出来的,不是记忆(`#734` 那一族,这是第四次)。
#   性题   1=always wrong → 4=not wrong at all   ⇒ 高 = **宽容**
#   非性题 1=not wrong    → 4=seriously wrong    ⇒ 高 = **严**
#   两套方向相反。统一成「高 = 严」:性题取负。
D["S"]=-z(D[SEX]).mean(axis=1); D["G"]=z(D[NON]).mean(axis=1)
print("\n  ⚠ 值标签实读:性题 高=宽容 · 非性题 高=严 —— 方向相反,性题已取负统一成「高=严」")

def rho(a,b): return sp(pd.Series(a),pd.Series(b))
def resid(y,X_):
    X_=np.c_[np.ones(len(X_)),X_]; return y-X_@np.linalg.lstsq(X_,y,rcond=None)[0]
def prho(a,b,c):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])

# --- 信度:两题的 Spearman-Brown,四题的 Cronbach α ---
def sb2(a,b):
    r=rho(a,b); return 2*r/(1+r)
def alpha(df):
    k=df.shape[1]; v=df.var(ddof=1); t=df.sum(axis=1).var(ddof=1)
    return k/(k-1)*(1-v.sum()/t)
relG=sb2(D.taxcheat,D.govcheat); relS=alpha(z(D[SEX]))
print(f"\n=== 信度(决定了偏掉 G 会校正不足多少)===")
print(f"  G(两题非性)Spearman-Brown = {relG:.4f}   ⚠ 越低 = 偏掉它校正得越不够 = 越偏向我想要的答案")
print(f"  S(四题性)Cronbach α      = {relS:.4f}")

# --- ① 收敛那一半 ---
rSG=rho(D.S,D.G)
print(f"\n=== ① 收敛:ρ(性的严, 一般道德的严) = {rSG:+.4f} ===")
print(f"  去衰减后 = {rSG/np.sqrt(relG*relS):+.4f}(两个合成量各自的信度都扣掉)")

# --- ② 判别:偏掉 G 之后四题性还剩多少一致性 ---
import itertools
pairs=list(itertools.combinations(SEX,2))
before={f"{a}×{b}":rho(D[a],D[b]) for a,b in pairs}
after ={f"{a}×{b}":prho(D[a].to_numpy(),D[b].to_numpy(),D.G.to_numpy()) for a,b in pairs}
print(f"\n=== ② 判别:偏掉 G 之后,四题性彼此还剩多少 ===")
print(f"  {'对':26s} {'偏前':>8s} {'偏后':>8s} {'保留':>7s}")
for k in before:
    print(f"  {k:26s} {before[k]:+8.4f} {after[k]:+8.4f} {after[k]/before[k]*100:6.1f}%")
mb=float(np.median(list(before.values()))); ma=float(np.median(list(after.values())))
print(f"  {'中位数':26s} {mb:+8.4f} {ma:+8.4f} {ma/mb*100:6.1f}%")

# ⚠ 信度校正的上界:若 G 完美可靠,偏掉它之后还剩多少(用去衰减后的 G 做控制)
#   做法:把 G 的秩残差按 1/sqrt(relG) 放大不可行(偏相关对控制量的线性缩放不变),
#   所以走公式:完美可靠控制下的偏相关 = (r_ab - r_aG*r_bG/relG) / sqrt((1-r_aG^2/relG)(1-r_bG^2/relG))
def prho_corrected(a,b,g,rel):
    rab,rag,rbg=rho(a,b),rho(a,g),rho(b,g)
    rag2,rbg2=rag**2/rel,rbg**2/rel
    den=np.sqrt(max(1-rag2,1e-12)*max(1-rbg2,1e-12))
    return float((rab-rag*rbg/rel)/den)
corr={f"{a}×{b}":prho_corrected(D[a],D[b],D.G,relG) for a,b in pairs}
mc=float(np.median(list(corr.values())))
print(f"\n=== ②b 若 G 完美可靠(信度校正上界)—— **W2 判在这一侧** ===")
for k in corr: print(f"  {k:26s} {corr[k]:+8.4f}")
print(f"  {'中位数':26s} {mc:+8.4f}  ⇒ 相对偏前保留 {mc/mb*100:.1f}%")

# --- 零:打乱 G 与人的配对,重算偏相关(保住 G 的边际与四题性之间的一切) ---
NP=2000
nul=[]
for _ in range(NP):
    gp=RNG.permutation(D.G.to_numpy())
    nul.append(float(np.median([prho(D[a].to_numpy(),D[b].to_numpy(),gp) for a,b in pairs])))
nq=(float(np.quantile(nul,.025)),float(np.quantile(nul,.975)))
print(f"\n=== 零:打乱「谁的一般道德配谁的性态度」,{NP} 次,重算偏相关 ===")
print(f"  零的 95% 区间 {nq[0]:+.4f} … {nq[1]:+.4f}(中位数 {np.median(nul):+.4f})")
print(f"  ⚠ 这个零检验的是**偏掉 G 有没有起作用**:若偏后 ≈ 零区间中心,说明 G 根本没在动它")

# --- 正控:G=0 时必须不通过(闸必须能失败) ---
gz=np.zeros(len(D))
pc=float(np.median([prho(D[a].to_numpy(),D[b].to_numpy(),gz+RNG.normal(0,1e-9,len(D))) for a,b in pairs]))
print(f"\n=== 正控:控制量为常数(g=0)时,偏相关必须退回偏前 {mb:+.4f} ===")
print(f"  实测 {pc:+.4f} · 差 {abs(pc-mb):.5f} ⇒ {'通过(闸能失败)' if abs(pc-mb)<0.01 else '不通过'}")

# --- 预注册判词 ---
print("\n"+"="*64)
if abs(pc-mb)>=0.01: v="**UNVERIFIED:正控没过,偏相关这具仪器本身不可信**"
elif not (nq[0]<mb<nq[1]) and mc<=nq[1]: v=f"**W2:信度校正后中位数 {mc:+.4f} 落进零区间 ⇒ 人这个单位上,性只是一般道德严厉的一个侧面**"
elif mc>=0.5*mb: v=(f"**W1:即使假设 G 完美可靠,四题性彼此仍保留 {mc/mb*100:.1f}%(中位数 {mc:+.4f},零上界 {nq[1]:+.4f})"
                    f" ⇒ 人这个单位上性也自成一条线**")
else: v=f"**W3:偏后 {ma:+.4f}、校正上界 {mc:+.4f} 介于两者之间 —— 这个设计判不了**"
print(v)
json.dump(dict(n=len(D),years=sorted(D.year.unique().astype(int).tolist()),rel_G=relG,rel_S=relS,
               rho_SG=rSG,rho_SG_disatt=rSG/np.sqrt(relG*relS),before=before,after=after,
               corrected=corr,med_before=mb,med_after=ma,med_corrected=mc,null_ci=nq,
               pos_control=pc,verdict=v),
          open(OUT/"person_level.json","w"),ensure_ascii=False,indent=1)
