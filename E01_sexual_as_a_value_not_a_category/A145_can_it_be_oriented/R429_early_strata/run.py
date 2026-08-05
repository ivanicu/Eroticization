import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A145 R429 -- 那条「越羞耻越少实践」的关系,是不是起始时间造出来的

`#384b`:`羞耻 → ACTED` = **−0.0833**,−6.85 sd。**但 `#384f④` 说方向不可识别。**
这份数据里有一个**时间上在前**的东西:`EARLY`(起始年龄)——
**一个人 9 岁时的起始,不可能由他今天的实践量造成。**

⚠⚠ **这不是工具变量,而且我不会把它写成工具变量。**
`EARLY` 与羞耻**直接相关**(`#332`:−0.102),排除性限制**不成立**。
**本轮做的是「部分定向」:问那条关系在 `EARLY` 的分层内部是否仍然成立、大小是否相近。**
若成立且相近 -> **它不是「起始时间造出来的」**。这仍然是关联,只是**少了一个候选解释**。
(§2 的 η:把这写成因果就是步子迈过了梯度。)

ESTIMAND        按 `EARLY` 分四层,层内 `ACTED ~ 羞耻 + S + c3⁻ + 类别数`;
                主量 ① 每层的羞耻系数;② **四层的全距**(是否随起始时间变)。
判据(**先标支**,`#379c`)
                【两支】guard 22(先证明是曲线)· guard 26(**用 MDE 扫描当正对照**,`#384d`)· 负对照。
                【非零支】四层**全部同号且各自越过零** -> 那条关系在早晚两端都成立。
                【零支】仅当某层未越阈时启用 MDE。
                另判:**全距是否越过随机分层的零**(= 关系随起始时间变)。
⚠ 零的种类     `offset_control`:**全距的零绝不是零**(任意四层都有差)->
                零 = **随机等大小分层**(层大小照旧,人打乱)后的全距分布。
IMPOSSIBLE      ① `EARLY` 与羞耻相关 -> 分层会**顺带**在羞耻上分层,层内羞耻的方差被压窄 ->
                   **层内系数天然比全样本小**,所以**只比符号与是否越阈,不比与 −0.0833 的大小**;
                ② 四层各 ~1,700 人 -> 功率低,先算 MDE;
                ③ 「不是起始时间造出来的」≠「是羞耻造成的」——**只少了一个候选解释**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
AC=next(c for c in d.columns if '41kpfir' in c)
ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ACTED)&np.isfinite(ncat)&np.isfinite(EARLY)
n=int(M.sum()); idx=np.flatnonzero(M)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def coef(y,g):
    k=int(g.sum())
    X=np.column_stack([np.ones(k),z(sh,g),z(S,g),z(C3,g),z(ncat,g)])
    yy=z(y,g); b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-5); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1])
QS=np.quantile(EARLY[M],[0,.25,.5,.75,1.]); QS[-1]+=1e-9
LAY=[idx[(EARLY[idx]>=QS[i])&(EARLY[idx]<QS[i+1])] for i in range(4)]
def gof(j):
    g=np.zeros(NN,bool); g[j]=True; return g
print(f"n=**{n:,}** · 按 `EARLY` 四分层(低 = 来得早)")
print(f"全样本(参照,`#384b`):**{coef(ACTED,M)[0]:+.4f}**\n")
rows=[]
for i,j in enumerate(LAY):
    g=gof(j); b,se=coef(ACTED,g)
    rows.append(dict(v_layer=i+1,v_n=len(j),v_early=float(np.median(EARLY[j])),v_b=b,v_se=se))
    print(f"   层{i+1}  n={len(j):>5,} · 起始中位 **{np.median(EARLY[j]):.1f}岁** · "
          f"羞耻→ACTED **{b:+.4f}** (se {se:.4f})")
T=pd.DataFrame(rows); check_columns(T,'R429')
RNG=float(T.v_b.max()-T.v_b.min())
print(f"\n四层全距 = **{RNG:.4f}**")

# ⚠ **第一版这里坏了,而两个门同时 FAIL 是同一个原因**:
# 我写 `ACTED[rg.permutation(NN)]` —— **打乱了整个数组(含 NaN)**,于是掩码内出现 NaN,
# `z()` 返回 NaN,阈变成 **NaN**,所有比较恒为 False,MDE 扫描于是 0/30。
# **一个 NaN 阈会让每一个门都「失败」,看起来像四个问题,其实是一个。**
# 修法:项目已有的 `perm_finite` 模式 —— **只在有限项内部打乱,保住缺失格局**;
# 而且层内的零必须**在层内**打乱,不能借全样本的。
def perm_in(v,g,seed):
    z_=v.copy(); jj=np.flatnonzero(g&np.isfinite(v))
    z_[jj]=v[jj][np.random.default_rng(seed).permutation(len(jj))]; return z_
NP_=400
nul_b=[]; nul_rng=[]
sizes=[len(j) for j in LAY]
gsm0=gof(LAY[0])
for s_ in range(NP_):
    rg=np.random.default_rng(8800+s_); p=rg.permutation(idx); c=0; bs=[]
    for k in sizes:
        bs.append(coef(ACTED,gof(p[c:c+k]))[0]); c+=k
    nul_rng.append(max(bs)-min(bs))
    nul_b.append(coef(perm_in(ACTED,gsm0,9900+s_),gsm0)[0])   # ★ 层内打乱,层内大小
nul_rng=np.array(nul_rng); nul_b=np.array(nul_b)
THR_B=float(np.percentile(np.abs(nul_b),95)); THR_R=float(np.percentile(nul_rng,95))
print(f"⚠ offset 零 ①(**层内**打乱 `ACTED`,层内 n={int(gsm0.sum()):,})|值| 95 分位 = **{THR_B:.4f}**")
print(f"⚠ offset 零 ②(**随机等大小分层**的全距;**任意四层都有差,所以零不是零**)= "
      f"**{nul_rng.mean():.4f} ± {nul_rng.std():.4f}** · 95 分位 **{THR_R:.4f}**")
OVER=[abs(b)>THR_B for b in T.v_b]; SAMES=bool((T.v_b<0).all())
print(f"   -> 四层越阈 **{sum(OVER)}/4** · 全部同号(负)**{SAMES}** · "
      f"全距 {RNG:.4f} vs 阈 {THR_R:.4f} -> "
      f"{'**越阈:关系随起始时间变**' if RNG>THR_R else '**未越阈:四层读作同一个数**'}")

print(f"\nguard 26 的正对照 = **MDE 扫描**(`#384d` 的规矩),每级 30 次,种在**最小的一层**上:")
gsm=gof(LAY[int(np.argmin(sizes))]); nsm=int(gsm.sum()); MDE=None
for gg in (0.04,0.06,0.08,0.12):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(5150+int(gg*100)*23+s_)
        y=np.full(NN,np.nan); y[gsm]=-gg*z(sh,gsm)+rg.standard_normal(nsm)
        if abs(coef(y,gsm)[0])>THR_B: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.15
print(f"   **MDE(最小层)= {MDE_:.2f}** · 争议幅度 = 全样本的 |−0.0833| = **0.083**")
bNEG,_=coef(perm_in(ACTED,gsm0,77),gsm0)   # ★ 同样在层内
T['v_thr']=THR_B; T.to_csv(pathlib.Path(__file__).parent/'results'/'strata.csv',index=False)

g=Gate('那条「越羞耻越少实践」的关系是不是起始时间造出来的')
g.asserted('★【两支】guard 22:先证明它是一条曲线(4 个不同的起始中位)',
           len(set(T.v_early.round(3)))>=3,f"{len(set(T.v_early.round(3)))} 个不同的 x",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:正对照用 **MDE 扫描**,不是单次种植',MDE_,0.0833,True,
    what='最小层的 80% 检出幅度 vs 全样本效应')
g.asserted('★【两支】负对照:打乱 `ACTED` -> 必须落回零',abs(bNEG)<=THR_B,
           f"{bNEG:+.5f} vs {THR_B:.4f}",kind='control')
g.asserted('★【两支】offset 零 ② 非退化(任意四层都有差)',nul_rng.std()>0,
           f"{nul_rng.mean():.4f} ± {nul_rng.std():.4f}",kind='control')
if MDE_<=0.0833 and abs(bNEG)<=THR_B:
    g.asserted('★【非零支】四层全部同号(负)且各自越过零 -> 不是起始时间造出来的',
               SAMES and all(OVER),f"同号 {SAMES} · 越阈 {sum(OVER)}/4 · {[round(x,4) for x in T.v_b]}")
    g.asserted('★【非零支/另判】全距未越过随机分层的零 -> 关系不随起始时间变',
               RNG<=THR_R,f"全距 {RNG:.4f} vs 阈 {THR_R:.4f}")
else:
    g.asserted('★ 正对照标定不足或负对照未过 -> 不判',False,
               f"MDE {MDE_:.2f} vs 争议 0.0833")
print(g)
print(f"\n⚠ **不是工具变量**:`EARLY` 与羞耻直接相关,排除性限制不成立。"
      f"本轮只做**部分定向** —— 少了一个候选解释,不是因果。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
