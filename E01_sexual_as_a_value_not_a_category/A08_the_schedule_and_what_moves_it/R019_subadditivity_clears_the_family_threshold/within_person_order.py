import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A137 R414 -- 对一个人来说,冷门的东西是先到的还是后到的,这件事本身跟羞耻有关吗

⛔ **先记一条:`#369` 的 NEXT 跑不了,而它自己的证据早就写在项目里。**
我写的是「人×类别格子 + **人固定效应**」。但 `R382` 的 IMPOSSIBLE 一栏白纸黑字:
**羞耻是人层的一个数,格层回归的因变量在人内是常数。**
**人固定效应会把一个人内恒定的结局整片吸收 —— 那个设计的因变量方差恰好是 0。**
一个跑不了的 NEXT 是一块指向墙的路牌,所以它先被记下来,再被换掉。

**换法**:不要人内**回归**,要人内**摘要** —— 每个人一个数,然后回到人层。

ESTIMAND        ① `ord_i` = 该人**自己**的类别里,「起始年龄」与「类别稀有度」的秩相关。
                   **正 = 冷门的东西对他来说是后到的。**
                ② 主量 = `corr(ord_i, 羞耻)`,**控制类别数与 `S`**(两者都与羞耻有关)。
两个活着的世界    A:顺序本身有内容 —— 冷门先到(「还没有框架的时候它就来了」)与羞耻有关。
                B:顺序只是人口时间表的投影 —— 个体偏差是噪声,与羞耻无关。
KILL(条件式)  仅当对照都过**且 MDE < 0.05** -> 判:`corr(ord_i, 羞耻)` 是否越过零的 95 分位。
⚠ 零的种类     ① 用 `offset_control`:**`ord_i` 的零绝不是零** —— 常见的东西在人口层就是先到的
                   (`#358` 的共同发育时间表),所以**任何人**都会有正的 `ord_i`。
                   零 = **保住人口时间表、打断个人偏差**:把每个人的每个起始年龄
                   换成**该类别自己的观测分布**里的一个抽样(缺失格局不动)。
                ② 用 `negative_control`:**这个零该是零** —— 打乱人,`ord_i` 与羞耻的配对被打断。
POSITIVE CTRL   种一个真的 `ord ↔ 羞耻` 关系 -> 必须被抓到。
⚠ 混淆(跑前写下) `ord_i` 依赖这个人报了多少类别;而勾选数与羞耻有关(`#357b`)。
                -> **类别数与 `S` 同一轮进模型**,不事后加。
⚠ 顺序          **先算 MDE 再看系数**(`#369a` 已成惯例)。
IMPOSSIBLE      起始年龄是回溯自报;`ord_i` 是一个秩相关,对类别数少的人极不稳 -> 设 ≥8 的下限。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R410=(ROOT/'E01_sexual_as_a_value_not_a_category/A136_is_c3_shame_its_own/R410_commonality_vs_person_variables/run.py').read_text()
exec(_R410.split('"""',2)[2].split('base=ok&np.isfinite(C3)')[0])

# ⚠ splice 会把一堆名字塞进我的命名空间,而我**从来看不见那张表**。
# 第一版把人内秩相关命名为 `ORD`,撞上了 R410 带进来的「普通块索引表 ORD」,
# 于是 `fit_apply` 里的 `hs(Ra, ORD)` 拿到一串 float 去做索引。变量遮蔽第四次(`B`·`diff`·`drop`·`ORD`)。
# 便宜的修法:**声明我要新建的名字,跑前和已有命名空间对一次**。
_MINE=['ORDV','RAR','O','J','prev','ncat','onsc','BINo','ordvec','offset_O','coef','MINC']
_CLASH=[x for x in _MINE if x in dir()]
print(f"⚠ 命名空间碰撞检查(splice 之后):{'**撞了:'+str(_CLASH)+'**' if _CLASH else '无碰撞'}")
assert not _CLASH, f"名字撞了:{_CLASH}"

BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
inv=pd.read_csv('data/derived/inventory.csv')
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
O=np.column_stack([d[c].map(BINo).values.astype(float) for c in onsc])   # (NN, J)
J=O.shape[1]
prev=np.isfinite(O).mean(0)                       # 该类别被报的比例 = 流行度
RAR=-np.log(np.clip(prev,1e-4,1.))                # 稀有度
ncat=np.isfinite(O).sum(1).astype(float)
print(f"{J} 个起始类别 · 稀有度范围 [{RAR.min():.3f}, {RAR.max():.3f}] · "
      f"人口层 corr(类别中位起始, 稀有度) = "
      f"**{np.corrcoef(np.nanmedian(np.where(np.isfinite(O),O,np.nan),0),RAR)[0,1]:+.4f}**")
print(f"   ⇒ 人口层「冷门的后到」本来就成立 —— **所以 `ord_i` 的零绝不是零。**\n")

MINC=8
def ordvec(Om):
    """每人一个数:自己的类别里 rank(起始年龄) 与 rank(稀有度) 的相关。"""
    out=np.full(NN,np.nan)
    F=np.isfinite(Om)
    for i in np.flatnonzero(F.sum(1)>=MINC):
        j=np.flatnonzero(F[i]); a=Om[i,j]; b=RAR[j]
        if np.std(a)<1e-9 or np.std(b)<1e-9: continue
        out[i]=np.corrcoef(rankdata(a),rankdata(b))[0,1]
    return out
ORDV=ordvec(O)
mA=ok&np.isfinite(ORDV)&np.isfinite(sh)&np.isfinite(ncat)
ALLR=np.flatnonzero(ok); S=fit_apply(ALLR,ALLR)[0]
mA&=np.isfinite(S); n=int(mA.sum())
print(f"① `ord_i`(≥{MINC} 个类别的人,n={n:,}):均值 **{np.nanmean(ORDV[mA]):+.4f}** · "
      f"中位 {np.nanmedian(ORDV[mA]):+.4f} · sd {np.nanstd(ORDV[mA]):.4f} · "
      f"为正 **{int((ORDV[mA]>0).sum()):,}/{n:,} = {(ORDV[mA]>0).mean():.1%}**")

# ---- offset 零:保住人口时间表,打断个人偏差 ----
def offset_O(seed):
    rg=np.random.default_rng(seed); P=np.full_like(O,np.nan)
    for j in range(J):
        v=O[np.isfinite(O[:,j]),j]; k=np.isfinite(O[:,j])
        P[k,j]=rg.choice(v,size=int(k.sum()),replace=True)   # 该类别自己的分布
    return P
NOFF=40
off=np.array([np.nanmean(ordvec(offset_O(2000+s))[mA]) for s in range(NOFF)])
print(f"⚠ offset 零(**保住人口时间表、打断个人偏差**,{NOFF} 次):"
      f"**{off.mean():+.4f} ± {off.std():.4f}** · 95 分位 {np.percentile(off,95):+.4f}")
print(f"   实测 {np.nanmean(ORDV[mA]):+.4f} -> "
      f"**{(np.nanmean(ORDV[mA])-off.mean())/max(off.std(),1e-12):+.2f} sd** · "
      f"{'**越阈:个人偏差之上还有内容**' if np.nanmean(ORDV[mA])>np.percentile(off,95) else '**未越阈:就是人口时间表**'}")

# ---- ② ord ↔ 羞耻,控制类别数与 S ----
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def coef(y,o=None,g=None):
    g=mA if g is None else g; o=ORDV if o is None else o
    gg=g&np.isfinite(o)&np.isfinite(y); k=int(gg.sum())
    X=np.column_stack([np.ones(k),z(o,gg),z(ncat,gg),z(S,gg)])
    b,*_=np.linalg.lstsq(X,z(y,gg),rcond=None); r=z(y,gg)-X@b
    s2=float(r@r)/(k-4); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1]),k
NUL=np.array([coef(sh[np.random.default_rng(3000+s).permutation(NN)])[0] for s in range(400)])
THR=float(np.percentile(np.abs(NUL),95))
print(f"\n② 零(**打乱人** —— 这个零该是零):**{NUL.mean():+.5f} ± {NUL.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   先算 MDE 再看系数(`#369a`),每级 40 次:")
MDE=None
for g_ in (0.02,0.03,0.04,0.05,0.08):
    hit=0
    for s_ in range(40):
        rg=np.random.default_rng(9000+int(g_*1000)*7+s_)
        y=np.full(NN,np.nan); y[mA]=g_*z(ORDV,mA)+rg.standard_normal(int(mA.sum()))
        if abs(coef(y)[0])>THR: hit+=1
    print(f"      种植 {g_:+.3f} -> 检出 **{hit}/40 = {hit*2.5:>5.1f}%**")
    if MDE is None and hit>=32: MDE=g_
MDE_=MDE if MDE else 0.10
print(f"   **MDE = {MDE_:.3f}** · 有意义 0.05 -> "
      f"{'**可以看系数**' if MDE_<0.05 else '**不要看系数**'}")
b,se,k=coef(sh)
print(f"\n   实测(n={k:,},已控制类别数与 `S`):**{b:+.5f}** · se {se:.5f} · "
      f"95% CI [{b-1.96*se:+.5f}, {b+1.96*se:+.5f}] · "
      f"距零 **{(b-NUL.mean())/max(NUL.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(b)>THR else '**未越阈**'}")
rg=np.random.default_rng(55)
ypos=np.full(NN,np.nan); ypos[mA]=0.12*z(ORDV,mA)+rg.standard_normal(int(mA.sum()))
b_pos=coef(ypos)[0]; b_neg=coef(sh[np.random.default_rng(99).permutation(NN)])[0]
print(f"   正对照(种植 0.12):**{b_pos:+.5f}** vs 阈 {THR:.5f} · "
      f"负对照(打乱人):**{b_neg:+.5f}**")
pd.DataFrame([dict(v_b=b,v_se=se,v_n=k,v_thr=THR,v_mde=MDE_,
                   v_ordmean=float(np.nanmean(ORDV[mA])),v_offmean=float(off.mean()),
                   v_offsd=float(off.std()))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'order.csv',index=False)

g=Gate('冷门的东西先到还是后到,这件事跟羞耻有关吗')
CP=abs(b_pos)>THR; CN=abs(b_neg)<=THR; CO=off.std()>0
g.asserted('★ 正对照:种植 0.12 -> 必须越阈',CP,f"{b_pos:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★ 负对照:打乱人 -> 必须是零',CN,f"{b_neg:+.5f} vs {THR:.5f}",kind='control')
g.asserted('★ offset 零非退化(人口时间表本身就造正的 ord)',CO,
           f"{off.mean():+.4f} ± {off.std():.4f}",kind='control')
if CP and CN and CO and MDE_<0.05:
    g.asserted('★ 注册的 kill:`ord ↔ 羞耻` 越过零的 95 分位',abs(b)>THR,f"{b:+.5f} vs {THR:.5f}")
    g.null_claim_uses_null_criteria('★ guard 21:这个零可发布吗','NULL',
        perm_quantile=float((np.abs(NUL)>abs(b)).mean()),mde=MDE_,
        sensitivity_shown=True,meaningful=0.05)
else:
    g.asserted('★ 注册的 kill(MDE ≥ 0.05 或对照未过 -> 不看系数)',False,
               f"MDE {MDE_:.3f} · 对照 {CP}/{CN}/{CO}")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
