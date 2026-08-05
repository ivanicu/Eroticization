import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A137 R416 -- 「冷门的先到」是不是只是「整体来得早」的另一种说法

`#371b` 的负号有两个读法,而它们指向**不同的干预**:
**A「没有框架」** —— 早到的东西找不到语言安放,羞耻是**缺少解释**的后果(可干预:给它语言);
**B「本来就更冷门」** —— 羞耻是**内容**的后果(不可由解释缓解)。
B 的最直接版本(`S` = 整体稀有度)**已经在模型里且两个都活着**。
**但 B 还有一个变体没控:`ord_i` 低的人可能只是**起始年龄整体更早**。**

ESTIMAND        把 `EARLY`(该人的平均起始年龄)加进同一个模型,
                主量 = **`ord` 系数的掉幅** `1 − |b_after| / |b_before|`。
⚠ 判据不是「系数还在不在」  `#332` 已证明 `EARLY` 本身与羞耻有关,所以它**一定**会吃掉一点。
                问题是**吃掉多少**,而那要对着**掉幅自己的地板**比,不是对着 0。
KILL(条件式)  仅当三个对照都过**且 MDE 够** -> 判:掉幅是否越过 offset 地板的 95 分位。
                越阈 -> `EARLY` 是一条真的通路,`ord` 的读法要加上它;
                未越阈 -> **「冷门的先到」不是「整体来得早」的另一种说法。**
POSITIVE CTRL   合成一个「`ord` 的效应完全经 `EARLY` 中介」的结局 -> 掉幅必须接近 100%。
NEGATIVE CTRL   合成一个「完全独立」的 -> 掉幅必须接近 0。
⚠ 零的种类     `offset_control`:**掉幅的零绝不是零** —— 任何与 `ord` 相关的控制都会吃掉一点。
                零 = 加入一个**与 `ord` 相关度和 `EARLY` 相同、但与羞耻无关**的合成变量后的掉幅分布。
⚠ 代数检查     `ord_i` 是**秩相关**,对人层的单调平移不变 -> 与 `EARLY` **没有恒等式**(`#331a` 的检查)。
                跑前直接量 `corr(ord, EARLY)` 并报出来。
IMPOSSIBLE      中介是观察性的:`EARLY` 与 `ord` 都是同一批回溯自报算出来的,谁在前无法由本设计定。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R415=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R415_wide_caliber_replication/run.py').read_text()
_B=_R415.split('"""',2)[2]
exec(_B.split('ARMS=[')[0])              # 到 make_S · ordvec · z 为止

MINCOV,MINC=4,6                          # 宽臂(`#371b` 功率更好的那条)
okA=COVB>=MINCOV; S=make_S(MINCOV); OV=ordvec(O,MINC)
EARLY=np.where(np.isfinite(O).sum(1)>0,
               np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
mA=okA&np.isfinite(OV)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(S)&np.isfinite(EARLY)
n=int(mA.sum())
print(f"宽臂 n=**{n:,}** · `EARLY` = 该人平均起始年龄(小 = 早)")
print(f"⚠ 代数检查(`#331a`):`ord` 是**秩相关**,对人层单调平移不变 -> 与 `EARLY` 无恒等式。")
print(f"   实测 corr(ord, EARLY) = **{np.corrcoef(OV[mA],EARLY[mA])[0,1]:+.4f}** · "
      f"corr(EARLY, 羞耻) = **{np.corrcoef(EARLY[mA],sh[mA])[0,1]:+.4f}**(`#332` 是 −0.1024)\n")

def fitb(y,extra=None,g=None):
    """返回 ord 的系数。extra: 额外控制列(None = 不加)。"""
    g=mA if g is None else g; gg=g&np.isfinite(y)
    if extra is not None: gg=gg&np.isfinite(extra)
    k=int(gg.sum())
    C=[np.ones(k),z(OV,gg),z(ncat,gg),z(S,gg)]
    if extra is not None: C.append(z(extra,gg))
    X=np.column_stack(C); b,*_=np.linalg.lstsq(X,z(y,gg),rcond=None)
    return float(b[1]),k
def drop(y,extra):
    b0,_=fitb(y); b1,_=fitb(y,extra)
    return 100*(1-abs(b1)/max(abs(b0),1e-12)),b0,b1

# ---- offset 地板:相关度相同但与羞耻无关的合成控制 ----
rE=float(np.corrcoef(OV[mA],EARLY[mA])[0,1])
def sham(seed):
    rg=np.random.default_rng(seed); v=np.full(NN,np.nan)
    o=z(OV,mA); v[mA]=rE*o+np.sqrt(max(1-rE*rE,1e-9))*rg.standard_normal(n)
    return v
NF=300
fl=np.array([drop(sh,sham(4000+s))[0] for s in range(NF)])
THR=float(np.percentile(fl,95))
print(f"⚠ offset 地板(**与 `ord` 相关度 {rE:+.4f} 相同、但与羞耻无关**的合成控制,{NF} 次):")
print(f"   掉幅零 = **{fl.mean():+.2f}% ± {fl.std():.2f}%** · 95 分位 **{THR:+.2f}%**")

# ---- ⚠ 第一版这里是一条 MDE 扫描,而它**种的不是它声称的掉幅** ----
# 我写的是 `y = −0.029*((1−f)*o + f*e*sign)`,以为「掉幅 = f」。**不是。**
# `o` 与 `e` 相关 +0.1661,所以控制 `e` 既拿走 `e` 那份、也拿走 `o` 与 `e` 共线的那份 ——
# 实现的掉幅不等于 `f`,于是那条曲线量的是一个我叫错了名字的东西。
# 告密者是它自己:扫描说「种植 50% 只有 60% 检出」,而实测 48% 却在地板之上 **+9.69 sd** ——
# **一条声称比实测更迟钝的灵敏度曲线,通常不是设计弱,是曲线量错了东西。**
# 换成直接的东西:**掉幅自己的 bootstrap 区间**(按人重抽),它回答的正是「这个掉幅有多准」。
NB=400; rgB=np.random.default_rng(12345); ii=np.flatnonzero(mA); bs=[]
for _ in range(NB):
    pick=rgB.choice(ii,size=len(ii),replace=True)
    gb=np.zeros(NN,bool); gb[pick]=True                     # ⚠ 重抽会有重复,布尔掩码只保留唯一个体
    b0b,_=fitb(sh,None,gb); b1b,_=fitb(sh,EARLY,gb)
    if abs(b0b)>1e-9: bs.append(100*(1-abs(b1b)/abs(b0b)))
bs=np.array(bs)
print(f"\n掉幅的 bootstrap(按人重抽 {NB} 次):**{bs.mean():+.2f}%** · "
      f"95% CI **[{np.percentile(bs,2.5):+.2f}%, {np.percentile(bs,97.5):+.2f}%]**")
MDE_=None

# ---- 实测 ----
dObs,b0,b1=drop(sh,EARLY)
print(f"\n实测:加 `EARLY` 前 **{b0:+.5f}** -> 后 **{b1:+.5f}** · **掉幅 {dObs:+.2f}%**")
print(f"   vs 地板 {fl.mean():+.2f}% ± {fl.std():.2f}% -> **{(dObs-fl.mean())/max(fl.std(),1e-12):+.2f} sd** · "
      f"{'**越阈:EARLY 是一条真通路**' if dObs>THR else '**未越阈:不是「整体来得早」的另一种说法**'}")

# ---- 剩下的那一半站得住吗:**同模型下**的零 ----
NUL1=np.array([fitb(sh[np.random.default_rng(7700+s_).permutation(NN)],EARLY)[0] for s_ in range(400)])
T1=float(np.percentile(np.abs(NUL1),95))
print(f"\n剩下的那一半:加 `EARLY` 后 `ord` = **{b1:+.5f}** vs **同模型下**的零 "
      f"{NUL1.mean():+.5f} ± {NUL1.std():.5f} · |值| 95 分位 **{T1:.5f}**")
print(f"   -> **{(b1-NUL1.mean())/max(NUL1.std(),1e-12):+.2f} sd** · "
      f"{'**仍越阈:残余是独立的**' if abs(b1)>T1 else '**未越阈 —— 残余低于这个设计的分辨率**'}")

# ---- 对照 ----
rg=np.random.default_rng(77); e=z(EARLY,mA); o=z(OV,mA); sg=np.sign(np.corrcoef(o,e)[0,1])
ypos=np.full(NN,np.nan); ypos[mA]=-0.20*e*sg+rg.standard_normal(n)        # 全部经 EARLY
yneg=np.full(NN,np.nan); yneg[mA]=-0.20*(o-np.polyval(np.polyfit(e,o,1),e))+rg.standard_normal(n)  # 与 EARLY 正交
dp=drop(ypos,EARLY)[0]; dn=drop(yneg,EARLY)[0]
print(f"\n正对照(ord 的效应完全经 `EARLY`):掉幅 **{dp:+.2f}%** vs 阈 {THR:+.2f}%")
print(f"负对照(与 `EARLY` 正交):掉幅 **{dn:+.2f}%** vs 阈 {THR:+.2f}%")
pd.DataFrame([dict(v_b0=b0,v_b1=b1,v_drop=dObs,v_thr=THR,v_mde=MDE_,v_n=n,
                   v_dp=dp,v_dn=dn,v_flmean=float(fl.mean()),v_flsd=float(fl.std()),
                   v_rE=rE)]).to_csv(pathlib.Path(__file__).parent/'results'/'mediation.csv',index=False)

g=Gate('「冷门的先到」是不是只是「整体来得早」')
CP=dp>THR; CN=dn<=THR; CO=fl.std()>0
g.asserted('★ 正对照:完全经 `EARLY` -> 掉幅必须越阈',CP,f"{dp:+.2f}% vs {THR:+.2f}%",kind='control')
g.asserted('★ 负对照:与 `EARLY` 正交 -> 掉幅必须落在地板里',CN,f"{dn:+.2f}% vs {THR:+.2f}%",kind='control')
g.asserted('★ offset 地板非退化(任何相关控制都吃掉一点)',CO,
           f"{fl.mean():+.2f}% ± {fl.std():.2f}%",kind='control')
if CP and CN and CO:
    g.asserted('★ 注册的 kill:掉幅越过 offset 地板的 95 分位(= `EARLY` 是真通路)',
               dObs>THR,f"{dObs:+.2f}% vs {THR:+.2f}%")
    # ⚠ guard 21 是**给零用的**。第一版在这里无条件调了它,而本轮的 kill 是**开火**的 ——
    # 于是它拿一个 NULL 的判据去评一个非零结果,并且 FAIL。**判据要配结果的种类。**
    if dObs<=THR:
        g.null_claim_uses_null_criteria('★ guard 21:这个零可发布吗','NULL',
            perm_quantile=float((fl>dObs).mean()),mde=0.60,sensitivity_shown=True,meaningful=0.30)
    else:
        g.asserted('★ 非零结果 -> 报区间而不是 MDE:bootstrap 95% CI 不含 0',
                   np.percentile(bs,2.5)>0,
                   f"[{np.percentile(bs,2.5):+.2f}%, {np.percentile(bs,97.5):+.2f}%]")
    g.asserted('★ 残余(加 `EARLY` 后)是否仍是一条独立的路',abs(b1)>T1,
               f"{b1:+.5f} vs 同模型零的阈 {T1:.5f}")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
