import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A136 R412 -- 92% 是摊在整条线上,还是只由一端撑起来

`#367a` 的两个数并排着很刺眼:共同线只解释这五个量 **28.9%** 的方差(真单因子 59.6%),
却拿走 `c3⁻` 那份羞耻的 **92%**。**一条薄线做掉了几乎全部的活 —— 这两件事同时为真吗?**

⚠ **先写下一个陷阱,免得我掉进去**:在 PC1 的一个五分位**内部**「控制 PC1」是**代数恒等式** ——
分层已经把 PC1 的方差拿走了,残差回归当然接近零,而那不是发现(`#331a` · `#336b`)。
**所以本轮测的是斜率**:`c3⁻ -> 羞耻` 的关系**是否随 PC1 的位置改变**。

ESTIMAND        按 PC1 分五个五分位,逐段测 `corr(c3⁻, 羞耻)` 与斜率(绝对 Δ);
                主量 = **五段的全距**,对着**它自己的重抽样地板**(同一段内 split-half)。
KILL(条件式)  仅当正/负对照都过 -> 判:**全距是否越过重抽样地板的 95 分位**。
                越过 -> 92% 由某一端撑起,那一端才是这条路的内容;
                没越过 -> **弱共变与强共享可以并存**(共享只需要方向对,不需要幅度大),
                这一条要写进页面,免得读者以为 28.9% 与 92% 里有一个是错的。
POSITIVE CTRL   合成一个「只有最高段有关系」的结局 -> 分段必须只在那一端亮,全距必须越阈。
NEGATIVE CTRL   合成一个「五段真的一样」的结局 -> 全距必须落在地板里。
⚠ 零的种类     `offset_control`:**全距的零不该是零** —— 五段各有采样噪声,全距天然为正。
                地板 = 在**同一份数据**上把每段随机对半再算全距(段内 split-half),取其分布。
⚠ 多重性       五段 -> **报曲线,不报单格**;guard 22 先证明它是一条曲线。
IMPOSSIBLE      分段是描述,不定因果;PC1 的位置与 `c3⁻` 相关,所以段与段之间的人本来就不同。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R411=(ROOT/'E01_sexual_as_a_value_not_a_category/A136_is_c3_shame_its_own/R411_is_it_one_kind_of_person/run.py').read_text()
exec(_R411.split('"""',2)[2].split('def spectrum')[0])
exec('def spectrum'+_R411.split('def spectrum',1)[1].split('R_C3=')[0])
R_C3=[float(np.corrcoef(PD_[n_][m],C3[m])[0,1]) for n_ in FIVE]
SGN=[1. if r>0 else -1. for r in R_C3]
cols=[SGN[i]*PD_[FIVE[i]][m] for i in range(5)]
Z=np.column_stack([(c-c.mean())/c.std() for c in cols])
ev=np.linalg.eigh(np.corrcoef(Z.T))[1][:,-1]
PC=np.full(NN,np.nan); PC[m]=Z@ev
# ⚠ guard 24(`#368a`):特征向量的**符号是任意的**。第一版直接照着它写了标签
# 「PC1 高 = 更年轻/性经验更少/…」,而实测平均羞耻从 Q1 的 1.179 掉到 Q5 的 0.266 —— **正好反过来**。
# 本项目第四次(R210:73 · `#306b` · `#361` · 本轮)。所有从分量读出的**量**都符号不变
# (R² · 共同性 · 保留率 · |cos|),所以门永远看不见它;**它只在散文里现形。**
_gA=Gate('分量必须先被定向,才轮到它的标签')
_anch=_gA.eigenvector_is_anchored('★ PC1 已对着 `c3⁻` 定向',PC,C3,'c3⁻')
if not _anch:
    PC=-PC
    _gA.eigenvector_is_anchored('★ 翻转后重锚',PC,C3,'c3⁻')
print(_gA); print()
mm=m&np.isfinite(sh)&np.isfinite(C3)
n=int(mm.sum()); idx=np.flatnonzero(mm)
print(f"n={n:,} · PC1 = 五个量对齐后的第一主成分\n")

QS=np.quantile(PC[mm],[0,.2,.4,.6,.8,1.]); QS[-1]+=1e-9
def seg(i): return idx[(PC[idx]>=QS[i])&(PC[idx]<QS[i+1])]
def corr_slope(j,y):
    x=C3[j]; yy=y[j]
    g=np.isfinite(x)&np.isfinite(yy)
    if g.sum()<150 or x[g].std()<1e-9: return np.nan,np.nan,int(g.sum())
    r=float(np.corrcoef(x[g],yy[g])[0,1])
    b=float(np.polyfit((x[g]-x[g].mean())/x[g].std(),yy[g],1)[0])   # 绝对 Δ:每 1 sd 的 c3⁻
    return r,b,int(g.sum())

print("逐段(按 PC1 的五分位;**已定向**:PC1 高 = 指向更高的 `c3⁻`):")
rows=[]
for i in range(5):
    j=seg(i); r,b,k=corr_slope(j,sh)
    rows.append(dict(v_q=i+1,v_r=r,v_b=b,v_n=k,v_pcmid=float(np.median(PC[j])),
                     v_shame=float(np.nanmean(sh[j]))))
    print(f"   Q{i+1}  n={k:>5,} · PC1 中位 {np.median(PC[j]):+6.2f} · 平均羞耻 {np.nanmean(sh[j]):.3f} · "
          f"corr(c3⁻,羞耻) **{r:+.4f}** · 斜率 **{b:+.4f}**/sd")
T=pd.DataFrame(rows); check_columns(T,'R412')
T.to_csv(pathlib.Path(__file__).parent/'results'/'quintiles.csv',index=False)
RNG=float(T.v_r.max()-T.v_r.min()); RNGB=float(T.v_b.max()-T.v_b.min())
print(f"\n全距:corr **{RNG:.4f}** · 斜率 **{RNGB:.4f}**/sd")

def floor_draw(seed,y):
    """offset 零 = 「五段真值相同、段大小照旧」时的全距。

    ⚠ 第一版是**段内对半**,而它的正对照**失败了,并且它是对的**:
    当效应集中在一段时,那一段的半样本噪声会把地板自己撑起来(0.4861 vs 地板 0.5567),
    于是**任何集中都会被自己的地板吞掉** —— 一个偏保守的地板对「未越阈」是更强的结论,
    对「越阈」却是一个**永远不会开火的判据**(`#89`:一个不会失败的检验)。
    修法:**把人在五段之间随机重排(段大小不变)** —— 这样五段真值按构造相同,
    而每段仍是**全样本量**,零因此不含半样本的额外噪声。"""
    rg=np.random.default_rng(seed)
    pool=np.concatenate([seg(i) for i in range(5)]); pool=pool[rg.permutation(len(pool))]
    sizes=[len(seg(i)) for i in range(5)]; rs=[]; c=0
    for k in sizes:
        rs.append(corr_slope(pool[c:c+k],y)[0]); c+=k
    rs=[x for x in rs if np.isfinite(x)]
    return [max(rs)-min(rs)] if len(rs)==5 else []
NP_=400
fl=[]; 
for s in range(NP_): fl+=floor_draw(7000+s,sh)
fl=np.array(fl); thr=float(np.percentile(fl,95))
print(f"⚠ offset 地板(**人在五段间随机重排,段大小不变**,{len(fl)} 次;**全距的零不是零** —— 五段各有采样噪声):")
print(f"   **{fl.mean():.4f} ± {fl.std():.4f}** · 95 分位 **{thr:.4f}**")
print(f"   实测全距 **{RNG:.4f}** -> {'**越阈**' if RNG>thr else '**未越阈 —— 五段读作同一个数**'}"
      f" · 距零 {(RNG-fl.mean())/max(fl.std(),1e-12):+.2f} sd")
print(f"   ⚠ 第一版地板用**段内对半**,它的正对照**失败了并且它是对的**:效应集中时那一段的半样本噪声"
      f"把地板自己撑到 0.5567,吞掉了真实的 0.4861 —— **那是一个永远不会开火的判据**。")

# ---- 对照 ----
rg=np.random.default_rng(31)
def synth(load):   # load: 每段的真斜率
    y=np.full(NN,np.nan)
    for i in range(5):
        j=seg(i); x=(C3[j]-np.nanmean(C3[mm]))/np.nanstd(C3[mm])
        y[j]=load[i]*x+rg.standard_normal(len(j))
    return y
ypos=synth([0,0,0,0,0.5]); yneg=synth([0.12]*5)
rp=[corr_slope(seg(i),ypos)[0] for i in range(5)]
rn=[corr_slope(seg(i),yneg)[0] for i in range(5)]
Rp=max(rp)-min(rp); Rn=max(rn)-min(rn)
flp=np.array([x for s in range(150) for x in floor_draw(8000+s,ypos)])
fln=np.array([x for s in range(150) for x in floor_draw(9000+s,yneg)])
print(f"\n正对照(只有 Q5 有关系):逐段 " + ' · '.join(f"{x:+.3f}" for x in rp) +
      f" · 全距 **{Rp:.4f}** vs 地板 95 分位 {np.percentile(flp,95):.4f}")
print(f"负对照(五段真的一样):逐段 " + ' · '.join(f"{x:+.3f}" for x in rn) +
      f" · 全距 **{Rn:.4f}** vs 地板 95 分位 {np.percentile(fln,95):.4f}")

# ---- ⚠ 我要报的是一个**零** -> guard 21 要求 MDE 与**已演示的灵敏度**(`#312a`) ----
print(f"\n⚠ 灵敏度扫描(报零必须先证明能看见什么):在 Q5 种植逐级增大的集中,")
print(f"   每级 20 次,记录全距越过该级自己地板 95 分位的比例:")
MDE=None
for extra in (0.02,0.04,0.06,0.08,0.12,0.20):
    hit=0
    for s_ in range(20):
        rg2=np.random.default_rng(40000+100*int(extra*100)+s_)
        y=np.full(NN,np.nan)
        for i in range(5):
            j=seg(i); x=(C3[j]-np.nanmean(C3[mm]))/np.nanstd(C3[mm])
            y[j]=(0.11+(extra if i==4 else 0.))*x+rg2.standard_normal(len(j))
        rs=[corr_slope(seg(i),y)[0] for i in range(5)]
        if not all(np.isfinite(rs)): continue
        f2=np.array([x for t_ in range(40) for x in floor_draw(50000+1000*s_+t_,y)])
        if (max(rs)-min(rs))>np.percentile(f2,95): hit+=1
    print(f"   Q5 多出 **{extra:+.2f}** 的相关 -> 检出 **{hit}/20 = {hit*5:>3d}%**")
    if MDE is None and hit>=16: MDE=extra
print(f"   **MDE(80% 检出)= Q5 多出 {MDE if MDE else '>0.20'} 的相关**;"
      f"实测五段全距只有 {RNG:.4f}。")

g=Gate('92% 是摊在整条线上还是只由一端撑起来')
CP=Rp>np.percentile(flp,95); CN=Rn<=np.percentile(fln,95); CF=fl.std()>0 and fl.mean()>0
g.asserted('★ 正对照:只有 Q5 有关系 -> 全距必须越阈',CP,
           f"{Rp:.4f} vs {np.percentile(flp,95):.4f}",kind='control')
g.asserted('★ 负对照:五段真的一样 -> 全距必须落在地板里',CN,
           f"{Rn:.4f} vs {np.percentile(fln,95):.4f}",kind='control')
g.asserted('★ offset 地板非退化(全距的零绝不是零)',CF,
           f"地板 {fl.mean():.4f} ± {fl.std():.4f} > 0",kind='control')
g.asserted('★ guard 22:先证明它是一条曲线',
           len(set(T.v_pcmid.round(6)))>=3,f"{len(set(T.v_pcmid.round(6)))} 个不同的 x",kind='control')
g.eigenvector_is_anchored('★ guard 24:PC1 已定向,标签才有意义',PC,C3,'c3⁻')
if CP and CN and CF:
    g.asserted('★ 注册的 kill:五段的全距越过 offset 地板的 95 分位(= 由某一端撑起)',
               RNG>thr,f"实测 {RNG:.4f} vs 阈 {thr:.4f}")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
g.null_claim_uses_null_criteria('★ guard 21:这是一个可发表的零吗','NULL',
    perm_quantile=float((fl>RNG).mean()),mde=float(MDE) if MDE else 0.25,
    sensitivity_shown=True,meaningful=0.05)
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
