import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A165 R465 -- 页面上说「广度」的时候,说的是哪一种

`#420a`:**「性行为计数」与「恋物类别计数」可以反号**(−0.076 vs +0.049,同一个预测量)。
**⇒ 那么页面上每一处用「广度」的地方,都必须问:是哪一种?**

⚠ **这一步是**列举**,不是猜**(`#407b`:不要 grep 一个没被标记的属性;
而「哪些量是广度」**是我可以显式列出的** —— 它们都在模型里当过控制项或主量)。

**页面/轮次里用过的广度量(逐条列举,每条注明出处):**
`ncat` **起始类别数**(几乎每个模型的控制项)·
`totalfetishcategory` **恋物类别计数**(`#366` panel · `#419` 锚)·
`Totalsexacts` **认可的性行为计数**(`#393` · `#419` · `#462` 锚)·
`COVB` **块覆盖数**(口径旋钮 ⑩,`#346`/`#371`)·
`PICKS` **所有多选块里勾选的总项数**(`#357b` 的「勾选数」)。

两个活着的世界:
**A 一个东西** -> 两两相关都高(**> 0.6**)-> 「广度」这个词在页面上是安全的;
**B 不是** -> 存在低相关甚至反号的一对 -> **页面上每一句「广度」都必须点名是哪一种**。

ESTIMAND        五个量两两相关(10 对);主量 = **最小的那一对**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**。
                【非零支】**最小的一对 < 0.6** -> 世界 B(必须点名);
                【零支】全都 ≥ 0.6 -> 世界 A。
⚠ 零的种类     `offset_control`:**几个计数之间相关的零绝不是零** ——
                它们都是**计数**,共享「这个人答了多少题」这个成分。
                零 = **在人层打乱其中一个**(`lib.nulls.perm_in`)后重算相关的分布 ——
                **这保住了每个量自己的分布与缺失格局,只打断配对。**
⚠ 多重性       10 对 -> **报分布**,不报单格。
IMPOSSIBLE      ① 五个量的**定义域不同**(有的需要块覆盖 ≥8)-> 用**共同掩码**,同轮报 n;
                ② 「相关低」不等于「其中一个是错的」——只等于**它们不是一个东西**;
                ③ 本轮不问哪一种「更对」。
"""
import numpy as np, pandas as pd, warnings, hashlib, itertools
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
def num(c): return pd.to_numeric(d[c],errors='coerce').values.astype(float)
PICKS=np.zeros(NN); seen=np.zeros(NN,bool)
for Mb,ppl in MB:
    PICKS[ppl]+=Mb.sum(1); seen[ppl]=True
PICKS=np.where(seen,PICKS,np.nan)
BR={'ncat 起始类别数':ncat,'恋物类别计数':num('totalfetishcategory'),
    '性行为计数':num('Totalsexacts'),'COVB 块覆盖数':COVB.astype(float),
    'PICKS 总勾选项数':PICKS}
M=np.ones(NN,bool)
for v in BR.values(): M&=np.isfinite(v)
M&=ok
n=int(M.sum())
print(f"⚠ **列举**(不是猜):**{len(BR)}** 个广度量 · 共同掩码 n=**{n:,}**")
print(f"⚠ **`#392e`:每个量先看它自己**")
for nm,v in BR.items():
    print(f"   {nm:<16} 范围 [{np.nanmin(v[M]):.0f}, {np.nanmax(v[M]):.0f}] · "
          f"众数 **{float(pd.Series(v[M]).mode().iloc[0]):g}** · sd {np.nanstd(v[M]):.2f}")
K=list(BR)
R=np.zeros((5,5))
for i in range(5):
    for j in range(5):
        R[i,j]=np.corrcoef(BR[K[i]][M],BR[K[j]][M])[0,1]
print(f"\n两两相关(10 对):")
print("                  "+"".join(f"{k[:6]:>9}" for k in K))
for i,k in enumerate(K):
    print(f"   {k[:14]:<16}"+"".join(f"{R[i,j]:+9.4f}" for j in range(5)))
pairs=[(K[i],K[j],R[i,j]) for i,j in itertools.combinations(range(5),2)]
pairs.sort(key=lambda t:t[2])
MINR=float(pairs[0][2])
print(f"\n   **最小的一对:{pairs[0][0]} ↔ {pairs[0][1]} = {MINR:+.4f}**")
print(f"   **最大的一对:{pairs[-1][0]} ↔ {pairs[-1][1]} = {pairs[-1][2]:+.4f}**")
print(f"   10 对里 ≥ 0.6 的:**{sum(1 for _,_,r in pairs if r>=0.6)}/10**")
NP_=400; nul=[]
for s_ in range(NP_):
    v2=perm_in(BR[K[0]],M,9200+s_)
    nul.append(float(np.corrcoef(v2[M],BR[K[1]][M])[0,1]))
nul=np.array(nul); THR=float(np.percentile(np.abs(nul),95))
print(f"\n⚠ offset 零(**在人层打乱其中一个**,保住每个量自己的分布与缺失格局,只打断配对;"
      f"**它们都是计数,共享「这个人答了多少题」-> 零不该是零**):")
print(f"   **{nul.mean():+.4f} ± {nul.std():.4f}** · |值| 95 分位 **{THR:.4f}**")
print(f"   -> 所有 10 对里,**{sum(1 for _,_,r in pairs if abs(r)>THR)}** 对越过这个零")
negs=np.array([float(np.corrcoef(perm_in(BR[K[0]],M,99940+s)[M],BR[K[1]][M])[0,1]) for s in range(200)])
rate=float((np.abs(negs)>THR).mean())
print(f"\n负对照(**越阈率**,200 次):**{100*rate:.1f}%**")
T=pd.DataFrame([dict(v_a=a,v_b=b,v_r=r) for a,b,r in pairs]); check_columns(T,'R465')
T.to_csv(pathlib.Path(__file__).parent/'results'/'breadth_pairs.csv',index=False)
NONNULL=MINR<0.6
# ⚠ **第一版把判据本身(0.6)当争议幅度种进去,而那按构造只有 ~50% 检出。**
# 判据是「估计的相关 ≥ 0.6」;**在真值恰好 = 0.6 处种植,估计落在阈两侧各半** ——
# **那不是功率不足,那是阈值的定义。**(与 `#406b` 的饱和同族:**主量的形状决定了扫描能不能量它**。)
# 正确的问法:**真值要低于 0.6 多少,我才判得准「不是一个东西」?**
print(f"\nguard 26 = **MDE 扫描**(⚠ 改问「真值低于 0.6 多少才判得准」),每级 30 次:")
MDE=None; det=[]
for gg in (0.58,0.55,0.50,0.45,0.40):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(8+int(gg*100)*167+s_)
        base=(BR[K[0]][M]-BR[K[0]][M].mean())/BR[K[0]][M].std()
        y=gg*base+np.sqrt(max(1-gg*gg,1e-9))*rg.standard_normal(n)
        if float(np.corrcoef(base,y)[0,1])<0.6: hit+=1     # 正确判为「不是一个东西」
    det.append(hit/30)
    print(f"   真值 **{gg:.2f}** -> 正确判为「不是一个东西」**{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=0.60-gg
MDE_=MDE if MDE else 0.25
print(f"   **MDE(相对 0.6 的偏离)= {MDE_:.2f}** · 而实测最小的一对偏离 0.6 达 **{0.60-MINR:.2f}**")

g=Gate('页面上说「广度」的时候,说的是哪一种')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:真值要低于 0.6 多少才判得准(⚠ 不是在阈值上种)',MDE_,0.60-MINR,True,
    what='偏离扫描',branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(它们都是计数,共享作答量)',nul.std()>0,
           f"{nul.mean():+.4f} ± {nul.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】**最小的一对 < 0.6** -> 世界 B(页面上每句「广度」都必须点名)',NONNULL,
               f"最小 {pairs[0][0]} ↔ {pairs[0][1]} = {MINR:+.4f} · ≥0.6 的 "
               f"{sum(1 for _,_,r in pairs if r>=0.6)}/10")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
