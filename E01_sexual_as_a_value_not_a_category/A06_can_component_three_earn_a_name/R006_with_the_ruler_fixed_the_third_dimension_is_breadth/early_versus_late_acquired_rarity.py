import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A64 R283 -- 羞耻贴的是「你现在敞开在哪」,还是「你什么时候开始敞开在那儿」

`#235b`:在与块不相交的起始仪器上,「在更少人去的领域敞开」↔ 羞耻 = **+0.1384**(n=9,944),
比块仪器强,且几乎不含勾选数混杂。而起始仪器带着块仪器**没有的东西:时间**。
所以可以问一个块仪器**结构上问不了**的问题。

WORLDS          ① **当下的地图** —— 羞耻贴的是你**现在**敞开在哪,与何时获得无关 ->
                   早半与晚半的罕见度**同等**预测羞耻
                ② **后来走进去的地方** —— 羞耻贴的是你**后来**才走进去的地方 ->
                   **只有晚半**预测羞耻,那是一个发展性事实
ESTIMAND        对每个人,把他报告了起始年龄的类别按**他自己的**起始年龄中位数劈成早/晚两半,
                各算一次**平均罕见度**(−log 流行度),跨人 z 化后分别与同时对羞耻回归。
KILL            **若晚半的 beta 明显大于早半(差 > 2× 展布)-> 世界②,发展性事实;
                若两者相当 -> 世界①,羞耻贴的是当下的地图。**
⚠ 最强混杂(跑之前写下)
                **晚获得的类别在人群层面天然更罕见**(`#128`–`#212` 整条线)。
                控制:两个量都**跨人 z 化**,人群层面的早/晚罕见度差因此被吸收进各自的均值;
                并**同时报**两半各自的分半信度(若晚半信度更高,beta 天然更大)与解衰减 beta。
                第二个:每半的类别数不同 -> 均值精度不同。要求**每半 ≥4 个类别**并报出类别数。
NEGATIVE CTRL   置换羞耻。
POSITIVE CTRL   两端:只贴晚半的合成结局必须被分开;与两半等相关的必须不被分开。
IMPOSSIBLE      起始年龄是**回忆**,而 `#114` 已测到人把最爱的兴趣记得更早
                (−0.2000 年/评分 sd)。所以「早/晚」带着回忆偏差;
                能判的是**在这份自报时间线上**羞耻贴哪一半,不是真实的发展顺序。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
V=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
OBS=np.isfinite(V); NC=V.shape[1]
RAR=-np.log(np.clip(OBS.mean(0),1e-4,1.))
print(f"起始仪器:{NC} 个类别;罕见度 {RAR.min():.2f}–{RAR.max():.2f}")
print(f"⚠ 人群层面「晚获得更罕见」:corr(类别中位起始年龄, 罕见度) = "
      f"{np.corrcoef(np.nanmedian(np.where(OBS,V,np.nan),0),RAR)[0,1]:+.4f}")

def halves(cats=None):
    """每人按自己的起始年龄中位数劈早/晚,各算平均罕见度。cats 限定用哪些类别(算信度用)。"""
    sel=np.ones(NC,bool) if cats is None else np.isin(np.arange(NC),cats)
    Vm=np.where(OBS&sel[None,:],V,np.nan); Rm=np.where(OBS&sel[None,:],RAR[None,:],np.nan)
    med=np.nanmedian(Vm,1)
    E=Vm<=med[:,None]; L=Vm>med[:,None]
    ne,nl=E.sum(1),L.sum(1)
    re_=np.where(ne>=4,np.nansum(np.where(E,Rm,0),1)/np.maximum(ne,1),np.nan)
    rl_=np.where(nl>=4,np.nansum(np.where(L,Rm,0),1)/np.maximum(nl,1),np.nan)
    return re_,rl_,ne,nl
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
RE,RL,NE,NL=halves()
print(f"每半 ≥4 个类别的人:早 {int(np.isfinite(RE).sum()):,} · 晚 {int(np.isfinite(RL).sum()):,};"
      f"类别数中位 早 {int(np.nanmedian(NE[np.isfinite(RE)]))} · 晚 {int(np.nanmedian(NL[np.isfinite(RL)]))}")
print(f"原始均值(未 z 化):早 {np.nanmean(RE):.4f} · 晚 {np.nanmean(RL):.4f}"
      f"  <- 混杂在这里,z 化后被吸收")
rng=np.random.default_rng(20260804)
rels=[]
for s in range(4):
    p=rng.permutation(NC); a=halves(p[:NC//2]); b=halves(p[NC//2:])
    rr=[]
    for k in (0,1):
        m=np.isfinite(a[k])&np.isfinite(b[k])
        r=float(np.corrcoef(a[k][m],b[k][m])[0,1]); rr.append(2*r/(1+r))
    rels.append(rr)
REL=np.mean(rels,0)
print(f"分半信度(类别劈半):早 **{REL[0]:+.4f}** · 晚 **{REL[1]:+.4f}**")

SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
zE,zL=z(RE),z(RL)
def fit(cols,yy=None):
    yy=y if yy is None else yy
    m=np.isfinite(yy)&np.all(np.isfinite(np.array(cols)),0)
    X=np.column_stack([np.ones(m.sum())]+[c[m] for c in cols])
    zy=(yy[m]-yy[m].mean())/yy[m].std()
    b=np.linalg.lstsq(X,zy,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i],zy[i],rcond=None)[0][j]
        for i in (rng.choice(m.sum(),m.sum(),True) for _ in range(200))])) for j in range(1,len(cols)+1)]
    return b,sd,int(m.sum())
bE,sE,nE=fit([zE]); bL,sL,nL_=fit([zL]); bB,sB,nB=fit([zE,zL])
print(f"\n对羞耻:")
print(f"  早半单独 **{bE[1]:+.4f}** ± {sE[0]:.4f}(n = {nE:,})· 解衰减 {bE[1]/np.sqrt(REL[0]):+.4f}")
print(f"  晚半单独 **{bL[1]:+.4f}** ± {sL[0]:.4f}(n = {nL_:,})· 解衰减 {bL[1]/np.sqrt(REL[1]):+.4f}")
print(f"  **同时放进去:早 {bB[1]:+.4f} ± {sB[0]:.4f} · 晚 {bB[2]:+.4f} ± {sB[1]:.4f}**(n = {nB:,})")
gap=bB[2]-bB[1]; gsd=float(np.hypot(sB[0],sB[1]))
print(f"  晚−早 = **{gap:+.4f}** vs 2×展布 {2*gsd:.4f}")
nul=[fit([zL],rng.permutation(y))[0][1] for _ in range(20)]
print(f"  置换羞耻的零(晚半):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")

m0=np.isfinite(zE)&np.isfinite(zL); n_=rng.standard_normal(NN)
y_late=np.where(m0,0.30*zL+n_,np.nan); y_both=np.where(m0,0.21*(zE+zL)+n_,np.nan)
pl=fit([zE,zL],y_late)[0]; pb=fit([zE,zL],y_both)[0]
print(f"\n正对照两端:只贴晚半 -> 晚−早 **{pl[2]-pl[1]:+.4f}** · 与两半等相关 -> **{pb[2]-pb[1]:+.4f}**")

T=pd.DataFrame([dict(half='早',rel=REL[0],beta_alone=bE[1],sd_alone=sE[0],beta_joint=bB[1],sd_joint=sB[0]),
                dict(half='晚',rel=REL[1],beta_alone=bL[1],sd_alone=sL[0],beta_joint=bB[2],sd_joint=sB[1])])
check_columns(T,'R283'); T.to_csv(pathlib.Path(__file__).parent/'results'/'early_vs_late.csv',index=False)

g=Gate('羞耻贴的是「在哪」还是「什么时候」')
g.asserted('正对照两端:只贴晚半必须分开,两半等相关必须不分开',
           (pl[2]-pl[1])>2*gsd and abs(pb[2]-pb[1])<(pl[2]-pl[1])/2,
           f"只贴晚半 {pl[2]-pl[1]:+.4f} · 等相关 {pb[2]-pb[1]:+.4f} · 2×展布 {2*gsd:.4f}")
g.asserted('⚠ 最强混杂已处理并报出:两半 z 化吸收人群层面的早/晚罕见度差,且两半信度同报',
           True, f"原始均值 早 {np.nanmean(RE):.4f} / 晚 {np.nanmean(RL):.4f};"
                 f"信度 早 {REL[0]:+.4f} / 晚 {REL[1]:+.4f}")
g.negative_control('置换羞耻(晚半)',abs(float(np.mean(nul))),abs(bL[1]),
                   null_spread=float(np.std(nul)),null_kind='跨人置换结局 —— 只打掉配对')
g.offset_control('★ 晚半 vs 早半(同时放进模型)',float(bB[2]),float(bB[1]),gsd,
                 null_kind='同一模型里早半的 beta —— 不是零假设,是「若羞耻贴的是当下的地图,晚半该落在哪」')
g.asserted('★ 注册的 kill:晚半明显大于早半 -> 发展性事实;相当 -> 当下的地图',
           gap>2*gsd, f"早 {bB[1]:+.4f} · 晚 {bB[2]:+.4f};差 {gap:+.4f} vs 2×展布 {2*gsd:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
