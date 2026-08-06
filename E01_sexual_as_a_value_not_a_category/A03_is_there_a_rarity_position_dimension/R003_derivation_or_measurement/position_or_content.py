import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A22 R208 -- 「稀有亲和」这个名字对吗?稳的是"偏爱冷门",还是"喜欢某一类东西"?

#162 的 NEXT:`#100` 的 S(跨不相交块分半 **0.4611**)是"勾选项的平均冷门程度",
它同时含两件事:你**挑了哪些**(内容)与你**挑得多冷门**(位置)。
`#160`/`#161`/`#162` 把「何时」那一侧四个量全部钉在 0.05–0.16;
**「什么」那 0.62 现在是本项目唯一一个大的人层量,而它的名字还没被检验过。**

分离器很直接:**把块按内容相似度分半。**

    POSITION  内容最不像的两半之间,信度仍然在 -> 稳的是**位置**:
              「你偏爱冷门的东西,不管话题是什么」。`#100` 的名字对
    CONTENT   内容不像的两半之间信度塌掉,而内容像的两半之间高 -> 稳的是**内容**:
              「你喜欢某一类东西,而那一类恰好冷门」。**「稀有亲和」这个名字要改**

ESTIMAND        S 的分半信度,在**两种分半**下:随机分半 vs **内容最不相似**分半;
                同一批人、同一个每半块数 k、同一个 Spearman-Brown。
IDENTIFICATION  块间相似度由**块层分数的跨人相关**给出(与 S 的定义无关的另一条信息),
                然后用它把块排成最不相似的两半。
SCOPE           有 >=2k 个可算块的人。
WORLDS          POSITION / CONTENT
KILL            条件式:**两个种植必须分开开火** ——
                种植一个**纯位置**倾向(每人在**所有**块上都挑更冷门的)必须在**不相似**分半上存活;
                种植一个**纯内容**偏好(每人只在**一族**块上挑更冷门的)必须在不相似分半上**塌掉**。
                两个都开火,这个设计才能分辨位置与内容。
POSITIVE CTRL   见上,两个,方向相反。
NEGATIVE CTRL   跨人置换块层分数。
NOISE FLOOR     5 个分半种子。
MULTIPLICITY    k ∈ {4,5,6} x {随机, 不相似} x {真实, 置换, 位置种植, 内容种植},整格发表。
IMPOSSIBLE      "内容相似度"由共现给出,而共现本身与冷门程度相关(#159b 量过 −0.89),
                所以"最不相似"的两半在稀有度构成上也不同 —— 这条写进范围,由种植对照来兜底。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BLK={}; RAR={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    if s.person.nunique()<1200 or s.option.nunique()<8: continue
    br=s.option.map(s.option.value_counts()/s.person.nunique())
    sp=-np.log(np.clip(br,1e-4,1.))
    g=pd.DataFrame({'person':s.person.values,'sp':sp.values}).groupby('person').sp.mean()
    v=np.full(len(df),np.nan); v[g.index.values]=g.values
    BLK[q.qi]=v
    RAR[q.qi]=float(np.mean(-np.log(np.clip(s.option.value_counts()/s.person.nunique(),1e-4,1.))))
B=np.c_[tuple(BLK.values())]; qs=list(BLK)
Bc=B-np.nanmean(B,axis=0,keepdims=True)
ok=np.isfinite(B); nb=ok.sum(1)
print(f"块 {B.shape[1]}  人(>=8 块)= {(nb>=8).sum():,}",flush=True)

# 块间相似度:块层分数的跨人相关(与"冷门程度"的定义无关的另一条信息)
C=pd.DataFrame(Bc).corr(min_periods=500).values.copy(); np.fill_diagonal(C,np.nan)
print(f"块间相关:中位 {np.nanmedian(C):+.3f}  范围 {np.nanmin(C):+.3f}..{np.nanmax(C):+.3f}",flush=True)
w,v=np.linalg.eigh(np.nan_to_num(C,nan=0.)); pc1=v[:,-1]
order=np.argsort(pc1)                                    # 沿最大对比轴排,两端最不相似

def make_perm(seed):
    """每个块独立跨人置换:保留每块的边际分布,只摧毁"同一个人的各块分数属于同一个人"。"""
    rg=np.random.default_rng(seed); P=Bc.copy()
    for j in range(P.shape[1]):
        idx=np.flatnonzero(ok[:,j]); P[idx,j]=Bc[rg.permutation(idx),j]
    return P

def rel(kind,k,seed,mode='random',plant=0.):
    global PERM
    PERM=make_perm(seed*7+11)
    """kind: 'real' | 'perm';mode: 'random' | 'dissimilar';plant: 位置或内容种植强度。"""
    rg=np.random.default_rng(seed); u=rg.standard_normal(len(B))
    fam=set(order[:len(order)//2].tolist())               # "一族"= 对比轴的一端(全局定义,与分半一致)
    A=[];Bb=[]
    for i in np.flatnonzero(nb>=2*k):
        av=np.flatnonzero(ok[i])
        if len(av)<2*k: continue
        if mode=='random':
            p=rg.permutation(av); h1,h2=p[:k],p[k:2*k]
        else:
            # ⚠ 第一版取"这个人自己块里 pc1 最低/最高的 k 个" —— 块少的人两端可能都落在
            #   同一族里,于是内容种植塌不下去(81%)。改成**显式跨族**:h1 全部取自
            #   对比轴的下半族,h2 全部取自上半族。
            lo=[j for j in av if j in fam]; hi=[j for j in av if j not in fam]
            if len(lo)<k or len(hi)<k: continue
            lo=sorted(lo,key=lambda j:pc1[j]); hi=sorted(hi,key=lambda j:-pc1[j])
            h1,h2=np.array(lo[:k]),np.array(hi[:k])
        # ⚠ 第一版写成 `x=Bc[rg.integers(len(Bc))]` —— **把整个人换成另一个人的整条向量**,
        #   两半仍来自同一个人,所以置换什么也没破坏,只是给人换了标签(零给到 0.38 = 效应的 71%)。
        #   正确的零:**每个块独立地跨人置换**,这样一个人的两半来自不同的人。
        x=Bc[i].copy() if kind=='real' else PERM[i].copy()
        if plant:
            # ⚠ 种植必须建在**零背景**上,否则真实信号垫在下面,内容种植塌不下去(第一版 82%)。
            x=PERM[i].copy()
            if plant>0:  x=x+plant*u[i]                                  # 纯位置:所有块同向
            else:        x=x+abs(plant)*u[i]*np.array([1. if j in fam else 0. for j in range(len(x))])
        a,b=np.nanmean(x[h1]),np.nanmean(x[h2])
        if np.isfinite(a) and np.isfinite(b): A.append(a); Bb.append(b)
    A=np.array(A); Bb=np.array(Bb)
    if len(A)<300: return np.nan,len(A)
    r=float(np.corrcoef(A,Bb)[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), len(A)

KS=[4,5,6]; rows=[]
print(f"\n{'k':<4}{'分半':<10}{'n':>7}{'真实 SB':>10}{'置换零':>10}{'位置种植':>10}{'内容种植':>10}")
for k in KS:
    for mode in ['random','dissimilar']:
        sd_=zlib.crc32(f'{mode}{k}'.encode())%9973
        sb,n=rel('real',k,sd_,mode)
        nul=np.nanmean([rel('perm',k,sd_+100+t,mode)[0] for t in range(3)])
        pp,_=rel('real',k,sd_+1,mode,plant=0.8)
        pc,_=rel('real',k,sd_+2,mode,plant=-0.8)
        rows.append(dict(k=k,mode=mode,n=n,sb=sb,null=nul,plant_pos=pp,plant_con=pc))
        print(f"{k:<4}{mode:<10}{n:>7,}{sb:>+10.4f}{nul:>+10.4f}{pp:>+10.4f}{pc:>+10.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'pos_vs_content.csv',index=False)
R=T[T['mode']=='random'].set_index('k'); Dz=T[T['mode']=='dissimilar'].set_index('k')
print(f"\n  随机分半均值 {R.sb.mean():+.4f}   不相似分半均值 {Dz.sb.mean():+.4f}   "
      f"保留 {100*Dz.sb.mean()/R.sb.mean():.0f}%")
print(f"  对照:位置种植 随机 {R.plant_pos.mean():+.4f} / 不相似 {Dz.plant_pos.mean():+.4f}"
      f"(保留 {100*Dz.plant_pos.mean()/R.plant_pos.mean():.0f}%)")
print(f"        内容种植 随机 {R.plant_con.mean():+.4f} / 不相似 {Dz.plant_con.mean():+.4f}"
      f"(保留 {100*Dz.plant_con.mean()/R.plant_con.mean():.0f}%)")

sd=float((R.sb-Dz.sb).std())
g=Gate('「稀有亲和」这个名字对吗')
g.asserted('正对照一:纯位置种植在不相似分半上存活',
           Dz.plant_pos.mean()/max(R.plant_pos.mean(),1e-9)>0.8,
           f"保留 {100*Dz.plant_pos.mean()/R.plant_pos.mean():.0f}%")
g.asserted('正对照二:纯内容种植在不相似分半上塌掉',
           Dz.plant_con.mean()/max(R.plant_con.mean(),1e-9)<0.8,
           f"保留 {100*Dz.plant_con.mean()/R.plant_con.mean():.0f}%")
g.negative_control('跨人置换块层分数',float(abs(Dz.null).mean()),float(Dz.sb.mean()))
g.require_resolvable_first('不相似分半上的信度本身可分辨',float(Dz.sb.mean()),sd)
g.offset_control('不相似分半 vs 随机分半',float(Dz.sb.mean()),float(R.sb.mean()),sd,
                 null_kind='同一批人、同一个 k 下随机分半的 SB 信度(不是零假设,是被比较的基准)')
g.no_sign_crossing('两种分半在每个 k 上同号',list((R.sb-Dz.sb).values))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
