import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A58 R276 -- 「对越轨类别宽」这个名字,挣不挣得到

`#230b`:成分 3 对羞耻给 +0.1286,比位置分的 +0.1185 还大,而两者只相关 +0.2036。
`#230c`:**得分明确,名字未挣**。读作"对越轨类别宽而不是对普通性行为宽"只是载荷的读法(`D5`)。

⚠ **诚实前置**:`c3` 的载荷我**已经看过**(`R274` 打印了高低端各 3 块)。
所以「载荷与块勾选率负相关」这个检验**不是独立的** —— 它只是把我已经看见的图样量化一遍。
`#201` `#202` 两次命名失败的教训是:**名字必须先写成一个能被杀死的预测,再去测。**
所以这一轮的主检验**不是**那个相关,而是一个**从名字推出来、且我没看过**的东西:

WORLDS          ① **名字对**:羞耻贴的是"在越轨类别上敞开" ->
                   把块按**自身平均勾选率**切成「普通」/「越轨」两半,
                   **只有越轨半的宽度分预测羞耻**
                ② **名字错**:羞耻贴的是宽度本身(或别的什么)->
                   两半的宽度分**同样**预测羞耻
ESTIMAND        两半各自的人层残差宽度分对同一道羞耻题的相关(`#179`/`#230` 同一道题)。
KILL            **若越轨半明显高于普通半(差 > 2× 展布,且解衰减后仍成立)-> 名字挣到第一份证据;
                若两半相当 -> 名字被杀掉,和 `#201` `#202` 一样,而 c3 仍是一个没有名字的真实维度。**
⚠ 最强混杂(跑之前写下)
                **越轨半的块可能只是信度更高** -> 相关天然更大。
                控制:**同一轮内报两半各自的分半信度,并给出解衰减后的相关**。
NEGATIVE CTRL   置换羞耻。
POSITIVE CTRL   构造一个只与越轨半宽度分相关的合成结局 -> 这条管道必须把差别测出来;
                再构造一个与两半等相关的合成结局 -> 差别必须消失。**两端都要过。**
IMPOSSIBLE      "越轨"在这里被操作化为**块的平均勾选率低**,那是一个粗代理;
                它不区分"社会禁忌"与"单纯少见"。名字挣到的只是这个代理下的第一份证据。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]; LAB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl)); LAB.append(str(q.col)[:52])
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
RATE=np.array([M.mean() for M,_ in MB])            # 块自身的平均勾选率 —— 外部量
ORD=np.argsort(-RATE)[:NB//2]; TRG=np.argsort(-RATE)[NB//2:]
print(f"块 {NB};n = {int(ok.sum()):,}")
print(f"普通半({len(ORD)} 块)平均勾选率 {RATE[ORD].mean():.3f} · "
      f"越轨半({len(TRG)} 块){RATE[TRG].mean():.3f}")
print("  越轨半最低 4 块:"+' | '.join(f"{LAB[i][:38]} {RATE[i]:.3f}" for i in np.argsort(RATE)[:4]))

def halves(seed):
    rg=np.random.default_rng(seed); H=np.full((2,NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        H[0,b,ppl]=M[:,o[:k]].mean(1); H[1,b,ppl]=M[:,o[k:2*k]].mean(1)
    return H
def profile(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
def halfscore(R,cols):
    sub=R[cols]; F=np.isfinite(sub)
    return np.where(F.sum(0)>=4,np.nansum(np.where(F,sub,0.0),axis=0)/np.maximum(F.sum(0),1),np.nan)
H=halves(500); Ra,Rb=profile(H[0]),profile(H[1])
SC={}; REL={}
for nm,cols in (('普通',ORD),('越轨',TRG)):
    a,b_=halfscore(Ra,cols),halfscore(Rb,cols)
    m=np.isfinite(a)&np.isfinite(b_)&ok
    r=float(np.corrcoef(a[m],b_[m])[0,1]); REL[nm]=2*r/(1+r) if r<0.999 else np.nan
    SC[nm]=(a+b_)/2
print(f"两半宽度分的分半信度:普通 **{REL['普通']:+.4f}** · 越轨 **{REL['越轨']:+.4f}**"
      f"(⚠ 最强混杂就在这两个数上)")

SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
def cr(x,yy=None):
    yy=y if yy is None else yy; m=np.isfinite(x)&np.isfinite(yy)&ok
    return float(np.corrcoef(x[m],yy[m])[0,1]), int(m.sum())
res={nm:cr(SC[nm]) for nm in SC}
dis={nm:res[nm][0]/np.sqrt(max(REL[nm],1e-6)) for nm in SC}
rng=np.random.default_rng(20260804)
boot={}
for nm in SC:
    m=np.flatnonzero(np.isfinite(SC[nm])&np.isfinite(y)&ok)
    boot[nm]=float(np.std([np.corrcoef(SC[nm][s_],y[s_])[0,1]
                           for s_ in (rng.choice(m,len(m),replace=True) for _ in range(200))]))
print(f"\n对羞耻题「{SHN[:46]}」:")
for nm in ('普通','越轨'):
    print(f"  {nm}半宽度分 **{res[nm][0]:+.4f}** ± {boot[nm]:.4f}(n = {res[nm][1]:,})· "
          f"解衰减后 {dis[nm]:+.4f}")
gap=res['越轨'][0]-res['普通'][0]; gsd=float(np.hypot(boot['越轨'],boot['普通']))
dgap=dis['越轨']-dis['普通']
print(f"  **差 {gap:+.4f} vs 2×展布 {2*gsd:.4f}** · 解衰减后的差 {dgap:+.4f}")
nullr=[cr(SC['越轨'],rng.permutation(y))[0] for _ in range(30)]
print(f"  置换羞耻的零:{np.mean(nullr):+.4f} ± {np.std(nullr):.4f}")

# 正对照两端
zt=(SC['越轨']-np.nanmean(SC['越轨']))/np.nanstd(SC['越轨'])
zo=(SC['普通']-np.nanmean(SC['普通']))/np.nanstd(SC['普通'])
n_=rng.standard_normal(NN)
y_only=np.where(np.isfinite(zt),0.35*zt+n_,np.nan)
y_both=np.where(np.isfinite(zt)&np.isfinite(zo),0.25*(zt+zo)+n_,np.nan)
po=cr(SC['越轨'],y_only)[0]-cr(SC['普通'],y_only)[0]
pb=cr(SC['越轨'],y_both)[0]-cr(SC['普通'],y_both)[0]
print(f"\n正对照两端:只贴越轨半的合成结局 -> 差 **{po:+.4f}** · 与两半等相关 -> 差 **{pb:+.4f}**")
# 已看过的那个相关(报告用,标注为非独立)
V3=None
try:
    L3=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A57_is_breadth_a_response_style_or_a_real_dimension'
                        /'R274_is_the_breadth_profile_one_dimension_or_several/results/block_loadings.csv')
    V3=L3.pc3.values
    print(f"⚠(非独立,载荷已看过)corr(c3 块载荷, 块平均勾选率) = "
          f"{np.corrcoef(V3,RATE)[0,1]:+.4f}")
except Exception as e: print(f"  载荷读取失败:{e}")

T=pd.DataFrame(dict(half=['普通','越轨'],n_blocks=[len(ORD),len(TRG)],
                    mean_rate=[RATE[ORD].mean(),RATE[TRG].mean()],
                    rel=[REL['普通'],REL['越轨']],
                    r_shame=[res['普通'][0],res['越轨'][0]],
                    boot_sd=[boot['普通'],boot['越轨']],
                    r_disattenuated=[dis['普通'],dis['越轨']]))
check_columns(T,'R276'); T.to_csv(pathlib.Path(__file__).parent/'results'/'halves_vs_shame.csv',index=False)

g=Gate('「对越轨类别宽」挣不挣得到这个名字')
g.asserted('正对照两端:只贴越轨半 -> 差必须出现;两半等相关 -> 差必须消失',
           po>2*gsd and abs(pb)<po/2, f"只贴越轨 {po:+.4f} · 两半等相关 {pb:+.4f} · 2×展布 {2*gsd:.4f}")
g.asserted('⚠ 最强混杂已量出:两半的信度差,解衰减后的差同时报告',
           True, f"信度 普通 {REL['普通']:+.4f} · 越轨 {REL['越轨']:+.4f};"
                 f"原始差 {gap:+.4f} · 解衰减差 {dgap:+.4f}")
g.negative_control('置换羞耻',abs(float(np.mean(nullr))),abs(res['越轨'][0]),
                   null_spread=float(np.std(nullr)),null_kind='跨人置换结局 —— 保留两个分数,只打掉配对')
g.asserted('★ 注册的 kill:越轨半明显高于普通半(差 > 2× 展布,且解衰减后仍成立)-> 名字挣到第一份证据',
           gap>2*gsd and dgap>0, f"差 {gap:+.4f} vs 2×展布 {2*gsd:.4f};解衰减差 {dgap:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
