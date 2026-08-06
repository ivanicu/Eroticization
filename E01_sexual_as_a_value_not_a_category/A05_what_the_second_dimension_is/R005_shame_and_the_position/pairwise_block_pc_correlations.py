import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A45 R250 -- 「内容维度」是一个维度,还是 32 个块内成分的和

`#204b` 只证明了**每一块自己**的载荷可复现(|r| = 0.97),
**没有证明 32 块的 PC1 指向同一个人层构念。**
**这一条本该在 `#188` 就做** —— `Cres` 一直是 32 个块 PC 分数的**和**,
而"和"只有在它们彼此相关时才是一个维度。

ESTIMAND        32 块两两之间(共 496 对)PC 人分数在**共同人群**上的相关。
                ⚠ 每块 PC 符号任意 -> 判 **|r|**,并对着**同样取绝对值的置换零**比(`#204b` 的教训)。
KILL            **若真实平均 |r| 不明显高于零 -> "内容维度"是 32 个互不相干成分的和,不是一个维度;
                `#188` 起所有"内容侧"结论都要按这个重新表述。**
NEGATIVE CTRL   题内跨人置换后走同一条管道 -> 零的 |r| 分布(它也取绝对值)。
POSITIVE CTRL   人为把一个共同人层因子种进**全部 32 块** -> 平均 |r| 必须明显上升。
SECOND READ     以第 1 块为**任意锚**对齐符号(锚不是"和",所以不会人为造出正相关),
                判其余 496 对里同号的比例。
IMPOSSIBLE      两块的共同人群因块而异(1,200–13,000),所以每对的 n 不同;
                本轮报的是**未加权**平均,不做跨对排名。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df_raw); RAW=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW.append(dict(qi=q.qi,M=M,ppl=ppl))
print(f"块 {len(RAW)}")

def scores(perm, plant, rng):
    S=[]
    u=rng.standard_normal(NN)
    for b in RAW:
        M=b['M'].copy()
        if perm:
            for j in range(M.shape[1]): M[:,j]=M[rng.permutation(len(M)),j]
        if plant:
            sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
            M=M+plant*np.outer(u[b['ppl']],sub)
        Z=M-M.mean(0,keepdims=True)
        w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
        f=np.full(NN,np.nan); f[b['ppl']]=Z@v[:,-1]
        S.append(f)
    return np.array(S)

def pairwise(S):
    K=len(S); out=[]
    for a in range(K):
        for b in range(a+1,K):
            m=np.isfinite(S[a])&np.isfinite(S[b])
            if m.sum()<300: continue
            out.append((a,b,float(np.corrcoef(S[a][m],S[b][m])[0,1]),int(m.sum())))
    return out

rng=np.random.default_rng(20260803)
P_real=pairwise(scores(False,0.0,rng))
P_null=[pairwise(scores(True,0.0,np.random.default_rng(9000+s))) for s in range(5)]
P_plant=pairwise(scores(False,0.6,np.random.default_rng(11)))
ar=np.mean([abs(x[2]) for x in P_real])
an=np.mean([[abs(x[2]) for x in p] for p in P_null])
sdn=np.std([np.mean([abs(x[2]) for x in p]) for p in P_null])
ap=np.mean([abs(x[2]) for x in P_plant])
print(f"\n对数 {len(P_real)}(共同人群 ≥300)")
print(f"真实 平均 |r| = **{ar:.4f}**  中位 {np.median([abs(x[2]) for x in P_real]):.4f}")
print(f"置换零 平均 |r| = {an:.4f} ± {sdn:.4f}(5 个种子)")
print(f"正对照(全 32 块种同一因子)平均 |r| = {ap:.4f}")

# 以第 1 块为任意锚对齐符号,看其余对的同号率
S=scores(False,0.0,np.random.default_rng(3))
anchor=0; sg=np.ones(len(S))
for k in range(len(S)):
    m=np.isfinite(S[anchor])&np.isfinite(S[k])
    if m.sum()>=300 and np.corrcoef(S[anchor][m],S[k][m])[0,1]<0: sg[k]=-1
Sa=np.array([sg[k]*S[k] for k in range(len(S))])
PA=pairwise(Sa); pos=np.mean([x[2]>0 for x in PA])
print(f"以第 1 块为锚对齐后:{len(PA)} 对里同号(正)的占 **{100*pos:.1f}%**")

T=pd.DataFrame([dict(a=x[0],b=x[1],r=x[2],n=x[3]) for x in P_real])
check_columns(T,'R250'); T.to_csv(pathlib.Path(__file__).parent/'results'/'pairwise.csv',index=False)
g=Gate('内容维度是一个维度还是 32 个成分的和')
g.asserted('正对照:把同一因子种进全部 32 块,平均 |r| 必须明显上升',ap>ar+0.05,
           f"真实 {ar:.4f} -> 种植 {ap:.4f}")
g.negative_control('题内跨人置换的平均 |r|',float(an),float(ar),null_spread=float(sdn))
g.resolvable('真实 |r| 减去零',float(ar-an),float(sdn))
g.asserted('⚠ 零也取绝对值(`#204b` 的教训)',True,f"零的 |r| = {an:.4f},不是 0")
g.asserted('第二读:以任意块为锚对齐后同号率明显高于 50%',pos>0.65,f"{100*pos:.1f}%")
g.asserted('注册的 kill:真实平均 |r| 不明显高于零 -> 不是一个维度',ar-an>2*sdn,
           f"差 {ar-an:+.4f} ± {sdn:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
