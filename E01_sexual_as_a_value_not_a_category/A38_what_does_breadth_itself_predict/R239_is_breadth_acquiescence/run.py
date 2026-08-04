import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A38 R239 -- 勾选数是作答风格,还是真实内容

`#193c` 我把勾选数称作「近乎普适的作答风格因子」—— **那是一个关于测量的断言,而我没检验过它。**
`#184` 已证这 20 道**没有反向计分题**(全是"我…"的正向陈述),
所以"作答风格"在本 release 上**无法用常规办法证伪**。

**但可以用语义天然相反的一对来判**:
    ACQUIESCENCE  K 对语义相反的一对**同号** -> 它测的是"倾向于同意",与内容无关
    CONTENT       K 对语义相反的一对**反号** -> 它测的是真实的偏好广度

ESTIMAND        先用**数据**找出彼此相关最负的题对(而不是我指定),再判 K 对每一对的符号。
KILL            **若 K 对全部强反向对都同号 -> `#193c` 的作答风格读法成立,
                20 道结局的全部相关都要按这个因子重新打折,包括 `#179` 的 +0.1155。**
POSITIVE CTRL   性别镜像对(想象自己是女性/男性地存在)必须在题间相关上**强负** ——
                否则"语义相反"这个筛子失效,整轮不可读。
NEGATIVE CTRL   随机配对的题:K 的符号应当**随机**,不集中。
IMPOSSIBLE      没有反向计分题 -> 本轮只能用**语义**相反,而语义相反 ≠ 心理学相反
                (一个人可以既被支配唤起也被支配唤起的反面唤起)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); cnt=np.zeros(NN); K=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}
    M=np.zeros((len(ppl),len(opt)))
    oi={o:i for i,o in enumerate(opt)}
    M[s.person.map(pi).values,s.option.map(oi).values]=1
    K[ppl]+=M.sum(1); cnt[ppl]+=1
K=np.where(cnt>=8,K,np.nan); bi=np.flatnonzero(np.isfinite(K))
Y=df[lik].values.astype(float)
CM=pd.DataFrame(Y[bi],columns=lik).corr()
print(f"n = {len(bi):,}")

# 数据自己挑出反向对(不是我指定)
pairs=[]
for i in range(len(lik)):
    for j in range(i+1,len(lik)):
        pairs.append((lik[i],lik[j],float(CM.iloc[i,j])))
pairs.sort(key=lambda t:t[2])
strong=[p for p in pairs if p[2]<-0.15]
print(f"题间相关 < −0.15 的强反向对:{len(strong)} 组")

def r_(y,ii):
    m=np.isfinite(y[ii]); jj=ii[m]
    return float(np.corrcoef(y[jj],K[jj])[0,1])

rows=[]
for a,b,rab in strong:
    ra=r_(df[a].values.astype(float),bi); rb_=r_(df[b].values.astype(float),bi)
    rows.append(dict(item_a=a[:46],item_b=b[:46],r_ab=rab,rK_a=ra,rK_b=rb_,
                     same_sign=bool(np.sign(ra)==np.sign(rb_))))
T=pd.DataFrame(rows); check_columns(T,'R239'); T.to_csv(pathlib.Path(__file__).parent/'results'/'pairs.csv',index=False)
print(f"\n{'题间 r':>9}{'r(K,a)':>9}{'r(K,b)':>9}{'同号':>6}  对")
for _,r in T.iterrows():
    print(f"{r.r_ab:>+9.3f}{r.rK_a:>+9.4f}{r.rK_b:>+9.4f}{'YES' if r.same_sign else 'no':>6}  "
          f"{r.item_a[:34]} | {r.item_b[:34]}")
n_same=int(T.same_sign.sum())
print(f"\n强反向对里 K 同号的:{n_same}/{len(T)}")

# 负对照:随机配对
rb=np.random.default_rng(20260803)
rK={c:r_(df[c].values.astype(float),bi) for c in lik}
rand_same=[]
for _ in range(2000):
    a,b=rb.choice(len(lik),2,replace=False)
    rand_same.append(np.sign(rK[lik[a]])==np.sign(rK[lik[b]]))
print(f"负对照(随机配对 2000 次)同号率 {100*np.mean(rand_same):.0f}%")
# 正对照:性别镜像对必须强负
mir=[p for p in pairs if 'existing' in p[0] and 'existing' in p[1]]
print(f"正对照(性别镜像对「想象自己是女性/男性地存在」)题间 r = {mir[0][2]:+.3f}" if mir else "正对照:未找到镜像对")

g=Gate('勾选数是作答风格还是真实内容')
g.asserted('正对照:性别镜像对在题间相关上强负',bool(mir and mir[0][2]<-0.15),
           f"{mir[0][2]:+.3f}" if mir else '未找到')
g.asserted('可判前提:确实存在强反向对',len(T)>=3,f"{len(T)} 组")
g.asserted('负对照:随机配对的同号率不是 100%',np.mean(rand_same)<0.95,
           f"{100*np.mean(rand_same):.0f}%")
g.asserted('注册的 kill:K 对全部强反向对都同号 -> 作答风格读法成立',n_same==len(T),
           f"同号 {n_same}/{len(T)}")
print(g)
print(f"\n  => {'ACQUIESCENCE —— 20 道结局的全部相关要按这个因子打折' if n_same==len(T) else 'CONTENT —— K 不是纯作答风格,`#193c` 的措辞要收窄'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
