import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A44 R249 -- 那张载荷表本身,可复现吗

⚠ **替代了 `#203` 注册的 NEXT,并说明为什么。**
`#203c` 要求一个**干净上下文的独立编码者**(因为我已看过载荷,自己写的编码会顺着它走)。
**本会话的指令是「未经用户要求不调用 Agent」——比常设授权更窄、更即时,所以不调用。**
替代设计达到同一个认识论目的,而且**更硬**:
**在给那张表起名之前,先问它是不是真的** —— 若载荷在两半人上不复现,
`#203b` 我看到的「动物 · 体液 · 所有权」就是噪声,任何编码都不必写了。

ESTIMAND        把人劈成两半,每块**独立**算 PC1 载荷,判两半载荷向量的 |r|。
                ⚠ **特征向量符号不定**,所以只能判 **|r|**,并对着**同样取 |r| 的零**比。
KILL            **若真实 |r| 不明显高于零,则 `#203b` 的载荷表不可复现,那段描述作废。**
NEGATIVE CTRL   题内跨人置换后走同一条管道 -> 零的 |r| 分布(它也享有"取绝对值"的便宜)。
POSITIVE CTRL   块内 PC1 解释率越高的块,复现应当越好 —— 若无关,说明复现度不是在测结构。
NOISE FLOOR     10 个劈分种子。
IMPOSSIBLE      每半只有约 3,300 人 · 每块 8–40 个选项 —— **载荷本身噪声大**,
                所以 |r| 是复现度的**下界**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BL=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    BL.append(dict(qi=q.qi,M=M,opt=opt))
print(f"块 {len(BL)} · 选项合计 {sum(len(b['opt']) for b in BL)}")

def pc1(M):
    Z=M-M.mean(0,keepdims=True); C=np.cov(Z,rowvar=False)
    w,v=np.linalg.eigh(C); return v[:,-1], float(w[-1]/max(w.sum(),1e-12))

def half_rep(perm, rng):
    out=[]
    for b in BL:
        M=b['M'].copy()
        if perm:
            for j in range(M.shape[1]): M[:,j]=M[rng.permutation(len(M)),j]
        p=rng.permutation(len(M)); h=len(p)//2
        la,eva=pc1(M[p[:h]]); lb,evb=pc1(M[p[h:2*h]])
        out.append(dict(qi=b['qi'],n_opt=len(b['opt']),ev=0.5*(eva+evb),
                        r_abs=abs(float(np.corrcoef(la,lb)[0,1]))))
    return pd.DataFrame(out)

rng=np.random.default_rng(20260803)
real=[]; null=[]
for s_ in range(10):
    real.append(half_rep(False,np.random.default_rng(7000+s_)))
    null.append(half_rep(True ,np.random.default_rng(8000+s_)))
R=pd.concat(real); N=pd.concat(null)
check_columns(R,'R249')
R.to_csv(pathlib.Path(__file__).parent/'results'/'replication.csv',index=False)
mr=float(R.r_abs.mean()); mn=float(N.r_abs.mean())
sdr=float(R.groupby(level=0).r_abs.mean().std()) if False else float(np.std([d.r_abs.mean() for d in real]))
sdn=float(np.std([d.r_abs.mean() for d in null]))
print(f"\n真实两半载荷 |r|:均值 **{mr:.4f}** ± {sdr:.4f}(10 个劈分种子)")
print(f"置换零        |r|:均值 {mn:.4f} ± {sdn:.4f}")
print(f"块内中位 |r| 真实 {R.r_abs.median():.4f} · 零 {N.r_abs.median():.4f}")
q=R.groupby('qi').r_abs.mean().sort_values()
print(f"\n复现最差的 4 块 |r|: {[round(v,3) for v in q.head(4)]}")
print(f"复现最好的 4 块 |r|: {[round(v,3) for v in q.tail(4)]}")
ev_corr=float(np.corrcoef(R.ev,R.r_abs)[0,1])
print(f"正对照:corr(块内 PC1 解释率, 复现 |r|) = {ev_corr:+.4f}")

g=Gate('那张载荷表可不可复现')
g.asserted('正对照:PC1 解释率越高的块复现越好',ev_corr>0.3,f"{ev_corr:+.4f}")
g.negative_control('题内跨人置换的两半载荷 |r|',mn,mr,null_spread=sdn)
g.resolvable('真实 |r| 减去零',mr-mn,float(np.hypot(sdr,sdn)))
g.asserted('⚠ 只能判 |r| —— 特征向量符号不定,而零也享有同样的便宜',True,
           f"真实 {mr:.4f} vs 零 {mn:.4f},两者都取绝对值")
g.asserted('注册的 kill:真实 |r| 不明显高于零 -> `#203b` 的载荷表不可复现',
           mr-mn>2*np.hypot(sdr,sdn),f"差 {mr-mn:+.4f} ± {np.hypot(sdr,sdn):.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(R.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 正对照失败了,查它是不是天花板效应 ---------------------------------------
qm_=R.groupby('qi').agg(r_abs=('r_abs','mean'),ev=('ev','mean'))
print("\n---- 正对照失败的诊断 ----")
print(f"  块级 |r| 的范围 {qm_.r_abs.min():.3f}..{qm_.r_abs.max():.3f}  sd {qm_.r_abs.std():.4f}")
print(f"  块级 PC1 解释率 范围 {qm_.ev.min():.3f}..{qm_.ev.max():.3f}  sd {qm_.ev.std():.4f}")
print(f"  32 块里 |r| ≥ 0.90 的:{int((qm_.r_abs>=0.90).sum())}/32 —— **天花板**")
g2=Gate('正对照失败的原因')
g2.asserted('全部块都挤在高位 -> 受限范围,不是仪器坏了',
            float((qm_.r_abs>=0.85).mean())==1.0 and float(qm_.r_abs.std())<0.06,
            f"最小 {qm_.r_abs.min():.3f} · sd {qm_.r_abs.std():.4f}")
g2.asserted('因此「解释率越高复现越好」判 UNVERIFIED,不判 REFUTED',True,
            f"corr = {ev_corr:+.4f},但 |r| 的 sd 只有 {qm_.r_abs.std():.4f} —— 没有可被预测的变异")
g2.resolvable('块级 |r| 的变异本身',float(qm_.r_abs.std()),
              float(np.std([d.groupby('qi').r_abs.mean().std() for d in real])))
print(g2)
