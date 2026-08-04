import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A22 R210 -- 「内容」这一侧到底有多大?

#163b:位置倾向里有一个小的**内容**成分(15%,5.7×)。
#159:「何时」那点结构几乎全在**类别**上(题目侧 0.96 vs 人侧 0.105)。
#164:位置倾向不是性别的影子(去性别后保留 102%)。

那么还没被问过的第三件事:**"这个人勾了哪些具体选项"作为一个人层性质,有多稳?**

    CONTENT>POSITION  内容分半信度明显高于位置的 0.60 -> 本项目一直在测的"位置"
                      只是更大的一个"内容"的投影,而 `#100` 的框架要重排
    COMPARABLE        两者同量级 -> 位置与内容是两个并列的人层维度
    POSITION>CONTENT  内容更低 -> "偏爱冷门"确实是这套数据里最稳的人层量

ESTIMAND        内容分数的分半信度:把块劈成不相交两半,在**每半的选项矩阵**上取主成分给人打分,
                跨人相关 + Spearman-Brown;并在**把位置分数残差化掉之后**再量一次。
IDENTIFICATION  两半在**块**上不相交,所以没有共享 item。位置(冷门程度)在同一半上算,
                作为协变量回归掉 —— 剩下的就是"与冷门程度无关的、你挑了哪些"。
SCOPE           两半各 >=k 个块的人。
WORLDS          见上三个。
KILL            条件式:**两个种植分开开火** —— 种一个纯**内容**信号(每人在一个固定的
                选项子集上多勾)必须被内容分数测到;种一个纯**位置**信号(每人处处挑更冷门)
                在**残差化后**的内容分数上必须**测不到**。
POSITIVE CTRL   见上,两个。
NEGATIVE CTRL   每块独立跨人置换(`#163c` 修好的那个零)。
NOISE FLOOR     3 个分半种子。
MULTIPLICITY    k ∈ {6,8} x {内容, 去位置后的内容, 位置} x {真实, 置换, 两种种植},整格发表。
IMPOSSIBLE      "内容"由主成分定义,所以它测的是**共享的**内容维度,不是"这个人独有的组合"。
                一个完全特异的组合(只有他一个人这么勾)对任何跨人成分都是不可见的。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BL={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    BL[q.qi]=dict(M=M,ppl=ppl,rar=-np.log(np.clip(M.mean(0),1e-4,1.)))
qs=list(BL); N=len(df)
print(f"块 {len(qs)}  选项合计 {sum(BL[q]['M'].shape[1] for q in qs)}",flush=True)

# 人 x 选项:按块给人打两个分 —— 内容(主成分)与位置(冷门程度)
# ⚠ #165:第一版的种植向量 u 在两半里各抽一次 -> 种下去的不是人层信号而是两半各自的噪声,
#        只会**稀释**真信号(正对照反而下降)。u 必须是**全局人层的同一个向量**。
#        同理四个臂必须共用**同一个块劈分**,否则 n 从 682 跳到 1346,臂之间不可比。
def build(half,u=None,plant=0.,kind='content',permute=False,pseed=0):
    con=np.zeros(N); pos=np.zeros(N); cnt=np.zeros(N)
    for bi,q in enumerate(half):
        M=BL[q]['M'].copy(); ppl=BL[q]['ppl']; rar=BL[q]['rar']
        if permute:                      # #163c:每块独立跨人置换
            rg2=np.random.default_rng(zlib.crc32(f'{q}{pseed}{bi}'.encode())%(1<<30))
            M=M[rg2.permutation(len(M))]
        if plant:
            up=u[ppl]                    # 同一个人在两半里拿到同一个 u
            if kind=='content':
                sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
                M=M+plant*np.outer(up,sub)
            else:                        # 纯位置:处处按冷门程度多勾
                M=M+plant*np.outer(up,(rar-rar.mean())/max(rar.std(),1e-9))
        Z=M-M.mean(0,keepdims=True)
        w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); pc=v[:,-1]
        con[ppl]+=Z@pc
        pos[ppl]+=(M@rar)/np.maximum(M.sum(1),1)
        cnt[ppl]+=1
    ok=cnt>=len(half)//2
    return (np.where(ok,con/np.maximum(cnt,1),np.nan),
            np.where(ok,pos/np.maximum(cnt,1),np.nan))

def sb(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<300: return np.nan,int(m.sum())
    r=float(np.corrcoef(a[m],b[m])[0,1])
    return (2*abs(r)/(1+abs(r)) if abs(r)<0.999 else np.nan),int(m.sum())

def resid(x,z):
    m=np.isfinite(x)&np.isfinite(z); out=np.full_like(x,np.nan)
    X=np.c_[np.ones(m.sum()),z[m]]
    out[m]=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]
    return out

ARMS=[('real',0.,'content',False),('perm',0.,'content',True),
      ('plant_content',0.6,'content',False),('plant_position',0.6,'position',False)]
KS=[6,8]; SEEDS=[0,1,2]; rows=[]
print(f"\n{'k':<4}{'seed':<6}{'臂':<15}{'n':>7}{'内容 SB':>10}{'去位置后':>10}{'位置 SB':>10}")
for k in KS:
    for sd_ in SEEDS:
        rg=np.random.default_rng(zlib.crc32(f'split{k}{sd_}'.encode())%(1<<30))
        p=rg.permutation(len(qs)); hA=[qs[i] for i in p[:k]]; hB=[qs[i] for i in p[k:2*k]]
        u=rg.standard_normal(N)                      # 全局人层的种植向量,两半共用
        for tag,pl,kind,perm in ARMS:
            cA,pA=build(hA,u,pl,kind,perm,sd_); cB,pB=build(hB,u,pl,kind,perm,sd_)
            s_c,n=sb(cA,cB); s_p,_=sb(pA,pB); s_cr,_=sb(resid(cA,pA),resid(cB,pB))
            rows.append(dict(k=k,seed=sd_,arm=tag,n=n,content=s_c,content_resid=s_cr,position=s_p))
            print(f"{k:<4}{sd_:<6}{tag:<15}{n:>7,}{s_c:>+10.4f}{s_cr:>+10.4f}{s_p:>+10.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'content_vs_position.csv',index=False)
R,P=T[T.arm=='real'],T[T.arm=='perm']
PC,PP=T[T.arm=='plant_content'],T[T.arm=='plant_position']
sd=float(R.content_resid.std())
print(f"\n  真实   :内容 {R.content.mean():+.4f}   去位置后 {R.content_resid.mean():+.4f} (sd {sd:.4f})   位置 {R.position.mean():+.4f}")
print(f"  置换零 :内容 {P.content.mean():+.4f}   去位置后 {P.content_resid.mean():+.4f}   位置 {P.position.mean():+.4f}")
print(f"  内容种植:内容 {PC.content.mean():+.4f} 去位置后 {PC.content_resid.mean():+.4f}")
print(f"  位置种植:内容 {PP.content.mean():+.4f} **去位置后 {PP.content_resid.mean():+.4f}**  位置 {PP.position.mean():+.4f}")

g=Gate('内容这一侧有多大')
g.asserted('正对照一:纯内容种植抬高内容分数',PC.content.mean()>R.content.mean()+0.05,
           f"{R.content.mean():+.4f} -> {PC.content.mean():+.4f}")
g.asserted('正对照二:纯位置种植在去位置后的内容分数上测不到',
           PP.content_resid.mean()<R.content_resid.mean()+2*sd,
           f"位置种植去位置后 {PP.content_resid.mean():+.4f} vs 真实 {R.content_resid.mean():+.4f} ± {sd:.4f}")
g.asserted('正对照二之二:位置种植确实抬高了**位置**分数',PP.position.mean()>R.position.mean(),
           f"位置 {R.position.mean():+.4f} -> {PP.position.mean():+.4f}")
g.negative_control('每块独立跨人置换(内容)',float(abs(P.content.mean())),float(R.content.mean()),
                   null_spread=float(P.content.std()))
g.require_resolvable_first('去位置后的内容信度可分辨',float(R.content_resid.mean()),sd)
g.offset_control('去位置后的内容 vs 位置',float(R.content_resid.mean()),float(R.position.mean()),sd,
                 null_kind='同一批人、同一块劈分上的位置分数分半信度 —— 不是零假设,是被比较的对象')
g.threshold_outside_noise('预注册阈值 0.60',float(R.content_resid.mean()),0.60,sd)
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 泄漏底噪校准(#165 的核心) ----------------------------------------
# 「去位置后的内容」这个零**不该是 0**:位置分数自带测量误差,回归掉一个有噪声的协变量
# 会留下衰减残留,而这个残留会被主成分读成"内容"。所以要造一个**只有位置、没有内容**的
# 合成世界(置换背景 + 纯位置种植),把种植强度调到它的**位置信度等于真实的 0.60**,
# 再读它的"去位置后内容" —— 那才是真实值 +0.26 该被比较的零。
print("\n---- 泄漏底噪:置换背景 + 纯位置种植,按位置信度匹配 ----")
cal=[]
for k in KS:
    for sd_ in SEEDS:
        rg=np.random.default_rng(zlib.crc32(f'split{k}{sd_}'.encode())%(1<<30))
        p=rg.permutation(len(qs)); hA=[qs[i] for i in p[:k]]; hB=[qs[i] for i in p[k:2*k]]
        u=rg.standard_normal(N)
        for st in [0.0,0.05,0.10,0.15,0.25,0.40,0.60]:
            cA,pA=build(hA,u,st,'position',True,sd_); cB,pB=build(hB,u,st,'position',True,sd_)
            s_p,n=sb(pA,pB); s_cr,_=sb(resid(cA,pA),resid(cB,pB))
            cal.append(dict(k=k,seed=sd_,strength=st,n=n,position=s_p,content_resid=s_cr))
C=pd.DataFrame(cal); C.to_csv(pathlib.Path(__file__).parent/'results'/'leak_calibration.csv',index=False)
tab=C.groupby('strength')[['position','content_resid']].agg(['mean','std'])
print(tab.round(4).to_string())

obs=float(R.position.mean())
gm=C.groupby('strength')[['position','content_resid']].mean()
i=(gm.position-obs).abs().idxmin()
leak=float(gm.loc[i,'content_resid']); leak_sd=float(C[C.strength==i].content_resid.std())
print(f"\n  真实位置信度 {obs:+.4f} -> 匹配到种植强度 {i} (位置 {gm.loc[i,'position']:+.4f})")
print(f"  该点的**泄漏底噪** = {leak:+.4f} ± {leak_sd:.4f}")
print(f"  真实的去位置后内容  = {R.content_resid.mean():+.4f} ± {sd:.4f}")

g2=Gate('去位置后的内容,是不是只是位置的泄漏')
g2.asserted('校准是单调的:种植越强,位置信度越高',
            bool((gm.position.diff().dropna()>0).all()),
            ' -> '.join(f'{v:.3f}' for v in gm.position))
g2.asserted('强度 0 时位置信度确实塌到零(置换背景生效)',abs(float(gm.loc[0.0,'position']))<0.15,
            f"{gm.loc[0.0,'position']:+.4f}")
g2.offset_control('真实内容 vs 位置泄漏底噪',float(R.content_resid.mean()),leak,
                  float(np.hypot(sd,leak_sd)),
                  null_kind='位置信度匹配到 0.60 的纯位置合成世界 —— 构造上内容为零,读到的全是回归衰减残留')
g2.resolvable('真实内容减去泄漏底噪',float(R.content_resid.mean()-leak),float(np.hypot(sd,leak_sd)))
print(g2)
print(f"\ncal sha1 {hashlib.sha1(C.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 两个人层分数彼此相关吗 ----------------------------------------------
print("\n---- 内容分 vs 位置分,在人层上的相关 ----")
cc=[]
for k in KS:
    for sd_ in SEEDS:
        rg=np.random.default_rng(zlib.crc32(f'split{k}{sd_}'.encode())%(1<<30))
        p=rg.permutation(len(qs)); hA=[qs[i] for i in p[:k]]; hB=[qs[i] for i in p[k:2*k]]
        for nm,h in (('A',hA),('B',hB)):
            c,pz=build(h)
            m=np.isfinite(c)&np.isfinite(pz)
            cc.append(dict(k=k,seed=sd_,half=nm,n=int(m.sum()),
                           r=float(np.corrcoef(c[m],pz[m])[0,1])))
X=pd.DataFrame(cc); X.to_csv(pathlib.Path(__file__).parent/'results'/'content_position_corr.csv',index=False)
print(X.groupby('k').r.agg(['mean','std','min','max']).round(4).to_string())
rbar=float(X.r.abs().mean()); rsd=float(X.r.std())
print(f"\n  |r| 平均 {rbar:.4f} ± {rsd:.4f}   共享方差 {100*rbar**2:.1f}%")
g3=Gate('内容与位置是同一个维度吗')
g3.resolvable('两分数的相关本身可分辨',rbar,rsd)
g3.equivalent_within('内容与位置近似正交',rbar,rsd,margin=0.30)
print(g3)
