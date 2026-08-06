import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A133 R405 -- 把「谁被纳入」和「`S` 被测得多好」分开

`#360` 的边界写着这两者分不开。**但有一个设计可以:
在**同一批人**(固定为 `cov>=16` 的那 950 人)上,用**不同数量的块**重估 `S`。**
样本不变 -> 只有估计质量在变。

⚠ **先报 MDE**:n=950 -> 2.8/√950 ≈ **0.091**。
**若 MDE 大过要测的差,如实报「这个子样本答不了」并停。**

ESTIMAND        固定 `cov>=16` 的人;k ∈ {4,6,8,12,16}:随机抽 k 个**他们都答过**的块重估 `S`,
                每档 ≥20 次抽样,报 `corr(S_k, 羞耻)` 的均值与展布。
KILL            **若曲线随 k 明显变化 -> 那 0.10 里有「估计质量」的成分;
                若曲线平 -> 那 0.10 是「谁被纳入」,是一条关于**人**的结论。**
POSITIVE CTRL   在**同一批人**上合成一个真值固定的关系 -> 曲线必须平。
NEGATIVE CTRL   `perm_finite`。
IMPOSSIBLE      这 950 人是**最极端**的子样本(答了 ≥16 块);结论只对他们成立,
                不能外推到全样本 —— **写在设计里。**
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN)
IN=np.zeros((NB,NN),bool)
for b,(M,ppl) in enumerate(MB): cov[ppl]+=1; IN[b,ppl]=True
FIX=np.flatnonzero((cov>=16)&np.isfinite(sh))
n=len(FIX); mde=2.8/np.sqrt(n)
print(f"固定样本:`cov>=16` 且羞耻有值 -> **n = {n:,}** · **MDE ≈ {mde:.4f}**")
print(f"⚠ 这 {n:,} 人是最极端的子样本(答了 ≥16 块);结论只对他们成立,不外推。\n")
COMMON=[b for b in range(NB) if IN[b,FIX].mean()>0.95]
print(f"他们几乎都答过的块:**{len(COMMON)}/{NB}**")
def S_from(blocks,rows):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for b in blocks:
        M,ppl=MB[b]
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    out=np.full(NN,np.nan); k=cv>=max(1,len(blocks)//2)
    out[k]=ps[k]/cv[k]; return out
def curve(y,seed=0,T=24):
    rg=np.random.default_rng(seed); res={}
    for k in (4,6,8,12,min(16,len(COMMON))):
        if k>len(COMMON): continue
        vals=[]
        for _ in range(T):
            bl=rg.choice(COMMON,k,replace=False)
            Sk=S_from(bl,FIX)
            m=np.isfinite(Sk)&np.isfinite(y); m2=np.zeros(NN,bool); m2[FIX]=True; m&=m2
            if m.sum()>=200: vals.append(float(np.corrcoef(Sk[m],y[m])[0,1]))
        if vals: res[k]=(float(np.mean(vals)),float(np.std(vals)),len(vals))
    return res
R=curve(sh,seed=11)
print(f"\n{'k 块':>6}{'corr(S_k, 羞耻)':>18}{'跨抽样展布':>12}")
for k,(m_,s_,c_) in R.items(): print(f"{k:>6}{m_:>+18.4f}{s_:>12.4f}")
ks=sorted(R); span=max(R[k][0] for k in ks)-min(R[k][0] for k in ks)
print(f"   曲线全距 **{span:.4f}** vs **MDE {mde:.4f}**")
rg=np.random.default_rng(5)
Sfull=S_from(COMMON,FIX); mf=np.isfinite(Sfull)
zz=np.full(NN,np.nan); zz[mf]=(Sfull[mf]-np.nanmean(Sfull[mf]))/np.nanstd(Sfull[mf])
ysyn=np.full(NN,np.nan); ysyn[mf]=0.20*zz[mf]+rg.standard_normal(int(mf.sum()))
P=curve(ysyn,seed=7)
psp=max(P[k][0] for k in P)-min(P[k][0] for k in P)
print(f"\n正对照(同一批人,真值固定 0.20):" + ' · '.join(f"k={k} **{P[k][0]:+.4f}**" for k in sorted(P)))
print(f"   全距 **{psp:.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
NGv=curve(perm_finite(sh,909),seed=3)
print(f"负对照(打乱人):" + ' · '.join(f"k={k} **{NGv[k][0]:+.4f}**" for k in sorted(NGv)))
T=pd.DataFrame([dict(v_k=k,v_r=R[k][0],v_sd=R[k][1]) for k in ks])
check_columns(T,'R405'); T.to_csv(pathlib.Path(__file__).parent/'results'/'fixed.csv',index=False)
gg=Gate('把「谁被纳入」和「测得多好」分开')
gg.asserted('★ 正对照:同一批人上真值固定的关系 -> 曲线必须平(全距 < MDE)',psp<mde,
            f"全距 {psp:.4f} vs MDE {mde:.4f}")
gg.asserted('★ 负对照:打乱人后全部 ≈ 0',max(abs(NGv[k][0]) for k in NGv)<0.08,
            ' · '.join(f"k={k} {NGv[k][0]:+.4f}" for k in sorted(NGv)))
gg.asserted('★ 注册的 kill:曲线随 k 是否明显变化(全距 > MDE)',span>mde,
            f"全距 **{span:.4f}** vs MDE **{mde:.4f}** —— "
            f"{'有估计质量的成分' if span>mde else '**平:那 0.10 是「谁被纳入」**'}")
gg.null_claim_uses_null_criteria('★ guard 21:若判平(零),三件套在不在',
    'NULL' if span<=mde else 'EFFECT',perm_quantile=0.5,mde=mde,
    sensitivity_shown=f"正对照全距 {psp:.4f};负对照全部 ≈ 0",meaningful=0.10)
gg.asserted('⚠ 边界:这 950 人是最极端的子样本',True,'结论只对他们成立,不外推到全样本')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
