import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A135 R409 -- 「单调」这个词,在嵌套点上和互斥层上不是同一件事

`#364b`:`#360a` 报的「单调」是**嵌套子样本**(阈 ≥k)上的;
本轮的**互斥分层**在低端**不单调**。**而页面上写的是「单调」。**

⚠ **嵌套点之间高度相关**(每个都是下一个的超集)-> **嵌套版本的「单调」比互斥版本更容易出现。**
**这一点必须写进读法,不是写进脚注。**

ESTIMAND        同一批人、同一个 `S`、同一个结局,两种切法**并排**:
                ① **嵌套**(`cov >= k`,k = 4,6,8,10,12,16)· ② **互斥**(六等分层)。
                各自报:序列 · 是否单调 · 两端差 · 以及**真值固定合成结局**在两种切法下的形状。
KILL            **若互斥版本的形状明显不同 -> 页面上「单调」要改成「随纳入阈单调」;
                若一致 -> 措辞不动。**
POSITIVE CTRL   真值固定的合成结局 -> 两种切法下都应当平(而它们各自的残余下降就是各自的基线)。
NEGATIVE CTRL   `perm_finite`。
⚠ guard 22     两条曲线各 6 点。
IMPOSSIBLE      嵌套与互斥是同一份数据的两种读法,不是两个独立证据。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
THS=[4,6,8,10,12,16]
def nested(y):
    out=[]
    for TH in THS:
        S=Spos(cov>=TH); m=(cov>=TH)&np.isfinite(S)&np.isfinite(y)
        out.append(float(np.corrcoef(S[m],y[m])[0,1]) if m.sum()>=200 else np.nan)
    return np.array(out)
S4=Spos(cov>=4); BASE=(cov>=4)&np.isfinite(S4)&np.isfinite(sh)
EDG=np.quantile(cov[BASE],[0,1/6,2/6,3/6,4/6,5/6,1.0])
def disjoint(y):
    out=[]
    for i in range(6):
        lo,hi=EDG[i],EDG[i+1]
        m=BASE&(cov>=lo)&((cov<hi) if i<5 else (cov<=hi))&np.isfinite(y)
        out.append(float(np.corrcoef(S4[m],y[m])[0,1]) if m.sum()>=200 else np.nan)
    return np.array(out)
N_=nested(sh); D_=disjoint(sh)
mono=lambda a: all(a[i]>=a[i+1]-1e-9 for i in range(len(a)-1))
print(f"{'':<8}" + ''.join(f"{k:>10}" for k in THS))
print(f"{'嵌套':<8}" + ''.join(f"{x:>+10.4f}" for x in N_) + f"   单调 **{'是' if mono(N_) else '否'}**")
print(f"{'互斥':<8}" + ''.join(f"{x:>+10.4f}" for x in D_) + f"   单调 **{'是' if mono(D_) else '否'}**")
print(f"\n两端差:嵌套 **{N_[0]-N_[-1]:+.4f}** · 互斥 **{D_[0]-D_[-1]:+.4f}**")
rg=np.random.default_rng(77)
z4=np.full(NN,np.nan); z4[BASE]=(S4[BASE]-S4[BASE].mean())/S4[BASE].std()
ys=np.full(NN,np.nan); ys[BASE]=0.15*z4[BASE]+rg.standard_normal(int(BASE.sum()))
PN,PD=nested(ys),disjoint(ys)
print(f"\n正对照(真值固定 0.15):")
print(f"   嵌套 " + ''.join(f"{x:>+9.4f}" for x in PN) + f"  单调 **{'是' if mono(PN) else '否'}** · "
      f"两端差 **{PN[0]-PN[-1]:+.4f}**")
print(f"   互斥 " + ''.join(f"{x:>+9.4f}" for x in PD) + f"  单调 **{'是' if mono(PD) else '否'}** · "
      f"两端差 **{PD[0]-PD[-1]:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
MN=np.mean([mono(nested(perm_finite(sh,600+i))) for i in range(30)])
MD=np.mean([mono(disjoint(perm_finite(sh,600+i))) for i in range(30)])
print(f"\n★ 负对照(打乱人 30 次)**单调出现的频率**:嵌套 **{100*MN:.0f}%** · 互斥 **{100*MD:.0f}%**")
print(f"   —— 这就是「嵌套的单调更容易出现」的直接测量。")
T=pd.DataFrame([dict(v_k=THS[i],v_nested=float(N_[i]),v_disjoint=float(D_[i])) for i in range(6)])
check_columns(T,'R409'); T.to_csv(pathlib.Path(__file__).parent/'results'/'nd.csv',index=False)
gg=Gate('「单调」在两种切法上不是同一件事')
gg.curve_has_enough_points('⚠ guard 22:两条曲线各几点',THS,min_points=3,what='嵌套/互斥曲线')
gg.asserted('★ 正对照:真值固定的合成结局在两种切法下的形状',True,
            f"嵌套两端差 {PN[0]-PN[-1]:+.4f}(单调 {'是' if mono(PN) else '否'})· "
            f"互斥 {PD[0]-PD[-1]:+.4f}(单调 {'是' if mono(PD) else '否'})")
gg.asserted('★ 负对照:打乱人后**单调**出现的频率(这是「嵌套更容易单调」的直接测量)',
            MN>MD,f"嵌套 {100*MN:.0f}% · 互斥 {100*MD:.0f}%")
gg.asserted('★ 注册的 kill:互斥版本的形状是否明显不同(单调性翻转)',
            mono(N_)!=mono(D_),
            f"嵌套单调 {'是' if mono(N_) else '否'} · 互斥单调 {'是' if mono(D_) else '否'} · "
            f"两端差 {N_[0]-N_[-1]:+.4f} vs {D_[0]-D_[-1]:+.4f}")
gg.asserted('⚠ 嵌套与互斥是同一份数据的两种读法',True,'不是两个独立证据')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
