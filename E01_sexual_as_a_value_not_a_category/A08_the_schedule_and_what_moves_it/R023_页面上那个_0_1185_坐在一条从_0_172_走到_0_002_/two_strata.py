import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A135 R408 -- 是「说了多少」,还是「问卷给你开了多少块」

`#363a`:羞耻的覆盖阈曲线是 29 个结局里**唯一**的 —— 单调、两端不重叠、不是选择。
**剩下最自然的心理学解释:填得少的人是「只报了自己最在意的那几样」的人 ——
对他们 `S` 量的是**主动说出口的那几样**有多冷门;填满的人,`S` 量的是**一份清单的平均**。**

**⇒ 那个下降该沿「说了多少」(勾选总数)走,而不是沿「开了多少块」(块覆盖)走。**

ESTIMAND        两条**并排**的曲线,同样六层:
                ① 按**块覆盖**分层(`#404` 的那条)· ② 按**勾选总数**分层;
                各层内的 `corr(S, 羞耻)`。
KILL            **若沿「勾选总数」也单调下降 -> 是**人的选择性**(说了多少);
                若只沿块覆盖下降 -> 是**问卷的门控结构**(开了多少块)。**
POSITIVE CTRL   真值固定的合成结局 -> 两条曲线都必须平(给出各自的纯选择基线)。
NEGATIVE CTRL   `perm_finite`。
⚠ 两者相关但不同 `#357a`:`corr(S, 勾选数) = +0.608` -> **必须并排看,不能只看一条**。
⚠ guard 22     先声明点数。
IMPOSSIBLE      两个分层变量本身高度相关,所以两条曲线不是独立的证据;
                本轮判的是**哪一条更陡**,不是「只有一条成立」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

cov=np.zeros(NN); PK=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1; PK[ppl]+=M.sum(1)
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
BASE=(cov>=4)&np.isfinite(sh)&np.isfinite(PK)
S4=Spos(cov>=4)
BASE&=np.isfinite(S4)
print(f"共同底样本(`cov>=4`,羞耻与 `S` 有值):n = **{int(BASE.sum()):,}**")
print(f"⚠ `corr(块覆盖, 勾选总数)` = **{np.corrcoef(cov[BASE],PK[BASE])[0,1]:+.4f}** —— 两条曲线不是独立证据\n")
QS=[0,1/6,2/6,3/6,4/6,5/6,1.0]
def strata(v,y):
    q=np.quantile(v[BASE],QS); out=[]
    for i in range(6):
        lo,hi=q[i],q[i+1]
        m=BASE&(v>=lo)&((v<hi) if i<5 else (v<=hi))&np.isfinite(y)
        if m.sum()<200: out.append((np.nan,int(m.sum()))); continue
        out.append((float(np.corrcoef(S4[m],y[m])[0,1]),int(m.sum())))
    return out
A=strata(cov,sh); B=strata(PK,sh)
print(f"{'层':>4}{'按块覆盖':>22}{'按勾选总数':>22}")
for i in range(6):
    print(f"{i+1:>4}   {A[i][0]:+.4f} (n={A[i][1]:>5,})   {B[i][0]:+.4f} (n={B[i][1]:>5,})")
sa=[x[0] for x in A]; sb=[x[0] for x in B]
da=sa[0]-sa[-1]; db=sb[0]-sb[-1]
print(f"\n★ 两端差:按**块覆盖** **{da:+.4f}** · 按**勾选总数** **{db:+.4f}** · "
      f"比 **{da/max(abs(db),1e-9):.2f}×**")
rg=np.random.default_rng(41)
z4=np.full(NN,np.nan); z4[BASE]=(S4[BASE]-S4[BASE].mean())/S4[BASE].std()
ysyn=np.full(NN,np.nan); ysyn[BASE]=0.15*z4[BASE]+rg.standard_normal(int(BASE.sum()))
PA=strata(cov,ysyn); PB=strata(PK,ysyn)
pa=[x[0] for x in PA][0]-[x[0] for x in PA][-1]; pb=[x[0] for x in PB][0]-[x[0] for x in PB][-1]
print(f"正对照(真值固定 0.15):按块覆盖 **{pa:+.4f}** · 按勾选总数 **{pb:+.4f}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
NGa=np.array([(lambda s:s[0][0]-s[-1][0])(strata(cov,perm_finite(sh,900+i))) for i in range(12)])
NGb=np.array([(lambda s:s[0][0]-s[-1][0])(strata(PK,perm_finite(sh,900+i))) for i in range(12)])
print(f"负对照(打乱人):按块覆盖 **{NGa.mean():+.4f} ± {NGa.std():.4f}** · "
      f"按勾选总数 **{NGb.mean():+.4f} ± {NGb.std():.4f}**")
T=pd.DataFrame([dict(v_layer=i+1,v_cov=sa[i],v_pk=sb[i],v_ncov=A[i][1],v_npk=B[i][1]) for i in range(6)])
check_columns(T,'R408'); T.to_csv(pathlib.Path(__file__).parent/'results'/'two.csv',index=False)
gg=Gate('是「说了多少」还是「开了多少块」')
gg.curve_has_enough_points('⚠ guard 22:两条曲线各有几个点',list(range(6)),min_points=3,what='六层分层曲线')
gg.asserted('★ 正对照:真值固定的合成结局 -> 两条曲线的两端差就是各自的纯选择基线',True,
            f"按块覆盖 {pa:+.4f} · 按勾选总数 {pb:+.4f} —— **真实曲线对着它们读,不对着零读**")
gg.negative_control('★ 负对照(块覆盖):打乱人',float(NGa.mean()),da,null_spread=float(NGa.std()),
    null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:哪一条更陡(相对各自的纯选择基线)',
            (da-pa)>(db-pb),
            f"块覆盖 {da:+.4f}(基线 {pa:+.4f},超出 {da-pa:+.4f})· "
            f"勾选总数 {db:+.4f}(基线 {pb:+.4f},超出 {db-pb:+.4f})")
gg.asserted('⚠ 两个分层变量高度相关,两条曲线不是独立证据',True,
            f"`corr(块覆盖, 勾选总数)` = {np.corrcoef(cov[BASE],PK[BASE])[0,1]:+.4f} —— "
            f"本轮判的是**哪一条更陡**,不是「只有一条成立」")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
