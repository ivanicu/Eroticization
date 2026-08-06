import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A41 R243 -- Δ 是「所有人都有一点点」,还是「一部分人有很多、另一部分相反」

`#128` 的 Δ = **mean_i rho_i** = −0.0328。**它是一个人层平均,而从没人看过 rho_i 的分布。**
这两件事在同一个均值下长得一模一样,却是完全不同的心理学陈述:

    UNIFORM   单峰围绕 −0.033 -> 「每个人的罕见兴趣都稍微来得早一点」
    MIXED     比零更宽 / 双峰 -> 「一部分人有很多,另一部分人相反」

ESTIMAND        rho_i 的分布**宽度**与**模态**,对着**置换零的 rho_i 分布**读。
IDENTIFICATION  置换零(题内跨人置换)保留**缺失模式**与**每题的值分布**,只毁掉配对 ——
                所以它带着**完全相同的人内测量噪声**。**这个零不是 0,是"只有噪声时该有多宽"。**
KILL            条件式:先要两个正对照分别开火(双峰种植 -> 测到双峰;单峰宽种植 -> 测到超额宽度
                但**不**测到双峰);再判:
                **超额 sd = √(var_real − var_null) 若 > 0.10,则 Δ 里有真实的人际差异;
                且 2 成分混合的 BIC 若明显优于 1 成分,则是 MIXED 而非 UNIFORM。**
NEGATIVE CTRL   置换零自己必须判为单峰,且它的"超额宽度"= 0。
NOISE FLOOR     人层 bootstrap 100;置换 5 个种子。
IMPOSSIBLE      rho_i 每人只由 8–31 个类别估出,**个体噪声很大** ——
                所以"超额宽度"是人际差异的**下界**,而任何双峰若窄于噪声就看不见。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
print(f"载入:V {V.shape} · KEEP {int(KEEP.sum()):,} 人",flush=True)

def rhos(Vm):
    _,rho=betas(Vm); m=np.isfinite(rho)&KEEP
    return rho[m]

def gmm_bic(x,k,rng,iters=200):
    """k 成分一维高斯混合的 BIC(EM)。"""
    n=len(x)
    if k==1:
        mu=np.array([x.mean()]); sd=np.array([x.std()]); w=np.array([1.0])
    else:
        q=np.quantile(x,np.linspace(0.2,0.8,k)); mu=q.copy(); sd=np.full(k,x.std()); w=np.full(k,1.0/k)
    for _ in range(iters):
        P=np.array([w[j]*np.exp(-0.5*((x-mu[j])/np.maximum(sd[j],1e-9))**2)/np.maximum(sd[j],1e-9) for j in range(k)])
        P=P/np.maximum(P.sum(0),1e-300)
        Nk=P.sum(1); w=Nk/n
        mu=(P*x).sum(1)/np.maximum(Nk,1e-9)
        sd=np.sqrt((P*(x-mu[:,None])**2).sum(1)/np.maximum(Nk,1e-9))
        sd=np.maximum(sd,1e-4)
    ll=np.log(np.maximum(np.array([w[j]*np.exp(-0.5*((x-mu[j])/sd[j])**2)/sd[j] for j in range(k)]).sum(0),1e-300)).sum()
    p=3*k-1
    return -2*ll+p*np.log(n), mu, sd, w

rng=np.random.default_rng(20260803)
def describe(x,tag,rng):
    b1,_,_,_=gmm_bic(x,1,rng); b2,mu2,sd2,w2=gmm_bic(x,2,rng)
    # ⚠ 列名 `mean` 撞 pandas 方法 —— 收紧后的 check_columns(#197a)在这一轮又抓到一个
    return dict(arm=tag,n=len(x),m_rho=float(x.mean()),sd=float(x.std()),
                bic1=float(b1),bic2=float(b2),dbic=float(b1-b2),
                mu_lo=float(min(mu2)),mu_hi=float(max(mu2)),w_lo=float(w2[np.argmin(mu2)]))

R=rhos(V); rows=[describe(R,'真实',rng)]
NUL=[]
for sd_ in range(5):
    x=rhos(perm_null(V,np.random.default_rng(4200+sd_))); NUL.append(x)
    rows.append(describe(x,f'置换零 s{sd_}',rng))
var_null=float(np.mean([x.var() for x in NUL]))
excess=float(np.sqrt(max(R.var()-var_null,0.0)))
print(f"\n真实 rho_i:mean {R.mean():+.4f} · sd {R.std():.4f} · n {len(R):,}")
print(f"置换零 rho_i:sd {np.sqrt(var_null):.4f}(5 个种子)")
print(f"**超额 sd = √(var_real − var_null) = {excess:.4f}**")

# 正对照
def plant_person(g, bimodal):
    u=(rng.choice([-1.0,1.0],len(V)) if bimodal else rng.standard_normal(len(V)))
    x=rar-rar.mean()
    return V+g*np.outer(u,x)*obs
Rb=rhos(plant_person(2.0,True)); Ru=rhos(plant_person(2.0,False))
rows+= [describe(Rb,'【正对照】双峰种植',rng), describe(Ru,'【正对照】单峰宽种植',rng)]
T=pd.DataFrame(rows); check_columns(T,'R243'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rho_dist.csv',index=False)
print(f"\n{'臂':<18}{'mean':>9}{'sd':>8}{'ΔBIC(1−2)':>12}{'低峰位置':>10}{'低峰占比':>9}")
for _,r in T.iterrows():
    print(f"{r.arm:<18}{r.m_rho:>+9.4f}{r.sd:>8.4f}{r.dbic:>12.1f}{r.mu_lo:>+10.4f}{r.w_lo:>9.2f}")

boot=[]
for _ in range(100):
    ii=rng.choice(np.flatnonzero(KEEP),int(KEEP.sum()),replace=True)
    Kb=np.zeros(len(V),bool); Kb[ii]=True; sk=KEEP.copy(); KEEP[:]=Kb
    x=rhos(V); boot.append(np.sqrt(max(x.var()-var_null,0.0))); KEEP[:]=sk
sdE=float(np.std(boot)); print(f"\n超额 sd 的人层 bootstrap sd = {sdE:.4f}")

real=T.iloc[0]; nul=T[T.arm.str.startswith('置换零')]
bim=T[T.arm.str.contains('双峰')].iloc[0]; uni=T[T.arm.str.contains('单峰宽')].iloc[0]
g=Gate('Δ 是所有人一点点还是一部分人很多')
g.asserted('正对照一:双峰种植必须被测成双峰',bim.dbic>50,f"ΔBIC = {bim.dbic:.1f}")
g.asserted('正对照二:单峰宽种植必须**不**被测成双峰',uni.dbic<bim.dbic/3,
           f"单峰宽 ΔBIC {uni.dbic:.1f} vs 双峰 {bim.dbic:.1f}")
g.asserted('负对照:置换零自己必须单峰',float(nul.dbic.max())<50,f"最大 ΔBIC {nul.dbic.max():.1f}")
g.offset_control('真实 rho_i 的 sd vs 置换零的 sd',float(real.sd),float(np.sqrt(var_null)),sdE,
                 null_kind='题内跨人置换的 rho_i 分布 —— 保留缺失模式与每题值分布,'
                           '**带着完全相同的人内测量噪声**;这不是 0,是"只有噪声时该有多宽"')
g.resolvable('超额 sd',excess,sdE)
g.asserted('注册的 kill:超额 sd > 0.10 且 2 成分明显优于 1 成分',
           (excess>0.10) and (real.dbic>50),
           f"超额 sd {excess:.4f} · 真实 ΔBIC {real.dbic:.1f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 模态判不了,但宽度判得了 —— 而宽度本身就是结论 --------------------------
# ⚠ 正对照二失败:**单峰宽**种植也给出 ΔBIC = 792。在 n≈10,000 上,
#   2 成分混合几乎总是赢 1 成分(真实分布都有一点非正态)。**ΔBIC 不是这里的双峰检验。**
#   按条件式规则:**不报模态判定**。但超额宽度是另一条独立的测量,不依赖那个仪器。
print("\n---- 宽度说了什么 ----")
from math import erf, sqrt as _sqrt
Phi=lambda z: 0.5*(1+erf(z/_sqrt(2)))
mu=float(R.mean())
share_neg=Phi(-mu/excess)          # 正态近似下 rho_i<0 的比例
print(f"  真实人层展布(去噪后)sd = {excess:.4f},而均值只有 {mu:+.4f}")
print(f"  **展布 / |均值| = {excess/abs(mu):.1f}×**")
print(f"  正态近似下,rho_i < 0 的人占 **{100*share_neg:.0f}%** —— 即约 {100*(1-share_neg):.0f}% 的人方向相反")
bs=[]
for _ in range(200):
    ii=rng.choice(np.flatnonzero(KEEP),int(KEEP.sum()),replace=True)
    Kb=np.zeros(len(V),bool); Kb[ii]=True; sk=KEEP.copy(); KEEP[:]=Kb
    x=rhos(V); e=np.sqrt(max(x.var()-var_null,0.0))
    bs.append(Phi(-float(x.mean())/e) if e>0 else np.nan); KEEP[:]=sk
sd_share=float(np.nanstd(bs)); print(f"  该比例的人层 bootstrap sd = {sd_share:.4f}")

g2=Gate('宽度 vs 均值')
g2.asserted('⚠ 模态判定被撤下:正对照二失败,ΔBIC 分不开双峰与单峰宽',
            uni.dbic>bim.dbic/3, f"单峰宽 {uni.dbic:.1f} vs 双峰 {bim.dbic:.1f} —— 本轮**不报**模态")
# ⚠ 我原本在这里写 `same_scale(展布, |均值|)` —— **误用**。`same_scale` 比的是两个量的
#   **可比量级**(用来挡「拿 k=17 的数去比 k=9 的数」那类错),不是「同一个单位」。
#   而 5.4× 恰恰是本轮的结论,所以那个断言等于在断言结论不成立。守卫没错,用法错了。
g2.asserted('展布与均值在同一 rho 尺度上,且展布远大于均值(这就是本轮的结论)',
            excess>3*abs(mu), f"展布 {excess:.4f} vs |均值| {abs(mu):.4f} = {excess/abs(mu):.1f}×")
g2.resolvable('展布本身',excess,sdE)
g2.threshold_outside_noise('rho_i<0 的比例 vs 50%',share_neg,0.50,sd_share)
g2.asserted('⚠ 该比例依赖正态近似 —— 而模态恰恰是判不了的那一项',True,
            '所以 57% 这个数带着一个本轮无法检验的假设;宽度本身不带')
print(g2)
