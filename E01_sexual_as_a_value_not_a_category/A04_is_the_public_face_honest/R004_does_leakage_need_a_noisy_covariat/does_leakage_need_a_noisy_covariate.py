import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A23 R211 -- 泄漏需要的是"协变量有噪声",还是"泛函不匹配"?

`#165` 的校准表其实在说一件我当时没读出来的事:
    强度 0.00(纯置换)-> 去位置后内容 +0.079
    强度 0.05(位置信度匹配真实 0.60)-> +0.076      <- 没有泄漏
    强度 0.25 -> +0.263    强度 0.60 -> +0.255      <- 泄漏在这里才出现
**在真实数据所在的那个点上,泄漏是零。** 那强度大时那 0.26 是什么?两个世界:

    NOISE      衰减泄漏:位置分数自带测量误差,回归掉一个有噪声的协变量必然留残差。
               -> 用**无噪声的真值**(种植向量 u 本身)去回归,泄漏应当消失。
    FUNCTIONAL 泛函不匹配:被读的是选项矩阵的**主成分投影**,被减的是**平均冷门程度**——
               同一个信号的两个不同汇总。强度大时种植主导了 PC1,而减掉平均值减不掉投影。
               -> 用无噪声真值回归,泄漏**仍在**,因为差的不是精度是形状。

ESTIMAND        在同一个"内容为零"的合成世界里,固定种植强度,只换**回归时用哪个协变量**:
                (a) 无噪声真值 u  vs  (b) 半边估计出来的位置分数。读残差化后的内容信度。
IDENTIFICATION  背景是每块独立跨人置换 -> 构造上没有共享人层内容。任何非零读数都是泄漏。
SCOPE           k ∈ {6,8} x 3 个劈分种子 x 强度 {0.05, 0.25, 0.60} x 协变量 {真值, 估计}。
KILL            条件式:先要**在强度 0.60 处读到明显的泄漏**(否则没有可判的对象);
                再要**强度 0.05 处两种协变量都塌到置换基线**(否则仪器在低强度就有偏)。
                两者都满足时才判 NOISE vs FUNCTIONAL。
NEGATIVE CTRL   强度 0:置换背景本身。
NOISE FLOOR     3 个劈分种子。
IMPOSSIBLE      判不了真实数据里"该用哪个泛函" —— 真实数据没有真值 u。
                本轮只能判**机制**,机制再决定 `#104`/`#164` 要不要重算。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A22_is_rare_affinity_the_right_name'
          /'R210_how_big_is_the_content_side'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- 泄漏底噪校准')[0].split('ARMS=')[0])

def resid_on(x,z):
    m=np.isfinite(x)&np.isfinite(z); out=np.full_like(x,np.nan)
    X=np.c_[np.ones(m.sum()),z[m]]
    out[m]=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]
    return out

rows=[]
print(f"\n{'k':<4}{'seed':<6}{'强度':<8}{'协变量':<10}{'n':>7}{'残差化后内容':>14}{'位置 SB':>10}")
for k in [6,8]:
    for sd_ in range(6):
        rg=np.random.default_rng(zlib.crc32(f'split{k}{sd_}'.encode())%(1<<30))
        p=rg.permutation(len(qs)); hA=[qs[i] for i in p[:k]]; hB=[qs[i] for i in p[k:2*k]]
        u=rg.standard_normal(N)
        for st in [0.05,0.25,0.60]:
            cA,pA=build(hA,u,st,'position',True,sd_); cB,pB=build(hB,u,st,'position',True,sd_)
            cov_ok=np.where(np.isfinite(cA),u,np.nan)        # 无噪声真值,只在有分的人上
            covB  =np.where(np.isfinite(cB),u,np.nan)
            s_p,n=sb(pA,pB)
            for cname,zA,zB in (('真值 u',cov_ok,covB),('估计位置',pA,pB)):
                s,_=sb(resid_on(cA,zA),resid_on(cB,zB))
                rows.append(dict(k=k,seed=sd_,strength=st,covariate=cname,n=n,resid=s,position=s_p))
                print(f"{k:<4}{sd_:<6}{st:<8}{cname:<10}{n:>7,}{s:>+14.4f}{s_p:>+10.4f}",flush=True)
        # 强度 0 的基线
        cA,pA=build(hA,u,0.,'position',True,sd_); cB,pB=build(hB,u,0.,'position',True,sd_)
        s,n=sb(cA,cB)
        rows.append(dict(k=k,seed=sd_,strength=0.0,covariate='未回归',n=n,resid=s,position=sb(pA,pB)[0]))

T=pd.DataFrame(rows)
# ⚠ #166:第一版这列叫 `cov` —— DataFrame 自带 `.cov()` 方法,`T.cov=='未回归'`
#        比较的是方法对象。这正是 #117e 记录、check_columns 为之而写的 bug,而我没调用它。
check_columns(T,'R211 结果表')
T.to_csv(pathlib.Path(__file__).parent/'results'/'leak_mechanism.csv',index=False)
piv=T.pivot_table(index='strength',columns='covariate',values='resid',aggfunc='mean')
sdv=T.pivot_table(index='strength',columns='covariate',values='resid',aggfunc='std')
print("\n均值:"); print(piv.round(4).to_string())
print("展布:"); print(sdv.round(4).to_string())

base=float(T[T.covariate=='未回归'].resid.mean()); bsd=float(T[T.covariate=='未回归'].resid.std())
hi=T[(T.strength==0.60)]; lo=T[(T.strength==0.05)]
hi_t=float(hi[hi.covariate=='真值 u'].resid.mean()); hi_e=float(hi[hi.covariate=='估计位置'].resid.mean())
lo_t=float(lo[lo.covariate=='真值 u'].resid.mean()); lo_e=float(lo[lo.covariate=='估计位置'].resid.mean())
spread=float(hi.resid.std())
print(f"\n  置换基线(未回归) {base:+.4f} ± {bsd:.4f}")
print(f"  强度 0.60:真值 u {hi_t:+.4f}   估计位置 {hi_e:+.4f}")
print(f"  强度 0.05:真值 u {lo_t:+.4f}   估计位置 {lo_e:+.4f}")

g=Gate('泄漏的机制:噪声还是泛函不匹配')
g.asserted('可判前提一:强度 0.60 处确有泄漏可判',hi_e>base+2*bsd,
           f"估计位置 {hi_e:+.4f} vs 基线 {base:+.4f} ± {bsd:.4f}")
g.asserted('可判前提二:强度 0.05 处两种协变量都回到基线附近',
           abs(lo_e-base)<2*bsd and abs(lo_t-base)<2*bsd,
           f"真值 {lo_t:+.4f} / 估计 {lo_e:+.4f} vs 基线 {base:+.4f} ± {bsd:.4f}")
g.offset_control('无噪声真值能不能杀掉泄漏',hi_t,base,spread,
                 null_kind='同一合成世界的置换基线(未做任何回归)—— 构造上内容为零')
g.negative_control('强度 0 的置换基线',abs(base),max(hi_e,hi_t))
print(g)
if hi_t>base+2*bsd:
    print("\n  => FUNCTIONAL:无噪声真值也杀不掉泄漏,差的是**形状**不是精度")
else:
    print("\n  => NOISE:无噪声真值杀掉了泄漏,机制是协变量的测量误差")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 配对:两种协变量共用同一批块、同一个 u -------------------------------
# 边际比较的展布来自 PC1 方向在小 k 上的不稳定,而那份不稳定**两个臂共有**。
# 分岔的正确统计量是 (估计位置 - 真值 u) 的**配对差**:若机制是 NOISE,
# 无噪声真值应当明显多去掉一截;若是 FUNCTIONAL,两者应当**等价**。
print("\n---- 配对差:估计位置 - 真值 u ----")
W=T[T.covariate!='未回归'].pivot_table(index=['k','seed','strength'],columns='covariate',values='resid')
W['diff']=W['估计位置']-W['真值 u']
print(W.round(4).to_string())
for st in [0.05,0.25,0.60]:
    d=W.xs(st,level='strength')['diff']
    print(f"  强度 {st}: 配对差 {d.mean():+.4f} ± {d.std():.4f}  (n={len(d)})")

d60=W.xs(0.60,level='strength')['diff']; m60,s60=float(d60.mean()),float(d60.std())
u60=float(W.xs(0.60,level='strength')['真值 u'].mean())
g2=Gate('无噪声真值比有噪声估计多去掉了多少')
g2.resolvable('配对差本身可分辨',m60,s60)
g2.equivalent_within('两种协变量去掉的一样多',m60,s60,margin=0.10)
g2.asserted('而真值 u 留下的泄漏是它自己的量级',u60>0.15,f"真值 u 残留 {u60:+.4f}")
print(g2)
verdict=('FUNCTIONAL' if (abs(m60)+2*s60<=0.10 and u60>0.15) else
         'NOISE' if m60>2*s60 else 'UNVERIFIED')
print(f"\n  => {verdict}")

# ---- 弧真正要的那个判断:真实数据所在的强度上,有没有泄漏 --------------------
# 机制(NOISE vs FUNCTIONAL)两次都不可分辨 -> 按 frontier §3 不追第三轮。
# 但决定 `#104`/`#164` 去留的不是机制,是**真实数据所在的那个工作点上泄漏是多少**。
# `#165` 已定:真实位置信度 0.60 对应种植强度 0.05。这里做同劈分的配对比较。
print("\n---- 工作点判断:强度 0.05 vs 置换基线,同劈分配对 ----")
B=T[T.covariate=='未回归'].set_index(['k','seed'])['resid']
out=[]
for st in [0.05,0.60]:
    for cn in ['真值 u','估计位置']:
        v=W.xs(st,level='strength')[cn]
        d=(v-B).dropna()
        out.append(dict(strength=st,covariate=cn,paired_mean=float(d.mean()),
                        paired_sd=float(d.std()),n=int(len(d))))
P=pd.DataFrame(out); print(P.round(4).to_string(index=False))
P.to_csv(pathlib.Path(__file__).parent/'results'/'paired_vs_baseline.csv',index=False)

r=lambda st,cn:P[(P.strength==st)&(P.covariate==cn)].iloc[0]
g3=Gate('真实工作点上有没有泄漏')
a=r(0.05,'估计位置'); b=r(0.60,'估计位置'); c=r(0.60,'真值 u')
g3.asserted('可判前提:强度 0.60 处泄漏确实可分辨(否则无从谈"没有")',
            b.paired_mean>2*b.paired_sd,f"{b.paired_mean:+.4f} ± {b.paired_sd:.4f}")
g3.equivalent_within('工作点(强度 0.05)上泄漏 = 0',a.paired_mean,a.paired_sd,margin=0.05)
g3.resolvable('强度 0.60 处无噪声真值仍留下的泄漏',c.paired_mean,c.paired_sd)
print(g3)
print(f"\n  工作点泄漏 {a.paired_mean:+.4f} ± {a.paired_sd:.4f} —— `#104`/`#164` 的保留率"
      f"{'不需要' if abs(a.paired_mean)+2*a.paired_sd<=0.05 else '需要'}泄漏校正")

# ---- 两个不同的问题,两个不同的分母(#166b) --------------------------------
# Q1「**这一次**残差化引入了多少泄漏?」-> `#104`/`#164` 各自是一次实现 -> 分母 = 实现展布。
# Q2「残差化这个动作**平均**会不会引入泄漏?」-> 分母 = 均值的标准误 = sd/sqrt(n)。
# 我原来只发表了 Q1 的答案,而弧要判的方法学问题是 Q2。两个都发表。
print("\n---- Q1 单次 vs Q2 平均 ----")
print(f"{'强度':<8}{'协变量':<10}{'配对均值':>10}{'实现展布':>10}{'均值SE':>9}{'Q1':>12}{'Q2':>12}")
for _,r_ in P.iterrows():
    se=r_.paired_sd/np.sqrt(r_.n)
    q1='可分辨' if abs(r_.paired_mean)>2*r_.paired_sd else '不可分辨'
    q2='可分辨' if abs(r_.paired_mean)>2*se else '不可分辨'
    print(f"{r_.strength:<8}{r_.covariate:<10}{r_.paired_mean:>+10.4f}{r_.paired_sd:>10.4f}"
          f"{se:>9.4f}{q1:>12}{q2:>12}")
a2=r(0.05,'估计位置'); se2=a2.paired_sd/np.sqrt(a2.n)
g4=Gate('Q2:残差化平均会不会引入泄漏')
g4.asserted('可判前提:强度 0.60 处在均值尺度上可分辨',
            r(0.60,'估计位置').paired_mean>2*r(0.60,'估计位置').paired_sd/np.sqrt(r(0.60,'估计位置').n),
            f"{r(0.60,'估计位置').paired_mean:+.4f}, SE {r(0.60,'估计位置').paired_sd/np.sqrt(r(0.60,'估计位置').n):.4f}")
g4.equivalent_within('工作点上平均泄漏 = 0',float(a2.paired_mean),float(se2),margin=0.05)
print(g4)

# ---- 机制,在均值尺度上重判(#166c) ---------------------------------------
# Q1 尺度上两次都不可分辨;但机制问的是"这个动作平均怎么表现",属于 Q2。
# ⚠ 边界 0.10 是本脚本在**三种子阶段就写下的**注册值。六种子后我一度想改成 0.05
#   (那样会变成 UNVERIFIED)—— 那是事后收紧,不算数。两个都写在这里。
d60=W.xs(0.60,level='strength')['diff']; se60=float(d60.std()/np.sqrt(len(d60)))
u60m=float(W.xs(0.60,level='strength')['真值 u'].mean())
print(f"\n---- 机制(均值尺度)----\n  强度 0.60:配对差 {d60.mean():+.4f}  SE {se60:.4f}  n={len(d60)}"
      f"\n  无噪声真值残留 {u60m:+.4f}")
g5=Gate('机制:NOISE 还是 FUNCTIONAL(均值尺度)')
g5.asserted('可判前提:无噪声真值仍留下量级可观的泄漏',u60m>0.15,f"真值 u 残留 {u60m:+.4f}")
g5.equivalent_within('两种协变量去掉的一样多(注册边界 0.10)',float(d60.mean()),se60,margin=0.10)
print(g5)
print(f"  [事后收紧到 0.05 会是:界 {abs(d60.mean())+2*se60:.4f} > 0.05 -> UNVERIFIED。不采纳。]")
print(f"\n  => {'FUNCTIONAL' if (abs(d60.mean())+2*se60<=0.10 and u60m>0.15) else 'UNVERIFIED'}"
      f":一个**完美的**协变量去掉的并不比有噪声的多,而它自己仍留下 {u60m:+.3f} —— 差的是形状不是精度")
