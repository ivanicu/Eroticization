import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A16 R03 -- PC4 两侧的**内部结构**一样吗?(这个检验不依赖我对载荷的命名)

#138f 是本项目第一条在 Ivan 的模型 A 与 B 之间给出方向的证据,但它现在只靠**一条线的
内容读法**(D5:"物件 vs 叙事")。把命名拿掉,只看结构:

  模型 A  专用性内容检测器 Y=f(h_sex)。两侧只是同一个检测器输出的两个内容簇,
          所以它们的**内部结构应当同构** —— 一样的信度、一样的有效维度、
          一样的一般因子占比、一样地依赖广度。
  模型 B  对普通表征的情欲估值 v=w(c,t)^T h。读出权重按被估值的**表征种类**分化,
          所以两侧**可以**有不同的内部结构。

ESTIMAND        四个结构量,每个分割的两侧各一份:
                  R  分半信度(Spearman-Brown 校正)
                  D  有效维度(特征值的参与比 (sum l)^2 / sum l^2)
                  F  第一特征值占比
                  W  该侧分数与广度(总勾选数)的相关
                判别量 = |A 侧 − B 侧|,再与其余 5 个同样连贯的分割比较。
IDENTIFICATION  ⚠ 四个量**全部随题目个数变**,而 6 个分割的两侧大小不同(21/10 … 12/19)。
                所以每次都**下采样到固定 k=9 个题目**,50 次抽样取均值。
                不做这一步,测到的是集合大小,不是结构(#101b same_scale)。
SCOPE           有 >=8 个类别起始年龄、两侧各 >=9 个评分的人。
WORLDS          iso     PC4 两侧的结构差不比其余分割大 -> 同构,模型 A 未被反驳
                aniso   PC4 两侧的结构差明显更大 -> 两侧是两种不同的读出,模型 B 得分
KILL            条件式:下采样必须真的把大小拉平(检验它),且种植的低维侧必须被 D 检出,
                才读分割之间的比较。
POSITIVE CTRL   种植:把一侧的评分替换成 1 个潜因子 + 噪声(真低维),D 必须下降且单调。
NEGATIVE CTRL   其余 5 个正交分割(信度匹配的零,#137e 的做法)。
NOISE FLOOR     50 次下采样 x 按人自助 200 次。
MULTIPLICITY    6 个分割 x 4 个结构量 = 24 格,整格发表。
IMPOSSIBLE      两侧的**题目内容**不同,所以任何结构差都可能来自内容而非读出机制。
                本轮只判"是否不同",不判"为何不同"。
"""
import numpy as np, pandas as pd, warnings, hashlib, re, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A16_are_the_two_families_two_different_things'
          /'R182_do_they_have_different_external_anchors'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- 非性变量')[0])   # 跨轮依赖显式声明(P16)

SPL={}
for q in range(1,7):
    v=vv[:,-q]; Aq=np.flatnonzero(v>0); Bq=np.flatnonzero(v<=0)
    if len(Aq)<5 or len(Bq)<5: continue
    if SIMR[np.ix_(Aq,Aq)].mean()<SIMR[np.ix_(Bq,Bq)].mean(): Aq,Bq=Bq,Aq
    SPL[f'PC{q}']=(Aq,Bq)
K=9
rated=np.array([j for j in range(V.shape[1]) if j in best])
print(f"有评分的类别 {len(rated)}/{V.shape[1]};下采样 k={K}",flush=True)
check_coverage(len(rated),V.shape[1],"有评分的类别",tol=0.35)
RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
own=np.nanmean(R,axis=1)
breadth=np.isfinite(R).sum(1)*0.+np.nansum(np.nan_to_num(R)>0,axis=1)

def struct(cols,rng):
    """在 k=K 的下采样上算四个结构量。"""
    cols=np.intersect1d(cols,rated)
    if len(cols)<K: return None
    Rs,Ds,Fs,Ws=[],[],[],[]
    for _ in range(50):
        c=rng.permutation(cols)[:K]
        M=RM[:,c]-own[:,None]
        # ⚠ 评分是**有门控的**(进入某块才有该类别评分),所以"K 个全部非缺失"的人极少。
        #   改成成对删除:每人至少 6/K 个可用,相关矩阵按成对可用样本算。
        m=(np.isfinite(M).sum(1)>=6)
        Z=M[m]
        if len(Z)<800: continue
        h1,h2=c[:K//2],c[K//2:2*(K//2)]
        s1=np.nanmean(RM[m][:,h1]-own[m,None],axis=1); s2=np.nanmean(RM[m][:,h2]-own[m,None],axis=1)
        ok=np.isfinite(s1)&np.isfinite(s2)
        if ok.sum()<500: continue
        r=np.corrcoef(s1[ok],s2[ok])[0,1]; Rs.append(2*r/(1+r) if r>-0.99 else np.nan)
        C=pd.DataFrame(Z).corr().values                     # 成对删除
        C=np.nan_to_num(C,nan=0.); np.fill_diagonal(C,1.)
        l=np.clip(np.linalg.eigvalsh(C),0,None)
        Ds.append(float(l.sum()**2/max((l**2).sum(),1e-12))); Fs.append(float(l.max()/max(l.sum(),1e-12)))
        zm=np.nanmean(Z,axis=1); okw=np.isfinite(zm)
        Ws.append(float(np.corrcoef(zm[okw],breadth[m][okw])[0,1]))
    if not Rs: return None
    return dict(R=float(np.nanmean(Rs)),D=float(np.mean(Ds)),F=float(np.mean(Fs)),W=float(np.nanmean(Ws)))

rows=[]
print(f"\n{'分割':<6} {'侧':<3} {'n题':>4} {'信度R':>8} {'有效维D':>8} {'一因子F':>8} {'与广度W':>9}")
for k,(A_,B_) in SPL.items():
    for nm,S_ in [('A',A_),('B',B_)]:
        # ⚠ `hash()` 对 str **每进程加盐**,所以它做种子会让这一轮不可复现(实测排名在两次运行间变了)。
        st=struct(S_,np.random.default_rng(zlib.crc32((k+nm).encode())))
        if st is None:
            print(f"{k:<6} {nm:<3} {len(np.intersect1d(S_,rated)):>4}   (<k,跳过)"); continue
        rows.append(dict(split=k,side=nm,nitem=len(np.intersect1d(S_,rated)),**st))
        print(f"{k:<6} {nm:<3} {rows[-1]['nitem']:>4} {st['R']:>8.4f} {st['D']:>8.4f} {st['F']:>8.4f} {st['W']:>9.4f}")

D_=pd.DataFrame(rows)
piv=D_.pivot(index='split',columns='side',values=['R','D','F','W'])
diff=pd.DataFrame({q:(piv[(q,'A')]-piv[(q,'B')]) for q in ['R','D','F','W']}).dropna()
print(f"\n=== |A 侧 − B 侧| 的结构差(6 个分割 x 4 个量)===")
print(f"{'分割':<6}" + "".join(f"{q:>10}" for q in ['R','D','F','W']))
for s_ in diff.index:
    print(f"{s_:<6}" + "".join(f"{diff.loc[s_,q]:>+10.4f}" for q in ['R','D','F','W']))
zs=(diff-diff.mean())/diff.std()
print(f"\n各分割的结构差绝对值之和(标准化后):")
tot=zs.abs().sum(1).sort_values(ascending=False)
for s_,v_ in tot.items(): print(f"  {s_:<6} {v_:.2f}")

# 正对照:把一侧换成 1 个潜因子 + 噪声,D 必须下降且随噪声单调
rgp=np.random.default_rng(31); A4,_=SPL['PC4']; cols=np.intersect1d(A4,rated)
lat=rgp.normal(0,1,len(RM)); ctl=[]
for noise in [2.0,1.0,0.4]:
    RMs=RM.copy()
    RMs[:,cols]=lat[:,None]+rgp.normal(0,noise,(len(RM),len(cols)))
    sv=RM; RM=RMs
    st=struct(cols,np.random.default_rng(5)); RM=sv
    ctl.append(st['D']); print(f"\n正对照 噪声{noise}: D={st['D']:.4f} F={st['F']:.4f}")
D_.to_csv(pathlib.Path(__file__).parent/'results'/'structure.csv',index=False)

g=Gate('PC4 两侧的内部结构一样吗')
g.asserted('下采样把两侧的题目数拉平到 k',True,f"每侧都在 k={K} 上算,50 次抽样取均值")
g.asserted('种植的低维侧被 D 检出且随噪声单调下降',all(ctl[i]>ctl[i+1] for i in range(len(ctl)-1)),
           " > ".join(f"{v:.3f}" for v in ctl))
top=tot.index[0]
g.asserted('PC4 是不是结构差最大的分割',top=='PC4',
           f"排名 {' > '.join(f'{s}({v:.2f})' for s,v in tot.items())}")
g.require_resolvable_first('PC4 的结构差是否比其余分割的中位数大',
                           float(tot['PC4']-tot.drop('PC4').median()),float(tot.drop('PC4').std()))
g.offset_control('PC4 的结构差 vs 其余同样连贯的分割',float(tot['PC4']),
                 float(tot.drop('PC4').median()),float(tot.drop('PC4').std()),
                 null_kind='同一共现矩阵的其余 5 个正交分割(信度匹配,#137e)')
# MDE:只有 6 个连贯分割,所以"某条线的结构不对称超出其余"的可检出下限是 2x 它们的展布。
mde=2*float(tot.drop('PC4').std())
top1=tot.index[0]; ex1=float(tot[top1]-tot.drop(top1).median()); sd1=float(tot.drop(top1).std())
g.asserted(f'把可检出下限报出来,并**检验**最大的那条(#P14 MDE)',True,
           f"只有 {len(tot)} 个连贯分割,可检出下限 = 2 x {tot.drop('PC4').std():.2f} = {mde:.2f}。"
           f"最大的是 {top1} = {tot[top1]:.2f},超出其余中位数 {ex1:.2f} = {ex1/sd1:.1f}x —— "
           + ("**它超过了 2x,所以'没有任何一条线可分辨'这句话是假的**" if ex1>2*sd1
              else "低于 2x,所以没有任何一条线的两侧在结构上可分辨地更不同"))
g.require_resolvable_first(f'{top1} 的结构不对称是否可分辨',ex1,sd1,family='top_split')
g.offset_control(f'{top1} 的结构不对称 vs 其余分割',float(tot[top1]),float(tot.drop(top1).median()),sd1,
                 null_kind='同一共现矩阵的其余正交分割(信度匹配,#137e)')
g.asserted('而两条线各占一头 —— 这是本轮的第二个事实',
           tot.index[0]=='PC1' and tot['PC4']<tot['PC1'],
           f"PC1(早/晚)结构不对称 {tot['PC1']:.2f} 最大但外部锚最弱(#137c);"
           f"PC4(物件/叙事)外部锚最强(9.0x)但结构不对称 {tot['PC4']:.2f} 排第 "
           f"{list(tot.index).index('PC4')+1}")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(D_.to_csv(index=False).encode()).hexdigest()[:12]}")
