import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A24 R212 -- 三条头部声明,各自的两个分母

`#166b`:Q1(这一次)用实现展布,Q2(平均)用均值 SE。README 上的声明读起来是 Q1,
判它们时用的却是 Q2 的分母。挑三条头部声明,**用各轮自己的代码**跨 8 个种子重跑,
只加分母、不重写统计量(`#143`/`#154`:轮次自己的代码才是那台仪器)。

    #100  跨不相交块集的分半信度 +0.432        (E01·A11·R147, rel_resid)
    #159  「何时」的结构:题目侧 0.96 vs 人侧 0.105  (E01·A21·R204, edge_item/edge_person)
    #163  内容最不相似分半保留 85%             (E01·A22·R208, rel dissimilar/random)

ESTIMAND        每条声明的点值,加上**两个**分母:实现展布 sd(Q1)与均值 SE = sd/sqrt(n)(Q2)。
KILL            条件式:先要**每条声明在 Q2 尺度上仍可分辨**(否则连平均都不成立,是另一个问题);
                再看有几条在 Q1 尺度上掉出来。**三条里 >=2 条掉出 Q1 -> README 加一列。**
IMPOSSIBLE      判不了"哪个分母才对" —— 那取决于声明想说什么。本轮只把两个都摆出来。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

def both(name,vals,floor=0.0,contrast=''):
    v=np.asarray([x for x in vals if np.isfinite(x)],float)
    sd=float(v.std(ddof=1)); se=sd/np.sqrt(len(v)); m=float(v.mean())
    # ⚠ 列名 `mean` 与 DataFrame.mean 撞名 —— check_columns 抓到了(#166c 的下一轮)。
    return dict(claim=name,n_seeds=len(v),point=m,floor=floor,effect=m-floor,
                sd_realization=sd,se_mean=se,
                q1_ratio=abs(m-floor)/sd if sd>0 else np.inf,
                q2_ratio=abs(m-floor)/se if se>0 else np.inf,contrast=contrast)

rows=[]

# ---- #163 (R208):内容最不相似分半保留率 -----------------------------------
print("== #163  E01·A21·R208 ==",flush=True)
S=(ROOT/'E01_sexual_as_a_value_not_a_category/A22_is_rare_affinity_the_right_name'
       /'R208_position_or_content'/'run.py').read_text()
exec(S.split('"""',2)[2].split('KS=')[0])
ret=[]
# ⚠ 账本里的 85% 是 k∈{4,5,6} 的**平均**(82.2/83.1/88.4),不是单个 k。
#   我第一版只跑了 k=6(88.1%),那是另一个对象 —— #147/#151/#161 的 same_scale 家族。
for sd_ in range(8):
    per=[]
    for kk in [4,5,6]:
        a=rel('real',kk,900+sd_,'dissimilar'); b=rel('real',kk,900+sd_,'random')
        per.append((a[0] if isinstance(a,tuple) else a)/(b[0] if isinstance(b,tuple) else b))
    ret.append(float(np.mean(per))); print(f"  seed {sd_}: 保留 {ret[-1]:.4f}  (k4/5/6 {per[0]:.3f}/{per[1]:.3f}/{per[2]:.3f})",flush=True)
rows.append(both('#163 内容不相似分半保留率',ret,floor=0.0,contrast='对 0;另见下方对 1.0'))
pd.DataFrame(dict(seed=range(8),retention=ret)).to_csv(
    pathlib.Path(__file__).parent/'results'/'c163_retention.csv',index=False)

# ---- #159 (R204):「何时」的结构在题目侧还是人侧 ----------------------------
print("\n== #159  E01·A21·R204 ==",flush=True)
S2=(ROOT/'E01_sexual_as_a_value_not_a_category/A21_who_owns_the_personal_20_percent'
        /'R204_person_or_item'/'run.py').read_text()
g=dict(globals()); exec(S2.split('"""',2)[2].split('for tag,Dx,perm in')[0],g)
ei=[];ep=[]
for sd_ in range(8):
    ri=g['edge_item'](g['D0'],4000+sd_); rp=g['edge_person'](g['D0'],5000+sd_)
    ri=ri[0] if isinstance(ri,tuple) else ri; rp=rp[0] if isinstance(rp,tuple) else rp
    ei.append(float(ri)); ep.append(float(rp)); print(f"  seed {sd_}: 题目侧 {ri:+.4f}  人侧 {rp:+.4f}",flush=True)
rows.append(both('#159 题目侧 − 人侧(配对差)',list(np.array(ei)-np.array(ep)),floor=0.0,
                 contrast='声明说的是两侧的差,不是各自对 0'))
rows.append(both('#159 人侧相关',ep,floor=0.0,contrast='对 0'))
pd.DataFrame(dict(seed=range(8),item=ei,person=ep)).to_csv(
    pathlib.Path(__file__).parent/'results'/'c159_edges.csv',index=False)

# ---- #100 (R147):跨不相交块集的分半信度 -----------------------------------
# 用 R147 **自己的源码**跑,只把 seed 范围放宽;并把 __file__ 指向本轮,
# 这样它的 to_csv 写进本轮的 results/,**不覆盖 R147 的 artifact**(L81)。
print("\n== #100  E01·A11·R147 ==",flush=True)
S3=(ROOT/'E01_sexual_as_a_value_not_a_category/A11_can_a_minority_structure_be_seen'
        /'R147_is_it_a_trait'/'run.py').read_text()
# ⚠ R147 的内层劈分种子是 `default_rng(700+rep)` —— **与外层 sd 无关**,所以真实臂
#   在 8 个种子上逐字节相同(sd=0.0000)。它的 −0.022 是 curveball **零假设的地板**,
#   不是 0.432 自己的误差。**这条声明从来没有过精度估计。**
#   两处文本改动(改的是驱动,不是统计量):把内层种子挂到 sd 上;把 real 臂的 a/b/picks 存出来。
S3=S3.replace('for sd in range(1,4):','for sd in range(1,9):')
S3=S3.replace('rr=np.random.default_rng(700+rep)','rr=np.random.default_rng(700+rep+1000*sd)')
S3=S3.replace("        R=np.array(rs)","        R=np.array(rs)\n        if w=='real': CAP.append((a.copy(),b.copy(),picks.copy()))")
CAP=[]
g3={'__file__':str(pathlib.Path(__file__)),'__name__':'r147',
    'pathlib':pathlib,'os':os,'sys':sys,'CAP':CAP}   # 源码的 import 在被切掉的文档字符串之前
exec(S3.split('"""',2)[2],g3)
G=g3['D']; G=G[G.world=='real']
cb=float(g3['D'].query("world=='cb'").rel_resid.mean())
rows.append(both('#100 跨块分半信度(去勾选数)',list(G.rel_resid),floor=cb,
                 contrast=f'curveball 地板 {cb:+.4f}'))

# 人层自助:用最后一个 sd 捕获的 a/b/picks(是 R147 自己算出来的量,不是我重算的)
a_,b_,pk=CAP[-1]
def rel_resid_of(ix):
    A,B,P=a_[ix],b_[ix],pk[ix]
    X=np.c_[np.ones(len(P)),P,np.log(P)]
    ra=A-X@np.linalg.lstsq(X,A,rcond=None)[0]; rb=B-X@np.linalg.lstsq(X,B,rcond=None)[0]
    r=np.corrcoef(ra,rb)[0,1]; return 2*r/(1+r)
rgb=np.random.default_rng(20260803); n_=len(a_)
boot=[rel_resid_of(rgb.integers(0,n_,n_)) for _ in range(400)]
print(f"  人层自助(400): {np.mean(boot):+.4f} ± {np.std(boot):.4f}   (n={n_:,} 人)")
rows.append(dict(claim='#100 同上,人层自助',n_seeds=400,point=float(np.mean(boot)),floor=cb,
                 effect=float(np.mean(boot))-cb,sd_realization=float(np.std(boot)),
                 se_mean=float(np.std(boot)),   # 自助分布的 sd **就是** 单次估计的 SE
                 q1_ratio=abs(np.mean(boot)-cb)/np.std(boot),
                 q2_ratio=abs(np.mean(boot)-cb)/np.std(boot),
                 contrast=f'curveball 地板 {cb:+.4f};自助 sd 同时是 Q1 与 Q2'))
G.to_csv(pathlib.Path(__file__).parent/'results'/'c100_grid.csv',index=False)

T=pd.DataFrame(rows); check_columns(T,'R212 汇总')
T.to_csv(pathlib.Path(__file__).parent/'results'/'two_denominators.csv',index=False)
print("\n"+T.round(4).to_string(index=False))
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 注册的 kill,以及它没抓到的那件事 --------------------------------------
q1=T[T.claim.str.contains('#163|配对差|自助')]
drop=int((q1.q1_ratio<2).sum())
print(f"\n三条头部声明在 Q1 尺度上掉出来的数量:{drop}/3  (注册 kill:>=2 -> README 加一列分母)")
g=Gate('README 那张声明表要不要加一列分母')
g.asserted('可判前提:每条声明在 Q2 尺度上仍可分辨',bool((q1.q2_ratio>2).all()),
           ' / '.join(f'{v:.0f}x' for v in q1.q2_ratio))
g.asserted('注册的 kill 不开火',drop<2,f"{drop}/3 掉出 Q1;最差 {q1.q1_ratio.min():.1f}x")
g.resolvable('#100 在人层自助尺度上',float(T[T.claim.str.contains('自助')].effect.iloc[0]),
             float(T[T.claim.str.contains('自助')].sd_realization.iloc[0]))
g.no_sign_crossing('#100 三种算法给出的点值',
                   [float(T[T.claim.str.contains('#100')].point.iloc[0]),
                    float(T[T.claim.str.contains('自助')].point.iloc[0]),0.432])
print(g)
pub,bs=0.432,float(T[T.claim.str.contains('自助')].point.iloc[0])
bsd=float(T[T.claim.str.contains('自助')].sd_realization.iloc[0])
print(f"\n  已发表点值 {pub:.4f} · 自助均值 {bs:.4f} · 差 {pub-bs:+.4f} = {(pub-bs)/bsd:.1f} 个自助 sd")
print(f"  => `#100` 的正确写法是 **0.432 ± 0.016(人层自助)**,而不是 0.432(地板 −0.022)")
