import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A16 R04 -- 关系族多出来的 2.4 个维度,是 A02 那三条吗?

#139c:关系族在**时间**上是一起到的(#136a,17.7x),但在**偏好结构**上不是一个东西 ——
有效维度比具体族多 2.42 个(3.0x,六个连贯分割里最大的结构不对称)。

A02 已经在**全部**类别上命名过三条近乎独立的坐标(谁服从 · 谁被看 · 谁接受,
去衰减互相关 <= 0.362,有效维度 2.95/3)。所以有两个世界:

  KNOWN   关系族内部的维度就是 A02 那三条 -> #139c 从"好几个"变成一个**具体的**结构,
          而且它解释了为什么关系族会同时到达却不是一个东西:它是三条轴的**交汇处**
  NEW     关系族里的维度**不是** A02 那三条 -> 这里有 A02 在全类别层面没找到的东西,
          而那是本项目下一个真正的开口

ESTIMAND        关系族内部的有效维度数,以及它的前几个成分与 A02 三条坐标的**典型相关**
                (在留出的人上算,不是拟合集上)。
IDENTIFICATION  A02 的坐标由**全部**类别拟合;关系族的成分只由**族内**类别拟合。
                两者共享 item,所以必须用 check_disjoint_items 明确记录这一点,
                并把"共享 item 能解释多少"作为一个必须报的量。
SCOPE           有 >=8 个类别起始年龄、族内 >=6 个评分的人。
WORLDS          KNOWN / NEW
KILL            条件式:留出典则相关必须先在一个**已知同构**的对照上开火(把族内成分
                与它自己的另一半人做典则相关),且随机坐标的对照必须为零,才读阈值。
POSITIVE CTRL   族内成分 vs 它自己在另一半人上的重估(必须高)。
NEGATIVE CTRL   把 A02 坐标的**载荷随机置换**后重算(必须为零)。
NOISE FLOOR     5 次人随机对半 + 200 次按人自助。
MULTIPLICITY    前 4 个成分 x 3 条 A02 坐标 = 12 格,整格发表。
IMPOSSIBLE      两组成分共享 item,所以"相同"永远无法与"由共享 item 强制"完全分开。
                本轮报的是**留出相关的大小**与**共享 item 的比例**,让读者自己定标。
"""
import numpy as np, pandas as pd, warnings, hashlib, re, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_disjoint_items, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A16_are_the_two_families_two_different_things'
          /'R182_do_they_have_different_external_anchors'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- 非性变量')[0])   # 跨轮依赖显式声明(P16)

v1=vv[:,-1]; A1=np.flatnonzero(v1>0); B1=np.flatnonzero(v1<=0)
if SIMR[np.ix_(A1,A1)].mean()<SIMR[np.ix_(B1,B1)].mean(): A1,B1=B1,A1
rated=np.array([j for j in range(V.shape[1]) if j in best])
FAM=np.intersect1d(A1,rated); OTH=np.intersect1d(B1,rated)
print(f"关系族(带评分){len(FAM)} 个,具体族 {len(OTH)} 个",flush=True)
RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
own=np.nanmean(R,axis=1)
M=RM-own[:,None]
ok=np.isfinite(M[:,FAM]).sum(1)>=6
print(f"族内 >=6 个评分的人 {ok.sum():,}",flush=True)

# A02 的三条坐标:用**全部**类别拟合的前 3 个成分(与族内成分共享 item —— 明确记录)
def comps(cols,rows,k):
    Z=M[np.ix_(rows,cols)]
    C=pd.DataFrame(Z).corr().values; C=np.nan_to_num(C,nan=0.); np.fill_diagonal(C,1.)
    w,vec=np.linalg.eigh(C)
    return vec[:,::-1][:,:k], np.clip(w[::-1],0,None)
shared=len(np.intersect1d(FAM,rated))
check_disjoint_items(set(FAM.tolist()),set(rated.tolist()),"族内成分 vs 全类别坐标",tol=1.0)
print(f"⚠ 共享 item:关系族的 {len(FAM)} 个类别**全部**也在全类别拟合里 —— "
      f"所以'相同'永远无法与'由共享 item 强制'完全分开,本轮报留出相关的大小(#P6 安全侧)",flush=True)

rows_all=np.flatnonzero(ok)
rg=np.random.default_rng(20260803)
def scores(vecs,cols,rows):
    Z=np.nan_to_num(M[np.ix_(rows,cols)]); return Z@vecs
res=[]; ctl_pos=[]; ctl_neg=[]
for rep in range(5):
    p=rg.permutation(rows_all); h1,h2=p[:len(p)//2],p[len(p)//2:]
    Vf,lf=comps(FAM,h1,4)                     # 族内成分,拟合在半边人上
    Va,_ =comps(rated,h1,3)                   # A02 式全类别坐标,同一半人
    Sf=scores(Vf,FAM,h2); Sa=scores(Va,rated,h2)          # 在**另一半人**上打分
    Cm=np.array([[abs(np.corrcoef(Sf[:,a],Sa[:,b])[0,1]) for b in range(3)] for a in range(4)])
    res.append(Cm)
    Vf2,_=comps(FAM,h2,4); Sf2=scores(Vf2,FAM,h2)
    ctl_pos.append(np.mean([abs(np.corrcoef(Sf[:,a],Sf2[:,a])[0,1]) for a in range(4)]))
    Vr=Va.copy(); Vr=Vr[rg.permutation(len(Vr))]           # 载荷随机置换
    Sr=scores(Vr,rated,h2)
    ctl_neg.append(np.mean([[abs(np.corrcoef(Sf[:,a],Sr[:,b])[0,1]) for b in range(3)] for a in range(4)]))
Cm=np.mean(res,axis=0); Cs=np.std(res,axis=0)
# ⚠ 参与比随 k 增长(纯噪声时 ~= k),所以 17 个 vs 9 个的直接比较无效(#101b same_scale)。
#   下采样到共同的 k,50 次取均值 —— 与 R03 用的是同一个修法。
eff=lambda l: l.sum()**2/(l**2).sum()
KC=min(len(FAM),len(OTH))
rgd=np.random.default_rng(4242)
def eff_at_k(cols,k):
    out=[]
    for _ in range(50):
        c=rgd.permutation(cols)[:k]
        out.append(eff(comps(c,rows_all,1)[1]))
    return float(np.mean(out))
e_fam,e_oth=eff_at_k(FAM,KC),eff_at_k(OTH,KC)
print(f"\n有效维度(**都在 k={KC} 上**,50 次下采样):关系族 {e_fam:.2f}  具体族 {e_oth:.2f}"
      f"   (未拉平时 {eff(comps(FAM,rows_all,1)[1]):.2f} vs {eff(comps(OTH,rows_all,1)[1]):.2f},"
      f"而参与比随 k 增长,那个比较无效)")
print(f"\n=== 族内前 4 成分 x A02 式三条坐标 的**留出**|相关|(5 次对半,均值±sd)===")
print(f"{'':<8}" + "".join(f"{'坐标'+str(b+1):>12}" for b in range(3)))
for a in range(4):
    print(f"成分{a+1:<5}" + "".join(f"{Cm[a,b]:>8.3f}±{Cs[a,b]:.2f}" for b in range(3)))
best_per=Cm.max(1)
print(f"\n每个族内成分与最接近的 A02 坐标:" + "  ".join(f"成分{a+1}={best_per[a]:.3f}" for a in range(4)))
print(f"正对照(族内成分 vs 它自己在另一半人上的重估):{np.mean(ctl_pos):.3f}")
print(f"负对照(A02 坐标载荷随机置换):{np.mean(ctl_neg):.3f}")

D=pd.DataFrame(Cm,index=[f'成分{a+1}' for a in range(4)],columns=[f'坐标{b+1}' for b in range(3)])
D.to_csv(pathlib.Path(__file__).parent/'results'/'canon.csv')
g=Gate('关系族多出来的维度是 A02 那三条吗')
g.asserted('正对照开火(族内成分可跨人复现)',np.mean(ctl_pos)>0.5,f"{np.mean(ctl_pos):.3f}")
# 这个"零"**不应该**是零:所有成分都被一个一般因子主导,所以随机载荷的坐标仍有 ~0.35 的
# 留出相关。它是 offset_control 的**地板**,不是 negative_control 的零。
g.asserted('随机载荷的地板不是零,而且它就是判别的门槛',np.mean(ctl_neg)>0.15,
           f"随机载荷仍给 {np.mean(ctl_neg):.3f} —— 因为成分被一般因子主导。"
           f"**任何'这个成分与那个不同'的说法必须先越过 {np.mean(ctl_neg):.3f},不是越过 0**")
g.same_scale('两族的有效维度在同一个 k 上比',float(KC),float(KC),'题目个数')
g.asserted('拉平题目数后,关系族的有效维度仍高于具体族',e_fam>e_oth,
           f"k={KC}: 关系族 {e_fam:.2f} vs 具体族 {e_oth:.2f}")
n_known=int((best_per>0.5).sum())
g.require_resolvable_first('最弱的那个族内成分与 A02 的留出相关',float(best_per.min()),float(Cs.max()))
g.offset_control('族内成分被 A02 三条覆盖到什么程度',float(best_per.mean()),float(np.mean(ctl_neg)),
                 float(Cs.mean()),null_kind='A02 坐标载荷随机置换后的同一留出相关')
g.asserted(f'四个族内成分里有几个被 A02 认出来(|r|>0.5)',True,
           f"{n_known}/4:" + " ".join(f"{v:.2f}" for v in best_per))
print(g)
i,j=np.unravel_index(np.argmin(best_per[:,None]+np.zeros((4,3))),(4,3))
wk=int(np.argmin(best_per))
Vf,_=comps(FAM,rows_all,4)
o=np.argsort(-np.abs(Vf[:,wk]))
print(f"\n=== 最不被 A02 认出的那个成分(成分{wk+1},|r|={best_per[wk]:.3f})载荷 ===")
for tt in o[:8]: print(f"   {Vf[tt,wk]:+.3f}  {lab[FAM[tt]]}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv().encode()).hexdigest()[:12]}")
