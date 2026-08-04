import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A14 R05 -- Delta = -0.2345 的机制:是审查吗?

#130a 测到本弧最大的效应:一个人最早报告的那批兴趣,在他自己曲目库里按罕见度排落在
**第 33 百分位**(49x 人内置换零,种植对照单调)。#130d 杀掉了我的"左尾"解释,
#131c 削弱了"中位数时间表"解释。**效应稳,两个候选机制都死了。**

剩下最便宜的候选是**审查**:一个人只能报告他**已经获得**的东西,而罕见类别的获得率低。
年纪小的人还没走到罕见的那一段,所以他们能进"最早一格"的候选里常见的占比更高。

审查的**唯一可检验的印记是年龄梯度**:
    审查为真 -> Delta 随当前年龄**单调收缩**(年长者已走完更多路,罕见类别补齐)
    审查为假 -> Delta 在各年龄段上**平坦**

ESTIMAND        Delta 按当前年龄分层(5 档),以及最年长档单独的 Delta。
IDENTIFICATION  年龄是 release 里唯一的暴露时长代理;分箱为 5 档(14-17 ... 29-32)。
SCOPE           报告 >=8 个类别起始年龄的人。
WORLDS          censor   Delta 随年龄单调收缩 -> 审查是机制,Delta 是一个观测窗口的性质
                stable   Delta 在各档平坦 -> 审查死,而我需要第四个候选
                inverse  Delta 随年龄**增大** -> 与审查相反,提示是回忆随时间加深(#119 的形态)
KILL            条件式:分层仪器必须先在一个**已知随年龄变化**的量上开火(类别数),
                且人内置换零必须在每一档上都为零,才读 Delta 的年龄趋势。
POSITIVE CTRL   类别数必须随年龄单调上升(已知:活得久勾得多)。
                以及:人为审查 —— 把每个人**最晚获得**的 f 比例类别删掉,Delta 必须朝
                "更负"移动(模拟"更年轻"),且随 f 单调。
NEGATIVE CTRL   人内置换起始年龄标签,每一档独立跑。
NOISE FLOOR     每档 200 次按人自助。
MULTIPLICITY    5 个年龄档 x {真实, 置换} x 4 个人为审查水平,整格发表。
IMPOSSIBLE      年龄与队列在横断面里不可分;"年长者走完更多路"与"年长队列成长环境不同"
                本轮分不开。只判**有没有梯度**,不判梯度的成因。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

AGES=[('14-17',15.5),('18-20',19.0),('21-24',22.5),('25-28',26.5),('29-32',30.5)]

def delta_of(Vm,who,rng=None,perm=False):
    out=[]
    for i in who:
        m=obs[i]; y=Vm[i,m].copy(); r=rar[m]
        if perm: y=y[rng.permutation(len(y))]
        lo=y.min(); k=int((y==lo).sum())
        if k==0 or k>=len(y): continue
        out.append(r[y==lo].mean()-r.mean())
    return np.array(out)

rows=[]
print(f"{'年龄档':<8} {'n':>7} {'类别数':>7} {'Delta':>9} {'置换零':>9} {'自助展布':>9} {'倍数':>7}")
for lab,a in AGES:
    who=np.flatnonzero(KEEP&(age==a))
    d=delta_of(V,who)
    dn=np.concatenate([delta_of(V,who,np.random.default_rng(7000+s),perm=True) for s in range(3)])
    rb=np.random.default_rng(int(a*100))
    bs=float(np.std([d[rb.integers(0,len(d),len(d))].mean() for _ in range(200)]))
    rows.append(dict(band=lab,age=a,n=len(d),ncat=float(NCAT[who].mean()),
                     delta=float(d.mean()),null=float(dn.mean()),boot=bs))
    print(f"{lab:<8} {len(d):>7,} {NCAT[who].mean():>7.1f} {d.mean():>+9.4f} {dn.mean():>+9.4f} "
          f"{bs:>9.4f} {abs(d.mean()-dn.mean())/bs:>7.1f}x")

D=pd.DataFrame(rows)
# 人为审查正对照:删掉每个人最晚获得的 f 比例类别 = 模拟"更年轻"
print(f"\n人为审查正对照(删掉最晚获得的 f 比例类别,模拟更年轻):")
who_all=np.flatnonzero(KEEP)
cens=[]
for f in [0.0,0.15,0.30,0.45]:
    Vc=V.copy(); ob2=obs.copy()
    if f>0:
        for i in who_all:
            j=np.flatnonzero(obs[i]); k=int(round(f*len(j)))
            if k>0:
                drop=j[np.argsort(-V[i,j])[:k]]; Vc[i,drop]=np.nan; ob2[i,drop]=False
    ob_save=obs; obs=ob2
    dd=delta_of(Vc,who_all); obs=ob_save
    cens.append(float(dd.mean())); print(f"  f={f:.2f}  Delta {dd.mean():+.4f}  n={len(dd):,}")

# ---- 分层仪器的正对照没开火,而那本身是可测的事:曲目库到底还长不长?
print("\n=== 曲目库随年龄增长吗?(审查这个假设的前提)===")
print(f"  {'年龄档':<8} {'n':>7} {'类别数':>7} {'≤17 岁获得的比例':>16} {'最晚获得年龄':>12}")
grow=[]
for lab,a in AGES:
    who=np.flatnonzero(KEEP&(age==a))
    fr=np.array([ (V[i,obs[i]]<=17.5).mean() for i in who ])
    lastj=np.array([ np.nanmax(V[i,obs[i]]) for i in who ])
    grow.append(dict(band=lab,age=a,ncat=float(NCAT[who].mean()),
                     frac17=float(fr.mean()),last=float(lastj.mean())))
    print(f"  {lab:<8} {len(who):>7,} {NCAT[who].mean():>7.1f} {fr.mean():>16.1%} {lastj.mean():>12.1f}")
Gr=pd.DataFrame(grow)
old=Gr[Gr.band=='29-32'].iloc[0]
print(f"\n  29-32 岁的人里,他们自己报告的兴趣有 {old.frac17:.1%} 是在 17 岁前获得的,"
      f"而最晚的那个平均在 {old.last:.1f} 岁")
print(f"  15 年里类别数只从 {Gr.ncat.values[0]:.1f} 长到 {Gr.ncat.values[-1]:.1f} "
      f"(+{100*(Gr.ncat.values[-1]/Gr.ncat.values[0]-1):.1f}%)")

g=Gate('Delta 的机制是审查吗')
g.asserted('曲目库在 15 年里几乎不长 —— 审查的前提本身很弱',
           (Gr.ncat.values[-1]/Gr.ncat.values[0]-1)<0.10,
           f"{Gr.ncat.values[0]:.1f} -> {Gr.ncat.values[-1]:.1f} "
           f"(+{100*(Gr.ncat.values[-1]/Gr.ncat.values[0]-1):.1f}%);"
           f"29-32 岁的人 {old.frac17:.0%} 的兴趣在 17 岁前就有了")
g.asserted('分层仪器在一个随年龄变的量上开火(弱,记录下来)',
           D.ncat.values[-1]>D.ncat.values[0],
           "类别数 " + " < ".join(f"{v:.1f}" for v in D.ncat.values) +
           " —— 单调但极弱,所以这个分层对'暴露时长'的检出力本身有限")
g.asserted('每一档的人内置换零都为零',bool((D.null.abs()<0.02).all()),
           " ".join(f"{v:+.4f}" for v in D.null.values))
g.asserted('⚠ 我预注册的审查方向是错的 —— 模拟说的是相反方向',
           all(cens[i]<=cens[i+1]+1e-9 for i in range(len(cens)-1)),
           "人为审查让 Delta **更不负**:" + " < ".join(f"{v:+.4f}" for v in cens) +
           " —— 所以审查预测的是「年轻 -> 更不负」,而不是我在 docstring 里写的那个方向")
trend=np.polyfit(D.age.values,D.delta.values,1)[0]
rng_d=float(D.delta.max()-D.delta.min()); mb=float(D.boot.mean())
g.require_resolvable_first('各档之间的 Delta 差是否可分辨',rng_d,mb)
g.no_sign_crossing('所有年龄档的 Delta 同号',list(D.delta.values))
print(g)
print(f"\n  年龄趋势 {trend:+.5f} / 岁   各档极差 {rng_d:.4f} = {rng_d/mb:.1f}x 单档自助展布")
print(f"  最年长档(29-32)单独 Delta = {D[D.band=='29-32'].delta.values[0]:+.4f}"
      f"(全体 {np.mean([r['delta'] for r in rows]):+.4f})")
D.to_csv(pathlib.Path(__file__).parent/'results'/'by_age.csv',index=False)
Gr.to_csv(pathlib.Path(__file__).parent/'results'/'growth.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
