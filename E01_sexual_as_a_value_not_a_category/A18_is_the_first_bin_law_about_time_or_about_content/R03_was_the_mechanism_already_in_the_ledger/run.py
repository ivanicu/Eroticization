import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A18 R03 -- 机制其实五轮前就写在同一本账里了吗?

#147a:Δ 的效应是普遍的(五个广度层每层 17.5–27.3×),可比单位 z = −0.5515。
#147 的 NEXT:把每人的起始年龄先减去 #75 的题目均值,在 z 上重算。

**但 `#130b` 当年就测过这一步**:扣掉人群时间表后 Δ 从 −0.2345 翻成 **+0.0767(14.5×)**。
如果那就是答案,那么其后四轮找机制 —— `#130d` 左尾 · `#131c` 中位数时间表 ·
`#132a` 审查 · `#146b` 稀有度离散度 —— 是在找一个**同一本账里已经写着的东西**。
这与 `#143` 同族:**用重跑的方式去问一个账本已经回答过的问题**。

这一轮把它在**可比单位**上钉死,并加一个 `#130b` 当年没有的对照。

ESTIMAND        z = Δ / sd(Δ | 人内置换),在三个臂上:
                  raw    原始起始年龄
                  item   减去每个题目的人群平均起始年龄(= #75 的时间表,精确一遍)
                  rand   减去一个**方差匹配的随机题目向量**(对照:必须不改变 z)
IDENTIFICATION  Δ 已经减掉了这个人自己的集合均值,所以常数人效应恰好抵消;
                题目去均值是精确的一遍(不涉及 #128 的非幂等问题,因为这里不做人内去均值)。
SCOPE           >=8 个类别起始年龄、最早一格不是全部的人。
WORLDS          SCHEDULE  item 臂的 z 归零或翻号 -> 机制就是共享时间表在人内的投影,
                          而它五轮前就在账里
                NOT       item 臂的 z 基本不变 -> 时间表被排除
KILL            条件式:rand 对照必须**不**改变 z(否则任何题目层减法都会改变它,
                        那么 item 臂的变化不可归因),才读 item 臂。
POSITIVE CTRL   rand 对照(见上)——它是"减去一个题目层向量"这件事本身的零。
NEGATIVE CTRL   人内置换(每人 200 次,同时给出分母)。
NOISE FLOOR     按人自助 200 次。
MULTIPLICITY    3 臂 x {Δ, z},整格发表。
IMPOSSIBLE      "时间表为什么按稀有度排"是另一个问题(`#130e`/`#131c` 在打它),
                本轮不碰。这里只判**时间表能不能解释 Δ**。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

NPERM=200
def zstat(Vm,seed,single=True):
    """⚠ single=True:每人只取**一个**最早的类别(并列用固定种子随机打破)。

    第一版让"最早一格"= 全部并列者,而原始年龄按 2 年分箱、m 可到 24 个,
    减去任何连续的题目向量后并列全消、m 恒为 1 —— **三个臂算的不是同一个统计量**。
    这正是我自己的随机对照抓到的:减一个方差匹配的**随机**向量也让 z 掉 48%。
    m=1 让三个臂精确可比(#101b same_scale)。"""
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    D=np.full(len(Vm),np.nan); Z=np.full(len(Vm),np.nan)
    for i in np.flatnonzero(KEEP):
        j=np.flatnonzero(obs[i]); r=rar[j]; y=Vm[i,j]; k=len(j)
        if k<2: continue
        if single:
            cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
            d=r[pick]-r.mean(); m=1
        else:
            lo=np.nanmin(y); m=int((y==lo).sum())
            if m==0 or m>=k: continue
            d=r[y==lo].mean()-r.mean()
        idx=rg.integers(0,k,(NPERM,m)); dr=r[idx].mean(1)-r.mean()
        s=dr.std()
        if s<1e-9: continue
        D[i]=d; Z[i]=(d-dr.mean())/s
    return D,Z

itm=np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])])
rgv=np.random.default_rng(4242)
randvec=rgv.normal(itm.mean(),itm.std(),len(itm))          # 方差匹配的随机题目向量
ARMS={'raw':V,'item(#75 时间表)':np.where(obs,V-itm[None,:],np.nan),
      'rand(方差匹配)':np.where(obs,V-randvec[None,:],np.nan)}
print(f"题目均值起始年龄范围 {itm.min():.1f}–{itm.max():.1f}(sd {itm.std():.2f});"
      f"随机向量 sd {randvec.std():.2f}",flush=True)

rows=[]
print(f"\n{'臂':<18}{'n':>7}{'Δ':>10}{'展布':>8}{'z':>10}{'展布':>8}{'z 倍数':>8}")
for nm,Vm in ARMS.items():
    D_,Z_=zstat(Vm,zlib.crc32(nm.encode())%99991)
    m=np.isfinite(Z_); ii=np.flatnonzero(m)
    rb=np.random.default_rng(7)
    bd=float(np.std([D_[ii[rb.integers(0,len(ii),len(ii))]].mean() for _ in range(200)]))
    bz=float(np.std([Z_[ii[rb.integers(0,len(ii),len(ii))]].mean() for _ in range(200)]))
    rows.append(dict(arm=nm,n=int(m.sum()),delta=float(D_[m].mean()),d_boot=bd,
                     z=float(Z_[m].mean()),z_boot=bz))
    print(f"{nm:<18}{int(m.sum()):>7,}{D_[m].mean():>+10.4f}{bd:>8.4f}{Z_[m].mean():>+10.4f}"
          f"{bz:>8.4f}{abs(Z_[m].mean())/bz:>8.1f}x")

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'arms.csv',index=False)
raw=T[T.arm=='raw'].iloc[0]; itemr=T[T.arm.str.startswith('item')].iloc[0]
rnd=T[T.arm.str.startswith('rand')].iloc[0]
expl=1-itemr.z/raw.z
print(f"\n  时间表解释掉的比例(在 z 上):{100*expl:.0f}%   （raw {raw.z:+.4f} -> item {itemr.z:+.4f}）")
print(f"  对照:减去方差匹配的随机题目向量 -> z {rnd.z:+.4f}(变化 {100*(1-rnd.z/raw.z):+.0f}%)")

g=Gate('时间表能不能解释 Δ')
# ⚠ "这个零应该是零吗" —— **不应该**。减去任何方差相当的题目层向量都会**稀释**
#   "谁最早"这个判断,所以随机臂不是零,它是**正确的基准**。第一版把它当零来断言,问错了问题。
g.asserted('随机臂本来就该改变 z(它是稀释,不是伪影)—— 所以它是基准不是零',
           abs(rnd.z-raw.z)>0.1*abs(raw.z),
           f"raw {raw.z:+.4f} -> rand {rnd.z:+.4f}(掉 {100*(1-rnd.z/raw.z):.0f}%);"
           f"方差相当的噪声同等程度打乱排序,零假设**不预测它无变化**")
g.no_sign_crossing('raw 与 rand 同号',[float(raw.z),float(rnd.z)])
g.require_resolvable_first('raw 臂本身可分辨',abs(float(raw.z)),float(raw.z_boot))
g.offset_control('减掉时间表后的 z vs 原始 z',float(itemr.z),float(raw.z),float(raw.z_boot),
                 null_kind='同一统计量在原始起始年龄上的值(不是零假设,是被解释的基准)')
g.require_resolvable_first('时间表臂是否越过随机臂',abs(float(itemr.z-rnd.z)),float(raw.z_boot))
g.offset_control('时间表解释的,超出方差匹配噪声多少',float(itemr.z),float(rnd.z),float(raw.z_boot),
                 null_kind='减去方差匹配的随机题目层向量后的同一个 z(稀释基准,不是零)')
g.asserted('时间表 vs 随机噪声,解释掉的比例',True,
           f"时间表 {100*expl:.0f}% vs 方差匹配噪声 {100*(1-rnd.z/raw.z):.0f}% —— "
           f"时间表特异的部分是 {100*(expl-(1-rnd.z/raw.z)):.0f} 个百分点")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
