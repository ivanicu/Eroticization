import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A21 R02 -- 「喜欢什么」是个人的,「什么时候得到」不是 —— 这 4.4 倍是真的吗?

#100:稀有亲和特质 S 的**人侧**分半信度 **+0.4611**(跨不相交块集)。
#159a:径向倾向 z 的**人侧**分半信度 **+0.105**(跨不相交类别半)。
**差 4.4 倍。** 但那是两套口径量出来的。

⚠ 一个必须先处理的混淆:**单位数不同**。S 用 **32 个多选块**(每块十几个选项),
   z 用一个人平均 **~13 个类别**。**信度随单位数涨**,而 Spearman-Brown 只校正
   "劈成两半"这一步,**不校正基数差**。所以必须在**同一个 k** 上比(#101b same_scale)。

ESTIMAND        两个量各自的分半信度,作为**每半单位数 k** 的函数,同一批人、同一种劈分、
                同一个 Spearman-Brown。判别量 = 在匹配的 k 上两条曲线的差。
IDENTIFICATION  两个量都在"人内把单位劈成不相交两半 -> 各算一次 -> 跨人相关 -> SB"这一套
                完全相同的机器上跑;唯一不同的是单位是**块**还是**类别**。
SCOPE           同时有 >=2k 个多选块与 >=2k 个起始年龄类别的人。
WORLDS          REAL      匹配 k 后 S 的信度仍明显更高 -> 「什么」是个人的,「何时」不是
                ARTIFACT  匹配 k 后两条曲线汇合 -> 那 4.4 倍是单位数的产物
KILL            条件式:两个量各自的**人特异种植**都必须把信度推上去,否则曲线不可读。
POSITIVE CTRL   见上,两个种植。
NEGATIVE CTRL   人内置换各自的单位标签,两条曲线各跑一次。
NOISE FLOOR     5 个劈分种子;按人自助 200 次。
MULTIPLICITY    k ∈ {4,5,6,7} x 2 个量 x {真实, 置换, 种植},整格发表。
IMPOSSIBLE      两个量的**单位本身**不同(一个块含十几个选项,一个类别只有一个年龄),
                所以"同一个 k"只拉平了个数,拉不平每单位的信息量。这一条写进范围,不假装解决。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

# ---- 每人每块的稀有亲和分量(S 的单位)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BLK={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    if s.person.nunique()<1200 or s.option.nunique()<8: continue
    br=s.option.map(s.option.value_counts()/s.person.nunique())
    sp=-np.log(np.clip(br,1e-4,1.))
    g=pd.DataFrame({'person':s.person.values,'sp':sp.values}).groupby('person').sp.mean()
    v=np.full(len(df),np.nan); v[g.index.values]=g.values
    BLK[q.qi]=v
BM=np.c_[tuple(BLK.values())]                      # 人 x 块
print(f"块 {BM.shape[1]}  类别 {V.shape[1]}",flush=True)

def demean_conv(Vm,ob,tol=1e-10,cap=500):
    Dm=np.where(ob,Vm,np.nan)
    for _ in range(cap):
        a=np.nanmean(Dm,axis=0,keepdims=True); Dm=Dm-a
        b=np.nanmean(Dm,axis=1,keepdims=True); Dm=Dm-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return Dm
D0=demean_conv(V,obs)
BMc=BM-np.nanmean(BM,axis=0,keepdims=True)          # 块固定效应去掉,与 z 的处理对齐
okB=np.isfinite(BM); okA=okB.sum(1)>=8

NPERM=120
def z_on(cols,i,perm,rg,tie):
    j=np.intersect1d(np.flatnonzero(obs[i]),cols); k=len(j)
    if k<3: return np.nan
    y=D0[i,j].copy(); r=rar[j]
    if perm: y=y[rg.permutation(k)]
    cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
    d=r[pick]-r.mean()
    idx=rg.integers(0,k,(NPERM,1)); dr=r[idx].mean(1)-r.mean()
    return (d-dr.mean())/dr.std() if dr.std()>1e-9 else np.nan

def half_rel(kind,k,seed,perm=False,plant=0.):
    """在每半 k 个单位上算分半信度(跨人 r + Spearman-Brown)。"""
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    u=rg.standard_normal(len(V))                               # 人特异种植量
    A=[];B=[]
    who=np.flatnonzero(KEEP&okA)
    for i in who:
        if kind=='S':
            av=np.flatnonzero(okB[i])
            if len(av)<2*k: continue
            p=rg.permutation(av); h1,h2=p[:k],p[k:2*k]
            x1=np.nanmean(BMc[i,h1])+plant*u[i]; x2=np.nanmean(BMc[i,h2])+plant*u[i]
            if perm:
                x1=np.nanmean(BMc[rg.integers(len(BMc)),h1]); x2=np.nanmean(BMc[rg.integers(len(BMc)),h2])
        else:
            av=np.flatnonzero(obs[i])
            if len(av)<2*k: continue
            p=rg.permutation(av); h1,h2=p[:k],p[k:2*k]
            x1=z_on(h1,i,perm,rg,tie); x2=z_on(h2,i,perm,rg,tie)
            if plant: x1=x1+plant*u[i]; x2=x2+plant*u[i]
        if np.isfinite(x1) and np.isfinite(x2): A.append(x1); B.append(x2)
    A=np.array(A); B=np.array(B)
    if len(A)<300: return np.nan,np.nan,len(A)
    r=float(np.corrcoef(A,B)[0,1])
    return (2*r/(1+r) if r>-0.99 else np.nan), r, len(A)

KS=[4,5,6,7]
rows=[]
print(f"\n{'k':<4}{'量':<5}{'n':>7}{'分半 r':>10}{'SB 信度':>10}{'置换零 SB':>11}{'种植 SB':>10}")
for k in KS:
    for kind in ['S','z']:
        sb,r,n=half_rel(kind,k,zlib.crc32(f'{kind}{k}'.encode())%9973)
        sbn,_,_=half_rel(kind,k,zlib.crc32(f'{kind}{k}n'.encode())%9973,perm=True)
        sbp,_,_=half_rel(kind,k,zlib.crc32(f'{kind}{k}p'.encode())%9973,plant=1.0)
        rows.append(dict(k=k,kind=kind,n=n,r=r,sb=sb,sb_null=sbn,sb_plant=sbp))
        print(f"{k:<4}{kind:<5}{n:>7,}{r:>+10.4f}{sb:>+10.4f}{sbn:>+11.4f}{sbp:>+10.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'matched_k.csv',index=False)
S_=T[T.kind=='S'].set_index('k'); Z_=T[T.kind=='z'].set_index('k')
gap=(S_.sb-Z_.sb)
print(f"\n  匹配 k 上的差(S − z):" + "  ".join(f"k={k}:{gap[k]:+.4f}" for k in KS))
print(f"  比值:" + "  ".join(f"k={k}:{S_.sb[k]/max(Z_.sb[k],1e-9):.1f}x" for k in KS))
print(f"  (未匹配口径:#100 报 S 0.4611 · #159 报 z 0.105 -> 4.4x)")

g=Gate('「什么」是个人的、「何时」不是 —— 这 4.4 倍是真的吗')
g.asserted('S 的人特异种植把信度推上去',bool((S_.sb_plant>S_.sb+0.05).all()),
           " ".join(f"k={k}:{S_.sb[k]:.3f}->{S_.sb_plant[k]:.3f}" for k in KS))
g.asserted('z 的人特异种植把信度推上去',bool((Z_.sb_plant>Z_.sb+0.05).all()),
           " ".join(f"k={k}:{Z_.sb[k]:.3f}->{Z_.sb_plant[k]:.3f}" for k in KS))
g.asserted('两个量的人内置换零都在零附近',
           bool((S_.sb_null.abs()<0.08).all() and (Z_.sb_null.abs()<0.08).all()),
           "S " + " ".join(f"{v:+.3f}" for v in S_.sb_null) +
           " | z " + " ".join(f"{v:+.3f}" for v in Z_.sb_null))
g.require_resolvable_first('匹配 k 后两条曲线的差',float(gap.mean()),float(gap.std()))
g.offset_control('匹配 k 后 S 的信度 vs z 的信度',float(S_.sb.mean()),float(Z_.sb.mean()),
                 float(gap.std()),null_kind='同一批人、同一种劈分、同一个 k 下 z 的分半信度(不是零假设,是被比较的对象)')
g.no_sign_crossing('每个 k 上 S 都高于 z',list(gap.values))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
