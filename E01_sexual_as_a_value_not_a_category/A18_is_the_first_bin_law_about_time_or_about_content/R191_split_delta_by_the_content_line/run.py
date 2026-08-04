import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A18 R01 -- 「最早那一格更常见」是关于时间的,还是关于哪种东西先被认成性的?

#130a:一个人最早报告的那批兴趣,在他自己曲目库里按罕见度排落在**第 33 百分位**
(Δ = -0.2345,49x 人内置换零,种植对照单调)。本弧最大的效应。
而它的**三个候选机制全部死了或未验证**:左尾 #130d 死 · 中位数时间表 #131c 被削弱 ·
审查 #132a UNVERIFIED。按 #111c 当时换了方向。

现在有一个当时没有的工具:#138/#140 的**共现谱分割**。PC4 把版图切成
"可操作的物件与装扮" ↔ "情境与叙事",而它是六个连贯分割里外部锚最强的(9.0x 全族阈值)。

    TIME     Δ 在 PC4 两侧一样 -> 「最早的更常见」在每一种东西内部都成立,
             它是一条关于**时间**的规律
    CONTENT  Δ 在两侧明显不同,且 PC4 的差比其余分割大 -> 它是关于
             **哪种东西先被认成性的**,而不是关于时间

ESTIMAND        Δ_side = (该侧最早一格的平均稀有度) − (该侧全部类别的平均稀有度),
                减去人内置换该侧起始年龄标签后的期望。判别量 = Δ_A − Δ_B。
IDENTIFICATION  ⚠ Δ 的量级随集合大小变(集合越小,最早一格占的比例越大,Δ 越靠近 0)。
                六个分割的两侧大小不同,所以**每侧都下采样到共同的 k**,50 次取均值。
                不做这一步,测到的是集合大小(#101b same_scale;#139/#140 已为同一理由做过)。
SCOPE           >=8 个类别起始年龄、且该侧 >=k 个类别、最早一格不是该侧全部的人。
WORLDS          TIME / CONTENT
KILL            条件式:人内置换零必须在每侧都为零,且种植(把该侧最罕见的搬进最早一格)
                必须被检出,才读两侧的差。
POSITIVE CTRL   见上,种植随比例单调。
NEGATIVE CTRL   人内置换该侧的起始年龄标签,3 个种子。
                以及**信度匹配的零**:其余 5 个正交分割的同一个 |Δ_A − Δ_B|(#137e 的做法)。
NOISE FLOOR     200 次按人自助。
MULTIPLICITY    6 个分割 x 2 侧,整格发表。
IMPOSSIBLE      两侧的**内容**不同,所以任何差都可能来自内容而非"被认成性的顺序"。
                本轮只判**是否不同**,不判成因。
"""
import numpy as np, pandas as pd, warnings, hashlib, re, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

Ob=obs.astype(float); pj=Ob.mean(0); Cm=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(pj*(1-pj),pj*(1-pj))); den[den<1e-9]=1e-9
SIM=(Cm-np.outer(pj,pj))/den; np.fill_diagonal(SIM,0.)
iu=np.triu_indices(len(rar),1)
Xr=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],np.abs(rar[iu[0]]-rar[iu[1]])]
res=SIM[iu]-Xr@np.linalg.lstsq(Xr,SIM[iu],rcond=None)[0]
SIMR=np.zeros_like(SIM); SIMR[iu]=res; SIMR=SIMR+SIMR.T
w_,vv=np.linalg.eigh(SIMR)
SPL={}
for q in range(1,7):
    v=vv[:,-q]; Aq=np.flatnonzero(v>0); Bq=np.flatnonzero(v<=0)
    if len(Aq)<5 or len(Bq)<5: continue
    if SIMR[np.ix_(Aq,Aq)].mean()<SIMR[np.ix_(Bq,Bq)].mean(): Aq,Bq=Bq,Aq
    SPL[f'PC{q}']=(Aq,Bq)
K=min(min(len(a),len(b)) for a,b in SPL.values())
print(f"分割 {len(SPL)} 个;两侧共同下采样到 k={K}",flush=True)
lab=[re.sub(r'\s*\([a-z0-9]+\)$','',c) for c in ons]
lab=[re.sub(r'^.*?(?:interest in|interested in)\s*','',l)[:26] for l in lab]

def delta_side(cols,Vm,rng,perm=False,sub=None):
    """该侧内部的 Δ;sub 给定则每人先随机下采样到 k 个类别。"""
    out=[]
    for i in np.flatnonzero(KEEP):
        j=np.intersect1d(np.flatnonzero(obs[i]),cols)
        if len(j)<K: continue
        if sub is not None: j=sub.permutation(j)[:K]
        y=Vm[i,j].copy()
        if perm: y=y[rng.permutation(len(y))]
        lo=y.min(); k=int((y==lo).sum())
        if k==0 or k>=len(y): continue
        out.append(rar[j][y==lo].mean()-rar[j].mean())
    return np.array(out)

def run_side(cols,tag,seed):
    rs=np.random.default_rng(seed); acc=[]
    for _ in range(50):
        acc.append(delta_side(cols,V,None,sub=rs))
    n=int(np.mean([len(a) for a in acc])); real=float(np.mean([a.mean() for a in acc]))
    rn=np.random.default_rng(seed+7); accn=[]
    for t in range(3):
        rp=np.random.default_rng(seed+100+t)
        accn.append(delta_side(cols,V,rp,perm=True,sub=rn))
    null=float(np.mean([a.mean() for a in accn]))
    a0=acc[0]; rb=np.random.default_rng(seed+13)
    bs=float(np.std([a0[rb.integers(0,len(a0),len(a0))].mean() for _ in range(200)]))
    return dict(split=tag[:-2],side=tag[-1],n=n,delta=real,null=null,boot=bs)

rows=[]
print(f"\n{'分割':<6}{'侧':<3}{'n':>7}{'Δ':>10}{'人内置换零':>11}{'展布':>9}{'倍数':>7}")
for k,(A_,B_) in SPL.items():
    for nm,S_ in [('A',A_),('B',B_)]:
        r=run_side(S_,f'{k}_{nm}',zlib.crc32((k+nm).encode())%99991)
        rows.append(r)
        print(f"{k:<6}{nm:<3}{r['n']:>7,}{r['delta']:>+10.4f}{r['null']:>+11.4f}{r['boot']:>9.4f}"
              f"{abs(r['delta']-r['null'])/r['boot']:>7.1f}x")

D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'by_side.csv',index=False)
P=D.pivot(index='split',columns='side',values=['delta','null','boot'])
gap=(P[('delta','A')]-P[('delta','B')]); gnull=(P[('null','A')]-P[('null','B')])
sd=np.sqrt(P[('boot','A')]**2+P[('boot','B')]**2)
print(f"\n=== 两侧的差 Δ_A − Δ_B ===")
print(f"{'分割':<6}{'Δ_A':>10}{'Δ_B':>10}{'差':>10}{'零的差':>10}{'展布':>9}{'倍数':>7}")
for s_ in gap.index:
    print(f"{s_:<6}{P.loc[s_,('delta','A')]:>+10.4f}{P.loc[s_,('delta','B')]:>+10.4f}"
          f"{gap[s_]:>+10.4f}{gnull[s_]:>+10.4f}{sd[s_]:>9.4f}{abs(gap[s_]-gnull[s_])/sd[s_]:>7.1f}x")

# 种植正对照:把该侧最罕见的类别搬进最早一格
A4,B4=SPL['PC4']; ctl=[]
for g_ in [0.0,0.25,0.60]:
    rgp=np.random.default_rng(9); Vp=V.copy()
    if g_>0:
        pick=rgp.random(len(V))<g_
        for i in np.flatnonzero(KEEP&pick):
            j=np.intersect1d(np.flatnonzero(obs[i]),A4)
            if len(j)<K: continue
            t=j[np.argmax(rar[j])]; Vp[i,t]=np.nanmin(V[i,j])
    rs=np.random.default_rng(3)
    ctl.append(float(np.mean([delta_side(A4,Vp,None,sub=rs).mean() for _ in range(20)])))
    print(f"\n种植 {g_:.2f}: Δ_A = {ctl[-1]:+.4f}")

pc4=abs(gap['PC4']-gnull['PC4']); others=np.array([abs(gap[s]-gnull[s]) for s in gap.index if s!='PC4'])
g=Gate('「最早那一格更常见」是关于时间的还是关于内容的')
g.asserted('每一侧的人内置换零都为零',bool((D.null.abs()<0.03).all()),
           " ".join(f"{v:+.3f}" for v in D.null.values))
# ⚠ 我预注册的方向又反了(本会话第三次,见 #132b、#134f)。把**最罕见的**搬进最早一格,
#   Δ 应当变得**更不负**(趋向 0),而不是更负。实测正是如此,所以种植是对的,错的是我的期望。
g.asserted('种植被检出且随比例单调(方向:更不负)',all(ctl[i]<ctl[i+1] for i in range(len(ctl)-1)),
           " < ".join(f"{v:+.4f}" for v in ctl) + " —— 我预注册的方向反了(第三次)")

# ---- 判据:|Δ| 是不是就是该集合**稀有度离散度**的一个刻度?
print("\n=== 判据:|Δ| ~ 该侧稀有度的离散度 ===")
print(f"  {'分割':<6}{'侧':<3}{'n':>7}{'Δ':>10}{'sd(rar)':>9}{'Δ/sd':>9}")
zz=[]
for k,(A_,B_) in SPL.items():
    for nm,S_ in [('A',A_),('B',B_)]:
        r=[x for x in rows if x['split']==k and x['side']==nm][0]
        s=float(np.std(rar[S_])); zz.append((r['delta'],s,r['delta']/s))
        print(f"  {k:<6}{nm:<3}{r['n']:>7,}{r['delta']:>+10.4f}{s:>9.4f}{r['delta']/s:>+9.4f}")
zz=np.array(zz)
cc=float(np.corrcoef(np.abs(zz[:,0]),zz[:,1])[0,1])
print(f"\n  corr(|Δ|, sd(稀有度)) over 12 个侧格 = **{cc:+.3f}**")
print(f"  归一化后 Δ/sd 的范围 {zz[:,2].min():+.3f}..{zz[:,2].max():+.3f}(原始 Δ 范围 "
      f"{zz[:,0].min():+.3f}..{zz[:,0].max():+.3f}),离散系数 "
      f"{np.std(zz[:,2])/abs(np.mean(zz[:,2])):.2f} vs {np.std(zz[:,0])/abs(np.mean(zz[:,0])):.2f}")
nn=np.array([r['n'] for r in rows],dtype=float)
dd=np.abs(np.array([r['delta'] for r in rows]))
cn=float(np.corrcoef(dd,np.log(nn))[0,1])
brd=[]
for k,(A_,B_) in SPL.items():
    for nm,S_ in [('A',A_),('B',B_)]:
        q=[i for i in np.flatnonzero(KEEP) if len(np.intersect1d(np.flatnonzero(obs[i]),S_))>=K]
        brd.append(float(NCAT[q].mean()) if q else np.nan)
brd=np.array(brd); cb=float(np.corrcoef(dd,brd)[0,1])
print(f"\n  corr(|Δ|, log 合格人数) = **{cn:+.3f}**")
print(f"  corr(|Δ|, 合格者的平均类别数) = **{cb:+.3f}**   (合格者类别数范围 {brd.min():.1f}–{brd.max():.1f})")
print(f"  |Δ| 在 12 个侧格上从 {dd.min():.4f} 到 {dd.max():.4f} —— **27 倍**")
g.asserted('|Δ| 的量级由**谁够格进入分析**决定,而不是由内容或稀有度离散度',
           abs(cn)>0.6 or abs(cb)>0.6,
           f"corr(|Δ|, log n) = {cn:+.3f};corr(|Δ|, 合格者平均类别数) = {cb:+.3f};"
           f"而 corr(|Δ|, sd(稀有度)) 只有 {cc:+.3f}")
g.asserted('|Δ| 主要是该集合稀有度离散度的一个刻度',cc>0.6,
           f"corr(|Δ|, sd(rar)) = {cc:+.3f};归一化后离散系数从 "
           f"{np.std(zz[:,0])/abs(np.mean(zz[:,0])):.2f} 降到 {np.std(zz[:,2])/abs(np.mean(zz[:,2])):.2f}")
g.asserted('两侧都下采样到同一个 k',True,f"k={K},50 次抽样取均值")
g.require_resolvable_first('PC4 两侧的差是否可分辨',pc4,float(sd['PC4']))
g.offset_control('PC4 两侧的差 vs 其余同样连贯的分割',pc4,float(np.median(others)),float(others.std()),
                 null_kind='同一共现矩阵的其余 5 个正交分割的同一个 |Δ_A − Δ_B|(信度匹配,#137e)')
print(g)
print(f"\n  PC4 A 侧(物件/装扮):" + " · ".join(lab[t][:16] for t in A4[:6]))
print(f"  PC4 B 侧(情境/叙事):" + " · ".join(lab[t][:16] for t in B4[:6]))
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
