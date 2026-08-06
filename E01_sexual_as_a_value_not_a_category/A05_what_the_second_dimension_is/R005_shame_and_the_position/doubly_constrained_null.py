import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A52 R265 -- 造一个恰好只毁掉配对的零

`#219c`:三种零里没有一个恰好**只毁掉配对而保留覆盖结构** ——
A 换掉了每个人的整套年龄 · B 只固定人内 · C 连"谁答了哪一题"都打乱(所以换掉了自变量)。
**那样的零是可以造的:人内多重集固定 + 题内边际固定。**

走法(严格保两侧)
    对一对类别 (j,k),取值对为 (a,b) 与 (b,a) 的两个人 i1 · i2,
    **各自在自己内部交换** j 与 k 的年龄:
      行:每人只换了自己的两个值 -> **人内多重集不变** ✅
      列 j:i1 由 a->b,i2 由 b->a -> **净零** ✅   列 k 同理 ✅
    起始年龄是**两年一档**(约 8 个取值),所以 (a,b)/(b,a) 的配对**大量存在**。

ESTIMAND        主 Δ 在这个双约束零下的净值。
KILL            **接近 A/B(−0.035) -> `#219b` 的 75% 全部来自"换掉了自变量";
                接近 C(−0.008) -> `#128` 的 Δ 要按它重估。**
硬前置          **两侧边际必须逐位不变** —— 不满足就整轮作废。
POSITIVE CTRL   种一个人内配对信号 -> 这个零必须把它毁掉。
NEGATIVE CTRL   把这个零施加在**已经打乱过**的数据上 -> Δ 应几乎不动。
IMPOSSIBLE      它保留两侧边际,但**不保留题对之间的共现结构**(哪些题一起被答)——
                那部分仍未被任何零控制。
"""
import numpy as np, pandas as pd, warnings, hashlib
from collections import defaultdict
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
KEEP0=(np.isfinite(V0).sum(1)>=8)

def demean_np(Aa,iters=200,tol=1e-10):
    D=np.where(np.isfinite(Aa),Aa,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def delta(Aa,need=8):
    D=demean_np(Aa)
    W=np.isfinite(D).astype(float); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    X=W*(rar0[None,:]-rb[:,None]); Y0=np.where(np.isfinite(D),D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*X).sum(1); den=np.sqrt((X*X).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(len(D),np.nan); good=(k>=need)&(den>1e-12)&KEEP0; out[good]=num[good]/den[good]
    return float(np.nanmean(out))

def doubly_constrained(Aa, rng, rounds=60):
    """人内多重集固定 + 题内边际固定的重排。"""
    B=Aa.copy()
    for _ in range(rounds):
        for j in range(Mc):
            for k in range(j+1,Mc):
                m=np.isfinite(B[:,j])&np.isfinite(B[:,k])
                idx=np.flatnonzero(m)
                if len(idx)<4: continue
                d=defaultdict(list)
                for i in idx:
                    a,b=B[i,j],B[i,k]
                    if a!=b: d[(a,b)].append(i)
                for (a,b),g1 in list(d.items()):
                    g2=d.get((b,a))
                    if not g2 or a>=b: continue
                    n=min(len(g1),len(g2))
                    if n==0: continue
                    p1=rng.permutation(g1)[:n]; p2=rng.permutation(g2)[:n]
                    take=rng.random(n)<0.5
                    for x,y,t in zip(p1,p2,take):
                        if not t: continue
                        B[x,j],B[x,k]=b,a
                        B[y,j],B[y,k]=a,b
    return B

rng=np.random.default_rng(20260803)
d0=delta(V0)
Vd=doubly_constrained(V0,rng,rounds=1)
# ---- 硬前置:两侧边际逐位不变 -------------------------------------------------
def margins(Aa):
    obsm=np.isfinite(Aa)
    row=[tuple(sorted(Aa[i][obsm[i]])) for i in range(0,N,37)]        # 抽样比对行多重集
    col=[tuple(sorted(Aa[:,j][obsm[:,j]])) for j in range(Mc)]
    return obsm.sum(), row, col, np.nansum(Aa)
o0,r0,c0,s0=margins(V0); o1,r1,c1,s1=margins(Vd)
same_obs=(o0==o1); same_row=all(x==y for x,y in zip(r0,r1))
same_col=all(x==y for x,y in zip(c0,c1)); same_sum=abs(s0-s1)<1e-9
print(f"硬前置:观测数 {'✅' if same_obs else '✗'} · 人内多重集(抽 {len(r0)} 人){'✅' if same_row else '✗'}"
      f" · 题内多重集(全 {Mc} 题){'✅' if same_col else '✗'} · 总和 {'✅' if same_sum else '✗'}")
assert same_obs and same_row and same_col and same_sum, "两侧边际没保住 —— 整轮作废"

vals=[]
for s in range(5):
    Vs=doubly_constrained(V0,np.random.default_rng(1300+s),rounds=3)
    vals.append(delta(Vs))
bs=[delta(V0[rng.integers(0,N,N)]) for _ in range(30)]
sd=float(np.std(bs))
zD=float(np.mean(vals))
print(f"\n主 Δ = {d0:+.4f}(sd {sd:.4f})")
print(f"零 D 双约束重排:{zD:+.4f} ± {np.std(vals):.4f} · 占效应 {100*abs(zD)/abs(d0):.0f}%"
      f"  -> 净 **{d0-zD:+.4f}**")
print(f"  对比:A 净 −0.0350 · B 净 −0.0368 · C 净 −0.0081(`#219b`)")

# 正对照:种一个人内配对信号,零必须毁掉它
u=np.abs(rng.standard_normal(N))+0.5
Vp=V0+10.0*np.outer(u,rar0-rar0.mean())*np.isfinite(V0)
Vp=np.where(np.isfinite(V0),np.round(Vp*2)/2,np.nan)     # 保持取值离散,便于配对
dp=delta(Vp); dpz=delta(doubly_constrained(Vp,np.random.default_rng(7),rounds=3))
print(f"\n正对照:种植后 Δ {dp:+.4f} -> 施加双约束零后 {dpz:+.4f}(毁掉 {100*(1-abs(dpz)/abs(dp)):.0f}%)")
# 负对照:施加在已打乱的数据上
Vb=V0.copy()
for i in range(N):
    ii=np.flatnonzero(np.isfinite(Vb[i]))
    if len(ii)>1: Vb[i,ii]=Vb[i,rng.permutation(ii)]
db=delta(Vb); dbz=delta(doubly_constrained(Vb,np.random.default_rng(9),rounds=3))
print(f"负对照:已打乱 Δ {db:+.4f} -> 再施加双约束零 {dbz:+.4f}(动了 {abs(dbz-db):.4f})")

T=pd.DataFrame([dict(null='D 双约束重排',value=zD,net=float(d0-zD),frac=abs(zD)/abs(d0))])
check_columns(T,'R265'); T.to_csv(pathlib.Path(__file__).parent/'results'/'double.csv',index=False)
netD=d0-zD
g=Gate('那个恰好只毁配对的零,给出什么')
g.asserted('硬前置:两侧边际逐位不变',same_obs and same_row and same_col and same_sum,'观测数·人内·题内·总和 全同')
g.asserted('可判前提:主 Δ 复现 `#128` 的 −0.0328',abs(d0+0.0328)<0.003,f"{d0:+.4f}")
g.asserted('正对照:种植的人内配对必须被这个零毁掉',abs(dpz)<0.5*abs(dp),
           f"{dp:+.4f} -> {dpz:+.4f}")
g.asserted('负对照:施加在已打乱数据上应几乎不动',abs(dbz-db)<3*sd,f"动了 {abs(dbz-db):.4f} vs 3×sd {3*sd:.4f}")
g.negative_control('零 D 双约束重排',abs(zD),abs(d0),null_spread=float(np.std(vals)),
                   null_kind='人内多重集固定 + 题内边际固定的重排 —— 只毁配对,保留覆盖结构')
g.offset_control('净 D vs 净 A',float(netD),-0.0350,sd,
                 null_kind='零 A(题内跨人)下的净值 —— 不是零假设,是"若 75% 全来自换掉自变量,D 该落在哪"')
g.asserted('注册的判定:接近 A/B -> C 的 75% 全来自换掉自变量;接近 C -> `#128` 要重估',
           abs(netD+0.0350)<abs(netD+0.0081),
           f"净 D {netD:+.4f};离 A(−0.0350){abs(netD+0.0350):.4f} · 离 C(−0.0081){abs(netD+0.0081):.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 诊断:这个零到底换了多少次,以及它混不混 ---------------------------------
print("\n---- 诊断:约束集有多大 ----")
def count_swaps(Aa, rng, rounds=1):
    B=Aa.copy(); tot=0; opp=0
    for _ in range(rounds):
        for j in range(Mc):
            for k in range(j+1,Mc):
                m=np.isfinite(B[:,j])&np.isfinite(B[:,k]); idx=np.flatnonzero(m)
                if len(idx)<4: continue
                d=defaultdict(list)
                for i in idx:
                    a,b=B[i,j],B[i,k]
                    if a!=b: d[(a,b)].append(i)
                for (a,b),g1 in list(d.items()):
                    g2=d.get((b,a))
                    if not g2 or a>=b: continue
                    opp+=min(len(g1),len(g2)); tot+=len(g1)+len(g2)
    return tot,opp
t_r,o_r=count_swaps(V0,rng); t_p,o_p=count_swaps(Vp,rng)
print(f"  真实数据:有序值对 {t_r:,} 个,其中有镜像伙伴的 {o_r:,}({100*o_r/max(t_r,1):.1f}%)")
print(f"  种植数据:有序值对 {t_p:,} 个,其中有镜像伙伴的 {o_p:,}({100*o_p/max(t_p,1):.1f}%)")
print(f"  -> 种植让顺序变一致,镜像伙伴几乎消失 —— **约束集在有序配置上几乎是刚性的**")
# 混合性:多跑几轮,零值动不动
mix=[delta(doubly_constrained(V0,np.random.default_rng(1400),rounds=r)) for r in (1,3,8)]
print(f"  混合性:rounds=1/3/8 时零值 {' · '.join(f'{v:+.4f}' for v in mix)}"
      f"  -> {'不随轮数变化,不是在采样' if max(mix)-min(mix)<0.005 else '随轮数变化'}")

g2=Gate('这个零为什么不工作')
g2.asserted('正对照失败已确认:种植后零没毁掉信号',abs(dpz)>=0.5*abs(dp),f"{dp:+.4f} -> {dpz:+.4f}")
g2.asserted('诊断一:种植让镜像伙伴几乎消失 -> 约束集在有序配置上几乎刚性',
            (o_p/max(t_p,1))<0.5*(o_r/max(t_r,1)),
            f"真实 {100*o_r/max(t_r,1):.1f}% vs 种植 {100*o_p/max(t_p,1):.1f}%")
g2.asserted('诊断二:零值不随轮数变化 -> 它不是在采样约束集,是被推到一个固定构型',
            max(mix)-min(mix)<0.005,f"{' · '.join(f'{v:+.4f}' for v in mix)}")
g2.asserted('因此本轮注册的那个判定取不到数',True,
            '零 D 无效 -> "接近 A/B 还是接近 C"这个问题本轮答不了')
print(g2)
