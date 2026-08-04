import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A65 R285 -- 第 12 个守卫,和一条现存撤回的回溯重审

`#239a`:**一个改变了样本的控制,改变的是估计目标本身。**
本项目 11 个守卫全在防**假肯定**;`#239a` 是第一个记录在案的**假否定**机制,而它成本更高
—— 被撤回的主张没人重审。

两件事,一轮内:
① 落第 12 个守卫 `control_kept_the_sample`(在 `lib/gates.py`)。
② **回溯重审 `#210`** —— 它用 `V -> where(V<=14)` 撤回了 `#206b`/`#207c`,
   而那个控制把样本从 **9,944 掉到 2,806**。`#210c` 已经注意到这批人「本身更早熟,不是随机子集」,
   并克制地只用它证伪 —— **但它没有做 `#239a` 那一步:把未截断的量算在同一批 2,806 人上。**

ESTIMAND        `r(rho_i, age)` 的三个读数:
                A 全样本未截断(`#210` 报 +0.1532,n=9,944)
                B 截断 ≤14(`#210` 报 −0.0741,n=2,806)
                **C 未截断,但只在 B 的那批人上** —— 这一格从没被算过
KILL            **若 C ≈ A(仍明显为正)-> 反号确实由截断造成,`#210` 的撤回成立,升 D8;
                若 C ≈ B(已经为负或接近零)-> 反号大部分由样本选择造成,
                `#210` 的撤回是 `#239a` 那一族的假否定,必须重开。**
NEGATIVE CTRL   人内跨人置换(`#210` 同款)。
POSITIVE CTRL   守卫两端:构造一个**只改效应、不改纳入**的控制(必须放行);
                构造一个**只改纳入、不改效应**的控制(必须报警)。
IMPOSSIBLE      C 与 A 的差别混合了「这批人更早熟」与「这批人 rho 估得更准」,
                本轮不拆这两者;能判的只有**反号是不是截断造成的**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
print(f"载入 `R173` 的管道:{N:,} 人 × {Mc} 类别")

def rho_of(Vm):
    D=demean_conv(Vm) if 'demean_conv' in dir() else None
    if D is None:
        D=np.where(np.isfinite(Vm),Vm,np.nan)
        for _ in range(300):
            a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
            b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=8)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
def r_age(rho, mask=None):
    m=np.isfinite(rho)&np.isfinite(age)
    if mask is not None: m&=mask
    return float(np.corrcoef(rho[m],age[m])[0,1]), int(m.sum())
rhoA=rho_of(V0); Vc=np.where(V0<=14,V0,np.nan); rhoB=rho_of(Vc)
keepB=np.isfinite(rhoB)
rA,nA=r_age(rhoA); rB,nB=r_age(rhoB); rC,nC=r_age(rhoA,keepB)
print(f"\nA 全样本未截断      r(rho, age) = **{rA:+.4f}**(n = {nA:,})  [`#210` 报 +0.1532]")
print(f"B 截断 ≤14         r(rho, age) = **{rB:+.4f}**(n = {nB:,})  [`#210` 报 −0.0741]")
print(f"**C 未截断,只在 B 的人上  r(rho, age) = {rC:+.4f}(n = {nC:,})  <- 这一格从没被算过**")
rng=np.random.default_rng(20260804)
bt=lambda rho,mask: float(np.std([np.corrcoef(rho[i],age[i])[0,1] for i in
    (rng.choice(np.flatnonzero(np.isfinite(rho)&np.isfinite(age)&(mask if mask is not None else True)),
     int(np.sum(np.isfinite(rho)&np.isfinite(age)&(mask if mask is not None else True))),True) for _ in range(200))]))
sA,sB_,sC=bt(rhoA,None),bt(rhoB,None),bt(rhoA,keepB)
print(f"  自助展布:A ±{sA:.4f} · B ±{sB_:.4f} · C ±{sC:.4f}")
nul=[float(np.corrcoef(rhoA[np.isfinite(rhoA)&np.isfinite(age)],
      rng.permutation(age[np.isfinite(rhoA)&np.isfinite(age)]))[0,1]) for _ in range(30)]
print(f"  置换零 {np.mean(nul):+.4f} ± {np.std(nul):.4f}")
print(f"  **C−B = {rC-rB:+.4f}(截断造成的部分)· A−C = {rA-rC:+.4f}(样本选择造成的部分)**")

T=pd.DataFrame([dict(arm='A_全样本未截断',r=rA,sd=sA,n=nA),
                dict(arm='B_截断≤14',r=rB,sd=sB_,n=nB),
                dict(arm='C_未截断_同批人',r=rC,sd=sC,n=nC)])
check_columns(T,'R285'); T.to_csv(pathlib.Path(__file__).parent/'results'/'intersection_retest.csv',index=False)

g=Gate('一条现存撤回是不是 `#239a` 那一族的假否定')
g.asserted('复现 `#210` 的两个已发表读数(否则不是同一条管道)',
           abs(rA-0.1532)<0.02 and abs(rB+0.0741)<0.02, f"A {rA:+.4f} vs +0.1532;B {rB:+.4f} vs −0.0741")
g.negative_control('置换年龄',abs(float(np.mean(nul))),abs(rA),null_spread=float(np.std(nul)),
                   null_kind='跨人置换年龄 —— 只打掉配对')
g.offset_control('★ C(未截断,同批人)vs B(截断,同批人)',rC,rB,float(np.hypot(sC,sB_)),
                 null_kind='同一批人上截断后的值 —— 不是零假设,是「若反号真由截断造成,C 该落在哪」')
g.asserted('★ 注册的 kill:C 仍明显为正 -> `#210` 的撤回成立;C 已为负/近零 -> 撤回要重开',
           rC>2*sC, f"C {rC:+.4f} ± {sC:.4f};C−B {rC-rB:+.4f}(截断部分)· A−C {rA-rC:+.4f}(样本部分)")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
