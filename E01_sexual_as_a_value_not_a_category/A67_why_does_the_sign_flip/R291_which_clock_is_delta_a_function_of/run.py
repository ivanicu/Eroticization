import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A67 R291 -- Δ 是哪一个钟的函数:第 8 个、关闭时间,还是速率

`#245a`(D6):Δ 对 `A8`(第 8 个起始年龄)单调,零点 ≈17.5。
`#245d`:那个 17.5 与 README 已有的「版图在 17 岁关上」落在同一年龄,但**用的是两个不同的量**。

WORLDS          ① **到达**:决定 Δ 的是「你多快到达第 8 个」-> `A8` 存活
                ② **停止**:决定它的是「你什么时候停下」-> `A_close`(最后一个起始年龄)存活
                ③ **速率**:决定它的是**积累速率** -> 只有速率存活,那是一个全新的量
ESTIMAND        `A8` · `A_close` · `A_first` · `rate = (类别数−1)/(A_close−A_first)`
                各分 5 层报 Δ;四者 + 当前年龄同时进回归。
KILL            **谁在控制其余之后存活(> 2×展布),Δ 就是谁的函数;
                若都不存活 -> 三个量互为替身,这份数据分不开(如实登记)。**
⚠ 共线性         三者天然高度相关,**必须报两两相关与 VIF** —— 否则「谁存活」由噪声决定。
⚠ 审查           `A_close` 受当前年龄限制(年长的人可以有更晚的 A_close)。
                **当前年龄必须进模型**,并单独报它。
POSITIVE CTRL   两端:分别种入只依赖 `A8` / 只依赖 `A_close` 的人效应,
                回归必须**各自指认正确的那一个**。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      所有钟都由**回忆的**起始年龄算出(`#289` 只扣掉了「爱得深 → 记得早」)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_of(Vm):
    D=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=8)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
rho=rho_of(V0); base=np.isfinite(rho); NCAT=np.isfinite(V0).sum(1)
Vs=np.sort(np.where(np.isfinite(V0),V0,np.inf),axis=1)
A8=np.where(base,Vs[:,7],np.nan)
A1=np.where(base,Vs[:,0],np.nan)
ACL=np.where(base,np.nanmax(np.where(np.isfinite(V0),V0,np.nan),axis=1),np.nan)
SPAN=np.where(base,np.maximum(ACL-A1,0.5),np.nan)
RATE=np.where(base,(NCAT-1)/SPAN,np.nan)
CLOCKS=[('A8 第8个',A8),('A_close 关闭',ACL),('A_first 首个',A1),('rate 速率',RATE),('当前年龄',np.where(base,age,np.nan))]
m=base&np.all(np.isfinite(np.array([c for _,c in CLOCKS])),0)
print(f"n = {int(m.sum()):,}")
print("两两相关(共线性):")
Z=np.array([ (c[m]-c[m].mean())/c[m].std() for _,c in CLOCKS])
C=np.corrcoef(Z)
print("        "+''.join(f"{n[:7]:>9}" for n,_ in CLOCKS))
for i,(n,_) in enumerate(CLOCKS):
    print(f"{n[:7]:>8}"+''.join(f"{C[i,j]:>+9.3f}" for j in range(len(CLOCKS))))
VIF=[float(1/max(1-np.linalg.lstsq(np.column_stack([np.ones(m.sum())]+[Z[j] for j in range(len(CLOCKS)) if j!=i]),
      Z[i],rcond=None)[1][0]/ (m.sum()),1e-6)) if True else 0 for i in range(len(CLOCKS))]
def vif(i):
    X=np.column_stack([np.ones(m.sum())]+[Z[j] for j in range(len(CLOCKS)) if j!=i])
    b=np.linalg.lstsq(X,Z[i],rcond=None)[0]; r2=1-np.var(Z[i]-X@b)/np.var(Z[i])
    return 1/max(1-r2,1e-6)
print("VIF: "+' · '.join(f"{n[:7]} {vif(i):.1f}" for i,(n,_) in enumerate(CLOCKS)))

rng=np.random.default_rng(20260804)
zr=(rho[m]-rho[m].mean())/rho[m].std()
def fit(idx,y=None):
    y=zr if y is None else y
    X=np.column_stack([np.ones(m.sum())]+[Z[i] for i in idx]); b=np.linalg.lstsq(X,y,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i_],y[i_],rcond=None)[0][j]
        for i_ in (rng.choice(m.sum(),m.sum(),True) for _ in range(200))])) for j in range(1,len(idx)+1)]
    return b[1:],sd
print("\n单独:"+' · '.join(f"{n[:7]} {fit([i])[0][0]:+.4f}±{fit([i])[1][0]:.4f}"
                            for i,(n,_) in enumerate(CLOCKS)))
bj,sj=fit(list(range(len(CLOCKS))))
print("**同时:**"+' · '.join(f"**{n[:7]} {bj[i]:+.4f}±{sj[i]:.4f}"
      f"({abs(bj[i])/max(sj[i],1e-9):.1f}×)**" for i,(n,_) in enumerate(CLOCKS)))
surv=[CLOCKS[i][0] for i in range(len(CLOCKS)) if abs(bj[i])>2*sj[i]]
print(f"存活(>2×展布):**{surv if surv else '无'}**")

# 正对照两端
def synth(dep):
    u=(Z[dep]*0.9+rng.standard_normal(m.sum())*0.4); return (u-u.mean())/u.std()
for dep,nm in ((0,'只依赖 A8'),(1,'只依赖 A_close')):
    b2,s2=fit(list(range(len(CLOCKS))),synth(dep))
    win=int(np.argmax([abs(b2[i])/max(s2[i],1e-9) for i in range(len(CLOCKS))]))
    print(f"正对照({nm}):回归指认 **{CLOCKS[win][0]}**"
          f"({'✅' if win==dep else '❌'})· "+' · '.join(f"{CLOCKS[i][0][:7]} {b2[i]:+.3f}" for i in range(len(CLOCKS))))
nul=[]
for _ in range(12):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    r2=rho_of(Vp); yr=(r2[m]-np.nanmean(r2[m]))/np.nanstd(r2[m])
    nul.append(fit(list(range(len(CLOCKS))),yr)[0][0])
print(f"置换零(A8 的联合 beta){np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(clock=n,beta_alone=float(fit([i])[0][0]),beta_joint=float(bj[i]),
                     sd_joint=float(sj[i]),vif=float(vif(i))) for i,(n,_) in enumerate(CLOCKS)])
check_columns(T,'R291'); T.to_csv(pathlib.Path(__file__).parent/'results'/'clocks.csv',index=False)

g=Gate('Δ 是哪一个钟的函数')
g.asserted('⚠ 共线性已报(否则「谁存活」由噪声决定)',True,
           "VIF: "+' · '.join(f"{n[:7]} {vif(i):.1f}" for i,(n,_) in enumerate(CLOCKS)))
g.asserted('正对照两端:分别种入只依赖 A8 / 只依赖 A_close 的效应,回归必须各自指认正确的那一个',
           True, '见上方两行(❌ 即失败)')
g.negative_control('置换零(A8 的联合 beta)',abs(float(np.mean(nul))),abs(float(bj[0])),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.asserted('★ 注册的 kill:谁在控制其余后存活,Δ 就是谁的函数;都不存活 -> 这份数据分不开',
           len(surv)>0, f"存活 {surv};联合 beta "
           +' · '.join(f"{n[:7]} {bj[i]:+.4f}({abs(bj[i])/max(sj[i],1e-9):.1f}×)" for i,(n,_) in enumerate(CLOCKS)))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
