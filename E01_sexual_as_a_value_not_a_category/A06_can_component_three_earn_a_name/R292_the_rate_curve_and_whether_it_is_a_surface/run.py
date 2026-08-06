import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A67 R292 -- 速率曲线,以及这是一条线还是一张曲面

`#246a`:`rate`(积累速率)是五个钟里最可分辨的(18.0×),但它现在只是一个联合 beta ——
**没人知道这条关系长什么样,也不知道它在哪里饱和或反转。**

ESTIMAND        ① 按 `rate` 分 8 层,各层报 Δ、它自己的置换零、**以及 n**;
                ② `rate × A_close` 的 4×3 曲面,每格报 Δ 与 n;
                ③ 交互项 `rate × A_close` 的联合 beta。
KILL            **若 Δ 沿 rate 单调且各 `A_close` 层内平行(交互项 < 2×展布)-> 两个量可加;
                若曲线反转或各层不平行 -> 存在交互,「形状」要具体化成一个交互项。**
POSITIVE CTRL   两端:① 只依赖 rate 的种入 -> 曲线必须出现且交互项不动;
                ② 种一个真交互 -> 交互项必须开火。
NEGATIVE CTRL   人内跨人置换。
⚠ 偏分布         `rate` 很可能极偏。**每层的 n 必须报**;
                n < 400 的层标注为**不可读**,不参与单调性判定。
IMPOSSIBLE      `rate = (类别数−1)/(最后−最初)`,分母是回忆的年龄跨度;
                跨度很小的人(集中报告)会得到极大的 rate,那既可能是真的快,也可能是回忆压缩。
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
A1=np.where(base,Vs[:,0],np.nan)
ACL=np.where(base,np.nanmax(np.where(np.isfinite(V0),V0,np.nan),axis=1),np.nan)
RATE=np.where(base,(NCAT-1)/np.maximum(ACL-A1,0.5),np.nan)
m0=base&np.isfinite(RATE)&np.isfinite(ACL)
print(f"n = {int(m0.sum()):,};rate 分位 "
      + ' · '.join(f"{p}%:{np.nanpercentile(RATE[m0],p):.2f}" for p in (1,10,50,90,99)))
rng=np.random.default_rng(20260804)
def D_of(mask,r=None):
    r=rho if r is None else r; mm=mask&np.isfinite(r); v=r[mm]
    if len(v)<50: return np.nan,np.nan,int(mm.sum())
    return (float(np.mean(v)),
            float(np.std([np.mean(v[i]) for i in (rng.choice(len(v),len(v),True) for _ in range(200))])),
            int(mm.sum()))
def curve(nb=8,r=None):
    q=np.nanpercentile(RATE[m0],np.linspace(0,100,nb+1)); out=[]
    for lo,hi in zip(q[:-1],q[1:]):
        mm=m0&(RATE>lo)&(RATE<=hi) if lo>q[0] else m0&(RATE<=hi)
        d=D_of(mm,r); out.append((float(np.nanmedian(RATE[mm])),)+d)
    return out
C=curve()
print(f"\nΔ 沿 rate(8 层):")
for rt,d,s,n in C:
    tag='' if n>=400 else '  ⚠ n<400 不可读'
    print(f"  rate 中位 {rt:>6.2f}  Δ = {d:+.4f} ± {s:.4f}  (n={n:,}){tag}")
ok_lv=[(rt,d) for rt,d,s,n in C if n>=400]
mono=all(ok_lv[i][1]>=ok_lv[i+1][1]-0.005 for i in range(len(ok_lv)-1))
print(f"  可读层单调不升 = **{mono}**;首末差 **{ok_lv[-1][1]-ok_lv[0][1]:+.4f}**")

qa=np.nanpercentile(ACL[m0],[33,67]); qr=np.nanpercentile(RATE[m0],[25,50,75])
print(f"\nrate × A_close 曲面(4×3,格内为 Δ,括号为 n):")
hdr=['A_close 低','A_close 中','A_close 高']
print(f"{'rate':>10}"+''.join(f"{h:>20}" for h in hdr))
grid=[]
for i,(lo,hi) in enumerate(zip([-np.inf]+list(qr),list(qr)+[np.inf])):
    row=[]
    for lo2,hi2 in zip([-np.inf]+list(qa),list(qa)+[np.inf]):
        mm=m0&(RATE>lo)&(RATE<=hi)&(ACL>lo2)&(ACL<=hi2); d=D_of(mm); row.append(d)
    grid.append(row)
    print(f"{('Q%d'%(i+1)):>10}"+''.join(f"{d:>+12.4f}({n:>5,})" for d,s,n in row))
zr=(rho[m0]-rho[m0].mean())/rho[m0].std()
zR=(RATE[m0]-RATE[m0].mean())/RATE[m0].std(); zA=(ACL[m0]-ACL[m0].mean())/ACL[m0].std()
def fit(cols,y=None):
    y=zr if y is None else y
    X=np.column_stack([np.ones(len(y))]+cols); b=np.linalg.lstsq(X,y,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i],y[i],rcond=None)[0][j]
        for i in (rng.choice(len(y),len(y),True) for _ in range(200))])) for j in range(1,len(cols)+1)]
    return b[1:],sd
b,s=fit([zR,zA,zR*zA])
print(f"\n交互模型:rate {b[0]:+.4f}±{s[0]:.4f}({abs(b[0])/s[0]:.1f}×)· "
      f"A_close {b[1]:+.4f}±{s[1]:.4f}({abs(b[1])/s[1]:.1f}×)· "
      f"**交互 {b[2]:+.4f}±{s[2]:.4f}({abs(b[2])/s[2]:.1f}×)**")
u1=(zR*0.9+rng.standard_normal(len(zr))*0.4); u1=(u1-u1.mean())/u1.std()
u2=(zR*zA*0.9+rng.standard_normal(len(zr))*0.4); u2=(u2-u2.mean())/u2.std()
p1=fit([zR,zA,zR*zA],u1); p2=fit([zR,zA,zR*zA],u2)
print(f"正对照①(只依赖 rate):交互 {p1[0][2]:+.4f}±{p1[1][2]:.4f}(必须 ≈0)· rate {p1[0][0]:+.4f}")
print(f"正对照②(真交互):交互 **{p2[0][2]:+.4f}**±{p2[1][2]:.4f}(必须开火)")
nul=[]
for _ in range(12):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    r2=rho_of(Vp); yy=(r2[m0]-np.nanmean(r2[m0]))/np.nanstd(r2[m0]); nul.append(fit([zR,zA,zR*zA],yy)[0][2])
print(f"置换零(交互项){np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(rate_median=rt,delta=d,sd=s_,n=n,readable=n>=400) for rt,d,s_,n in C])
check_columns(T,'R292'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rate_curve.csv',index=False)

g=Gate('速率曲线是一条线还是一张曲面')
g.asserted('正对照两端:只依赖 rate 的种入交互必须 ≈0;真交互必须开火',
           abs(p1[0][2])<2*p1[1][2] and abs(p2[0][2])>4*p2[1][2],
           f"① 交互 {p1[0][2]:+.4f}±{p1[1][2]:.4f} · ② 交互 {p2[0][2]:+.4f}±{p2[1][2]:.4f}")
g.asserted('⚠ 每层 n 已报,n<400 的层标注不可读且不参与单调性判定',
           True, ' · '.join(f"{n:,}{'' if n>=400 else '⚠'}" for _,_,_,n in C))
g.negative_control('置换零(交互项)',abs(float(np.mean(nul))),abs(float(b[2])),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.asserted('★ 注册的 kill:曲线单调且交互项 < 2×展布 -> 可加;否则存在交互',
           mono and abs(b[2])<2*s[2],
           f"单调={mono};交互 {b[2]:+.4f} ± {s[2]:.4f}({abs(b[2])/s[2]:.1f}×)")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
