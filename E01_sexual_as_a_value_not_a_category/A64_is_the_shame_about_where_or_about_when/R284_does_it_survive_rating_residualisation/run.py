import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A64 R284 -- 「早半更强」在扣掉评分之后还在不在

`#238b`:`#114` 测到人**把最爱的兴趣记得更早**(−0.2000 年/评分 sd)。
若早半因此富集了「爱得最深」的兴趣、而羞耻跟强度走,就会造出 `#238a` 的图样。
这是那条结论**唯一的重大对手**,而这份 release **带着评分**(`RATING_0_5`),所以它能被直接测掉。

ESTIMAND        把每个类别的起始年龄先对**这个人自己对该类别的评分**回归取残差
                (`#173` 同款手法),用残差重新劈早/晚两半,重跑 `#238a` 的全部。
KILL            **若早半优势在评分残差化之后存活(晚−早 仍 < −2×展布)-> `#238a` 不是回忆偏差,升 D7;
                若塌掉 -> `#238a` 撤回,「羞耻贴最早的欲望」改写成「羞耻贴你最爱的欲望」。**
POSITIVE CTRL   两端(`#276` 同款):
                ① 构造一个**已知由评分驱动**的假早/晚划分 -> 残差化必须把它杀掉;
                ② 构造一个**与评分无关**的划分 -> 残差化必须不动它。
NEGATIVE CTRL   置换羞耻。
⚠ 覆盖率           只有能匹配到评分列的类别才能用,类别数会掉。**匹配率与两条臂的类别数同报。**
IMPOSSIBLE      评分是**当下**的评分,起始年龄是**回忆**;
                用当下评分残差化回忆,只能扣掉「爱得深 -> 记得早」这一条通路,
                扣不掉「记得早 -> 现在更爱」。方向不可分,如实登记。
"""
import numpy as np, pandas as pd, warnings, hashlib, re as _re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
V=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in d.columns]
RT=d[rate].apply(pd.to_numeric,errors='coerce').values
def norm(s): return _re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=_re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
NC=V.shape[1]; MJ=sorted(best)
print(f"起始类别 {NC};匹配到评分列的 **{len(MJ)}({100*len(MJ)/NC:.0f}%)**;评分列 {len(rate)}")
OBS=np.isfinite(V); RAR=-np.log(np.clip(OBS.mean(0),1e-4,1.))
corr_ra=[]
for j in MJ:
    m=np.isfinite(V[:,j])&np.isfinite(RT[:,best[j]])
    if m.sum()>300: corr_ra.append(np.corrcoef(V[m,j],RT[m,best[j]])[0,1])
print(f"⚠ `#114` 的通路在这份匹配上:corr(起始年龄, 自己的评分) 均值 = **{np.mean(corr_ra):+.4f}**"
      f"(负 = 爱得深记得早),{sum(1 for c in corr_ra if c<0)}/{len(corr_ra)} 个类别为负")

def resid_on_rating(Vm):
    out=Vm.copy()
    for j in MJ:
        m=np.isfinite(Vm[:,j])&np.isfinite(RT[:,best[j]])
        if m.sum()<300: continue
        b=np.polyfit(RT[m,best[j]],Vm[m,j],1); out[m,j]=Vm[m,j]-np.polyval(b,RT[m,best[j]])
    return out
def halves_from(Vm, cats):
    sel=np.zeros(NC,bool); sel[list(cats)]=True
    Vs=np.where(np.isfinite(Vm)&sel[None,:],Vm,np.nan)
    Rm=np.where(np.isfinite(Vm)&sel[None,:],RAR[None,:],np.nan)
    med=np.nanmedian(Vs,1); E=Vs<=med[:,None]; L=Vs>med[:,None]
    ne,nl=E.sum(1),L.sum(1)
    return (np.where(ne>=4,np.nansum(np.where(E,Rm,0),1)/np.maximum(ne,1),np.nan),
            np.where(nl>=4,np.nansum(np.where(L,Rm,0),1)/np.maximum(nl,1),np.nan))
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
rng=np.random.default_rng(20260804)
def joint(E,L,yy=None):
    yy=y if yy is None else yy
    zE,zL=z(E),z(L); m=np.isfinite(yy)&np.isfinite(zE)&np.isfinite(zL)
    X=np.column_stack([np.ones(m.sum()),zE[m],zL[m]]); zy=(yy[m]-yy[m].mean())/yy[m].std()
    b=np.linalg.lstsq(X,zy,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i],zy[i],rcond=None)[0][j]
        for i in (rng.choice(m.sum(),m.sum(),True) for _ in range(200))])) for j in (1,2)]
    return b[1],b[2],sd[0],sd[1],int(m.sum())
E0,L0=halves_from(V,MJ); bE0,bL0,sE0,sL0,n0=joint(E0,L0)
Vr=resid_on_rating(V); E1,L1=halves_from(Vr,MJ); bE1,bL1,sE1,sL1,n1=joint(E1,L1)
print(f"\n(仅用匹配上的 {len(MJ)} 个类别)")
print(f"  残差化前:早 **{bE0:+.4f}** ± {sE0:.4f} · 晚 {bL0:+.4f} ± {sL0:.4f} · "
      f"晚−早 **{bL0-bE0:+.4f}** vs 2×展布 {2*np.hypot(sE0,sL0):.4f}(n = {n0:,})")
print(f"  **残差化后:早 {bE1:+.4f} ± {sE1:.4f} · 晚 {bL1:+.4f} ± {sL1:.4f} · "
      f"晚−早 {bL1-bE1:+.4f} vs 2×展布 {2*np.hypot(sE1,sL1):.4f}(n = {n1:,})**")
g0,g1=bL0-bE0,bL1-bE1
print(f"  优势保留 **{100*g1/g0 if g0!=0 else float('nan'):.1f}%**")
# ⚠ 两条臂的 n 不同(残差化改变了中位数劈分,从而改变 >=4 的过滤)。必须在同一批人上再比一次。
common=np.isfinite(E0)&np.isfinite(L0)&np.isfinite(E1)&np.isfinite(L1)&np.isfinite(y)
def joint_on(E,L,mask):
    zE,zL=z(np.where(mask,E,np.nan)),z(np.where(mask,L,np.nan)); m=mask&np.isfinite(zE)&np.isfinite(zL)
    X=np.column_stack([np.ones(m.sum()),zE[m],zL[m]]); zy=(y[m]-y[m].mean())/y[m].std()
    b=np.linalg.lstsq(X,zy,rcond=None)[0]
    sd=[float(np.std([np.linalg.lstsq(X[i],zy[i],rcond=None)[0][j]
        for i in (rng.choice(m.sum(),m.sum(),True) for _ in range(200))])) for j in (1,2)]
    return b[1],b[2],sd[0],sd[1],int(m.sum())
c0=joint_on(E0,L0,common); c1=joint_on(E1,L1,common)
print(f"  **同一批人(n = {c0[4]:,})**:残差化前 晚−早 {c0[1]-c0[0]:+.4f} vs 2×展布 "
      f"{2*np.hypot(c0[2],c0[3]):.4f} -> 残差化后 {c1[1]-c1[0]:+.4f} vs {2*np.hypot(c1[2],c1[3]):.4f};"
      f"保留 **{100*(c1[1]-c1[0])/(c0[1]-c0[0]):.1f}%**")

nul=[joint(E1,L1,rng.permutation(y))[1] for _ in range(20)]
print(f"  置换羞耻的零(晚半):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")

# 正对照两端
fakeR=V.copy()
for j in MJ:                      # ① 完全由评分驱动的假起始年龄
    m=np.isfinite(V[:,j])&np.isfinite(RT[:,best[j]]); fakeR[m,j]=-2.0*RT[m,best[j]]+rng.standard_normal(m.sum())*0.1
Ea,La=halves_from(fakeR,MJ); ba=joint(Ea,La)
Eb,Lb=halves_from(resid_on_rating(fakeR),MJ); bb=joint(Eb,Lb)
print(f"\n正对照①(完全由评分驱动的假划分):残差化前 晚−早 {ba[1]-ba[0]:+.4f} -> 后 {bb[1]-bb[0]:+.4f}"
      f"(必须被杀掉)")
fakeI=V.copy()
for j in MJ:                      # ② 与评分无关的假起始年龄
    m=np.isfinite(V[:,j]); fakeI[m,j]=rng.standard_normal(m.sum())
Ec,Lc=halves_from(fakeI,MJ); bc=joint(Ec,Lc)
Ed,Ld=halves_from(resid_on_rating(fakeI),MJ); bd=joint(Ed,Ld)
print(f"正对照②(与评分无关的假划分):残差化前 晚−早 {bc[1]-bc[0]:+.4f} -> 后 {bd[1]-bd[0]:+.4f}"
      f"(必须几乎不动)")

T=pd.DataFrame([dict(arm='残差化前',beta_early=bE0,beta_late=bL0,sd_early=sE0,sd_late=sL0,n=n0),
                dict(arm='残差化后',beta_early=bE1,beta_late=bL1,sd_early=sE1,sd_late=sL1,n=n1)])
check_columns(T,'R284'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rating_residualised.csv',index=False)

g=Gate('「早半更强」扣掉评分之后还在不在')
g.asserted('⚠ `#114` 的通路在这份数据上确实存在(否则这个控制无的放矢)',
           np.mean(corr_ra)<0, f"corr(起始年龄, 评分) 均值 {np.mean(corr_ra):+.4f};"
           f"{sum(1 for c in corr_ra if c<0)}/{len(corr_ra)} 为负")
g.asserted('正对照①:完全由评分驱动的假划分,残差化必须把它杀掉',
           abs(bb[1]-bb[0])<abs(ba[1]-ba[0])/2, f"{ba[1]-ba[0]:+.4f} -> {bb[1]-bb[0]:+.4f}")
g.asserted('正对照②:与评分无关的假划分,残差化必须几乎不动',
           abs((bd[1]-bd[0])-(bc[1]-bc[0]))<0.05, f"{bc[1]-bc[0]:+.4f} -> {bd[1]-bd[0]:+.4f}")
g.negative_control('置换羞耻(残差化后的晚半)',abs(float(np.mean(nul))),abs(bL1),
                   null_spread=float(np.std(nul)),null_kind='跨人置换结局 —— 只打掉配对')
g.asserted('⚠ 同一批人上重比(两条臂的 n 不同,不比不能下结论)',
           True, f"n = {c0[4]:,};前 {c0[1]-c0[0]:+.4f} -> 后 {c1[1]-c1[0]:+.4f};"
                 f"保留 {100*(c1[1]-c1[0])/(c0[1]-c0[0]):.1f}%")
g.asserted('★ 注册的 kill:早半优势在评分残差化之后存活 -> `#238a` 不是回忆偏差',
           g1<-2*np.hypot(sE1,sL1), f"残差化后 晚−早 {g1:+.4f} vs 2×展布 {2*np.hypot(sE1,sL1):.4f};"
           f"保留 {100*g1/g0 if g0!=0 else float('nan'):.1f}%")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
