import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A19 R02 -- 接缝检验:「ρ 与 z_resid 是同一件事」这句话,跨到 S 上还成立吗?

#149a:同一批类别上 corr(z_resid, ρ) = **−0.5741**(91.6× 地板)—— 两个量读的是同一件事。
#129e:`corr(ρ, S) = −0.0417`(按类别数卡钳匹配后,3.1×),是一条**现存声明**。

**预测(跑之前写死):** 既然 corr(z,ρ) 是**负**的,那么
    corr(z_resid, S)  应当与  corr(ρ, S)  **符号相反、量级相当**,
即 corr(z_resid,S) ≈ **+0.04**,且 `corr(z_resid,S)` 与 `−corr(ρ,S)` 应在各自自助展布内相符。

    HOLDS  两者在展布内相符 -> 接缝是紧的,#129e 的读法不变,而本项目两条主线正式并成一条
    BREAKS 不相符 -> 「同一件事」这句话跨到 S 时破了。**那么 #129e 说的到底是 ρ 的性质
           还是那个共同成分的性质,就不再清楚**,而这条现存声明的读法必须改

ESTIMAND        corr(z_resid, S) 与 corr(ρ, S),都按**答题类别数**卡钳 1:1 匹配(#129e 的做法);
                判别量 = corr(z_resid,S) + corr(ρ,S)(若"同一件事"成立,它应当 ≈ 0)。
IDENTIFICATION  S 来自**多选题选项**,与起始年龄题目零重叠;S 对勾选数的残差化用
                `check_residualized` 断言(#129 的教训:那一行错过一次,整轮作废)。
SCOPE           >=8 个类别起始年龄、S 可算的人。
WORLDS          HOLDS / BREAKS
KILL            条件式:匹配必须真的把类别数差压到 <0.1 sd,且两个量在同一批人上算,
                才读判别量。
POSITIVE CTRL   种植一个同时驱动 z 与 ρ 的人特异径向信号,并让 S 与它相关 ——
                两个相关必须都被拉动且符号相反。
NEGATIVE CTRL   按人置换 S 的标签(打掉 S 与任何人内量的联系)。
NOISE FLOOR     按人自助 200 次;5 个匹配种子。
MULTIPLICITY    2 个量 x {未匹配, 匹配},整格发表。
IMPOSSIBLE      S 与 z/ρ 的因果方向;以及 S 本身的仪器成分(#5 的覆盖度定律)。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

def demean_conv(Vm,tol=1e-10,cap=500):
    D=np.where(obs,Vm,np.nan)
    for _ in range(cap):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-a
        b=np.nanmean(D,axis=1,keepdims=True); D=D-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return D
NPERM=200
def two_stats(Dres,seed):
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    Z=np.full(len(Dres),np.nan); Rho=np.full(len(Dres),np.nan)
    for i in np.flatnonzero(KEEP):
        j=np.flatnonzero(obs[i]); y=Dres[i,j]; r=rar[j]; k=len(j)
        if k<4: continue
        cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
        d=r[pick]-r.mean()
        idx=rg.integers(0,k,(NPERM,1)); dr=r[idx].mean(1)-r.mean()
        if dr.std()<1e-9: continue
        Z[i]=(d-dr.mean())/dr.std()
        s=np.nanstd(y)
        if s<1e-9: continue
        Rho[i]=np.corrcoef(y,r)[0,1]
    return Z,Rho

Dres=demean_conv(V)
Z,Rho=two_stats(Dres,zlib.crc32(b'A19R02'))
base=KEEP&np.isfinite(Z)&np.isfinite(Rho)&np.isfinite(S)
ii=np.flatnonzero(base); NC=NCAT.astype(float)
print(f"同一批人 {len(ii):,};corr(z,ρ) = {np.corrcoef(Z[ii],Rho[ii])[0,1]:+.4f}",flush=True)

def matched(vec,seed):
    rg=np.random.default_rng(seed); med=np.median(S[ii])
    hi=ii[S[ii]>med]; lo=ii[S[ii]<=med]
    c=(NC-NC[ii].mean())/NC[ii].std(); used=np.zeros(len(S),bool); P=[]
    for a in hi[rg.permutation(len(hi))]:
        d_=np.abs(c[lo]-c[a]); d_[used[lo]]=np.inf; k=int(np.argmin(d_))
        if d_[k]<0.25: used[lo[k]]=True; P.append((a,lo[k]))
    P=np.array(P); sel=np.r_[P[:,0],P[:,1]]
    return (float(np.corrcoef(vec[sel],S[sel])[0,1]),
            abs(NC[P[:,0]].mean()-NC[P[:,1]].mean())/NC[ii].std(), len(P))

rb=np.random.default_rng(5)
def boot(vec):
    return float(np.std([ (lambda s_: np.corrcoef(vec[s_],S[s_])[0,1])(
        ii[rb.integers(0,len(ii),len(ii))]) for _ in range(200)]))
rows=[]
for nm,vec in [('z_resid',Z),('rho',Rho)]:
    raw=float(np.corrcoef(vec[ii],S[ii])[0,1])
    mv=[matched(vec,600+s) for s in range(5)]
    mc=float(np.mean([x[0] for x in mv])); bal=float(np.mean([x[1] for x in mv]))
    rows.append(dict(stat=nm,raw=raw,matched=mc,bal=bal,boot=boot(vec),n=int(np.mean([x[2] for x in mv]))))
    print(f"  {nm:<9} 未匹配 {raw:+.4f}   匹配后 {mc:+.4f}   类别数残差 {bal:.3f} sd   "
          f"自助展布 {rows[-1]['boot']:.4f}",flush=True)

T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'seam.csv',index=False)
cz=T[T.stat=='z_resid'].iloc[0]; cr=T[T.stat=='rho'].iloc[0]
disc=float(cz.matched+cr.matched); sd=float(np.sqrt(cz.boot**2+cr.boot**2))
print(f"\n  预测:corr(z_resid,S) ≈ −corr(ρ,S) = {-cr.matched:+.4f}")
print(f"  实测:corr(z_resid,S) = {cz.matched:+.4f}")
print(f"  判别量 corr(z,S)+corr(ρ,S) = **{disc:+.4f}**(若两者是同一件事应 ≈ 0),"
      f"展布 {sd:.4f} -> {abs(disc)/sd:.1f}x")

# 正对照:种植一个同时驱动 z 与 ρ 且与 S 相关的信号
rgp=np.random.default_rng(31); x=rar-rar.mean()
u=np.where(np.isfinite(S),S,0.)+rgp.normal(0,0.5,len(S))
Dp=demean_conv(np.where(obs,V+1.2*np.outer(u,x),np.nan))
Zp,Rp=two_stats(Dp,555); mp=KEEP&np.isfinite(Zp)&np.isfinite(Rp)&np.isfinite(S); jp=np.flatnonzero(mp)
pz=float(np.corrcoef(Zp[jp],S[jp])[0,1]); pr=float(np.corrcoef(Rp[jp],S[jp])[0,1])
print(f"\n  正对照(种植与 S 相关的径向信号):corr(z,S) {pz:+.4f}  corr(ρ,S) {pr:+.4f}  "
      f"判别量 {pz+pr:+.4f}")
rgn=np.random.default_rng(77); Sp=S.copy(); Sp[ii]=S[ii][rgn.permutation(len(ii))]
nz=float(np.corrcoef(Z[ii],Sp[ii])[0,1]); nr=float(np.corrcoef(Rho[ii],Sp[ii])[0,1])
print(f"  负对照(按人置换 S):corr(z,S) {nz:+.4f}  corr(ρ,S) {nr:+.4f}")

g=Gate('「ρ 与 z_resid 是同一件事」跨到 S 上还成立吗')
g.asserted('匹配把类别数差压下去了',float(cz.bal)<0.1 and float(cr.bal)<0.1,
           f"z {cz.bal:.3f} sd / ρ {cr.bal:.3f} sd")
g.asserted('正对照:种植的信号让两个相关符号相反且都被拉动',pz*pr<0 and abs(pz)>0.1 and abs(pr)>0.1,
           f"corr(z,S) {pz:+.4f} · corr(ρ,S) {pr:+.4f}")
g.negative_control('按人置换 S 后两个相关都归零',max(abs(nz),abs(nr)),max(abs(cz.matched),abs(cr.matched)))
g.no_sign_crossing('两个相关符号相反(预测的方向)',[float(cz.matched),-float(cr.matched)])
# ⚠ 假设本身是**"两者相同"**,所以 require_resolvable_first 用反了 ——
#   它为"我要它非零"设计,会把想要的结果报成 FAIL 并把整族标 MOOT。正确的是**等价界**。
MARGIN=0.5*abs(float(cr.matched))          # 预设边界:效应的一半
g.equivalent_within('判别量落在预设等价边界内(效应的一半)',disc,sd,MARGIN)
g.asserted('把这个设计能排除的差异报出来(#P14 MDE)',True,
           f"判别量 {disc:+.4f},95% 上界 {abs(disc)+2*sd:.4f};效应本身 {abs(cz.matched):.4f} —— "
           f"**本设计只能排除大于效应 {100*(abs(disc)+2*sd)/abs(cz.matched):.0f}% 的差异**")
g.asserted('预注册的点预测是否命中',abs(cz.matched-(-cr.matched))<2*sd,
           f"预测 {-cr.matched:+.4f},实测 {cz.matched:+.4f},差 {cz.matched+cr.matched:+.4f}"
           f"(展布 {sd:.4f})")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
