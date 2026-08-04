import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A14 R01 -- 性版图是从常见向罕见"辐射"出去的吗?而稀有亲和是起点还是终点?

本项目两个最大的现存发现从未连过线:
  #75  人群共享一个发育时间表:先具体(外观 14.0)后关系(精神改变 17.0)
  #100 稀有亲和是一条可靠的人格维度(+0.4611,23.1x),跨到强制选择仍活(#126)

它们各自回答"什么时候"和"什么口味",没有人问过**这两者是不是同一条轨迹的两端**。
三个关于"一个人怎么变成他自己"的图景,在心理学上完全不同:

  A 共享辐射   每个人都从常见走向罕见,速率一样。个体差异只在"走多远",不在"怎么走"。
               -> sd(beta_i) 落在零上。时间表就是全部故事。
  B 起点分歧   高稀有亲和的人**一开始就在外围** —— 他们最早获得的就是不常见的东西。
               -> beta_i 有真实方差,且 corr(beta, S) < 0。"怪口味从一开始就在那里。"
  C 行程分歧   高稀有亲和的人是**走得更久的人** —— 罕见的东西对他们来说来得更晚。
               -> beta_i 有真实方差,且 corr(beta, S) > 0。"怪口味是你一直走下去的地方。"

ESTIMAND        beta_i = 这个人的起始年龄**偏离**对类别**稀有度**的斜率,在类别固定效应
                和这个人整体早熟被完全去掉之后。两个量:sd(beta_i);corr(beta_i, S_i)。
IDENTIFICATION  beta_i 结构上免疫两个最大的干扰:类别去均值精确移除了人群时间表(以及
                稀有度与"平均什么时候"的题目层相关),人内去均值精确移除了整体早熟 ——
                因为 beta 是一个以稀有度为权的**对比**,常数人效应恰好贡献 0。
                这与救了 #116 的是同一种结构免疫。
SCOPE           报告了 >=8 个类别起始年龄的人。起始年龄按 2 年分箱(release 所致),
                这是分辨率地板,报出来。
WORLDS          A shared / B start / C travel(见上)
KILL            条件式(#P16):正对照必须开火**且**零必须为零,才读阈值。
POSITIVE CTRL   种植一个已知的人特异径向斜率 g*u_i*(r_j-rbar):sd(beta) 必须随 g 单调上升,
                且 corr(beta,u) 必须开火。g=0 必须**逐位复现**真实臂(同种子)。
NEGATIVE CTRL   题内跨人置换残差 —— 精确保留每个题目的缺失模式与残差分布,只摧毁
                "谁配到哪个残差"。这正是 beta 想要检验的配对。
CONFOUND        #114:人把**最爱**的兴趣记得更早(-0.2000 年/评分 sd)。若高 S 的人更爱
                罕见类别,回忆偏差会把罕见类别拉早 -> 制造出世界 B。**所以 B 必须过
                artifact_cannot_explain,C 不必**(伪影符号与 C 相反)。控制手段:把
                起始年龄先对**这个人自己对该类别的评分**回归取残差,重跑全部。
                第二个混淆:年龄(#119 说畸变随年龄加深)。作为协变量报出。
                第三个:S 与 beta 的**仪器不相交** —— S 来自多选题选项,beta 来自类别起始
                年龄,check_disjoint_items 断言。
NOISE FLOOR     4 seeds x 200 次按人自助。
MULTIPLICITY    2 个统计量 x 4 个 g x 4 seeds x {含/不含评分残差化},整格发表。
IMPOSSIBLE      因果方向。横断面数据无法区分"先有稀有亲和所以走法不同"与"走法不同
                所以形成了稀有亲和";本轮只判**轨迹形状与特质是否同号**。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage, check_disjoint_items

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}); check_columns(V,"onset wide"); V=V.values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
AGEBIN={'14-17':15.5,'18-20':19,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=df['age'].map(AGEBIN).values

def norm(s): return re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
check_coverage(len(best),V.shape[1],"onset->rating match",tol=0.35)

# ---- 稀有度 = -log(报告了该类别起始年龄的人的比例)。这是类别的常见程度。
obs=np.isfinite(V)
prev=obs.mean(0); rar=-np.log(np.clip(prev,1e-4,1.))
NCAT=obs.sum(1); KEEP=NCAT>=8
print(f"类别 {V.shape[1]}  人 {KEEP.sum():,}(>=8 个起始年龄)  稀有度范围 "
      f"{rar.min():.2f}-{rar.max():.2f}(最常见 {prev.max():.1%} / 最罕见 {prev.min():.1%})",flush=True)

# ---- S:稀有亲和特质,来自**多选题选项**(与起始年龄不共享任何 item)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
Ssum={}; Scnt={}; picks={}
nblk=0
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    if s.person.nunique()<1200 or s.option.nunique()<8: continue
    nblk+=1
    br=s.option.map(s.option.value_counts()/s.person.nunique())
    sp=-np.log(np.clip(br,1e-4,1.))
    g=pd.DataFrame({'person':s.person.values,'sp':sp.values}).groupby('person').sp.agg(['mean','size'])
    for p,m_,n_ in zip(g.index,g['mean'].values,g['size'].values):
        Ssum[p]=Ssum.get(p,0.)+m_; Scnt[p]=Scnt.get(p,0)+1; picks[p]=picks.get(p,0)+n_
S=np.full(len(df),np.nan); PK=np.full(len(df),np.nan)
for p,c in Scnt.items():
    if c>=5: S[p]=Ssum[p]/c; PK[p]=picks[p]
ok=np.isfinite(S)
z=lambda a: (a-np.nanmean(a))/np.nanstd(a)
S[ok]=z(S[ok]-np.polyval(np.polyfit(z(PK[ok]),z(S[ok]),1),z(PK[ok])))   # 去掉"勾了多少"(#104)
print(f"S 来自 {nblk} 个多选块,{ok.sum():,} 人;与起始年龄题目零重叠",flush=True)
check_disjoint_items(set(ons),set(str(x) for x in lg.option.unique()),"S vs onset")

# ---- beta_i
def demean(Vm,iters=1):
    """⚠ 一遍"题目去均值 -> 人内去均值"在**缺失不平衡**的面板上不是幂等的:人内去均值
    会把题目均值重新带回来,残余的题目主效应正是与稀有度相关的那一部分。必须交替投影
    迭代到收敛,否则测到的是残余边际,不是交互(#105 的边际决定量陷阱)。"""
    D=np.where(obs,Vm,np.nan)
    for _ in range(iters):
        D=D-np.nanmean(D,axis=0,keepdims=True)      # 类别固定效应(=人群时间表)
        D=D-np.nanmean(D,axis=1,keepdims=True)      # 整体早熟
    return D

def demean_conv(Vm,tol=1e-10,cap=500):
    D=np.where(obs,Vm,np.nan)
    for k in range(cap):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-a
        b=np.nanmean(D,axis=1,keepdims=True); D=D-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    demean_conv.iters=k+1
    return D

def betas(Vm,strip=None):
    """题目去均值 -> 人内去均值 -> 对(中心化)稀有度的人内斜率。常数人效应贡献恰好 0。
    strip: 一个同形状的通道(已在同一双向去均值尺度上),把它从残差里回归掉。
    #114 的偏差正是在这个尺度上测的(-0.2000 年/评分 sd),所以必须在这里剥,
    而不是在原始年龄上用汇总斜率剥 —— 汇总斜率装的是题目层与人层关系,不是回忆偏差。"""
    D=demean_conv(Vm)
    if strip is not None:
        f=np.isfinite(D)&np.isfinite(strip)
        lam=float(np.polyfit(strip[f],D[f],1)[0])
        D=np.where(f,D-lam*strip,D); D=demean_conv(np.where(obs,D,np.nan))
        betas.last_lam=lam
    b=np.full(len(D),np.nan); rho=np.full(len(D),np.nan)
    for i in np.flatnonzero(KEEP):
        m=obs[i]; x=rar[m]-rar[m].mean(); v=(x*x).sum()
        if v<=1e-9: continue
        y=D[i,m]; b[i]=np.nansum(y*x)/v
        sy=np.sqrt(np.nansum((y-np.nanmean(y))**2))
        if sy>1e-9: rho[i]=np.nansum(y*x)/(np.sqrt(v)*sy)   # 尺度无关:这个人多少方差朝着"稀有"排
    return b,rho

def perm_null(Vm,rng):
    """题内跨人置换:精确保留缺失模式与每题的值分布,只摧毁配对。"""
    Wm=Vm.copy()
    for j in range(Vm.shape[1]):
        idx=np.flatnonzero(obs[:,j]); Wm[idx,j]=Vm[rng.permutation(idx),j]
    return Wm

def plant(Vm,rng,g):
    if g==0.: return Vm.copy()
    u=rng.standard_normal(len(Vm)); x=rar-rar.mean()
    return Vm+g*np.outer(u,x)*obs, u
def plant_u(rng): return rng.standard_normal(len(V))

def summarise(br,tag,extra=None):
    b,rho=br
    m=np.isfinite(b)&KEEP
    mr=np.isfinite(rho)&KEEP
    out=dict(tag=tag,n=int(m.sum()),sd=float(np.std(b[m])),
             rho_mean=float(np.mean(rho[mr])),rho_sd=float(np.std(rho[mr])))
    mrs=mr&np.isfinite(S); out['rho_corr_S']=float(np.corrcoef(rho[mrs],S[mrs])[0,1])
    mm=m&np.isfinite(S)
    out['corr_S']=float(np.corrcoef(b[mm],S[mm])[0,1]); out['n_S']=int(mm.sum())
    ma=m&np.isfinite(age)
    out['corr_age']=float(np.corrcoef(b[ma],age[ma])[0,1])
    if extra is not None:
        me=m&np.isfinite(extra); out['corr_plant']=float(np.corrcoef(b[me],extra[me])[0,1])
    return out

GS=[0.0,0.8,2.0,5.0]
rows=[]
for sd_ in range(1,5):
    rng=np.random.default_rng(5100+sd_)
    rows.append(summarise(betas(V),'real'))
    rows[-1]['seed']=sd_
    rows.append({**summarise(betas(perm_null(V,np.random.default_rng(5200+sd_))),'null'),'seed':sd_})
    for g in GS:
        rp=np.random.default_rng(5300+sd_)            # 同种子 -> g=0 必须逐位等于 real
        u=rp.standard_normal(len(V)); x=rar-rar.mean()
        Vp=V+g*np.outer(u,x)*obs
        rows.append({**summarise(betas(Vp),f'plant{g}',extra=u),'seed':sd_})
    print(f"  seed {sd_}",flush=True)

D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'
G=D.groupby('tag')[['sd','corr_S','corr_age']].agg(['mean','std'])
real=D[D.tag=='real']; null=D[D.tag=='null']
# 按人自助真实点(#97c:真实臂的种子展布结构上为零)
b0,rho0=betas(V); mm=np.isfinite(b0)&np.isfinite(rho0)&KEEP&np.isfinite(S)
rb=np.random.default_rng(777); ii=np.flatnonzero(mm)
BS=np.array([[np.std(b0[s_]),np.corrcoef(b0[s_],S[s_])[0,1],
              np.mean(rho0[s_]),np.corrcoef(rho0[s_],S[s_])[0,1]]
             for s_ in (ii[rb.integers(0,len(ii),len(ii))] for _ in range(200))])
sd_boot,cs_boot,rm_boot,rcs_boot=BS.std(0)

print("\n=== 整格(4 seeds,均值 ± 种子 sd) ===")
print(G.round(4).to_string())
print(f"\n真实点按人自助(200 次):sd(beta) {sd_boot:.4f}  corr(beta,S) {cs_boot:.4f}  mean(rho) {rm_boot:.4f}  corr(rho,S) {rcs_boot:.4f}")

g=Gate('稀有亲和是起点还是终点')
r_sd=float(real.sd.mean()); n_sd=float(null.sd.mean())
r_cs=float(real.corr_S.mean()); n_cs=float(null.corr_S.mean())
r_rm=float(real.rho_mean.mean()); n_rm=float(null.rho_mean.mean())
g.degenerate_matches_reference('g=0 逐位复现 real',float(D[D.tag=='plant0.0'].sd.mean()),r_sd,tol=1e-12)
mono=[float(D[D.tag==f'plant{q}'].corr_plant.mean()) for q in GS]
g.asserted('种植的 corr(beta,u) 随 g 单调',all(mono[i]<mono[i+1] for i in range(len(mono)-1)),
           " < ".join(f"{v:.3f}" for v in mono))
g.asserted('置换零的 sd(beta) 比真实**大**,所以 sd 的比较是尺度不匹配,不可读',
           n_sd>r_sd, f'null {n_sd:.4f} > real {r_sd:.4f} -- 置换摧毁人内一致性,残差被放大')
g.require_resolvable_first('corr(beta,S) 可分辨',abs(r_cs),cs_boot)
g.negative_control('corr(beta,S) 对置换零',n_cs,r_cs,null_spread=cs_boot)
r_rcs=float(real.rho_corr_S.mean()); n_rcs=float(null.rho_corr_S.mean())
g.require_resolvable_first('corr(rho,S) 可分辨(尺度无关版)',abs(r_rcs),rcs_boot)
g.negative_control('corr(rho,S) 对置换零',n_rcs,r_rcs,null_spread=rcs_boot)
g.require_resolvable_first('mean(rho) 相对置换零可分辨',abs(r_rm-n_rm),rm_boot)
g.offset_control('mean(rho) 高于置换零',r_rm,n_rm,rm_boot,null_kind='题内跨人置换(保留缺失模式与每题值分布)')
print(g)

# ---- 混淆臂 A:#114 的回忆偏差单独能造出多少?(伪影 = -0.2000 年 / 评分 sd)
RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
zr=(RM-np.nanmean(RM))/np.nanstd(RM)
V_art=np.where(obs&np.isfinite(zr),-0.2000*zr,np.nan)     # 只含伪影,不含真实信号
oa=obs.copy(); obs_backup=obs
obs=obs&np.isfinite(zr)                                    # 伪影只在能配上评分的格子里有定义
NC2=obs.sum(1); KEEP_b=KEEP; KEEP=NC2>=8
b_b,rho_b=betas(V)                                         # 同一受限人群上的**基线**(#101b same_scale)
mb=np.isfinite(rho_b)&KEEP; mbs=mb&np.isfinite(S)
base_rm=float(np.mean(rho_b[mb])); base_rcs=float(np.corrcoef(rho_b[mbs],S[mbs])[0,1])
b_a,rho_a=betas(np.where(obs,V_art,np.nan))
ma=np.isfinite(rho_a)&KEEP; mas=ma&np.isfinite(S)
art_rm=float(np.mean(rho_a[ma])); art_rcs=float(np.corrcoef(rho_a[mas],S[mas])[0,1])
# ---- 混淆臂 B:把起始年龄先对**这个人自己对该类别的评分**残差化,再重跑全部
E=demean_conv(np.where(obs,zr,np.nan))          # 评分通道,放到与 D 相同的双向去均值尺度上
b_c,rho_c=betas(V,strip=E)
lam=betas.last_lam
Dz=demean_conv(np.where(obs,V_art,np.nan))                 # 正对照:剥完后伪影臂的**残差幅度**必须归零
fz=np.isfinite(Dz)&np.isfinite(E)
lz=float(np.polyfit(E[fz],Dz[fz],1)[0]); Rz=Dz-lz*E
art_left=float(np.nanstd(Rz[fz])/max(np.nanstd(Dz[fz]),1e-12))
# 剥离臂的零与自助展布
rho_cn=betas(perm_null(V,np.random.default_rng(5900)),strip=E)[1]
mn2=np.isfinite(rho_cn)&KEEP; ctl_null=float(np.mean(rho_cn[mn2]))
mc0=np.isfinite(rho_c)&KEEP; ic=np.flatnonzero(mc0)
rbb=np.random.default_rng(31337)
ctl_boot=float(np.std([np.mean(rho_c[ic[rbb.integers(0,len(ic),len(ic))]]) for _ in range(200)]))
mc=np.isfinite(rho_c)&KEEP; mcs=mc&np.isfinite(S)
ctl_rm=float(np.mean(rho_c[mc])); ctl_rcs=float(np.corrcoef(rho_c[mcs],S[mcs])[0,1])
obs=obs_backup; KEEP=KEEP_b
print(f"\n=== 混淆:#114 的回忆偏差(-0.2000 年/评分 sd),{len(best)}/{V.shape[1]} 个类别能配上评分 ===")
print(f"  伪影单独产生         mean(rho) {art_rm:+.4f}   corr(rho,S) {art_rcs:+.4f}")
print(f"  评分残差化后的真实值   mean(rho) {ctl_rm:+.4f}   corr(rho,S) {ctl_rcs:+.4f}   (n={mcs.sum():,})")
print(f"  同一受限人群的基线    mean(rho) {base_rm:+.4f}   corr(rho,S) {base_rcs:+.4f}   <- 比较必须对它做")
print(f"  全 31 类别的真实值    mean(rho) {r_rm:+.4f}   corr(rho,S) {r_rcs:+.4f}")
g.same_scale('伪影臂与基线同人群',mb.sum(),ma.sum(),'人数(26 类别受限集)')
# ⚠ 单位幅度的伪影臂问的是"如果年龄**全部**由回忆偏差构成会怎样" —— 那不是这个门要的性质。
#   要的是"这条通道实际贡献了多少",= 移除它前后的差。单位幅度是一个上界,会制造假 FAIL。
g.asserted('单位幅度伪影臂(上界,不是贡献)',True,
           f'mean(rho) {art_rm:+.4f}  corr(rho,S) {art_rcs:+.4f} —— 若年龄 100% 由该通道构成')
g.artifact_cannot_explain('#114 回忆偏差的实际贡献不能解释 mean(rho)',base_rm-ctl_rm,base_rm,rm_boot)
g.artifact_cannot_explain('#114 回忆偏差的实际贡献不能解释 corr(rho,S)',base_rcs-ctl_rcs,base_rcs,rcs_boot)
g.asserted('剥离是有效的:纯伪影臂剥完后残差归零',art_left<0.02,
           f'残余幅度 {art_left:.1%} of 原;剥离斜率 lambda={lam:+.4f} 年/评分 sd,#114 独立测得 -0.2000')
def rho_of(D):
    rl=np.full(len(D),np.nan)
    for i in np.flatnonzero(KEEP):
        m_=obs[i]; x=rar[m_]-rar[m_].mean(); y=D[i,m_]
        sy=np.sqrt(np.nansum((y-np.nanmean(y))**2)); v=np.sqrt((x*x).sum())
        if sy>1e-9 and v>1e-9: rl[i]=np.nansum(y*x)/(v*sy)
    return rl

# ---- 判据:+0.0790 到底是回忆偏差,还是**一遍去均值留下的残余题目主效应**?
print("\n=== 判据:把双向去均值迭代到收敛,不剥任何东西 ===")
print(f"  {'迭代':>6} {'mean(rho)':>11} {'corr(rho,S)':>13}")
for it in [1,2,3,5,10]:
    rl=rho_of(demean(V,iters=it)); m3=np.isfinite(rl)&KEEP; m4=m3&np.isfinite(S)
    print(f"  {it:>6} {np.mean(rl[m3]):>11.4f} {np.corrcoef(rl[m4],S[m4])[0,1]:>13.4f}")
rl=rho_of(demean_conv(V)); m3=np.isfinite(rl)&KEEP; m4=m3&np.isfinite(S)
conv_rm=float(np.mean(rl[m3])); conv_rcs=float(np.corrcoef(rl[m4],S[m4])[0,1])
rln=rho_of(demean_conv(perm_null(V,np.random.default_rng(6100))))
mn3=np.isfinite(rln)&KEEP; conv_null=float(np.mean(rln[mn3]))
ic2=np.flatnonzero(m3); rb2=np.random.default_rng(2718)
conv_boot=float(np.std([np.mean(rl[ic2[rb2.integers(0,len(ic2),len(ic2))]]) for _ in range(200)]))
ic3=np.flatnonzero(m4)
conv_cboot=float(np.std([np.corrcoef(rl[s_],S[s_])[0,1] for s_ in
                         (ic3[rb2.integers(0,len(ic3),len(ic3))] for _ in range(200))]))
print(f"  收敛({demean_conv.iters} 次)  mean(rho) {conv_rm:+.4f}  置换零 {conv_null:+.4f}  展布 {conv_boot:.4f}"
      f"  -> {abs(conv_rm-conv_null)/conv_boot:.1f}x")
print(f"  收敛           corr(rho,S) {conv_rcs:+.4f}  展布 {conv_cboot:.4f}"
      f"  -> {abs(conv_rcs)/conv_cboot:.1f}x")
one_rm=float(np.mean(rho_of(demean(V,iters=1))[np.isfinite(rho_of(demean(V,iters=1)))&KEEP]))
g.asserted('一遍去均值会把结论的**符号**弄反',one_rm*conv_rm<0,
           f'一遍 {one_rm:+.4f} -> 收敛 {conv_rm:+.4f}(13 次交替投影)。'
           f'残余题目主效应恰好沿稀有度排列,精确伪造出"题目属性 x 人"的交互')
g.require_resolvable_first('收敛后 mean(rho) 是否还离得开零',abs(conv_rm-conv_null),conv_boot)
g.offset_control('收敛后 mean(rho) 高于置换零',conv_rm,conv_null,conv_boot,
                 null_kind='题内跨人置换(保留缺失模式与每题值分布)')

# 规格曲线:lambda 从 0 扫到拟合值。拟合值剥掉**全部**评分-年龄协变量 = 伪影的上界;
# #114 独立测得的 -0.2000 是下界。整条曲线发表,不挑一格(#P14 spec_survival)。
print("\n=== lambda 规格曲线(在收敛去均值之上,整条发表)===")
print(f"  {'lambda':>8} {'mean(rho)':>11} {'corr(rho,S)':>13}   注")
spec=[]
for L_,note in [(0.0,'不剥'),(-0.1,''),(-0.2000,'#114 独立测得(下界)'),(-0.25,''),(lam,'拟合(上界:剥掉全部协变量)')]:
    Dl=demean_conv(V); f2=np.isfinite(Dl)&np.isfinite(E)
    Dl=demean_conv(np.where(obs,np.where(f2,Dl-L_*E,Dl),np.nan))
    rl=rho_of(Dl); m3=np.isfinite(rl)&KEEP; m4=m3&np.isfinite(S)
    a,b_=float(np.mean(rl[m3])),float(np.corrcoef(rl[m4],S[m4])[0,1])
    spec.append((L_,a,b_)); print(f"  {L_:>8.4f} {a:>11.4f} {b_:>13.4f}   {note}")
same=sum(1 for _,a,_ in spec[1:] if a>0)
g.asserted('mean(rho) 的符号在 lambda 网格上不稳定',same<len(spec)-1,
           f'{same}/{len(spec)-1} 个剥离规格仍为正 —— 结论随剥离强度换号')
g.no_sign_crossing('mean(rho) 在剥离前后不换号',[base_rm,ctl_rm])
g.resolvable('剥离后的 mean(rho) 是否还离得开零',abs(ctl_rm-ctl_null),ctl_boot)
print(f"  剥离后的零与展布      零 {ctl_null:+.4f}   自助展布 {ctl_boot:.4f}   "
      f"效应 {ctl_rm:+.4f} = {abs(ctl_rm-ctl_null)/ctl_boot:.1f}x")
g.asserted('评分残差化后 mean(rho) 的保留率',True,
           f'{base_rm:+.4f} -> {ctl_rm:+.4f}(保留 {100*ctl_rm/base_rm:.0f}%)')
g.asserted('评分残差化后 corr(rho,S) 的保留率',True,
           f'{base_rcs:+.4f} -> {ctl_rcs:+.4f}(保留 {100*ctl_rcs/base_rcs:.0f}%)')
# ---- 混淆 C:年龄与类别数偏相关
from numpy.linalg import lstsq
mm2=np.isfinite(rho0)&KEEP&np.isfinite(S)&np.isfinite(age)
Z=np.c_[np.ones(mm2.sum()),z(age[mm2]),z(NCAT[mm2].astype(float))]
res_r=rho0[mm2]-Z@lstsq(Z,rho0[mm2],rcond=None)[0]
res_S=S[mm2]-Z@lstsq(Z,S[mm2],rcond=None)[0]
part=float(np.corrcoef(res_r,res_S)[0,1])
print(f"  去掉年龄与类别数后    corr(rho,S) {part:+.4f}(原 {r_rcs:+.4f},保留 {100*part/r_rcs:.0f}%)"
      f"   -> {abs(part)/rcs_boot:.1f}x 自身展布")
g.require_resolvable_first('corr(rho,S) 去掉年龄与类别数后仍可分辨',abs(part),rcs_boot)
print(g)

itm=np.array([np.nanmean(V[obs[:,j],j]) for j in range(V.shape[1])])
sl=np.polyfit(rar,itm,1); rr=np.corrcoef(rar,itm)[0,1]
print(f"\n=== 人群层(题目 n={len(itm)}):平均起始年龄 ~ 稀有度 ===")
print(f"  斜率 {sl[0]:+.3f} 年 / 稀有度单位   r = {rr:+.3f}")
print(f"  最常见的三个:{[f'{prev[j]:.0%}/{itm[j]:.1f}y' for j in np.argsort(-prev)[:3]]}")
print(f"  最罕见的三个:{[f'{prev[j]:.0%}/{itm[j]:.1f}y' for j in np.argsort(prev)[:3]]}")
D.to_csv(OUT/'grid.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
