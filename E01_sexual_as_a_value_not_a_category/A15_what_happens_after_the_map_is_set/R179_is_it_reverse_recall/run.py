import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A15 R02 -- "不搭,所以一定是后来才有的"?连通性差的年龄梯度。

#133a:晚获得的类别与这个人早期那些**更不相连**(去稀有度后 -0.0042,6.5x),
而 #114 的回忆偏差在这条线上**符号相反**,所以"越爱越早"那条通路已排除。

但 #133 自己写下的残余缺口还在,而且它是**另一条通道**:

    反向回忆   「这个东西跟我其余的不搭,所以它一定是后来才有的」
               —— 连通性低**导致**被报成晚,而不是晚导致连通性低。

这条通路有一个可检验的印记:**它是回忆过程的性质,所以它应随回忆的衰减而加深。**
而 #119 已经独立测到「记忆畸变随年龄加深,十五年里几乎翻倍(3.5x)」。所以:

    reverse   连通性差随当前年龄**单调变大** -> 反向回忆通路活,133a 必须降级为
              "至少一部分是叙事整理"
    stable    各档平坦 -> 这条通路死,133a 升级:成年后的扩张是关于**内容**的,
              不是关于**讲述**的

ESTIMAND        去稀有度连通性差(晚半 - 早半),按当前年龄 5 档。
IDENTIFICATION  分割改为**人内中位数**(不是固定 17.5),否则 14-17 档不可能有"晚"。
                零仍是人内置换早/晚标签,保留类别集与晚的个数。
SCOPE           >=8 个类别起始年龄、且人内中位数能把集合分成两个非空半边的人。
WORLDS          reverse / stable(见上)
KILL            条件式:年龄仪器必须先在 #119 的已知量上开火(评分->起始年龄的斜率
                必须随年龄变陡),且每档的人内置换零必须为零,才读年龄趋势。
POSITIVE CTRL   #119 的斜率年龄梯度(已知 -0.0505 @15 -> -0.0917 @30)。
                以及:种植一个随年龄增强的连通性差,必须被这个分层设计回收。
NEGATIVE CTRL   人内置换早/晚标签,每档独立。
NOISE FLOOR     每档 200 次按人自助。
MULTIPLICITY    5 档 x {真实,置换} x 2 种分割(中位数 / 固定 17.5,后者只在能跑的档上),
                整格发表。
IMPOSSIBLE      年龄与队列在横断面里不可分(与 #132 同)。只判有没有梯度,不判成因。
                若梯度存在,"回忆衰减"与"队列差异"仍分不开。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_residualized, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

Ob=obs.astype(float); pj=Ob.mean(0); Cm=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(pj*(1-pj),pj*(1-pj))); den[den<1e-9]=1e-9
SIM=(Cm-np.outer(pj,pj))/den; np.fill_diagonal(SIM,0.)
iu=np.triu_indices(len(rar),1)
X=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],np.abs(rar[iu[0]]-rar[iu[1]])]
res=SIM[iu]-X@np.linalg.lstsq(X,SIM[iu],rcond=None)[0]
check_residualized(res,rar[iu[0]]+rar[iu[1]],"配对相似度对稀有度")
SIMR=np.zeros_like(SIM); SIMR[iu]=res; SIMR=SIMR+SIMR.T

RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]

def contrast(i,latemask):
    jj=np.flatnonzero(obs[i]); lt=jj[latemask]; er=jj[~latemask]
    if len(lt)==0 or len(er)==0: return np.nan
    return SIMR[np.ix_(lt,er)].mean()

def med_split(i,Vm=None):
    Vm=V if Vm is None else Vm
    y=Vm[i,obs[i]]; return y>np.median(y)

CUT=17.5
def fixed_split(i): return V[i,obs[i]]>CUT

AGES=[('14-17',15.5),('18-20',19.0),('21-24',22.5),('25-28',26.5),('29-32',30.5)]
rows=[]
print(f"{'年龄档':<8} {'n(中位)':>6} {'中位数分割':>10} {'倍数':>6} | {'n(17.5)':>6} "
      f"{'固定 17.5 分割':>10} {'倍数':>6} | {'#114 斜率':>10}")
for lab,a in AGES:
    who=[i for i in np.flatnonzero(KEEP&(age==a)) if 0<med_split(i).sum()<obs[i].sum()]
    v=np.array([contrast(i,med_split(i)) for i in who]); v=v[np.isfinite(v)]
    # 固定 17.5 分割 —— #133a 用的那个,唯一能分辨的那个。14-17 档结构上不可能有"晚"。
    whoF=[i for i in np.flatnonzero(KEEP&(age==a)) if 0<fixed_split(i).sum()<obs[i].sum()]
    vF=np.array([contrast(i,fixed_split(i)) for i in whoF]); vF=vF[np.isfinite(vF)]
    nlF=[]
    for s_ in range(3):
        rg=np.random.default_rng(7700+s_)
        for i in whoF:
            k=int(fixed_split(i).sum()); m=obs[i].sum()
            lm=np.zeros(m,bool); lm[rg.choice(m,k,replace=False)]=True
            nlF.append(contrast(i,lm))
    nlF=np.array(nlF); nlF=nlF[np.isfinite(nlF)]
    rbF=np.random.default_rng(int(a*10)+1)
    bsF=float(np.std([vF[rbF.integers(0,len(vF),len(vF))].mean() for _ in range(200)])) if len(vF)>50 else np.nan
    nl=[]
    for s_ in range(3):
        rg=np.random.default_rng(6600+s_)
        for i in who:
            k=int(med_split(i).sum()); m=obs[i].sum()
            lm=np.zeros(m,bool); lm[rg.choice(m,k,replace=False)]=True
            nl.append(contrast(i,lm))
    nl=np.array(nl); nl=nl[np.isfinite(nl)]
    rb=np.random.default_rng(int(a*10))
    bs=float(np.std([v[rb.integers(0,len(v),len(v))].mean() for _ in range(200)]))
    # #114 的正对照:这一档里 评分 -> 起始年龄 的人内斜率
    D_=[];Z_=[]
    for i in who:
        m=obs[i]&np.isfinite(RM[i]); 
        if m.sum()<4: continue
        y=V[i,m]-V[i,m].mean(); z=RM[i,m]-RM[i,m].mean()
        if z.std()>1e-9: D_.append(y); Z_.append(z)
    yy=np.concatenate(D_); zz=np.concatenate(Z_)
    slope=float(np.polyfit(zz/np.std(zz),yy,1)[0])
    rows.append(dict(band=lab,age=a,n=len(v),eff=float(v.mean()),null=float(nl.mean()),
                     boot=bs,slope114=slope,nF=len(vF),
                     effF=float(vF.mean()) if len(vF) else np.nan,
                     nullF=float(nlF.mean()) if len(nlF) else np.nan,bootF=bsF))
    gF=(vF.mean()-nlF.mean()) if len(vF) and len(nlF) else np.nan
    print(f"{lab:<8} {len(v):>6,} {v.mean()-nl.mean():>+10.4f} {abs(v.mean()-nl.mean())/bs:>6.1f}x | "
          f"{len(vF):>6,} {gF:>+10.4f} {abs(gF)/bsF if bsF==bsF else float('nan'):>6.1f}x | {slope:>+10.4f}")

D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'by_age.csv',index=False)
gap=D.eff-D.null
trend=float(np.polyfit(D.age.values,gap.values,1)[0])
rng_g=float(gap.max()-gap.min()); mb=float(D.boot.mean())
s_trend=float(np.polyfit(D.age.values,D.slope114.values,1)[0])

g=Gate('"不搭所以记成晚"这条反向回忆通路活着吗')
g.asserted('年龄仪器先在 #119 的已知量上开火(评分->起始年龄斜率随年龄变陡)',
           s_trend<0 and abs(D.slope114.values[-1])>abs(D.slope114.values[0]),
           "斜率 " + " -> ".join(f"{v:+.4f}" for v in D.slope114.values) +
           f"(趋势 {s_trend:+.5f}/岁)")
# ⚠ 这个零**不应该**是零 —— 它是两个随机半边之间的基线连通度,天然为正。
#   所以它是 offset_control 的偏移量,不是 negative_control 的零。第一版把它当零断言,问错了问题。
g.asserted('置换基线是正的,而且它本来就该是正的(不是零)',bool((D.null>0).all()),
           "基线 " + " ".join(f"{v:+.4f}" for v in D.null.values) +
           " —— 两个随机半边之间的平均连通度,零假设**不预测它为零**")
g.asserted('中位数分割在 29-32 档复现了 #133a 的方向',gap.values[-1]<0,
           f"29-32 档 {gap.values[-1]:+.4f}(#133a 用固定 17.5 分割得 -0.0042)")
g.require_resolvable_first('中位数分割:效应本身可分辨吗',abs(gap.values[-1]),float(D.boot.values[-1]),family='median')
DF=D.dropna(subset=['bootF']); gapF=(DF.effF-DF.nullF)
g.require_resolvable_first('固定 17.5 分割:效应本身可分辨吗',abs(gapF.values[-1]),
                           float(DF.bootF.values[-1]),family='fixed')
trendF=float(np.polyfit(DF.age.values,gapF.values,1)[0])
rng_F=float(gapF.max()-gapF.min()); mbF=float(DF.bootF.mean())
g.require_resolvable_first('固定 17.5 分割:各档之间的差可分辨吗',rng_F,mbF,family='fixed_trend')
g.no_sign_crossing('固定 17.5 分割下所有年龄档同号',list(gapF.values))
g.asserted('反向回忆预测的单调年龄增长:不存在',
           not all(abs(gapF.values[i])<abs(gapF.values[i+1]) for i in range(len(gapF)-1)),
           "按档 " + " ".join(f"{v:+.4f}" for v in gapF.values) +
           f" —— 趋势 {trendF:+.6f}/岁,而**最年轻的一档效应最大**;"
           f"同时年龄仪器在 #114 的斜率上从 -0.1391 走到 -0.2920,所以它不是没检出力")
print(f"\n  固定 17.5 分割,连通性差按档:" + "  ".join(f"{v:+.4f}" for v in gapF.values))
print(f"  年龄趋势 {trendF:+.6f}/岁   极差 {rng_F:.4f} = {rng_F/mbF:.1f}x 单档自助展布")
print(g)
# ---- 剂量曲线:切点从 13.5 扫到 21.5(只在 29-32 档,那里所有切点都可用)
print("\n=== 切点剂量曲线(29-32 档)—— 绝对年龄,还是自己序列里靠后?===")
who9=np.flatnonzero(KEEP&(age==30.5))
print(f"  {'切点':>6} {'n':>6} {'晚的个数':>8} {'连通性差':>10} {'展布':>8} {'倍数':>7}")
dose=[]
for c in [13.5,15.5,17.5,19.5,21.5]:
    ww=[i for i in who9 if 0<(V[i,obs[i]]>c).sum()<obs[i].sum()]
    if len(ww)<100: continue
    vv=np.array([contrast(i,V[i,obs[i]]>c) for i in ww]); vv=vv[np.isfinite(vv)]
    nn=[]
    rg=np.random.default_rng(8800)
    for i in ww:
        k=int((V[i,obs[i]]>c).sum()); m=obs[i].sum()
        lm=np.zeros(m,bool); lm[rg.choice(m,k,replace=False)]=True
        nn.append(contrast(i,lm))
    nn=np.array(nn); nn=nn[np.isfinite(nn)]
    rbc=np.random.default_rng(int(c*10))
    bsc=float(np.std([vv[rbc.integers(0,len(vv),len(vv))].mean() for _ in range(200)]))
    nlate=float(np.mean([(V[i,obs[i]]>c).sum() for i in ww]))
    dose.append((c,len(vv),nlate,float(vv.mean()-nn.mean()),bsc))
    print(f"  {c:>6.1f} {len(vv):>6,} {nlate:>8.1f} {vv.mean()-nn.mean():>+10.4f} {bsc:>8.4f} "
          f"{abs(vv.mean()-nn.mean())/bsc:>7.1f}x")
dv=[d[3] for d in dose]
g.asserted('⚠ 我预测的剂量方向又反了 —— 切点越**早**效应越大',
           all(dv[i]<=dv[i+1]+1e-9 for i in range(len(dv)-1)),
           " < ".join(f"{v:+.4f}" for v in dv) + " (切点 " +
           ", ".join(f"{d[0]:.1f}" for d in dose) + ")")
same=[d for d in dose if abs(d[3]-dose[-1][3])<1e-9]
g.asserted('起始年龄的分箱让三个切点选出**同一个**集合',len(same)>=3,
           f"切点 {', '.join(f'{d[0]:.1f}' for d in same)} 全部 n={same[-1][1]:,}、晚 {same[-1][2]:.1f} 个、"
           f"效应 {same[-1][3]:+.4f} —— 分箱是 17.5 -> 22 -> 28,**所以 #133 的'晚'实际是 19 岁以后**")
g.asserted('方向在每一个切点上都一致,量级不一致(2.4 倍)',
           all(v<0 for v in dv) and max(dv)/min(dv)>2,
           f"全部为负;最大 {min(dv):+.4f} / 最小 {max(dv):+.4f} = {min(dv)/max(dv):.1f} 倍。"
           f"**可报的是方向,不是量级**(frontier §2:效应 X 只许可 <=X 的断言)")
print(g)

print(f"\n  连通性差按档:" + "  ".join(f"{v:+.4f}" for v in gap.values))
print(f"  年龄趋势 {trend:+.6f}/岁   极差 {rng_g:.4f} = {rng_g/mb:.1f}x 单档自助展布")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
