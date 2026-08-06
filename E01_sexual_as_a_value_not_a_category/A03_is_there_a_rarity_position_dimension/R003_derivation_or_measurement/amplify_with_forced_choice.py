import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A19 R03 -- 用强制选择放大效应,把接缝的等价边界压到有意义的水平

#150b:判别量 `corr(z,S)+corr(ρ,S)` = −0.0034(0.2×),但它的 95% 上界是 0.0317,
而效应本身只有 0.0370 —— **只能排除大于效应 86% 的差异**。CONSISTENT,不是 PROVEN。
而限制**不是 n**(已近万,两个相关的展布都在 0.010),是**效应小**。

#150 的 NEXT:用 `#126` 的**强制选择**块重建 S —— 它在构造上移除作答水平。

⚠ **但强制选择的 10 个块里有 8 个就是起始年龄的类别**(束缚 · 温柔 · 非自愿 · 权力动态 ·
   怀孕 · 玩具 · 羞辱 · 精神改变)。若不处理,S_fc 与 z/ρ 共享内容,相关会被共享 item 抬高
   —— 那正是 `#126c` 在设计时漏掉、`#127` 回头补的那个洞。
   **所以 z 与 ρ 必须在剔除那 8 个类别后的 23 个上重算(留块法),`check_disjoint_items` 断言。**

ESTIMAND        corr(z_resid, S_fc) 与 corr(ρ, S_fc),z/ρ 在留块后的类别上算,
                按类别数卡钳匹配;判别量 = 两者之和(若"同一件事"成立应 ≈ 0)。
IDENTIFICATION  S_fc 来自强制选择(每人每块只能选一个,作答水平在构造上被移除);
                z/ρ 来自剩下 23 个类别的起始年龄。两者**零重叠**,断言。
SCOPE           >=6 个留块后类别的起始年龄、且 >=5 个强制选择块有作答的人。
WORLDS          AMPLIFIED  |corr| 明显大于 0.037,等价边界随之收紧 -> 接缝可以焊死
                SAME_SIZE  |corr| 与 0.037 同量级 -> 强制选择没有放大它,
                           而"作答水平"本来就不是这条相关的主要成分
                KILLED     |corr| 归零 -> 之前的 0.037 里有共享 item 的成分,#150 要降级
KILL            条件式:留块必须真的把重叠清零(断言),匹配必须把类别数差压到 <0.1 sd,
                且负对照(按人置换 S_fc)必须归零,才读效应。
POSITIVE CTRL   种植一个与 S_fc 相关的人特异径向信号,两个相关必须被拉动且符号相反。
NEGATIVE CTRL   按人置换 S_fc 的标签。
NOISE FLOOR     按人自助 200 次;5 个匹配种子。
MULTIPLICITY    2 个量 x {多选 S, 强制选择 S_fc} x {未匹配, 匹配},整格发表。
IMPOSSIBLE      强制选择只有 10 块,所以 S_fc 的信度天然低于多选 S;
                若它没放大,分不清是"作答水平不重要"还是"S_fc 太吵"。噪声地板报出。
"""
import numpy as np, pandas as pd, warnings, hashlib, re, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_disjoint_items

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

FC=inv[inv['kind']=='FORCED_CHOICE_MOST']['col'].tolist()
print(f"强制选择块 {len(FC)}:{FC}",flush=True)
# 把 FC 块名映射到起始年龄类别
KEYS={'bondage':'bondage','gentleness':'gentleness','nonconsent':'nonconsent',
      'powerdynamic':'power dynamics','pregnancy':'pregnancy','toys':'toys',
      'humiliation':'humiliation','mentalalteration':'mental alteration'}
drop=set()
for c in FC:
    for k,v in KEYS.items():
        if k in c.lower():
            for j,o in enumerate(ons):
                if v.split()[0] in o.lower(): drop.add(j)
KEEPC=np.array([j for j in range(V.shape[1]) if j not in drop])
print(f"剔除 {len(drop)} 个与强制选择重叠的类别,剩 {len(KEEPC)} 个",flush=True)
check_disjoint_items(set(int(j) for j in drop),set(int(j) for j in KEEPC),"留块后 z/ρ vs 强制选择块")

# S_fc:每人每块所选项的 -log(基率),跨块平均(作答水平在构造上被移除)
sm=np.zeros(len(df)); sc=np.zeros(len(df))
for c in FC:
    s=df[c]; m=s.notna()
    v=s[m]; pop=v.map(v.value_counts()/len(v)).astype(float).values
    sm[np.flatnonzero(m)]+=-np.log(np.clip(pop,1e-4,1.)); sc[np.flatnonzero(m)]+=1
Sfc=np.full(len(df),np.nan); okf=sc>=5; Sfc[okf]=sm[okf]/sc[okf]
Sfc[okf]=(Sfc[okf]-Sfc[okf].mean())/Sfc[okf].std()
print(f"S_fc 可算的人 {int(okf.sum()):,};与多选 S 的相关 "
      f"{np.corrcoef(Sfc[okf&np.isfinite(S)],S[okf&np.isfinite(S)])[0,1]:+.4f}",flush=True)

obs2=obs.copy(); obs2[:,list(drop)]=False
def demean_conv(Vm,ob,tol=1e-10,cap=500):
    Dm=np.where(ob,Vm,np.nan)
    for _ in range(cap):
        a=np.nanmean(Dm,axis=0,keepdims=True); Dm=Dm-a
        b=np.nanmean(Dm,axis=1,keepdims=True); Dm=Dm-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return Dm
NPERM=200
def two_stats(Dres,ob,seed):
    rg=np.random.default_rng(seed); tie=np.random.default_rng(20260803)
    Z=np.full(len(Dres),np.nan); Rho=np.full(len(Dres),np.nan)
    for i in range(len(Dres)):
        j=np.flatnonzero(ob[i]); k=len(j)
        if k<6: continue
        y=Dres[i,j]; r=rar[j]
        cand=np.flatnonzero(y==np.nanmin(y)); pick=cand[tie.integers(len(cand))]
        d=r[pick]-r.mean()
        idx=rg.integers(0,k,(NPERM,1)); dr=r[idx].mean(1)-r.mean()
        if dr.std()<1e-9 or np.nanstd(y)<1e-9: continue
        Z[i]=(d-dr.mean())/dr.std(); Rho[i]=np.corrcoef(y,r)[0,1]
    return Z,Rho

Dres=demean_conv(V,obs2)
Z,Rho=two_stats(Dres,obs2,zlib.crc32(b'A19R03'))
NC2=obs2.sum(1).astype(float)
def analyse(Svec,tag):
    base=np.isfinite(Z)&np.isfinite(Rho)&np.isfinite(Svec)
    ii=np.flatnonzero(base); rb=np.random.default_rng(5)
    out=[]
    for nm,vec in [('z_resid',Z),('rho',Rho)]:
        raw=float(np.corrcoef(vec[ii],Svec[ii])[0,1])
        mv=[]
        for sd_ in range(5):
            rg=np.random.default_rng(600+sd_); med=np.median(Svec[ii])
            hi=ii[Svec[ii]>med]; lo=ii[Svec[ii]<=med]
            c=(NC2-NC2[ii].mean())/NC2[ii].std(); used=np.zeros(len(Svec),bool); P=[]
            for a in hi[rg.permutation(len(hi))]:
                d_=np.abs(c[lo]-c[a]); d_[used[lo]]=np.inf; kk=int(np.argmin(d_))
                if d_[kk]<0.25: used[lo[kk]]=True; P.append((a,lo[kk]))
            P=np.array(P); sel=np.r_[P[:,0],P[:,1]]
            mv.append((float(np.corrcoef(vec[sel],Svec[sel])[0,1]),
                       abs(NC2[P[:,0]].mean()-NC2[P[:,1]].mean())/NC2[ii].std()))
        bo=float(np.std([(lambda s_: np.corrcoef(vec[s_],Svec[s_])[0,1])(
            ii[rb.integers(0,len(ii),len(ii))]) for _ in range(200)]))
        out.append(dict(S=tag,stat=nm,n=len(ii),raw=raw,matched=float(np.mean([x[0] for x in mv])),
                        bal=float(np.mean([x[1] for x in mv])),boot=bo))
    return out

rows=analyse(S,'多选 S')+analyse(Sfc,'强制选择 S_fc')
T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'amplify.csv',index=False)
print(f"\n{'S 的来源':<14}{'量':<9}{'n':>7}{'未匹配':>10}{'匹配后':>10}{'类别数残差':>10}{'展布':>9}")
for _,r in T.iterrows():
    print(f"{r.S:<14}{r.stat:<9}{int(r.n):>7,}{r.raw:>+10.4f}{r.matched:>+10.4f}{r.bal:>10.3f}{r.boot:>9.4f}")

g=Gate('强制选择能不能把接缝的等价边界压下去')
g.asserted('留块把重叠清零',len(set(KEEPC.tolist())&drop)==0,f"剔除 {len(drop)} 个,剩 {len(KEEPC)} 个")
for tag in ['多选 S','强制选择 S_fc']:
    cz=T[(T.S==tag)&(T.stat=='z_resid')].iloc[0]; cr=T[(T.S==tag)&(T.stat=='rho')].iloc[0]
    disc=float(cz.matched+cr.matched); sdv=float(np.sqrt(cz.boot**2+cr.boot**2))
    marg=0.5*abs(float(cr.matched))
    print(f"\n  【{tag}】效应 |corr(ρ,S)| = {abs(cr.matched):.4f};判别量 {disc:+.4f} ± {sdv:.4f};"
          f"95% 上界 {abs(disc)+2*sdv:.4f};只能排除 > 效应 {100*(abs(disc)+2*sdv)/abs(cz.matched):.0f}% 的差异")
    g.asserted(f'{tag}:匹配把类别数差压下去了',float(cz.bal)<0.1,f"{cz.bal:.3f} sd")
    g.equivalent_within(f'{tag}:判别量落在等价边界内(效应一半)',disc,sdv,marg)
fcz=T[(T.S=='强制选择 S_fc')&(T.stat=='rho')].iloc[0]
msz=T[(T.S=='多选 S')&(T.stat=='rho')].iloc[0]
# ---- 分歧义:是"作答水平不重要",还是"S_fc 太吵"?独立量 S_fc 与 S 的信度。
def sfc_half(cols):
    sm2=np.zeros(len(df)); sc2=np.zeros(len(df))
    for c in cols:
        s=df[c]; m=s.notna(); v=s[m]
        pop=v.map(v.value_counts()/len(v)).astype(float).values
        sm2[np.flatnonzero(m)]+=-np.log(np.clip(pop,1e-4,1.)); sc2[np.flatnonzero(m)]+=1
    o=sc2>=len(cols)//2; out=np.full(len(df),np.nan); out[o]=sm2[o]/sc2[o]; return out
rgh=np.random.default_rng(2026); pm=rgh.permutation(len(FC))
h1,h2=[FC[i] for i in pm[:5]],[FC[i] for i in pm[5:]]
A,B=sfc_half(h1),sfc_half(h2); mh=np.isfinite(A)&np.isfinite(B)
r_half=float(np.corrcoef(A[mh],B[mh])[0,1]); rel_fc=2*r_half/(1+r_half)
REL_S=0.4611                                        # #100/#127b 的跨不相交块分半
pred=abs(float(msz.matched))*np.sqrt(max(rel_fc,1e-9)/REL_S)
print("\n=== 分歧义:S_fc 是更干净,还是更吵?===")
print(f"  S_fc 分半(5+5 块)r = {r_half:+.4f} -> Spearman-Brown 信度 **{rel_fc:.4f}**"
      f"(而多选 S 的是 {REL_S:.4f},#100/#127b)")
print(f"  衰减预测:|corr(ρ,S_fc)| ≈ |corr(ρ,S)| x sqrt(rel_fc/rel_S) = "
      f"{abs(msz.matched):.4f} x {np.sqrt(rel_fc/REL_S):.3f} = **{pred:.4f}**")
print(f"  实测 |corr(ρ,S_fc)| = **{abs(fcz.matched):.4f}**   差 {abs(fcz.matched)-pred:+.4f}")
g.asserted('S_fc 的信度独立量出来了',np.isfinite(rel_fc),
           f"分半 r={r_half:+.4f} -> SB {rel_fc:.4f} vs 多选 S 的 {REL_S:.4f}")
g.equivalent_within('缩小量是否就是衰减(实测 vs 衰减预测)',
                    abs(fcz.matched)-pred,float(fcz.boot),0.5*abs(msz.matched))
g.asserted('强制选择是否放大了效应',abs(fcz.matched)>abs(msz.matched),
           f"多选 {abs(msz.matched):.4f} -> 强制选择 {abs(fcz.matched):.4f}"
           f"({'放大' if abs(fcz.matched)>abs(msz.matched) else '没放大'} "
           f"{abs(fcz.matched)/abs(msz.matched):.2f}×)")
print(g)
print(f"\nartifact sha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
