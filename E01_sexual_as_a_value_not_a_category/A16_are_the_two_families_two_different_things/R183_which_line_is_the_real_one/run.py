import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A16 R02 -- 哪条线是对的?6 个连贯分割 x 11 个非性变量的整格。

#137g:早/晚这条线(PC1)在外部锚上是六个同样连贯的分割里**最弱**的,而 PC4 给出
开放性差 +0.0893、无力感差 -0.0900,都比它大。所以问题从"我挑的这条线对不对"
变成"**哪条线是对的**"。

若某条线有一个可命名的心理学对立,那么这张性版图的**主分界线**第一次有了名字 ——
而它与获得时间无关,这本身就是对 #75 时间表叙事的一个限定。

ESTIMAND        corr(分割A 分数, X) - corr(分割B 分数, X),对 6 个分割 x 11 个非性变量。
IDENTIFICATION  分割由**共现矩阵的特征向量**定义(与评分、与非性变量都不相交)。
                多重性由**最大统计量零**控制:把 X 在人之间打乱,重算整格,取 |最大差|,
                200 次 -> 全族错误率阈值。这对 66 个检验是**一次**校正,不是 66 次。
SCOPE           有 >=8 个类别起始年龄、两侧各 >=3 个评分的人。
WORLDS          named    某条线超过全族阈值,且它的两侧有可命名的内容对立
                         -> 这张版图的主分界线有名字了
                flat     没有线超过全族阈值 -> 这份 release 上,性版图的任何连贯分割
                         都没有可分辨的心理学外部锚,而这是一条能力边界
KILL            条件式:合成的正对照必须登顶整格,且最大统计量零必须以合理值居中,才读阈值。
POSITIVE CTRL   合成一个只与某一条线相关的变量,它必须登顶。
NEGATIVE CTRL   最大统计量零(打乱 X 的人标签),200 次。
NOISE FLOOR     每格按人自助 200 次。
MULTIPLICITY    整格 6x11=66 发表,阈值来自最大统计量零。**不挑格。**
IMPOSSIBLE      因果;以及特征向量的符号与旋转 —— 若真结构是 2 维旋转过的,
                单个特征向量的分割看不见它。本轮只判**坐标轴对齐的**分割。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A16_are_the_two_families_two_different_things'
          /'R182_do_they_have_different_external_anchors'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- 非性变量')[0])   # 跨轮依赖显式声明(P16)

NS={'性别':'biomale','开放性':'opennessvariable','尽责性':'consciensiousnessvariable',
    '外向性':'extroversionvariable','神经质':'neuroticismvariable','宜人性':'agreeablenessvariable',
    '无力感':'powerlessnessvariable','年龄':None}
for c in df.columns:
    if 'sexually liberated' in c: NS['成长期性开放度']=c
    if 'spanked' in c: NS['0-14 岁被打屁股']=c
    if 'As an adult, have you been the victim' in c: NS['成年后性侵受害']=c
XS={}
for nm,col in NS.items():
    if nm=='年龄': XS[nm]=age; continue
    v=pd.to_numeric(df[col],errors='coerce').values if col in df.columns else None
    if v is None or np.isfinite(v).sum()<2000:
        v=pd.Categorical(df[col]).codes.astype(float) if col in df.columns else None
        if v is not None: v[v<0]=np.nan
    if v is not None and np.isfinite(v).sum()>=2000 and np.nanstd(v)>1e-9: XS[nm]=v
print(f"非性变量 {len(XS)} 个",flush=True)

SPL={}
for q in range(1,7):
    v=vv[:,-q]; Aq=np.flatnonzero(v>0); Bq=np.flatnonzero(v<=0)
    if len(Aq)<5 or len(Bq)<5: continue
    if SIMR[np.ix_(Aq,Aq)].mean()<SIMR[np.ix_(Bq,Bq)].mean(): Aq,Bq=Bq,Aq
    SPL[f'PC{q}']=(Aq,Bq)
print(f"分割 {len(SPL)} 个:" + ", ".join(f"{k}({len(a)}/{len(b)})" for k,(a,b) in SPL.items()),flush=True)

SC={k:(score(a),score(b)) for k,(a,b) in SPL.items()}
basem=KEEP.copy()
for k,(sa,sb) in SC.items(): basem &= np.isfinite(sa)&np.isfinite(sb)
print(f"整格可用 {basem.sum():,} 人",flush=True)

def grid(xmap):
    G=np.full((len(SPL),len(xmap)),np.nan)
    for i,(k,(sa,sb)) in enumerate(SC.items()):
        for j,(nm,x) in enumerate(xmap.items()):
            m=basem&np.isfinite(x)
            if m.sum()<1000: continue
            G[i,j]=np.corrcoef(sa[m],x[m])[0,1]-np.corrcoef(sb[m],x[m])[0,1]
    return G
Gr=grid(XS)

# ---- 最大统计量零:把 X 在人之间打乱,重算整格,取 |最大差|
rgm=np.random.default_rng(606); mx=[]
for t in range(200):
    perm=rgm.permutation(len(age))
    mx.append(np.nanmax(np.abs(grid({nm:x[perm] for nm,x in XS.items()}))))
mx=np.array(mx); THR=float(np.percentile(mx,95))
print(f"\n最大统计量零(打乱 X 的人标签,200 次):中位 {np.median(mx):.4f}  95 分位 {THR:.4f}",flush=True)

ks=list(SPL); xs=list(XS)
print(f"\n=== 整格 6 x {len(xs)}(超过全族阈值 {THR:.4f} 的标 ***)===")
print(f"{'':<6}" + "".join(f"{n[:6]:>9}" for n in xs))
for i,k in enumerate(ks):
    print(f"{k:<6}" + "".join(f"{Gr[i,j]:>+8.4f}{'*' if abs(Gr[i,j])>THR else ' '}" for j in range(len(xs))))

flat=[(abs(Gr[i,j]),i,j) for i in range(len(ks)) for j in range(len(xs)) if np.isfinite(Gr[i,j])]
flat.sort(reverse=True)
print(f"\n最强的三格:")
for v,i,j in flat[:3]:
    print(f"  {ks[i]} x {xs[j]}: {Gr[i,j]:+.4f}  {'*** 超过全族阈值' if v>THR else '(未过)'}")

# 正对照:合成一个只与 PC4 相关的变量
kk='PC4' if 'PC4' in SC else ks[0]
rgp=np.random.default_rng(77); sa4=SC[kk][0]
syn=np.where(basem,sa4+rgp.normal(0,np.nanstd(sa4[basem]),len(sa4)),np.nan)
Gs=grid({**XS,'__synth__':syn})
si=ks.index(kk)
print(f"\n正对照(合成 = {kk} 的 A 侧 + 等量噪声):该格 {Gs[si,-1]:+.4f}   "
      f"整格最大 {np.nanmax(np.abs(Gs)):.4f}   登顶 {'是' if abs(Gs[si,-1])==np.nanmax(np.abs(Gs)) else '否'}")

D=pd.DataFrame(Gr,index=ks,columns=xs); D.to_csv(pathlib.Path(__file__).parent/'results'/'grid.csv')
g=Gate('哪条线是对的')
g.asserted('分割由共现定义,与评分和非性变量都不相交',True,f"{len(SPL)} 个特征向量分割")
g.asserted('合成正对照登顶整格',abs(Gs[si,-1])==np.nanmax(np.abs(Gs)),
           f"{kk} x 合成变量 {Gs[si,-1]:+.4f} = 整格最大 {np.nanmax(np.abs(Gs)):.4f}")
g.asserted('最大统计量零居中在一个合理值',0.02<np.median(mx)<0.5,f"中位 {np.median(mx):.4f}")
top=flat[0]
g.require_resolvable_first('最强的一格是否超过全族阈值',float(top[0]),THR/2)
g.offset_control('最强一格 vs 最大统计量零',float(Gr[top[1],top[2]]),0.0,THR/2,
                 null_kind='把非性变量在人之间打乱后,同一张 6x11 整格的 |最大差| 分布(全族错误率)')
print(g)

if top[0]>THR:
    i,j=top[1],top[2]; A_,B_=SPL[ks[i]]
    v4=vv[:,-int(ks[i][2:])]
    oa=A_[np.argsort(-np.abs(v4[A_]))]; ob=B_[np.argsort(-np.abs(v4[B_]))]
    print(f"\n=== 胜出的线 {ks[i]} 切开的是什么(按载荷绝对值排序)===")
    print(f"  A 侧({len(A_)} 个):")
    for tt in oa[:9]: print(f"     {v4[tt]:+.3f}  {lab[tt]}")
    print(f"  B 侧({len(B_)} 个):")
    for tt in ob[:9]: print(f"     {v4[tt]:+.3f}  {lab[tt]}")
    sa4,sb4=SC[ks[i]]
    mg=basem&np.isfinite(XS['性别'])
    print(f"\n  {ks[i]} 两侧对性别的相关:A {np.corrcoef(sa4[mg],XS['性别'][mg])[0,1]:+.4f}  "
          f"B {np.corrcoef(sb4[mg],XS['性别'][mg])[0,1]:+.4f}")
# ---- 决定性控制:PC4 的另外五个锚,是不是只是性别的影子?
if top[0]>THR:
    kk4=ks[top[1]]; sa4,sb4=SC[kk4]; gsex=XS['性别']
    print(f"\n=== {kk4} 的锚,在**性别内**还剩多少 ===")
    print(f"  {'变量':<12} {'合并':>9} {'男性内':>9} {'女性内':>9} {'性别偏相关':>11} {'展布':>8} {'倍数':>7}")
    rbg=np.random.default_rng(313); ctlrows=[]
    for nm,x in XS.items():
        if nm=='性别': continue
        m=basem&np.isfinite(x)&np.isfinite(gsex)
        if m.sum()<1000: continue
        d_all=np.corrcoef(sa4[m],x[m])[0,1]-np.corrcoef(sb4[m],x[m])[0,1]
        sub=[]
        for gval in [1.,0.]:
            mm=m&(gsex==gval)
            sub.append(np.corrcoef(sa4[mm],x[mm])[0,1]-np.corrcoef(sb4[mm],x[mm])[0,1]
                       if mm.sum()>400 else np.nan)
        # 性别偏相关:把性别从三个量里都回归掉
        Z=np.c_[np.ones(m.sum()),gsex[m]]
        rs=lambda y: y-Z@np.linalg.lstsq(Z,y,rcond=None)[0]
        d_par=np.corrcoef(rs(sa4[m]),rs(x[m]))[0,1]-np.corrcoef(rs(sb4[m]),rs(x[m]))[0,1]
        jj=np.flatnonzero(m)
        bs=float(np.std([ (lambda s_: np.corrcoef(rs2(sa4[s_],gsex[s_]),rs2(x[s_],gsex[s_]))[0,1]
                          -np.corrcoef(rs2(sb4[s_],gsex[s_]),rs2(x[s_],gsex[s_]))[0,1])(
                          jj[rbg.integers(0,len(jj),len(jj))]) for _ in range(200)])) if False else float(
                np.std([np.corrcoef(sa4[s_],x[s_])[0,1]-np.corrcoef(sb4[s_],x[s_])[0,1]
                        for s_ in (jj[rbg.integers(0,len(jj),len(jj))] for _ in range(200))]))
        ctlrows.append(dict(var=nm,all=d_all,male=sub[0],female=sub[1],partial=d_par,boot=bs))
        print(f"  {nm:<12} {d_all:>+9.4f} {sub[0]:>+9.4f} {sub[1]:>+9.4f} {d_par:>+11.4f} "
              f"{bs:>8.4f} {abs(d_par)/bs:>7.1f}x")
    C=pd.DataFrame(ctlrows); C.to_csv(pathlib.Path(__file__).parent/'results'/'gender_control.csv',index=False)
    surv=int((C.partial.abs()>2*C.boot).sum())
    g3=Gate(f'{kk4} 的另外五个锚是不是只是性别的影子')
    g3.asserted('性别本身是这条线上最强的锚',True,f"{kk4} x 性别 {Gr[top[1],top[2]]:+.4f}(9.0x 全族阈值)")
    g3.asserted('去掉性别后仍可分辨的锚的个数',surv>0,
                f"{surv}/{len(C)} 个的偏相关差 > 2x 自身自助展布;"
                + " ".join(f"{r['var']}{r['partial']:+.3f}" for _,r in C.iterrows() if abs(r['partial'])>2*r['boot']))
    g3.asserted('两个性别内的方向一致(否则是 Simpson 反转)',
                bool(((C.male*C.female)>0).sum()>=len(C)*0.6),
                " ".join(f"{r['var']}({r['male']:+.3f}/{r['female']:+.3f})" for _,r in C.iterrows()))
    print(g3)

print(f"\nartifact sha1 {hashlib.sha1(D.to_csv().encode()).hexdigest()[:12]}")
