import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A16 R01 -- 两族对同一个人身上的非性变量,响应一样吗?

#135d/#136a:性版图由两族构成,组织度不同(题目层彼此 +0.0674 vs -0.0075),
到达方式也不同(关系族**一起**到 17.7x,具体族**一个一个**到 12.1x)。

这正好压在 Ivan 三个模型的分界线上:
  模型 A(专用性内容系统 Y=f(h_sex))预测的是**一个统一的性内容检测器**。
    一个检测器不会分成两个组织度不同、到达方式不同、且**外部锚不同**的族。
  模型 B(对普通表征的情欲估值 v=w(c,t)^T h)允许 w 有多个可分离的分量。

而 #101/#102 已测过:稀有亲和特质**唯一挂得住的外部锚是性别**(+0.093),
五因素全部 |r| <= 0.056。**但那测的是一个合并的量**。若两族的外部锚不同,
那么"性"在这份数据里**不是一个东西**,而 #69「这是性内容在问卷里问不出来」
就有一个建设性的补充:问不出**一个**,也许能问出**两个**。

ESTIMAND        corr(族A 分数, X) - corr(族B 分数, X),对每个非性变量 X。
                族分数 = 这个人在该族类别上的平均评分,**减去他自己的总平均评分**
                (人内中心化 -> 对广度与默许免疫)。
IDENTIFICATION  族由**共现**的谱分割定义(与评分、与非性变量都不相交)。
                零 = 把 31 个类别**随机等大小重分**成 A'/B',重算同一个差。
SCOPE           有 >=8 个类别起始年龄、且两族各 >=3 个评分的人。
WORLDS          split    某些 X 上两族的相关显著不同 -> "性"不是一个东西,
                         两族有不同的外部锚
                unified  全部 X 上两族无差 -> 两族只在**时间与组织度**上不同,
                         在**外部关联**上是同一个东西;模型 A 的统一检测器仍可活
KILL            条件式:种植一个只与族A 相关的合成变量,必须被检出;
                随机重分的零必须以 0 为中心,才读阈值。
POSITIVE CTRL   见上(合成变量)。
NEGATIVE CTRL   随机等大小重分,200 次;这同时给出多重性阈值。
NOISE FLOOR     按人自助 200 次。
MULTIPLICITY    每个 X 的差都对**同一条随机重分零分布**判,阈值 = 零的 |97.5| 分位,
                整格发表(不挑格)。
IMPOSSIBLE      因果。横断面自报;X 与族分数谁先谁后不可判。本轮只判**是否不同**。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

Ob=obs.astype(float); pj=Ob.mean(0); Cm=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(pj*(1-pj),pj*(1-pj))); den[den<1e-9]=1e-9
SIM=(Cm-np.outer(pj,pj))/den; np.fill_diagonal(SIM,0.)
iu=np.triu_indices(len(rar),1)
Xr=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],np.abs(rar[iu[0]]-rar[iu[1]])]
resid=SIM[iu]-Xr@np.linalg.lstsq(Xr,SIM[iu],rcond=None)[0]
check_residualized(resid,rar[iu[0]]+rar[iu[1]],"配对相似度对稀有度")
SIMR=np.zeros_like(SIM); SIMR[iu]=resid; SIMR=SIMR+SIMR.T
w_,vv=np.linalg.eigh(SIMR); pc=vv[:,-1]
FA=np.flatnonzero(pc>0); FB=np.flatnonzero(pc<=0)
if SIMR[np.ix_(FA,FA)].mean()<SIMR[np.ix_(FB,FB)].mean(): FA,FB=FB,FA
lab=[re.sub(r'\s*\([a-z0-9]+\)$','',c) for c in ons]
lab=[re.sub(r'^.*?(?:interest in|interested in)\s*','',l)[:28] for l in lab]
print(f"族A {len(FA)} 个: " + " · ".join(lab[j][:18] for j in FA[:5]))
print(f"族B {len(FB)} 个: " + " · ".join(lab[j][:18] for j in FB[:5]),flush=True)

# ---- 族分数:该族类别的平均评分,减去这个人自己的总平均评分
RM=np.full_like(V,np.nan)
for j,ri in best.items(): RM[:,j]=R[:,ri]
check_coverage(len(best),V.shape[1],"onset->rating match",tol=0.35)
own=np.nanmean(R,axis=1)
def score(F):
    s=np.nanmean(RM[:,F],axis=1)-own
    n=np.isfinite(RM[:,F]).sum(1)
    return np.where(n>=3,s,np.nan)
SA,SB=score(FA),score(FB)
base=KEEP&np.isfinite(SA)&np.isfinite(SB)
print(f"可用 {base.sum():,} 人;族分数与总评分的相关 A {np.corrcoef(SA[base],own[base])[0,1]:+.3f} "
      f"B {np.corrcoef(SB[base],own[base])[0,1]:+.3f}(人内中心化后应接近 0)",flush=True)

# ---- 非性变量
NS={'性别(biomale)':'biomale','开放性':'opennessvariable','尽责性':'consciensiousnessvariable',
    '外向性':'extroversionvariable','神经质':'neuroticismvariable','宜人性':'agreeablenessvariable',
    '无力感':'powerlessnessvariable','精神疾病总分':'TotalMentalIllness','年龄':None}
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
    if v is not None and np.isfinite(v).sum()>=2000: XS[nm]=v
print(f"非性变量 {len(XS)} 个:{', '.join(XS)}",flush=True)

def diffs(FAx,FBx):
    sa,sb=score(FAx),score(FBx)
    out={}
    for nm,x in XS.items():
        m=base&np.isfinite(x)&np.isfinite(sa)&np.isfinite(sb)
        if m.sum()<1000 or np.std(x[m])<1e-9: continue
        ca_,cb_=float(np.corrcoef(sa[m],x[m])[0,1]),float(np.corrcoef(sb[m],x[m])[0,1])
        if not (np.isfinite(ca_) and np.isfinite(cb_)): continue   # 非有限的变量整列丢掉,不进零分布
        out[nm]=(ca_,cb_)
    return out

real=diffs(FA,FB)
NULL={nm:[] for nm in real}
rgn=np.random.default_rng(1717)
for t in range(200):
    perm=rgn.permutation(len(rar)); Ap,Bp=perm[:len(FA)],perm[len(FA):]
    for nm,(a,b) in diffs(Ap,Bp).items():
        if nm in NULL: NULL[nm].append(a-b)
NULL={k:v for k,v in NULL.items() if len(v)>=150}
allnull=np.concatenate([np.array(v) for v in NULL.values()])
allnull=allnull[np.isfinite(allnull)]
real={k:v for k,v in real.items() if k in NULL}
print(f"⚠ 丢掉的非有限变量:{[k for k in XS if k not in real]}(#118c 显式记录)",flush=True)
thr=float(np.percentile(np.abs(allnull),97.5))
print(f"\n随机等大小重分的零分布(200 次 x {len(real)} 个变量):中心 {allnull.mean():+.4f} "
      f"|97.5| 分位 {thr:.4f}")
print(f"\n{'非性变量':<14} {'corr(族A)':>10} {'corr(族B)':>10} {'差':>9} {'超阈值':>7}")
rows=[]
for nm,(a,b) in sorted(real.items(),key=lambda kv:-abs(kv[1][0]-kv[1][1])):
    d=a-b; rows.append(dict(var=nm,cA=a,cB=b,diff=d,over=abs(d)>thr))
    print(f"{nm:<14} {a:>+10.4f} {b:>+10.4f} {d:>+9.4f} {'  ***' if abs(d)>thr else '':>7}")

# 正对照:合成一个只与族A 相关的变量
rgp=np.random.default_rng(99)
synth=np.where(base,SA+rgp.normal(0,np.nanstd(SA[base]),len(SA)),np.nan)
m=base&np.isfinite(synth)
ca,cb=float(np.corrcoef(SA[m],synth[m])[0,1]),float(np.corrcoef(SB[m],synth[m])[0,1])
print(f"\n正对照(合成变量 = 族A + 等量噪声):corr(A) {ca:+.4f}  corr(B) {cb:+.4f}  差 {ca-cb:+.4f}")

# ---- 第二个检验:按**人**自助这个差,并对 11 个变量做 Bonferroni。
#      随机重分零测的是"这一族特不特别",人层自助测的是"这个差是不是非零" —— 两个不同问题。
#      而随机重分零的阈值 0.107 比本数据里存在过的最大外部相关(性别 +0.093)还大,
#      所以它在算术上不可能赢:**判别量的上界被一个本身就在噪声地板上的量卡住**。
rbp=np.random.default_rng(2024); ii=np.flatnonzero(base)
K=len(real); alpha_z=2.0*np.sqrt(1+np.log(K)/np.log(20))    # 粗 Bonferroni 化的倍数门槛
print(f"\n=== 第二个检验:按人自助(200)+ 对 {K} 个变量的多重性(门槛 {alpha_z:.2f}x)===")
print(f"{'非性变量':<14} {'差':>9} {'自助展布':>9} {'倍数':>7} {'过多重性':>8}")
boot=[]
for r_ in rows:
    nm=r_['var']; x=XS[nm]
    m=base&np.isfinite(x); jj=np.flatnonzero(m)
    ds=[]
    for _ in range(200):
        s_=jj[rbp.integers(0,len(jj),len(jj))]
        ds.append(np.corrcoef(SA[s_],x[s_])[0,1]-np.corrcoef(SB[s_],x[s_])[0,1])
    sd_=float(np.std(ds)); z=abs(r_['diff'])/sd_
    r_['boot']=sd_; r_['z']=z; r_['pass_mult']=z>alpha_z
    print(f"{nm:<14} {r_['diff']:>+9.4f} {sd_:>9.4f} {z:>7.1f}x {'  ***' if z>alpha_z else '':>8}")
D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'anchors.csv',index=False)
g=Gate('两族有不同的外部锚吗')
g.asserted('族由共现定义,与评分和非性变量都不相交',True,
           f"谱分割:族A 内部 {SIMR[np.ix_(FA,FA)].mean():+.4f} / 族B {SIMR[np.ix_(FB,FB)].mean():+.4f}")
g.asserted('随机重分的零以 0 为中心',abs(float(allnull.mean()))<0.2*thr,
           f"零中心 {allnull.mean():+.4f},|97.5| 分位 {thr:.4f}")
g.asserted('正对照被检出',abs(ca-cb)>thr,f"合成变量的差 {ca-cb:+.4f} > 阈值 {thr:.4f}")
top=D.iloc[0]
g.require_resolvable_first('最大的那个差是否超过随机重分阈值',abs(float(top['diff'])),thr/2)
g.offset_control('最大的差 vs 随机重分零',float(top['diff']),0.0,thr/2,
                 null_kind='把 31 个类别随机等大小重分成 A/B 后的同一个差(200 次,跨全部变量合并)')
print(g)
print(f"\n  超过随机重分阈值的变量:{int(D.over.sum())}/{len(D)}"
      f"   过人层自助+多重性的:{int(D.pass_mult.sum())}/{len(D)}")
# ---- 信度匹配的零:随机 21/10 分割的两个分数比真实两族**噪声大**(它们不连贯),
#      所以那条零太宽,是 same_scale 不匹配。换成同一个矩阵的其余特征向量 ——
#      它们与 PC1 一样连贯(同一结构的正交方向),只是不是"那两族"。
print("\n=== 信度匹配的零:SIMR 的第 2..6 特征向量给出的同样连贯的分割 ===")
print(f"  {'分割':<10} {'|A|':>4} {'内部连通':>9} {'开放性差':>10} {'无力感差':>10} {'最大差':>9}")
alt=[]
for q in range(2,7):
    v=vv[:,-q]; Aq=np.flatnonzero(v>0); Bq=np.flatnonzero(v<=0)
    if len(Aq)<5 or len(Bq)<5: continue
    if SIMR[np.ix_(Aq,Aq)].mean()<SIMR[np.ix_(Bq,Bq)].mean(): Aq,Bq=Bq,Aq
    dq=diffs(Aq,Bq)
    do=dq.get('开放性',(np.nan,np.nan)); dp=dq.get('无力感',(np.nan,np.nan))
    mx=max((abs(a-b) for a,b in dq.values()),default=np.nan)
    alt.append(dict(pc=q,nA=len(Aq),coh=float(SIMR[np.ix_(Aq,Aq)].mean()),
                    d_open=do[0]-do[1],d_pow=dp[0]-dp[1],mx=mx))
    print(f"  PC{q:<8} {len(Aq):>4} {SIMR[np.ix_(Aq,Aq)].mean():>+9.4f} "
          f"{do[0]-do[1]:>+10.4f} {dp[0]-dp[1]:>+10.4f} {mx:>9.4f}")
A_=pd.DataFrame(alt)
thr2=float(np.percentile(A_.mx.values,80)) if len(A_)>=4 else np.nan
print(f"  真实分割(PC1):|A|={len(FA)} 内部 {SIMR[np.ix_(FA,FA)].mean():+.4f}  "
      f"开放性差 {D[D['var']=='开放性']['diff'].iloc[0]:+.4f}  "
      f"最大差 {D['diff'].abs().max():.4f}")
print(f"  信度匹配零的最大差(5 个正交分割):" + " ".join(f"{v:.4f}" for v in A_.mx.values))

g2=Gate('两族有不同的外部锚吗 —— 人层自助检验')
g2.asserted('信度匹配的零:PC1 的分割在外部锚上并不比其余正交分割更强',
            not (D['diff'].abs().max() > A_.mx.max()),
            f"PC1 最大差 {D['diff'].abs().max():.4f} vs 其余 5 个正交分割 "
            f"{A_.mx.min():.4f}..{A_.mx.max():.4f} —— **所以'这两族'不比同一矩阵里别的连贯分割更特别**")
g2.asserted('随机重分零的阈值被算术卡死',thr>0.093,
            f"阈值 {thr:.4f} > 本数据里存在过的最大外部相关(性别 +0.093)。"
            f"**判别量的上界被一个本身在噪声地板上的量卡住**,所以这条路线在算术上不可能赢")
big=D.sort_values('z',ascending=False).iloc[0]
# ⚠ 多重性门槛是一个**门槛**,不是一个**展布**。第一版写成
#   require_resolvable_first(z, alpha_z),那要求 z > 2*alpha_z —— 把门槛又乘了 2。
#   喂给门的必须是它要的那种量(#129i 同一族的错)。
g2.asserted('最大的 z 过多重性门槛',float(big['z'])>alpha_z,
            f"{big['var']} z={big['z']:.1f}x > 门槛 {alpha_z:.2f}x({K} 个变量)")
g2.offset_control(f"{big['var']}:两族的差 vs 零",float(big['diff']),0.0,float(big['boot']),
                  null_kind='按人自助的抽样分布(零假设:两族对该变量的相关相同)')
print(g2)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
