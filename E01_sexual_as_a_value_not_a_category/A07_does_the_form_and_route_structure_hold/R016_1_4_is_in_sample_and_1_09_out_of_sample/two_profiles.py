import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A103 R354 -- 这两份羞耻,是不是同一种羞耻

`#308a`:两条路的贡献**相加**。**但「相加」不区分「同一种感受的两个来源」和「两种不同的感受」。**

⚠ 本轮做**两个**检验,而第二个比第一个有力:
- **① 结局剖面**(注册的那个):两个正交残差各自与其余 28 个结局的相关剖面,判两条剖面的相关。
  **guard 15 必须挂上** —— 而按构造两者分数层相关是 0,所以 guard 15 会 FAIL,
  **那正是要写清楚的地方**:两个正交量的剖面相似**只能**读成「它们落在同一片结局上」。
- **② 调节剖面**(更有力):对 10 个调节变量各按中位切两层,在每层里分别估两条路 ↔ 羞耻,
  得两条 20 维**条件效应**剖面,判它们的相关。
  **同一种羞耻应当对同一批东西同样地敏感** —— 这比「落在同一片结局上」强,
  因为它问的是**同一个结局在不同条件下怎么变**,而不是**哪些结局被碰到**。

KILL            **② 的剖面相关明显高于它自己的零 -> 两条路对条件的反应一致,
                「同一种羞耻」这个读法可以活着(但仍不是证明);
                不高 -> 两条路是两种不同的感受,而公开页面上「两条路」要改写成「两种羞耻」。**
POSITIVE CTRL   合成**同一潜变量驱动**的两个量 -> ② 必须高;
                合成**两个对调节反应相反**的量 -> ② 必须低。两个都要,否则读不出量程。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ 预先说明     本轮诚实的产出多半是一个**上界**,不是一个命名 —— 写在跑之前。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
S=Q[0]; C3=-Q[4]
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&ok
def resid(a,b,m):
    out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
    out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
RS=resid(S,C3,m0); RC=resid(C3,S,m0)
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(m0 if m is None else m)
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>150 else np.nan
print(f"两个正交残差:corr(RS, RC) = **{cor(RS,RC):+.6f}**(按构造 ≈ 0)")
print(f"各自 ↔ 羞耻:RS **{cor(RS,sh):+.4f}** · RC **{cor(RC,sh):+.4f}**")

# ---- ① 结局剖面(注册的) ----
OTH=[(nm,y.astype(float)) for nm,y in OUT if str(nm)!=SHAME]
pa=np.array([cor(RS,y) for _,y in OTH]); pb=np.array([cor(RC,y) for _,y in OTH])
k=np.isfinite(pa)&np.isfinite(pb)
prof=float(np.corrcoef(pa[k],pb[k])[0,1])
print(f"\n① **结局剖面**({int(k.sum())} 个结局):两条剖面的相关 **{prof:+.4f}**")

# ---- ② 调节剖面(更有力) ----
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
MOD={'年龄':d['age'].map(AGE).values.astype(float),
 '开放性':pd.to_numeric(d['opennessvariable'],errors='coerce').values.astype(float),
 '尽责性':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce').values.astype(float),
 '外向性':pd.to_numeric(d['extroversionvariable'],errors='coerce').values.astype(float),
 '神经质':pd.to_numeric(d['neuroticismvariable'],errors='coerce').values.astype(float),
 '宜人性':pd.to_numeric(d['agreeablenessvariable'],errors='coerce').values.astype(float),
 '无力感':pd.to_numeric(d['powerlessnessvariable'],errors='coerce').values.astype(float),
 '性别':pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}).values.astype(float),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1}).values.astype(float)}
def modprof(y):
    a=[];b=[];lab=[]
    for nm,v in MOD.items():
        f=np.isfinite(v)&m0
        if f.sum()<600: continue
        th=np.median(v[f])
        for side,sel in (('低',v<=th),('高',v>th)):
            mm=f&sel
            if mm.sum()<250: continue
            a.append(cor(RS,y,mm)); b.append(cor(RC,y,mm)); lab.append(f"{nm}·{side}")
    return np.array(a),np.array(b),lab
ma,mb,lab=modprof(sh)
kk=np.isfinite(ma)&np.isfinite(mb)
mprof=float(np.corrcoef(ma[kk],mb[kk])[0,1])
print(f"② **调节剖面**({int(kk.sum())} 个条件格):两条条件效应剖面的相关 **{mprof:+.4f}**")
for i in np.argsort(-np.abs(ma-mb))[:4]:
    print(f"     分歧最大:{lab[i]:<12} RS {ma[i]:+.4f} · RC {mb[i]:+.4f}(差 {ma[i]-mb[i]:+.4f})")
def perm_finite(v,seed):
    z=v.copy(); j=np.flatnonzero(np.isfinite(z))
    z[j]=z[np.random.default_rng(seed).permutation(j)]; return z
nul=[]
for i in range(12):
    a2,b2,_=modprof(perm_finite(sh,300+i)); k2=np.isfinite(a2)&np.isfinite(b2)
    if k2.sum()>=10: nul.append(float(np.corrcoef(a2[k2],b2[k2])[0,1]))
print(f"   打乱人的零:**{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
rg=np.random.default_rng(21); zS=np.where(m0,(RS-np.nanmean(RS[m0]))/np.nanstd(RS[m0]),np.nan)
zC=np.where(m0,(RC-np.nanmean(RC[m0]))/np.nanstd(RC[m0]),np.nan)
AGEv=MOD['年龄']; hi=np.isfinite(AGEv)&(AGEv>np.nanmedian(AGEv))
print(f"\n正对照:")
CT={}
for tag,w in (('同一潜变量驱动',None),('对调节反应相反',hi)):
    y=np.full(NN,np.nan)
    base=0.25*(zS[m0]+zC[m0])+rg.standard_normal(int(m0.sum()))
    if w is None: y[m0]=base
    else:
        g=np.where(w[m0],1.0,-1.0)
        y[m0]=0.25*(zS[m0]*g - zC[m0]*g)+rg.standard_normal(int(m0.sum()))
    a3,b3,_=modprof(y); k3=np.isfinite(a3)&np.isfinite(b3)
    CT[tag]=float(np.corrcoef(a3[k3],b3[k3])[0,1])
    print(f"   {tag:<16} ② 调节剖面相关 **{CT[tag]:+.4f}**")
# ⚠ 剖面为什么没功效:**两条路的条件效应在 19 个格里几乎是常数** —— 没有变化就没有剖面可相关。
#    把这一点直接测出来,比继续修剖面有用。
def cell_ns():
    ns=[]
    for nm,v in MOD.items():
        f=np.isfinite(v)&m0
        if f.sum()<600: continue
        th=np.median(v[f])
        for sel in (v<=th,v>th):
            mm=f&sel
            if mm.sum()>=250: ns.append(int(mm.sum()))
    return np.array(ns)
NS=cell_ns(); se=1/np.sqrt(NS)
print(f"\n⚠ 剖面为什么没功效 —— 两条路的条件效应几乎是常数:")
for tag,v in (('RS(位置侧)',ma[kk]),('RC(广度型侧)',mb[kk])):
    print(f"   {tag:<14} 均值 **{v.mean():+.4f}** · 跨 {len(v)} 格的 sd **{v.std():.4f}** · "
          f"各格自身 se 中位 **{np.median(se):.4f}** -> **{v.std()/np.median(se):.2f}×**")
flat=max(ma[kk].std(),mb[kk].std())/np.median(se)
print(f"   -> 两条路都**没有被这 10 个调节变量调节**(最大 {flat:.2f}× 自身 se);"
      f"**没有变化,就没有剖面可以比**。")
LB=[l for l,g in zip(lab,kk) if g]
print(f"\n   全部 {len(LB)} 个条件格(RS / RC):")
for l,x,y2 in zip(LB,ma[kk],mb[kk]): print(f"     {l:<16} {x:+.4f} / {y2:+.4f}")

# ⚠ #300a:给「没有被调节」这条新结论发明一个旋钮 —— 换成**极端三分位**切,而不是中位切。
def modprof_q(y,lo,hi):
    a=[];b=[];ns=[]
    for nm,v in MOD.items():
        f=np.isfinite(v)&m0
        if f.sum()<600: continue
        ql,qh=np.quantile(v[f],lo),np.quantile(v[f],hi)
        for sel in (v<=ql,v>=qh):
            mm=f&sel
            if mm.sum()>=250: a.append(cor(RS,y,mm)); b.append(cor(RC,y,mm)); ns.append(int(mm.sum()))
    a,b,ns=np.array(a),np.array(b),np.array(ns); k=np.isfinite(a)&np.isfinite(b)
    return a[k],b[k],1/np.sqrt(ns[k])
KN=[]
for tag,(lo,hi) in (('中位(0.5/0.5)',(0.5,0.5)),('三分位(1/3,2/3)',(1/3,2/3)),('极端四分位(.25/.75)',(0.25,0.75))):
    a4,b4,se4=modprof_q(sh,lo,hi)
    KN.append((tag,a4.std()/np.median(se4),b4.std()/np.median(se4),len(a4)))
print(f"\n发明的旋钮(切分规则)· 跨格 sd 相对各格自身 se:")
for tag,x,y2,n_ in KN: print(f"   {tag:<20} RS **{x:.2f}×** · RC **{y2:.2f}×**  ({n_} 格)")
KMAX=max(max(x,y2) for _,x,y2,_ in KN)

T=pd.DataFrame([dict(v_test='①结局剖面',v_val=prof),dict(v_test='②调节剖面',v_val=mprof),
                dict(v_test='②的零',v_val=float(np.mean(nul)))])
check_columns(T,'R354'); T.to_csv(pathlib.Path(__file__).parent/'results'/'two_profiles.csv',index=False)
gg=Gate('这两份羞耻是不是同一种羞耻')
gg.asserted('★ 正对照:同一潜变量驱动 -> ② 必须高;对调节反应相反 -> ② 必须低',
            CT['同一潜变量驱动']>0.5 and CT['对调节反应相反']<0.2,
            f"同一潜变量 **{CT['同一潜变量驱动']:+.4f}** · 反应相反 **{CT['对调节反应相反']:+.4f}** —— 量程")
gg.negative_control('★ 负对照:打乱人后的 ② 剖面相关',float(np.mean(nul)),mprof,
    null_spread=float(np.std(nul)),null_kind='`perm_finite` 题内跨人打乱 —— 保住缺失格局(#264b)')
gg.profile_similarity_is_not_identity('⚠ guard 15:① 结局剖面 vs 分数层相关(按构造为 0)',
                                      prof,cor(RS,RC))
gg.asserted('★ 发明的旋钮:三种切分下「没有被调节」还成不成立(全部 < 2× 自身 se)',
            KMAX<2.0,' · '.join(f"{t} RS {x:.2f}× RC {y2:.2f}×" for t,x,y2,_ in KN))
gg.asserted('★ 为什么 ② 没功效:两条路的条件效应跨 19 格的 sd,相对各格自身 se',
            flat<2.0,
            f"RS sd {ma[kk].std():.4f} · RC sd {mb[kk].std():.4f} · 各格 se 中位 {np.median(se):.4f}"
            f" -> 最大 **{flat:.2f}×** —— **两条路都没有被这 10 个变量调节,所以没有剖面可比**")
gg.asserted('★ 注册的 kill:② 调节剖面相关是否明显高于它自己的零',
            (mprof-np.mean(nul))>2*np.std(nul),
            f"② **{mprof:+.4f}** vs 零 {np.mean(nul):+.4f} ± {np.std(nul):.4f} "
            f"({abs(mprof-np.mean(nul))/max(2*np.std(nul),1e-9):.1f}× 的 2×展布)")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
