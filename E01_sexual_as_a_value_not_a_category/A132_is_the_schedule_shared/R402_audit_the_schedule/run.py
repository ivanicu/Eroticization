import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A132 R402 -- 那张「人群共享的时间表」,证据够不够它被引用的次数

`#120`:性兴趣按一张**人群共享**的时间表到来 —— 外观 14.0 · 身体部位 14.4 · 衣物 14.7 →
权力动态 16.8 · 束缚 16.9 · 精神改变 17.0(内容类早、关系类晚,差 2–3 年)。
**这条被 `#130` `#333c` `#336` 等十几轮引用,而它自己的证据只有一组均值。**

ESTIMAND        ① 31 个类别起始年龄均值的**自助 95% 区间**;
                ② **不依赖「内容/关系」这个二分**的量:起始年龄与该类别**稀有度**的相关
                   (⚠ `#384` 已报「块层 `c3⁻` 载荷对应」不可算 —— **不再试**);
                ③ **顺序**在跨性别与跨当前年龄档下稳不稳(「共享」要求顺序稳)。
KILL            **若区间彼此大量重叠 -> 那张表是一个排序噪声,「2–3 年」的说法要收窄;
                若顺序跨组稳 -> 「共享」这个词站得住,而那是这条声明真正的内容。**
POSITIVE CTRL   合成一个**已知顺序**的时间表 -> 跨组顺序相关必须 ≈ 1。
NEGATIVE CTRL   打乱人 -> 跨组顺序相关必须 ≈ 0。
⚠ 多重性       31 个类别 -> 报**分布与顺序统计量**,不逐格宣称。
⚠ 二分归属     「内容类 / 关系类」**是谁分的**必须先说清 —— 若无来源,**本轮不使用它**。
IMPOSSIBLE      起始年龄是**回溯自报**;顺序稳只说明回忆一致,不说明事件一致。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')

BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
NC=ONS.shape[1]; HAS=np.isfinite(ONS)
NAM=[re.sub(r'^How old were you when you first experienced (sexual )?interest in ','',str(c))[:44]
     if (re:=__import__('re')) else str(c)[:44] for c in ons]
PREV=HAS.mean(0); RAR=-np.log(np.clip(PREV,1e-4,1.))
mu=np.array([np.nanmean(ONS[HAS[:,j],j]) for j in range(NC)])
n_j=HAS.sum(0)
rg=np.random.default_rng(1202); B=600
CI=np.zeros((NC,2))
for j in range(NC):
    v=ONS[HAS[:,j],j]
    bs=np.array([v[rg.integers(0,len(v),len(v))].mean() for _ in range(B)])
    CI[j]=np.percentile(bs,[2.5,97.5])
o=np.argsort(mu)
print(f"① 31 个类别的起始年龄均值 + 自助 95% 区间(按均值排序,前 5 / 后 5):")
for j in o[:5]: print(f"   {mu[j]:5.2f} [{CI[j,0]:5.2f}, {CI[j,1]:5.2f}]  n={n_j[j]:>5,}  {NAM[j]}")
print("   ...")
for j in o[::-1][:5]: print(f"   {mu[j]:5.2f} [{CI[j,0]:5.2f}, {CI[j,1]:5.2f}]  n={n_j[j]:>5,}  {NAM[j]}")
span=mu.max()-mu.min(); wid=float(np.mean(CI[:,1]-CI[:,0]))
ov=sum(1 for a in range(NC) for b in range(a+1,NC) if CI[a,0]<CI[b,1] and CI[b,0]<CI[a,1])
print(f"\n   全距 **{span:.2f} 年** · 区间平均宽度 **{wid:.2f} 年** · "
      f"**{ov}/{NC*(NC-1)//2}** 对区间重叠({100*ov/(NC*(NC-1)//2):.0f}%)")
print(f"\n② 不依赖任何二分的量:`corr(类别起始均值, 类别稀有度)` = "
      f"**{np.corrcoef(mu,RAR)[0,1]:+.4f}**(31 点)")
sex=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ag=d['age'].map(AGE).values.astype(float)
def means(mask):
    return np.array([np.nanmean(ONS[mask&HAS[:,j],j]) if (mask&HAS[:,j]).sum()>=60 else np.nan
                     for j in range(NC)])
from scipy.stats import spearmanr
def rho(a,b):
    k=np.isfinite(a)&np.isfinite(b); return float(spearmanr(a[k],b[k]).statistic) if k.sum()>=10 else np.nan
gm,gf=means(sex==1),means(sex==0)
print(f"\n③ 顺序稳不稳(Spearman ρ):")
print(f"   男 vs 女 **{rho(gm,gf):+.4f}**")
BANDS=[(0,18),(18,24),(24,99)]
mb=[means((ag>=a)&(ag<b)) for a,b in BANDS]
print(f"   年龄档两两:" + ' · '.join(f"{rho(mb[i],mb[j]):+.4f}" for i in range(3) for j in range(i+1,3)))
rgp=np.random.default_rng(5)
p1=means(rgp.random(NN)<0.5); p2=means(~(rgp.random(NN)<0.5))
print(f"   **随机劈半(参照臂)** **{rho(p1,p2):+.4f}**")
TRUE=np.linspace(12,20,NC)
SYN=np.where(HAS,TRUE[None,:]+rgp.standard_normal((NN,NC))*3,np.nan)
sm=lambda mask: np.array([np.nanmean(SYN[mask&HAS[:,j],j]) if (mask&HAS[:,j]).sum()>=60 else np.nan
                          for j in range(NC)])
print(f"\n正对照(合成已知顺序的时间表):男 vs 女 **{rho(sm(sex==1),sm(sex==0)):+.4f}**")
PERM=ONS[rgp.permutation(NN)]
pmw=np.array([np.nanmean(PERM[(sex==1)&np.isfinite(PERM[:,j]),j]) for j in range(NC)])
pmf=np.array([np.nanmean(PERM[(sex==0)&np.isfinite(PERM[:,j]),j]) for j in range(NC)])
print(f"⚠ 负对照说明:打乱**人**不破坏类别顺序(每列还是那一列)——"
      f"正确的负对照是**打乱每个人的类别标签**:")
SHF=np.array([ONS[i][rgp.permutation(NC)] for i in range(NN)])
sf=lambda mask: np.array([np.nanmean(SHF[mask&np.isfinite(SHF[:,j]),j]) for j in range(NC)])
print(f"   打乱类别标签后:男 vs 女 **{rho(sf(sex==1),sf(sex==0)):+.4f}**")
# ⚠ #300a:上页面前发明一个旋钮 —— **换掉我编的分档中点**(`#332c` 同款,这次作用在顺序上)。
RANKMAP={k:i+1.0 for i,k in enumerate(BIN)}
ALTMAP={'0-4yo':2.5,'5-6yo':6,'7-8yo':8,'9-10yo':10,'11-12yo':12,'13-14yo':14,
        '15-16yo':16,'17-18yo':18,'19-25yo':21,'26yo+':32}
print(f"\n发明的旋钮 · 起始年龄的数值编码:")
KN=[]
for tag,MP in (('我的中点',BIN),('序号 1..10',RANKMAP),('另一套中点(26+→32)',ALTMAP)):
    O2=np.column_stack([d[c].map(MP).values.astype(float) for c in ons])
    H2=np.isfinite(O2)
    m2=lambda mask: np.array([np.nanmean(O2[mask&H2[:,j],j]) if (mask&H2[:,j]).sum()>=60 else np.nan
                              for j in range(NC)])
    mu2=np.array([np.nanmean(O2[H2[:,j],j]) for j in range(NC)])
    KN.append((tag,rho(m2(sex==1),m2(sex==0)),float(np.corrcoef(mu2,RAR)[0,1])))
    print(f"   {tag:<20} 男女顺序 ρ **{KN[-1][1]:+.4f}** · `corr(均值, 稀有度)` **{KN[-1][2]:+.4f}**")
kok=(max(x[1] for x in KN)-min(x[1] for x in KN))<0.10 and max(abs(x[2]) for x in KN)<0.25

T=pd.DataFrame([dict(v_cat=NAM[j],v_mu=float(mu[j]),v_lo=float(CI[j,0]),v_hi=float(CI[j,1]),
                     v_n=int(n_j[j]),v_rar=float(RAR[j])) for j in range(NC)])
check_columns(T,'R402'); T.to_csv(pathlib.Path(__file__).parent/'results'/'schedule.csv',index=False)
gg=Gate('那张共享时间表的证据')
gg.asserted('★ 发明的旋钮:换三套年龄编码,顺序稳定性与「不是稀有度」两条都不动',kok,
            ' · '.join(f"{t} ρ={r_:+.4f} rar={c_:+.4f}" for t,r_,c_ in KN))
gg.asserted('★ 正对照:合成已知顺序 -> 跨性别顺序相关必须 ≈ 1',
            rho(sm(sex==1),sm(sex==0))>0.9,f"{rho(sm(sex==1),sm(sex==0)):+.4f}")
gg.asserted('★ 负对照:打乱**类别标签**后跨性别顺序相关必须 ≈ 0',
            abs(rho(sf(sex==1),sf(sex==0)))<0.3,f"{rho(sf(sex==1),sf(sex==0)):+.4f}")
gg.asserted('★ 注册的 kill ①:区间是否彼此大量重叠',ov/(NC*(NC-1)//2)>0.5,
            f"{ov}/{NC*(NC-1)//2} 对重叠({100*ov/(NC*(NC-1)//2):.0f}%)· "
            f"全距 {span:.2f} 年 vs 区间宽 {wid:.2f} 年")
gg.asserted('★ 注册的 kill ②:顺序跨组稳不稳(男女 ρ > 0.8)',rho(gm,gf)>0.8,
            f"男vs女 {rho(gm,gf):+.4f} · 随机劈半参照臂 {rho(p1,p2):+.4f}")
gg.asserted('⚠ 二分归属:「内容类/关系类」是谁分的',True,
            '**本轮不使用它** —— ② 用的是与稀有度的相关,不依赖任何二分')
gg.asserted('⚠ 边界:起始年龄是回溯自报',True,'顺序稳只说明**回忆**一致,不说明**事件**一致')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
