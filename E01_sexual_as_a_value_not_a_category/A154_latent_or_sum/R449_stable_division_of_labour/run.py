import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A154 R449 -- 两半的分工比例,在四个结局上稳不稳

`#404d`:两半在**区分那群人**上对等(`#403a`),在**治疗性**上却几乎全是「常规也管用」那一半在做工。
**⇒ 若「两边都有」真是一件事,两半的分工比例应当在多个结局上**稳定**;
若它只是两件事的加法,**比例会随结局漂移**。**

四个结局(**都已在页面上,没有一个是为这个构造造的**):
**羞耻** · **能不能改**(`BELIEF`)· **治疗性** · **实践了多少**(`ACTED`)。

ESTIMAND        每个结局上拟合 `结局 ~ z(S) + (−z(五题)) + c3⁻ + 类别数`
                (**控制集对四个结局相同,且从不控制结局自己**);
                分工 = **角度** `atan2(b_五题, b_S)`(以及**比值** `b_S / b_五题`);
                主量 = **四个角度的离散度**(圆离散度)。
判据(**先标支**,`#379c`)
                【两支】**两种度量(角度 / 比值)结论一致** —— 否则是度量在说话(`#444`);
                        负对照用**越阈率**;guard 26 **显式传 branch**。
                【非零支】离散度**低于** offset 零的下侧 -> 分工稳定 -> **一件事**;
                【零支】落在零里 -> 与随机无异 -> **只是加法**,启用 MDE。
⚠ 零的种类     `offset_control`:**四个角度的离散度的零绝不是零** ——
                四个系数各有采样噪声,即使真值完全相同也会有离散度。
                零 = **`lib.nulls.perm_in` 打乱结局**后重算四个角度的离散度分布。
IMPOSSIBLE      ① 四个结局彼此相关 -> 四个角度**不是四次独立观测**,离散度会被人为压低,
                   **而这偏向「稳定」** -> 结论方向上**不保守**,必须明说;
                ② 角度对**两个系数都接近零**的结局极不稳;
                ③ 「分工稳定」不等于「存在一个潜变量」(`#404c` 那条边界仍在)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
THC=next(c for c in d.columns if 'vmq8jqw' in str(c)); TH=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
BC=next(c for c in d.columns if '7lgg41e' in str(c))
BMAP={'Impossible':0.,'With an extreme amount of effort, maybe':1.,'With a lot of effort, yes':2.,
      'With some effort, yes':3.,'With little effort, yes':4.}
BEL=d[BC].map(BMAP).values.astype(float)
AC=next(c for c in d.columns if '41kpfir' in str(c)); ACT=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
OUT={'羞耻':sh,'能不能改':BEL,'治疗性':TH,'实践了多少':ACT}
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(ncat)&np.isfinite(INT)
for v in OUT.values(): M&=np.isfinite(v)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
print(f"⚠ **`#392e`:四个结局各自先看清楚**")
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
for nm,v in OUT.items():
    g=M&np.isfinite(anc)
    print(f"   {nm:<8} 取值 {np.unique(v[M]).tolist()[:7]} · 众数 **{float(pd.Series(v[M]).mode().iloc[0]):g}** · "
          f"与锚相关 **{np.corrcoef(v[g],anc[g])[0,1]:+.4f}**")
A=np.full(NN,np.nan); Bv=np.full(NN,np.nan); A[M]=z(S,M); Bv[M]=-z(INT,M)
CTRL=[C3,ncat]
def coefs(y):
    X=np.column_stack([np.ones(n),z(A,M),z(Bv,M)]+[z(v,M) for v in CTRL]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); return float(b[1]),float(b[2])
def disp(angles):
    """圆离散度:1 − |平均单位向量|。0 = 完全一致,1 = 完全散开。"""
    a=np.asarray(angles); return float(1-abs(np.mean(np.exp(1j*a))))
print(f"\nn=**{n:,}** · 控制集对四个结局相同(`c3⁻` · 类别数),**从不控制结局自己**\n")
rows=[]; ang=[]; rat=[]
for nm,y in OUT.items():
    bS,bI=coefs(y); th_=float(np.arctan2(bI,bS)); r_=bS/bI if abs(bI)>1e-9 else np.nan
    rows.append(dict(v_out=nm,v_bS=bS,v_bI=bI,v_angle=th_,v_ratio=r_))
    ang.append(th_); rat.append(r_)
    print(f"   {nm:<8} z(S) **{bS:+.4f}** · −z(五题) **{bI:+.4f}** · "
          f"角度 **{np.degrees(th_):+7.2f}°** · 比值 **{r_:+.3f}**")
T=pd.DataFrame(rows); check_columns(T,'R449')
T.to_csv(pathlib.Path(__file__).parent/'results'/'labour.csv',index=False)
D_ANG=disp(ang); D_RAT=float(np.std(rat))
print(f"\n主量:四个角度的**圆离散度** = **{D_ANG:.5f}**(0=完全一致)· 比值的 sd = **{D_RAT:.4f}**")
NP_=400; nulA=[]; nulR=[]
for s_ in range(NP_):
    aa=[];rr=[]
    for k,(nm,y) in enumerate(OUT.items()):
        bS,bI=coefs(perm_in(y,M,8600+s_*7+k)); aa.append(np.arctan2(bI,bS))
        rr.append(bS/bI if abs(bI)>1e-9 else np.nan)
    nulA.append(disp(aa)); nulR.append(float(np.nanstd(rr)))
nulA=np.array(nulA); nulR=np.array(nulR)
LOA=float(np.percentile(nulA,5)); LOR=float(np.percentile(nulR,5))
print(f"\n⚠ offset 零(**`lib.nulls.perm_in` 打乱结局**;**四个系数各有采样噪声,"
      f"即使真值相同也有离散度 -> 这个零绝不是零**):")
print(f"   角度离散度零 **{nulA.mean():.5f} ± {nulA.std():.5f}** · 5 分位(下侧)**{LOA:.5f}**")
print(f"   比值 sd 零   **{nulR.mean():.4f} ± {nulR.std():.4f}** · 5 分位(下侧)**{LOR:.4f}**")
OKA=D_ANG<LOA; OKR=D_RAT<LOR
print(f"   -> 角度 {'**低于零的下侧 -> 稳定**' if OKA else '**落在零里 -> 与随机无异**'} · "
      f"比值 {'**低于**' if OKR else '**落在零里**'}")
negs=[]
for s_ in range(200):
    aa=[]
    for k,(nm,y) in enumerate(OUT.items()):
        bS,bI=coefs(perm_in(y,M,98000+s_*7+k)); aa.append(np.arctan2(bI,bS))
    negs.append(disp(aa))
negs=np.array(negs); rate=float((negs<LOA).mean())
print(f"\n负对照(**越界率**,打乱结局 200 次落到下侧 5% 之外的比例):**{100*rate:.1f}%**(合格 1–12%)")
CONS=(OKA==OKR)
MDE=None
print(f"\nguard 26 = **MDE 扫描**(种一个真实的共同分工),每级 20 次:")
for gg in (0.05,0.10,0.20,0.35):
    hit=0
    for s_ in range(20):
        rg=np.random.default_rng(150+int(gg*100)*97+s_); aa=[]
        for k,(nm,y) in enumerate(OUT.items()):
            yy=np.full(NN,np.nan)
            yy[M]=gg*(0.5*z(A,M)+0.5*z(Bv,M))+rg.standard_normal(n)
            bS,bI=coefs(yy); aa.append(np.arctan2(bI,bS))
        if disp(aa)<LOA: hit+=1
    print(f"   共同分工强度 **{gg:.2f}** -> 检出 **{hit}/20 = {hit*5:>3d}%**")
    if MDE is None and hit>=16: MDE=gg
MDE_=MDE if MDE else 0.50
g=Gate('两半的分工比例在四个结局上稳不稳')
g.asserted('★【两支】两种度量(角度 / 比值)结论一致 —— 否则是度量在说话',CONS,
           f"角度 {OKA} · 比值 {OKR}",kind='control')
g.asserted('★【两支】负对照:**越界率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 有意义的共同分工强度 0.20',MDE_,0.20,True,
    what='MDE 扫描 80% 检出',branch='non_null' if OKA else 'null')
g.asserted('★【两支】offset 零非退化(四个系数各有采样噪声)',nulA.std()>0,
           f"{nulA.mean():.5f} ± {nulA.std():.5f}",kind='control')
if CONS and 0.01<=rate<=0.12:
    g.asserted('★【非零支】离散度**低于** offset 零的下侧 -> 分工稳定 -> 一件事',OKA,
               f"角度 {D_ANG:.5f} vs 下侧 {LOA:.5f} · 比值 {D_RAT:.4f} vs {LOR:.4f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **IMPOSSIBLE ①**:四个结局彼此相关 -> 四个角度**不是四次独立观测**,"
      f"离散度被人为压低,**而那偏向「稳定」** —— 结论方向上**不保守**。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
