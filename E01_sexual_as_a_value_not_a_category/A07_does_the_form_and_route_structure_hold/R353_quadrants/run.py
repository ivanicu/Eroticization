import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A102 R353 -- 羞耻的两条路,在同一个人身上是相加还是相乘

`#231` 起羞耻有两条路(`S` +0.1185 · `c3` +0.1286,彼此相关 0.2036,联合 +0.0963/+0.1090)。
**它们在人群层近独立 —— 但那不等于它们在同一个人身上是两件事。**
`#286`–`#307` 这一整段都在查仪器。**回到问题本身。**

⚠ **符号约定,先写下来**(`#306b` 的教训):`c3` 与羞耻的相关是**负**的(−0.1278),
所以本轮全程用 **`c3⁻ = −c3`**,让「两条路都指向更多羞耻」——
这只是一个**约定**,不是一个发现。

ESTIMAND        按 (`S`, `c3⁻`) 中位切成四象限,报四格的羞耻均值;
                **交互 = 观测(高,高) − 加性预测**,加性预测 = 低低 + (高低−低低) + (低高−低低)。
                另报**连续版**:羞耻 ~ zS + zc3⁻ + zS·zc3⁻ 的交互系数。
KILL            **加性 -> 四格落在加性预测上,连续交互项落在它自己的置换零里;
                相乘 -> 高×高 那一格明显高于加性预测。**
POSITIVE CTRL   合成一个**只有加性**的结局 -> 四格必须落回加性预测(否则读到的是切分伪影);
                再合成一个**带已知交互**的结局 -> 必须被抓到(灵敏度)。
NEGATIVE CTRL   `perm_finite` 题内跨人打乱。
⚠ 报           **每格的 n 与自身展布** —— 四格切分让每格 n 掉到约四分之一(`#348` 的 MDE 纪律)。
⚠ KNOB         切分规则本身是旋钮(中位 / 极端三分位 / 连续交互项)-> **报规格曲线,不报一格**。
IMPOSSIBLE      交互是**统计的**,不是机制的:即使有交互,它也不能说「一条路放大了另一条」,
                只能说「两者同高时比相加更高」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok)
Q=fit_apply(ALLR,ALLR)
S=Q[0]; C3=-Q[4]                                    # ⚠ 约定:c3⁻ = −c3
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&ok
print(f"n={int(m0.sum()):,} · corr(S, c3⁻) = **{np.corrcoef(S[m0],C3[m0])[0,1]:+.4f}**")
print(f"corr(S, 羞耻) = **{np.corrcoef(S[m0],sh[m0])[0,1]:+.4f}** · "
      f"corr(c3⁻, 羞耻) = **{np.corrcoef(C3[m0],sh[m0])[0,1]:+.4f}**")
def quad(y,lo_q,hi_q,mask=None):
    m=m0 if mask is None else (m0&mask)
    a,b=S[m],C3[m],; yy=y[m]
    qa=(np.quantile(a,lo_q),np.quantile(a,hi_q)); qb=(np.quantile(b,lo_q),np.quantile(b,hi_q))
    cells={}
    for i,(la,ha) in enumerate([(-np.inf,qa[0]),(qa[1],np.inf)]):
        for j,(lb,hb) in enumerate([(-np.inf,qb[0]),(qb[1],np.inf)]):
            k=(a>=la)&(a<=ha)&(b>=lb)&(b<=hb)
            cells[(i,j)]=(float(yy[k].mean()),int(k.sum()),float(yy[k].std()/max(np.sqrt(k.sum()),1)))
    add=cells[(0,0)][0]+(cells[(1,0)][0]-cells[(0,0)][0])+(cells[(0,1)][0]-cells[(0,0)][0])
    inter=cells[(1,1)][0]-add
    se=np.sqrt(sum(cells[c][2]**2 for c in cells))
    return cells,add,inter,se
def cont(y,mask=None):
    m=m0 if mask is None else (m0&mask)
    z=lambda v:(v[m]-v[m].mean())/v[m].std()
    X=np.column_stack([np.ones(m.sum()),z(S),z(C3),z(S)*z(C3)])
    yy=(y[m]-y[m].mean())/y[m].std()
    return float(np.linalg.lstsq(X,yy,rcond=None)[0][3])
print(f"\n{'切分':<14}{'低低':>9}{'高低':>9}{'低高':>9}{'高高':>9}{'加性预测':>10}{'交互':>10}{'±se':>8}  最小 n")
SPEC=[]
for tag,(lo,hi) in (('中位',(0.5,0.5)),('三分位',(1/3,2/3)),('四分位',(0.25,0.75)),('极端五分位',(0.2,0.8))):
    cells,add,inter,se=quad(sh,lo,hi)
    nmin=min(c[1] for c in cells.values())
    print(f"{tag:<14}{cells[(0,0)][0]:>9.4f}{cells[(1,0)][0]:>9.4f}{cells[(0,1)][0]:>9.4f}"
          f"{cells[(1,1)][0]:>9.4f}{add:>10.4f}{inter:>+10.4f}{se:>8.4f}  {nmin:,}")
    SPEC.append(dict(v_cut=tag,inter=inter,se=se,nmin=nmin))
ci=cont(sh)
rgp=np.random.default_rng(909)
def perm_finite(v,seed):
    z=v.copy(); j=np.flatnonzero(np.isfinite(z))
    z[j]=z[np.random.default_rng(seed).permutation(j)]; return z
nul=[cont(perm_finite(sh,800+i)) for i in range(20)]
print(f"\n连续版交互系数 **{ci:+.4f}** · 置换零 **{np.mean(nul):+.4f} ± {np.std(nul):.4f}** "
      f"-> **{abs(ci-np.mean(nul))/max(2*np.std(nul),1e-9):.1f}× 的 2×展布**")
T=pd.DataFrame(SPEC); check_columns(T,'R353')
T.to_csv(pathlib.Path(__file__).parent/'results'/'quadrants.csv',index=False)
zS=np.where(m0,(S-np.nanmean(S[m0]))/np.nanstd(S[m0]),np.nan)
zC=np.where(m0,(C3-np.nanmean(C3[m0]))/np.nanstd(C3[m0]),np.nan)
rg=np.random.default_rng(6)
print(f"\n正对照:")
for g,lab in ((0.0,'只有加性'),(0.15,'加性 + 已知交互 0.15')):
    y=np.full(NN,np.nan)
    y[m0]=0.2*zS[m0]+0.2*zC[m0]+g*zS[m0]*zC[m0]+rg.standard_normal(int(m0.sum()))
    _,ad,it,se_=quad(y,0.5,0.5); cc=cont(y)
    print(f"   {lab:<22} 四象限交互 **{it:+.4f} ± {se_:.4f}** · 连续交互 **{cc:+.4f}**")
    if g==0.0: c0,i0,s0=cc,it,se_
    else: c1_,i1=cc,it
# ⚠ #300a 的规矩:上页面前先**发明能弄坏它的旋钮**。这里试两个最可能的:
#   ① 控制**勾选数**(广度)及其与两条路的交互 —— 广度若同时驱动两者,可能掩盖交互;
#   ② **分性别** —— #286–#295 说结构按性别不不变,交互可能只存在于一侧而在合并时相消。
NQ=np.zeros(NN); CV=np.zeros(NN)
for M,ppl in MB: CV[ppl]+=1; NQ[ppl]+=M.sum(1)
NQ=np.where(CV>=8,NQ/np.maximum(CV,1),np.nan)
def cont_ctrl(y,mask=None,ctrl=False):
    m=(m0 if mask is None else (m0&mask))&(np.isfinite(NQ) if ctrl else True)
    z=lambda v:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
    cols=[np.ones(m.sum()),z(S),z(C3),z(S)*z(C3)]
    if ctrl: cols+= [z(NQ),z(NQ)*z(S),z(NQ)*z(C3)]
    X=np.column_stack(cols); yy=(y[m]-y[m].mean())/y[m].std()
    return float(np.linalg.lstsq(X,yy,rcond=None)[0][3]),int(m.sum())
SEXV=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
KN=[('原样',cont_ctrl(sh)),('控制勾选数+其交互',cont_ctrl(sh,ctrl=True)),
    ('仅 biomale=0',cont_ctrl(sh,mask=(SEXV==0))),('仅 biomale=1',cont_ctrl(sh,mask=(SEXV==1)))]
print(f"\n发明的旋钮(上页面前的破坏尝试):")
for tag,(v,nn) in KN: print(f"   {tag:<22} 连续交互 **{v:+.4f}**  (n={nn:,})")
kv=[v for _,(v,_) in KN]
print(f"   -> 跨四个口径:**{min(kv):+.4f} … {max(kv):+.4f}**(极差 {max(kv)-min(kv):.4f};"
      f"置换零展布 {np.std(nul):.4f})")

gg=Gate('羞耻的两条路:相加还是相乘')
gg.asserted('★ 发明的旋钮:控制广度 / 分性别,交互有没有冒出来',
            max(abs(v) for v in kv)<2*np.std(nul),
            ' · '.join(f"{t} {v:+.4f}" for t,(v,_) in KN)+f" —— 全部 vs 2×置换零展布 {2*np.std(nul):.4f}")
gg.asserted('★ 正对照①:只有加性的合成结局,四象限交互必须落回零',abs(i0)<2*s0,
            f"交互 {i0:+.4f} vs 2×se {2*s0:.4f} · 连续 {c0:+.4f}")
gg.asserted('★ 正对照②:带已知交互 0.15 的合成结局必须被抓到',abs(c1_)>0.08,
            f"连续交互 {c1_:+.4f}(真值 0.15)· 四象限 {i1:+.4f}")
gg.negative_control('★ 负对照:`perm_finite` 打乱人后的连续交互',float(np.mean(nul)),ci,
    null_spread=float(np.std(nul)),
    null_kind='`perm_finite` 只在有限项内打乱人 —— 题内跨人的零,保住缺失格局(#264b/#278b)')
gg.asserted('★ 注册的 kill:加性 还是 相乘',
            abs(ci-np.mean(nul))<2*np.std(nul),
            f"连续交互 **{ci:+.4f}** vs 置换零 {np.mean(nul):+.4f} ± {np.std(nul):.4f}")
gg.asserted('⚠ 规格曲线:四种切分下交互的符号是否一致(⚠ 真零的签名就是符号乱走 + 全部小于自身 se)',
            len(set(np.sign(T.inter)))==1,
            ' · '.join(f"{r.v_cut} {r.inter:+.4f}±{r.se:.4f}" for _,r in T.iterrows()))
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
