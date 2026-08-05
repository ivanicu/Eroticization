import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A151 R442 -- 换一个结局去扫:「治疗性」

`#397b`:页面上前三条羞耻路径都是我**去找羞耻**时找到的,第四条是**在找别的东西时**掉出来的。
**⇒ 换一个结局去扫,等于换一组假设** —— 而这个项目至今只在**两个**结局上扫过。

第三个结局:`"Engaging with or fantasizing about what arouses me feels therapeutic or healing to me"`
**「治疗性」** —— 页面上**从没被当作结局用过**,而它是羞耻的**反面候选**,不是它的同义词。

两个活着的世界:
**A 治疗性 = 羞耻取反** -> 去掉羞耻后,残差里**没有**东西越过族内零;
**B 治疗性是自己的一条线** -> **有**,而且**不是**已知的那几条(`S` · `c3⁻` · `EARLY` · 五题)。

ESTIMAND        `治疗性` 的残差(去掉 **羞耻** · `S` · `c3⁻` · 类别数)对**所有**数值型人层量求相关;
                主量 = **最大 |r|**(对**族内取最大零**)+ 逐个的名字。
判据(**先标支**,`#379c`)
                【两支】**阳性参照必须先开火**(`#419`):用**羞耻自己**,在**未去掉它**的残差上,
                        对**单变量零**(`#388b` 的类别错误不再犯;`#391c` 的恒等式陷阱不再犯);
                        负对照用**越阈率**(`#395b`);guard 26 用 **MDE 扫描**。
                【非零支】最大 |r| 越过族内零 -> 世界 B,报名字;
                          并**逐个核对**它是不是已知的那四条之一(`#357b`:相关高可能只是它已在页面上)。
                【零支】未越阈 -> 世界 A,启用 MDE。
⚠ 零的种类     `offset_control`:**最大相关的零绝不是零**(数百候选取最大,最大值天然为正)
                -> 零 = **同样多的合成噪声列取最大**(`#419`/`#432` 已验证的方法)。
⚠ `#392e`      结局进模型前先打印它的**取值集合、众数、与方向已知锚的相关**。
IMPOSSIBLE      ① 残差去掉四个量 -> 与它们共线的候选被系统性压低,**这不是「它不相关」**;
                ② 全自报,同源方差;③ 扫描是**发现**不是**检验** —— 命中的要另一轮独立复制。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
TC=next(c for c in d.columns if 'vmq8jqw' in c)
THER=pd.to_numeric(d[TC],errors='coerce').values.astype(float)
# ---- ⚠ `#392e`:结局进模型前先看它自己 ----
u=np.unique(THER[np.isfinite(THER)])
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
ga=np.isfinite(THER)&np.isfinite(anc)
print(f"⚠ **`#392e`:结局进模型前先看它自己**")
print(f"   `{str(TC)[:60]}`")
print(f"   取值集合 = {u.tolist()} · 众数 = **{float(pd.Series(THER[np.isfinite(THER)]).mode().iloc[0]):g}** · "
      f"n = {int(np.isfinite(THER).sum()):,}")
print(f"   与方向已知锚 `Totalsexacts` 相关 = **{np.corrcoef(THER[ga],anc[ga])[0,1]:+.4f}**")
print(f"   与羞耻的相关 = **{np.corrcoef(THER[np.isfinite(THER)&np.isfinite(sh)],sh[np.isfinite(THER)&np.isfinite(sh)])[0,1]:+.4f}** "
      f"-> ⚠ **若这个数接近 −1,世界 A 就成立而不必扫**\n")
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(THER)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
CTRL=[sh,S,C3,ncat]
X0=np.column_stack([np.ones(n)]+[z(v,M) for v in CTRL]); yT=z(THER,M)
b0,*_=np.linalg.lstsq(X0,yT,rcond=None); RES=yT-X0@b0
X1=np.column_stack([np.ones(n)]+[z(v,M) for v in [S,C3,ncat]])
b1,*_=np.linalg.lstsq(X1,yT,rcond=None); RES_NS=yT-X1@b1     # ⚠ **未去掉羞耻**的残差(`#391c` 的教训)
print(f"n=**{n:,}**")
CAND=[]
for c in d.columns:
    if c==TC: continue
    v=pd.to_numeric(d[c],errors='coerce').values.astype(float)
    if np.isfinite(v[M]).sum()>=2000 and np.nanstd(v[M])>1e-9: CAND.append((c,v))
K=len(CAND)
def cr(res,v):
    g=np.isfinite(v[M])
    return float(np.corrcoef(res[g],z(v,M)[g])[0,1]) if g.sum()>=2000 else np.nan
cor=sorted(((abs(cr(RES,v)),cr(RES,v),c) for c,v in CAND if np.isfinite(cr(RES,v))),reverse=True)
print(f"候选 **{K}** 个 · 最大 |r| = **{cor[0][0]:.4f}**")
for k,(a,r_,c) in enumerate(cor[:8]): print(f"   {k+1}. {r_:+.4f}  {str(c)[:66]}")
offmax=[]
for s in range(200):
    rg=np.random.default_rng(3300+s); mx=0.
    for _ in range(K): mx=max(mx,abs(float(np.corrcoef(RES,rg.standard_normal(n))[0,1])))
    offmax.append(mx)
offmax=np.array(offmax); OTHR=float(np.percentile(offmax,95))
print(f"⚠ offset 零(**同样多({K})的合成噪声列取最大**,200 次)= "
      f"**{offmax.mean():.4f} ± {offmax.std():.4f}** · 95 分位 **{OTHR:.4f}**")
print(f"   -> 实测最大 {cor[0][0]:.4f} {'**越阈**' if cor[0][0]>OTHR else '**未越阈**'} · "
      f"越阈候选 **{sum(1 for a,_,_ in cor if a>OTHR)}/{len(cor)}**")
nulS=np.array([abs(float(np.corrcoef(RES_NS,perm_in(sh,M,7100+s)[M])[0,1])) for s in range(400)])
STHR=float(np.percentile(nulS,95)); rS=abs(cr(RES_NS,sh))
print(f"\n★ 阳性参照(**羞耻自己**,在**未去掉它**的残差上 —— `#391c` 的恒等式陷阱不再犯):"
      f"|r| = **{rS:.4f}** vs **单变量零** {STHR:.4f} -> {'**开火**' if rS>STHR else '**不开火**'}")
negs=np.array([abs(float(np.corrcoef(RES,np.random.default_rng(9200+s).standard_normal(n))[0,1]))
               for s in range(200)])
rate=float((negs>OTHR).mean())
print(f"  负对照(**越阈率**,单列噪声 vs 族内零):**{100*rate:.1f}%** —— "
      f"⚠ 族内零是**最大值**的零,单列本就该几乎不越,合格 ≤2%")
print(f"\nguard 26 的正对照 = **MDE 扫描**,每级 30 次:")
MDE=None
for gg in (0.03,0.05,0.08,0.12):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(1800+int(gg*100)*61+s_)
        v=gg*RES+np.sqrt(max(1-gg*gg,1e-9))*rg.standard_normal(n)
        if abs(float(np.corrcoef(RES,v)[0,1]))>OTHR: hit+=1
    print(f"   种植 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.15
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 = 实测最大 **{cor[0][0]:.4f}**")
T=pd.DataFrame([dict(v_col=str(c)[:70],v_r=r_,v_absr=a) for a,r_,c in cor[:30]])
check_columns(T,'R442'); T.to_csv(pathlib.Path(__file__).parent/'results'/'ther_corr.csv',index=False)
g=Gate('「治疗性」是羞耻取反,还是自己的一条线')
g.asserted('★【两支】阳性参照:**羞耻自己**必须在单变量零上开火',rS>STHR,
           f"{rS:.4f} vs {STHR:.4f}",kind='control')
g.asserted('★【两支】负对照:单列噪声对**族内**零的越阈率 ≤2%',rate<=0.02,
           f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 实测最大',MDE_,float(cor[0][0]),True,what='MDE 扫描 80% 检出')
g.asserted('★【两支】offset 零非退化(取最大天然为正)',offmax.std()>0,
           f"{offmax.mean():.4f} ± {offmax.std():.4f}",kind='control')
if rS>STHR and rate<=0.02:
    g.asserted('★【非零支】最大 |r| 越过族内零 -> 世界 B(治疗性是自己的一条线)',cor[0][0]>OTHR,
               f"{cor[0][0]:.4f} vs {OTHR:.4f} · 头名 {str(cor[0][2])[:44]}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
