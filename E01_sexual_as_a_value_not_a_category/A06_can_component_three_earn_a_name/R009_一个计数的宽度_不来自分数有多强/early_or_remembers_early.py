import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A66 R289 -- 「早期积累」还是「把一切都记得更早」

`#243d`:`#243a` 的分组是「**报告**在 14 岁前有 ≥8 个起始年龄」,
它同时携带「真的早」与「记得早」(`#114`:最爱的记得早 −0.2000 年/评分 sd)。

ESTIMAND        **只调整分组,不动 Δ**:把每个类别的起始年龄对**这个人自己对该类别的评分**
                回归,扣掉评分的贡献但**保留年龄尺度**
                (`V_adj = V − b·(rating − mean_rating)`),再用同一条 `≤14 且 ≥8` 规则重新划组,
                重跑 `#243a` 的三行。
KILL            **若反号仍在(早期组明显为负、匹配组明显为正)-> `#243a` 升 D7,
                「早熟」可正式改称「早期积累」;
                若塌掉 -> `#243a` 降为「把一切记得更早的人身上这个关系更强」,
                那是一条关于记忆而不是关于发展的事实。**
⚠ 守卫 12        残差化会改变 `≥8` 的纳入 -> **必须在交集样本上再报一次**(`#239a`)。
POSITIVE CTRL   沿用 `#284` 两端:① 完全由评分驱动的假分组必须被杀;
                ② 与评分无关的假分组必须不动。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      评分是**当下**的,起始年龄是**回忆**;
                只能扣掉「爱得深 -> 记得早」,扣不掉「记得早 -> 现在更爱」。方向不可分。
"""
import numpy as np, pandas as pd, warnings, hashlib, re as _re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
RT=df[rate].apply(pd.to_numeric,errors='coerce').values
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(
    {'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}).notna().sum()>300]
def norm(s): return _re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=_re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
print(f"类别 {Mc};匹配到评分的 {len(best)}({100*len(best)/Mc:.0f}%)")
def adj(Vm, drive=None):
    """扣掉评分贡献但保留年龄尺度。drive 不为 None 时,先造一个由评分驱动的假起始年龄。"""
    out=Vm.copy()
    for j,i in best.items():
        m=np.isfinite(out[:,j])&np.isfinite(RT[:,i])
        if m.sum()<300: continue
        b=np.polyfit(RT[m,i],out[m,j],1)[0]
        out[m,j]=out[m,j]-b*(RT[m,i]-np.nanmean(RT[m,i]))
    return out
def rho_of(Vm):
    D=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=8)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
rho=rho_of(V0); base=np.isfinite(rho); NCAT=np.isfinite(V0).sum(1)
rng=np.random.default_rng(20260804)
def D_of(mask):
    m=mask&base; v=rho[m]
    return (float(np.mean(v)),
            float(np.std([np.mean(v[i]) for i in (rng.choice(len(v),len(v),True) for _ in range(300))])),
            int(m.sum()))
def three_rows(Vgroup, tag, restrict=None):
    E=base&((np.isfinite(Vgroup)&(Vgroup<=14)).sum(1)>=8)
    if restrict is not None: E=E&restrict
    pool=base&(~E)
    if restrict is not None: pool=pool&restrict
    sel=[]
    for c in np.unique(NCAT[E]):
        need=int((NCAT[E]==c).sum()); cand=np.flatnonzero(pool&(NCAT==c))
        if len(cand): sel.append(rng.choice(cand,min(need,len(cand)),replace=False))
    M=np.zeros(N,bool)
    if sel: M[np.concatenate(sel)]=True
    de,dm=D_of(E),D_of(M)
    print(f"  {tag:<22} 早 **{de[0]:+.4f}**±{de[1]:.4f}(n={de[2]:,})· "
          f"匹配 **{dm[0]:+.4f}**±{dm[1]:.4f}(n={dm[2]:,})· 差 **{de[0]-dm[0]:+.4f}**")
    return de,dm,E,M
print(f"\n三行:")
d0=three_rows(V0,'原始分组(`#243a`)')
Va=adj(V0); d1=three_rows(Va,'评分调整后的分组')
common=d0[2]|d0[3]|d1[2]|d1[3]
inter=(d0[2]|d0[3])&(d1[2]|d1[3])
print(f"  ⚠ 守卫 12:交集样本 n = {int((inter&base).sum()):,}")
d0c=three_rows(V0,'原始(交集样本上)',restrict=inter)
d1c=three_rows(Va,'调整后(交集样本上)',restrict=inter)
# 正对照两端
fakeR=V0.copy()
for j,i in best.items():
    m=np.isfinite(V0[:,j])&np.isfinite(RT[:,i]); fakeR[m,j]=20-3.0*RT[m,i]+rng.standard_normal(m.sum())*0.5
pa=three_rows(fakeR,'正对照①假(评分驱动)'); pb=three_rows(adj(fakeR),'正对照①调整后')
fakeI=V0.copy()
for j in range(Mc):
    m=np.isfinite(V0[:,j]); fakeI[m,j]=rng.uniform(2,28,m.sum())
pc=three_rows(fakeI,'正对照②假(与评分无关)'); pd_=three_rows(adj(fakeI),'正对照②调整后')
nul=[]
for _ in range(15):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    r2=rho_of(Vp); nul.append(float(np.nanmean(r2[np.isfinite(r2)&d1[2]])))
print(f"\n置换零(调整后的早期组上){np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(arm='原始',d_early=d0[0][0],d_match=d0[1][0],gap=d0[0][0]-d0[1][0],n_early=d0[0][2]),
                dict(arm='评分调整',d_early=d1[0][0],d_match=d1[1][0],gap=d1[0][0]-d1[1][0],n_early=d1[0][2]),
                dict(arm='原始_交集',d_early=d0c[0][0],d_match=d0c[1][0],gap=d0c[0][0]-d0c[1][0],n_early=d0c[0][2]),
                dict(arm='调整_交集',d_early=d1c[0][0],d_match=d1c[1][0],gap=d1c[0][0]-d1c[1][0],n_early=d1c[0][2])])
check_columns(T,'R289'); T.to_csv(pathlib.Path(__file__).parent/'results'/'adjusted_grouping.csv',index=False)

g=Gate('「早期积累」还是「记得更早」')
g.asserted('正对照①:完全由评分驱动的假分组,调整必须把它杀掉',
           abs(pb[0][0]-pb[1][0])<abs(pa[0][0]-pa[1][0])/2,
           f"{pa[0][0]-pa[1][0]:+.4f} -> {pb[0][0]-pb[1][0]:+.4f}")
g.asserted('正对照②:与评分无关的假分组,调整必须几乎不动',
           abs((pd_[0][0]-pd_[1][0])-(pc[0][0]-pc[1][0]))<0.03,
           f"{pc[0][0]-pc[1][0]:+.4f} -> {pd_[0][0]-pd_[1][0]:+.4f}")
g.negative_control('置换零(调整后的早期组)',abs(float(np.mean(nul))),abs(d1[0][0]),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.control_kept_the_sample('★ 评分调整这个控制',before=d0[0][0]-d0[1][0],after=d1[0][0]-d1[1][0],
                          n_before=d0[0][2],n_after=d1[0][2],
                          before_common=d0c[0][0]-d0c[1][0],after_common=d1c[0][0]-d1c[1][0],
                          n_common=int((inter&base).sum()))
g.asserted('★ 注册的 kill:反号在评分调整后仍在 -> 升 D7,改称「早期积累」',
           d1c[0][0]<-0.05 and d1c[1][0]>-0.02,
           f"交集样本上 调整后 早 {d1c[0][0]:+.4f} · 匹配 {d1c[1][0]:+.4f} · 差 {d1c[0][0]-d1c[1][0]:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
