import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A123 R382 -- 到人×类别那一层问:重叠那段是不是「在冷门的东西上比别人更早」

`#336b`:人层扣不掉时间表(代数上是减一个常数)。**必须在人×类别层做。**

⚠⚠ **本轮最大的陷阱,先写下来**:一个人贡献**多行**而羞耻只有**一个**值 ->
**标准误必须按人聚类**,否则 n 被虚增约 30 倍。
**本轮不用解析标准误,只用「打乱人」的置换零** —— 那自动是人层的。

ESTIMAND        长表 `(人 i, 类别 j)`:因变量 = 该人的羞耻(人层,重复);
                预测量 = `dev_ij`(该格起始相对该类别中位的偏差,**正 = 比别人早**)·
                `rar_j`(该类别的稀有度 = −log 流行度)· **交互 `dev × rar`**。
KILL            **若 `dev × rar` 在格层出现 -> 重叠那段是「在**冷门**的东西上比别人更早」;
                若不出现 -> 重叠是人层的,格层看不到它。**
POSITIVE CTRL   在格层种入 `dev × rar` -> 必须被抓到。
NEGATIVE CTRL   **打乱人**(保持每个人的格结构)—— 这是本轮唯一合法的零。
IMPOSSIBLE      羞耻是**人层**的一个数,所以格层回归的因变量在人内是常数;
                这限制了格层能解释的方差上限。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')

BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
NC=ONS.shape[1]
MED=np.nanmedian(ONS,0); PREV=np.isfinite(ONS).mean(0)
RAR=-np.log(np.clip(PREV,1e-4,1.))
nc=np.isfinite(ONS).sum(1)
keep_p=np.flatnonzero((nc>=5)&np.isfinite(sh))
DEV=MED[None,:]-ONS                                   # 正 = 比人群典型更早
ii,jj=np.nonzero(np.isfinite(DEV[keep_p]))
pid=keep_p[ii]; dev=DEV[keep_p][ii,jj]; rar=RAR[jj]; ycell=sh[pid]
NROW=len(pid); NPERS=len(np.unique(pid))
print(f"长表 **{NROW:,}** 行 · **{NPERS:,}** 人 · {NC} 类别 · 平均每人 **{NROW/NPERS:.1f}** 行")
print(f"⚠⚠ 一个人贡献多行而羞耻只有一个值 -> **只用「打乱人」的置换零**,不用解析标准误")
z=lambda v:(v-v.mean())/max(v.std(),1e-12)
zd,zr=z(dev),z(rar)
X=np.column_stack([np.ones(NROW),zd,zr,zd*zr])
def beta(yy):
    b,*_=np.linalg.lstsq(X,(yy-yy.mean())/max(yy.std(),1e-12),rcond=None); return b[1:]
b0=beta(ycell)
print(f"\n格层系数:`dev` **{b0[0]:+.5f}** · `rar` **{b0[1]:+.5f}** · "
      f"**`dev × rar` {b0[2]:+.5f}**")
uniq=np.unique(pid); pos={p:k for k,p in enumerate(uniq)}
pidx=np.array([pos[p] for p in pid])
shp=sh[uniq]
rg=np.random.default_rng(606); NP=400
NUL=np.array([beta(shp[rg.permutation(len(uniq))][pidx]) for _ in range(NP)])
for k,nm in enumerate(['dev','rar','**dev × rar**']):
    q=float((np.abs(NUL[:,k])>=abs(b0[k])).mean())
    print(f"   {nm:<14} 零 **{NUL[:,k].mean():+.5f} ± {NUL[:,k].std():.5f}** · "
          f"|零| ≥ |观测| **{q:.3f}** · **{abs(b0[k]-NUL[:,k].mean())/max(2*NUL[:,k].std(),1e-12):.2f}× 的 2×展布**")
qint=float((np.abs(NUL[:,2])>=abs(b0[2])).mean())
rgp=np.random.default_rng(77)
for g in (0.0,0.02,0.05):
    yp=shp+0.0
    cell=g*(zd*zr)+rgp.standard_normal(NROW)
    per=np.zeros(len(uniq))
    np.add.at(per,pidx,cell); cnt=np.bincount(pidx,minlength=len(uniq))
    yy=(per/np.maximum(cnt,1))[pidx]
    bb=beta(yy)
    print(f"正对照 g={g:.2f}(格层种入 `dev × rar`,再聚合回人层):"
          f"`dev × rar` **{bb[2]:+.5f}** · "
          f"{abs(bb[2]-NUL[:,2].mean())/max(2*NUL[:,2].std(),1e-12):.2f}× 的 2×展布")
    if g==0.05: pc=bb[2]
    if g==0.0: pc0=bb[2]
T=pd.DataFrame([dict(v_term='dev',v_b=float(b0[0])),dict(v_term='rar',v_b=float(b0[1])),
                dict(v_term='dev×rar',v_b=float(b0[2])),dict(v_term='零sd',v_b=float(NUL[:,2].std()))])
check_columns(T,'R382'); T.to_csv(pathlib.Path(__file__).parent/'results'/'cell.csv',index=False)
gg=Gate('格层:重叠那段是不是「在冷门的东西上比别人更早」')
gg.asserted('★ 正对照:格层种入 `dev × rar` 必须被抓到,且 g=0 落零',
            abs(pc-NUL[:,2].mean())>2*NUL[:,2].std() and abs(pc0-NUL[:,2].mean())<2*NUL[:,2].std(),
            f"g=0.05 -> {pc:+.5f} · g=0 -> {pc0:+.5f} · 零 sd {NUL[:,2].std():.5f}")
gg.negative_control('★ 负对照:**打乱人**(保持格结构)后的 `dev × rar`',
    float(NUL[:,2].mean()),float(b0[2]),null_spread=float(NUL[:,2].std()),
    null_kind='打乱**人**、保持每个人的格结构 —— 唯一合法的人层零')
gg.asserted('★ 注册的 kill:`dev × rar` 在格层出不出现',
            abs(b0[2]-NUL[:,2].mean())>2*NUL[:,2].std(),
            f"**{b0[2]:+.5f}** vs 零 {NUL[:,2].mean():+.5f} ± {NUL[:,2].std():.5f} "
            f"(|零| ≥ |观测| **{qint:.3f}**)")
gg.asserted('⚠⚠ 聚类:只用打乱人的置换零,不用解析标准误',True,
            f"{NROW:,} 行来自 {NPERS:,} 人 —— 解析 se 会把 n 虚增 {NROW/NPERS:.0f} 倍")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
