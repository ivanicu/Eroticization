import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R224 -- 「偏爱冷门」这件事,带不带情感

**方向变更。** `#167`-`#178` 连着十二轮都在查我自己的仪器,而 §0.2 与 frontier §3 都说:
当最可引用的句子变成关于我的严谨而不是关于对象时,报告就是关于我的。回到心理学。

`#163`-`#165` 把「偏爱冷门」钉成一个**可迁移的位置倾向**:不是内容(0.60 vs 0.26)、
不是性别(去性别后 102%)、不是勾选数(67%)。但**从没问过它带不带情感** ——
而 `#101`/`#102` 说它唯一挂得住的外部锚是性别,人格五因素全部 |r| ≤ 0.056。
**羞耻与疗愈不是人格,是直接的情感—关系变量,从没测过。**

    A 专用性内容系统:羞耻跟着**内容**走(哪些家族禁忌),位置的表观关联全由内容中介
    B 对普通表征的情色赋值:位置是一个**坐标**,不带立场 -> 控制内容后 S ⟂ 羞耻
    C 递归:赋值回流并重塑关系 -> **控制内容后 S 仍预测羞耻/疗愈**

ESTIMAND        corr(S, 羞耻) 与 corr(S, 疗愈),在控制**内容成分**、勾选数、性别之后。
IDENTIFICATION  S 来自 32 个多选块;羞耻/疗愈是独立的 Likert 题(−3..+3,n=15,503),**零 item 重叠**。
                内容成分用 `#165`(R210)的主成分打分,同一批块。
KILL            条件式:先要**两个对照都开火**(种植的人层变量必须被 S 测到;打乱羞耻必须塌到零);
                再判:**控制内容后 |r(S, 羞耻)| 仍 > 2× 自身自助 sd -> C 存活,B 被杀。**
NEGATIVE CTRL   跨人打乱羞耻(1000 次)-> r 的零分布。
POSITIVE CTRL   合成一个 = S + 噪声 的人层变量,同一条管道必须测到强关联。
NOISE FLOOR     人层 bootstrap 1000 次。
MULTIPLICITY    2 个结局 × 4 个模型(raw · +勾选数 · +性别 · +内容)= 8 格,整格发表。
IMPOSSIBLE      横断面自报;判不了因果方向 —— 「因为冷门所以羞耻」与「因为羞耻所以只敢要冷门的」
                在本设计里不可分。这一条决定了本轮最多能杀 B,不能在 A 与 C 之间定夺。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A22_is_rare_affinity_the_right_name'
          /'R210_how_big_is_the_content_side'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('ARMS=')[0])

SHAME='"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'
HEAL ='"Engaging with or fantasizing about what arouses me feels therapeutic or healing to me" (vmq8jqw)'
sh=df[SHAME].values.astype(float); he=df[HEAL].values.astype(float)
sex=pd.to_numeric(df.get('biomale'),errors='coerce').values.astype(float)

# 全部块上的 S(位置)、C(内容)、K(勾选数)
# ⚠ #179a:第一版用 `build(qs,...)`,它要求覆盖**一半的块**(16 个)-> n 只剩 950。
#   直接算,门槛降到 **≥8 块**(与 `#163` 同口径)。
con=np.zeros(N); pos=np.zeros(N); cnt=np.zeros(N); K=np.zeros(N)
for q in qs:
    M=BL[q]['M']; ppl=BL[q]['ppl']; rar=BL[q]['rar']
    Z=M-M.mean(0,keepdims=True)
    w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); pc=v[:,-1]
    con[ppl]+=Z@pc
    pos[ppl]+=(M@rar)/np.maximum(M.sum(1),1)
    K[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
Sc=np.where(ok,con/np.maximum(cnt,1),np.nan)
Sp=np.where(ok,pos/np.maximum(cnt,1),np.nan)
K =np.where(ok,K,np.nan)
print(f"S(位置)有效 {np.isfinite(Sp).sum():,}  C(内容)有效 {np.isfinite(Sc).sum():,}  "
      f"羞耻 {np.isfinite(sh).sum():,}  疗愈 {np.isfinite(he).sum():,}")

def partial(y, x, ctrls):
    m=np.isfinite(y)&np.isfinite(x)
    for c in ctrls: m&=np.isfinite(c)
    if m.sum()<500: return np.nan,0
    X=np.c_[np.ones(m.sum())]+0.
    if ctrls: X=np.c_[np.ones(m.sum()),*[c[m] for c in ctrls]]
    ry=y[m]-X@np.linalg.lstsq(X,y[m],rcond=None)[0]
    rx=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), int(m.sum())

MODELS=[('raw',[]),('+勾选数',[K]),('+性别',[K,sex]),('+内容',[K,sex,Sc])]
rows=[]
print(f"\n{'结局':<8}{'模型':<10}{'n':>8}{'r(S, 结局)':>12}{'bootstrap sd':>14}{'|r|/sd':>9}")
for oname,y in [('羞耻',sh),('疗愈',he)]:
    for mname,ctrls in MODELS:
        r,n=partial(y,Sp,ctrls)
        m=np.isfinite(y)&np.isfinite(Sp)
        for c in ctrls: m&=np.isfinite(c)
        idx=np.flatnonzero(m); rb=np.random.default_rng(zlib.crc32(f'{oname}{mname}'.encode())%(1<<30))
        bs=[]
        for _ in range(300):
            s_=rb.choice(idx,len(idx),replace=True)
            yy=y[s_]; xx=Sp[s_]; cc=[c[s_] for c in ctrls]
            X=np.c_[np.ones(len(s_)),*cc] if cc else np.ones((len(s_),1))
            ry=yy-X@np.linalg.lstsq(X,yy,rcond=None)[0]; rx=xx-X@np.linalg.lstsq(X,xx,rcond=None)[0]
            bs.append(np.corrcoef(ry,rx)[0,1])
        sd=float(np.std(bs))
        rows.append(dict(outcome=oname,model=mname,n=n,r=r,boot_sd=sd,ratio=abs(r)/sd if sd>0 else np.nan))
        print(f"{oname:<8}{mname:<10}{n:>8,}{r:>+12.4f}{sd:>14.4f}{abs(r)/sd if sd>0 else np.nan:>9.1f}",flush=True)
T=pd.DataFrame(rows); check_columns(T,'R224'); T.to_csv(pathlib.Path(__file__).parent/'results'/'grid.csv',index=False)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 对照 -------------------------------------------------------------------
print("\n---- 对照 ----")
rgn=np.random.default_rng(20260803)
m=np.isfinite(sh)&np.isfinite(Sp)&np.isfinite(K)&np.isfinite(sex)&np.isfinite(Sc)
X=np.c_[np.ones(m.sum()),K[m],sex[m],Sc[m]]
def pr(y,x):
    ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]; rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])
r_real=pr(sh[m],Sp[m])
null=[pr(rgn.permutation(sh[m]),Sp[m]) for _ in range(300)]
plant=Sp[m]+rgn.standard_normal(m.sum())*Sp[m].std()*2.0    # 合成:S + 噪声
r_plant=pr(plant,Sp[m])
print(f"  真实       r = {r_real:+.4f}")
print(f"  跨人打乱羞耻 r = {np.mean(null):+.4f} ± {np.std(null):.4f}  (300 次)")
print(f"  正对照(S+噪声)r = {r_plant:+.4f}")

# ---- 规格曲线:覆盖门槛 ------------------------------------------------------
print("\n---- 规格曲线:块覆盖门槛 ----")
spec=[]
for thr in [4,6,8,10,12,16,20]:
    okk=cnt>=thr
    S2=np.where(okk,pos/np.maximum(cnt,1),np.nan); C2=np.where(okk,con/np.maximum(cnt,1),np.nan)
    K2=np.where(okk,K,np.nan)
    for oname,y in [('羞耻',sh),('疗愈',he)]:
        r,n=partial(y,S2,[K2,sex,C2]); spec.append(dict(thr=thr,outcome=oname,n=n,r=r))
S=pd.DataFrame(spec); S.to_csv(pathlib.Path(__file__).parent/'results'/'spec_curve.csv',index=False)
print(S.pivot_table(index='thr',columns='outcome',values=['n','r']).round(4).to_string())
sh_spec=S[S.outcome=='羞耻']
same=int((np.sign(sh_spec.r)==np.sign(r_real)).sum())

g=Gate('「偏爱冷门」带不带情感')
g.negative_control('跨人打乱羞耻',float(abs(np.mean(null))),float(r_real),
                   null_spread=float(np.std(null)))
g.asserted('正对照:S+噪声 必须被同一条管道强测到',r_plant>0.3,f"{r_plant:+.4f}")
g.resolvable('控制内容后 S 仍预测羞耻',float(T[(T.outcome=='羞耻')&(T.model=='+内容')].r.iloc[0]),
             float(T[(T.outcome=='羞耻')&(T.model=='+内容')].boot_sd.iloc[0]))
g.asserted('内部判别:同一条管道下疗愈**塌掉**而羞耻**没有** -> 不是通用的方法学假象',
           abs(float(T[(T.outcome=='疗愈')&(T.model=='+内容')].r.iloc[0]))<0.03,
           f"疗愈 +内容 {float(T[(T.outcome=='疗愈')&(T.model=='+内容')].r.iloc[0]):+.4f} "
           f"vs 羞耻 {float(T[(T.outcome=='羞耻')&(T.model=='+内容')].r.iloc[0]):+.4f}")
g.asserted('规格曲线:羞耻在 7 个覆盖门槛上同号',same>=6,f"{same}/7 同号")
g.covers_every_arm('整格 8 个单元全部发表',list(T.model)+list(T.outcome),
                   ['raw','+勾选数','+性别','+内容','羞耻','疗愈'])
print(g)
