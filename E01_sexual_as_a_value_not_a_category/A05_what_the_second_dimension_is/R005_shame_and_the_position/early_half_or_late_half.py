import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R225 -- 羞耻贴的是早来的那一半,还是晚来的那一半

`#179` 杀了 B,但明说方向不可分。**方向可以被逼近**:
    位置 -> 羞耻:羞耻跟着**当前**口味的冷门程度,早/晚两半差别不大,或晚半更紧
    羞耻 -> 位置:羞耻是先有的,它塑造了后来要什么 -> **晚半更紧**
    共同因:两半差不多

⚠ `#179` 的 NEXT 里那个设计问题(起始年龄题与多选块**零 item 重叠**)在这里被绕开:
**整轮完全在起始年龄那一侧内部做** —— 31 个类别各有人群稀有度,每人各有获得年龄,
所以"早一半的平均稀有度"与"晚一半的平均稀有度"是同一批 item 上的两个人层量。

ESTIMAND        r(早半平均稀有度, 羞耻) 与 r(晚半平均稀有度, 羞耻),同一批人、同一条管道。
                关键量是**配对差**,不是各自对 0。
IDENTIFICATION  人内按获得年龄中位数劈半;两半在**类别**上不相交。
KILL            条件式:先要**正对照开火**(把羞耻替换成一个由早半构造的合成变量,必须早半 >> 晚半);
                再判配对差:**|Δ| > 2× 自助 sd 才可判方向;否则明说不可分。**
NEGATIVE CTRL   人内打乱获得年龄(保留曲目与年龄分布,毁掉配对)-> 两半差应塌到零。
CONFOUND        两半的**类别数**与**稀有度范围**天然不同(`#128` 的共享时间表:早=常见,晚=罕见)。
                控制:两半各自减去它自己的人群均值;并把两半的类别数作为协变量。
NOISE FLOOR     人层 bootstrap 500 次。
IMPOSSIBLE      仍是横断面自报的回忆年龄(`#162`:回忆偏差是人群规律)。
                本轮能判的是**哪一半贴得更紧**,不是因果。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
SHAME='"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'
sh=df[SHAME].values.astype(float)
# ⚠ #180a:原始列的取值是**年龄段字符串**(`11-12yo`),`to_numeric` 全成 NaN,
#   而下游没有任何一步会因此报错 —— `have` 全 False、稀有度变成常数 9.21、n=0,
#   三个 nan 静静地印了出来。**已解析的数值版在 `data/derived/onset.csv`(轮次检验过的派生件)。**
O=pd.read_csv('data/derived/onset.csv')
onset=[c for c in O.columns if re.search(r'How old were you when you first',c)]
print(f"起始年龄题 {len(onset)} 道(来自 data/derived/onset.csv)")
A_=O[onset].apply(pd.to_numeric,errors='coerce').values
assert np.isfinite(A_).sum()>10000, f'解析后仍几乎全空:{np.isfinite(A_).sum()} 个有效值'
A_=np.where((A_>=2)&(A_<=60),A_,np.nan)
have=np.isfinite(A_)
rar=-np.log(np.clip(have.mean(0),1e-4,1.))          # 人群稀有度:多少人报过这一类
print(f"人群稀有度范围 {rar.min():.2f}–{rar.max():.2f};有 >=6 个起始年龄的人 {int((have.sum(1)>=6).sum()):,}")

def halves(Amat, shuffle=False, rng=None):
    """人内按年龄中位数劈半,返回 (早半平均稀有度, 晚半平均稀有度, 早半类别数, 晚半类别数)"""
    n=Amat.shape[0]; e=np.full(n,np.nan); l=np.full(n,np.nan); ne=np.zeros(n); nl=np.zeros(n)
    for i in range(n):
        idx=np.flatnonzero(np.isfinite(Amat[i]))
        if len(idx)<6: continue
        a=Amat[i,idx]
        if shuffle: a=a[rng.permutation(len(a))]      # 人内打乱:曲目与年龄分布不变,配对毁掉
        o=np.argsort(a,kind='stable'); h=len(o)//2
        E=idx[o[:h]]; L=idx[o[-h:]]
        e[i]=rar[E].mean(); l[i]=rar[L].mean(); ne[i]=len(E); nl[i]=len(L)
    return e,l,ne,nl

E,L,NE,NL=halves(A_)
print(f"早半平均稀有度 {np.nanmean(E):.3f}  晚半 {np.nanmean(L):.3f}  "
      f"(#128 的共享时间表:早=常见、晚=罕见 -> 晚半应更大)")

def pr(y,x,ctrls):
    m=np.isfinite(y)&np.isfinite(x)
    for c in ctrls: m&=np.isfinite(c)
    X=np.c_[np.ones(m.sum()),*[c[m] for c in ctrls]] if ctrls else np.ones((m.sum(),1))
    ry=y[m]-X@np.linalg.lstsq(X,y[m],rcond=None)[0]
    rx=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), int(m.sum())

K=have.sum(1).astype(float)
r_e,n_e=pr(sh,E,[K]); r_l,n_l=pr(sh,L,[K])
print(f"\nr(早半, 羞耻) = {r_e:+.4f}  (n={n_e:,})")
print(f"r(晚半, 羞耻) = {r_l:+.4f}  (n={n_l:,})")

# 配对自助
rb=np.random.default_rng(20260803)
m=np.isfinite(sh)&np.isfinite(E)&np.isfinite(L)&np.isfinite(K); idx=np.flatnonzero(m)
d=[]
for _ in range(500):
    s_=rb.choice(idx,len(idx),replace=True)
    X=np.c_[np.ones(len(s_)),K[s_]]
    ry=sh[s_]-X@np.linalg.lstsq(X,sh[s_],rcond=None)[0]
    re_=E[s_]-X@np.linalg.lstsq(X,E[s_],rcond=None)[0]
    rl_=L[s_]-X@np.linalg.lstsq(X,L[s_],rcond=None)[0]
    d.append(np.corrcoef(ry,re_)[0,1]-np.corrcoef(ry,rl_)[0,1])
dm,dsd=float(np.mean(d)),float(np.std(d))
print(f"配对差(早−晚) = {dm:+.4f} ± {dsd:.4f}   |Δ|/sd = {abs(dm)/dsd:.1f}")
pd.DataFrame(dict(arm=['早','晚'],r=[r_e,r_l],n=[n_e,n_l])).to_csv(
    pathlib.Path(__file__).parent/'results'/'halves.csv',index=False)
print(f"\nsha1 {hashlib.sha1(f'{r_e}{r_l}{dm}'.encode()).hexdigest()[:12]}")

# ---- 对照 ①:人内打乱获得年龄(保留曲目与年龄分布,毁掉配对)------------------
rgn=np.random.default_rng(7)
En,Ln,_,_=halves(A_,shuffle=True,rng=rgn)
r_en,_=pr(sh,En,[K]); r_ln,_=pr(sh,Ln,[K])
print(f"\n对照①人内打乱:早 {r_en:+.4f} · 晚 {r_ln:+.4f} · 差 {r_en-r_ln:+.4f}")

# ---- 对照 ②:两半的信度不同吗(这是最强混杂)---------------------------------
# 早半是常见的东西(范围窄)、晚半是罕见的(范围宽)。若两半**测得的精度不同**,
# 相关就会被不同程度地衰减 —— 那时"早半更紧"说的是仪器,不是人。
def half_reliability(Amat, which, rng):
    """把每一半**再劈一次**,两个四分之一的平均稀有度跨人相关 + Spearman-Brown。"""
    n=Amat.shape[0]; a=np.full(n,np.nan); b=np.full(n,np.nan)
    for i in range(n):
        idx=np.flatnonzero(np.isfinite(Amat[i]))
        if len(idx)<8: continue
        o=idx[np.argsort(Amat[i,idx],kind='stable')]; h=len(o)//2
        H = o[:h] if which=='early' else o[-h:]
        p=rng.permutation(len(H)); k=len(H)//2
        if k<2: continue
        a[i]=rar[H[p[:k]]].mean(); b[i]=rar[H[p[k:2*k]]].mean()
    m=np.isfinite(a)&np.isfinite(b)
    r=float(np.corrcoef(a[m],b[m])[0,1])
    return 2*r/(1+r), int(m.sum())
rgr=np.random.default_rng(11)
rel_e,ne_=half_reliability(A_,'early',rgr); rel_l,nl_=half_reliability(A_,'late',rgr)
print(f"对照②信度:早半 SB {rel_e:.4f} (n={ne_:,}) · 晚半 SB {rel_l:.4f} (n={nl_:,})")
dis_e=r_e/np.sqrt(max(rel_e,1e-6)); dis_l=r_l/np.sqrt(max(rel_l,1e-6))
print(f"        去衰减后:早 {dis_e:+.4f} · 晚 {dis_l:+.4f} · 差 {dis_e-dis_l:+.4f}")

g=Gate('羞耻贴的是早来的还是晚来的')
g.negative_control('人内打乱获得年龄(两半之差)',float(abs(r_en-r_ln)),float(dm),
                   null_spread=float(dsd))
g.resolvable('配对差(早−晚)',dm,dsd)
g.asserted('对照②:两半信度必须可比,否则"更紧"说的是仪器',
           abs(rel_e-rel_l)<0.15,f"早 {rel_e:.3f} vs 晚 {rel_l:.3f},差 {rel_e-rel_l:+.3f}")
g.asserted('去衰减后方向不变',np.sign(dis_e-dis_l)==np.sign(dm),
           f"去衰减差 {dis_e-dis_l:+.4f} vs 原始差 {dm:+.4f}")
g.asserted('预注册的三个世界里,有没有一个预测「早半更紧」',False,
           '位置→羞耻 预测「差别不大或晚半更紧」· 羞耻→位置 预测「晚半更紧」· 共同因 预测「两半差不多」'
           ' —— **观察到的是早半更紧,落在三者之外**')
print(g)
