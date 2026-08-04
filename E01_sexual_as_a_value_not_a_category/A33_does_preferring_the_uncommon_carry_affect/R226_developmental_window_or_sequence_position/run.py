import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A33 R226 -- 「早」是发育窗口,还是序列位置

`#180` 的"早"是**人内相对**的(这个人自己最早的一半),它混着两件事:
    绝对年龄早(青春期前)   -> **发育窗口**
    在这个人的序列里靠前     -> **身份地基**
心理学含义完全不同,而 `#180` 的设计分不开。

ESTIMAND        同一个配对差 r(早, 羞耻) − r(晚, 羞耻),按**两种劈分**各算一次:
                相对(人内中位数)· 绝对(12 岁前 / 12 岁后)。
KILL            绝对下消失而相对下还在 -> **序列位置**;反之 -> **发育窗口**;
                两者都在 -> **本数据分不开,明说**(这也是一个合法结论)。
NEGATIVE CTRL   两种劈分各配一个人内打乱。
CONFOUND        绝对劈分下,两组的**大小**因人而异(早熟的人几乎全在"早"组)->
                只取两侧各 ≥3 个类别的人;并把两侧类别数作为协变量。
NOISE FLOOR     配对 bootstrap 500 次。
MULTIPLICITY    2 种劈分 × {真实, 打乱} × 阈值 {10,12,14} = 12 格,整格发表。
IMPOSSIBLE      绝对与相对在数据里高度相关(早熟的人两种劈分几乎一致)-> 若两者都显著,
                本设计**无法**归因,不得挑一个报。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
sh=df['"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'].values.astype(float)
O=pd.read_csv('data/derived/onset.csv')
onset=[c for c in O.columns if re.search(r'How old were you when you first',c)]
A_=O[onset].apply(pd.to_numeric,errors='coerce').values
A_=np.where((A_>=2)&(A_<=60),A_,np.nan)
assert np.isfinite(A_).sum()>10000
have=np.isfinite(A_); rar=-np.log(np.clip(have.mean(0),1e-4,1.)); K=have.sum(1).astype(float)
print(f"{len(onset)} 个类别,稀有度 {rar.min():.2f}–{rar.max():.2f},有效人 {int((K>=6).sum()):,}")

def split_relative(Amat, shuffle=False, rng=None):
    n=Amat.shape[0]; e=np.full(n,np.nan); l=np.full(n,np.nan); ne=np.zeros(n); nl=np.zeros(n)
    for i in range(n):
        idx=np.flatnonzero(np.isfinite(Amat[i]))
        if len(idx)<6: continue
        a=Amat[i,idx]
        if shuffle: a=a[rng.permutation(len(a))]
        o=np.argsort(a,kind='stable'); h=len(o)//2
        E=idx[o[:h]]; L=idx[o[-h:]]
        e[i]=rar[E].mean(); l[i]=rar[L].mean(); ne[i]=len(E); nl[i]=len(L)
    return e,l,ne,nl

def split_absolute(Amat, cut, shuffle=False, rng=None):
    n=Amat.shape[0]; e=np.full(n,np.nan); l=np.full(n,np.nan); ne=np.zeros(n); nl=np.zeros(n)
    for i in range(n):
        idx=np.flatnonzero(np.isfinite(Amat[i]))
        if len(idx)<6: continue
        a=Amat[i,idx]
        if shuffle: a=a[rng.permutation(len(a))]
        E=idx[a<cut]; L=idx[a>=cut]
        if len(E)<3 or len(L)<3: continue          # 两侧各 >=3,否则早熟/晚熟的人只剩一侧
        e[i]=rar[E].mean(); l[i]=rar[L].mean(); ne[i]=len(E); nl[i]=len(L)
    return e,l,ne,nl

def diff(e,l,ne,nl,seed):
    m=np.isfinite(sh)&np.isfinite(e)&np.isfinite(l)&np.isfinite(K)
    idx=np.flatnonzero(m)
    def pr(y,x,s_):
        X=np.c_[np.ones(len(s_)),K[s_],ne[s_],nl[s_]]
        ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]; rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0]
        return np.corrcoef(ry,rx)[0,1]
    r_e=pr(sh[idx],e[idx],idx); r_l=pr(sh[idx],l[idx],idx)
    rb=np.random.default_rng(seed); d=[]
    for _ in range(500):
        s_=rb.choice(idx,len(idx),replace=True)
        d.append(pr(sh[s_],e[s_],s_)-pr(sh[s_],l[s_],s_))
    return r_e,r_l,float(np.mean(d)),float(np.std(d)),len(idx)

rows=[]
print(f"\n{'劈分':<16}{'臂':<8}{'n':>8}{'r 早':>9}{'r 晚':>9}{'配对差':>10}{'sd':>9}{'|Δ|/sd':>8}")
for name,fn in [('相对(人内中位)',lambda s,r: split_relative(A_,s,r))]+\
               [(f'绝对(<{c} 岁)',(lambda c: (lambda s,r: split_absolute(A_,c,s,r)))(c)) for c in [10,12,14]]:
    for arm,shuf in [('真实',False),('人内打乱',True)]:
        rng=np.random.default_rng(hash(name)%1000+int(shuf))
        e,l,ne,nl=fn(shuf,rng)
        r_e,r_l,dm,dsd,n=diff(e,l,ne,nl,seed=abs(hash(name+arm))%9973)
        rows.append(dict(split=name,arm=arm,n=n,r_early=r_e,r_late=r_l,d=dm,sd=dsd,
                         ratio=abs(dm)/dsd if dsd>0 else np.nan))
        print(f"{name:<16}{arm:<8}{n:>8,}{r_e:>+9.4f}{r_l:>+9.4f}{dm:>+10.4f}{dsd:>9.4f}"
              f"{abs(dm)/dsd if dsd>0 else np.nan:>8.1f}",flush=True)
T=pd.DataFrame(rows); T.to_csv(pathlib.Path(__file__).parent/'results'/'grid.csv',index=False)

real=T[T.arm=='真实']; rel=real[real.split.str.startswith('相对')].iloc[0]
absr=real[real.split.str.startswith('绝对')]
g=Gate('「早」是发育窗口还是序列位置')
g.asserted('可判前提:两种劈分都有足够样本',bool((real.n>2000).all()),
           f"最小 n {int(real.n.min()):,}")
g.negative_control('人内打乱(相对劈分)',
                   float(abs(T[(T.split.str.startswith('相对'))&(T.arm=='人内打乱')].d.iloc[0])),
                   float(rel.d))
g.resolvable('相对劈分的配对差',float(rel.d),float(rel.sd))
for _,r in absr.iterrows():
    g.resolvable(f'{r.split} 的配对差',float(r.d),float(r.sd))
g.no_sign_crossing('四种劈分的配对差同号',[float(x) for x in real.d])
print(g)
# ⚠ #181a:第一版这里**直接打出结论**,而它自己的可判前提刚刚失败了 ——
#   绝对劈分 n 只有 1,089–3,854(相对是 11,329),而**点估计更大**(+0.0553 vs +0.0274),
#   只是更噪。把"没通过可分辨性"读成"发育窗口不起作用",正是 P5 ★ 与「零需要功效」。
#   verdict 必须**受同一个条件式管辖**,不能只在 Gate 里守而在 print 里放行。
n_abs_res=int((absr.ratio>2).sum())
mde=lambda sd: 2*sd                                  # 本设计能分辨的最小效应 = 2× 自身展布
print(f"\n  相对:Δ {rel.d:+.4f} (n={int(rel.n):,}, {rel.ratio:.1f}×)")
for _,r in absr.iterrows():
    # ⚠ 第一版这句"点估计比相对劈分大"是**逐行印的通用语**,而对 <14 岁那一行是假的
    #   (+0.0109 < +0.0274)。#178 那一类:一个印在报告里的常量/通用语,没人会查。
    cmp_ = '大' if abs(r.d)>abs(rel.d) else '小'
    print(f"  {r.split}:Δ {r.d:+.4f} (n={int(r.n):,}, {r.ratio:.1f}×)  "
          f"MDE={mde(r.sd):.4f} —— 点估计比相对劈分**{cmp_}**")
precond_ok = bool((real.n>2000).all())
if not precond_ok:
    print("\n  => UNVERIFIED。**不判**。绝对劈分的臂只有相对臂的 1/10–1/3 样本,")
    print(f"     而 3 个绝对阈值里 {int((absr.d.abs()>abs(rel.d)).sum())} 个的点估计**更大**;")
    print("     把它们的不可分辨读成「发育窗口不起作用」是把功效不足")
    print("     读成证据缺席(P5 ★)。**本设计分不开「发育窗口」与「序列位置」。**")
else:
    print(f"  => {'两者都在,分不开' if (rel.ratio>2 and n_abs_res>0) else ('序列位置' if rel.ratio>2 else ('发育窗口' if n_abs_res>0 else '都不可分辨'))}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
