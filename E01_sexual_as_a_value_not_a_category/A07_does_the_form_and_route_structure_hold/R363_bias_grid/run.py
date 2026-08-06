import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A109 R363 -- 去衰减的上偏,是一张表还是一个数

`#317b`:去衰减带约 **+8%** 的上偏 —— 但那是在**一个**真相关(0.20)与**一组**信度上测的,
而项目里的去衰减数字散在 **0.1–0.6** 之间。**上偏很可能随真相关与信度而变。**

⚠ **这是 Closure(工具),如实标注** —— 但它服务于页面上**三个**已发表的去衰减数字。

ESTIMAND        合成网格:真相关 ∈ {0.05, 0.10, 0.20, 0.40, 0.60} × 复合信度 ∈ {0.2, 0.4, 0.6, 0.8, 1.0},
                每格 12 个种子,n 与真实数据同量级(6,478);
                报**绝对上偏**(读数 − 真值)与**相对上偏**(/真值)。
KILL            **若相对上偏在网格上近似恒定 -> 一个乘性系数就够,写进 `CALIBER.md`;
                若随信度剧烈变化 -> 必须查表,而「上偏约一成」这句话不成立。**
POSITIVE CTRL   **信度 1.0 那一列的上偏必须 ≈ 0** —— 否则是流程本身有偏,不是去衰减的问题。
NEGATIVE CTRL   真相关 0 那一行:去衰减读数必须 ≈ 0(上偏在零上没有比例可言,只看绝对)。
IMPOSSIBLE      合成数据是**正态**且指标**平行**;真实指标不是,所以这张表是**下界式的**参考,
                不是对真实上偏的点估计。
"""
import numpy as np, pandas as pd, hashlib
from lib.gates import Gate, check_columns

N=6478; NS=12
def cell(true_r,R,seed):
    rg=np.random.default_rng(seed)
    L=rg.standard_normal(N)
    y=true_r*L+np.sqrt(max(1-true_r**2,1e-12))*rg.standard_normal(N)
    if R>=0.999: a=L.copy(); w=L.copy()
    else:
        r_item=R/(2-R); s=np.sqrt(1/r_item-1)
        a=L+s*rg.standard_normal(N); w=L+s*rg.standard_normal(N)
    za=(a-a.mean())/a.std(); zw=(w-w.mean())/w.std()
    rr=float(np.corrcoef(za,zw)[0,1]); rel=2*rr/(1+rr) if rr>0 else np.nan
    rel=min(rel,1.0)
    f=za+zw; obs=float(np.corrcoef(f,y)[0,1])
    return obs/np.sqrt(max(rel,1e-9)),rel,obs
TRUE=(0.05,0.10,0.20,0.40,0.60); RELS=(0.2,0.4,0.6,0.8,1.0)
rows=[]
print(f"n={N:,} · 每格 {NS} 种子\n")
print("绝对上偏(去衰减读数 − 真值):")
print("真相关\\信度 " + ''.join(f"{R:>10.1f}" for R in RELS))
for t in TRUE:
    line=f"{t:>9.2f}  "
    for R in RELS:
        v=np.array([cell(t,R,4000+1000*int(100*t)+100*int(10*R)+i) for i in range(NS)])
        bias=float(v[:,0].mean()-t)
        rows.append(dict(true_r=t,rel_target=R,rel_meas=float(v[:,1].mean()),
                         read=float(v[:,0].mean()),bias=bias,
                         rel_bias=bias/t if t>0 else np.nan,sd=float(v[:,0].std())))
        line+=f"{bias:>+10.4f}"
    print(line)
T=pd.DataFrame(rows); check_columns(T,'R363')
T.to_csv(pathlib.Path(__file__).parent/'results'/'bias_grid.csv',index=False)
print("\n相对上偏(绝对上偏 / 真值):")
print("真相关\\信度 " + ''.join(f"{R:>10.1f}" for R in RELS))
for t in TRUE:
    line=f"{t:>9.2f}  "
    for R in RELS:
        rb=float(T[(T.true_r==t)&(T.rel_target==R)].rel_bias.iloc[0])
        line+=f"{100*rb:>+9.1f}%"
    print(line)
sub=T[T.rel_target<0.999]
rb_sd=float(sub.rel_bias.std()); rb_mean=float(sub.rel_bias.mean())
ab_sd=float(sub.bias.std())
print(f"\n★ 信度<1 的 {len(sub)} 格:相对上偏 **{100*rb_mean:+.1f}% ± {100*rb_sd:.1f}pp** · "
      f"绝对上偏 sd **{ab_sd:.4f}**")
print(f"   -> {'**相对**' if rb_sd<ab_sd/max(abs(sub.true_r.mean()),1e-9) else '**绝对**'}上偏更稳")
pc=T[T.rel_target>=0.999]
print(f"\n正对照(信度 1.0 那一列):绝对上偏 " + ' · '.join(f"{b:+.4f}" for b in pc.bias) +
      f" -> max |bias| **{float(pc.bias.abs().max()):.4f}**")
ng=np.array([cell(0.0,R,9000+i)[0] for R in RELS for i in range(NS)])
print(f"负对照(真相关 0):去衰减读数 **{ng.mean():+.4f} ± {ng.std():.4f}**")
gg=Gate('去衰减的上偏,是一张表还是一个数')
gg.asserted('★ 正对照:信度 1.0 那一列的上偏必须 ≈ 0',float(pc.bias.abs().max())<0.01,
            f"max |bias| {float(pc.bias.abs().max()):.4f}")
gg.asserted('★ 负对照:真相关 0 时去衰减读数 ≈ 0',abs(ng.mean())<0.01,
            f"{ng.mean():+.4f} ± {ng.std():.4f}")
gg.asserted('★ 注册的 kill:相对上偏在网格上是否近似恒定(sd < 5pp)',rb_sd<0.05,
            f"相对上偏 {100*rb_mean:+.1f}% ± {100*rb_sd:.1f}pp,跨 {len(sub)} 格")
gg.asserted('⚠ 边界:合成数据是正态且指标平行,真实指标不是',True,
            '所以这张表是**参考**,不是对真实上偏的点估计')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
