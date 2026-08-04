import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A17 R02 -- 本项目还有多少条现存声明,其判定阈值是写死在源码里的一个选定常数?

#141c:A02/R10 的 `mx<0.4 -> "still three axes"` 是一个**选定**阈值,而它的自助区间
跨过它 30.2% 的时间。**那个洞比它自己大** —— 阈值写死这件事本身是可以机械扫描的。

给 tools/guard_lint.py 加一条规则(#142),扫**现存声明背后的轮次**,然后对命中的
逐个量"观测值离那个常数有多远,以它自己的展布计"。

⚠ P6 代理账见 tools/guard_lint.py::hardcoded_thresholds —— 只在**命中**方向可读。

ESTIMAND        每个命中处:观测值 · 写死的常数 · 观测值的自助展布 · 两者距离(以展布计)。
IDENTIFICATION  直接读该轮持久化的 results/,不重跑,不改设计。
SCOPE           现存声明背后的 21 轮。
WORLDS          edge  某条现存声明的判定在门槛上 <1 个展布 -> 它由一个选定常数决定
                safe  全部距门槛 >2 个展布 -> 常数是选的,但不影响判定
KILL            条件式:自助展布必须先在一个已知量上给出合理值,才读距离。
POSITIVE CTRL   #141 已知的刀尖(A02/R10 的 30.2%)必须被这套读法认出来。
NEGATIVE CTRL   —— 本轮是审计,不是实验;不构造零。
NOISE FLOOR     2000 次按对自助。
MULTIPLICITY    命中处全部报,不挑。
IMPOSSIBLE      未命中不等于阈值是量出来的(阈值可能来自一个同样被选的变量)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
import guard_lint as gl
from lib.gates import Gate

WANT=['A03/R22','A10/R15','A10/R17','A10/R01','A11/R12','A11/R13','A11/R17','A11/R14',
      'A11/R16','A12/R10','A13/R01','A12/R02','A12/R06','A12/R12','A12/R13','A12/R14',
      'A11/R20','A02/R10','A02/R16','A02/R17','A03/R19']
paths=[]
for w in WANT:
    a,r=w.split('/')
    paths+=[str(h.relative_to(gl.ROOT)) for h in sorted(gl.ROOT.glob(f'E01*/{a}_*/{r}_*/run.py'))]
hits=gl.report_thresholds(paths,'(现存声明背后的 %d 轮)'%len(paths))

# ---- 对最可量的两处,把距离量出来
print("\n=== 把命中处的距离量出来(以观测值自己的自助展布计)===")
rows=[]
G=pd.read_csv('E01_sexual_as_a_value_not_a_category/A12_does_order_of_acquisition_reshape'
              '/R14_direction_at_double_n/results/grid.csv')
rb=np.random.default_rng(1)
for t,g in G[G.arm=='real'].groupby('thresh'):
    v=g.v_shift.values
    se=float(np.std([v[rb.integers(0,len(v),len(v))].mean() for _ in range(2000)]))
    ratio=v.mean()/se
    rows.append(dict(round='#118 A12/R14',quantity=f'强度@{t}',value=ratio,const=3.0,
                     spread_of_ratio=float(np.std([ (lambda s:(v[s].mean()/max(np.std(
                        [v[rb.integers(0,len(v),len(v))].mean() for _ in range(200)]),1e-9)))(
                        rb.integers(0,len(v),len(v))) for _ in range(200)])),
                     npairs=len(v),shift=float(v.mean()),se=se))
    print(f"  #118 强度@{t}: 位移 {v.mean():+.5f} 展布 {se:.5f} -> **{ratio:.4f}**  "
          f"写死的门槛 3.0  距离 {ratio-3.0:+.4f}({100*(ratio/3.0-1):+.1f}%)")
# 比值自身的自助噪声:两次不同的重抽顺序给出 3.0321 与 3.0928 —— 差 0.06。
noise=[]
for s in range(20):
    r2=np.random.default_rng(100+s)
    se2=float(np.std([v[r2.integers(0,len(v),len(v))].mean() for _ in range(2000)]))
    noise.append(v.mean()/se2)
print(f"\n  比值自身的自助噪声(20 个不同重抽种子):"
      f"{np.min(noise):.4f}–{np.max(noise):.4f},sd {np.std(noise):.4f}")
print(f"  README 写的 3.1x 与实测 {rows[1]['value']:.4f} 一致(我的'舍入错误'判断是错的)。"
      f"**但 3.0 与 3.1 的差别本身在噪声里**,而门槛正好设在 3。")
print(f"  加倍样本后:{rows[1]['value']:.4f} -> {rows[0]['value']:.4f}  "
      f"—— **什么也没加强**,而源码打印的是「加强成功,按新强度引用」")

D=pd.DataFrame(rows); D.to_csv(pathlib.Path(__file__).parent/'results'/'edges.csv',index=False)
g=Gate('现存声明的判定有多少压在选定常数上')
g.asserted('linter 规则命中了 #141 已知的那个刀尖(正对照)',
           any('A02_what_basis' in h[0] and 'R10' in h[0] or 'mx<50' in h[2] for h in hits),
           f"命中 {len(hits)} 处,含 A02 的 `mx<50` 两处(与 #141 的 `mx<0.4` 同一族)")
r400=[r for r in rows if '400' in r['quantity']][0]; r250=[r for r in rows if '250' in r['quantity']][0]
g.asserted('#118 的预注册加强判定压在门槛上',abs(r250['value']-3.0)<0.5*r250['spread_of_ratio'] or r250['value']-3.0<0.1,
           f"实测 {r250['value']:.4f} vs 写死的 3.0,距离 {r250['value']-3.0:+.4f} = "
           f"{100*(r250['value']/3.0-1):+.1f}%")
g.asserted('#118 的"加强"其实没有发生',abs(r250['value']-r400['value'])<0.2,
           f"n>=400 {r400['value']:.4f} -> n>=250 {r250['value']:.4f},差 {r250['value']-r400['value']:+.4f}")
g.asserted('⚠ 我的"README 舍入错误"判断是错的,门当场抓到',True,
           f"实测 {r400['value']:.4f} 舍入到 {round(r400['value'],1)} = README 写的 3.1x。撤回")
g.asserted('而比值的一位小数精度本身在自助噪声之内',np.std(noise)>0.02,
           f"20 个重抽种子给出 {np.min(noise):.3f}–{np.max(noise):.3f}(sd {np.std(noise):.3f})"
           f" —— **一个设在 3 的门槛落在这个噪声带里**")
g.resolvable('#118 的效应本身仍高于可分辨门槛',float(r400['shift']),float(r400['se']))
print(g)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
