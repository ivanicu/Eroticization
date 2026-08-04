import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A15 R03 -- 有边界吗?还是连通度从一开始就单调下降?

#134d:「晚获得的类别更不相连」绑在**绝对年龄**上(人内中位数分割不可分辨,1.0x),
不是"你自己序列里靠后的"。#134c:量级随切点变 2.4 倍,只能报方向。

但切点扫描被**分箱**卡死了:起始年龄的分箱是 ... 13.5 -> 15.5 -> 17.5 -> 22 -> 28,
所以只有三个真正不同的切点。**绕开切点**:直接对每个 (人, 类别) 求这个类别与这个人
其余类别的平均连通度,按类别的**获得年龄分箱**画一条曲线。

    KNEE   曲线有拐点 -> 边界存在,#132/#133 用"青春期"这个词是对的
    LINE   曲线是直线 -> 没有边界,"版图在 17 岁关上"必须改成
           "连通度从一开始就单调下降",而所有关于 17 岁的措辞都要撤

ESTIMAND        conn(bin) = 该获得年龄分箱内,类别与同一人其余类别的平均去稀有度连通度,
                减去人内置换起始年龄标签后的同分箱期望。
IDENTIFICATION  零 = 人内置换起始年龄标签,精确保留这个人的类别集与年龄分布。
                所以曲线的**形状**是配对结构的,不是边际的。
SCOPE           >=8 个类别起始年龄的人,全部年龄档合并(#134a:五档同号,无梯度)。
WORLDS          KNEE / LINE
KILL            条件式:拐点检验必须先通过**两个**对照才读 ——
                (a) 种植一个真拐点,必须被检出;
                (b) 种植一条**纯直线**,拐点检验必须**不**开火。
                (b) 是关键:一个总会开火的拐点检验是一个不能失败的检查(#96a)。
POSITIVE CTRL   见 (a)。
NEGATIVE CTRL   见 (b),以及人内置换零。
NOISE FLOOR     200 次按人自助;拐点增益的零分布由置换数据给出(200 次)。
MULTIPLICITY    分段点扫过全部内部分箱,整条增益曲线发表。
IMPOSSIBLE      分箱本身是 release 给的(2 年宽,尾部更宽),所以拐点的**位置**至多精确到
                一个分箱边界。本轮只判**有没有拐点**,不判它精确在哪一岁。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])   # 跨轮依赖显式声明(P16)

Ob=obs.astype(float); pj=Ob.mean(0); Cm=(Ob.T@Ob)/len(Ob)
den=np.sqrt(np.outer(pj*(1-pj),pj*(1-pj))); den[den<1e-9]=1e-9
SIM=(Cm-np.outer(pj,pj))/den; np.fill_diagonal(SIM,0.)
iu=np.triu_indices(len(rar),1)
X=np.c_[np.ones(len(iu[0])),rar[iu[0]]+rar[iu[1]],rar[iu[0]]*rar[iu[1]],np.abs(rar[iu[0]]-rar[iu[1]])]
res=SIM[iu]-X@np.linalg.lstsq(X,SIM[iu],rcond=None)[0]
check_residualized(res,rar[iu[0]]+rar[iu[1]],"配对相似度对稀有度")
SIMR=np.zeros_like(SIM); SIMR[iu]=res; SIMR=SIMR+SIMR.T

BINS=sorted(set(np.unique(V[obs]).tolist()))
who=np.flatnonzero(KEEP)
print(f"{len(who):,} 人  获得年龄分箱 {BINS}",flush=True)

def curve(Vm,rng=None,perm=False):
    """每个 (人,类别) 的连通度,按获得年龄分箱聚合。"""
    s={b:0. for b in BINS}; c={b:0 for b in BINS}
    for i in who:
        jj=np.flatnonzero(obs[i]); y=Vm[i,jj].copy()
        if perm: y=y[rng.permutation(len(y))]
        k=len(jj)
        conn=(SIMR[np.ix_(jj,jj)].sum(1))/(k-1)     # 与这个人其余类别的平均连通度
        for b,v in zip(y,conn):
            if b in s: s[b]+=v; c[b]+=1
    return np.array([s[b]/max(c[b],1) for b in BINS]), np.array([c[b] for b in BINS])

obsv,n_b=curve(V)
null=np.mean([curve(V,np.random.default_rng(4400+t),perm=True)[0] for t in range(3)],axis=0)
gap=obsv-null
KEEPB=n_b>=300                        # 太稀的分箱不读
print(f"\n{'年龄分箱':>8} {'n':>8} {'连通度':>10} {'置换零':>10} {'差':>10}")
for b,nn,o,z,gp,kp in zip(BINS,n_b,obsv,null,gap,KEEPB):
    print(f"{b:>8.1f} {nn:>8,} {o:>+10.5f} {z:>+10.5f} {gp:>+10.5f}{'' if kp else '   (n<300,不读)'}")

xb=np.array(BINS)[KEEPB]; yb=gap[KEEPB]; wb=n_b[KEEPB].astype(float)
def fit_lin(x,y,w): 
    Z=np.c_[np.ones(len(x)),x]; b=np.linalg.lstsq(Z*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
    return float(np.sum(w*(y-Z@b)**2))
def fit_knee(x,y,w):
    best=(np.inf,None)
    for k in x[1:-1]:
        Z=np.c_[np.ones(len(x)),x,np.maximum(x-k,0)]
        b=np.linalg.lstsq(Z*np.sqrt(w)[:,None],y*np.sqrt(w),rcond=None)[0]
        r=float(np.sum(w*(y-Z@b)**2))
        if r<best[0]: best=(r,float(k))
    return best
sse_l=fit_lin(xb,yb,wb); sse_k,knee=fit_knee(xb,yb,wb)
gain=1-sse_k/max(sse_l,1e-18)
print(f"\n直线 SSE {sse_l:.3e}   拐点 SSE {sse_k:.3e}   增益 {gain:.3f}   最佳拐点 {knee}")

# 零分布:同样的拐点检验跑在置换曲线上
gn=[]
for t in range(200):
    yz=curve(V,np.random.default_rng(9000+t),perm=True)[0][KEEPB]-null[KEEPB]
    sl=fit_lin(xb,yz,wb); sk,_=fit_knee(xb,yz,wb); gn.append(1-sk/max(sl,1e-18))
gn=np.array(gn); thr=float(np.percentile(gn,95))
print(f"拐点增益的置换零分布:中位 {np.median(gn):.3f}  95 分位 {thr:.3f}")

# 对照 (a) 种植真拐点 / (b) 种植纯直线 —— 检验必须只在 (a) 开火
def synth(shape):
    yy=np.zeros(len(xb))
    if shape=='knee': yy=-0.0006*np.maximum(xb-15.5,0)
    if shape=='line': yy=-0.00012*(xb-xb.mean())
    return yy+np.random.default_rng(7).normal(0,float(np.std(yb))*0.25,len(xb))
res_ctl={}
for sh in ['knee','line']:
    ys=synth(sh); sl=fit_lin(xb,ys,wb); sk,kk=fit_knee(xb,ys,wb); res_ctl[sh]=(1-sk/max(sl,1e-18),kk)
    print(f"  对照 {sh}: 增益 {res_ctl[sh][0]:.3f}  拐点 {kk}")

# ---- 判据:R03 的曲线(晚的连通度**更高**)与 #133a(晚→早的连通度**更低**)方向相反。
#      两者可以同时为真,当且仅当**晚获得的类别彼此抱团**。直接做块分解。
CUT=17.5
def blocks(Vm,rng=None,perm=False):
    ee=[];ll=[];el=[]
    for i in who:
        jj=np.flatnonzero(obs[i]); y=Vm[i,jj].copy()
        if perm: y=y[rng.permutation(len(y))]
        lt=jj[y>CUT]; er=jj[y<=CUT]
        if len(lt)<2 or len(er)<2: continue
        ee.append(SIMR[np.ix_(er,er)][np.triu_indices(len(er),1)].mean())
        ll.append(SIMR[np.ix_(lt,lt)][np.triu_indices(len(lt),1)].mean())
        el.append(SIMR[np.ix_(lt,er)].mean())
    return np.array(ee),np.array(ll),np.array(el)
EE,LL,EL=blocks(V)
rgn=np.random.default_rng(3300); nEE,nLL,nEL=blocks(V,rgn,perm=True)
rbk=np.random.default_rng(101)
def bs(a): return float(np.std([a[rbk.integers(0,len(a),len(a))].mean() for _ in range(200)]))
print(f"\n=== 块分解(切点 {CUT},n={len(EE):,} 人)===")
print(f"  {'块':<14} {'真实':>10} {'人内置换零':>11} {'差':>10} {'展布':>9} {'倍数':>7}")
for nm,a,b_ in [('早 x 早',EE,nEE),('晚 x 晚',LL,nLL),('晚 x 早',EL,nEL)]:
    d=a.mean()-b_.mean(); s=bs(a)
    print(f"  {nm:<14} {a.mean():>+10.5f} {b_.mean():>+11.5f} {d:>+10.5f} {s:>9.5f} {abs(d)/s:>7.1f}x")
d_ll=LL.mean()-nLL.mean(); d_el=EL.mean()-nEL.mean(); d_ee=EE.mean()-nEE.mean()

# ---- 晚到的那一族**是什么**?(把数字变成一句关于内容的话)
import re as _re
lab=[_re.sub(r'\s*\([a-z0-9]+\)$','',c) for c in ons]
lab=[_re.sub(r'^.*?(?:interest in|interested in)\s*','',l)[:34] for l in lab]
latefrac=np.array([np.nanmean((V[obs[:,j],j]>CUT)) for j in range(V.shape[1])])
ordr=np.argsort(-latefrac)
print(f"\n=== 最常晚到的类别 vs 最常早到的(晚 = {CUT} 岁后,即 19+,见 #134e)===")
print(f"  {'最常晚到':<38}{'晚到比例':>8}   |  {'最常早到':<34}{'晚到比例':>8}")
for a,b_ in zip(ordr[:6],ordr[::-1][:6]):
    print(f"  {lab[a][:36]:<38}{latefrac[a]:>8.1%}   |  {lab[b_][:32]:<34}{latefrac[b_]:>8.1%}")
top=ordr[:8]; bot=ordr[-8:]
mt=SIMR[np.ix_(top,top)][np.triu_indices(8,1)].mean()
mb=SIMR[np.ix_(bot,bot)][np.triu_indices(8,1)].mean()
mx=SIMR[np.ix_(top,bot)].mean()
print(f"\n  这 8 个最晚到的彼此去稀有度连通度 {mt:+.5f};8 个最早到的彼此 {mb:+.5f};两组之间 {mx:+.5f}")

g=Gate('有边界吗,还是一条直线')
g.asserted('晚到的那一族在**题目层**也彼此更连通 —— 所以它有内容,不只是人层的配对结构',
           mt>mb and mt>mx,
           f"最晚到的 8 个彼此 {mt:+.5f} > 最早到的 8 个彼此 {mb:+.5f},两组之间 {mx:+.5f}")
g.asserted('矛盾解开:晚获得的类别**彼此抱团**,但与早期那批**断开**',
           d_ll>0 and d_el<0,
           f"晚x晚 {d_ll:+.5f}(高于置换零) 而 晚x早 {d_el:+.5f}(低于置换零);早x早 {d_ee:+.5f}")
g.asserted('(a) 种植的真拐点被检出',res_ctl['knee'][0]>thr,
           f"增益 {res_ctl['knee'][0]:.3f} > 零的 95 分位 {thr:.3f}")
g.asserted('(b) 种植的纯直线**不**触发拐点检验 —— 否则这是个不能失败的检查(#96a)',
           res_ctl['line'][0]<=thr,
           f"增益 {res_ctl['line'][0]:.3f} vs 零的 95 分位 {thr:.3f}")
g.offset_control('真实曲线的拐点增益是否高于置换零',gain,float(np.median(gn)),float(gn.std()),
                 null_kind='人内置换起始年龄标签后,同一个拐点检验的增益分布')
g.no_sign_crossing('可读分箱上的差全部同号',list(yb))
print(g)
D=pd.DataFrame(dict(bin=BINS,n=n_b,conn=obsv,null=null,gap=gap,readable=KEEPB))
D.to_csv(pathlib.Path(__file__).parent/'results'/'curve.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
