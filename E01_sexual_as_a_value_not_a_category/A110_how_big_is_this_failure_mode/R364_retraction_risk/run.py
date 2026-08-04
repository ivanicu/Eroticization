import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A110 R364 -- 「无网格」这个特征的撤回风险比

`#318d`:本次会话四次「单格结果写成无格的话」,四次都被下一轮抓到。
**那是我事后数出来的。数一遍全账本,看这个失败模式有多大。**

⚠⚠ **这是自评,而自评是空的**(door ③)。所以本轮**只报计数与比值,不报「所以我该怎样」** ——
后者需要一个干净上下文的对抗者,而本会话不能派。

ESTIMAND        全部账本条目:① 它**当时有没有网格**(规格曲线 / 多种子 / 多切法 / 多口径)
                ② 它**有没有被后续条目撤回或降级**;
                2×2 -> **「无网格」的撤回风险比 RR**。
KILL            **若 RR 明显 > 1 -> 「无网格」是一个可查的撤回风险因子;
                若 ≈ 1 -> 我这四次是巧合,而 `#318d` 那句自述没有支持。**
POSITIVE CTRL   `#294`/`#299`/`#310`/`#318` 四条必须被识别为「撤回/降级了某条」。
NEGATIVE CTRL   随机重排「被撤回」标签 -> RR 必须回到 1(报 200 次重排的分布)。
IMPOSSIBLE      「被撤回」只统计**账本里明写**的;一条错而没人发现的条目在这里记为「未被撤回」——
                **所以 RR 是一个下界**,而且它对「我更愿意去复查哪一类条目」这个选择敏感。
"""
import re,pandas as pd,numpy as np,hashlib
from lib.gates import Gate, check_columns

LED=pathlib.Path('RETRACTIONS.md').read_text()
ENT={};cur=None
for l in LED.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); ENT[cur]=[]
    elif cur is not None: ENT[cur].append(l)
ENT={k:'\n'.join(v) for k,v in ENT.items()}
GRID=re.compile(r'规格曲线|specification|跨种子|多种子|两向|扫描|sweep|网格|grid|旋钮|knob|'
                r'三种|四种|五点|逐变量|每一格|各口径|口径',re.I)
RETR=re.compile(r'撤回|降级|纠正|推翻|收窄|retract|downgrad|overturn',re.I)
REF=re.compile(r'`?#(\d+)[a-z]?`?')
retracted=set()
for k,t in ENT.items():
    for line in t.split('\n'):
        if not RETR.search(line): continue
        for r in REF.findall(line):
            r=int(r)
            if r in ENT and r<k: retracted.add(r)
rows=[dict(v_entry=k,grid=bool(GRID.search(t)),retr=(k in retracted)) for k,t in ENT.items()]
T=pd.DataFrame(rows); check_columns(T,'R364')
T.to_csv(pathlib.Path(__file__).parent/'results'/'risk.csv',index=False)
a=int(((~T.grid)&T.retr).sum()); b=int(((~T.grid)&(~T.retr)).sum())
c=int((T.grid&T.retr).sum());  e=int((T.grid&(~T.retr)).sum())
p_no=a/max(a+b,1); p_yes=c/max(c+e,1); RR=p_no/max(p_yes,1e-9)
print(f"账本条目 **{len(T)}** 条 · 被后续条目明写撤回/降级/纠正的 **{int(T.retr.sum())}** 条 "
      f"({100*T.retr.mean():.1f}%)")
print(f"检出网格的 **{int(T.grid.sum())}** 条({100*T.grid.mean():.0f}%)\n")
print(f"{'':<12}{'被撤回':>10}{'未被撤回':>12}{'撤回率':>10}")
print(f"{'无网格':<12}{a:>10}{b:>12}{100*p_no:>9.1f}%")
print(f"{'有网格':<12}{c:>10}{e:>12}{100*p_yes:>9.1f}%")
print(f"\n★ **风险比 RR = {RR:.2f}**")
rg=np.random.default_rng(2024); lab=T.retr.values.copy(); NPERM=200
nul=[]
for _ in range(NPERM):
    p=rg.permutation(lab)
    a2=int(((~T.grid.values)&p).sum()); b2=int(((~T.grid.values)&(~p)).sum())
    c2=int((T.grid.values&p).sum());   e2=int((T.grid.values&(~p)).sum())
    nul.append((a2/max(a2+b2,1))/max(c2/max(c2+e2,1),1e-9))
nul=np.array(nul); q=float((nul>=RR).mean())
print(f"负对照(随机重排「被撤回」标签 {NPERM} 次):RR **{nul.mean():.2f} ± {nul.std():.2f}** · "
      f"零里 ≥ 观测的比例 **{q:.3f}**")
# ⚠ **事后切片,明确标注:不是检验。** RR≈1 有三个可能原因,其中一个可以便宜地看一眼:
#    这四次是不是**近期**行为,而非全账本的性质。
def rr_of(sub):
    a2=int(((~sub.grid)&sub.retr).sum()); b2=int(((~sub.grid)&(~sub.retr)).sum())
    c2=int((sub.grid&sub.retr).sum());   e2=int((sub.grid&(~sub.retr)).sum())
    if min(a2+b2,c2+e2)<5: return np.nan,(a2,b2,c2,e2)
    return (a2/max(a2+b2,1))/max(c2/max(c2+e2,1),1e-9),(a2,b2,c2,e2)
print(f"\n⚠ **事后切片(标注:不是检验,是描述)**:")
for lo in (0,200,260,280):
    sub=T[T.v_entry>=lo]; r_,cnt=rr_of(sub)
    print(f"   条目 ≥ #{lo:<4} n={len(sub):>3} · RR = **{r_:.2f}**" if np.isfinite(r_)
          else f"   条目 ≥ #{lo:<4} n={len(sub):>3} · RR 不可算(格子太小 {cnt})")
print(f"   ⚠ n 小时 RR 的抽样噪声极大(全账本的置换零就有 ±0.23),**这一行不作为证据**。")

PC=[294,299,310,318]
found=[k for k in PC if any(RETR.search(l) and REF.findall(l) for l in ENT.get(k,'').split('\n'))]
print(f"\n正对照:`#294`/`#299`/`#310`/`#318` 被识别为「撤回了某条」的:**{found}**")
print(f"   它们各自撤回的目标是否被标记:" +
      ' · '.join(f"#{k}->{sorted(int(x) for l in ENT[k].split(chr(10)) if RETR.search(l) for x in REF.findall(l) if int(x)<k)[:2]}"
                 for k in found))
gg=Gate('「无网格」的撤回风险比')
gg.asserted('★ 正对照:四条已知的撤回条目必须被识别',len(found)==4,f"识别到 {found}")
gg.asserted('★ 负对照:随机重排后 RR 回到 1',abs(nul.mean()-1)<0.25,
            f"重排 RR {nul.mean():.2f} ± {nul.std():.2f}")
gg.asserted('★ 注册的 kill:RR 是否明显 > 1(置换分位数 < 0.05)',q<0.05,
            f"RR **{RR:.2f}** · 零里 ≥ 观测 **{q:.3f}**(零 {nul.mean():.2f} ± {nul.std():.2f})")
gg.asserted('⚠⚠ 边界:这是自评,而自评是空的(door ③)',True,
            '只报计数与比值,**不报「所以我该怎样」** —— 后者需要干净上下文的对抗者')
gg.asserted('⚠ 边界:「被撤回」只统计账本明写的 -> RR 是一个下界',True,
            '一条错而没人发现的条目在这里记为「未被撤回」;它也对「我更愿意复查哪一类」敏感')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
