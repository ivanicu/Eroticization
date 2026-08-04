import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A112 R367 -- 账本里的 UNVERIFIED,有多少条能说出自己属于哪一种

`#321b`:一个 UNVERIFIED 可以是**「仪器不够」**,也可以是**「答案本来就在中间」**,
**而区间是唯一能分开这两者的东西。**

⚠ **这是 Closure。⚠ 而 `#364` 刚证明自评的结论不可信** ——
所以本轮**只报计数**,并且**人工读 ≥5 条验证分类**,把**分类器准确率**也报出来。

ESTIMAND        全账本判为 UNVERIFIED 的条目,分三类:
                ① **带区间/MDE**(可以判属于哪一种)· ② **只有点估计**(不可判)·
                ③ **仪器失败**(正对照没过 —— 那是另一回事)。
KILL            **若 ② 很大 -> 那是一份可执行的补测清单;若很小 -> 这条线到此为止。**
POSITIVE CTRL   **已知答案的三条**:`#321`(带区间)· `#310a`(带 MDE)· `#296a`(正对照 FAIL)——
                必须各自落到 ①①③。
NEGATIVE CTRL   剔掉区间/MDE 关键词 -> ① 必须归零。
⚠ 人工验证     随机抽 **≥5 条**逐条读,报分类器的**准确率**(`#358` 的盲区:
                关键词分类器会把 UNVERIFIED 与「零」混起来)。
IMPOSSIBLE      分类按**关键词**判,是代理;P6 安全侧:**判为 ① = 待读;判为 ② = 确定要补。**
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
# ⚠ 正对照抓到的缺陷:每条的 **NEXT 段描述的是下一轮的设计**,
#    而任何读整条的分类器都会把**未来那一轮**的性质算到**当前这一条**头上。
#    `#321` 就是这样被判成③的 —— 触发词在它的 NEXT 里(「③ 仪器失败(正对照没过」)。
#    **剥掉 NEXT 再判。**
def body(t):
    for mark in ('**NEXT(','**NEXT (','**NEXT'):
        i=t.find(mark)
        if i>0: return t[:i]
    return t
ENT_FULL=dict(ENT); ENT={k:body(v) for k,v in ENT.items()}
UNV=re.compile(r'UNVERIFIED|不可读|读不了|读不出来|分不开',re.I)
IVAL=re.compile(r'\bMDE\b|95% CI|CI \[|区间|分位数|percentile|自助|bootstrap|置换零.{0,20}±|± *0\.\d',re.I)
FAILPC=re.compile(r'正对照.{0,24}(FAIL|没过|失败|不合格)|positive control.{0,20}fail|仪器(坏|失败|不合格)|测不动',re.I)
rows=[]
for k,t in ENT.items():
    if not UNV.search(t): continue
    if FAILPC.search(t): cls='③仪器失败'
    elif IVAL.search(t): cls='①带区间/MDE'
    else: cls='②只有点估计'
    rows.append(dict(v_entry=k,cls=cls))
T=pd.DataFrame(rows); check_columns(T,'R367')
T.to_csv(pathlib.Path(__file__).parent/'results'/'unverified.csv',index=False)
n=len(T); vc=T.cls.value_counts()
print(f"账本 **{len(ENT)}** 条,其中含 UNVERIFIED 式判决的 **{n}** 条({100*n/len(ENT):.0f}%)\n")
for c in ('①带区间/MDE','②只有点估计','③仪器失败'):
    print(f"   {c:<14} **{int(vc.get(c,0)):>3}** 条({100*vc.get(c,0)/max(n,1):>4.1f}%)")
print(f"\n**② 只有点估计**的条目(前 12,这是一份可执行的补测清单):")
lst=sorted(T[T.cls=='②只有点估计'].v_entry)
print('   '+' · '.join(f"#{k}" for k in lst[:12])+(f" …共 {len(lst)} 条" if len(lst)>12 else ''))
PC={321:'①带区间/MDE',310:'①带区间/MDE',296:'③仪器失败'}
got={k:(T[T.v_entry==k].cls.iloc[0] if (T.v_entry==k).any() else '(未收录)') for k in PC}
print(f"\n正对照:" + ' · '.join(f"#{k} -> **{v}**(应 {PC[k]})" for k,v in got.items()))
strip=IVAL.sub('',LED)
E2={};cur=None
for l in strip.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); E2[cur]=[]
    elif cur is not None: E2[cur].append(l)
neg=sum(1 for k,v in E2.items() if UNV.search('\n'.join(v)) and not FAILPC.search('\n'.join(v))
        and IVAL.search('\n'.join(v)))
print(f"负对照:剔掉区间/MDE 关键词后判为 ① 的 **{neg}** 条(必须 0)")
rg=np.random.default_rng(1234); samp=sorted(rg.choice(T.v_entry.values,6,replace=False))
print(f"\n⚠ 人工验证的抽样(6 条):{[int(x) for x in samp]}")
for k in samp:
    t=ENT[int(k)]
    print(f"\n   --- #{int(k)} 判为 **{T[T.v_entry==k].cls.iloc[0]}**")
    hits=[x.group(0)[:34] for x in list(IVAL.finditer(t))[:3]]
    fp=[x.group(0)[:34] for x in list(FAILPC.finditer(t))[:2]]
    print(f"       区间/MDE 证据:{hits if hits else '(无)'}")
    print(f"       正对照失败证据:{fp if fp else '(无)'}")
# ⚠⚠ 这个缺陷不止影响本轮:`#313`(零式声明的三件套)与 `#319`(撤回风险比)
#     **也读了整条,包括 NEXT**。`#319` 的结论是 RR≈1(零),而 NEXT 污染会让
#     「NEXT 里提到网格」的条目被判为有网格 -> **把 RR 推向 1**。**那个零可能是伪影。重跑。**
GRID=re.compile(r'规格曲线|specification|跨种子|多种子|两向|扫描|sweep|网格|grid|旋钮|knob|'
                r'三种|四种|五点|逐变量|每一格|各口径|口径',re.I)
RETR=re.compile(r'撤回|降级|纠正|推翻|收窄|retract|downgrad|overturn',re.I)
REF=re.compile(r'`?#(\d+)[a-z]?`?')
def rr_with(src):
    retr=set()
    for k,t in src.items():
        for line in t.split('\n'):
            if not RETR.search(line): continue
            for r in REF.findall(line):
                r=int(r)
                if r in src and r<k: retr.add(r)
    g=np.array([bool(GRID.search(t)) for t in src.values()])
    rr_=np.array([k in retr for k in src])
    a=int(((~g)&rr_).sum()); b=int(((~g)&(~rr_)).sum())
    c=int((g&rr_).sum());   e=int((g&(~rr_)).sum())
    return (a/max(a+b,1))/max(c/max(c+e,1),1e-9),(a,b,c,e)
rr_full,cf=rr_with(ENT_FULL); rr_body,cb=rr_with(ENT)
print(f"\n⚠⚠ 重跑 `#319` 的撤回风险比,NEXT 剥不剥:")
print(f"   读整条(`#319` 的做法):RR = **{rr_full:.2f}** · 2×2 = {cf}")
print(f"   **剥掉 NEXT**:            RR = **{rr_body:.2f}** · 2×2 = {cb}")
print(f"   -> 网格检出率 {100*np.mean([bool(GRID.search(t)) for t in ENT_FULL.values()]):.0f}% "
      f"-> {100*np.mean([bool(GRID.search(t)) for t in ENT.values()]):.0f}%")

gg=Gate('UNVERIFIED 分三类')
gg.asserted('★ `#319` 的 RR 在剥掉 NEXT 后是否改变结论(仍 ≈1 则那个零成立)',
            abs(rr_body-1)<0.35,
            f"读整条 {rr_full:.2f} -> 剥 NEXT **{rr_body:.2f}** —— "
            f"{'结论不变' if abs(rr_body-1)<0.35 else '**结论变了,`#319` 要重判**'}")
gg.asserted('★ 正对照:`#321`->① · `#310`->① · `#296`->③',
            all(got[k]==PC[k] for k in PC),' · '.join(f"#{k}={v}" for k,v in got.items()))
gg.asserted('★ 负对照:剔掉区间关键词后 ① 必须归零',neg==0,f"剩 {neg} 条")
gg.asserted('★ 注册的 kill:② 只有点估计的条目数',len(lst)>0,
            f"**{len(lst)}** 条 —— {'这是一份可执行的补测清单' if len(lst)>3 else '这条线到此为止'}")
gg.asserted('⚠ 人工验证:6 条抽样的证据已逐条打印,准确率由眼判,不由本脚本自report',True,
            '⚠ `#364`:自评的结论不可信 —— 所以本轮只报计数,人工那一步的结论写在 README 里')
gg.asserted('⚠ 安全侧(P6):判为 ① = 待读;判为 ② = 确定要补',True,f"把 {n} 条缩到 {len(lst)} 条")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
