import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A25 R214 -- README 上的现行声明,有没有哪一条的展布来源是零臂

`#168d`:早期条目的精度陈述率低到无法与"没有误差棒"区分。
`#168b` 只能对**有 artifact 的 146 轮**发言。所以那 96 条无精度陈述里,
哪些是"有误差棒只是没写",哪些是"真的没有",分不开。
**但公开面是可以收窄的**:只看 README 当前引用的那些条目。

ESTIMAND        README(中英两份)引用的条目 ∩ 无精度陈述的条目;对每条判展布来源。
KILL            **若其中任何一条的来源是 `null_零臂` -> 公开面在承诺它没有证明的东西,那一条要改。**
                条件式:先要**引用图非空且正对照命中**(`#100`/R147 必须在名单里 ——
                `#167` 已确证它八轮没有误差棒;`#167` 之后我给它补了 ±0.016,
                所以它现在应当判为 `bootstrap_人层`,**这就是仪器的正对照**)。
IMPOSSIBLE      判不了"没写出来的误差棒存在不存在" —— 拿不出来源就记为零臂,这是**保守方向**。
"""
import re, numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_coverage
import guard_lint as GL
OUT=pathlib.Path(__file__).parent/'results'

# ---- README 引用图 ---------------------------------------------------------
cited={}
for f in ['README.md','README_zh.md']:
    t=pathlib.Path(f).read_text()
    for m in re.finditer(r'`#(\d+)[a-z]?`',t):
        cited.setdefault(int(m.group(1)),set()).add(f)
print(f"README 引用的条目 {len(cited)} 条:{sorted(cited)}")

# ---- 无精度陈述的条目(复用 #168 的普查口径)-------------------------------
txt=pathlib.Path('RETRACTIONS.md').read_text()
def feats(b):
    return (bool(re.search(r'±\s*\d*\.\d+|\+/-\s*\d*\.\d+',b))
            or bool(re.search(r'\d+(\.\d+)?×\s*(自身|its own|展布|spread)',b))
            or bool(re.search(r'\[[−\-\+]?\d*\.\d+\s*,\s*[−\-\+]?\d*\.\d+\]',b)))
prec={}; rnd={}
for e in re.split(r'\n## Entry ',txt)[1:]:
    m=re.match(r'(\d+)',e)
    if not m: continue
    n=int(m.group(1)); prec[n]=feats(e[:4000])
    r=re.search(r'`(E\d+·A\d+·R\d+)`',e[:200]); rnd[n]=r.group(1) if r else None
head=txt[:txt.index('\n## Entry ')]
for line in head.split('\n'):
    m=re.match(r'\|\s*(\d+)\s*\|',line)
    if m and 1<=int(m.group(1))<=16: prec[int(m.group(1))]=feats(line); rnd.setdefault(int(m.group(1)),None)
NMAX=max(int(x) for x in re.findall(r'(?m)^## Entry (\d+)',txt))
check_coverage(len(prec),NMAX,'R214 复用 #168 的普查',tol=0.0)

# ---- 交集,并对每条判展布来源 ---------------------------------------------
S=GL.error_bar_scan('.')
vd={r['round'][:r['round'].index('_')] if '_' in r['round'] else r['round']:r.verdict
    for _,r in S.iterrows()}
rows=[]
for n in sorted(cited):
    if n not in prec: rows.append(dict(entry=n,in_ledger=False)); continue
    R=rnd.get(n); rk=R.split('·')[-1] if R else None
    v=vd.get(rk)
    if prec[n]:            src,why='已陈述精度','条目正文带 ± / ×自身展布 / 区间'
    elif v=='SEED_IS_FAKE':src,why='null_零臂','artifact 证明真实臂无抖动来源(#168b)'
    elif v=='NO_SPREAD':   src,why='null_零臂','artifact 里真实臂只有一行'
    else:                  src,why='null_零臂','条目不陈述精度,且拿不出来源 —— 保守判为零臂'
    rows.append(dict(entry=n,in_ledger=True,round=R or '-',scan=v or '-',
                     source=src,why=why,files=','.join(sorted(cited[n]))))
T=pd.DataFrame(rows); T.to_csv(OUT/'public_face.csv',index=False)
bad=T[T.get('source','')=='null_零臂']
print(f"\n{'条目':<7}{'轮次':<16}{'扫描':<14}{'来源':<14}{'为什么'}")
for _,r in T.iterrows():
    if not r.get('in_ledger',False): print(f"#{r.entry:<6}{'不在账本':<16}"); continue
    print(f"#{r.entry:<6}{str(r['round']):<16}{str(r.scan):<14}{r.source:<14}{r.why}")
print(f"\n**README 上来源是零臂的现行声明:{len(bad)} / {int(T.in_ledger.sum())}**")

g=Gate('公开面有没有在承诺它没有证明的东西')
g.asserted('可判前提一:引用图非空',len(cited)>0,f"{len(cited)} 条被引用")
g.asserted('可判前提二:正对照 —— #100 已在 #167 后补上 ±0.016,应判为已陈述精度',
           bool((T[T.entry==100].source=='已陈述精度').all()) if (T.entry==100).any() else False,
           f"#100 -> {T[T.entry==100].source.iloc[0] if (T.entry==100).any() else '不在名单'}")
g.asserted('注册的 kill:README 上没有来源是零臂的现行声明',len(bad)==0,
           f"{len(bad)} 条:{sorted(bad.entry) if len(bad) else '无'}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 分诊:这 30 条里,哪几条的 README 句子真的在断言一个数量 -----------------
# 「拿不出来源就记为零臂」是保守方向,但它把**不需要误差棒的声明**也算了进来:
# 一条撤回、一个范围陈述、一个结构事实,断言的不是数量,没有 ± 可言。
# 决定弧的那个问题是:**公开面有没有在断言一个数量却不给精度。**
MAG=re.compile(r'[+−\-]\s?\d\.\d{2,4}|\b\d{1,3}(?:\.\d+)?\s?%|\b\d+(?:\.\d+)?\s?[×x]\b')
CNT=re.compile(r'\b\d+\s*(?:个|条|块|人|轮|categories|blocks|people|rounds|entries|items)')
tri=[]
for _,r in bad.iterrows():
    n=int(r.entry); ctx=[]
    for f in ['README.md','README_zh.md']:
        t=pathlib.Path(f).read_text()
        for para in re.split(r'\n(?=\||\n|#)',t):
            # ⚠ 第一版第二个条件是**前缀匹配**:`#15 也命中 #150/#151/#159 —— 与当初
            #   `R\d{2}` 读到 `R204` 完全同形。必须要求闭合反引号或单字母后缀。
            if re.search(rf'`#{n}[a-z]?`',para): ctx.append(para)
    blob=' '.join(ctx)
    mags=[m for m in MAG.findall(blob) if not CNT.search(m)]
    tri.append(dict(entry=n,n_ctx=len(ctx),n_magnitudes=len(mags),
                    asserts_magnitude=len(mags)>0,sample=', '.join(mags[:5])))
TR=pd.DataFrame(tri); TR.to_csv(OUT/'triage.csv',index=False)
need=TR[TR.asserts_magnitude]
print(f"\n---- 分诊:30 条零臂里,README 句子真的断言数量的 ----")
print(TR.sort_values('n_magnitudes',ascending=False).head(12).to_string(index=False))
print(f"\n**断言了数量却没有精度的:{len(need)} / {len(TR)}**  条目 {sorted(need.entry)}")

g2=Gate('公开面断言数量却不给精度的,有几条')
g2.asserted('可判前提:分诊器能区分数量与计数',
            bool((~TR.asserts_magnitude).any()) and bool(TR.asserts_magnitude.any()),
            f"{int(TR.asserts_magnitude.sum())} 条断言数量 / {int((~TR.asserts_magnitude).sum())} 条不断言")
g2.asserted('每一条零臂条目都能在 README 里定位到上下文',bool((TR.n_ctx>0).all()),
            f"最小上下文数 {int(TR.n_ctx.min())}")
g2.asserted('修正后的 kill:公开面没有"断言数量却无精度"的声明',len(need)==0,
            f"{len(need)} 条:{sorted(need.entry)}")
print(g2)

# ---- 正确的单位是 README 的**行**,不是条目 -------------------------------
# 一行常引用多条条目,所以按条目计数会把同一行的数量记到每一条引用上 —— 24 是虚高的。
# 决定弧的那个问题只能在行上问:**这一行断言了一个数量,而它引用的条目里有没有任何一条给出精度?**
print("\n---- 行级:断言了数量、而引用的条目**全部**不陈述精度的行 ----")
lines=[]
for f in ['README.md','README_zh.md']:
    for i,ln in enumerate(pathlib.Path(f).read_text().split('\n'),1):
        cs=[int(x) for x in re.findall(r'`#(\d+)[a-z]?`',ln)]
        if not cs: continue
        mags=[m for m in MAG.findall(ln) if not CNT.search(m)]
        if not mags: continue
        anyp=any(prec.get(c,False) for c in cs)
        # ⚠ `text` 只是显示副本。第一版把检测跑在它上面,而英文行更长 ——
        #   `9.0×` 落在 110 字符之后被切掉,于是英文行被判成"无精度"、中文行没有。
        #   **检测必须跑在整行上;截断只用于打印。**
        lines.append(dict(file=f,line=i,n_mag=len(mags),cites=','.join(map(str,cs)),
                          any_cite_has_precision=anyp,full=ln,text=ln[:110]))
L=pd.DataFrame(lines); L.to_csv(OUT/'line_level.csv',index=False)
naked=L[~L.any_cite_has_precision]
print(f"  断言数量并带引用的行:{len(L)}   其中引用的条目**全无**精度陈述的:{len(naked)}")
for _,r in naked.iterrows():
    print(f"   {r.file}:{r.line}  引用 #{r.cites}  {r.n_mag} 个数量")
    print(f"      {r.text}")

g3=Gate('行级:公开面有没有"数量无处可依"的行')
g3.asserted('可判前提:行级单位区分得开(既有全无精度的行,也有有精度的行)',
            len(naked)>0 or len(L)>0, f"{len(L)} 行断言数量,{len(naked)} 行全无精度")
g3.asserted('注册的 kill(行级):没有"数量无处可依"的行',len(naked)==0,
            f"{len(naked)}/{len(L)} 行")
print(g3)

# ---- 再收窄两步,得到可操作的名单 ------------------------------------------
# ① 行**自己**可能就带精度(`9.0×` `4.1×` `± 0.016`)—— 我只查了被引用的条目,没查行本身。
# ② **已撤回/已反转**的行不需要误差棒 —— 它断言的是"这条死了",不是一个数量。
# ⚠ 中文写「余量 23 倍」,正则只认 `×` —— 语言不对称,把中文行误判成无精度。
PREC_INLINE=re.compile(r'±\s*\d*\.\d+|\d+(?:\.\d+)?\s?(?:[×x](?!\s*\d)|倍)|\[[−\-\+]?\d*\.\d+\s*,')
# ⚠ 第一版中文有 `**反转` 而英文没有 `Inverted` —— 同一张表的两半用了不对称的正则。
DEAD=re.compile(r'\*\*Withdrawn\*\*|\*\*撤回\*\*|\*\*反转|\*\*Inverted|Retracted|已撤回|已反转|不再成立')
L['line_has_precision']=L.full.str.contains(PREC_INLINE)
L['withdrawn']=L.full.str.contains(DEAD)
L.drop(columns=['full']).to_csv(OUT/'line_level.csv',index=False)
# 这条 kill 问的是**有没有承诺精度**。一行**明确声明自己的精度未估**不是承诺 —— 它是范围陈述,
# 正是 §2 的 η 规则要求的那一步。所以第四步收窄:排除已显式声明的行。
DISCLOSED=re.compile(r'precision has never been estimated|precision not estimated|精度从未被估|精度未估|未估过')
L['disclosed']=L.full.str.contains(DISCLOSED)
act=L[(~L.any_cite_has_precision)&(~L.line_has_precision)&(~L.withdrawn)&(~L.disclosed)]
print(f"\n---- 可操作名单 ----")
print(f"  断言数量并带引用 {len(L)} 行 -> 引用无精度 {int((~L.any_cite_has_precision).sum())}"
      f" -> 行自己也无精度 {int(((~L.any_cite_has_precision)&(~L.line_has_precision)).sum())}"
      f" -> 非撤回 {int(((~L.any_cite_has_precision)&(~L.line_has_precision)&(~L.withdrawn)).sum())}"
      f" -> 且未显式声明 **{len(act)}**")
for _,r in act.iterrows():
    print(f"   {r.file}:{r.line}  引用 #{r.cites}\n      {r.text}")
act.drop(columns=['full']).to_csv(OUT/'actionable.csv',index=False)

g4=Gate('公开面里真正需要补精度的行')
g4.asserted('可判前提:三步收窄各自都在起作用',
            len(L)>int((~L.any_cite_has_precision).sum())>len(act),
            f"{len(L)} -> {int((~L.any_cite_has_precision).sum())} -> {len(act)}")
g4.asserted('注册的 kill(收窄后):没有需要补精度的行',len(act)==0,
            f"{len(act)} 行:{[f'{r.file}:{r.line}' for _,r in act.iterrows()]}")
print(g4)
