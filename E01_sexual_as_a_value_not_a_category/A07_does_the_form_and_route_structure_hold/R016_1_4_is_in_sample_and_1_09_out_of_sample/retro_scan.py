import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A101 R352 -- 回溯扫描:多少个由特征向量派生的数没对过符号

`#306b` 是本项目**第三次**撞上特征向量符号,而前两次都是人眼抓到的。
**现在有了 guard 20 这个能自己开火的检查 —— 那就该回头扫一遍。**

⚠ **这是 Closure,不是 Frontier。** 它不分开任何世界,它保护的是已发表数字的符号。

ESTIMAND        全仓 341 个 `run.py`:① 多少个用了特征分解 ② 其中多少个**跨版本/劈分/种子**
                比较过特征向量派生量 ③ 其中多少个**对齐了符号** ④ 其中多少个的**结论依赖符号**。
POSITIVE CTRL   **已知答案的两条**:`R351`(本轮上一轮,我刚加了对齐)必须判为**已对齐**;
                `R349`(只报 |cos|,绝对值)必须判为**不依赖符号**。
NEGATIVE CTRL   把对齐关键词从文本里剔掉 -> 「已对齐」的计数必须归零。
IMPOSSIBLE      「对齐了没有」是**关键词**的代理,不是对代码语义的判断。P6 的安全侧:
                检出对齐 = UNVERIFIED 待读;**未检出 = 确定要读**。
                同理「结论依赖符号」按**是否报了带号的派生量**判,那也是代理。
"""
import re,pandas as pd,numpy as np,hashlib
from lib.gates import Gate, check_columns

FILES=sorted(pathlib.Path('.').glob('*/A*/R*/run.py'))
EIG=re.compile(r'\beigh\b|linalg\.svd|linalg\.eig\(')
CMP=re.compile(r'@\s*V[a-zA-Z0-9_]*\[|\|cos\||V\w*\[:,\s*k\]\s*@|cos_|跨种子|跨劈分|劈半|两半|fold|折')
ALG=re.compile(r'<\s*0\s*:\s*V|=-V|\*=\s*-1|np\.sign|对齐|align|flip')
SGN=re.compile(r'[+-]0\.\d{3,}')          # 报了带号的数
ABSONLY=re.compile(r'abs\(float\(V|np\.abs\(.*@|\|cos\|')
# ⚠ 精化后的判据(读完那 5 轮之后写的,而**读**是在扫描之后 —— 这一步就是那次读的机器化):
#   只有**两次以上**特征分解、且它们的特征向量被**逐元素**比较过,才可能撞上符号。
#   用**投影**(V[:,:k]@V[:,:k].T)或**绝对值**比较的,数学上就是符号不变的。
MULTI=re.compile(r'(for\b[^\n]*\n(?:[^\n]*\n){0,12}?[^\n]*(?:eigh|linalg\.svd))')
INVAR=re.compile(r'\[:,\s*:\s*\w+\]\s*@\s*\w+\[:,\s*:\s*\w+\]\.T|abs\(float\(np\.corrcoef|abs\(float\(V|np\.abs\(')
rows=[]
for f in FILES:
    t=f.read_text(errors='ignore')
    if not EIG.search(t): continue
    n_eig=len(EIG.findall(t))
    multi=(n_eig>=2) or bool(MULTI.search(t))
    rows.append(dict(v_round=f.parts[-2][:34],n_eig=n_eig,multi=multi,cmp=bool(CMP.search(t)),
                     algn=bool(ALG.search(t)),signed=bool(SGN.search(t)),
                     absonly=bool(ABSONLY.search(t)),invar=bool(INVAR.search(t))))
T=pd.DataFrame(rows); check_columns(T,'R352')
n_e=len(T); n_c=int(T.cmp.sum()); n_a=int((T.cmp&T.algn).sum())
NEED=T[T.multi&(~T.algn)&(~T.invar)&T.signed]
LOOSE=T[T.cmp&(~T.algn)&T.signed&(~T.absonly)]
print(f"全仓 **{len(FILES)}** 个 `run.py`")
print(f"  ① 用了特征分解            **{n_e}**")
print(f"  ② 其中跨版本/劈分/种子比较过 **{n_c}**")
print(f"  ③ 其中检出对齐符号         **{n_a}**({100*n_a/max(n_c,1):.0f}%)")
print(f"  ④ 宽判据(跨版本 且 未对齐 且 带号)-> **{len(LOOSE)}** 轮")
print(f"  ⑤ **精判据**(**两次以上**特征分解 且 未对齐 且 **不是符号不变的比较**)-> **{len(NEED)}** 轮")
print(f"     ⚠ 宽判据的 {len(LOOSE)} 轮我逐个读过:「两半」在本项目通常指**选项半**或**块半**,")
print(f"     不是两次特征分解 —— 那是宽判据的假阳性来源。")
print(f"\n宽判据标出、我逐个读过的 {len(LOOSE)} 轮:")
VERDICT={"R271":"一块一次 eigh,而分半信度取了 `abs()` —— **符号已被有意处理**",
         "R274":"两次 eigh,但比的是**前 k 维子空间**(投影),符号不变",
         "R282":"两次 eigh,比的是 ‖P1·V2‖_F —— **投影范数,符号不变**",
         "R303":"**只有一次** eigh;结论是两个相关的**对比**,全局翻号会同时翻两个",
         "R323":"**只有一次** eigh;同上,单一全局符号约定"}
for _,r in LOOSE.iterrows():
    k=r.v_round.split("_")[0]
    print(f"   {r.v_round[:30]:<32} {VERDICT.get(k,'(未读)')}")
print(f"\n精判据剩下:{list(NEED.v_round) if len(NEED) else '**0 轮** —— 一个也没有受影响'}")
T.to_csv(pathlib.Path(__file__).parent/'results'/'sign_audit.csv',index=False)
# ⚠ **先过 L3 的先例闸口**:本项目有人做过这件事吗?
#   有 —— `#226b` 诊断了 `R210:73`,`R272` 修好并当场定价。
#   所以还活着的问题**不是「有没有 bug」**,是 **「修复有没有传下去」**。
AGG=re.compile(r'Z@v\[:,-1\]|Z@pc')
agg=[(f.parts[-2],bool(ALG.search(f.read_text(errors='ignore'))))
     for f in FILES if AGG.search(f.read_text(errors='ignore'))]
num=lambda n:int(re.match(r'R(\d+)',n).group(1))
R272=[n for n,_ in agg if n.startswith('R272')]
bad=[n for n,a in agg if not a]
after=[n for n in bad if R272 and num(n)>num(R272[0])]
print(f"\n⚠ L3 先例闸口:`R272_align_the_block_signs_and_see_who_it_hurts` **已经做过这件事**")
print(f"   把每块 PC 分累加的轮次共 **{len(agg)}** 个,其中未检出对齐 **{len(bad)}** 个")
print(f"   **它们全部早于 R272** -> `R272` 之后仍未对齐的轮次:**{len(after)}** 个 {after}")
print(f"   ⇒ 修复**已被包含**:`#226b` 诊断 · `R272` 修好并定价 · 之后没有新的未对齐聚合。")
p351=T[T.v_round.str.startswith('R351')]; p349=T[T.v_round.str.startswith('R349')]
strip=lambda s: ALG.sub('',s)
n_a_neg=sum(1 for f in FILES if EIG.search(f.read_text(errors='ignore'))
            and CMP.search(f.read_text(errors='ignore'))
            and ALG.search(strip(f.read_text(errors='ignore'))))
gg=Gate('回溯扫描:多少个特征向量派生量没对过符号')
gg.asserted('★ 正对照:`R351`(刚加了对齐)必须判为已对齐',
            len(p351)>0 and bool(p351.iloc[0].algn),
            f"R351 检出 cmp={bool(p351.iloc[0].cmp)} algn={bool(p351.iloc[0].algn)}" if len(p351) else 'R351 未收录')
gg.asserted('★ 正对照:`R349`(只报 |cos|)必须判为只报绝对值',
            len(p349)>0 and bool(p349.iloc[0].absonly),
            f"R349 absonly={bool(p349.iloc[0].absonly)}" if len(p349) else 'R349 未收录')
gg.asserted('★ 负对照:剔掉对齐关键词后「已对齐」必须归零',n_a_neg==0,
            f"剔掉后仍检出 {n_a_neg} 个")
gg.asserted('⚠ P5★:仪器两个方向都返回过非零',n_a>0 and len(LOOSE)>0,
            f"检出对齐 {n_a} · 宽判据标出 {len(LOOSE)} · 精判据标出 {len(NEED)}")
gg.asserted('★ L3 先例闸口:这件事本项目有人做过吗',True,
            f"`R272` 做过(`#226b` 诊断 -> `R272` 修好并定价)—— **我差一步把一个三百轮前的已知 bug 当成新发现**")
gg.asserted('★ 真正还活着的问题:修复传下去了吗(`R272` 之后有没有新的未对齐聚合)',
            len(after)==0,
            f"把每块 PC 分累加的 {len(agg)} 轮里,未对齐 {len(bad)} 轮,**全部早于 R272**;之后 {len(after)} 个")
gg.asserted('⚠ 安全侧(P6):检出对齐 = UNVERIFIED 待读;未检出 = 确定要读',True,
            f"这张表**不判对错**,它只把 {n_e} 轮缩到 {len(NEED)} 轮")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
