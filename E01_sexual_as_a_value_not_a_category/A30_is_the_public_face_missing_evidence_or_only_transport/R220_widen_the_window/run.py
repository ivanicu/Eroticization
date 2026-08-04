import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A30 R220 -- ② 证据缺口,有多少是我的窗口太窄造出来的

`#174` 的「结构上做不到」:限定语按**同一行**判,账本里写在**下一段**的看不见 -> 会高估 ②。

ESTIMAND        同一批裸数,把匹配窗口从**行**放宽到**条目**(`## Entry N` 到下一个 `## Entry`),
                重跑分类,量 ② 的下降。
KILL            **② 掉到 5 行以下 -> "证据缺口"基本上是窗口造出来的,剩下的全部工作就是搬字。**
POSITIVE CTRL   窗口放宽必须让**至少一个**已知的 ② 变成 ① —— 否则放宽没生效,整轮不可读。
NEGATIVE CTRL   注入 `+9.9991`(账本里不存在):无论窗口多宽都必须留在 ②。
GAUGE           放宽窗口同时会**放宽假阳性** —— 一条条目里任何地方出现的限定语都会算数。
                所以同时报告**窗口内的距离**:限定语与那个数隔了几行。
                距离中位数若很大,说明"同一条条目"这个窗口已经太宽,匹配到的不是同一件事。
IMPOSSIBLE      仍然按**串**匹配不按**所指**(`#170b`)。窗口越宽,这个问题越重。
"""
import re, pandas as pd, hashlib
import readme_ledger_audit as A
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'
norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')

led=pathlib.Path('RETRACTIONS.md').read_text()
# 条目切块:`## Entry N` 到下一个;顶部表格算一块
parts=re.split(r'(?m)^## Entry ',led)
entries=[parts[0]]+['## Entry '+p for p in parts[1:]]
print(f"账本切成 {len(entries)} 块(含顶部表格块)")

naked=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A29_is_it_the_numbers_or_the_qualifiers_that_leak'
                       '/R218_qualifiers_stripped_from_surviving_numbers/results/naked_number_lines.csv')

def classify(tok, window):
    """window='line' | 'entry'。返回 (类, 证据, 距离行数)"""
    if window=='line':
        homes=[l for l in led.split('\n') if tok in norm(l) or tok in l]
        qual=[l for l in homes if A._QUAL.search(l)]
        return ('①搬运缺口' if qual else '②证据缺口'), (qual[0][:150] if qual else ''), 0
    best=None
    for e in entries:
        L=e.split('\n')
        hit=[i for i,l in enumerate(L) if tok in norm(l) or tok in l]
        if not hit: continue
        q=[i for i,l in enumerate(L) if A._QUAL.search(l)]
        if not q: continue
        d=min(abs(i-j) for i in hit for j in q)
        if best is None or d<best[0]: best=(d,L[min(q,key=lambda j:min(abs(i-j) for i in hit))][:150])
    return ('①搬运缺口' if best else '②证据缺口'), (best[1] if best else ''), (best[0] if best else -1)

rows=[]
for _,r in naked.iterrows():
    ln=str(r.line)
    for t in set(norm(m) for m in A._MAGNUM.findall(ln)):
        c1,_,_ = classify(t,'line'); c2,ev,d = classify(t,'entry')
        rows.append(dict(token=t,line=ln[:100],cls_line=c1,cls_entry=c2,dist=d,evidence=ev[:120]))
T=pd.DataFrame(rows).drop_duplicates(['token','line']); T.to_csv(OUT/'widened.csv',index=False)

def byline(col):
    return T.groupby('line')[col].agg(lambda s:'①搬运缺口' if '①搬运缺口' in set(s) else '②证据缺口')
b1,b2=byline('cls_line'),byline('cls_entry')
n2_line=int((b1=='②证据缺口').sum()); n2_entry=int((b2=='②证据缺口').sum())
print(f"\n② 证据缺口:窗口=行 **{n2_line}** -> 窗口=条目 **{n2_entry}**")
print(f"① 搬运缺口:{int((b1=='①搬运缺口').sum())} -> {int((b2=='①搬运缺口').sum())}")
moved=set(b1[b1=='②证据缺口'].index)&set(b2[b2=='①搬运缺口'].index)
print(f"\n放宽后从 ② 变成 ① 的行:{len(moved)}")
for l in list(moved)[:6]: print(f"  {l[:100]}")

d=T[(T.cls_entry=='①搬运缺口')&(T.dist>=0)].dist
print(f"\nGAUGE 窗口内距离:中位 {d.median():.0f} 行 · 90 分位 {d.quantile(0.9):.0f} 行 · 最大 {d.max():.0f} 行")
print(f"  其中同一行(距离 0)的占 {100*(d==0).mean():.0f}%")

# ⚠ #175a:上一轮的哨兵 `+9.9991` 被 `#174` 的条目正文**写进了账本** ——
#   它在被写下来的那一刻就不再"不存在"了,于是负对照判成 ①。
#   **一个发表在账本里的哨兵值,就不再是缺席的。** 哨兵必须在运行时**自证缺席**。
NEG=None
for cand in ['+9.'+str(k)+'321' for k in range(10)]+['+8.'+str(k)+'764' for k in range(10)]:
    if cand not in led and norm(cand) not in norm(led): NEG=cand; break
assert NEG, '找不到一个账本里缺席的哨兵'
print(f"\n负对照哨兵(运行时自证缺席):{NEG}")
neg=classify(NEG,'entry')[0]
g=Gate('② 有多少是窗口造出来的')
g.asserted('正对照:放宽窗口至少让一个 ② 变成 ①',len(moved)>0,f"{len(moved)} 行")
g.asserted('负对照:账本里不存在的数,窗口再宽也留在 ②',neg=='②证据缺口',f"{NEG} -> {neg}")
g.asserted('GAUGE:放宽后匹配到的多数仍在近处(否则窗口太宽,匹配的不是同一件事)',
           float((d<=3).mean())>0.5,f"距离 ≤3 行的占 {100*(d<=3).mean():.0f}%")
check_coverage(len(b2),len(naked),'R220 逐行覆盖',tol=0.10)
g.asserted('注册的 kill:② 掉到 5 行以下',n2_entry<5,f"{n2_line} -> {n2_entry}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 产出:条目窗口下的完整搬运清单 -----------------------------------------
print("\n---- 搬运清单(条目窗口)----")
P=T[T.cls_entry=='①搬运缺口'][['token','line','dist','evidence']].drop_duplicates('token')
P=P.sort_values('dist')
P.to_csv(OUT/'transport_list_entry_window.csv',index=False)
for _,r in P.head(8).iterrows():
    print(f"  {r.token:<9} d={r.dist}  README: {r.line[:60]}")
    print(f"  {'':<9}         账本  : {r.evidence[:100]}")
print(f"\n  共 {len(P)} 个 token,全部有账本内的限定语原句 -> results/transport_list_entry_window.csv")
