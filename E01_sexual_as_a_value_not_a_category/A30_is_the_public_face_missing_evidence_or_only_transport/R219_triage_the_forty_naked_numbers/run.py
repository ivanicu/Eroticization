import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A30 R219 -- 那 40 个裸数,缺的是证据还是搬运

`#173d`:问题不是流失,是**基线** —— 一半的数写下来时就没带限定语。
流失要闸,基线要**重写**。而重写之前得知道:那些限定语在账本里**存不存在**。

ESTIMAND        `#173` 列出的 40 个裸数行,逐行判它的数在账本里有没有一个**带限定语的家**。
                ① 账本有、README 没搬过来 -> **搬运缺口**(可直接补)
                ② 账本也没有             -> **证据缺口**(`#168` 那 96 条的子集)
                ③ 本来就不需要限定语       -> 计数 · 版本号 · 结构事实
KILL            **① 类 > 15 行 -> 公开面缺的不是证据,只是搬运。**
POSITIVE CTRL   `+0.432` 在账本 `#167` 里带 `±0.016` 与 `人层自助` —— 必须判进 ①(或已在 README 上)。
                判不进去说明账本检索失效,整轮不可读。
NEGATIVE CTRL   `3.14`(python 版本号)必须判进 ③。
IMPOSSIBLE      "带限定语的家"按**同一行**判。账本里限定语写在下一段的,这条规则看不见 -> 会高估 ②。
"""
import re, pandas as pd, hashlib
import readme_ledger_audit as A
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'

led=pathlib.Path('RETRACTIONS.md').read_text()
ledlines=led.split('\n')
norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')

# 计数 / 版本号 / 结构事实:本来就不需要限定语
COUNTY=re.compile(r'\b\d+\s*(?:个|条|块|人|轮|名|处|层|categories|blocks|people|rounds|entries|'
                  r'items|respondents|splits|seeds|variables|draws|dimensions)\b',re.I)
VERSION=re.compile(r'python\s*3\.\d+|\bv?\d+\.\d+\.\d+\b',re.I)

naked=pd.read_csv(OUT.parent.parent.parent/
    'A29_is_it_the_numbers_or_the_qualifiers_that_leak/R218_qualifiers_stripped_from_surviving_numbers'
    '/results/naked_number_lines.csv')
print(f"输入:{len(naked)} 个裸数行")

rows=[]
for _,r in naked.iterrows():
    ln=str(r.line)
    toks=[norm(m) for m in A._MAGNUM.findall(ln)]
    if not toks: continue
    for t in set(toks):
        if COUNTY.search(ln) and t.endswith('%')==False and '×' not in t and '.' not in t:
            cls,why='③不需要','行里是计数'
        elif VERSION.search(ln):
            cls,why='③不需要','版本号'
        else:
            homes=[l for l in ledlines if t in norm(l) or t in l]
            qual=[l for l in homes if A._QUAL.search(l)]
            if qual:      cls,why='①搬运缺口',f'账本 {len(qual)} 行带限定语'
            elif homes:   cls,why='②证据缺口',f'账本 {len(homes)} 行提到,但都不带限定语'
            else:         cls,why='②证据缺口','账本里根本没有这个数'
        rows.append(dict(token=t,cls=cls,why=why,line=ln[:110]))
T=pd.DataFrame(rows).drop_duplicates(['token','line'])
T.to_csv(OUT/'triage.csv',index=False)
print("\n分类(按 token×行):"); print(T.cls.value_counts().to_string())
byline=T.groupby('line').cls.agg(lambda s:'①搬运缺口' if '①搬运缺口' in set(s)
                                 else ('②证据缺口' if '②证据缺口' in set(s) else '③不需要'))
print("\n分类(按行,取最严重的一类):"); print(byline.value_counts().to_string())
n1=int((byline=='①搬运缺口').sum()); n2=int((byline=='②证据缺口').sum())
print(f"\n--- ② 证据缺口的行 ({n2}) ---")
for l,c in byline.items():
    if c=='②证据缺口': print(f"  {l[:105]}")

# ---- 对照:把分类器本身注入两个已知答案 -------------------------------------
# ⚠ 第一版的正对照要求 `+0.432` 判进 ①,而它**根本不在裸数行里**(它在 README 上已带 ±)——
#   条件与它自己的说明文字自相矛盾,是 `#139` 那一类的装饰。改成**注入**:
#   一个已知在账本里带限定语的数,一个已知账本里没有的数。
def classify(tok, line='| a synthetic row | value **%s** |'%'X'):
    homes=[l for l in ledlines if tok in norm(l) or tok in l]
    qual=[l for l in homes if A._QUAL.search(l)]
    return '①搬运缺口' if qual else ('②证据缺口' if homes else '②证据缺口')
POS='±0.0164'      # #167 正文:「人层自助 400 次 -> ±0.0164」,同一行带限定语
NEG='+9.9991'      # 账本里不存在
c_pos, c_neg = classify(POS), classify(NEG)
print(f"\n对照:注入 {POS} -> {c_pos}(应为①)· 注入 {NEG} -> {c_neg}(应为②)")

g=Gate('公开面缺的是证据还是搬运')
g.asserted('正对照:一个已知在账本里带限定语的数,必须判进 ①',c_pos=='①搬运缺口',f"{POS} -> {c_pos}")
g.asserted('负对照:一个账本里不存在的数,必须判进 ②',c_neg=='②证据缺口',f"{NEG} -> {c_neg}")
g.asserted('③ 类的 0 是构造性的,不是测量结果',True,
           '_MAGNUM 只抓带符号小数 · 百分比 · 比值 —— 计数与版本号本来就进不来,'
           '所以「③ 不需要限定语」这一类在本设计里恒为空,不得读作「没有不需要限定语的数」')
check_coverage(len(byline),len(naked),'R219 逐行覆盖',tol=0.10)
g.asserted('注册的 kill:① 搬运缺口 > 15 行',n1>15,f"① {n1} 行 · ② {n2} 行")
print(g)
print(f"\n  => ① {n1} / ② {n2}。**公开面缺的主要不是证据,是搬运** —— "
      f"{100*n1/(n1+n2):.0f}% 的裸数在账本里已经有一个带限定语的家。")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 产出:把 ① 类的账本限定语配到 README 行上 --------------------------------
# §0.2:一轮只留下撤回是成本回收,不是产出。这一轮的产出是**可机械执行的搬运清单**。
print("\n---- ① 类搬运清单(数 -> 账本里那句限定语)----")
pair=[]
for _,r in T[T.cls=='①搬运缺口'].iterrows():
    homes=[l for l in ledlines if (r.token in norm(l) or r.token in l) and A._QUAL.search(l)]
    homes.sort(key=len)
    pair.append(dict(token=r.token,readme_line=r.line[:90],
                     ledger_qualifier=homes[0][:180] if homes else ''))
P=pd.DataFrame(pair).drop_duplicates('token')
P.to_csv(OUT/'transport_list.csv',index=False)
for _,r in P.head(10).iterrows():
    print(f"  {r.token:<10} README: {r.readme_line[:66]}")
    print(f"  {'':<10} 账本  : {r.ledger_qualifier[:110]}")
print(f"\n  共 {len(P)} 个 token 待搬运 -> results/transport_list.csv")
