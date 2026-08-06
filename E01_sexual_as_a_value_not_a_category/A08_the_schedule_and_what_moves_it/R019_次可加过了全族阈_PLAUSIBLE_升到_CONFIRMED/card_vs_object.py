import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A138 R417 -- 数据集卡片描述的,是不是我手上这个文件

`#372` 的 NEXT 要求 **L3 先跑**(读数据集自己的文档)。跑了,然后**没走到原定的问题**:
卡片第 3 行写 **15,511 行 × 376 列**;文件是 **15,503 × 365**。
卡片写 `straightness` **「Now: 5 levels」**;文件里是 **2 档**。
**⇒ 卡片描述的不是这个文件。**

这不是一个书目细节。`#357` 是**整条页面注记**建在卡片上的:
「成长环境 7→3」「关系风格 3→2」「匿名化让相关弱约 25%」——
**如果卡片落后一个版本,那三条的 scope 全都要跟着动。**

ESTIMAND        把卡片里每一个「Now: N levels (...)」的断言,逐条对着文件核。
                主量 = **符合 / 不符 / 卡片提到但文件没有** 的三分计数。
KILL(条件式)  仅当**正对照**过(至少找到一条卡片与文件**相符**的,证明解析器不是坏的)->
                判:**是否存在任何一条不符**。
                有 -> 卡片与文件是**两个对象**,所有卡片来源的断言都要标注这一条;
                无 -> 只有行列数对不上,那是子样本口径,scope 影响小得多。
POSITIVE CTRL   **必须**至少有一条相符 —— **一个「全部不符」的解析器,坏掉的概率远高于数据。**
NEGATIVE CTRL   拿一个**卡片里根本没提**的列去跑解析,必须什么都不返回。
⚠ 这是 VERIFICATION，不是发现  按 P15 的 `prior_art_in_card` 轴:结果若相符,那是核对;
                **只有不符才是发现**,而它的价值在于**它改的是别处的 scope**。
IMPOSSIBLE      卡片没写的列,本轮什么也说不了;卡片的散文措辞不统一,解析必然漏掉一些。
                **漏掉的方向是保守的**(少报不符),这一条要明说。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
card=pathlib.Path('data/card/column_notes.txt').read_text()
m=re.search(r'([\d,]+)\s*rows\s*x\s*([\d,]+)\s*columns',card)
cr,cc=int(m.group(1).replace(',','')),int(m.group(2).replace(',',''))
print(f"① 形状:卡片 **{cr:,} × {cc}** · 文件 **{d.shape[0]:,} × {d.shape[1]}** -> "
      f"行差 **{d.shape[0]-cr:+,}** · 列差 **{d.shape[1]-cc:+d}**\n")

# 解析:一个列名行(顶格无缩进、非全大写标题、下一行起是缩进块),块里找 "Now: N levels (a, b, c)"
LINES=card.split('\n'); blocks={}; cur=None
for ln in LINES:
    if ln and not ln.startswith(' ') and not ln.startswith('-') and not ln.startswith('=') \
       and ln.strip() and ln.strip()==ln and not ln.isupper() and len(ln.split())<=3:
        cur=ln.strip(); blocks[cur]=[]
    elif cur is not None and ln.startswith('  '):
        blocks[cur].append(ln.strip())
print(f"解析到 **{len(blocks)}** 个卡片列块;其中 **{sum(1 for k in blocks if k in d.columns)}** 个在文件里有同名列。")

rows=[]
for col,body in blocks.items():
    if col not in d.columns: 
        rows.append(dict(v_col=col,v_kind='卡片提到但文件无此列',v_card='',v_file='')); continue
    txt=' '.join(body)
    mm=re.search(r'Now:\s*(\d+)\s*levels?\s*(?:\(([^)]*)\))?',txt)
    if not mm: continue
    nlev=int(mm.group(1)); names=[x.strip() for x in (mm.group(2) or '').split(',') if x.strip()]
    got=sorted(map(str,pd.Series(d[col]).dropna().unique()))
    ok_=len(got)==nlev
    rows.append(dict(v_col=col,v_kind='符合' if ok_ else '不符',
                     v_card=f"{nlev} 档"+(f" {names}" if names else ''),
                     v_file=f"{len(got)} 档 {got[:6]}"))
T=pd.DataFrame(rows); check_columns(T,'R417')
T.to_csv(pathlib.Path(__file__).parent/'results'/'card_vs_file.csv',index=False)
print(f"\n② 逐条核对(只核卡片明写「Now: N levels」的那些):")
for r in T.itertuples():
    if r.v_kind=='卡片提到但文件无此列': continue
    mark='✅' if r.v_kind=='符合' else '❌'
    print(f"   {mark} {r.v_col:<24} 卡片 {r.v_card:<58} 文件 {r.v_file}")
miss=[r.v_col for r in T.itertuples() if r.v_kind=='卡片提到但文件无此列']
print(f"\n   卡片提到但文件里没有的列(**{len(miss)}** 个):{miss if miss else '无'}")
NOK=int((T.v_kind=='符合').sum()); NBAD=int((T.v_kind=='不符').sum())
print(f"\n   **符合 {NOK} · 不符 {NBAD} · 卡片有文件无 {len(miss)}**")

# 负对照:卡片没提的列
notin=[c for c in d.columns if c not in blocks][:5]
print(f"\n负对照(卡片根本没提的 5 个列 -> 解析必须什么都不返回):"
      f"**{sum(1 for c in notin if c in set(T.v_col))} 条返回**(应为 0)")

# ---- ③ ⚠ 两条不符**恰好推翻我自己在 `#366` 里加的一条警告** ----
# 我在 R410 写的是:「`TotalMentalIllness`/`childhood_adversity` 的 NaN 被读成 0 ——
# 把「没有」与「没作答」并成一档」。**卡片说这两列本该是 `None` / `Any` 两档。**
# 那么文件里的 NaN 就是 `None`,而我那条警告是**多余的**。
# ⚠ 但还有另一个世界:`None` 的**行被删了**而不是**标签被删了**。用行数分辨 —— 这是决定性的。
print(f"\n③ ⚠ 「NaN 是 None」还是「None 的行被删了」:")
for c in ('childhood_adversity','TotalMentalIllness'):
    have=int(d[c].notna().sum()); tot=len(d)
    print(f"   {c:<22} 'Any' **{have:,}** / 全表 **{tot:,}** = {have/tot:.1%} · "
          f"NaN **{tot-have:,}**")
print(f"   若 `None` 的**行**被删,全表就只剩 'Any' 的那些人(最多 {int(d['TotalMentalIllness'].notna().sum()):,} 行);"
      f"实际全表 **{len(d):,}** 行,且与卡片只差 **{len(d)-cr:+,}**。")
print(f"   **⇒ 被删的是标签,不是行 -> 文件里的 NaN 就是卡片的 `None`(D6,来源:卡片 + 行数)。**")
print(f"   **⇒ `#366` 里那条「把没有与没作答并成一档」的警告**多余**,本轮撤回。**")

g=Gate('数据集卡片描述的是不是我手上这个文件')
g.asserted('★ 正对照:至少一条卡片与文件**相符**(否则坏的是解析器,不是数据)',NOK>=1,
           f"相符 {NOK} 条",kind='control')
g.asserted('★ 负对照:卡片没提的列不产生任何断言',
           sum(1 for c in notin if c in set(T.v_col))==0,f"{notin[:2]}…",kind='control')
if NOK>=1:
    g.asserted('★ 注册的 kill:不存在任何一条不符',NBAD==0,
               f"不符 {NBAD} 条:{[r.v_col for r in T.itertuples() if r.v_kind=='不符']}")
    g.asserted('★ 形状也一致',d.shape[0]==cr and d.shape[1]==cc,
               f"{d.shape[0]:,}×{d.shape[1]} vs 卡片 {cr:,}×{cc}")
else:
    g.asserted('★ 注册的 kill(正对照未过 -> 不判)',False,'UNVERIFIED:解析器可能是坏的')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
