import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A138 R418 -- 「已删除」是一个可核的断言

`#373b` 只核了卡片**明写档数**的 5 个列 —— 卡片里最容易核的一小块。
卡片还有一整节 `DROPPED GT QUESTIONS (not in this dataset)`,**逐条点名**说某些列已被删除。
**删除是可核的**:那些列名若仍在文件里,卡片就不只是**过时**,而是**反的**。
而 `#373` 已经证明卡片与文件是两个对象 —— **所以这一节没有理由被默认相信。**

ESTIMAND        卡片点名的每一个「已删除」标识,在文件的列名里查(精确 · 大小写不敏感 · 子串)。
                主量 = **仍然存在的个数**。
KILL(条件式)  仅当正/负对照都过 -> 判:**仍然存在的个数是否为 0**。
                非 0 -> 卡片的删除清单是**假的**,而那是比过时严重得多的一件事;
                为 0 -> 删除清单**成立**,这一节可以继续被引用。
POSITIVE CTRL   **必须**先证明这段解析能在文件里**找到**东西 ——
                拿一批**确定还在**的列名走同一条匹配,必须全部命中(`#373b` 的教训照搬)。
NEGATIVE CTRL   拿一批**编造的**标识走同一条匹配,必须一个都不命中。
⚠ 匹配的方向    「没找到」是一个**否定**断言,而 `P5★` 说:一个从未返回过非零的仪器给出的 0 不是测量。
                所以正对照不是可选的,**它是这一轮成立的前提**。
IMPOSSIBLE      卡片用的是 GT 原始问题的名字,发布版可能**改过名** ——
                所以「没找到」也可能是「改名了」。本轮只能证否(找到 = 卡片错),不能证成。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
COLS=list(d.columns); LOW=[c.lower() for c in COLS]
# ⚠ **第一版是裸子串匹配,而它给了两个假阳性** —— 正是 §2 点名的 `personal`/`persona` 陷阱:
#   `race`  命中 `Totalsexyrace` · `You find raceplay to be:` —— 那是**情色内容**,不是人口学的族裔题;
#   `smost` 命中 `"I find {bodypartsmost}:"` —— 那是**引用** bodypartsmost 的派生题,不是那个列本身。
# **代理账(P6):属性 = 「该列还在文件里」· 代理 = 「列名里出现这个子串」·
#  蕴含方向 = 「不出现 ⇒ 不在」(健全)、「出现 ⇒ 在」(**不健全**)· 见证 = 上面两条。**
# 所以裸子串只能用来**证否**,不能用来**证在**。分成两个匹配器,把不健全的一侧留成 UNVERIFIED。
import re as _re
def find_sub(tok):
    t=tok.lower().strip()
    return [COLS[i] for i,c in enumerate(LOW) if t and t in c]
def find_token(tok):
    """把列名切成标识符片段,要求**整片**相等 —— 这才是「这个列名是它」的近似。"""
    t=tok.lower().strip()
    out=[]
    for i,c in enumerate(LOW):
        parts=set(_re.split(r'[^a-z0-9]+',c))
        if t==c or t in parts: out.append(COLS[i])
    return out
find=find_sub

# 卡片 `DROPPED GT QUESTIONS` 一节里**点名的标识**(逐字抄自 data/card/column_notes.txt)
DROPPED=[
 ('人口学','country'),('人口学','ethnicity'),('人口学','race'),('人口学','religion_importance'),
 ('人口学','sexage'),
 ('源列(被计算列取代)','attractedtomasculine'),('源列(被计算列取代)','attractedtofeminine'),
 ('源列(被计算列取代)','abuse'),('源列(被计算列取代)','sexualAssaultLevel'),
 ('源列(被计算列取代)','tkpfq6j'),('源列(被计算列取代)','gjt7igv'),
 ('内容','pedophil'),('内容','smost'),
]
rows=[]
for grp,tok in DROPPED:
    hs_=find_sub(tok); ht=find_token(tok)
    rows.append(dict(v_group=grp,v_token=tok,v_nsub=len(hs_),v_ntok=len(ht),
                     v_hits='; '.join(x[:52] for x in hs_[:4])))
T=pd.DataFrame(rows); check_columns(T,'R418')
T.to_csv(pathlib.Path(__file__).parent/'results'/'dropped.csv',index=False)

# ---- 正对照:先证明这条匹配能找到东西 ----
KNOWN=['biomale','age','politics','straightness','pornhabit','animated','written','powerlessnessvariable']
posn=[]
posn=[k for k in KNOWN if find(k)]
print(f"正对照(**必须先证明匹配器能命中**,`P5★`):{len(posn)}/{len(KNOWN)} 个确定还在的列名命中")
FAKE=['zzq_not_a_column','xyzzy_plugh','__nope__','qqqqq_absent']
negn=[k for k in FAKE if find(k)]
print(f"负对照(编造的标识):{len(negn)}/{len(FAKE)} 命中(应为 0)\n")

print("卡片点名「已删除」的标识,逐条在文件列名里查(**整片相等** = 健全的『仍在』判据):")
for r in T.itertuples():
    if r.v_ntok: mark='❌ **仍在**'
    elif r.v_nsub: mark='⚠ 仅子串命中(代理不健全 -> UNVERIFIED,不算仍在)'
    else: mark='✅ 确已不在'
    print(f"   {mark}  [{r.v_group}] `{r.v_token}`" + (f" -> {r.v_hits}" if r.v_nsub else ""))
STILL=int((T.v_ntok>0).sum()); SUBONLY=int(((T.v_nsub>0)&(T.v_ntok==0)).sum())
print(f"\n   **仍然存在(整片相等):{STILL} / {len(T)}** · "
      f"**仅子串命中(我自己的匹配器造的噪声):{SUBONLY}**")

g=Gate('卡片的「已删除」清单是不是真的')
CP=len(posn)==len(KNOWN); CN=len(negn)==0
g.asserted('★ 正对照:确定还在的列名必须全部命中(否则这一轮的 0 是沉默不是无罪)',CP,
           f"{len(posn)}/{len(KNOWN)}",kind='control')
g.asserted('★ 正对照(严格匹配器):同样 8 个必须整片命中',
           all(find_token(k) for k in KNOWN),
           f"{sum(1 for k in KNOWN if find_token(k))}/{len(KNOWN)}",kind='control')
g.asserted('★ 负对照:编造的标识一个都不能命中',CN,f"{len(negn)}/{len(FAKE)}",kind='control')
if CP and CN:
    g.asserted('★ 注册的 kill:卡片点名删除的,一个都不在文件里(整片相等)',STILL==0,
               f"仍在 {STILL} 个:{[r.v_token for r in T.itertuples() if r.v_ntok]} · "
               f"仅子串 {SUBONLY} 个(不算,代理不健全)")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
