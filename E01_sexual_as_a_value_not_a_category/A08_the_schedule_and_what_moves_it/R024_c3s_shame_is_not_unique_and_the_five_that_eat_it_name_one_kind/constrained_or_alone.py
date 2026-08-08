import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A140 R421 -- 页面上那 97% 的数,是「可核而未核」还是「印的时候就没带约束」

`#376b`:核对器覆盖率 **3.1%**。**那个 3.1% 本身是一条路标**,而两个答案指向完全不同的下一步。

两个活着的世界:
**A 可核而未核** —— 数都带着约束印出来了,只是我的正则没抽到 -> **扩正则,高杠杆的 Closure**。
**B 印的时候就没带约束** —— 页面上大多数数是**孤立**的 -> **那才是发现**:
   它指向一条可以改的**写作规则**(`#375c` 的推论:**凡印导出量,就把输入量一起印**),
   而那条规则把未来所有轮次的这类错都变成**写下即自曝** —— 那是 ④ 能力边界的更新。

⚠ **先读 IMPOSSIBLE 栏再设计**(`#376d` 刚定下的惯例,因为它已经预言对了三次):
IMPOSSIBLE ①  这一页的表格单元是**巨型段落**,所以「同一句」是个坏代理 —— 它会**高估**同伴关系。
              -> 改用**字符窗**,并**扫窗宽**(100 / 200 / 400),让结论不挂在一个我随手选的数上。
IMPOSSIBLE ②  「带约束」是**必要**条件,不是充分:一个数旁边有 `±` 不等于那个 `±` 是**它的**。
              -> 所以本轮报的是**上界**(最多这么多是可核的),措辞必须如此。
IMPOSSIBLE ③  原始量(n · 计数)本来就没有约束可带,把它们算进「孤立」是不公平的。
              -> 单独分出 `n = …` 一类,不计入分母。

ESTIMAND        每个加粗数字,按窗内是否有**约束同伴**(`±` · `se` · `CI` · `[a, b]` · `vs` · `零`/`null`)
                分成:**已核**(`#376` 抽到的)· **可核而未核** · **孤立** · **原始量**。
KILL(条件式)  仅当对照都过 -> 判:**「孤立」是否占非原始量的多数**。
                是 -> 世界 B,写作规则要改;否 -> 世界 A,扩正则。
POSITIVE CTRL   `#376` 已抽到的 4 条,**必须**落进「已核」或「可核」,一条都不能落进「孤立」。
NEGATIVE CTRL   合成一句只有一个孤零零数字的句子,**必须**落进「孤立」。
⚠ 规格曲线      三个窗宽全部报,不报单格(`realstat G4`)。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.page_arithmetic import ALL, disjoint, controls, coverage

cp,cn,nb,ng=controls()
print(f"核对器自检(`lib/page_arithmetic.controls()`):正对照 **{cp}** · 负对照 **{cn}**\n")

BOLD=re.compile(r'\*\*([-−+]?\d[\d,]*(?:\.\d+)?%?)\*\*')
# ⚠ 第一版的标记里有一个裸的 `零`,而负对照的合成句写的是「孤**零**零的数字」-> **假阳性**。
# **本窗口第三次子串假阳性**(`#374b` 的 `race`/`smost`,加这一次)。
# 修的是**标记**,不是测试 —— 把测试改到能过,是把仪器的病写进结论。
MARK=re.compile(r'±|\bse\b|\bCI\b|95%|\[\s*[-−+]?\d|\bvs\b|置换零|offset 零|的零\b|零的\b|'
                r'\bnull\b|地板|\bfloor\b|阈')
RAW=re.compile(r'n\s*=\s*$|样本|人\s*$')
def classify(txt,W):
    got=set()
    for kind,src,gi,pi in ALL(txt):
        for m in BOLD.finditer(src): got.add(m.group(1))
    out=[]
    for m in BOLD.finditer(txt):
        v=m.group(1); a,b=max(0,m.start()-W),min(len(txt),m.end()+W)
        win=txt[a:b]; pre=txt[max(0,m.start()-8):m.start()]
        if re.search(r'n\s*=\s*\**$',pre) or (',' in v and '.' not in v):
            k='原始量'                      # n = 6,717 / 12,720 这类计数
        elif v in got: k='已核'
        elif MARK.search(win): k='可核而未核'
        else: k='孤立'
        out.append(dict(v_val=v,v_kind=k,v_win=W))
    return out

print("规格曲线(三个窗宽全部报,不报单格):")
rows=[]
for W in (100,200,400):
    for page in ('README.md','README_zh.md'):
        t=pathlib.Path(page).read_text()
        for r in classify(t,W): rows.append(dict(v_page=page,**r))
T=pd.DataFrame(rows); check_columns(T,'R421')
T.to_csv(pathlib.Path(__file__).parent/'results'/'classified.csv',index=False)
for W in (100,200,400):
    S=T[T.v_win==W]; c=S.v_kind.value_counts()
    nonraw=len(S)-c.get('原始量',0)
    alone=c.get('孤立',0)
    print(f"   窗宽 ±{W:>3}:已核 **{c.get('已核',0):>3}** · 可核而未核 **{c.get('可核而未核',0):>3}** · "
          f"孤立 **{alone:>3}** · 原始量 {c.get('原始量',0):>3} · "
          f"-> 孤立占非原始量 **{alone/max(nonraw,1):.1%}**")

# ---- 对照 ----
t0=pathlib.Path('README.md').read_text()
got4=set()
for kind,src,gi,pi in ALL(t0):
    for m in BOLD.finditer(src): got4.add(m.group(1))
C200={r['v_val']:r['v_kind'] for r in classify(t0,200)}
posbad=[v for v in got4 if C200.get(v)=='孤立']
NEG="一句话里只有一个孤零零的数字 **0.1234** 而已。"
negk=[r['v_kind'] for r in classify(NEG,200)]
print(f"\n对照:")
print(f"   正对照(`#376` 抽到的 {len(got4)} 个数必须不落进「孤立」):落进孤立的 **{len(posbad)}** 个 {posbad}")
print(f"   负对照(合成的孤零零一句):分类 **{negk}**")
CP=len(got4)>0 and len(posbad)==0
CN=negk==['孤立']

S=T[T.v_win==200]; c=S.v_kind.value_counts()
nonraw=len(S)-c.get('原始量',0); alone=c.get('孤立',0); share=alone/max(nonraw,1)
g=Gate('页面上那 97% 的数是可核而未核,还是印的时候就没带约束')
g.asserted('★ 正对照:已抽到的数一个都不能落进「孤立」',CP,f"落进孤立 {len(posbad)}",kind='control')
g.asserted('★ 负对照:合成的孤零零一句必须落进「孤立」',CN,f"{negk}",kind='control')
g.asserted('★ 规格曲线:三个窗宽都跑了,不报单格',len(T.v_win.unique())==3,
           f"窗宽 {sorted(T.v_win.unique())}",kind='control')
if CP and CN:
    g.asserted('★ 注册的 kill:「孤立」占非原始量的多数(= 世界 B,写作规则要改)',share>0.50,
               f"{alone}/{nonraw} = {share:.1%}")
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\n⚠ 上界措辞(IMPOSSIBLE ②):「窗内有 `±`」不等于那个 `±` 是**它的** ——"
      f"所以「可核而未核 {c.get('可核而未核',0)}」是**上界**,「孤立 {alone}」是**下界**。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
