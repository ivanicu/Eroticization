import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A30 R221 -- 那个"家",是不是它自己的家

`#175b`:公开面上每一个裸数,账本里都有一个带限定语的家 -> 100% 搬运缺口。
`#175c`:但清单按**串**匹配不按**所指** —— `10%` 与 `−0.028` 配到的都是无关条目。
**所以 `#175b` 数的可能是"有没有一个家",不是"有没有它自己的家"。**

ESTIMAND        70 个 token 逐个判:证据句**所在条目号** ∈ README 那一行**引用的条目号集合**。
KILL            **通过率 < 70% -> 清单不能作为搬运依据,`#175b` 的 100% 要相应降级。**
POSITIVE CTRL   `0.3×` 的证据来自 `#120`,而 README 行正是引 `#120` —— 必须判通过。
                (这一条我已手工搬过,所以它的正确答案是已知的。)
NEGATIVE CTRL   `10%` 的证据来自 `#12`(一张无关表格的 ±10%),而 README 行引的不是 `#12` ——
                必须判**不通过**。这两条把仪器的两个方向都钉住。
IMPOSSIBLE      README 有些行**不带任何引用标记**。那些行判不了 -> 单列一类,不计入分母。
"""
import re, pandas as pd, hashlib
import readme_ledger_audit as A
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'
norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')

led=pathlib.Path('RETRACTIONS.md').read_text()
parts=re.split(r'(?m)^## Entry ',led)
entries=[(0,parts[0])]+[(int(re.match(r'(\d+)',p).group(1)),'## Entry '+p)
                        for p in parts[1:] if re.match(r'(\d+)',p)]

def home_entry(tok):
    """返回 (条目号, 距离, 证据句);条目号 0 = 顶部表格块。"""
    best=None
    for num,e in entries:
        L=e.split('\n')
        hit=[i for i,l in enumerate(L) if tok in norm(l) or tok in l]
        if not hit: continue
        q=[i for i,l in enumerate(L) if A._QUAL.search(l)]
        if not q: continue
        j=min(q,key=lambda j:min(abs(i-j) for i in hit))
        d=min(abs(i-j) for i in hit)
        if best is None or d<best[1]: best=(num,d,L[j][:150])
    return best

W=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A30_is_the_public_face_missing_evidence_or_only_transport'
                   '/R220_widen_the_window/results/widened.csv')
W=W[W.cls_entry=='①搬运缺口'].drop_duplicates('token')
# ⚠ #176a:`widened.csv` 的 `line` 是**截断到 100 字符的显示副本**,而引用标记 `#NNN`
#   写在行尾 —— 于是 70 个 token 里只有 5 个"可判",两个对照都拿不到答案。
#   **第三次**在截断副本上做检测(`#173c` 第二次,`#169d`② 第一次)。
#   修法:用截断前缀回到**当前 README 的整行**。
FULL=[]
for f in ['README.md','README_zh.md']:
    FULL += pathlib.Path(f).read_text().split('\n')
def full_line(prefix):
    p=str(prefix).rstrip()
    for l in FULL:
        if l.startswith(p[:60]): return l
    return p
rows=[]
for _,r in W.iterrows():
    fl=full_line(r.line)
    cited={int(x) for x in re.findall(r'`#(\d+)[a-z]?`',fl)}
    h=home_entry(str(r.token))
    if h is None: continue
    num,d,ev=h
    rows.append(dict(token=r.token,home_entry=num,cited=','.join(map(str,sorted(cited))) or '',
                     n_cited=len(cited),own_home=(num in cited) if cited else None,
                     dist=d,line=fl[:90],full_len=len(fl),evidence=ev[:110]))
T=pd.DataFrame(rows); T.to_csv(OUT/'own_home.csv',index=False)
judgeable=T[T.own_home.notna()]
rate=judgeable.own_home.mean() if len(judgeable) else float('nan')
print(f"token 总数 {len(T)}   可判(README 行带引用){len(judgeable)}   "
      f"不可判(行不带引用){int(T.own_home.isna().sum())}")
print(f"\n**通过率 {judgeable.own_home.sum():.0f}/{len(judgeable)} = {100*rate:.1f}%**  (注册阈值 70%)")
print("\n--- 判为「不是它自己的家」的 ---")
for _,r in judgeable[~judgeable.own_home.astype(bool)].head(12).iterrows():
    print(f"  {str(r.token):<9} 证据在 #{r.home_entry}, 而行引的是 #{r.cited}")
    print(f"  {'':<9} {r.line[:80]}")

pos=T[T.token.astype(str).str.startswith('0.3')]
neg=T[T.token.astype(str)=='10%']
g=Gate('那个家是不是它自己的家')
g.asserted('正对照:0.3× 的证据来自 #120,而 README 行正引 #120 -> 必须通过',
           bool(len(pos) and (pos.own_home==True).any()),
           f"{list(zip(pos.token,pos.home_entry,pos.cited,pos.own_home)) if len(pos) else '不在表里'}")
g.asserted('负对照:10% 的证据来自无关条目 -> 必须不通过',
           bool(len(neg) and (neg.own_home==False).all()),
           f"{list(zip(neg.token,neg.home_entry,neg.cited,neg.own_home)) if len(neg) else '不在表里'}")
check_coverage(len(T),len(W),'R221 token 覆盖',tol=0.10)
g.asserted('注册的 kill:通过率 >= 70%',rate>=0.70,f"{100*rate:.1f}%")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
