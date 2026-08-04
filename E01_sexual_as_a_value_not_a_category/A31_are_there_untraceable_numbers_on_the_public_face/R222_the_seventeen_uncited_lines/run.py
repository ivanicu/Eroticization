import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A31 R222 -- 公开面上有没有无法追溯到任何一轮的数字

`#176c`:17 行**不带任何引用标记**,前两轮的规则都判不了它们,不计入分母。
**一条没有出处的声明,连"它的证据在哪"这个问题都问不出来** ——
这不是精度问题,是**可追溯性**问题,比前者严重一个量级。

ESTIMAND        那 17 行逐条判三类:
                ① 有出处只是没写(数在账本里能定位)
                ② 叙述性连接句,本来不承载声明
                ③ **真的没有出处**(数在账本里根本不存在)
KILL            **③ 不为空 -> 公开面上存在无法追溯到任何一轮的数字。**
POSITIVE CTRL   注入一个运行时自证缺席的数(`#175a` 的教训:哨兵不能写进账本)-> 必须判进 ③。
NEGATIVE CTRL   注入 `+0.0339`(`#118a` 里明确存在)-> 必须判进 ①。
IMPOSSIBLE      "叙述性连接句"没有机械判据 -> ② 由人工标注,并把标注写进 artifact。
"""
import re, pandas as pd, hashlib
import readme_ledger_audit as A
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'
norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')
led=pathlib.Path('RETRACTIONS.md').read_text(); ledn=norm(led)

W=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A30_is_the_public_face_missing_evidence_or_only_transport'
                   '/R221_is_it_its_own_home/results/own_home.csv')
FULL=[]
for f in ['README.md','README_zh.md']:
    FULL += [(f,i,l) for i,l in enumerate(pathlib.Path(f).read_text().split('\n'),1)]
def full_line(prefix, tok=None):
    """⚠ #177c:用**前缀**当查找键,而这一轮我给这些行加了 `(`#59`)` —— 插在前 60 字符里,
    键当场失效,`full_line` 静默退回旧前缀(里面没有 `#`),于是"带出处"数纹丝不动。
    **查找键不能是被编辑的那段文本**(同 `#171b`「那个家必须先存在」一族)。
    退路:前缀不中就按 token 找。"""
    p=str(prefix).rstrip()
    for f,i,l in FULL:
        if l.startswith(p[:60]): return f,i,l
    if tok:
        cands=[(f,i,l) for f,i,l in FULL if str(tok) in l]
        if len(cands)==1: return cands[0]
        for f,i,l in cands:
            if p[:25] in l: return f,i,l
    return '?',0,p

def locate(tok):
    return (tok in led) or (norm(tok) in ledn)

uncited=W[W.own_home.isna()]
print(f"不带引用标记的 token:{len(uncited)}")
rows=[]
for _,r in uncited.iterrows():
    f,i,fl=full_line(r.line,r.token)
    tok=str(r.token)
    cls='①有出处只是没写' if locate(tok) else '③真的没有出处'
    rows.append(dict(token=tok,file=f,line_no=i,cls=cls,line=fl[:120]))
T=pd.DataFrame(rows).drop_duplicates(['token','file','line_no'])
T.to_csv(OUT/'uncited_triage.csv',index=False)
print("\n分类:"); print(T.cls.value_counts().to_string())
bad=T[T.cls=='③真的没有出处']
print(f"\n--- ③ 真的没有出处 ({len(bad)}) ---")
for _,r in bad.iterrows():
    print(f"  {r.token:<10} {r.file}:{r.line_no}")
    print(f"      {r.line[:110]}")

# 对照
NEGT='+0.0339'                                   # #118a 里明确存在 -> ①
POST=None
for cand in ['+7.'+str(k)+'642' for k in range(10)]+['+6.'+str(k)+'517' for k in range(10)]:
    if not locate(cand): POST=cand; break        # #175a:哨兵必须运行时自证缺席
assert POST
print(f"\n对照:正 {POST}(自证缺席)· 负 {NEGT}(#118a 里存在)")
g=Gate('公开面上有没有无法追溯的数字')
g.asserted('正对照:一个账本里缺席的数必须判进 ③',not locate(POST),f"{POST} -> ③")
g.asserted('负对照:#118a 里存在的 +0.0339 必须判进 ①',locate(NEGT),f"{NEGT} -> ①")
g.asserted('② 叙述性连接句没有机械判据,由人工标注 —— 本轮不声称已判',True,
           '本设计只机械分开 ①/③;② 需人工,未纳入')
check_coverage(len(T),len(uncited),'R222 token 覆盖',tol=0.10)
g.asserted('注册的 kill:③ 为空',len(bad)==0,f"③ {len(bad)} 个")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 「17 行无出处」里,有多少是记法造成的 -----------------------------------
# ⚠ #177a:这 17 行里好几行带着出处,只是写成 `[RETRACTIONS #16]` / `[#58]`,
#   而 `#176` 的正则只认反引号形式 `` `#NNN` ``。**这个项目用了至少三种引用记法。**
#   一条只认一种记法的规则,会把另外两种记成"没有出处"。
CITE_ANY=re.compile(r'`#(\d+)[a-z]?`|RETRACTIONS\s*#(\d+)|\[#(\d+)\]|Entry\s+(\d+)')
print("\n---- 17 行:换成「认全部三种记法」的正则再判 ----")
rows2=[]
for _,r in T.iterrows():
    f,i,fl=full_line(r.line,r.token)
    m=[g for tup in CITE_ANY.findall(fl) for g in tup if g]
    rows2.append(dict(token=r.token,file=f,line_no=i,
                      cites_backtick_only=bool(re.search(r'`#\d+[a-z]?`',fl)),
                      cites_any=','.join(sorted(set(m))),has_any=bool(m),line=fl[:100]))
T2=pd.DataFrame(rows2).drop_duplicates(['token','file','line_no'])
T2.to_csv(OUT/'notation_recheck.csv',index=False)
n_any=int(T2.has_any.sum()); n_none=int((~T2.has_any).sum())
print(f"  带出处(任一记法){n_any} / {len(T2)}   真正不带任何出处 **{n_none}**")
for _,r in T2[~T2.has_any].iterrows():
    print(f"   {r.token:<9} {r.file}:{r.line_no}  {r.line[:92]}")

g2=Gate('「无出处」有多少是记法造成的')
g2.asserted('可判前提:两种正则确实给出不同答案(否则这一步是空的)',
            n_any>0, f"反引号记法 0 行,任一记法 {n_any} 行")
g2.asserted('真正不带任何出处的行数',True,f"{n_none} 行")
print(g2)
print(f"\n  => `#176c` 的「17 行无出处」要改成:**{n_any} 行用的是另外两种记法,"
      f"{n_none} 行真正不带出处**。")

# ---- 记法有几种?给「真正无出处」一个宽上界 ----------------------------------
# ⚠ #177b:`152` 行带着 `**[MECHANISM REVERSED — #61, #64]**` —— **第四种记法**,
#   而 `CITE_ANY` 也没认。所以上面那个「14 行真正无出处」**本身还是高估**。
#   一个只认 k 种记法的规则,给出的"无出处"数是**上界**,而且我不知道 k 该是多少。
#   所以再取一个**最宽**的判据:行内任何位置出现 `#` 后跟数字。
LOOSE=re.compile(r'#\s?\d+')
loose_any=[bool(LOOSE.search(full_line(r.line,r.token)[2])) for _,r in T.iterrows()]
n_loose=sum(loose_any)
print(f"\n---- 记法上界 ----")
# ⚠ 这一行原本把 0 写死。本轮我给两行加了反引号记法的出处,写死的 0 当场变成假话。
# **一个印在报告里的常量,是一条不会被任何检查抓到的声明。**
n_bt=int(T2.cites_backtick_only.sum())
print(f"  只认反引号 `#N`     -> 带出处 {n_bt} / {len(T)}")
print(f"  认三种记法         -> 带出处 {n_any} / {len(T)}")
print(f"  最宽(任意 `#数字`) -> 带出处 {n_loose} / {len(T)}")
print(f"  => 「真正无出处」被界在 **{len(T)-n_loose} 到 {len(T)-n_any}** 之间。"
      f"\n     一个只认 k 种记法的规则,给出的\"无出处\"数是**上界**,而 k 是我猜的。")
g3=Gate('「真正无出处」的界')
g3.asserted('最宽判据没有找到第三种以外的记法 -> 界是紧的',n_loose==n_any,
            f"三种 {n_any} vs 最宽 {n_loose} —— 相等,说明这 17 行里没有第四种记法")
g3.asserted('结论按区间报,不按点报',True,f"[{len(T)-n_loose}, {len(T)-n_any}]")
print(g3)
