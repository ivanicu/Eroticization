#!/usr/bin/env python3
"""
tools/readme_ledger_audit.py -- README 上的数字,有没有被后来的账本条目撤回过?

由 #143 触发。#141 花了一整轮"重新定价"一条声明,而 **Entry 24 在 118 条之前就已经把它
杀掉了** —— 我用**重跑那一轮**的方式去审计它,却没读那条声明后来的账本条目。
一个轮次自己的输出**不是**那条声明的当前状态;账本才是。

⚠ P6 代理账:
  PROPERTY   README 上的这个数,是不是已经被后来的某条账本条目撤回或改写
  PROXY      该数字的字面串,出现在账本里一行**含撤回类词汇**的文本中
  IMPLICATION 只有一个方向可靠:**命中 -> 这个数确实出现在一条撤回语境里,必须人工读**。
             反过来"未命中 -> 这个数是当前的"**不可靠** —— 撤回可能换了写法或没写数字。
  WITNESS    大量命中是良性的:那个数字是撤回行里的**更正值**,不是被撤的值。必须人工分诊。
  SAFE SIDE  输出是**待读清单**,不是判决。
"""
import re,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
KILL=('撤回','withdrawn','WITHDRAWN','retracted','RETRACTED','killed','RETRACTION',
      'no longer follows','降级','INVALID','reduced')

def audit(readme='README.md',ledger='RETRACTIONS.md'):
    R=(ROOT/readme).read_text(); L=(ROOT/ledger).read_text().splitlines()
    nums=set(re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',R))
    nums|=set(re.findall(r'(?<![\w.])\d\.\d{2}\s*(?:of|/)\s*3',R))
    out=[]
    for n in sorted(nums):
        for i,l in enumerate(L,1):
            if n in l and any(k in l for k in KILL):
                out.append((n,i,l.strip())); break
    return len(nums),out

if __name__=='__main__':
    tot,hits=audit()
    print(f"README 里可辨识的数字 {tot} 个;{len(hits)} 个出现在账本的撤回语境里(**待读清单,非判决**)\n")
    for n,i,l in hits: print(f"  {n:<10} RETRACTIONS.md:{i}\n     {l[:170]}\n")
    print("⚠ SAFE SIDE(#P6):大量命中是良性的 —— 那个数字往往是撤回行里的**更正值**。")
    print("   未命中**不等于**这个数是当前的。输出必须人工分诊。")


# ---------------------------------------------------------------------------
# #144:README 内部一致性 —— 同一个事实被写了两遍,而修正只落在一遍上
#
# #143 找到的是"账本的撤回没走到前页"。更坏的一种是:**走到了前页的一半**。
# 本 README 有两套并行叙述(中文表 + 英文正文),同一条声明在两处各写一次:
#   色情诱导 82.7%    行 56(已修)  vs  行 131(仍写已被认定为"账本里根本没有"的 +0.2515)
#   覆盖度定律 0.815   行 87        vs  行 122
#   三条轴 2.95       行 54        vs  行 111
# 这是 §P16「一个事实一个家」的直接违反,而它的代价与 #143 相同:
# 读者拿到哪一份取决于他先读到哪一段。
#
# ⚠ P6 代理账:
#   PROPERTY   同一条声明的两处叙述,数字是否一致
#   PROXY      同一个**引用标记**(#NN / ANN / EntryNN)出现在两侧,而两侧的数字集合不同
#   IMPLICATION 只有一个方向可靠:**数字集合不同 -> 两处确实写了不同的数,必须人工读**。
#              反过来"相同 -> 一致"**不可靠** —— 两处可能都是旧的。
#   SAFE SIDE  输出是待读清单,不是判决。
import re as _re

def internal_consistency(readme='README.md'):
    lines=(ROOT/readme).read_text().splitlines()
    def nums(s): return set(_re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',s))
    def has_cjk(s): return bool(_re.search(r'[\u4e00-\u9fff]',s))
    cites={}
    for i,l in enumerate(lines,1):
        for c in set(_re.findall(r'#\d{1,3}\b|\bA\d{2}\b|Entry \d{1,3}',l)):
            cites.setdefault(c,[]).append((i,has_cjk(l),nums(l)))
    out=[]
    for c,rows in sorted(cites.items()):
        sides={True:set(),False:set()}
        for i,cj,ns in rows: sides[cj]|=ns
        if rows and len(set(cj for _,cj,_ in rows))==2 and sides[True]!=sides[False]:
            out.append((c,[(i,cj,sorted(ns)) for i,cj,ns in rows]))
    return out

def uncited_numbers(readme='README.md'):
    """#144b:带数字但**不带任何引用标记**的行 —— 它们对任何基于引用的审计都是不可见的。

    #143 的那个真错误(行 131 的 rho +0.2515 与 "85% surviving response-style control")
    正是这样一行:账本明说「账本里根本没有这个数」,而它活了下来,**因为没有任何东西
    把它连回账本**。一个不带出处的数字,是一个无法被撤回的数字。"""
    # ⚠ #145:第一版**按行**判,而引用往往在同一**段**的别处 —— 它把盲区高估了 14 倍
    #   (17 行 vs 实际 2 段,其中一段是 python 版本号)。检查的单位必须与"出处"的单位一致。
    txt=(ROOT/readme).read_text(); paras=_re.split(r'\n\s*\n',txt)
    off=1; out=[]
    for para in paras:
        n=para.count('\n')+2
        ns=_re.findall(r'(?<![\w.])\d+\.\d{2,4}(?![\w])',para)
        if ns and not _re.search(r'#\d{1,3}\b|RETRACTIONS|Entry \d|\bA\d{2}\b|\bR\d{2}\b',para):
            out.append((off,ns,para.strip().replace('\n',' ')[:110]))
        off+=n
    return out

if __name__=='__main__' and '--internal' in __import__('sys').argv:
    hits=internal_consistency()
    print(f"\n=== README 内部一致性:{len(hits)} 个引用标记在两套叙述里带着不同的数字 ===")
    for c,rows in hits:
        print(f"  {c}")
        for i,cj,ns in rows:
            print(f"     行{i:>4} {'中文表' if cj else '英文正文'}  {ns}")
    print("\n⚠ SAFE SIDE(#P6):数字集合相同**不等于**一致 —— 两处可能都是旧的。待读清单。")
    un=uncited_numbers()
    print(f"\n=== 带数字但**不带出处**的行:{len(un)} 行(对任何基于引用的审计都不可见)===")
    for i,ns,s in un: print(f"  行{i:>4}  {ns}\n        {s}")
    print("\n⚠ 一个不带出处的数字,是一个**无法被撤回**的数字。#143 的那个真错误就在这类行里。")


# ============================================================================
# #170:账本点名过的 README 缺陷,现在修好了没有
# ----------------------------------------------------------------------------
# #169c:边跑了、输出了、没有人读。`#144d` 早就记过「两套并行叙述」,而 `+0.815`
# 的第二处断言仍活到 149 条条目之后。**一条被记进账本的 README 缺陷,和一条被修好的
# README 缺陷,在账本里长得一模一样。**
#
# P6 代理账:
#   PROPERTY    账本点名的那处 README 缺陷已经修好
#   PROXY       被点名的字符串在当前 README 里**仍能原样命中**,且附近没有撤回标记
#   IMPLICATION 仍原样命中且无标记 => 未修(**命中方向可读**)
#   WITNESS     ① 一处缺陷可以在改写后仍保留同一个数字(`+0.815` 现在活在删除线里,
#                  带 `KILLED — #20`)—— 所以**必须查邻近的撤回标记**,否则修好的也被判未修。
#               ② **同一个字符串可以有两个所指**:`+0.023` 在行 50 是开放性相关,
#                  在行 111 是信度阶梯里的 sham 值。代理按**字符串**匹配,不按**所指**匹配 ——
#                  所以「同文件内重复」的判定必须人工确认所指相同(#170b)。
#   SAFE SIDE   只在**命中且无标记**方向判「未修」。未命中 != 已修(措辞可能变了)。
# ============================================================================
# ⚠ `#565`:页面上实际在用、而这张词表里**没有**的标记词 —— `DOWNGRADED` 与 `降级`。
#   `README.md:411` 写着 `[DOWNGRADED — RETRACTIONS #24, re-affirmed #143]`,
#   那是一句**带条目指针的范围声明**,而检查因为不认识这个词把它报成「未标记」。
#   ⚠ 补词表会让这道检查**变松**,所以补的同时必须报**前后差值**(见 `#565b`),
#   否则这就是「为了让闸门说是而放松它」——`#564d` 刚写过的那件事的另一面。
_MARK = _re.compile(r'~~|KILLED|CORRECTED|Withdrawn|Inverted|DOWNGRADED|撤回|已杀|反转|未修|降级|scope stated')
_NAMES = _re.compile(r'README[^\n]{0,400}')

def named_defects(readmes=('README.md','README_zh.md'), ledger='RETRACTIONS.md'):
    import pathlib
    L=pathlib.Path(ledger).read_text()
    cur={f:pathlib.Path(f).read_text() for f in readmes}
    out=[]
    for e in _re.split(r'\n## Entry ',L)[1:]:
        m=_re.match(r'(\d+)',e)
        if not m: continue
        n=int(m.group(1))
        for seg in _NAMES.findall(e):
            # 被点名的具体串:反引号里的东西,或一个带符号的数量
            toks=set(_re.findall(r'`([^`\n]{3,40})`',seg))|set(_re.findall(r'[+−-]\d\.\d{3,4}|\b\d{1,3}\.\d%',seg))
            for tk in toks:
                if tk.startswith('#') or tk.startswith('tools/') or len(tk)<4: continue
                # ⚠ #179b:`tk in ln` 是**无边界的子串匹配** —— `+0.093`(性别锚)
                #   会匹配进 `+0.0937`(疗愈),于是一条新增声明凭空造出两个"未修"。
                #   同一个串两个所指(`#170b`)的第二种形态:**前缀包含**。
                #   数字型 token 必须要求右侧不是数字。
                import re as _r2
                pat=(_r2.escape(tk)+r'(?!\d)') if _r2.search(r'\d$',tk) else _r2.escape(tk)
                for f,txt in cur.items():
                    for ln_no,ln in enumerate(txt.split('\n'),1):
                        if _r2.search(pat,ln):
                            out.append(dict(entry=n,token=tk,file=f,line=ln_no,
                                            marked=bool(_MARK.search(ln)),excerpt=ln[:100]))
    import pandas as pd
    D=pd.DataFrame(out)
    return D.drop_duplicates(['entry','token','file','line']) if len(D) else D

def report_named_defects(**kw):
    D=named_defects(**kw)
    if D.empty: print('没有可判的点名'); return D
    live=D[~D.marked]
    # #170a:同一个数在**同一份文件里出现两次** = 并行叙述(缺陷);
    #        中英各一次 = 镜像(正常)。这是 `#144d` 那条缺陷的精确形状。
    dup=live.groupby(['entry','token','file']).size().reset_index(name='n')
    dup=dup[dup.n>1]
    print(f"\n账本点名 × 当前 README 命中 {len(D)} 处,无撤回标记 {len(live)} 处,"
          f"**同一文件内重复(= 并行叙述){len(dup)} 处**")
    for _,r in dup.iterrows():
        where=live[(live.entry==r.entry)&(live.token==r.token)&(live.file==r.file)]
        print(f"  #{r.entry} 点名 `{r.token}` -> {r.file} 的 " +
              ', '.join(str(x) for x in where.line) + " 行各写了一遍")
    if dup.empty: print("  (无:每个被点名的数在每份文件里只有一个家)")
    print("\n⚠ SAFE SIDE(#P6):只在**命中且无标记**方向判「未修」。"
          "\n   未命中 ≠ 已修(措辞可能变了);带标记 = 已改写,不是未修。")
    return D


# ============================================================================
# #171:去重删掉的是重复,还是内容
# ----------------------------------------------------------------------------
# #170c 排出三层(守卫没调用 / 输出没人读 / 记了没人修),第四层没查:**修得对不对。**
# `#170a` 把 R06 正文段的数字删掉只留指针 —— 而没有任何检查能证明删掉的那些数
# 与声明表里留下的那些数是**同一批**。删错一个,读者会看到一个表里没有、正文也没有的洞。
#
# P6 代理账:
#   PROPERTY    一次去重只删掉了重复,没有删掉内容
#   PROXY       修前 README 里出现过的**数量**,修后在整份文件里**一次都不出现**
#   IMPLICATION 修前有、修后全无 => 那个数离开了这一页(**命中方向可读**)
#   WITNESS     一个数可以被**改写**而不是删除(0.2285 → "the correction" 这种改写),
#               所以命中项必须逐条人工判「它是被删了,还是被换了说法」
#   SAFE SIDE   只在**离开**方向判。没有离开 != 去重是对的(可能删的是内容、留的是重复)。
# ============================================================================
_MAGNUM = _re.compile(r'[+−-]\s?\d\.\d{2,4}|\b\d{1,3}(?:\.\d+)?\s?%|\b\d+(?:\.\d+)?\s?[×倍]')

def numbers_that_left(rev='HEAD~1', files=('README.md','README_zh.md')):
    import subprocess, pathlib, pandas as pd
    out=[]
    for f in files:
        try: old=subprocess.run(['git','show',f'{rev}:{f}'],capture_output=True,text=True,check=True).stdout
        except subprocess.CalledProcessError: continue
        new=pathlib.Path(f).read_text()
        # ⚠ #172a:英文写 `10.4×`、中文写 `10.4 倍` —— 字符串不同,于是**每一次翻译
        #   都会被判成删除**(`63e03b7` 一次假报 7 个)。差集之前必须归一化。
        norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')
        o=set(norm(m) for m in _MAGNUM.findall(old)); n=set(norm(m) for m in _MAGNUM.findall(new))
        for tk in sorted(o-n):
            ctx=[l for l in old.split('\n') if tk in l]
            out.append(dict(file=f,token=tk,n_old_lines=len(ctx),
                            old_excerpt=(ctx[0][:120] if ctx else '')))
    return pd.DataFrame(out)

def report_numbers_that_left(rev='HEAD~1', **kw):
    D=numbers_that_left(rev=rev, **kw)
    if D.empty:
        print(f"\n对比 {rev} → 工作区:没有数量离开这一页"); return D
    print(f"\n对比 {rev} → 工作区:**{len(D)} 个数量离开了这一页**")
    for _,r in D.iterrows():
        print(f"  {r.file}  `{r.token}`  (修前出现在 {r.n_old_lines} 行)\n      {r.old_excerpt}")
    print("\n⚠ SAFE SIDE(#P6):只在**离开**方向可读。一个数可以被**改写**而不是删除 ——"
          "\n   命中项必须逐条人工判「被删了,还是换了说法」。没有离开 ≠ 去重是对的。")
    return D


# ============================================================================
# #173:被删掉的是限定语,还是数字
# ----------------------------------------------------------------------------
# `#172` 自己的「结构上做不到」写着:一句被删掉的**限定语**不带数字,这道闸看不见。
# 而 §2 的 η 规则说,`#167`-`#171` 这一整串的病根就是**限定语** ——
# 「a number reported without the scope over which it holds」。
#
# P6 代理账:
#   PROPERTY    一次改动没有把某个数的限定语剥掉
#   PROXY       某个被删除的行同时含**限定语标记**与**一个数量**,而那个数量在改后**仍然存在**
#   IMPLICATION 命中 => 数活着、限定它的那句话死了(**命中方向可读**)
#   WITNESS     限定语可以被**改写到别处**(移到表头、移到脚注)——
#               所以命中项必须人工确认「它是没了,还是搬家了」
#   SAFE SIDE   只在**命中**方向判。没有命中 != 限定语都还在:
#               一句不含数字的限定语(「这是上界」「仪器无法逐轮判定」)这条规则也看不见。
# ============================================================================
_QUAL = _re.compile(
    r'UNVERIFIED|bounded|upper bound|at most|only|not resolvable|cannot|structurally|'
    r'未估|不可分辨|上界|下界|保守|单切分|假阳性|结构上|不能|无法|仅|只在|方向可读|'
    r'\bn\s?=|p95|±|scope|instrument|population|噪声带|自身展布', _re.I)

def qualifiers_stripped_texts(old, new):
    """纯文本版:同一段判定逻辑,不碰 git。**正对照就跑这一支** ——
    否则对照要穿过 git/cwd/文件名三层,任何一层出错都会让对照静默失败,
    而一个失败的正对照会把仪器的零读成无罪(P5 ★)。"""
    import pandas as pd
    norm=lambda m: m.strip().replace(' ','').replace('倍','×').replace('−','-')
    newset=new.split('\n'); out=[]
    for ln in set(old.split('\n'))-set(newset):
        if not _QUAL.search(ln): continue
        nums=[norm(m) for m in _MAGNUM.findall(ln)]
        naked=[]
        for t in nums:
            homes=[l for l in newset if t in norm(l) or t in l]
            if not homes: continue
            if not any(_QUAL.search(h) for h in homes): naked.append(t)
        if naked: out.append(dict(surviving=', '.join(sorted(set(naked))[:6]),deleted_line=ln[:150]))
    return pd.DataFrame(out)


def qualifiers_stripped(rev_from, rev_to='WORKTREE', files=('README.md','README_zh.md')):
    """git 版:只负责**取文本**,判定一律交给 `qualifiers_stripped_texts`。

    #173c:第一版两条路各写一份逻辑,于是正对照走的是 git 路、**失败在管道上**,
    而我差点把它读成「逻辑不成立」。**一份逻辑,一个家。**
    返回 (DataFrame, plumbing) —— plumbing 记每个文件读到多少字符,供管道对照用。
    """
    import subprocess, pathlib as _pl, pandas as pd
    def read(rev, f):
        if rev == 'WORKTREE':
            try: return _pl.Path(f).read_text()
            except FileNotFoundError: return ''
        r = subprocess.run(['git','show',f'{rev}:{f}'], capture_output=True, text=True)
        return r.stdout if r.returncode == 0 else ''
    # 「新家」跨两份 README 找(#173b):`63e03b7` 把 README.md 从中文换成英文。
    new_all = '\n'.join(read(rev_to, g) for g in files)
    frames = []; plumbing = {}
    for f in files:
        old = read(rev_from, f); plumbing[f] = (len(old), len(new_all))
        if not old: continue
        D = qualifiers_stripped_texts(old, new_all)
        if len(D): frames.append(D.assign(file=f, rev_from=rev_from))
    out = pd.concat(frames, ignore_index=True) if frames else \
          pd.DataFrame(columns=['surviving','deleted_line','file','rev_from'])
    return out, plumbing
