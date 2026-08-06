"""E03·A18·R692 —— 修单位重跑:段级并集下,重引还分得开随机吗

**类型:FRONTIER**。`#655` 是 UNVERIFIED:句级下真实重引与随机配对都是 1.000,**分不开**。
诊断是**单位错了** —— 限定词住在条目的表和上下文里,不住在带数字的那一句里。

⚠ **`#111c`:这是这一问题的第一次修正后重试。**
**预注册的硬停止**:段级下若真/随机**仍然相同**,⇒ **这个问题在账本这个载体上判不了,
写进「做不到」并停,不试第三次。**

修法(`#655` 写死的):
① 单位从「句」改成「段」——原侧 = 源条目里**所有提到该数字的句子**的限定词并集;
   引侧 = 引用条目里**那一段**的限定词并集。
② **offset 必须用同样的段级并集算**(`#655` 的⑤:段级并集会让分母变大,
   **而分母变大会把丢失率系统性压低**;若 offset 还用句级,我就是拿变宽的分母去比没变的基线)。
③ 正对照仍是 `#653`→`#542`,段级下必须仍被抓到。
④ **判据只留能分辨的那一刀**:真实丢失率中位比 offset **低 >= 0.15** -> 有分辨力,然后才读中位;
   **低于 0.15 -> 判不了 -> 停。**
IMPOSSIBLE(不写 planned):措辞性丢失查不到 · 只查账本不查页面 · 词表是我列的 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
TXT=pathlib.Path("RETRACTIONS.md").read_text()
ENT={}
for m in re.finditer(r'^## Entry (\d+)[^\n]*\n(.*?)(?=^## Entry |\Z)', TXT, re.M|re.S):
    ENT[int(m.group(1))]=m.group(2)
print(f"=== 硬规则①:先打印 ===\n  账本条目 = **{len(ENT)}**  编号 {min(ENT)}–{max(ENT)}")

QUAL = {
 "绝对值": r'\|ρ\||\|r\||\|rho\||绝对值|abs\(',
 "中位":   r'中位|median',
 "逐年":   r'逐年|per-year|年内|within.year',
 "归一":   r'归一|normali[sz]ed|/ *r_?max|天花板',
 "样本量": r'\bn\s*=|\bn\s*=\s*\d|样本量',
 "区间":   r'CI|区间|\[[-+−]?\d*\.\d+ *, *[-+−]?\d*\.\d+\]|95%',
 "留出":   r'留出|held.out|半样本|split.half',
 "逐对":   r'逐对|每一对|per-pair|每对',
 "块重抽": r'块 ?bootstrap|block bootstrap|bootstrap',
 "秩":     r'秩相关|Spearman|Kendall|rank',
 "对齐":   r'对齐|align|反向计分|反向措辞',
 "生/同号": r'未归一|生相关|生阈|同号|signed|raw',
 "控制":   r'正对照|安慰剂|置换零|positive control|placebo',
 "对数":   r'\d+ *对\b|npairs|\d+ 对',
}
NUM=re.compile(r'(?<![\w.])([0-9]+\.[0-9]{3,4})(?![\w])')
def quals(s): return {k for k,p in QUAL.items() if re.search(p,s,re.I)}
def sents(body):
    return [x.strip() for x in re.split(r'[。\n]|(?<=[.;])\s', body) if x.strip()]

# ── 预注册⑤:先在正对照上验证词表 ─────────────────────────────────────────
print("\n=== 预注册⑤:先验证词表能不能抓到 `#653`→`#542` 的 `|·|` ===")
src=[s for s in sents(ENT[542]) if "0.5789" in s]
dst=[s for s in sents(ENT[653])+sents(ENT[654]) if "0.579" in s or "0.5789" in s]
print(f"  `#542` 里含 0.5789 的句子 {len(src)} 句 · 后来条目里重引它的 {len(dst)} 句")
if src: print(f"    原句限定 = {sorted(quals(src[0]))}")
for s in dst[:3]: print(f"    引句限定 = {sorted(quals(s))}   « {s[:70]} »")
ok_vocab = bool(src) and ("绝对值" in quals(src[0])) and any("绝对值" not in quals(s) for s in dst)
print(f"  ⇒ 词表抓得到那次丢失? **{ok_vocab}**" + ("" if ok_vocab else "  ⚠ 抓不到 ⇒ 先补词表,不许跑全库"))

# ── 段级:全库配对(修法①)────────────────────────────────────────────────
def paras(body):
    return [p for p in re.split(r'\n\s*\n', body) if p.strip()]
def src_union(x, num):
    """原侧 = 源条目里**所有提到该数字**的句子的限定词并集"""
    u=set()
    for s in sents(ENT[x]):
        if num in s: u |= quals(s)
    return u
def dst_union(y, x, num):
    """引侧 = 引用条目里**那一段**(含该反向引用且含该数字)的限定词并集"""
    u=set(); found=False
    for p in paras(ENT[y]):
        if f"#{x}" in p and num in p:
            u |= quals(p); found=True
    return u, found

pairs=[]
for y,bd in ENT.items():
    refs={int(r) for r in re.findall(r'#(\d{3})', bd)} & set(ENT)
    for x in refs:
        if x>=y: continue
        ynums=set(NUM.findall(bd)); xnums=set(NUM.findall(ENT[x]))
        for num in sorted(ynums & xnums):
            q0=src_union(x,num)
            q1,found=dst_union(y,x,num)
            if not found: continue
            if not q0: pairs.append(dict(src=x,dst=y,num=num,lost=None)); continue
            pairs.append(dict(src=x,dst=y,num=num,n_src=len(q0),
                              lost=len(q0-q1)/len(q0),lost_set=sorted(q0-q1)))
            break
measur=[p for p in pairs if p["lost"] is not None]
noq=[p for p in pairs if p["lost"] is None]
print(f"\n=== 段级配对 ===\n  {len(pairs)} 对 · **可量 {len(measur)}** · 原侧零限定(量不了) {len(noq)}")
L=[p["lost"] for p in measur]
med=float(np.median(L)) if L else np.nan
print(f"  原侧限定数中位 = {np.median([p['n_src'] for p in measur]):.1f}(句级时多为 1–3)")
print(f"  **段级丢失率:中位 {med:.3f} · 均值 {np.mean(L):.3f} · 无损 {sum(1 for v in L if v==0)}/{len(L)}**")
from collections import Counter
C=Counter(w for p in measur for w in p.get("lost_set",[]))
print("  最常被丢掉的限定:"+" · ".join(f"{k} {v}" for k,v in C.most_common(8)))

# ── offset:同样的段级并集(修法②)───────────────────────────────────────
allparas=[p for b in ENT.values() for p in paras(b)]
allents=list(ENT)
def placebo(seed):
    rng=np.random.default_rng(seed); out=[]
    for p in measur:
        x2=allents[rng.integers(0,len(allents))]
        q0=src_union(x2,p["num"]) or quals(ENT[x2][:1200])
        q1=quals(allparas[rng.integers(0,len(allparas))])
        if q0: out.append(len(q0-q1)/len(q0))
    return float(np.median(out)) if out else np.nan
pl=float(np.median([placebo(s) for s in SEEDS]))
print(f"\n=== 控制 ===\n  offset(同样的段级并集,随机配)= **{pl:.3f}**")
print(f"  正对照 词表抓到 `#653`->`#542` = **{ok_vocab}**")
gap=pl-med
print(f"  **分辨力 = offset − 真实 = {gap:+.3f}**(判据 >= 0.15)")
from lib.gates import Gate
G=Gate("段级并集下,重引还分得开随机吗")
p1=G.positive_control("词表必须抓到 #653->#542 丢掉的 `|·|`",planted=float(1.0 if ok_vocab else 0.0),floor=0.5,spread=0.01)
p2=G.offset_control("真实重引的丢失率必须比随机配对低 >= 0.15",effect=med,offset=pl,spread=0.02,
                    null_kind="同样的段级并集,随机配一条源条目与一段引文 —— 系统性基线,不该趋零")
if p1 and p2 and gap>=0.15:
    verdict=f"**有分辨力(差 {gap:+.3f})⇒ 段级丢失率中位 = {med:.3f}**"
elif p1 and p2:
    verdict=f"**判不了 —— 分辨力 {gap:+.3f} < 0.15。按 `#655` 的硬停止:写进「做不到」并停,不试第三次。**"
else:
    verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · offset {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(entries=len(ENT),pairs=len(pairs),measurable=len(measur),no_qual=len(noq),
               median_loss=med,offset=pl,gap=gap,vocab_ok=bool(ok_vocab),verdict=verdict,
               most_lost=C.most_common(10),
               worst=sorted(measur,key=lambda p:-p["lost"])[:10],unchallenged=True),
          open(OUT/"lossy_requote_para.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'lossy_requote_para.json'}")
