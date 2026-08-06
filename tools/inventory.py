"""tools/inventory.py — 「文件在」与「语义在」是两列独立的布尔

造于 `E02·A219·R590`(`#544` 的 NEXT)。行动类型:**PRODUCTION**。

**为什么造它:** 本项目三次撞上同一件事 —— **文件完整、可跑、无缺失,而仍然不可用,
因为它缺的不是数据,是数据的含义**(`#541c` RWAS 无题目文本 · `#541b` MSSCQ 无非性条目 ·
`#544b` MFQ 无条目)。而这三次都是**我逐个手查出来的**,一次比一次晚。

**判据(先于使用写死,并带正/负对照):**
   `has_semantics` = 该来源里存在**长度 > 30 字符的变量标签或题目文本行**。
   正对照:`MSSCQ`(码本逐题列出全文)必须判 True;
   负对照:`RWAS`(码本只写 `Q1 - Q22`)必须判 False。
   **两个对照任一不过,判据本身作废,不得输出任何一行。**

用法:`python tools/inventory.py`(无参数),输出表格并写 `results/inventory.json`。
"""
import json, pathlib, re, sys, zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
LONG = 30            # 「一句话」的最短长度,先于使用写死


# ⚠ 第三次改判据,又是被对象改的,不是被我想到的(`#546b`):
# 只认 `1.` 与 `Qn` 时,HEXACO(`EFear8 我…`)· SD3(`M1<TAB>我…`)· EQSQ(`E1. 我…`)· MIES
# 四个包全被判「语义不在」—— 而它们的题目文本都在,只是编号是**字母+数字**码。
# 抽查 `False` 抓住了它:一个刚改过的判据,它的 False 必须对着对象看一眼。
# ⚠⚠ 第四次,而这次改的不是正则,是**判据的方向**(`#546c`)。
# 实测:「没找到条目行」预测「语义不在」错了 **7 次** —— gss(.dta 内嵌)· dplace(csv 定义)·
# HEXACO/SD3/EQSQ(字母数字编号)· MIES(**JSON** `"Q1" : "…"`)· YRBS(**SAS `label` 块**)。
# 每一次我都加一条路,而对象每次都拿出第八种编码方式。
# ⇒ **没有任何有限的路径集合能许可 False。** 按 P6 的安全侧:**只在可靠方向上判决** ——
# 本工具的词汇里**删掉 False**:只有 `True`(找到了)与 `None`(没找到 = 判不了)。
# 这条规则我在第一版的 docstring 里就写了,然后连续三版都在发 False。
ITEM = re.compile(r'^\s*(?:[A-Za-z]{0,8}\d{1,3}[.):\s\t]|_column\(\d+\))\s*.*?["\s\t](\S[^"]{%d,})' % LONG)


def long_lines(text, n=LONG):
    """⚠ 第一版数的是「长行」,负对照当场失败:RWAS 的码本有散文和 TIPI 题目,
    于是「有长行」判不了「有本量表的条目」。**改成数编号条目行** ——
    `1. <文本>` / `Q7 <文本>` / `_column(n) ... "<标签>"`,文本本身须 >LONG 字符。
    判据是被它自己的负对照改掉的,不是被我想到的。"""
    return [m.group(1).strip() for m in (ITEM.match(l) for l in text.splitlines()) if m]


EMBEDDED = (".dta", ".sav", ".XPT", ".xpt")


def embedded_labels(f: pathlib.Path):
    """⚠ P6 的代理账:旁挂码本**在** ⇒ 语义在(可靠);旁挂码本**不在** ⇏ 语义不在(**不可靠**)。
    第一版只找旁挂文件,于是 gss(.dta 内嵌标签)· dplace(variables.csv)· mfq(.sav)
    全被判成「语义不在」—— **三个假阴性,方向恰好是让我少干活的那个方向。**
    本函数读内嵌元数据;读不到返回 None(判不了),**永远不返回 False**。"""
    try:
        if f.suffix == ".dta":
            import pandas as pd
            vl = pd.read_stata(f, iterator=True).variable_labels()
            return sum(1 for v in vl.values() if v and len(str(v).strip()) > LONG)
        if f.suffix == ".sav":
            import pyreadstat
            _d, m = pyreadstat.read_sav(str(f), metadataonly=True)
            return sum(1 for v in m.column_names_to_labels.values() if v and len(str(v).strip()) > LONG)
    except Exception:
        return None
    return None


def csv_definitions(d: pathlib.Path):
    """有些来源把定义放在一个普通 csv 里(如 D-PLACE 的 variables.csv)。"""
    for c in list(d.rglob("variables.csv")) + list(d.rglob("*codebook*.csv")):
        try:
            import pandas as pd
            df = pd.read_csv(c, nrows=4000)
            for col in ("definition", "title", "description", "label"):
                if col in df.columns:
                    n = int((df[col].astype(str).str.len() > LONG).sum())
                    if n >= 5: return n, c.name
        except Exception:
            pass
    return 0, None


def probe(p: pathlib.Path):
    """返回 (file_present, semantics_present, evidence)。三值:True/False/None(判不了)。"""
    if not p.exists(): return False, None, "路径不存在"
    if p.is_dir():
        txts = list(p.rglob("codebook.txt")) + list(p.rglob("*.dct")) + list(p.rglob("*.sas"))
        data = list(p.rglob("*.csv")) + list(p.rglob("*.dat")) + list(p.rglob("*.sav")) \
             + list(p.rglob("*.dta")) + list(p.rglob("*.zip")) + list(p.rglob("*.XPT"))
        if not data: return False, None, "目录内无数据文件"
        if not txts:
            n_csv, where = csv_definitions(p)
            if n_csv >= 5: return True, True, f"{len(data)} 数据文件 · 无旁挂码本,但 {where} 含 {n_csv} 条 >{LONG} 字符的定义"
            for f in data[:6]:
                if f.suffix in EMBEDDED:
                    n = embedded_labels(f)
                    if n and n >= 5: return True, True, f"{len(data)} 数据文件 · 无旁挂码本,但 {f.name} 内嵌 {n} 个 >{LONG} 字符的标签"
            return True, None, f"{len(data)} 数据文件,0 个旁挂码本,内嵌元数据也未读到 -> **判不了,不判 False**(P6:旁挂不在 ⇏ 语义不在)"
        best, where = 0, None
        for t in txts[:12]:
            try: ll = long_lines(t.read_text(errors="replace"))
            except Exception: continue
            if len(ll) > best: best, where = len(ll), t.name
        if best >= 5: return True, True, f"{len(data)} 数据文件 · 码本 {len(txts)} 个 · 最长者含 {best} 条编号条目({where})"
        n_csv, w2 = csv_definitions(p)
        if n_csv >= 5: return True, True, f"{len(data)} 数据文件 · 旁挂码本无编号条目,但 {w2} 含 {n_csv} 条定义"
        for f in data[:6]:
            if f.suffix in EMBEDDED:
                n = embedded_labels(f)
                if n and n >= 5: return True, True, f"{len(data)} 数据文件 · 旁挂码本无编号条目,但 {f.name} 内嵌 {n} 个标签"
        return True, None, f"{len(data)} 数据文件 · 码本 {len(txts)} 个 · 0 条编号条目 · 无 csv 定义 · 无内嵌标签({where})"
    if p.suffix == ".zip":
        try:
            with zipfile.ZipFile(p) as z:
                names = z.namelist()
                txt = [n for n in names if n.lower().endswith((".txt", ".dct", ".sas", ".do", ".r"))]
                if not txt: return True, None, f"{len(names)} 成员,0 个文本成员 -> **判不了**(成员可能内嵌标签,未解压)"
                best = 0
                for n in txt[:6]:
                    try: best = max(best, len(long_lines(z.read(n).decode("utf8", "replace"))))
                    except Exception: pass
                return (True, True, f"{len(names)} 成员 · 文本成员 {len(txt)} 个 · 最长者 {best} 条编号条目") \
                    if best >= 5 else (True, None, f"{len(names)} 成员 · 文本成员 {len(txt)} 个 · 0 条编号条目 -> **判不了**")
        except Exception as e:
            return True, None, f"zip 读取失败:{e}"
    return True, None, "既非目录也非 zip"


SOURCES = {
    "gss": EXT / "gss", "nsfg": EXT / "nsfg", "dplace": EXT / "dplace", "yrbs": EXT / "yrbs",
    "brfss": EXT / "brfss", "openpsych/MSSCQ": EXT / "openpsych/MSSCQ",
    "openpsych/RWAS": EXT / "openpsych/RWAS", "openpsych/HEXACO.zip": EXT / "openpsych/HEXACO.zip",
    "openpsych/SD3.zip": EXT / "openpsych/SD3.zip", "openpsych/ECR.zip": EXT / "openpsych/ECR-data-1March2018.zip",
    "openpsych/SCS.zip": EXT / "openpsych/SCS.zip", "openpsych/EQSQ.zip": EXT / "openpsych/EQSQ.zip",
    "openpsych/MIES.zip": EXT / "openpsych/MIES_Dev_Data.zip",
    "dataverse/mfq": EXT / "dataverse/mfq", "ngram": EXT / "ngram",
}
for z in sorted((EXT / "dataverse").glob("*.zip")): SOURCES[f"dataverse/{z.stem[:22]}"] = z
# :zip 未解压是「判不了」的唯一来源。解压后的目录优先于 zip 本身。
for k in list(SOURCES):
    v = SOURCES[k]
    if isinstance(v, pathlib.Path) and v.suffix == ".zip":
        d = v.with_name(v.stem + "_x")
        if d.is_dir(): SOURCES[k] = d


def main():
    # --- 判据的正/负对照,先跑,不过就不输出
    pos = probe(EXT / "openpsych/MSSCQ")
    neg = probe(EXT / "openpsych/RWAS")
    print(f"判据对照 —— 正:MSSCQ semantics={pos[1]}(必须 True) · 负:RWAS semantics={neg[1]}(必须 None —— False 已从词汇中删除)")
    if not (pos[1] is True and neg[1] is None):
        print("⛔ 判据未通过对照 —— 不输出任何一行(`#P5★`:未经正对照的判定是沉默,不是无罪)")
        sys.exit(2)
    print("✅ 判据通过对照\n")
    rows = []
    for name, p in SOURCES.items():
        f, s, ev = probe(p)
        rows.append(dict(source=name, path=str(p.relative_to(ROOT)), file_present=f,
                         semantics_present=s, evidence=ev))
    w = max(len(r["source"]) for r in rows)
    print(f"{'来源'.ljust(w)}  文件在  语义在  证据")
    for r in rows:
        s = {True: "  是  ", False: "**否(不应出现)**", None: " 判不了"}[r["semantics_present"]]
        print(f"{r['source'].ljust(w)}   {'是' if r['file_present'] else '否'}    {s}  {r['evidence']}")
    n_f = sum(1 for r in rows if r["file_present"])
    n_s = sum(1 for r in rows if r["semantics_present"] is True)
    n_no = sum(1 for r in rows if r["semantics_present"] is False)   # 恒为 0:False 已从词汇中删除
    n_un = sum(1 for r in rows if r["semantics_present"] is None)
    print(f"\n共 {len(rows)} 个来源:文件在 {n_f} · **语义在 {n_s}** · **语义不在 {n_no}** · 判不了 {n_un}")
    out = pathlib.Path(__file__).resolve().parents[1] / \
        "E02_condemnation_is_not_rarity/A219_the_dataset_on_my_own_disk/R590_inventory_tool/results"
    out.mkdir(parents=True, exist_ok=True)
    json.dump(dict(rows=rows, criterion=f"存在 >={5} 行长度 >{LONG} 字符的标签/题目文本",
                   controls=dict(positive="MSSCQ -> True", negative="RWAS -> False", passed=True),
                   totals=dict(n=len(rows), file_present=n_f, semantics=n_s,
                               no_semantics=n_no, undecidable=n_un),
                   note="「文件在」与「语义在」是两列独立布尔,不合并成一个「可用」标签",
                   unchallenged=True), open(out / "inventory.json", "w"), indent=1)
    print(f"wrote {out/'inventory.json'}")


if __name__ == "__main__":
    main()
