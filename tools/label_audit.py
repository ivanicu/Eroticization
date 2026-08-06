"""tools/label_audit.py — 「语义在」不等于「语义对」

造于 `E02·A221·R592`(`#546e` 的第一条)。行动类型:**PRODUCTION**。
`tools/inventory.py` 只回答「这份数据有没有条目文本」。它**不检查那些文本是否描述了它所标的列**。
而本项目已经被这件事咬过:`samesexany` 的标签是
"Ever Had Sexual Contact with a Female (computed)",而它的取值是 **1=是 / 5=否 / 7=其他** ——
**标签一个字都没提这件事**,`#495b` 那次自伤就是这么来的。

**判据(先于使用写死,带正/负对照):**
   若标签里出现**取值词**(`yes`/`no`/`male`/`female`/`1=`/`0=`… 或 `(computed)`/`(recode)` 之类的加工标记),
   则该列的**实际取值集合**必须与之相容。三值:
   `INCOMPATIBLE`(标签暗示二元/命名取值,实际取值集合不是它)· `OK` · `N/A`(标签不含取值词)。
   ⚠ 与 `#546c` 同样的纪律:**只在可靠方向上判决** —— 本工具**不发「一定对」**,
   它只发 `INCOMPATIBLE`(找到了矛盾)与 `N/A`/`OK-weak`(没找到矛盾 ≠ 没有矛盾)。
   正对照:`nsfg.samesexany` 必须判 `INCOMPATIBLE`;
   负对照:`gss.year`(纯数值、标签无取值词)必须判 `N/A`,**不得**判 `INCOMPATIBLE`。
"""
import json, pathlib, re, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
BINWORD = re.compile(r'\b(yes|no|male|female|ever|whether|has|had|any)\b', re.I)
COMPUTED = re.compile(r'\(\s*(computed|recode|recoded|derived)\s*\)', re.I)


def verdict(label, values):
    """三值。只在可靠方向判决:找到矛盾才说话。"""
    vs = sorted({v for v in values if v == v})
    lab = str(label or "")
    if not lab.strip(): return "N/A", "无标签"
    binary_hint = bool(BINWORD.search(lab)) or bool(COMPUTED.search(lab))
    if not binary_hint: return "N/A", "标签不含取值词或加工标记"
    if len(vs) > 12: return "N/A", f"取值 {len(vs)} 个,非命名型,判据不适用"
    canon = {0.0, 1.0}
    if set(vs) <= canon or set(vs) <= {1.0, 2.0}:
        return "OK-weak", f"取值 {vs} 与「是/否」相容(**弱** —— 没找到矛盾 ≠ 没有矛盾)"
    return "INCOMPATIBLE", f"标签暗示是/否,实际取值 {vs} —— **不看码本就会读错**"


def main():
    import numpy as np, pandas as pd
    rows = []
    # --- NSFG(旁挂 .dct 提供标签)
    def parse(p):
        out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
        for line in open(p, errors="replace"):
            m = pat.search(line)
            if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
        return out
    NS = ROOT / "data/external/nsfg"
    LAY = parse(NS / "setup/2011_2013_FemRespSetup.dct")
    pick = ["samesexany", "oppsexany", "hadsex", "evrmarry", "cohever", "samesex", "sxok18"]
    cols = {n: LAY[n] for n in pick if n in LAY}
    buf = {n: [] for n in cols}
    for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
        for n, (s, w, _) in cols.items():
            v = line[s:s + w].strip()
            buf[n].append(float(v) if v not in ("", ".") else float("nan"))
    for n, (s, w, lab) in cols.items():
        v, why = verdict(lab, buf[n])
        rows.append(dict(source="nsfg", var=n, label=lab, verdict=v, why=why,
                         n=int(sum(1 for x in buf[n] if x == x)),
                         inclusion=["2011-2013 女性问卷", "取值取全样本非缺失", "标签来自 .dct"]))
    # --- GSS(.dta 内嵌标签)
    G = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
    vl = pd.read_stata(G, iterator=True).variable_labels()
    gpick = ["year", "sex", "homosex", "premarsx", "abany", "cappun", "grass", "fepol"]
    g = pd.read_stata(G, columns=[c for c in gpick if c in vl], convert_categoricals=False)
    for c in g.columns:
        v, why = verdict(vl.get(c, ""), g[c].dropna().unique().tolist())
        rows.append(dict(source="gss", var=c, label=str(vl.get(c, ""))[:60], verdict=v, why=why,
                         n=int(g[c].notna().sum()),
                         inclusion=["GSS 1972-2024 全样本", "取值取非缺失唯一值", "标签内嵌于 .dta"]))
    # --- 对照
    pos = next(r for r in rows if r["source"] == "nsfg" and r["var"] == "samesexany")
    neg = next(r for r in rows if r["source"] == "gss" and r["var"] == "year")
    print(f"判据对照 —— 正:nsfg.samesexany = {pos['verdict']}(必须 INCOMPATIBLE) · "
          f"负:gss.year = {neg['verdict']}(必须 N/A)")
    if not (pos["verdict"] == "INCOMPATIBLE" and neg["verdict"] == "N/A"):
        print("⛔ 判据未通过对照 —— 不输出任何一行"); sys.exit(2)
    print("✅ 判据通过对照\n")
    w = max(len(f"{r['source']}.{r['var']}") for r in rows)
    for r in rows:
        print(f"  {(r['source']+'.'+r['var']).ljust(w)}  {r['verdict']:13s} n={r['n']:6d}  "
              f"{r['label'][:40]:42s} {r['why']}")
    n_bad = sum(1 for r in rows if r["verdict"] == "INCOMPATIBLE")
    print(f"\n共 {len(rows)} 列:**不相容 {n_bad}** · 弱相容 "
          f"{sum(1 for r in rows if r['verdict']=='OK-weak')} · 判据不适用 "
          f"{sum(1 for r in rows if r['verdict']=='N/A')}")
    print("⚠ 本工具**不发「一定对」** —— `OK-weak` 只表示没找到矛盾(`#546c` 同一条纪律)。")
    out = ROOT / "E02_condemnation_is_not_rarity/A221_semantics_present_is_not_semantics_right/R592_label_value_compat/results"
    out.mkdir(parents=True, exist_ok=True)
    json.dump(dict(rows=rows, n_incompatible=n_bad,
                   controls=dict(positive="nsfg.samesexany -> INCOMPATIBLE",
                                 negative="gss.year -> N/A", passed=True),
                   note="不发「一定对」;OK-weak = 没找到矛盾 != 没有矛盾",
                   unchallenged=True), open(out / "label_audit.json", "w"), indent=1)
    print(f"wrote {out/'label_audit.json'}")


if __name__ == "__main__":
    main()
