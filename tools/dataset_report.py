"""tools/dataset_report.py — 四级闸接成一条路径,产出一份体检报告

造于 `E02·A223·R595`(`#549` 的 NEXT)。行动类型:**PRODUCTION**。
把三件已有工具接起来,并补上第四级,对**任何一个数据目录**产出三段式报告:
**能证明什么 · 不能证明什么 · 判不了什么**。

四级:
  ① `inventory`     有没有条目文本            —— 只能证明「有」(`#546c`)
  ② `label_audit`   标签暗示的取值与实际取值相容吗 —— 只能证明「不相容」(`#547b`)
  ③ 全库扫描        该来源的不相容率            —— 需先量特异度(`#548b`)
  ④ `value_range_guard` 二元列的基础率在预注册区间吗 —— 通过不证明映射对(`#549e`)
**四级都只在一个方向上可靠。报告的第一段永远比第二、三段短,这是设计,不是缺陷。**

⚠ `P2`:大文件**不整读**。`.XPT`/`.dta` 只读表头与一小块(`chunksize`)。
用法:`python tools/dataset_report.py <目录或文件>`
"""
import json, pathlib, sys, re
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools")); sys.path.insert(0, str(ROOT))
from inventory import probe
from label_audit import verdict

CHUNK = 20000          # 只读这么多行来看取值集合 —— P2


def head_of(f: pathlib.Path):
    """返回 (列名->标签, 前 CHUNK 行的 DataFrame)。只读表头与一小块。"""
    import pandas as pd
    if f.suffix.upper() == ".XPT":
        it = pd.read_sas(f, format="xport", chunksize=CHUNK, encoding="latin-1")
        d = next(iter(it))
        return {c: "" for c in d.columns}, d          # XPT 不带长标签
    if f.suffix == ".dta":
        vl = pd.read_stata(f, iterator=True).variable_labels()
        d = next(pd.read_stata(f, chunksize=CHUNK, convert_categoricals=False))
        return vl, d
    if f.suffix == ".csv":
        d = pd.read_csv(f, nrows=CHUNK)
        return {c: "" for c in d.columns}, d
    return None, None


def report(target: pathlib.Path):
    can, cannot, undec = [], [], []
    f_present, sem, ev = probe(target)
    (can if sem is True else undec).append(
        f"① 条目文本:{'**在**' if sem is True else '**判不了**'} —— {ev}")
    if sem is not True:
        undec.append("① 只能证明「在」;没找到 ≠ 不在(`#546c`),所以下面几级只在可读到标签时有意义")
    data = sorted(target.rglob("*.XPT")) + sorted(target.rglob("*.dta")) + sorted(target.rglob("*.csv")) \
        if target.is_dir() else [target]
    data = [d for d in data if d.stat().st_size > 0][:1]
    rows = []
    if not data:
        undec.append("②③④ 没有可读的 .XPT/.dta/.csv —— 全部判不了")
    else:
        f = data[0]
        print(f"  读表头与前 {CHUNK} 行:{f.name}({f.stat().st_size/1e6:.0f} MB,**不整读**)")
        vl, d = head_of(f)
        if d is None:
            undec.append(f"②③④ 不支持的格式 {f.suffix} —— 判不了")
        else:
            appl = bad = 0
            for c in d.columns:
                v, why = verdict(vl.get(c, ""), d[c].dropna().unique().tolist())
                if v == "N/A": continue
                appl += 1
                if v == "INCOMPATIBLE":
                    bad += 1
                    rows.append(dict(var=str(c), label=str(vl.get(c, ""))[:50], why=why))
            if appl == 0:
                undec.append(f"②③ 判据适用列 0 个(该来源的标签{'为空' if not any(vl.values()) else '不含开头二元标记'})"
                             f" —— **判不了,不判「干净」**")
            else:
                rate = bad / appl
                (cannot if rate > 0 else undec).append(
                    f"②③ 判据适用 {appl} 列,**不相容 {bad}**(率 {rate:.4f})"
                    + (f";例:{[r['var'] for r in rows[:5]]}" if rows else ""))
            undec.append(f"④ 基础率守卫需要**调用者声明哪些列是二元的**,本报告不猜 —— "
                         f"跑 `Gate.value_range_guard` 时必须写 `coding_note`(`#549b`)")
    print("\n【能证明什么】"); [print("  " + x) for x in can] or print("  (无)")
    print("【不能证明什么 / 已发现的问题】"); [print("  " + x) for x in cannot] or print("  (无)")
    print("【判不了什么】"); [print("  " + x) for x in undec]
    return dict(target=str(target), can=can, cannot=cannot, undecidable=undec,
                incompatible_examples=rows[:20], chunk=CHUNK,
                note="四级都只在一个方向上可靠;第一段比后两段短是设计,不是缺陷")


if __name__ == "__main__":
    t = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data/external/brfss"
    print(f"=== 数据体检:{t} ===")
    r = report(t)
    out = ROOT / "E02_condemnation_is_not_rarity/A223_a_report_for_any_dataset/R595_dataset_report/results"
    out.mkdir(parents=True, exist_ok=True)
    json.dump(r, open(out / f"report_{t.name}.json", "w"), indent=1, ensure_ascii=False)
    print(f"\nwrote {out/('report_'+t.name+'.json')}")
