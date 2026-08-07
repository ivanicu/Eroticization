r"""#869 · E03·A98·R308 —— 数一数还有几条守则只靠拷贝活着(PRODUCTION)

**本轮不问任何关于人的问题。如实标 Production,不假装 Frontier。**

逻辑住在两个可复用的地方,本文件只是把它们跑一遍并存下产物:
  · `tools/exception_audit.py` —— 审计本身(**每个形状自带正控/负控**,不受控的不计入)
  · `lib/bks_items.py`        —— 本轮搬走的那一条(BKS likert 题集 + `biomale` 例外)

⚠ 详细的设计、代理账、攻击向量与边界,写在那两个文件的 docstring 与账本 `#869` 里。

⚠⚠ **仪器:本项目自己的 862 个 `.py` 源文件。而这一轮换不了仪器,理由不是「还没做」:**
**本轮的对象就是这个代码库本身** —— 「项目里有多少条只靠拷贝活着的守则」这个问题,
**只此一具仪器**,第二具仪器在定义上不存在(别人的代码库回答的是别人的问题)。
⇒ 所以本轮**不做跨仪器复制,而这是结构性的,不是省略**。
**它换来的是另一种严格:审计的每一个形状各自带正控与负控,不受控的形状不计入。**
"""
import subprocess, sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.exit(subprocess.call([sys.executable, str(ROOT / "tools/exception_audit.py")], cwd=ROOT))
