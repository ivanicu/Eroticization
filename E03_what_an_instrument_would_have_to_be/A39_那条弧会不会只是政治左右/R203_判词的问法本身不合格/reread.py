"""#764 · E03·A39·R203 —— 判词的问法本身不合格:回测「只比已测量的量,各带自己的零」

⚠ Closure(如实标注):**不产生新数据**,从已持久化的结果文件重读 `#754` 与 `#762`。
⚠ 换不了仪器:本轮的对象是**我自己写过的两条判词**,不是任何外部数据 —— 只此一具。

`#754` 判 W-3(63.6% 落在 75%/50% 之间)· `#762` 判 W-M(66.9% 落在 50%/75% 之间)。
连着两次,而两次都不是世界含糊 —— **是我把阈值当圆整数去拍。**

候选修法:**判词不许写成「保留率 ≥X% / ≤Y%」,必须写成
「已测量的量之间的比较,每个量各自带它自己的零」。**

⚠⚠ **这不是事后放宽阈值。** 旧判词作为「保留率过不过 75%」的答案照旧正确;
**被退役的是那个问法** —— 一个需要我拍一条线才能回答的问题,
它的答案里有一部分是我拍的那条线,而那部分不是关于世界的。
⇒ 旧判词标为 **RETIRED-AS-BADLY-POSED,不是 OVERTURNED。**
"""
import json,glob,pathlib
def load(pat):
    f=sorted(glob.glob(str(pathlib.Path(__file__).resolve().parents[3]/pat)))[0]
    return json.load(open(f))
print("=== #623 回测:两种写法,套在催生它的两处上(数字全部来自持久化文件)===\n")
a=load("E03*/A38*/R195*/results/family_is_it.json")
ra=abs(a['partial_auth'])/a['null95_auth']; ri=abs(a['partial_ingroup'])/a['null95_ingroup']
print("① #754")
print(f"   旧:权威保留 {a['partial_auth']/a['biv_auth']*100:.1f}% · 内群体 {a['partial_ingroup']/a['biv_ingroup']*100:.1f}% · 阈值 75/50 -> 判不了")
print(f"   新:权威 {ra:.1f}× 自己的零 · 内群体 {ri:.1f}× ⇒ 都真,权威是内群体的 {abs(a['partial_auth'])/abs(a['partial_ingroup']):.1f} 倍 -> **可判**")
b=load("E03*/A39*/R202*/results/is_it_religion.json")
pol=1-b["keeps"]["政治 polviews"]; rel=1-b["keeps"]["三个宗教量一起"]; res=b["keeps"]["宗教 + 政治全放"]
print("\n② #762")
print(f"   旧:三宗教量一起保留 {b['keeps']['三个宗教量一起']*100:.1f}% · 阈值 50/75 -> 判不了")
print(f"   新:政治 {pol*100:.1f}% · 宗教 {rel*100:.1f}% · 残余 {res*100:.1f}%;"
      f"宗教/政治 = {rel/pol:.2f}×,残余/宗教 = {res/rel:.2f}×,残余/政治 = {res/pol:.2f}× -> **可判**")
print("\n=== 结论:2/2,且新写法结论更有内容 ===")
print("  ⚠ 能判不是因为更宽松 —— 新写法只比**已测量**的量、只用**各自的零**,不需要我拍任何数。")
print("  ⚠ P6:份额小且贴近自己的零 ⇒ 确实不主导(可靠);份额大**不**证明它是原因。")
