import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A134 R406 -- guard 22:一条只有一个点的曲线,会让下游每一个判据都通过

`#361b`:`R405` 的设计退化成一个点,而**注册的 kill 与 guard 21 都在那个点上通过了**。
**`#296b` 挡的是一个**数**越出定义域;守卫库里没有一个挡**设计**退化成一个点的。**

⚠ **Closure(工具),如实标注。**
"""
import pandas as pd,hashlib
from lib.gates import Gate

CASES=[('★ 正对照:`#404` 的六点覆盖曲线',[4,6,8,10,12,16],None,3,True),
       ('★ 正对照:`#390` 的三点节点扫描',[0.33,0.50,0.67],None,3,True),
       ('★ 负对照:`#405` 的一点曲线',[1],None,3,False),
       ('负对照:x 全相同',[8,8,8,8],None,3,False),
       ('负对照:非数',None,None,3,False),
       ('负对照:x/y 长度不一致',[1,2,3],[0.1,0.2],3,False)]
g=Gate('guard 22 的验收'); res=[]
for nm,xs,ys,mp,want in CASES:
    got=g.curve_has_enough_points(nm,xs,ys=ys,min_points=mp)
    res.append(dict(v_case=nm[:40],want=want,got=got,ok=(got==want)))
print(g)
T=pd.DataFrame(res)
print(f"\n六个用例:**{int(T.ok.sum())}/{len(T)}** 与预期一致")
T.to_csv(pathlib.Path(__file__).parent/'results'/'guard22.csv',index=False)
g2=Gate('guard 22 本身')
g2.asserted('★ 六个用例全部与预期一致',bool(T.ok.all()),f"{int(T.ok.sum())}/{len(T)}")
g2.asserted('⚠ 边界:它检的是**有没有曲线**,不检曲线**对不对**',True,
            '一条有六个点但每个点都算错的曲线,它照样放行')
print(g2)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
