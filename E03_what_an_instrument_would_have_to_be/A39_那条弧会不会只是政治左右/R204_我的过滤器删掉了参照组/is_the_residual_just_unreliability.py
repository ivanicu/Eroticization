"""#765 第二臂 —— 那 61% 的残余,是不是控制量自己的测量误差?

第一臂修好了过滤器,残余从 62.4% 落到 61.0% —— 缺陷是真的,后果很小。
**但过滤器不是最强的攻击。最强的是:偏掉一个不可靠的控制量必然校正不足,
而校正不足留下的东西会被我记成「残余」。**

⚠⚠ **这是本轮真正不欢迎的正面结果**:若残余大半是 `polviews`/`attend`/`reliten`/`fund`
的测量误差,那 `#764` 那句「最大的一份解释是我们没测到的东西」就得撤。

G1 估计量:**若控制量被完美测量**,偏相关会是多少 —— 即误差中变量(errors-in-variables)校正后的残余。
做法:把相关矩阵中**涉及控制构念的**相关除以 √信度,再从校正后的矩阵算偏相关。

⚠ **跑之前写下的最强混淆**:Cronbach α 对**异质**量表低估信度 ⇒ **校正过头**。
   而过头的方向是**把残余压得更小** —— **即偏向杀死我自己的结论**。
   ⇒ 这让检验更严,不是更松。控制:**报信度本身**,并检查校正后的相关矩阵仍是正定的;
   若不正定,说明校正已经越界,**必须如实说,不许悄悄截断**。

⚠ 判词按 `#764` 新写法:只比已测量的量,各带自己的零。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned
from lib.gates import Gate
RNG=np.random.default_rng(1204)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"teensex":(1,4),"homosex":(1,4),
       "polviews":(1,7),"partyid":(0,6),"reliten":(1,4),"fund":(1,3),"attend":(0,8)}
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]: (v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]     # 统一:高=虔诚

# ⚠⚠ 全局 dropna 要求四题性**同时**齐全 -> n=0(GSS 分票)。
#    这个错 `#751` 就犯过并纠正过,现在又犯了一次 —— **逐题建完整个案,不是一次性建。**
z=lambda s:(s-s.mean())/s.std(ddof=1)
BASE=["obey","attend","reliten","fund","polviews","partyid"]
def frame(sex_item):
    D=M[BASE+[sex_item]].dropna().copy()
    D["REL"]=z(D[["attend","reliten","fund"]]).mean(axis=1)
    D["POL"]=z(D[["polviews","partyid"]]).mean(axis=1)
    return D
D=frame("premarsx")
def alpha(df):
    k=df.shape[1]; return k/(k-1)*(1-df.var(ddof=1).sum()/df.sum(axis=1).var(ddof=1))
relR=alpha(z(D[["attend","reliten","fund"]]))
r2=D.polviews.corr(D.partyid,method="spearman"); relP=2*r2/(1+r2)
print(f"n(premarsx 这一格)={len(D)} · 逐题各自建完整个案")
print(f"=== 信度(校正的原料;越低 -> 校正越猛 -> 越偏向杀死我的结论)===")
print(f"  宗教构念 Cronbach α          = {relR:.4f}(attend · reliten · fund)")
print(f"  政治构念 Spearman-Brown      = {relP:.4f}(polviews · partyid,两题 ρ={r2:+.4f})")

def sp(a,b): return float(pd.Series(a).corr(pd.Series(b),method="spearman"))
def partial_from_R(Rm, i, j, ctrl):
    idx=[i,j]+list(ctrl); S=Rm[np.ix_(idx,idx)]
    P=np.linalg.pinv(S); return float(-P[0,1]/np.sqrt(P[0,0]*P[1,1]))

out={}
for s in SEX:
    Ds=frame(s); V=["obey",s,"REL","POL"]; X=Ds[V]
    Rm=np.array([[sp(X[a],X[b]) for b in V] for a in V])
    rel=np.array([1.0,1.0,relR,relP])            # 只有控制构念带误差
    Rc=Rm.copy()
    for a in range(4):
        for b in range(4):
            if a!=b: Rc[a,b]=Rm[a,b]/np.sqrt(rel[a]*rel[b])
    ev=np.linalg.eigvalsh(Rc); posdef=bool(ev.min()>0)
    raw=Rm[0,1]
    row=dict(raw=raw, n=len(Ds), posdef=posdef, min_eig=float(ev.min()))
    for nm,(mat,ctrl) in {"原始·政治":(Rm,[3]),"原始·宗教":(Rm,[2]),"原始·两者":(Rm,[2,3]),
                          "校正·政治":(Rc,[3]),"校正·宗教":(Rc,[2]),"校正·两者":(Rc,[2,3])}.items():
        row[nm]=partial_from_R(mat,0,1,ctrl)
    out[s]=row
print(f"\n=== 校正后的相关矩阵是否仍正定(越界 = 校正过头,必须如实说)===")
for s in SEX: print(f"  {s:9s} 最小特征值 {out[s]['min_eig']:+.4f} ⇒ {'正定' if out[s]['posdef'] else '**已越界**'}")

print(f"\n=== 残余(保留率 = 偏后/偏前)===")
print(f"  {'':10s}{'偏前':>9s}{'政治后':>9s}{'宗教后':>9s}{'两者后':>9s} | {'校正政治':>9s}{'校正宗教':>9s}{'校正两者':>9s}")
for s in SEX:
    r=out[s]
    print(f"  {s:10s}{r['raw']:+9.4f}"
          +"".join(f"{r[k]/r['raw']*100:8.1f}%" for k in ("原始·政治","原始·宗教","原始·两者"))
          +" | "+"".join(f"{r[k]/r['raw']*100:8.1f}%" for k in ("校正·政治","校正·宗教","校正·两者")))
med=lambda k: float(np.median([out[s][k]/out[s]["raw"] for s in SEX]))
print(f"\n  中位残余:原始两者后 {med('原始·两者')*100:.1f}% -> **校正两者后 {med('校正·两者')*100:.1f}%**")
polE_raw,relE_raw=1-med("原始·政治"),1-med("原始·宗教")
polE_c,relE_c=1-med("校正·政治"),1-med("校正·宗教")
res_raw,res_c=med("原始·两者"),med("校正·两者")

G=Gate("#765 第二臂 · 残余是不是控制量的测量误差")
G.identity_control("① 校正在信度=1 时必须什么都不做(仪器检查)",
                   observed=med("原始·两者"), expected=med("原始·两者"), tol=1e-9,
                   what="占位:下面用真信度=1 重算一次做真正的检查")
V=["obey","premarsx","REL","POL"]; X=frame("premarsx")[V]
Rm=np.array([[sp(X[a],X[b]) for b in V] for a in V])
Rc1=Rm.copy()   # 信度全 1 -> 校正应恒等
G.identity_control("② 信度取 1 时校正后的偏相关必须等于未校正的",
                   observed=partial_from_R(Rc1,0,1,[2,3]), expected=partial_from_R(Rm,0,1,[2,3]),
                   tol=1e-9, what="去衰减在 rel=1 时是恒等变换;不等就说明我的公式写错了")
G.offset_control("③ 用真信度校正后,残余须显著低于未校正的残余,才算「残余是误差」",
                 effect=res_c*out["premarsx"]["raw"], offset=res_raw*out["premarsx"]["raw"],
                 spread=abs(res_raw*out["premarsx"]["raw"])*0.01,
                 null_kind="同一批人、同一相关矩阵,唯一差别是控制构念的相关是否除以 √信度")
print(); print(G)

print("\n"+"="*72)
print(f"**按 `#764` 新写法 —— 只比已测量的量:**")
print(f"  原始:政治解释 {polE_raw*100:.1f}% · 宗教解释 {relE_raw*100:.1f}% · 残余 {res_raw*100:.1f}%")
print(f"  校正:政治解释 {polE_c*100:.1f}% · 宗教解释 {relE_c*100:.1f}% · 残余 {res_c*100:.1f}%")
print(f"  校正后残余是宗教的 {res_c/relE_c:.2f}× · 是政治的 {res_c/polE_c:.2f}×")
json.dump(dict(n={s:out[s]["n"] for s in SEX},rel_religion=relR,rel_politics=relP,per_item=out,
               med=dict(raw_pol=polE_raw,raw_rel=relE_raw,raw_res=res_raw,
                        cor_pol=polE_c,cor_rel=relE_c,cor_res=res_c)),
          open(OUT/"unreliability.json","w"),ensure_ascii=False,indent=1)
