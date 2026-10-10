needsPackage "PathSignatures"
needsPackage "FourTiTwo"
debug PathSignatures -- needed to use private commands
weight = word -> (sum(toList word));

wnextLyndonWord = method();
wnextLyndonWord(Array,List,ZZ) := Array => (ar,alph,k) -> (
    revalph:= hashTable apply(0..length(alph)-1,i->(alph#i => i));
    l := toList ar;
    nl := fold((ceiling(k/weight(l))):l, (i,j)->i|j);
    if(weight(nl)>k) then (nl = nl_{0..length(nl)-1});
    while(nl_(-1) == alph#-1 and length(nl)>1) do (
        nl = nl_{0..length(nl)-2};
    );
    if(nl != {alph#-1}) then nl = nl_{0..length(nl)-2} | {alph#(revalph#(nl#-1)+1)};

    new Array from nl
);

-- lyndonWords(d,k) returns a list of all Lyndon words of length at most k in d letters

wlyndonWords = method();
wlyndonWords (List,ZZ) := (alph,k) -> (
    if(k <= 0) then error("k must be a positive integer in lyndonWords(d,k).");
    l={[alph#0]};
    while(l_(-1) != [alph#(-1)]) do (
        l = l | {wnextLyndonWord(l_(-1),alph,k)};
    );
    select(l,i->weight(i)<=k)
)


stuffle = method() -- stuffle recursive definition
stuffle (List,List, NCRing) := (l1,l2,R) -> (
    if(length(l1) == 0) then return((new Array from l2)_R);
    if(length(l2) == 0) then return((new Array from l1)_R);
    l1f := l1_{0..length(l1)-2};
    l2f := l2_{0..length(l2)-2};
    a := l1#-1;
    b := l2#-1;
    return(stuffle(l1, l2f, R) * [b]_R + stuffle(l1f,l2, R) * [a]_R + stuffle(l1f, l2f, R) * ([a + b]_R))
)

stuffle(List,NCRingElement) := (l1, t) -> (
    R := ring t;
    linExt(i -> stuffle(l1, i, R), t)
)

stuffle(NCRingElement,NCRingElement) := (t1, t2) -> (
    linExt(i -> stuffle(i,t2), t1)
)

zeroOneToN = method(); -- $\pi_Y: \mathbb Q\langle \letter 0,\letter 1\rangle \to \mathbb Q\langle\mathbb N\rangle$

zeroOneToN(List) := (i) -> (
    zo := i;
    if zo == {} then return 1;
    if zo#-1 == 1 then return 0;
    n := {};
    zz := 1;
    scan(length(zo), j-> (if zo_j == 2 then (n = n | {zz}; zz = 1) else zz += 1));
    return(n)
);

zeroOneToN(List, NCRing) := (i, A) -> (
    zo := i;
    if zo == {} then return 1_A;
    if zo#-1 == 1 then return 0_A;
    n := {};
    zz := 1;
    scan(length(zo), j-> (if zo_j == 2 then (n = n | {zz}; zz = 1) else zz += 1));
    return(new Array from n)_A
);
zeroOneToN(NCRingElement, NCRing) := (t, Ak) -> linExt(i-> zeroOneToN(i,Ak), t) ;

hoffmanExpTerm = method();
hoffmanExpTerm(List, List, NCRing) := (w, c, A) -> (
    int := reverse accumulate(reverse c, 0, plus);
    int = {toList(0..(int#0 - 1))} | apply(length(int) - 1, i -> toList(int#i..(int#(i+1) - 1)));
    return 1 / (product apply(c, j -> j !)) * product(apply(int, i -> [sum(w_i)]_A))
)

hoffmanExp = method();
hoffmanExp(List, NCRing) := (w, A) -> (
    comp := keys tally apply(compositions(length(w)), i -> select(i, j -> j!=0));
    return sum(comp, i -> hoffmanExpTerm(w, i, A))
)

hoffmanExp(NCRingElement) := (f) -> (
    linExt(w -> hoffmanExp(w, ring f), f)
)

hoffmanLogTerm = method();
hoffmanLogTerm(List, List, NCRing) := (w, c, A) -> (
    int := reverse accumulate(reverse c, 0, plus);
    int = {toList(0..(int#0 - 1))} | apply(length(int) - 1, i -> toList(int#i..(int#(i+1) - 1)));
    return (-1)^(length(w) - length(c)) / (product c) * product(apply(int, i -> [sum(w_i)]_A))
)

hoffmanLog = method();
hoffmanLog(List, NCRing) := (w, A) -> (
    comp := keys tally apply(compositions(length(w)), i -> select(i, j -> j!=0));
    return sum(comp, i -> hoffmanLogTerm(w, i, A))
)

hoffmanLog(NCRingElement) := (f) -> (
    linExt(w -> hoffmanLog(w, ring f), f)
)


-- initialization

baseField = QQ; --weight grading --Question: do we need two different weights? should we implement depth?

-- create parametrization of universal variety

---- create variables for shuffle algebra


-- strrWords = select(rWords, i -> sum(toList i) == shLevel);
-- strstR = baseField new Array from (apply(strrWords, i-> v_i))

-- strRUnivZetaI = eliminate(runivZetaI, apply(select(rWords, i-> sum(toList i)!=shLevel), i-> v_i));
-- strRUnivZetaI = sub(strRUnivZetaI, strstR);
compute = wt -> (

print("Constructing rings and maps...");

shWords = select(toList apply(flatten(apply(1..wt, i->toList((i:1)..(i:2)))), i-> new Array from i), i -> length i >= 1);


shLWords := lyndonWords(2, wt);
--shLWords = shLWords_{1..length(shLWords)-2};
shLR := baseField new Array from apply(shLWords, i->x_i) | {Degrees => apply(shLWords, i->length(i))}; --commutative ring
A2 := wordAlgebra(2, CoefficientRing => shLR); -- NCRing

shLog := sum(shLWords, i -> x_i * lieBasis(i,A2)); --concatenation log in concatenation-deshuffle Hopf algebra
univShT := sum(0..wt, k -> tensorExp(shLog,k));

-- univShT // wordFormat

A2 = wordAlgebra(2, CoefficientRing => shLR);
Ak := wordAlgebra(wt, CoefficientRing => shLR);
logReg := sum(2..wt, n-> (-1)^(n-1)/n * inner(univShT, (new Array from (((n-1) : 1) | (1:2)))_A2) * ([1]_Ak)^n); --\eqref{eq:reg_comp}
reg := sum(0..wt, i -> tensorExp(logReg, i));

stMapT := sum(select(terms (reg * zeroOneToN(univShT, Ak)), i -> sum(toList value wordString (keys coefficientHTable i)#0) == wt)); --stuffle map tensor (parametrized tensor)

stWords := toList apply( flatten apply(wt..wt, l -> select(keys tally apply(compositions(l), i-> select(i, j-> j!=0)), i-> true)), i-> new Array from i);
stR := baseField new Array from (apply(stWords, i-> v_i));

print("Computing matrix...");

stLWords := wlyndonWords(toList(1..wt),wt);
stLR := baseField new Array from apply(stLWords, i-> y_i) | {Degrees => apply(stLWords, i-> sum(toList i))};

mons := flatten entries basis(wt,stLR);

stuffleOut := pol -> (polh := applyKeys(standardForm(pol - part(0,pol)), i-> applyKeys(i,j-> last baseName stLR_j)); sum(pairs polh, (term,coef) -> coef * fold(flatten apply(pairs term, (word,ex)-> toList (ex:(word_Ak))),(i,j)->stuffle(i,j))) + (sub(part(0,pol), QQ))_Ak);

stls := apply(mons, stuffleOut);
stlsIm := apply(stls, i-> inner(stMapT,i));

zeroOneSub := apply(gens shLR, i-> zeroOneToN(apply(reverse(toList(last baseName i)),j->3-j)));
zeroOneSubV := apply(zeroOneSub, i-> if(class i === List) then y_(reverse new Array from i) else 0);
zeroOneMap := map(stLR,shLR,zeroOneSubV);
auto := apply(stlsIm, i-> zeroOneMap(i));
autoM := sub(matrix apply(auto, i->flatten entries (coefficients(i,Monomials => mons))#1),QQ);

isZeroMon := i -> ((coefficients(i, Variables => {y_[1]}))#0)_0_0 != 1;

print("Computing toric variety...");

monvars := apply(mons,i->z_i);
monind = select(apply(length(mons),i->if(isZeroMon(mons#i)) then -1 else i),j->j!=-1);
monnind = select(apply(length(mons),i->if(isZeroMon(mons#i)) then i else -1),j->j!=-1);
monR := QQ new Array from monvars_monind;
monvars = gens monR;
exps := transpose matrix apply(monvars, i-> (exponents last baseName i)#0);
autoMsm := (autoM)_monind;

I := ideal mingens(toricMarkov(exps,monR));

print("Intersecting variety with image under isomorphism...");
autoMap := map(monR,monR,(entries (autoMsm * (vector monvars)))_monind);
imI := autoMap(I) + ideal ((transpose matrix (autoMsm * (vector monvars)))_monnind);
zetaI := imI + I;
return((zetaI,autoM));
);
