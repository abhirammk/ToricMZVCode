needsPackage "PathSignatures"
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

baseField = QQ;
shLevel = stLevel = 6; --weight grading --Question: do we need two different weights? should we implement depth?

-- create parametrization of universal variety

---- create variables for shuffle algebra


-- strrWords = select(rWords, i -> sum(toList i) == shLevel);
-- strstR = baseField new Array from (apply(strrWords, i-> v_i))

-- strRUnivZetaI = eliminate(runivZetaI, apply(select(rWords, i-> sum(toList i)!=shLevel), i-> v_i));
-- strRUnivZetaI = sub(strRUnivZetaI, strstR);

shWords = select(toList apply(flatten(apply(1..shLevel, i->toList((i:1)..(i:2)))), i-> new Array from i), i -> length i >= 1);
shR = baseField new Array from apply(shWords, i->u_i); --commutative ring
varWordTable = hashTable apply(shWords, i-> (u_i, i));

shLWords = lyndonWords(2, shLevel);
shLR = baseField new Array from apply(shLWords, i->x_i); --commutative ring
A2 = wordAlgebra(2, CoefficientRing => shLR); -- NCRing

shLog = sum(shLWords, i -> x_i * lieBasis(i,A2)); --concatenation log in concatenation-deshuffle Hopf algebra
univShT = sum(0..shLevel, k -> tensorExp(shLog,k));

univShParam = tensorParametrization(univShT, VarWordTable => varWordTable); -- parametrization of projected shuffle variety

-- create map to stuffle algebra

A2 = wordAlgebra(2, CoefficientRing => shR);
Ak = wordAlgebra(shLevel, CoefficientRing => shR)
shRT = sum(shWords, i-> u_i * i_A2) + 1; --generic element of weight $\leq k$ of concatenation-deshuffle Hopf algebra
logReg = sum(2..stLevel, n-> (-1)^(n-1)/n * inner(shRT, (new Array from (((n-1) : 1) | (1:2)))_A2) * ([1]_Ak)^n); --\eqref{eq:reg_comp}
reg = sum(0..stLevel, i -> tensorExp(logReg, i));

stMapT = sum(select(terms (reg * zeroOneToN(shRT, Ak)), i -> sum(toList value wordString (keys coefficientHTable i)#0) == stLevel)); --stuffle map tensor (parametrized tensor)

stWords = toList apply( flatten apply(stLevel..stLevel, l -> select(keys tally apply(compositions(l), i-> select(i, j-> j!=0)), i-> true)), i-> new Array from i);
stR = baseField new Array from (apply(stWords, i-> v_i)); --commutative ring
--| {Degrees => apply(stWords, i-> sum(toList i))});
varWordTable = hashTable apply(stWords, i -> (v_i, i));

stMapParam = tensorParametrization(stMapT, VarWordTable => varWordTable);

-- create parametrization of stuffle algebra

stLWords = wlyndonWords(toList(1..stLevel),stLevel) -- fix this!
stLR = baseField new Array from apply(stLWords, i-> y_i) | {Degrees => apply(stLWords, i-> sum(toList i))};

Ak = wordAlgebra(shLevel, CoefficientRing => stLR); -- NCRing
setWeights(Ak, toList(1..shLevel)); 

shLog = sum(stLWords, i -> y_i * lieBasis(i,Ak)); -- overwrites previous definition
univShT = sum(1..shLevel, k -> tensorExp(shLog,k));

univStT = sum(stWords, i -> inner(univShT, hoffmanLog(i_Ak)) * i_Ak);
univStParam = tensorParametrization(univStT, VarWordTable => varWordTable); --parametrization of projected stuffle variety

-- compute varieties and intersect

univStI = preimage(univStParam, ideal(y_[1]))

-- el1 = map(stLR, stLR, {0} | (gens stLR)_{1..length(gens stLR)-1});
-- elunivStParam = el1 * univStParam;
-- needsPackage "MultigradedImplicitization"
-- -- bnd = 2;
-- munivStI = ideal flatten values componentsOfKernel(, elunivStParam)

univShI = preimage(univShParam * stMapParam,ideal(x_[1],x_[2]))

-- el12 = map(shLR, shLR, {0,0} | (gens shLR)_{2..length(gens shLR)-1});
-- elunivShParam = el12 * univShParam * stMapParam;
-- munivShI = ideal flatten values componentsOfKernel(14, elunivShParam);

univZetaI = univStI + univShI; -- variety $\mathcal{Z}_k$
-- univZetaI = sub(munivStI, stR) + sub(munivShI, stR);
hunivZetaI = eliminate(univZetaI, apply(select(stWords, i-> i#0 == 1), i-> v_i)); -- honest zeta values
hWords = select(stWords, i-> i#0 != 1);
hstR = baseField new Array from (apply(hWords, i-> v_i)); -- smaller commutative ring with only admissable words
hunivZetaI = sub(hunivZetaI, hstR)


-- compute matrix

univShLin = select(flatten entries gens univShI, i->degree i == {1});
univShLin = matrix apply(univShLin, i-> flatten entries (coefficients(i,Monomials=>gens stR))#1);
univShLin = transpose sub(univShLin, QQ);
(PSh,LSh,USh) = LUdecomposition univShLin

univStLin = select(flatten entries gens univStI, i->degree i == {1});
univStLin = matrix apply(univStLin, i-> flatten entries (coefficients(i,Monomials=>gens stR))#1);
univStLin = transpose sub(univStLin, QQ)
(PSt,LSt,USt) = LUdecomposition univStLin

S = QQ[x_(1,1)..x_(#gens stR, #gens stR)];
A = genericMatrix(S,#gens stR, #gens stR);
myeq = A*univStLin - sub(univShLin,S);
eqI = ideal(myeq) + ideal apply(1..512, i-> (gens S)_i);


