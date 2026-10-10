needsPackage "PathSignatures"
debug PathSignatures

-- Implementation of Duval's algorithm. Given d and k, nextLyndon(w,d,k) creates the next Lyndon word of length at most k in d letters after w in lexicographical order
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

dimen = k -> (
lwords = wlyndonWords({2,3},k);
lwords = sort(lwords,i->length(i));
lwords = select(lwords,i->sum(toList i)<=k);
R = QQ new Array from apply(lwords, i-> x_i) | {Degrees => apply(lwords,i->sum(toList i))};
mons = flatten entries basis(k,R);

A = transpose matrix apply(mons,i->(exponents(i))#0);

rank A
)

apply(2..30,i->dimen(i))

needsPackage "FourTiTwo"

Q = QQ[y_1..y_(length mons)]; I = toricMarkov(A,Q);

degree I

dimen(18); transpose A
"mymatrix.txt" << toString transpose A << endl << close

dimen(21); Q = QQ[y_1..y_(length mons)]; I = toricMarkov(A,Q);
P = convexHull(A)
isNormal P

time volume(P) * (rank A-1)!

needsPackage "Polyhedra"

k=23
d = dimen(k);
Simplex = convexHull diagonalMatrix apply(select(lwords,i->weight(i)!=k-1), i-> k/(weight(i)));
latticeVolume (convexHull latticePoints Simplex)
round(volume Simplex * (d-1)!)


m = map(R,Q,mons);
I = kernel m;
dim I
degree I
