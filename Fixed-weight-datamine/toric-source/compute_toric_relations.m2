needsPackage "FourTiTwo"

-- Implementation of Duval's algorithm. Given d and k, nextLyndon(w,d,k) creates the next Lyndon word of length at most k in d letters after w in lexicographical order
weight = word -> (sum(toList word));

wnextLyndonWord = method();
wnextLyndonWord(Array,List,ZZ) := Array => (ar,alph,k) -> (
    l := toList ar;
    nl := fold((i,j)->i|j, (ceiling(k/weight(l))):l);
    if(weight(nl)>k) then (nl = nl_{0..length(nl)-1});
    while(nl_(-1) == alph#-1 and length(nl)>1) do (
        nl = nl_{0..length(nl)-2};
    );
    if(nl != {alph#-1}) then nl = nl_{0..length(nl)-2} | {alph#(position(alph, a -> a == nl#-1)+1)};

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

j = 10
lynds = sort(wlyndonWords({2,3},j),i->length(i))

zetaRing = QQ new Array from apply(lynds, i -> x_i) | {Degrees => apply(lynds, i->weight(i))}
mons = flatten entries basis(j,zetaRing)
A = transpose matrix apply(mons, i-> (exponents i)#0)

monRing = QQ new Array from apply(mons, i -> z_i)
relExps = entries toricMarkov(A)

posPart = v -> apply(v, i -> max(0,i))
negPart = v -> apply(v, i -> - min(0,i))


jsonString = "{\n";
jsonString = jsonString | "\"vars\": " | toString(new Array from lynds) | ",\n"
jsonString = jsonString | "\"mons\": " |toString(new Array from apply(entries transpose A, i-> new Array from i)) | ",\n"
jsonString = jsonString | "\"binoms\": [\n";
binomsJson = apply(relExps, i -> "{\n\"mon1\": " | toString(new Array from posPart(i)) | ",\n" | "\"mon2\": " | toString(new Array from negPart(i)) | "\n},\n")
binomsJson = concatenate(binomsJson)
jsonString = jsonString | (if #binomsJson == 0 then "" else substring(0, #binomsJson-2, binomsJson)) | "]\n}"

"toric-relations-weight-" | toString j | ".json" << jsonString << endl << close
