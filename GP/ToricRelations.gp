default(realprecision, 1200);
default(parisizemax, 7g);

checknat(n, name) = {
  \\ Require a nonnegative integer.
  if(type(n) != "t_INT" || n < 0, error(Str(name, " must be a nonnegative integer")));
  return(n);
};

varvec(var, n0, n1) = {
  \\ Construct quoted symbolic variables from a prefix and an index range.
  checknat(n0, "n0");
  if(type(n1) != "t_INT", error("n1 must be an integer"));
  if(n1 < n0, return([]));
  return(vector(n1-n0+1, i, eval(Str("'", var, n0+i-1))));
};

monomial(vec, exps) = {
  \\ Evaluate one monomial given its exponent vector.
  if(#vec != #exps, error("Incompatible numbers of variables and exponents"));
  my(term = 1);
  for(i = 1, #vec,
    if(exps[i], term *= vec[i]^exps[i]);
  );
  return(term);
};

monomialvec(vec, exps) = {
  \\ Evaluate a vector of monomials in a common set of variables.
  return(apply(e -> monomial(vec, e), exps));
};

monomiallabel(labels, exps) = {
  \\ Format a monomial from textual factor labels and integer exponents.
  if(#labels != #exps, error("Incompatible labels and exponents"));
  my(terms = List());
  for(i = 1, #labels,
    if(exps[i], listput(~terms, if(exps[i] == 1, labels[i], Str(labels[i],"^",exps[i]))));
  );
  return(if(#terms,strjoin(Vec(terms),"*"),"1"));
};

relrank(rels) = {
  \\ Compute the exact rational rank of integer coefficient columns.
  if(#rels == 0, return(0));
  return(matrank(Mat(apply(c -> Vec(c)~, rels))));
};

partition(n, a, b) = {
  \\ List multiplicity choices satisfying a*x+b*y=n.
  checknat(n, "n"); checknat(a, "a"); checknat(b, "b");
  if(a == 0 || b == 0 || a == b, error("Use distinct integers"));
  my(res = List(), rem, x);
  for(y = 0, n\b,
    rem = n-b*y;
    if(rem % a == 0,
      x = rem/a;
      listput(~res, concat(vector(x, i, a), vector(y, i, b)));
    );
  );
  return(Vec(res));
};

perms(x, y, a, b) = {
  \\ List distinct words with x copies of a and y copies of b.
  checknat(x, "x"); checknat(y, "y");
  if(a == b, return([vector(x+y, i, a)]));
  if(x == 0, return([vector(y, i, b)]));
  if(y == 0, return([vector(x, i, a)]));
  my(res = List(), taila = perms(x-1, y, a, b), tailb = perms(x, y-1, a, b));
  for(i = 1, #taila, listput(~res, concat([a], taila[i])));
  for(i = 1, #tailb, listput(~res, concat([b], tailb[i])));
  return(Vec(res));
};

parts(n, a = 2, b = 3) = {
  \\ List all ordered words of weight n over two distinct positive integers.
  checknat(n, "n"); checknat(a, "a"); checknat(b, "b");
  if(a == 0 || b == 0 || a == b, error("Use distinct integers"));
  my(res = List(), tail, rem);
  for(y = 0, n\b,
    rem = n-b*y;
    if(rem % a == 0,
      tail = perms(rem/a, y, a, b);
      for(i = 1, #tail, listput(~res, tail[i]));
    );
  );
  return(Vec(res));
};

compositions(n) = {
  \\ Enumerate compositions by increasing binary cut mask.
  checknat(n, "n");
  if(n == 0, return([[]]));
  my(res = List(), word, last);
  for(mask = 0, 2^(n-1)-1,
    word = List(); last = 0;
    for(j = 1, n-1,
      if(bittest(mask, j-1), listput(~word, j-last); last = j);
    );
    listput(~word, n-last);
    listput(~res, Vec(word));
  );
  return(Vec(res));
};

admissible(k) = {
  \\ List all convergent MZV indices of weight k in binary-cut order.
  checknat(k, "k");
  if(k < 2, return([]));
  my(res = compositions(k-1));
  for(i = 1, #res, res[i][1]++);
  return(res);
};

getexponents(n, deg) = {
  \\ Enumerate nonnegative exponent vectors of length n and sum deg.
  checknat(n, "n"); checknat(deg, "deg");
  if(deg == 0, return([vector(n)]));
  if(n == 0, return([]));
  if(n == 1, return([[deg]]));
  my(res = List(), tails);
  for(e = 0, deg,
    tails = getexponents(n-1, deg-e);
    for(j = 1, #tails, listput(~res, concat([e], tails[j])));
  );
  return(Vec(res));
};

nonlin(vec, deg) = {
  \\ Generate homogeneous monomials without identifying equal images.
  return(monomialvec(vec, getexponents(#vec, deg)));
};

veccmp(u, v) = {
  \\ Compare integer words lexicographically, including unequal lengths.
  for(i = 1, min(#u,#v),
    if(u[i] < v[i], return(-1));
    if(u[i] > v[i], return(1));
  );
  return(sign(#u-#v));
};

islyndon(w, flag = 0) = {
  \\ Test strict or non-strict minimality among cyclic rotations.
  if(flag != 0 && flag != 1, error("flag must be 0 or 1"));
  my(n = #w, rot, cmp);
  if(n == 0, return(0));
  for(i = 1, n-1,
    rot = concat(w[i+1..n], w[1..i]);
    cmp = veccmp(w, rot);
    if(cmp > 0 || (flag == 0 && cmp == 0), return(0));
  );
  return(1);
};

lyndonwords(k, flag = 0) = {
  \\ Enumerate weight-k Lyndon words over the ordered alphabet 2<3.
  checknat(k, "k");
  if(flag != 0 && flag != 1, error("flag must be 0 or 1"));
  return(vecsort(select(w -> islyndon(w,flag), parts(k,2,3)), veccmp));
};

lyndonarray(k) = {
  \\ Group strict Lyndon words by weights 1 through k.
  checknat(k, "k");
  return(vector(k, i, lyndonwords(i)));
};

lwto(k) = {
  \\ Flatten the strict Lyndon words of weights at most k.
  return(concat(lyndonarray(k)));
};

lwbasis(k, flag = 0) = {
  \\ Name the selected weight-k Lyndon coordinates without asserting independence.
  return(varvec('LW, 1, #lyndonwords(k,flag)));
};

totlyn(k) = {
  \\ Count strict Lyndon words of weights at most k.
  return(#lwto(k));
};

weightwalk(weights, pos, remaining) = {
  \\ Recursively distribute remaining weight among the remaining generators.
  if(remaining == 0, return([vector(#weights-pos+1)]));
  if(pos > #weights, return([]));
  my(res = List(), tails);
  for(e = 0, remaining\weights[pos],
    tails = weightwalk(weights, pos+1, remaining-e*weights[pos]);
    for(j = 1, #tails, listput(~res, concat([e], tails[j])));
  );
  return(Vec(res));
};

weightexponents(weights, k) = {
  \\ Enumerate monomials of weighted degree k with positive integer weights.
  checknat(k, "k");
  for(i = 1, #weights,
    checknat(weights[i], "weight");
    if(weights[i] == 0, error("Weights must be positive"));
  );
  return(weightwalk(weights, 1, k));
};

weightpartition(weights, exps) = {
  \\ Expand a weighted exponent vector into its sorted factor weight partition.
  if(#weights != #exps, error("Incompatible weights and exponents"));
  return(vecsort(concat(vector(#weights, i, vector(exps[i], j, weights[i])))));
};

lyndondata(k) = {
  \\ Build a duplicate-free fixed-weight coordinate model and its index dictionary.
  checknat(k, "k");
  if(k < 2, error("The coordinate weight must be at least 2"));
  my(arr = lyndonarray(k), words, vars, weights, exps, pats, groups);
  words = concat(arr);
  vars = concat(vector(k, i, varvec(Str("L",i,"n"),1,#arr[i])));
  weights = apply(vecsum, words);
  exps = weightexponents(weights, k);
  pats = Set(apply(e -> weightpartition(weights,e), exps));
  pats = vecsort(pats, (u,v) -> if(#u == #v, veccmp(u,v), sign(#u-#v)));
  groups = vector(#pats, i, select(e -> weightpartition(weights,e) == pats[i], exps));
  exps = concat(groups);
  return(Map([
    "weight", k;
    "words", words;
    "vars", vars;
    "weights", weights;
    "exponents", exps;
    "monomials", monomialvec(vars,exps);
    "partitions", pats;
    "groups", apply(g -> monomialvec(vars,g), groups)
  ]));
};

getlyndonarray(k) = {
  \\ Return the fixed-weight monomials grouped by factor weight partition.
  return(mapget(lyndondata(k), "groups"));
};

exponentimage(exps, powers) = {
  \\ Add coordinate exponent vectors with the specified multiplicities.
  if(#exps != #powers || #exps == 0, error("Invalid coordinate exponent data"));
  my(res = vector(#exps[1]));
  for(i = 1, #exps,
    if(powers[i], res += powers[i]*exps[i]);
  );
  return(res);
};

productdata(exps, deg) = {
  \\ Group degree-d monomials by exact image exponents and record a sparse kernel basis.
  checknat(deg, "deg");
  if(#exps == 0, error("At least one coordinate is required"));
  my(powers = getexponents(#exps,deg));
  my(images = apply(e -> exponentimage(exps,e), powers));
  my(seen = Map(), distinct = List(), reps = List(), pairs = List());
  my(ids = vector(#powers), pos);
  for(i = 1, #images,
    if(mapisdefined(seen, images[i], &pos),
      ids[i] = pos;
      listput(~pairs, [reps[pos],i]),
      listput(~distinct, images[i]);
      listput(~reps, i);
      ids[i] = #distinct;
      mapput(~seen, images[i], #distinct);
    );
  );
  return(Map([
    "coordinates", #exps;
    "degree", deg;
    "powers", powers;
    "images", images;
    "distinct", Vec(distinct);
    "representatives", Vec(reps);
    "ids", ids;
    "pairs", Vec(pairs);
    "products", #powers;
    "imagecount", #distinct;
    "rank", #pairs
  ]));
};

formalrels(k, deg) = {
  \\ Compute exact degree-d relations of the weight-k monomial parametrization.
  my(model = lyndondata(k));
  my(data = productdata(mapget(model,"exponents"),deg));
  mapput(~data, "weight", k);
  mapput(~data, "model", model);
  mapput(~data, "status", "Exact formal relations");
  return(data);
};

formalpolys(data) = {
  \\ Convert the sparse formal basis into polynomials in X1 through XN.
  my(vars = varvec('X,1,mapget(data,"coordinates")));
  my(slots = monomialvec(vars,mapget(data,"powers")));
  return(apply(p -> slots[p[1]]-slots[p[2]], mapget(data,"pairs")));
};

formalvectors(data) = {
  \\ Materialize dense integer coefficient columns only when explicitly requested.
  my(n = mapget(data,"products"));
  return(apply(p -> vector(n,i,(i == p[1])-(i == p[2]))~, mapget(data,"pairs")));
};

verifyformal(data) = {
  \\ Verify the stored exponent collisions and the sparse basis cardinality.
  my(images = mapget(data,"images"), pairs = mapget(data,"pairs"));
  if(#pairs != #images-#Set(images), return(0));
  if(#Set(apply(p -> p[2],pairs)) != #pairs, return(0));
  for(i = 1, #pairs,
    if(pairs[i][1] >= pairs[i][2], return(0));
    if(images[pairs[i][1]] != images[pairs[i][2]], return(0));
  );
  return(1);
};

formalsummary(data) = {
  \\ Print the exact coordinate, image, and homogeneous-kernel counts.
  printf("Weight %d, degree %d: %d coordinates, %d monomials, %d images, %d relations.\n",
    mapget(data,"weight"), mapget(data,"degree"), mapget(data,"coordinates"),
    mapget(data,"products"), mapget(data,"imagecount"), mapget(data,"rank"));
};

printformal(data, limit = 18) = {
  \\ Print a bounded number of exact formal binomial relations.
  checknat(limit, "limit");
  my(pols = formalpolys(data), n);
  n = if(limit == 0, #pols, min(limit,#pols));
  for(i = 1, n, printf("Formal relation %d: %s = 0\n", i, Str(pols[i])));
  if(n < #pols, printf("Showing %d of %d formal relations.\n", n, #pols));
};

getrels(wt,deg) = formalsummary(formalrels(wt,deg));


\\ Change this accordingly 

wtlow = 10;
wtup = 15;
deglow = 2;
degup = 3;
for(i = wtlow, wtup, for( j = deglow, degup, getrels(i,j)));

