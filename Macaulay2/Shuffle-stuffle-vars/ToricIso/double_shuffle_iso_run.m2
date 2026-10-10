load "double_shuffle_iso1.m2";

for i from 11 to 11 do (
wt = i;
(zetaI,autM) = compute(wt);
"./results/ideals/zetaIdealWt" | toString(wt) | ".txt" << toString(zetaI) << endl << close;
"./results/matrices/zetaMatrixWt" | toString(wt) | ".txt" << toString(autM) << endl << close;

print("Dimension in weight " | toString(wt) | ": " | toString(dim zetaI));
print("Generators in weight " | toString(wt) | ":");
print(betti mingens zetaI);

"./results/log/zetaLogWt" | toString(wt) | ".txt" << ("Dimension in weight " | toString(wt) | ": " | toString(dim zetaI) | "\n") << ("Generators in weight " | toString(wt) | ":") << betti mingens zetaI << endl << close;

);
