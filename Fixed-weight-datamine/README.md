# Toric relations in fixed-weight MZV bases

This code computes generators of the toric ideal associated with the Lyndon
monomial basis in a fixed weight, and expresses them in the Hoffmann or a
depth-minimal MZV basis. The resulting polynomial relations are reduced by
LLL and written to text and JSON files.

## Usage

Run with SageMath (tested with 10.9):

```sh
sage -python generate_relation_files.py
```

The default is weights 10–18, in both bases. Existing output files and
change-of-basis matrices are reused. To select a weight, basis or range:

```sh
sage -python generate_relation_files.py --weight 12
sage -python generate_relation_files.py --weight 16 --basis depth-minimal
sage -python generate_relation_files.py --min-weight 15 --max-weight 18
```

Add `--force` to recompute the relations using the saved matrices, or
`--output-dir path` to write the results elsewhere. To recompute a matrix,
remove its saved file; if using another output directory, also remove the
corresponding project copy, which is used as a fallback.

Macaulay2 with FourTiTwo/4ti2 is needed only when the toric source file is
missing. Set `M2_EXECUTABLE` if Macaulay2 is not on PATH. The driver accepts
weights 2–18. High-weight calculations can require substantial time and
memory; PARI's stack is allowed to grow to 8 GiB.

## Calculation

The source monomials are products of the MZVs indexed by Lyndon words on
2 and 3. Macaulay2 first treats the Lyndon generators as formal independent
variables and computes binomial generators of the monomial toric ideal.
The change-of-basis matrix expresses each weight-w source monomial as a
linear combination of the chosen weight-w MZVs. Substituting these linear
forms into the binomials gives the relations in the new basis. A polynomial
of degree d in the basis variables has MZV weight d times w.

The Hoffmann basis is Sage's `basis_brown`, with index words in 2 and 3.
For the depth-minimal basis, independent MZVs are selected in increasing
depth, exhausting each depth before proceeding. This gives as many
low-depth elements as possible; choices within a depth depend on the
enumeration order. Indices use Sage's increasing summation convention:
the last index is at least 2. Depth is the number of indices.

After substitution, denominators are cleared and each polynomial is divided
by its integer coefficient gcd. LLL is applied separately in each polynomial
degree, with `delta=0.99`. The nonzero output rows form an integer lattice
basis of these primitive input relations. No saturation or further gcd
division is performed. The output generates the same ideal over the rationals
and is linearly independent degree by degree. It need not be a minimal set
of ideal generators. Different runs of LLL can give different reduced bases.

Basis conversion uses Sage's map to the auxiliary F-algebra, including
numerical rational reconstruction of some coefficients. The subsequent
matrix and polynomial calculations are exact over the rationals. The code
checks invertibility, the unimodular LLL transformation and preservation
of the relation spans.

## Code

- `generate_relation_files.py`: the command-line driver; obtains the source relations and matrices, then runs conversion, reduction and export.
- `toric-source/compute_toric_relations.m2`: enumerates Lyndon words and fixed-weight monomials, and computes their toric ideal generators with FourTiTwo.
- `basis_coordinates.py`: selects the MZV basis and computes the change-of-basis matrix from the F-algebra images of the source products and target MZVs.
- `high_weight_coordinates.py`: computes the F-algebra images needed at weights 16–18. It uses Sage's coaction calculation and determines the remaining zeta(w) coefficient by period reconstruction, evaluating products factor by factor. Reconstructions must agree at 1024 and 1536 bits; this is a consistency check. At weight 18 it adds the multiplicative generators from [the MZV datamine paper, Appendix B, equation (B.15)](https://arxiv.org/pdf/0907.2557).
- `transform_relations.py`: substitutes the change-of-basis matrix into the toric generators and checks the rational relation spans after reduction.
- `reduce_and_export_relations.py`: normalizes the integer coefficients, applies LLL, checks the integer transformation, and writes the relation files.
- `product_basis_cache.py`: keeps Sage's repeatedly used expansions of multiplicative basis products in memory during a calculation.

## Files

- `toric-source/toric-relations-weight-<w>.json`: Lyndon words (`vars`), source monomial exponents (`mons`), and binomial exponents (`binoms`).
- `mzv-linear-bases/change-of-basis-weight-<w>-<basis>.json`: target MZV index lists and the rational change-of-basis matrix. These are local calculations using Sage and, at weight 18, the published datamine generators.
- `<basis>-basis/relations-weight-<w>-<basis>.txt`: definitions of the basis variables followed by the integer polynomial equations.
- `<basis>-basis/relations-weight-<w>-<basis>.json`: the same basis definitions and equations, encoded as `basis` and `relations`. Each relation contains `terms`, with integer strings `coefficient` and exponent vectors `exponents` in basis order.

## Other bases

```sh
sage -python generate_relation_files.py --weight 12 --basis-file my-basis.json
```

A basis file contains `name` and `basis` (a list of MZV index lists), with an
optional `coordinates` matrix. It has the same format as the saved
change-of-basis files. Row i of the matrix expresses source monomial i,
in the source `mons` order, in the listed basis. Entries must be integers
or rational strings, not floats.

Without a matrix, the code checks independence in the F-algebra and computes
the coordinates. A supplied matrix is checked for dimensions and
invertibility; its MZV identities are assumed correct.

## AI usage

Except for the Macaulay2 (M2) file, all the code in this folder was generated using GPT-6.1-Sol via Codex.
