# Fixed-weight MZV relations

## AI Usage
Except for the Macaulay2 (M2) code, all content in this folder was generated using GPT-6.1-Sol via Codex.

## How to use

```sh
sage -python generate_relation_files.py
```

Generates weights 10–18 in the Hoffmann and depth-minimal bases. Completed
files are kept and saved conversion matrices are reused automatically.

```sh
sage -python generate_relation_files.py --weight 12
sage -python generate_relation_files.py --weight 16 --basis depth-minimal
sage -python generate_relation_files.py --min-weight 15 --max-weight 18 --force
```

`--force` regenerates relations using saved coordinates. `--output-dir path`
changes the output directory. A custom basis can be supplied with
`--weight N --basis-file my-basis.json`.

## Files

- `generate_relation_files.py`: runs the complete generation workflow.
- `basis_coordinates.py`: selects MZV bases and computes change-of-basis matrices.
- `high_weight_coordinates.py`: supplies coordinates at weights 16–18.
- `transform_relations.py`: substitutes coordinates and checks relation spans.
- `reduce_and_export_relations.py`: applies integer LLL and writes text/JSON.
- `product_basis_cache.py`: caches Sage's product-basis expansions in memory.

- `toric-source/`: the Macaulay2 script and its generated Lyndon monomial bases and toric relations.
- `mzv-linear-bases/`: saved MZV linear basis definitions and conversion matrices for both bases, weights 10–18.
- `hoffmann-basis/`, `depth-minimal-basis/`: final text and JSON relation files.

These linear bases and conversion matrices are computed locally using Sage
and, where needed, the published MZV datamine generators. The folder is not
a copy of the full datamine.

The generator runs Macaulay2 inside `toric-source/`, keeping its generated
`toric-relations-weight-*.json` files there.

Coordinate JSON contains `name`, `basis` (MZV index lists) and `coordinates`
(exact rational strings). Row i expresses source monomial i in the target
basis, following the source JSON's `mons` order. A custom basis file uses
this format; omit `coordinates` to compute them. To recompute saved
coordinates, remove the relevant coordinate file and rerun.

Relation JSON contains only `basis` and `relations`. Each relation has
`terms`, each with an integer string `coefficient` and an `exponents` vector
in basis order. Text files contain basis definitions and polynomial equations.
Indices follow Sage's increasing summation convention: the last index is at least 2.

## Calculation

The conversion matrix is substituted into the toric polynomials. Each
polynomial is made primitive with integer coefficients, then LLL with
`delta=0.99` reduces each polynomial degree separately. Zero rows are removed.
The unimodular transformation and exact relation span are checked. The integer
lattice is preserved without saturation. LLL is not canonical or guaranteed
to find globally shortest vectors; regenerated equations may differ while
spanning the same lattice.

Depth-minimal bases greedily exhaust each depth before proceeding to the next.
Sage computes the MZV coordinates. Weights 16–18 use factored period evaluations;
weight 18 uses the generators in
[the MZV datamine paper, Appendix B, equation (B.15)](https://arxiv.org/pdf/0907.2557).
Period reconstruction is checked at 1024 and 1536 bits, but is not an independent
proof of MZV identities. Subsequent polynomial and lattice arithmetic is exact.

Requires SageMath (tested with 10.9). Missing toric source files also require
Macaulay2 with FourTiTwo/4ti2; set `M2_EXECUTABLE` if needed. Supported weights
are 2–18. PARI's stack may grow to 8 GiB; high weights can be expensive.
