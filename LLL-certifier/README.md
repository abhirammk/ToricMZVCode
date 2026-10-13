# Certified exclusion of bounded MZV relations

Create and independently verify certificates excluding nonzero integer
relations among MZV monomials. The default coefficient cap is `10^23 - 1`.

## Usage

Run from this folder. Verification requires Python 3.10 or newer and only
its standard library. Creation also requires PARI/GP (`gp` on PATH).

```sh
python3 certified.py verify examples/weight10-degree2-certificate.json
python3 certified.py create 10 2 --digits 23 --output weight10.json
python3 -m unittest discover -s tests -v
```

GP-dependent tests skip when GP is unavailable. Creation saves a certificate
only after exact verification succeeds, atomically replacing the output file.

Use `--digits h` for the coefficient cap `10^h - 1`, `--scale-digits N`
for the lattice scale `10^N`, and `--gp` or `GP_BIN` for the GP executable.
The default scale exponent is `n*(h+2)+80`, where `n` is the monomial count.
`GP_STACK` and `GP_STACK_MAX` control the initial GP stack and growth ceiling
(defaults: `256M` and `4G`).

An inconclusive result means the sufficient exclusion inequality failed;
it does not establish a relation. Increasing `--scale-digits` may help.
Exit status is 0 on success, 1 on errors or inconclusive results, and 130
on interruption.

## Monomials and proof

Generators are MZVs indexed by strict Lyndon words over `2 < 3`, with weight
at most `w`. Form all generator monomials of weight `w`, then products of
`d` such monomials, retaining each distinct exponent vector once. Degree
`d` gives total weight `d*w`; it does not specify the number of MZV factors.

The evaluator encloses exact values using integer arithmetic and Akhilesh's
double-tail recurrence. GP proposes a lattice transform. The verifier
regenerates the monomials and enclosures, checks unimodularity, and compares
exact Gram determinants against the bound for a hypothetical relation.
The certificate format and exact verification mathematics are unchanged.

See [the proof and references](docs/PROOF.md) and
[the mathematical PDF](docs/certificate-idea.pdf).

## Files

- `certified.py`: monomials, enclosures, exact verifier, creation, and CLI.
- `reducers.py`: PARI/GP lattice reduction and stack handling.
- `examples/`: a valid weight-10, degree-2 certificate.
- `tests/`: arithmetic, enclosure, serialization, and corruption checks.
- `docs/`: mathematical justification.
