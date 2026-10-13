#!/usr/bin/env python3
"""Exact bounded MZV relation certificates. Verification uses only the stdlib.

MZV enclosures: Akhilesh, arXiv:2105.12156, Theorems 3 and 5.
All endpoints are integers divided by 2**bits. No floating-point operations.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

# Lattice coordinates and certificates routinely exceed Python's default
# 4300-digit decimal conversion guard. Allow arbitrary-size exact integers
# for GP input/output and JSON, including when imported by other programs.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

FORMAT = "mzv-bounded-exclusion-v1"


class ExclusionInconclusive(ValueError):
    """The sufficient lattice exclusion test did not pass."""


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def ceildiv(a, b):
    return -(-a // b)


def words_at_weight(weight, prefix=()):
    if weight == 0:
        yield prefix
    for letter in (2, 3):
        if weight >= letter:
            yield from words_at_weight(weight-letter, prefix+(letter,))


def monomial_data(weight, degree):
    integer(weight, "weight")
    integer(degree, "degree")
    words = [w for k in range(1, weight+1) for w in words_at_weight(k)
             if all(w < w[j:]+w[:j] for j in range(1, len(w)))]
    weights = [sum(w) for w in words]

    def weighted(pos, remaining, prefix):
        if pos == len(words):
            if remaining == 0:
                yield prefix
            return
        for exponent in range(remaining//weights[pos]+1):
            yield from weighted(pos+1, remaining-exponent*weights[pos],
                                prefix+(exponent,))

    base = list(weighted(0, weight, ()))
    support = {tuple(0 for _ in words)}
    for _ in range(degree):
        support = {tuple(a+b for a, b in zip(u, v))
                   for u in support for v in base}
    return words, sorted(support)


def binary_word(composition):
    if not composition or composition[0] < 2 or any(s < 1 for s in composition):
        raise ValueError("MZV composition must be nonempty and admissible")
    return ''.join('0'*(s-1)+'1' for s in composition)


def split_word(w):
    """w = 0 1**(b-1) mid 0**(a-1) 1 (Akhilesh Theorem 5)."""
    if len(w) < 2 or w[0] != '0' or w[-1] != '1':
        raise ValueError("expected an admissible binary word")
    left = 1
    while left < len(w)-1 and w[left] == '1':
        left += 1
    right = len(w)-1
    while right > left and w[right-1] == '0':
        right -= 1
    a, b = len(w)-right, left
    return a, b, w[:right], w[left:], w[left:right]


def mzv_enclosures(compositions, bits, steps):
    """Return positive MZV interval endpoints in units of 2**(-bits).

    At n=N every admissible double tail lies in [0, 2*4**(-N)].
    Empty/one-letter tails are exact 1/binomial(2n,n) and that quantity/n.
    Descend with Theorem 5; all coefficients are positive. Integer division
    rounds each lower endpoint down and each upper endpoint up.
    """
    integer(bits, "bits", 1)
    integer(steps, "steps", 1)
    requested = [binary_word(c) for c in compositions]
    rules = {}

    def include(w):
        if len(w) < 2 or w in rules:
            return
        rule = split_word(w)
        rules[w] = rule
        for subword in rule[2:]:
            include(subword)

    for w in requested:
        include(w)
    D = 1 << bits
    initial_upper = ceildiv(2*D, 1 << (2*steps))
    values = {w: (0, initial_upper) for w in rules}
    central_binomial = math.comb(2*steps, steps)
    for n in range(steps, 0, -1):
        values[''] = (D//central_binomial, ceildiv(D, central_binomial))
        den = n*central_binomial
        values['0'] = values['1'] = (D//den, ceildiv(D, den))
        next_values = {}
        for w, (a, b, init, fin, mid) in rules.items():
            lo, hi = values[w]
            for subword, power in ((init, a), (fin, b), (mid, a+b)):
                denominator = n**power
                lower, upper = values[subword]
                lo += lower//denominator
                hi += ceildiv(upper, denominator)
            next_values[w] = (lo, hi)
        values = next_values
        # binomial(2(n-1), n-1) = binomial(2n,n)*n / (2*(2n-1)).
        central_binomial = central_binomial*n//(2*(2*n-1))
    return [values[w] for w in requested]


def monomial_enclosures(words, support, bits, steps):
    D = 1 << bits
    intervals = mzv_enclosures(words, bits, steps)
    maxexp = [max((e[j] for e in support), default=0) for j in range(len(words))]
    powers = []
    for (lo, hi), count in zip(intervals, maxexp):
        table = [(D, D)]
        for _ in range(count):
            lower, upper = table[-1]
            table.append((lower*lo//D, ceildiv(upper*hi, D)))
        powers.append(table)
    result = []
    for exponents in support:
        lo = hi = D
        for j, exponent in enumerate(exponents):
            lower, upper = powers[j][exponent]
            lo, hi = lo*lower//D, ceildiv(hi*upper, D)
        result.append((lo, hi))
    return result


def determinant(matrix):
    """Bareiss determinant with exact division and row pivoting."""
    A = [row[:] for row in matrix]
    size = len(A)
    previous, sign = 1, 1
    for k in range(size-1):
        pivot_row = next((i for i in range(k, size) if A[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            A[k], A[pivot_row] = A[pivot_row], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, size):
            for j in range(k+1, size):
                numerator = pivot*A[i][j]-A[i][k]*A[k][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("inexact Bareiss division")
                A[i][j] = quotient
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1] if size else 1


def gram_pivots(columns):
    """Leading Gram determinants; successive ratios are squared GS norms."""
    n = len(columns)
    G = [[sum(x*y for x, y in zip(u, v)) for v in columns] for u in columns]
    previous = 1
    pivots = [1]
    for k in range(n):
        pivot = G[k][k]
        if pivot <= 0:
            raise ValueError("basis is not independent")
        pivots.append(pivot)
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*G[i][j]-G[i][k]*G[k][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("inexact Gram Bareiss division")
                G[i][j] = G[j][i] = quotient
        previous = pivot
    return pivots


def exclusion_check(U, rounded, errors, D, multiplier, cap):
    n = len(rounded)
    if len(U) != n or any(type(row) is not list or len(row) != n for row in U):
        raise ValueError("wrong transformation dimensions")
    if any(type(entry) is not int for row in U for entry in row):
        raise ValueError("transformation entries must be integers")
    if abs(determinant(U)) != 1:
        raise ValueError("transformation is not unimodular")
    columns = [[multiplier*U[i][j] for i in range(n)] +
               [sum(rounded[i]*U[i][j] for i in range(n))] for j in range(n)]
    pivots = gram_pivots(columns)
    # Any relation in the box gives length^2 <= cap^2*(K^2*n+(sum errors/D)^2).
    rhs_numerator = cap*cap*(multiplier*multiplier*n*D*D+sum(errors)**2)
    if not all(pivots[i+1]*D*D > pivots[i]*rhs_numerator for i in range(n)):
        raise ExclusionInconclusive("inconclusive: exact lattice bound does not exclude the coefficient box")
    return pivots


def verify(cert):
    if type(cert) is not dict or cert.get('format') != FORMAT:
        raise ValueError("unsupported certificate format")
    weight = integer(cert['weight'], 'weight')
    degree = integer(cert['degree'], 'degree')
    digits = integer(cert['max_digits'], 'max_digits', 1)
    bits = integer(cert['bits'], 'bits', 1)
    steps = integer(cert['steps'], 'steps', 1)
    scale = integer(cert['scale'], 'scale', 1)
    multiplier = integer(cert['multiplier'], 'multiplier', 1)
    words, support = monomial_data(weight, degree)
    if not support:
        raise ValueError("empty monomial basis")
    if cert['words'] != [list(w) for w in words] or cert['exponents'] != [list(e) for e in support]:
        raise ValueError("basis metadata does not match independently regenerated supports")
    rounded = cert['rounded']
    if type(rounded) is not list or len(rounded) != len(support) or any(type(r) is not int for r in rounded):
        raise ValueError("invalid rounded coordinates")
    intervals = monomial_enclosures(words, support, bits, steps)
    D = 1 << bits
    errors = [max(abs(r*D-scale*lo), abs(r*D-scale*hi))
              for r, (lo, hi) in zip(rounded, intervals)]
    exclusion_check(cert['transform'], rounded, errors, D, multiplier, 10**digits-1)
    return len(support)


def create(weight, degree, digits=23, scale_digits=None, gp="gp"):
    """Propose and verify a certificate; never return an unchecked proposal."""
    from reducers import transform

    integer(digits, "digits", 1)
    def announce(message):
        print(message, file=sys.stderr)

    announce(f"Generating monomials for weight {weight}, degree {degree}...")
    words, support = monomial_data(weight, degree)
    n = len(support)
    if not n:
        raise ValueError("empty monomial basis")
    if scale_digits is None:
        scale_digits = n*(digits+2)+80
    integer(scale_digits, 'scale_digits', 1)
    scale = 10**scale_digits
    bits = (10**(scale_digits+50)).bit_length()
    steps = (bits+1)//2+32
    multiplier = 10*n
    announce(f"Enclosing {len(words)} MZVs for {n} monomials ({bits} bits, {steps} steps)...")
    intervals = monomial_enclosures(words, support, bits, steps)
    D = 1 << bits
    rounded = [(scale*(lo+hi)+D)//(2*D) for lo, hi in intervals]
    announce("Reducing the integer lattice with GP...")
    U = transform(rounded, multiplier, gp)
    cert = dict(format=FORMAT, weight=weight, degree=degree, max_digits=digits,
                bits=bits, steps=steps, scale=scale, multiplier=multiplier,
                words=[list(w) for w in words], exponents=[list(e) for e in support],
                rounded=rounded, transform=U)
    announce("Checking certificate using exact arithmetic...")
    verify(cert)
    return cert


def write_json(path, value):
    """Replace the output after writing the complete certificate."""
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    producer = commands.add_parser('create', help='create a certificate only if exact verification succeeds')
    producer.add_argument('weight', type=int)
    producer.add_argument('degree', type=int)
    producer.add_argument('--digits', type=int, default=23)
    producer.add_argument('--scale-digits', type=int)
    producer.add_argument('--gp', default=os.environ.get('GP_BIN', 'gp'))
    producer.add_argument('--output', type=Path, required=True)
    checker = commands.add_parser('verify', help='verify without GP or numerical libraries')
    checker.add_argument('certificate', type=Path)
    args = parser.parse_args()
    try:
        if args.command == 'create':
            integer(args.digits, 'digits', 1)
            cert = create(args.weight, args.degree, args.digits, args.scale_digits, args.gp)
            write_json(args.output, cert)
            print(f"PROVED: no nonzero relation with |a_i| < 10^{args.digits}; certificate: {args.output}")
        else:
            cert = json.loads(args.certificate.read_text())
            n = verify(cert)
            print(f"VERIFIED: weight={cert['weight']}, degree={cert['degree']}, dimension={n}, |a_i| <= {10**cert['max_digits']-1}")
    except ExclusionInconclusive as error:
        print(f"INCONCLUSIVE: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("INTERRUPTED: no new certificate saved", file=sys.stderr)
        return 130
    except (ValueError, KeyError, TypeError, OSError, ArithmeticError, subprocess.SubprocessError) as error:
        print(f"NO CERTIFICATE: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
