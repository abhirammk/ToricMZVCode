"""Extend Sage's motivic coordinates to weight 18 using a published basis.

Weight-18 generators: Blumlein, Broadhurst, Vermaseren, arXiv:0907.2557,
Appendix B, equation (B.15). Indices are reversed to Sage's convention.
The coaction is exact; the remaining period coefficient uses PARI rational
reconstruction, as in Sage, with an additional higher-precision check.
"""
from functools import lru_cache
from sage.all import QQ, Multizetas, Compositions, VectorSpace, matrix, pari, prod, vector
from sage.modular.multiple_zeta import (
    B_data, Values, compute_u_on_compo, extend_multiplicative_basis,
    phi_on_basis, phi_on_multiplicative_basis, rho_matrix_inverse,
)
from sage.modular.multiple_zeta_F_algebra import F_algebra

GENERATORS_18 = [
    (5,13), (3,15), (3,5,5,5), (5,3,3,7), (3,5,3,7),
    (3,3,3,9), (1,1,6,10), (1,1,3,3,4,6),
]


class Weight18Coordinates:
    def __init__(self, progress=None):
        if len(B_data) == 18:
            B_data.append(list(GENERATORS_18))
        elif len(B_data) < 19 or B_data[18] != GENERATORS_18:
            raise ValueError('unexpected Sage multiplicative basis data at weight 18')
        self.M = Multizetas(QQ)
        self.F = F_algebra(QQ)
        self.progress = progress
        P = matrix(QQ,[phi_on_basis(factors).homogeneous_to_vector()
                       for factors in extend_multiplicative_basis(B_data,18)])
        if P.nrows() != 65 or P.rank() != 65:
            raise ArithmeticError('published weight-18 generators do not give a full-rank coordinate system')
        self.inverse = P.inverse()
        if progress:
            progress('Published weight-18 multiplicative basis: rank 65 verified')

    @lru_cache(None)
    def period(self, c, precision):
        # Evaluate factors individually: avoid expanding rho_inverse into
        # thousands of MZVs and avoid allocating a weight-18 numerical cache.
        return pari.zetamult(list(reversed(c)), precision=precision)

    @lru_cache(None)
    def basis_periods(self, weight, precision):
        return tuple(prod((self.period(tuple(c),precision) for c in factors),pari(1))
                     for factors in extend_multiplicative_basis(B_data,weight))

    def period_coefficient(self, c, u, precision):
        weight = sum(c)
        rho_period = pari(0)
        if u:
            inv = self.inverse if weight == 18 else rho_matrix_inverse(weight)
            coefficients = u.homogeneous_to_vector()*inv
            rho_period = sum((cf*value for cf,value in zip(coefficients,self.basis_periods(weight,precision))),pari(0))
        ratio = (self.period(c,precision)-rho_period)/self.period((weight,),precision)
        return QQ(ratio.bestappr().sage())

    @lru_cache(None)
    def u(self, c):
        return compute_u_on_compo(tuple(c))

    @lru_cache(None)
    def phi(self, c):
        c = tuple(c)
        weight = sum(c)
        if weight < 16:
            return self.M(c).phi()
        u = self.u(c)
        if c in B_data[weight]:
            return phi_on_multiplicative_basis(c)
        a = self.period_coefficient(c,u,1024)
        b = self.period_coefficient(c,u,1536)
        if a != b:
            raise ArithmeticError(f'period rational reconstruction is unstable for {c}: {a} versus {b}')
        return u+a*self.F.gen(weight)

    def vector(self,c):
        return self.phi(tuple(c)).homogeneous_to_vector()

    def u_vector(self, c):
        u = self.u(tuple(c))
        return u.homogeneous_to_vector() if u else vector(QQ, len(self.M.basis_brown(sum(c))))

    def depth_basis(self,weight=18):
        dimension = len(self.M.basis_brown(weight))
        space = VectorSpace(QQ,dimension)
        span = space.subspace([])
        chosen = []
        # zeta(18) supplies the pure even direction. From now on, replacing
        # phi(c) by u(c) does not change independence or any depth cutoff.
        for depth in range(1,weight):
            tested = 0
            for c in Compositions(weight,length=depth):
                if c[0] < 2:
                    continue
                tested += 1
                c = tuple(reversed(c))
                v = self.F.gen(weight).homogeneous_to_vector() if depth == 1 else self.u_vector(c)
                if v not in span:
                    span = space.subspace(span.basis()+[v])
                    chosen.append(c)
                if len(chosen) == dimension:
                    if self.progress:
                        self.progress(f'Weight {weight}: depth {depth}, full rank {dimension} reached after {tested} candidates')
                    return chosen
            if self.progress:
                self.progress(f'Weight {weight}: depth {depth}, rank {len(chosen)}/{dimension}; {tested} candidates checked')
        raise ArithmeticError(f'depth selection failed to fill weight {weight}')
