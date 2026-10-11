"""MZV basis selection and reusable exact change-of-basis coordinates."""
import json
from functools import lru_cache
from pathlib import Path
from sage.all import QQ, Multizetas, Compositions, VectorSpace, matrix, prod
from sage.modular.multiple_zeta import Values
from product_basis_cache import with_cached_product_bases


@lru_cache(None)
def high_weight_backend():
    from high_weight_coordinates import Weight18Coordinates
    Values.update(15,1024)
    return Weight18Coordinates()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def exponents(v, n):
    if not isinstance(v, list) or len(v) != n or any(type(e) is not int or e < 0 for e in v):
        raise ValueError(f'expected {n} nonnegative integer exponents')
    return v


def word(w):
    if not isinstance(w, list) or not w or any(type(i) is not int or i < 1 for i in w) or w[-1] < 2:
        raise ValueError('MZV indices must be positive, with last index >= 2 (Sage convention)')
    return tuple(w)


@with_cached_product_bases
def depth_minimal_basis(M, weight, progress=None):
    """Sage's depth-filtration greedy selection, stopping once a basis is full.

    Enumerate in the same order as M.basis_filtration(weight). Exhaust each
    smaller depth before moving on; after reaching the full weight dimension,
    all remaining candidates are necessarily dependent and can be skipped.
    """
    Values.update(weight, 1024)
    dimension = len(M.basis_brown(weight))
    space = VectorSpace(QQ, dimension)
    span = space.subspace([])
    chosen = []
    for depth in range(1, weight):
        tested = 0
        for c in Compositions(weight, length=depth):
            if c[0] < 2:
                continue
            tested += 1
            b = M(tuple(reversed(c)))
            v = b.phi_as_vector()
            if v not in span:
                span = space.subspace(span.basis() + [v])
                chosen.append(b)
            if len(chosen) == dimension:
                if progress:
                    progress(f'Weight {weight}: depth {depth}, full rank {dimension} reached after {tested} candidates')
                return chosen
        if progress:
            progress(f'Weight {weight}: depth {depth}, rank {len(chosen)}/{dimension}; {tested} candidates checked')
    raise ArithmeticError('depth-filtration search did not find a full basis')


@with_cached_product_bases
def coordinate_system(data, name='hoffmann', custom=None, progress=None):
    """Matrix rows express source monomials in the columns' target basis."""
    if custom is None and name not in ('hoffmann','depth-minimal'):
        raise ValueError(f'unknown basis: {name}')
    words = [word(w) for w in data['vars']]
    mons = [exponents(m, len(words)) for m in data['mons']]
    weights = {sum(sum(w)*e for w, e in zip(words, m)) for m in mons}
    if len(weights) != 1 or next(iter(weights)) <= 0:
        raise ValueError('source monomials must have one positive weight')
    weight = weights.pop()
    for rel in data['binoms']:
        a, b = [exponents(rel[k], len(mons)) for k in ('mon1', 'mon2')]
        if any(sum((x-y)*m[j] for x, y, m in zip(a, b, mons)) for j in range(len(words))):
            raise ValueError('input binomial is not a toric relation')
    M = Multizetas(QQ)
    if custom:
        indices = [word(w) for w in custom['basis']]
    elif name == 'hoffmann':
        indices = [tuple(b.support()[0]) for b in M.basis_brown(weight)]
    elif name == 'depth-minimal':
        indices = (high_weight_backend().depth_basis(weight) if weight >= 16 else
                   [tuple(b.support()[0]) for b in depth_minimal_basis(M,weight,progress)])
    else:
        raise ValueError(f'unknown basis: {name}')
    if not indices or any(sum(w) != weight for w in indices):
        raise ValueError('target basis must consist of MZVs of the source weight')
    descriptors = [{'indices':list(map(int,w))} for w in indices]
    if custom and 'coordinates' in custom:
        rows = custom['coordinates']
        if len(rows) != len(mons) or any(len(r) != len(indices) for r in rows):
            raise ValueError('coordinates must have one row per source monomial and one column per target MZV')
        if any(type(c) not in (int,str) for r in rows for c in r):
            raise ValueError('coordinates must be integers or rational strings')
        C = matrix(QQ,rows)
    else:
        phi = high_weight_backend().phi if weight >= 16 else lambda w: M(w).phi()
        P = matrix(QQ,[phi(w).homogeneous_to_vector() for w in indices])
        if P.nrows() != P.ncols() or P.rank() != len(indices):
            raise ValueError('target MZVs must form a basis of the full weight component')
        unit = phi(words[0]).parent().one()
        generators = [phi(w) if any(m[i] for m in mons) else unit for i,w in enumerate(words)]
        source = matrix(QQ,[prod((g**e for g,e in zip(generators,m)),unit).homogeneous_to_vector() for m in mons])
        C = source*P.inverse()
    if C.nrows() != C.ncols() or C.rank() != C.nrows():
        raise ValueError('change of basis must be square and invertible')
    return weight, descriptors, C


def compute_coordinates(data, name='hoffmann', custom=None, progress=None):
    """Return the uniform name/basis/coordinates file format."""
    weight, descriptors, C = coordinate_system(data,name,custom,progress)
    return {'name':name,'basis':[d['indices'] for d in descriptors],
            'coordinates':[[str(c) for c in row] for row in C.rows()]}
