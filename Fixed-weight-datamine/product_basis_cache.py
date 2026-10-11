"""Cache Sage's expanded product bases within one calculation.

Sage's rho_inverse calls basis_data for every MZV it reduces. That iterator
re-expands the same products each time. The context below materializes each
fixed basis once and returns fresh iterators over it on subsequent calls.
No arithmetic, ordering, or precision settings are changed.
"""
from contextlib import contextmanager
from functools import wraps
from sage.modular.multiple_zeta import Multizetas


@contextmanager
def cached_product_bases():
    if getattr(Multizetas.basis_data, '_mzv_products_cached', False):
        yield
        return
    original = Multizetas.basis_data
    bases = {}

    def cached(self, basering, weight):
        key = (self, basering, weight)
        if key not in bases:
            bases[key] = tuple(original(self, basering, weight))
        return iter(bases[key])

    cached._mzv_products_cached = True
    Multizetas.basis_data = cached
    try:
        yield
    finally:
        Multizetas.basis_data = original


def with_cached_product_bases(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with cached_product_bases():
            return function(*args, **kwargs)
    return wrapped
