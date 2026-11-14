import sys
import time
from math import isqrt
import numpy as np


def mobius_table(n):
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    prod = np.ones(n + 1, dtype=np.int64)
    idx = np.arange(n + 1, dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    for p in np.nonzero(sieve)[0]:
        p = int(p)
        mu[p::p] *= -1
        prod[p::p] *= p
        if p * p <= n:
            mu[p * p::p * p] = 0
    flip = (prod != idx) & (mu != 0)
    mu[flip] *= -1
    return mu


def squarefree_count(x, D=None):
    r = isqrt(x)
    if D is None:
        D = max(2, int(0.5 * x ** 0.4))
    D = min(D, r)
    mu = mobius_table(D)
    M = np.cumsum(mu, dtype=np.int64)
    d = np.arange(1, D + 1, dtype=object if x >= 2 ** 62 else np.int64)
    if x >= 2 ** 62:
        S1 = sum(int(mu[i]) * (x // (i * i)) for i in range(1, D + 1) if mu[i])
    else:
        S1 = int(np.sum(mu[1:].astype(np.int64) * (x // (d * d))))
    T = x // ((D + 1) * (D + 1))
    Mb = np.zeros(T + 2, dtype=np.int64)
    for t in range(T, 0, -1):
        y = isqrt(x // t)
        s = isqrt(y)
        kb = isqrt(T // t)
        total = 0
        if kb >= 2:
            ks = np.arange(2, kb + 1, dtype=np.int64)
            total += int(Mb[t * ks * ks].sum())
        if kb < s:
            ks = np.arange(kb + 1, s + 1, dtype=np.int64)
            total += int(M[y // ks].sum())
        sk = max(s, kb)
        vmax = y // (sk + 1)
        if vmax >= 1:
            vs = np.arange(1, vmax + 1, dtype=np.int64)
            nx = y // (vs + 1)
            hi = np.empty_like(nx)
            hi[0] = y
            hi[1:] = nx[:-1]
            lo = np.maximum(nx, sk)
            total += int((M[1:vmax + 1] * (hi - lo)).sum())
        Mb[t] = 1 - total
    S2 = int(Mb[1:T + 1].sum()) - T * int(M[D])
    return S1 + S2


def naive(x):
    r = isqrt(x)
    mu = mobius_table(r)
    d = np.arange(1, r + 1, dtype=np.int64)
    return int(np.sum(mu[1:].astype(np.int64) * (x // (d * d))))


if __name__ == "__main__":
    x = 10 ** int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 12
    t0 = time.time()
    q = squarefree_count(x)
    print(q, round(time.time() - t0, 3))
