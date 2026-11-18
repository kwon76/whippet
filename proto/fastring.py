"""Exact arithmetic in R_q for large operands through float64 matrix products (BLAS).

The routines compute the same integers as the reference routines of ring.py and tail.py:
  matvec_small_batch(A, S)   == np.stack([ring.matvec_small(A, S[i]) for i in range(len(S))])
  lincomb_small(c, S)        == tail.mul_int(c[:, None, :], S).sum(axis=0)
  lincomb_uniform(al, A)     == ring.mul(al[:, None, :], A).sum(axis=0) % q
  scalar_dot(a, S)           == sum_j a_j * S[:, j, :] mod q
They are exact because every float64 product is a sum of integers below 2^53: operands with entries
in [0, q) are split into two 16-bit limbs, and small operands are checked against an explicit bound
before use (otherwise the reference routine is called).  test() compares them with the reference
routines on random inputs; run_tests.py calls it.

Setting the environment variable WHIPPET_SLOW=1 makes every routine call the reference routine.
"""
import os
import numpy as np
import ring as R

Q = R.Q
QU = R.QU
D = R.D
SLOW = os.environ.get("WHIPPET_SLOW") == "1"
EXACT = float(2 ** 53)
MEM = 1 << 28                      # bytes of float64 temporaries per chunk (about 256 MB)
_R32 = (1 << 32) % Q               # 2^32 mod q = 99

_K = np.arange(D)[:, None]
_L = np.arange(D)[None, :]
_U = (_K - _L) % D                 # M(c)[k, l] = c[(k - l) mod 64] * (-1 if l > k)
_S = np.where(_L > _K, -1, 1)


def negacyclic(c):
    """c: (..., 64) integers.  Returns M (..., 64, 64) with (c * s)[k] = sum_l M[k, l] s[l]."""
    c = np.asarray(c, dtype=np.int64)
    return c[..., _U] * _S


def _limbs(a):
    a = np.asarray(a, dtype=np.int64)
    return (a >> 16).astype(np.float64), (a & 0xFFFF).astype(np.float64)


def _absmax(s):
    """an upper bound on the absolute values of the entries (the type bound for 8- and 16-bit arrays,
    which avoids a full pass over large digit arrays)"""
    s = np.asarray(s)
    if s.dtype in (np.int8, np.int16):
        return int(np.iinfo(s.dtype).max) + 1
    return int(np.abs(s).max()) if s.size else 0


# ---------------------------------------------------------------- A s for small s

class Prepared:
    """A (k, n, 64) with entries in [0, q), prepared for repeated products with small vectors."""

    def __init__(self, A):
        self.A = A
        self.k, self.n = A.shape[0], A.shape[1]
        if not SLOW:
            A2 = np.ascontiguousarray(np.asarray(A, dtype=np.int64).transpose(0, 2, 1).reshape(self.k * D, self.n))
            self.hi, self.lo = _limbs(A2)
            del A2

    def apply(self, S):
        """S: (r, n, 64) small integers.  Returns (r, k, 64) uint64, row i = A S[i]."""
        S = np.asarray(S)
        r = S.shape[0]
        assert S.shape[1] == self.n
        if SLOW or _absmax(S) * self.n * 65536.0 >= EXACT:
            return np.stack([R.matvec_small(self.A, S[i].astype(np.int64)) for i in range(r)])
        out = np.empty((r, self.k, D), dtype=np.uint64)
        rc = max(1, min(r, MEM // (8 * self.n)))                       # one coefficient of rc blocks at a time
        for s0 in range(0, r, rc):
            blk = S[s0:s0 + rc]                                         # (rc, n, 64), any integer type
            nr = blk.shape[0]
            acc = np.zeros((nr, self.k, D), dtype=np.int64)
            for u in range(D):
                bu = blk[:, :, u]
                if not bu.any():
                    continue
                su = np.ascontiguousarray(bu.T, dtype=np.float64)       # (n, rc)
                ph = self.hi @ su
                pl = self.lo @ su
                part = ((ph.astype(np.int64) % Q) * 65536 + pl.astype(np.int64)) % Q     # (k*64, rc)
                part = part.reshape(self.k, D, nr).transpose(2, 0, 1)                    # (rc, k, 64)
                rr = part[:, :, R.RIDX[u]]
                if u:
                    rr[:, :, :u] = -rr[:, :, :u]
                acc += rr
            out[s0:s0 + nr] = (acc % Q).astype(np.uint64)
        return out


def matvec_small_batch(A, S):
    return Prepared(A).apply(S)


# ---------------------------------------------------------------- sum_i c_i S_i for small c and S

def lincomb_small(c, S):
    """c: (r, 64) small integers (challenges), S: (r, n, 64) small integers.  Exact sum_i c_i S_i, int64."""
    c = np.asarray(c, dtype=np.int64)
    S = np.asarray(S)
    r, n = S.shape[0], S.shape[1]
    bound = _absmax(c) * _absmax(S) * r * D
    if SLOW or bound >= EXACT:
        acc = np.zeros((n, D), dtype=np.int64)
        for i in range(r):
            acc += _mul_int(c[i], S[i])
        return acc
    M = negacyclic(c).transpose(0, 2, 1).astype(np.float64)            # (r, 64[l], 64[k]) = M(c_i)^T
    out = np.zeros((n, D), dtype=np.float64)
    rc = max(1, min(r, MEM // (8 * n * D)))
    for s0 in range(0, r, rc):
        blk = S[s0:s0 + rc].astype(np.float64)                          # (rc, n, 64)
        nr = blk.shape[0]
        out += blk.transpose(1, 0, 2).reshape(n, nr * D) @ M[s0:s0 + nr].reshape(nr * D, D)
    return np.rint(out).astype(np.int64)


def _mul_int(a, b):
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    acc = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.int64)
    for u in range(D):
        au = a[..., u:u + 1]
        if not au.any():
            continue
        rr = b[..., R.RIDX[u]]
        if u:
            rr = rr.copy()
            rr[..., :u] = -rr[..., :u]
        acc += au * rr
    return acc


# ---------------------------------------------------------------- sum_k alpha_k A_k for uniform alpha and A

def lincomb_uniform(alpha, A):
    """alpha: (k, 64), A: (k, n, 64), entries in [0, q).  Returns sum_k alpha_k A_k mod q, (n, 64) uint64."""
    alpha = np.asarray(alpha, dtype=np.uint64)
    A = np.asarray(A, dtype=np.uint64)
    k, n = A.shape[0], A.shape[1]
    if SLOW or k * D * float(1 << 32) >= EXACT:
        return R.mul(alpha[:, None, :], A).sum(axis=0) % QU
    M = negacyclic(alpha.astype(np.int64))                              # (k, 64, 64), entries in (-q, q)
    M = np.where(M < 0, M + Q, M)                                       # in [0, q)
    Mt = M.transpose(0, 2, 1).reshape(k * D, D)                         # rows (k, l), columns out index
    mh, ml = _limbs(Mt)
    out = np.empty((n, D), dtype=np.uint64)
    nc = max(1, min(n, MEM // (8 * k * D * 3)))
    for j0 in range(0, n, nc):
        blk = np.asarray(A[:, j0:j0 + nc], dtype=np.int64).transpose(1, 0, 2).reshape(-1, k * D)
        ah, al = _limbs(blk)
        hh = (ah @ mh).astype(np.int64) % Q
        mid = ((ah @ ml).astype(np.int64) + (al @ mh).astype(np.int64)) % Q
        ll = (al @ ml).astype(np.int64) % Q
        out[j0:j0 + blk.shape[0]] = ((hh * _R32 + mid * 65536 + ll) % Q).astype(np.uint64)
    return out


# ---------------------------------------------------------------- sum_j a_j S_j for scalars a_j in Z_q

def scalar_dot(a, S):
    """a: (n,) entries in [0, q); S: (r, n, 64) small integers.  Returns (r, 64) uint64."""
    a = np.asarray(a, dtype=np.uint64)
    S = np.asarray(S)
    r, n = S.shape[0], S.shape[1]
    if SLOW or _absmax(S) * n * 65536.0 >= EXACT:
        out = np.zeros((r, D), dtype=np.uint64)
        aq = a.astype(object)
        for i in range(r):
            out[i] = ((S[i].astype(object) * aq[:, None]).sum(axis=0) % Q).astype(np.uint64)
        return out
    ah, al = _limbs(a)
    out = np.empty((r, D), dtype=np.uint64)
    rc = max(1, min(r, MEM // (8 * n * D)))
    for s0 in range(0, r, rc):
        blk = S[s0:s0 + rc].astype(np.float64).transpose(0, 2, 1)       # (rc, 64, n)
        ph = (blk @ ah).astype(np.int64) % Q
        pl = (blk @ al).astype(np.int64)
        out[s0:s0 + blk.shape[0]] = ((ph * 65536 + pl) % Q).astype(np.uint64)
    return out


# ---------------------------------------------------------------- self-test against the reference routines

def test(trials=40, seed=11, log=print):
    rng = np.random.default_rng(seed)
    bad = 0
    for t in range(trials):
        k = int(rng.integers(1, 9))
        n = int(rng.integers(1, 300))
        r = int(rng.integers(1, 12))
        lim = int(rng.choice([1, 2, 9, 33, 129, 1 << 12, 1 << 19]))
        A = R.uniform(rng, (k, n, D))
        S = rng.integers(-lim, lim + 1, size=(r, n, D), dtype=np.int64)
        if t % 5 == 0:
            S[:, :, 1:] = 0                                              # sparse and structured inputs
        ref = np.stack([R.matvec_small(A, S[i]) for i in range(r)])
        bad += int(not (Prepared(A).apply(S) == ref).all())
        bad += int(not (Prepared(A).apply(S.astype(np.int16) if lim < 32000 else S) == ref).all())
        c = R.challenges(rng, r)
        ref = np.zeros((n, D), dtype=np.int64)
        for i in range(r):
            ref += _mul_int(c[i], S[i])
        bad += int(not (lincomb_small(c, S) == ref).all())
        al = R.uniform(rng, (k, D))
        bad += int(not (lincomb_uniform(al, A) == R.mul(al[:, None, :], A).sum(axis=0) % QU).all())
        a = R.uniform(rng, (n,))
        ref = np.stack([((S[i].astype(object) * a.astype(object)[:, None]).sum(axis=0) % Q).astype(np.uint64)
                        for i in range(r)])
        bad += int(not (scalar_dot(a, S) == ref).all())
    log("fastring.test: %d trials, %d mismatches" % (trials, bad))
    return bad == 0


if __name__ == "__main__":
    import sys
    sys.exit(0 if test(int(sys.argv[1]) if len(sys.argv) > 1 else 40) else 1)
