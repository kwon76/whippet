"""Arithmetic for the reference implementation.

R_q = Z_q[X]/(X^64+1), q = 2^32 - 99 (q = 5 mod 8).  Elements are numpy uint64 arrays (..., 64) with
entries in [0, q).  The field F = F_{q^4} is the image of S+ = Z[eta], eta = X^8 + X^-8, inside R_q;
an element of F is stored by its four coordinates in the basis (1, eta_1, eta_2, eta_3) with
eta_k = X^{8k} + X^{-8k} = X^{8k} - X^{64-8k}.
"""
import hashlib
import numpy as np

D = 64
Q = 2**32 - 99
QU = np.uint64(Q)
HALF = (Q + 1) // 2                      # inverse of 2 mod q

_IDX = np.arange(D)


def _rot_tables():
    idx = np.zeros((D, D), dtype=np.int64)
    neg = np.zeros((D, D), dtype=bool)
    for u in range(D):
        for c in range(D):
            j = c - u
            if j < 0:
                j += D
                neg[u, c] = True
            idx[u, c] = j
    return idx, neg


RIDX, RNEG = _rot_tables()


def rot(b, u):
    """X^u * b (negacyclic shift), entries stay in [0, q)."""
    r = b[..., RIDX[u]]
    if u:
        r = r.copy()
        r[..., :u] = (QU - r[..., :u]) % QU
    return r


def mul(a, b):
    """Product in R_q of arrays broadcastable against each other (both with entries in [0,q))."""
    a = np.asarray(a, dtype=np.uint64)
    b = np.asarray(b, dtype=np.uint64)
    shape = np.broadcast_shapes(a.shape, b.shape)
    acc = np.zeros(shape, dtype=np.uint64)
    for u in range(D):
        au = a[..., u:u + 1]
        if not au.any():
            continue
        acc += (au * rot(b, u)) % QU
    return acc % QU


def mul_small(a, s):
    """a in R_q (uint64), s with small signed coefficients (int64, |s| < 2^20); a.shape broadcastable.
    Returns a*s mod q.  Sums are kept in int64 without intermediate reduction."""
    a = np.asarray(a, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    shape = np.broadcast_shapes(a.shape, s.shape)
    acc = np.zeros(shape, dtype=np.int64)
    for u in range(D):
        su = s[..., u:u + 1]
        if not su.any():
            continue
        r = a[..., RIDX[u]]
        if u:
            r = r.copy()
            r[..., :u] = -r[..., :u]
        acc += su * r
        if u % 16 == 15:
            acc %= Q
    return (acc % Q).astype(np.uint64)


def matvec_small(A, s):
    """A: (k, n, 64) uint64, s: (n, 64) small int64.  Returns A s in R_q^k."""
    A = np.asarray(A, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    k = A.shape[0]
    acc = np.zeros((k, D), dtype=np.int64)
    for u in range(D):
        su = s[:, u]
        if not su.any():
            continue
        part = np.tensordot(A, su, axes=([1], [0]))          # (k, 64): sum_j A[k,j,:] * s[j,u]
        part %= Q
        r = part[:, RIDX[u]]
        if u:
            r = r.copy()
            r[:, :u] = -r[:, :u]
        acc += r
    return (acc % Q).astype(np.uint64)


def dot(a, b):
    """<a, b> for a, b of shape (n, 64) in R_q."""
    return mul(a, b).sum(axis=0) % QU


def center(a):
    """centred representatives in (-q/2, q/2] as int64"""
    a = np.asarray(a).astype(np.int64)
    return np.where(a > Q // 2, a - Q, a)


def to_q(s):
    return (np.asarray(s, dtype=np.int64) % Q).astype(np.uint64)


# ---------------------------------------------------------------- the field F inside R_q

def emb(c):
    """F -> R_q: coordinates (..., 4) to ring elements (..., 64)."""
    c = np.asarray(c, dtype=np.uint64)
    out = np.zeros(c.shape[:-1] + (D,), dtype=np.uint64)
    out[..., 0] = c[..., 0]
    for k in (1, 2, 3):
        out[..., 8 * k] = c[..., k]
        out[..., D - 8 * k] = (QU - c[..., k]) % QU
    return out


def theta(a):
    """Normalised trace R_q -> F: (..., 64) -> (..., 4)."""
    a = np.asarray(a, dtype=np.uint64)
    out = np.zeros(a.shape[:-1] + (4,), dtype=np.uint64)
    out[..., 0] = a[..., 0]
    for k in (1, 2, 3):
        d = (a[..., 8 * k] + QU - a[..., D - 8 * k]) % QU
        out[..., k] = (d * np.uint64(HALF)) % QU
    return out


_SPARSE = (0, 8, 16, 24, 40, 48, 56)


def fmul_ring(c, a):
    """c in F (coordinates (...,4)) times a in R_q (..., 64): seven shifted multiples."""
    e = emb(c)
    a = np.asarray(a, dtype=np.uint64)
    shape = np.broadcast_shapes(e.shape, a.shape)
    acc = np.zeros(shape, dtype=np.uint64)
    for u in _SPARSE:
        acc += (e[..., u:u + 1] * rot(a, u)) % QU
    return acc % QU


def fmul(c1, c2):
    """product in F"""
    return theta(fmul_ring(c1, emb(c2)))


def fadd(c1, c2):
    return (np.asarray(c1, dtype=np.uint64) + np.asarray(c2, dtype=np.uint64)) % QU


def fsub(c1, c2):
    return (np.asarray(c1, dtype=np.uint64) + QU - np.asarray(c2, dtype=np.uint64)) % QU


def fconst(x):
    return np.array([x % Q, 0, 0, 0], dtype=np.uint64)


FONE = fconst(1)
FZERO = fconst(0)


def theta_prod_sum(a, b):
    """sum_x Theta(a_x * b_x) for a, b of shape (M, 64): only the coefficients at multiples of 8 are
    computed.  Returns an element of F."""
    a = np.asarray(a, dtype=np.uint64)
    b = np.asarray(b, dtype=np.uint64)
    co = np.zeros(8, dtype=object)
    for k in range(8):
        tgt = 8 * k
        idx = (tgt - _IDX) % D
        pos = _IDX <= tgt
        prod = (a * b[:, idx]) % QU                      # (M, 64)
        sp = int(prod[:, pos].sum(dtype=object)) if pos.any() else 0
        sn = int(prod[:, ~pos].sum(dtype=object)) if (~pos).any() else 0
        co[k] = (sp - sn) % Q
    out = np.zeros(4, dtype=np.uint64)
    out[0] = co[0]
    for k in (1, 2, 3):
        out[k] = ((co[k] - co[8 - k]) * HALF) % Q
    return out


def eq_table(rho):
    """rho: (mu, 4) elements of F.  Returns the table eq(rho, k), k in [2^mu], as (2^mu, 4); the first
    variable is the least significant bit of k."""
    tab = FONE[None, :].copy()
    for j in range(rho.shape[0] - 1, -1, -1):
        hi = fmul(tab, rho[j][None, :])
        lo = fsub(tab, hi)
        new = np.empty((2 * tab.shape[0], 4), dtype=np.uint64)
        new[0::2] = lo
        new[1::2] = hi
        tab = new
    return tab


# ---------------------------------------------------------------- randomness from a transcript

class Transcript:
    def __init__(self, label=b"whippet-ref"):
        self.h = hashlib.shake_256()
        self.h.update(label)
        self.ctr = 0

    def absorb(self, data):
        if isinstance(data, np.ndarray):
            data = np.ascontiguousarray(data).tobytes()
        self.h.update(len(data).to_bytes(8, "little"))
        self.h.update(data)

    def seed(self, tag):
        h = self.h.copy()
        h.update(b"|" + tag + b"|" + self.ctr.to_bytes(8, "little"))
        self.ctr += 1
        s = h.digest(32)
        self.h.update(s)
        return s

    def rng(self, tag):
        return np.random.Generator(np.random.PCG64(int.from_bytes(self.seed(tag), "little")))


def uniform(rng, shape):
    return rng.integers(0, Q, size=shape, dtype=np.uint64)


def _ternary_ref(rng, shape):
    bits = rng.integers(0, 4, size=shape, dtype=np.int8)
    return (bits & 1) - (bits >> 1)


# Generator.integers(0, 4, dtype=int8) returns, for a PCG64 generator without a cached half word, the top two
# bits of the bytes of the raw 64-bit outputs in little-endian order (Lemire's method on buffered bytes).  The
# table maps a raw byte directly to the ternary value; ternary() uses it when this applies (same values).
_TLUT = {np.int8: np.array([((b >> 6) & 1) - ((b >> 6) >> 1) for b in range(256)], dtype=np.int8)}
_TLUT[np.float64] = _TLUT[np.int8].astype(np.float64)


def ternary(rng, shape, dtype=np.int8):
    """entries -1, 0, 1 with probabilities 1/4, 1/2, 1/4"""
    n = int(np.prod(shape))
    bg = rng.bit_generator
    if n % 8 == 0 and isinstance(bg, np.random.PCG64) and not bg.state['has_uint32']:
        return _TLUT[dtype][bg.random_raw(n // 8).view(np.uint8)].reshape(shape)
    return _ternary_ref(rng, shape).astype(dtype)


_ROOTS = np.exp(1j * np.pi * (2 * np.arange(D // 2) + 1) / D)
_VAND = _ROOTS[:, None] ** np.arange(D)[None, :]


def opnorm(c):
    return np.abs(_VAND @ np.asarray(c, dtype=np.float64)).max()


TAU, TOP = 71, 15.0


def challenge(rng):
    """element of the challenge set of [BS23]: 23 zeros, 31 +-1, 10 +-2, operator norm at most 15"""
    while True:
        c = np.zeros(D, dtype=np.int64)
        pos = rng.permutation(D)
        sg = rng.integers(0, 2, size=41) * 2 - 1
        c[pos[:31]] = sg[:31]
        c[pos[31:41]] = 2 * sg[31:]
        if opnorm(c) <= TOP:
            return c


def challenges(rng, r):
    return np.stack([challenge(rng) for _ in range(r)])
