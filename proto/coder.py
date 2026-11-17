"""Bit packing and a range coder with a discretised Gaussian model.

The range coder keeps a 64-bit range and 32-bit frequencies; the encoder accumulates `low` as a Python
integer, so no carry handling is needed.  Values outside [-K, K] are sent with an escape symbol followed
by 40 raw bits.
"""
from math import exp, ceil, log2
import numpy as np

PREC = 32
TOT = 1 << PREC
TOP = 1 << 64
BOT = 1 << 56


class Model:
    def __init__(self, sigma):
        self.sigma = max(float(sigma), 0.5)
        K = int(ceil(10 * self.sigma)) + 2
        self.K = K
        x = np.arange(-K, K + 1, dtype=np.float64)
        p = np.exp(-x * x / (2 * self.sigma ** 2))
        p /= p.sum()
        f = np.floor(p * (TOT - 4 * K - 16)).astype(np.int64) + 1
        esc = 1
        f[K] += TOT - int(f.sum()) - esc                 # put the rounding slack on the symbol 0
        self.freq = np.concatenate([f, [esc]])
        self.cum = np.concatenate([[0], np.cumsum(self.freq)])
        assert int(self.cum[-1]) == TOT and (self.freq > 0).all()

    def bits(self, v):
        """ideal code length of the vector v under the model"""
        v = np.asarray(v, dtype=np.int64)
        inside = np.abs(v) <= self.K
        idx = np.where(inside, v + self.K, 2 * self.K + 1)
        return float((PREC - np.log2(self.freq[idx].astype(np.float64))).sum() + 40 * (~inside).sum())


def encode(values, sigma):
    md = Model(sigma)
    K = md.K
    low, rng, nb = 0, TOP - 1, 64
    freq, cum = md.freq, md.cum
    for v in (int(t) for t in np.asarray(values).ravel()):
        syms = [(v + K, None)] if abs(v) <= K else [(2 * K + 1, v)]
        for sym, raw in syms:
            r = rng >> PREC
            low += r * int(cum[sym])
            rng = r * int(freq[sym])
            while rng < BOT:
                low <<= 8
                rng <<= 8
                nb += 8
            if raw is not None:
                for sh in (20, 0):                       # 40 raw bits, two 20-bit symbols
                    part = ((raw + (1 << 39)) >> sh) & ((1 << 20) - 1)
                    r = rng >> 20
                    low += r * part
                    rng = r
                    while rng < BOT:
                        low <<= 8
                        rng <<= 8
                        nb += 8
    t = rng.bit_length() - 2
    v = ((low >> t) + 1) << t
    nbits = nb - t
    nbytes = (nbits + 7) // 8
    code = (v >> t) << (nbytes * 8 - nbits)
    return code.to_bytes(nbytes, "big")


def decode(data, count, sigma):
    md = Model(sigma)
    K = md.K
    freq, cum = md.freq, md.cum
    pos = 0

    def nxt():
        nonlocal pos
        b = data[pos] if pos < len(data) else 0
        pos += 1
        return b

    rel = 0
    for _ in range(8):
        rel = (rel << 8) | nxt()
    rng = TOP - 1
    out = np.zeros(count, dtype=np.int64)
    for i in range(count):
        r = rng >> PREC
        cnt = min(rel // r, TOT - 1)
        sym = int(np.searchsorted(cum, cnt, side="right")) - 1
        rel -= r * int(cum[sym])
        rng = r * int(freq[sym])
        while rng < BOT:
            rel = (rel << 8) | nxt()
            rng <<= 8
        if sym == 2 * K + 1:
            raw = 0
            for _ in range(2):
                r = rng >> 20
                part = min(rel // r, (1 << 20) - 1)
                rel -= r * part
                rng = r
                while rng < BOT:
                    rel = (rel << 8) | nxt()
                    rng <<= 8
                raw = (raw << 20) | part
            out[i] = raw - (1 << 39)
        else:
            out[i] = sym - K
    return out


def pack_radix(values, base):
    """values in [0, base): one big integer in radix `base`, as bytes"""
    vals = [int(v) for v in np.asarray(values).ravel()]
    nbits = int(ceil(len(vals) * log2(base))) if vals else 0
    x = 0
    for v in reversed(vals):
        assert 0 <= v < base
        x = x * base + v
    return x.to_bytes((nbits + 7) // 8, "big")


def unpack_radix(data, count, base):
    x = int.from_bytes(data, "big")
    out = np.zeros(count, dtype=np.int64)
    for i in range(count):
        x, out[i] = divmod(x, base)
    if x:
        raise ValueError("non-canonical packing")
    return out


def radix_len(count, base):
    return (int(ceil(count * log2(base))) + 7) // 8
