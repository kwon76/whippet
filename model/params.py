#!/usr/bin/env python3
"""Proof-size model for LaBRADOR-style recursion.

The model follows the parameter selection described in [BS23, Sections
5.4-5.7] (decomposition bases, norm estimates, split into nu and mu parts,
last level without outer commitments) and adds the options studied in the
paper:

  E  entropy coding of Gaussian-like values (amortised opening, short
     garbage, projections),
  O  rounded outer commitments (low part joins the next witness),
  I  rounded inner commitments at the last level (low part recomputed by
     the verifier),
  G  derived garbage (one h polynomial at the last level and the constant
     coefficients of the lifting polynomials are not sent),
  P  norm bound gamma = sqrt(2 tau) beta for the amortised opening
     (completeness >= 1/2 for every witness, by Markov's inequality).

Two accountings of the Module-SIS norm bounds are supported:

  A  the bounds of [BS23, Theorem 5.1] together with Remark 5.2,
  B  bounds computed from the expected norms only: 8*T*gamma for the inner
     commitments and 2*beta' for the outer commitments.

Sizes are in bits internally; 1 KB = 1024 bytes.
"""
from math import sqrt, log2, ceil, pi, e, floor
import json
import sys

D = 64
LOGQ = 32.0
TAU = 71.0
T = 15.0
SL = sqrt(128 / 30)
K2 = ceil(128 / LOGQ)                 # number of lifted aggregates
KB = 8 * 1024.0
LOWDIV = 12.0                         # squared low-part coefficient: 4^rho/12 (expected) or 4^rho/4 (worst case)


def log2delta(b):
    return log2(((pi * b) ** (1 / b) * b / (2 * pi * e)) ** (1 / (2 * (b - 1))))


BKZ = 438                             # 0.292*b >= 128  <=>  b >= 438.36
LD = log2delta(BKZ)
UNIT = 4 * D * LOGQ * LD              # BKZ-438 fails iff log2(B)^2 < UNIT*kappa


def rank(B):
    """smallest Module-SIS rank over R_q that resists BKZ-438 at norm B"""
    lb = log2(B)
    if lb >= LOGQ - 0.01:
        return None
    k = max(1, ceil(lb * lb / UNIT))
    while 2 * sqrt(k * D * LOGQ * LD) <= lb:
        k += 1
    return k


def commitment_bits(B):
    """continuous version of rank(): bits of a binding commitment for norm B"""
    return log2(B) ** 2 / (4 * LD)


class Opt:
    def __init__(self, E=False, O=False, I=False, G=False, P=False,
                 eta=0.10, cap=1.10, acc='A', TL=False, LIN=False):
        self.TL = TL
        self.LIN = LIN      # linear relation: no quadratic garbage g_ij
        self.E, self.O, self.I, self.G, self.P = E, O, I, G, P
        self.eta, self.cap, self.acc = eta, cap, acc


def gbits(sigma, o):
    """bits per coefficient for a discretised Gaussian of standard deviation sigma"""
    sigma = max(sigma, 0.6)
    if o.E:
        # arithmetic coding with a discretised Gaussian model of width sigma: a vector of N
        # coefficients and norm at most sigma*sqrt(N) takes N*log2((sigma*sqrt(2 pi)+1)*sqrt(e)) bits
        return log2((sigma * sqrt(2 * pi) + 1) * sqrt(e))
    return log2(12 * sigma)


def hibits(dr):
    """bits per coefficient of the high part after dropping dr bits"""
    return log2(2 ** (LOGQ - dr) + 2) if dr > 0 else LOGQ


_LC = {}
_LA = {}


def level(n, r, beta, dec, o, sr, wide=True):
    """One level with outer commitments.  dec: z is decomposed in two parts."""
    key = (n, r, round(beta, 4), dec, round(sr, 6), wide, o.E, o.O, o.G, o.P, o.eta, o.acc, o.TL, o.LIN, BKZ)
    if key not in _LC:
        _LC[key] = _level(n, r, beta, dec, o, sr, wide)
    return _LC[key]


def _level(n, r, beta, dec, o, sr, wide):
    s = beta / sqrt(r * n * D)
    gam = beta * sqrt(TAU) * (sqrt(2) if o.P else 1.0)
    sigz = gam / sqrt(n * D)
    if dec:
        b = max(2, round(sqrt(sqrt(12) * sigz)))
        zn2 = 2 * gam * gam / (b * b)
        zf = b + 1
        base = b
    else:
        b = 0
        zn2 = gam * gam
        zf = 1
        base = max(2, round(sqrt(12) * sigz)) if wide else max(2, round(sqrt(sqrt(12) * sigz)))
    lb = log2(base)
    t1 = max(1, round(LOGQ / lb))
    b1 = ceil(2 ** (LOGQ / t1))
    gs = max(2.0, sqrt(24 * n * D) * s * s)
    t2 = max(1, round(log2(gs) / lb))
    b2 = max(2, round(gs ** (1 / t2)))
    ng = (r * r + r) // 2
    kap, k1, do = 8, 4, 0
    for _ in range(80):
        g1 = b1 * b1 * t1 / 12 * r * kap * D + (0.0 if o.LIN else b2 * b2 * t2 / 12 * ng * D)
        g2 = b1 * b1 * t1 / 12 * ng * D
        core = zn2 + g1 + g2
        do = 0
        if o.O:
            do = floor(0.5 * log2(max(1.0, o.eta * core * LOWDIV / (2 * k1 * D))))
            do = max(0, min(int(LOGQ) - 4, do))
        bp = sqrt(core + (2 * k1 * D * 4 ** do / LOWDIV if do > 0 else 0.0))
        if o.acc == 'A':
            nk1 = rank(2 * bp * sr)
            nkap = rank(8 * T * zf * bp * sr if o.LIN else max(8 * T * zf * bp * sr, 2 * zf * bp * sr + 4 * T * SL * beta))
        else:
            nk1 = rank(2 * bp)
            nkap = rank(8 * T * gam)
        if nk1 is None or nkap is None:
            return None
        if (nk1, nkap) == (k1, kap):
            break
        k1, kap = nk1, nkap
    outer = 2 * k1 * D * hibits(do)
    jl = 256 * gbits(beta / sqrt(2), o)
    lift = K2 * (3 if o.TL else D - (1 if o.G else 0)) * LOGQ
    chal = 4 * 128
    m = r * t1 * kap + (t1 + (0 if o.LIN else t2)) * ng + (2 * k1 if do > 0 else 0)
    return dict(kind='rec', n=n, r=r, beta=beta, s=s, b=b, t1=t1, t2=t2, b1=b1, b2=b2,
                kap=kap, k1=k1, Do=do, bp=bp, wide=wide, bits=outer + jl + lift + chal,
                outer=outer, jl=jl, lift=lift, chal=chal, m=m, dec=dec, base=base,
                wbits=n * r * D * log2(max(2.0, sqrt(12) * s)))


def last(n, r, beta, nu, o):
    """Last level: no outer commitments, interactive garbage, final message in the clear."""
    key = (n, r, round(beta, 4), nu, o.E, o.I, o.G, o.P, o.cap, o.acc == 'C', o.TL, o.LIN, BKZ)
    if key not in _LA:
        _LA[key] = _last(n, r, beta, nu, o)
    return _LA[key]


def _last(n, r, beta, nu, o):
    s = beta / sqrt(r * n * D)
    gam = beta * sqrt(TAU) * (sqrt(2) if o.P else 1.0)
    sigz = gam / sqrt(n * D)
    best = None
    for di in (range(0, 16) if o.I else (0,)):
        kap = 8
        ok = True
        for _ in range(80):
            l2 = r * kap * D * 4 ** di / LOWDIV if di > 0 else 0.0
            e2 = TAU * l2 * (2 if o.P else 1)
            gb = sqrt(gam * gam + e2)
            bb = sqrt(beta * beta + l2)
            nk = rank(2 * gb) if o.acc == 'C' else rank(8 * T * gb if o.LIN else max(8 * T * gb, 2 * gb + 4 * T * SL * bb))
            if nk is None:
                ok = False
                break
            if nk == kap:
                break
            kap = nk
        if not ok:
            continue
        ratio = bb / beta
        if ratio > o.cap + 1e-9:
            continue
        zb = n * D * gbits(sigz, o)
        tb = r * kap * D * hibits(di)
        sg = sqrt(2 * n * D) * s * s * TAU * r
        gbt = 0.0 if o.LIN else (2 * nu + 1) * D * min(LOGQ, gbits(sg, o))
        hb = (2 * r - 1 - (1 if o.G else 0)) * D * LOGQ
        jl = 256 * gbits(bb / sqrt(2), o)
        lift = K2 * (3 if o.TL else D - (1 if o.G else 0)) * LOGQ
        chal = 4 * 128
        bits = zb + tb + gbt + hb + jl + lift + chal
        cand = dict(kind='last', n=n, r=r, beta=beta, s=s, kap=kap, Di=di, ratio=ratio, bits=bits,
                    z=zb, t=tb, g=gbt, h=hb, jl=jl, lift=lift, chal=chal, nu=nu,
                    wbits=n * r * D * log2(max(2.0, sqrt(12) * s)))
        if best is None or cand['bits'] < best['bits']:
            best = cand
    return best


_SP = {}


def all_splits(n, m, dec, rmax):
    """All shapes (nu, mu, n', r') of the next witness with multiplicity at most rmax."""
    key = (n, m, dec, rmax)
    if key not in _SP:
        res = {}
        f = 2 if dec else 1
        for nu in range(1, n + 1):
            npz = ceil(n / nu)
            mu = max(1, ceil(m / npz))
            n2 = max(npz, ceil(m / mu))
            r2 = f * nu + mu
            if r2 > rmax or npz < 2:
                break
            res.setdefault((n2, r2), (nu, mu, n2, r2))
        _SP[key] = list(res.values())
    return _SP[key]


def ranked(n, lv, dec, o, rmax, k):
    """The k shapes whose next level produces the smallest witness (exact one-level look-ahead)."""
    out = []
    for nu, mu, n2, r2 in all_splits(n, lv['m'], dec, rmax):
        nx = level(n2, r2, lv['bp'], True, o, SL, True)
        if nx is None:
            continue
        score = (2 * n2 + nx['m']) * D * log2(max(2, nx['base'])) + nx['bits']
        out.append((score, nu, mu, n2, r2))
    out.sort()
    return [x[1:] for x in out[:k]]


def chain(W0, beta0, L, o, thorough=False):
    """Best total over the split choices for a recursion with L levels."""
    best = [None]

    def rec(n, r, beta, lev, acc, trace, nu_prev):
        if lev == L:
            la = last(n, r, beta, nu_prev, o)
            if la is None:
                return
            tot = acc + la['bits']
            if best[0] is None or tot < best[0][0]:
                best[0] = (tot, trace + [la])
            return
        final = (lev == L - 1)
        sr = SL * (o.cap if (o.I and final) else 1.0)
        variants = [(True, True), (False, True), (False, False)] if final else [(True, True)]
        for dec, wide in variants:
            lv = level(n, r, beta, dec, o, sr, wide)
            if lv is None:
                continue
            if best[0] is not None and acc + lv['bits'] > best[0][0]:
                continue
            if final:
                cands = all_splits(n, lv['m'], dec, 20 if thorough else 14)
            elif thorough:
                k = 5 if lev == L - 2 else (3 if lev == 1 else (2 if lev <= 3 else 1))
                cands = ranked(n, lv, dec, o, 400 if lev == 1 else 80, k)
            else:
                k = 4 if lev == L - 2 else (3 if lev == 1 else (2 if lev == 2 else 1))
                cands = ranked(n, lv, dec, o, 400 if lev == 1 else 80, k)
            for nu, mu, n2, r2 in cands:
                rec(n2, r2, lv['bp'], lev + 1, acc + lv['bits'],
                    trace + [dict(lv, nu=nu, mu=mu)], 2 * nu if dec else nu)

    ropt = max(4.0, (2 * W0 / 12) ** (1 / 3))
    facs = ((0.35, 0.45, 0.55, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0) if thorough
            else (0.4, 0.55, 0.75, 1.0, 1.3, 1.7, 2.2, 3.0))
    for fac in facs:
        r0 = max(2, round(ropt * fac))
        rec(ceil(W0 / r0), r0, beta0, 1, 0.0, [], 0)
    return best[0]


def optimise(W0, beta0, o, Ls=range(5, 9), thorough=False):
    best = None
    for L in Ls:
        c = chain(W0, beta0, L, o, thorough)
        if c is not None and (best is None or c[0] < best[0]):
            best = c
    return best


def grid(flags):
    etas = (0.1, 0.3) if flags.get('O') else (0.1,)
    caps = (1.10, 1.25) if flags.get('I') else (1.0,)
    return [(eta, cap) for eta in etas for cap in caps]


def tune(W0, beta0, thorough=False, **flags):
    """Optimise the recursion shape and the two rounding hyper-parameters."""
    best = None
    for eta, cap in grid(flags):
        o = Opt(eta=eta, cap=cap, **flags)
        c = optimise(W0, beta0, o, thorough=thorough)
        if c is not None and (best is None or c[0] < best[0]):
            best = (c[0], c[1], o)
    return best


def replay(tr, W0, beta0, **flags):
    """Evaluate a fixed recursion shape (taken from the trace tr) under other options."""
    best = None
    for eta, cap in grid(flags):
        o = Opt(eta=eta, cap=cap, **flags)
        r = tr[0]['r']
        n = ceil(W0 / r)
        beta = beta0
        out = []
        tot = 0.0
        nu_prev = 0
        ok = True
        for i, lv0 in enumerate(tr):
            if lv0['kind'] == 'rec':
                final = (i == len(tr) - 2)
                sr = SL * (o.cap if (o.I and final) else 1.0)
                lv = level(n, r, beta, lv0['dec'], o, sr, lv0['wide'])
                if lv is None:
                    ok = False
                    break
                nu, mu = lv0['nu'], lv0['mu']
                out.append(dict(lv, nu=nu, mu=mu))
                tot += lv['bits']
                f = 2 if lv0['dec'] else 1
                n = max(ceil(n / nu), ceil(lv['m'] / mu))
                r = f * nu + mu
                beta = lv['bp']
                nu_prev = f * nu
            else:
                la = last(n, r, beta, nu_prev, o)
                if la is None:
                    ok = False
                    break
                out.append(la)
                tot += la['bits']
        if ok and (best is None or tot < best[0]):
            best = (tot, out, o)
    return best


def show(tr):
    for i, lv in enumerate(tr):
        if lv['kind'] == 'rec':
            print(" L%d n=%6d r=%3d beta=2^%.2f s=%6.2f b=%3d t1=%d t2=%d kap=%2d k1=%d Do=%2d "
                  "out=%.2fKB wit=%.1fKB nu=%d mu=%d dec=%d" % (
                      i + 1, lv['n'], lv['r'], log2(lv['beta']), lv['s'], lv['b'], lv['t1'], lv['t2'],
                      lv['kap'], lv['k1'], lv['Do'], lv['bits'] / KB, lv['wbits'] / KB, lv['nu'],
                      lv['mu'], lv['dec']))
        else:
            print(" L%d n=%6d r=%3d beta=2^%.2f s=%6.2f kap=%2d Di=%d out=%.2fKB "
                  "(z %.2f t %.2f g %.2f h %.2f) wit=%.1fKB ratio=%.3f" % (
                      i + 1, lv['n'], lv['r'], log2(lv['beta']), lv['s'], lv['kap'], lv['Di'],
                      lv['bits'] / KB, lv['z'] / KB, lv['t'] / KB, lv['g'] / KB, lv['h'] / KB,
                      lv['wbits'] / KB, lv['ratio']))


def components(tr):
    """aggregate a trace into the components used in the breakdown figure"""
    c = dict(outer=0.0, jl=0.0, lift=0.0, chal=0.0, z=0.0, t=0.0, gh=0.0)
    for lv in tr:
        c['jl'] += lv['jl']
        c['lift'] += lv['lift']
        c['chal'] += lv['chal']
        if lv['kind'] == 'rec':
            c['outer'] += lv['outer']
        else:
            c['z'] += lv['z']
            c['t'] += lv['t']
            c['gh'] += lv['g'] + lv['h']
    return {k: v / KB for k, v in c.items()}


if __name__ == "__main__":
    print("log2 delta(%d) = %.6f, rank unit = %.3f" % (BKZ, LD, UNIT))
    W0 = 2 ** 22
    beta0 = 2.0 ** 14                    # binary R1CS with 2^25 constraints [BS23, Table 3]
    for acc in ('A', 'B'):
        tot, tr, o = tune(W0, beta0, thorough=True, acc=acc)
        print("baseline, accounting %s: %.2f KB (+0.5 KB R1CS reduction), %d levels" % (
            acc, tot / KB, len(tr)))
        show(tr)
        for nm, fl in (("E", dict(E=True)), ("G", dict(G=True)), ("O", dict(O=True)), ("I", dict(I=True)),
                       ("EOIG", dict(E=True, O=True, I=True, G=True))):
            t2, tr2, o2 = replay(tr, W0, beta0, acc=acc, **fl)
            print("  same shape, %-5s %.2f KB (%+.2f)" % (nm, t2 / KB, (t2 - tot) / KB))
        t3, tr3, o3 = tune(W0, beta0, acc=acc, E=True, O=True, I=True, G=True)
        print("  re-optimised, all: %.2f KB" % (t3 / KB))
