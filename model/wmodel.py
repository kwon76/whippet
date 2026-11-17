#!/usr/bin/env python3
"""Proof-size model for the linear tail ("Whippet").

One level of the tail takes a witness of m ring elements with Euclidean norm at most beta, splits it
into r blocks of rank n = ceil(m/r) and sends

  u1   outer commitment to the (rounded) inner commitments t_i = hi(A s_i),
  p    256 Johnson-Lindenstrauss projections,
  B_k  K2 lifted aggregates, each an element of the field F = F_{q^fF} with one coordinate fixed,
  g_j  the sum-check round polynomials over F (two elements of F per round, log2(r n) rounds),
  u2   outer commitment to the r partial evaluations W_i in R_q (or the W_i themselves),

and hands the amortised opening z = sum c_i s_i together with the openings of u1 and u2 to the next
level.  The last level sends the inner commitments, the W_i and z.

Accountings of the Module-SIS norms (same as in params.py):
  A  bounds required by the security proof, including the projection slack of the next level,
  B  expected norms, 8T*gamma for inner commitments, 2*beta' for outer commitments,
  C  as B, but the inner commitments of the last level only for 2*gamma.
"""
from math import sqrt, log2, ceil, floor
import params as P


class WOpt:
    def __init__(self, E=True, O=True, I=True, Pm=False, acc='A', eta=0.10, eta_in=0.10, cap=1.10,
                 fF=4, inner=True, clear=True, gc=False, theta=2.0, tmar=128 ** -0.5):
        # gc: norm bounds that hold for every witness (deterministic bounds for digits and low parts, a joint
        # Markov bound with factor theta for (z, e), projection threshold raised by (1 + tmar)); accounting A
        self.gc, self.theta, self.tmar = gc, theta, tmar
        self.E, self.O, self.I, self.P = E, O, I, Pm
        self.acc, self.eta, self.eta_in, self.cap, self.fF = acc, eta, eta_in, cap, fF
        self.inner, self.clear = inner, clear


def set_bkz(b):
    P.BKZ = b
    P.LD = P.log2delta(b)
    P.UNIT = 4 * P.D * P.LOGQ * P.LD
    _C.clear()
    _T.clear()


def felt():
    return P.LOGQ


_C = {}
_T = {}


def wlevel(m, beta, r, dec, wmode, ein, o):
    key = (m, round(beta, 3), r, dec, wmode, ein, o.E, o.O, o.P, o.acc, o.eta, o.cap, o.fF, P.BKZ,
           o.gc, o.theta, o.tmar)
    if key not in _C:
        _C[key] = (_wlevel_gc if o.gc else _wlevel)(m, beta, r, dec, wmode, ein, o)
    return _C[key]


def _wlevel(m, beta, r, dec, wmode, ein, o):
    D, LOGQ, TAU, T, SL = P.D, P.LOGQ, P.TAU, P.T, P.SL
    n = ceil(m / r)
    if n < 2:
        return None
    s = beta / sqrt(r * n * D)
    pm = 2.0 if o.P else 1.0
    gam2 = beta * beta * TAU * pm
    sigz = sqrt(gam2 / (n * D))
    if dec:
        b = max(2, round(sqrt(sqrt(12) * sigz)))
        zn2 = 2 * gam2 / (b * b)
        zf = b + 1
        base = b
    else:
        b = 0
        zn2 = gam2
        zf = 1
        base = max(2, round(sqrt(12) * sigz))
    lb = log2(base)
    commit = (wmode == 'commit')
    opn = (wmode == 'open')
    sr = SL
    kap, k1 = 8, 4
    for _ in range(100):
        di = 0
        if ein > 0:
            di = floor(0.5 * log2(max(1.0, ein * beta * beta * P.LOWDIV / (r * kap * D))))
            di = max(0, min(int(LOGQ) - 4, di, floor(log2(max(1.0, sqrt(12) * s)))))
        lam2 = r * kap * D * 4.0 ** di / P.LOWDIV if di > 0 else 0.0
        bbar2 = beta * beta + lam2
        game2 = TAU * pm * lam2
        gbar = sqrt(gam2 + game2)
        tb = LOGQ - di
        t1 = max(1, round(tb / lb))
        b1 = ceil(2 ** (tb / t1))
        nt2 = 0.0 if opn else b1 * b1 * t1 / 12 * r * kap * D
        tw = max(1, round(LOGQ / lb))
        bw = ceil(2 ** (LOGQ / tw))
        nw2 = bw * bw * tw / 12 * r * D if commit else 0.0
        te, ne2 = 0, 0.0
        if di > 0:
            if dec:
                te = 2
                ne2 = kap * D * (b * b / 12 + game2 / (kap * D) / (b * b))
            else:
                te = 1
                ne2 = game2
        core = zn2 + nt2 + nw2 + ne2
        nout = 0 if opn else k1 * (2 if commit else 1)
        do = 0
        if o.O and not opn:
            do = floor(0.5 * log2(max(1.0, o.eta * core * P.LOWDIV / (nout * D))))
            do = max(0, min(int(LOGQ) - 4, do))
        bp = sqrt(core + (nout * D * 4.0 ** do / P.LOWDIV if do > 0 else 0.0))
        if o.acc == 'A':
            nk1 = P.rank(2 * bp * sr)
            nkap = P.rank(8 * T * zf * bp * sr)
        else:
            nk1 = P.rank(2 * bp)
            nkap = P.rank(8 * T * gbar)
        if opn:
            nk1 = k1
        if nk1 is None or nkap is None:
            return None
        if (nk1, nkap) == (k1, kap):
            break
        k1, kap = nk1, nkap
    else:
        return None
    mu = max(1, ceil(log2(n)))
    f = felt()
    outer = r * kap * D * P.hibits(di) if opn else k1 * D * P.hibits(do)
    wbits = k1 * D * P.hibits(do) if commit else r * D * LOGQ
    jl = 256 * P.gbits(beta / sqrt(2), o)
    lift = P.K2 * (o.fF - 1) * f
    sc = 2 * o.fF * f * mu
    chal = 4 * 128
    m2 = (2 if dec else 1) * n + (0 if opn else r * kap * t1) + (r * tw if commit else 0) + kap * te \
        + (nout if do > 0 else 0)
    return dict(kind='rec', m=m, n=n, r=r, beta=beta, s=s, b=b, base=base, dec=dec, wmode=wmode,
                kap=kap, k1=k1, di=di, do=do, t1=t1, tw=tw, te=te, mu=mu, bp=bp, m2=m2,
                b1=b1, bw=bw, gam2=gam2, game2=game2, lam2=lam2, ein=ein,
                outer=outer, garb=wbits, jl=jl, lift=lift, sc=sc, chal=chal,
                bits=outer + wbits + jl + lift + sc + chal)


def wlast(m, beta, r, o):
    key = (m, round(beta, 3), r, o.E, o.I, o.P, o.acc, o.cap, o.fF, P.BKZ, o.gc, o.theta, o.tmar)
    if key not in _T:
        _T[key] = (_wlast_gc if o.gc else _wlast)(m, beta, r, o)
    return _T[key]


def _wlast(m, beta, r, o):
    D, LOGQ, TAU, T, SL = P.D, P.LOGQ, P.TAU, P.T, P.SL
    n = ceil(m / r)
    if n < 1:
        return None
    s = beta / sqrt(r * n * D)
    pm = 2.0 if o.P else 1.0
    gam2 = beta * beta * TAU * pm
    sigz = sqrt(gam2 / (n * D))
    mu = max(1, ceil(log2(n)))
    f = felt()
    best = None
    for di in (range(0, 17) if o.I else (0,)):
        kap = 8
        ok = False
        for _ in range(100):
            lam2 = r * kap * D * 4.0 ** di / P.LOWDIV if di > 0 else 0.0
            e2 = TAU * pm * lam2
            gb = sqrt(gam2 + e2)
            bb = sqrt(beta * beta + lam2)
            nk = P.rank(2 * gb) if o.acc == 'C' else P.rank(8 * T * gb)
            if nk is None:
                break
            if nk == kap:
                ok = True
                break
            kap = nk
        if not ok:
            continue
        zb = n * D * P.gbits(sigz, o)
        tb = r * kap * D * P.hibits(di)
        wb = r * D * LOGQ
        jl = 256 * P.gbits(beta / sqrt(2), o)
        lift = P.K2 * (o.fF - 1) * f
        sc = 2 * o.fF * f * mu
        chal = 4 * 128
        bits = zb + tb + wb + jl + lift + sc + chal
        c = dict(kind='last', m=m, n=n, r=r, beta=beta, s=s, kap=kap, di=di, mu=mu, ratio=bb / beta,
                 gam2=gam2, game2=e2, lam2=lam2,
                 z=zb, t=tb, garb=wb, jl=jl, lift=lift, sc=sc, chal=chal, bits=bits)
        if best is None or bits < best['bits']:
            best = c
    return best


def _wlevel_gc(m, B, r, dec, wmode, ein, o):
    """A recursive level with norm bounds that hold for every witness of norm at most B.
    Digits and low parts are counted with their maximal absolute values; the prover repeats the last
    challenge until Q = ||z||^2 + ||e||^2 <= theta*tau*(B^2 + Lambda^2), which happens with probability
    at least 1 - 1/theta (Markov); the verifier accepts projections of norm sqrt(128 (1 + tmar)) B, which
    the honest projection satisfies with probability at least 128 tmar^2 / (1 + 128 tmar^2) (Cantelli).
    Ranks as in accounting A, for the verifier's bound sqrt(1 + tmar) * B' of the next level."""
    D, LOGQ, TAU, T, SL = P.D, P.LOGQ, P.TAU, P.T, P.SL
    n = ceil(m / r)
    if n < 2:
        return None
    tm = sqrt(1.0 + o.tmar)
    s = B / sqrt(r * n * D)
    commit = (wmode == 'commit')
    kap, k1 = 8, 4
    for _ in range(100):
        di = 0
        if ein > 0:
            di = floor(0.5 * log2(max(1.0, ein * B * B * 4 / (r * kap * D))))
            di = max(0, min(int(LOGQ) - 4, di, floor(log2(max(1.0, sqrt(12) * s)))))
        lam2 = r * kap * D * 4.0 ** di / 4 if di > 0 else 0.0
        qmax = o.theta * TAU * (B * B + lam2)
        nz = n + (kap if di > 0 else 0)
        if dec:
            b = max(2, round(sqrt(2 * sqrt(qmax / (nz * D)))))
            det0 = nz * D * b * b / 4.0
            zn2 = det0 + ((sqrt(qmax) + sqrt(det0)) / b) ** 2
            zf, base = b + 1, b
        else:
            b, zn2, zf = 0, qmax, 1
            base = max(2, round(sqrt(12 * qmax / (nz * D))))
        lb = log2(base)
        tb = LOGQ - di
        t1 = max(1, round(tb / lb))
        b1 = ceil(2 ** (tb / t1))
        top1 = (2.0 ** (tb - 1) + 0.5) / float(b1) ** (t1 - 1) + 1
        nt2 = r * kap * D * ((t1 - 1) * b1 * b1 / 4.0 + top1 * top1)
        tw = max(1, round(LOGQ / lb))
        bw = ceil(2 ** (LOGQ / tw))
        topw = 2.0 ** (LOGQ - 1) / float(bw) ** (tw - 1) + 1
        nw2 = r * D * ((tw - 1) * bw * bw / 4.0 + topw * topw) if commit else 0.0
        te = (2 if dec else 1) if di > 0 else 0
        core = zn2 + nt2 + nw2
        nout = k1 * (2 if commit else 1)
        do = 0
        if o.O:
            do = floor(0.5 * log2(max(1.0, o.eta * core * 4 / (nout * D))))
            do = max(0, min(int(LOGQ) - 4, do))
        bp = sqrt(core + (nout * D * 4.0 ** do / 4 if do > 0 else 0.0))
        nk1 = P.rank(2 * bp * tm * SL)
        nkap = P.rank(8 * T * zf * bp * tm * SL)
        if nk1 is None or nkap is None:
            return None
        if (nk1, nkap) == (k1, kap):
            break
        k1, kap = nk1, nkap
    else:
        return None
    mu = max(1, ceil(log2(n)))
    f = felt()
    outer = k1 * D * P.hibits(do)
    wbits = k1 * D * P.hibits(do) if commit else r * D * LOGQ
    jl = 256 * P.gbits(B * tm / sqrt(2), o)
    lift = P.K2 * (o.fF - 1) * f
    sc = 2 * o.fF * f * mu
    chal = 4 * 128
    m2 = (2 if dec else 1) * n + r * kap * t1 + (r * tw if commit else 0) + kap * te + (nout if do > 0 else 0)
    return dict(kind='rec', m=m, n=n, r=r, beta=B, bjl=B * tm, s=s, b=b, base=base, dec=dec, wmode=wmode,
                kap=kap, k1=k1, di=di, do=do, t1=t1, tw=tw, te=te, mu=mu, bp=bp, m2=m2,
                b1=b1, bw=bw, gam2=qmax, game2=0.0, lam2=lam2, qmax=qmax, ein=ein, gc=True,
                outer=outer, garb=wbits, jl=jl, lift=lift, sc=sc, chal=chal,
                bits=outer + wbits + jl + lift + sc + chal)


def _wlast_gc(m, B, r, o):
    """Terminal level with the joint check ||z||^2 + ||e||^2 <= theta*tau*(B^2 + Lambda^2)."""
    D, LOGQ, TAU, T = P.D, P.LOGQ, P.TAU, P.T
    n = ceil(m / r)
    if n < 1:
        return None
    tm = sqrt(1.0 + o.tmar)
    s = B / sqrt(r * n * D)
    mu = max(1, ceil(log2(n)))
    f = felt()
    best = None
    for di in (range(0, 17) if o.I else (0,)):
        kap = 8
        ok = False
        for _ in range(100):
            lam2 = r * kap * D * 4.0 ** di / 4 if di > 0 else 0.0
            qmax = o.theta * TAU * (B * B + lam2)
            nk = P.rank(8 * T * sqrt(qmax))
            if nk is None:
                break
            if nk == kap:
                ok = True
                break
            kap = nk
        if not ok:
            continue
        zb = n * D * P.gbits(sqrt(qmax / (n * D)), o)
        tb = r * kap * D * P.hibits(di)
        wb = r * D * LOGQ
        jl = 256 * P.gbits(B * tm / sqrt(2), o)
        lift = P.K2 * (o.fF - 1) * f
        sc = 2 * o.fF * f * mu
        chal = 4 * 128
        bits = zb + tb + wb + jl + lift + sc + chal
        c = dict(kind='last', m=m, n=n, r=r, beta=B, bjl=B * tm, s=s, kap=kap, di=di, mu=mu,
                 ratio=sqrt(1 + lam2 / (B * B)), gam2=qmax, game2=0.0, lam2=lam2, qmax=qmax, joint=True, gc=True,
                 z=zb, t=tb, garb=wb, jl=jl, lift=lift, sc=sc, chal=chal, bits=bits)
        if best is None or bits < best['bits']:
            best = c
    return best


def best_last(m, beta, o, rs=(1, 2, 3, 4, 5, 6, 7, 8, 10, 12)):
    best = None
    for r in rs:
        if r > m:
            continue
        la = wlast(m, beta, r, o)
        if la is not None and (best is None or la['bits'] < best['bits']):
            best = la
    return best


def rcands(m):
    rs = set(range(2, 13))
    r0 = sqrt(max(4.0, m / 20.0))
    for fac in (0.25, 0.35, 0.5, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.2, 4.0):
        rs.add(max(2, round(r0 * fac)))
    return sorted(x for x in rs if 2 * x <= m)


def search(W0, beta0, o, beam=120, depth=8):
    """Beam search over the shapes of the recursion; returns (bits, trace)."""
    best = None
    la = best_last(W0, beta0, o)
    if la is not None:
        best = (la['bits'], [la])
    states = [(0.0, W0, beta0, [])]
    for _ in range(depth):
        new = {}
        for acc, m, beta, tr in states:
            for r in rcands(m):
                for dec in (True, False):
                    for wmode in (('commit', 'clear') if o.clear else ('commit',)):
                        if wmode == 'clear' and r > 16:
                            continue
                        for ein in ((0.0, o.eta_in) if o.inner else (0.0,)):
                            lv = wlevel(m, beta, r, dec, wmode, ein, o)
                            if lv is None or lv['m2'] >= m:
                                continue
                            a2 = acc + lv['bits']
                            if best is not None and a2 >= best[0]:
                                continue
                            key = (lv['m2'], round(log2(lv['bp']), 2))
                            if key not in new or a2 < new[key][0]:
                                new[key] = (a2, lv['m2'], lv['bp'], tr + [lv])
        cand = []
        for a2, m2, bp, tr in new.values():
            la = best_last(m2, bp, o)
            if la is None:
                continue
            tot = a2 + la['bits']
            if best is None or tot < best[0]:
                best = (tot, tr + [la])
            cand.append((tot, a2, m2, bp, tr))
        cand.sort(key=lambda x: x[0])
        states = [(a2, m2, bp, tr) for _, a2, m2, bp, tr in cand[:beam]]
        if not states:
            break
    return best


def comps(tr):
    c = dict(outer=0.0, garb=0.0, jl=0.0, lift=0.0, sc=0.0, chal=0.0, z=0.0, t=0.0)
    for lv in tr:
        for k in ('jl', 'lift', 'sc', 'chal'):
            c[k] += lv[k]
        if lv['kind'] == 'rec':
            if lv['wmode'] == 'open':
                c['t'] += lv['outer']
            else:
                c['outer'] += lv['outer'] + (lv['garb'] if lv['wmode'] == 'commit' else 0.0)
            c['garb'] += lv['garb'] if lv['wmode'] != 'commit' else 0.0
        else:
            c['z'] += lv['z']
            c['t'] += lv['t']
            c['garb'] += lv['garb']
    return {k: v / P.KB for k, v in c.items()}


def show(tr):
    for i, lv in enumerate(tr):
        if lv['kind'] == 'rec':
            print(" L%d m=%7d n=%6d r=%4d beta=2^%.2f b=%3d dec=%d W=%-6s kap=%2d k1=%d di=%2d do=%2d "
                  "t1=%d mu=%2d out=%.2fKB -> m'=%d" % (
                      i + 1, lv['m'], lv['n'], lv['r'], log2(lv['beta']), lv['b'], lv['dec'], lv['wmode'],
                      lv['kap'], lv['k1'], lv['di'], lv['do'], lv['t1'], lv['mu'], lv['bits'] / P.KB,
                      lv['m2']))
        else:
            print(" L%d m=%7d n=%6d r=%4d beta=2^%.2f kap=%2d di=%2d mu=%2d out=%.2fKB (z %.2f t %.2f W %.2f "
                  "jl %.2f sc %.2f)" % (
                      i + 1, lv['m'], lv['n'], lv['r'], log2(lv['beta']), lv['kap'], lv['di'], lv['mu'],
                      lv['bits'] / P.KB, lv['z'] / P.KB, lv['t'] / P.KB, lv['garb'] / P.KB, lv['jl'] / P.KB,
                      lv['sc'] / P.KB))


if __name__ == "__main__":
    import run_eval as RE
    import sys
    for bkz in (438,):
        set_bkz(bkz)
        for logN in (26, 28, 30):
            for acc in "CBA":
                g = RE.greyhound_instance(logN, True)
                o = WOpt(acc=acc)
                tot, tr = search(g['W'], g['beta'], o)
                own = (g['n1'] * P.D * P.hibits(g['Du']) + (o.fF - 1) * P.LOGQ + 2 * 128) / P.KB
                print("N=2^%d acc %s bkz %d: tail %.2f KB + own %.2f KB = %.2f KB (%d levels)" % (
                    logN, acc, bkz, tot / P.KB, own, tot / P.KB + own, len(tr)))
                if logN == 30:
                    show(tr)
                    print("  ", {k: round(v, 2) for k, v in comps(tr).items()})
                sys.stdout.flush()
