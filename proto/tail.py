"""Reference implementation of the linear tail: non-interactive prover and verifier.

The implementation follows the protocol of the paper step by step and is written for clarity, not for
speed.  Parameters of every level (ranks, bases, numbers of dropped bits, norm bounds) are taken from
the size model (../model/wmodel.py), so that measured proof sizes can be compared with the model.
"""
import hashlib
import os
import sys
import time
from math import sqrt, log2, ceil

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import params as P          # noqa: E402
import wmodel as W          # noqa: E402
import ring as R            # noqa: E402
import coder as C           # noqa: E402
import fastring as F        # noqa: E402

Q = R.Q
QU = R.QU
D = R.D
JL = 256
CHUNK = 1 << 16


# ------------------------------------------------------------------ small helpers

def crs(tag, shape):
    seed = int.from_bytes(hashlib.shake_256(b"crs|" + tag).digest(32), "little")
    return R.uniform(np.random.Generator(np.random.PCG64(seed)), shape)


def round_hi(t, rho):
    """t in R_q (uint64).  Returns (H, L) with centred(t) = 2^rho H + L, |L| <= 2^(rho-1)."""
    tc = R.center(t)
    if rho == 0:
        return tc, np.zeros_like(tc)
    H = (tc + (1 << (rho - 1))) >> rho
    return H, tc - (H << rho)


def hi_base(rho):
    """number of values of H + offset; offset = 2^(31-rho) + 1"""
    return (1 << (32 - rho)) + 3, (1 << (31 - rho)) + 1


def digits(x, base, k):
    """centred base-`base` digits of the integer array x; the last digit is the remaining quotient"""
    x = np.asarray(x, dtype=np.int64).copy()
    out = []
    hb = base // 2
    for _ in range(k - 1):
        dgt = ((x + hb) % base) - hb
        out.append(dgt)
        x = (x - dgt) // base
    out.append(x)
    return np.stack(out)


def recompose(dg, base):
    acc = np.zeros_like(dg[0])
    for j in range(dg.shape[0] - 1, -1, -1):
        acc = acc * base + dg[j]
    return acc


def mul_int(a, b):
    """exact product in R of small integer arrays (int64)"""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    acc = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.int64)
    for u in range(D):
        au = a[..., u:u + 1]
        if not au.any():
            continue
        r = b[..., R.RIDX[u]]
        if u:
            r = r.copy()
            r[..., :u] = -r[..., :u]
        acc += au * r
    return acc


def sigma_m1(v):
    """the automorphism X -> X^-1 on (..., 64) arrays with entries in [0, q)"""
    out = np.empty_like(v)
    out[..., 0] = v[..., 0]
    out[..., 1:] = (QU - v[..., :0:-1]) % QU
    return out


def fsum(terms):
    acc = np.zeros(4, dtype=object)
    for t in terms:
        acc = acc + np.asarray(t).astype(object)
    return (acc % Q).astype(np.uint64)


def theta_dot(a, b):
    """Theta(<a, b>) for a, b of shape (M, 64) with entries in [0, q)"""
    a = np.asarray(a, dtype=np.uint64)
    b = np.asarray(b, dtype=np.uint64)
    co = []
    idx_all = np.arange(D)
    for k in range(8):
        tgt = 8 * k
        idx = (tgt - idx_all) % D
        prod = (a * b[:, idx]) % QU
        col = prod.sum(axis=0, dtype=np.uint64)            # each column sum < 2^32 * M
        pos = idx_all <= tgt
        co.append((int(col[pos].sum(dtype=object)) - int(col[~pos].sum(dtype=object))) % Q)
    out = np.zeros(4, dtype=np.uint64)
    out[0] = co[0]
    for k in (1, 2, 3):
        out[k] = ((co[k] - co[8 - k]) * R.HALF) % Q
    return out


# ------------------------------------------------------------------ proof container

class Writer:
    def __init__(self):
        self.parts = []
        self.sizes = {}

    def put(self, kind, data):
        self.parts.append(bytes(data))
        self.sizes[kind] = self.sizes.get(kind, 0) + len(data)

    def put_var(self, kind, data):
        self.put(kind, len(data).to_bytes(3, "big") + bytes(data))

    def bytes(self):
        return b"".join(self.parts)


class Reader:
    def __init__(self, data):
        self.data = data
        self.pos = 0

    def get(self, n):
        if self.pos + n > len(self.data):
            raise ValueError("proof too short")
        out = self.data[self.pos:self.pos + n]
        self.pos += n
        return out

    def get_var(self):
        n = int.from_bytes(self.get(3), "big")
        return self.get(n)


def put_ring(wr, kind, a):
    wr.put(kind, np.ascontiguousarray(np.asarray(a, dtype="<u4")).tobytes())


def get_ring(rd, shape):
    n = int(np.prod(shape))
    a = np.frombuffer(rd.get(4 * n), dtype="<u4").astype(np.uint64).reshape(shape)
    if (a >= QU).any():
        raise ValueError("non-canonical element")
    return a


def put_hi(wr, kind, H, rho):
    base, off = hi_base(rho)
    wr.put(kind, C.pack_radix(H + off, base))


def get_hi(rd, shape, rho):
    base, off = hi_base(rho)
    n = int(np.prod(shape))
    return (C.unpack_radix(rd.get(C.radix_len(n, base)), n, base) - off).reshape(shape)


# ------------------------------------------------------------------ statements

class TopStatement:
    """K constraints <phi_k, w> = b_k over R_q and L constraints Theta(<phi'_l, w>) = b'_l over F"""

    def __init__(self, phis, bs, fphis, fbs):
        self.phis, self.bs, self.fphis, self.fbs = phis, bs, fphis, fbs
        self.K, self.L, self.m = len(phis), len(fphis), (phis[0] if len(phis) else fphis[0]).shape[0]

    def agg(self, alpha, mu):
        Phi = np.zeros((self.m, D), dtype=np.uint64)
        terms = []
        for k in range(self.K):
            Phi = (Phi + R.mul(alpha[k][None, :], self.phis[k])) % QU
            terms.append(R.theta(R.mul(alpha[k], self.bs[k])))
        for l in range(self.L):
            Phi = (Phi + R.fmul_ring(mu[l][None, :], self.fphis[l])) % QU
            terms.append(R.fmul(mu[l], self.fbs[l]))
        return Phi, fsum(terms)


class LevelStatement:
    """The relation on the last message of a level: openings of the commitments, the folded evaluation
    claim and (if the partial evaluations are committed) the final sum-check claim."""

    def __init__(self, lv, A, Bm, Cm, c, E, H, H1, H2, Wclear, PhiP, claim):
        self.lv, self.A, self.Bm, self.Cm, self.c, self.E = lv, A, Bm, Cm, c, E
        self.H, self.H1, self.H2, self.Wclear, self.PhiP, self.claim = H, H1, H2, Wclear, PhiP, claim
        mode = lv['wmode']
        self.commit, self.opn = (mode == 'commit'), (mode == 'open')
        self.K = lv['kap'] + 1 + (0 if self.opn else lv['k1']) + (lv['k1'] if self.commit else 0)
        self.L = 1 if self.commit else 0
        self.m = lv['m2']

    def agg(self, alpha, mu):
        lv = self.lv
        n, r, kap, k1, t1, tw = lv['n'], lv['r'], lv['kap'], lv['k1'], lv['t1'], lv['tw']
        b, di, do, b1, bw = lv['b'], lv['di'], lv['do'], lv['b1'], lv['bw']
        dec = lv['dec']
        aA, aE = alpha[:kap], alpha[kap]
        pos = kap + 1
        terms = []
        segs = []
        zc = F.lincomb_uniform(aA, self.A)
        zc = (zc + R.fmul_ring(self.E, aE[None, :])) % QU
        segs.append(zc)
        if dec:
            segs.append((zc * np.uint64(b)) % QU)
        cq = R.to_q(self.c)
        if not self.opn:
            be = alpha[pos:pos + k1]
            pos += k1
            ac = R.mul(aA[None, :, :], cq[:, None, :])                       # (r, kap, 64)
            seg = np.zeros((r, t1, kap, D), dtype=np.uint64)
            for j in range(t1):
                f = np.uint64((-(1 << di) * pow(b1, j, Q)) % Q)
                seg[:, j] = (ac * f) % QU
            seg = (seg.reshape(-1, D) + F.lincomb_uniform(be, self.Bm)) % QU
            segs.append(seg)
            terms += [R.theta(R.mul(be[k], R.to_q(self.H1[k] << do))) for k in range(k1)]
        else:
            rhs = (R.mul(cq[:, None, :], R.to_q(self.H << di)).sum(axis=0)) % QU     # (kap, 64)
            terms += [R.theta(R.mul(aA[k], rhs[k])) for k in range(kap)]
        if self.commit:
            ga = alpha[pos:pos + k1]
            pos += k1
            ec = R.mul(aE[None, :], cq)                                      # (r, 64)
            mp = R.fmul_ring(mu[0][None, :], self.PhiP)                      # (r, 64)
            seg = np.zeros((r, tw, D), dtype=np.uint64)
            for j in range(tw):
                f = np.uint64(pow(bw, j, Q))
                seg[:, j] = (((QU - ec) % QU + mp) % QU * f) % QU
            seg = (seg.reshape(-1, D) + F.lincomb_uniform(ga, self.Cm)) % QU
            segs.append(seg)
            terms += [R.theta(R.mul(ga[k], R.to_q(self.H2[k] << do))) for k in range(k1)]
            terms.append(R.fmul(mu[0], self.claim))
        else:
            cw = R.mul(cq, self.Wclear).sum(axis=0) % QU
            terms.append(R.theta(R.mul(aE, cw)))
        if di > 0:
            neg = (QU - aA) % QU
            segs.append(neg)
            if dec:
                segs.append((neg * np.uint64(b)) % QU)
        if do > 0 and not self.opn:
            segs.append((QU - be) % QU)
            if self.commit:
                segs.append((QU - ga) % QU)
        Phi = np.concatenate(segs, axis=0)
        assert Phi.shape[0] == self.m, (Phi.shape, self.m)
        return Phi, fsum(terms)


# ------------------------------------------------------------------ steps shared by prover and verifier

def jl_chunks(seed, ncols):
    rng = np.random.Generator(np.random.PCG64(int.from_bytes(seed, "little")))
    for start in range(0, ncols, CHUNK):
        w = min(CHUNK, ncols - start)
        yield start, R.ternary(rng, (JL, w), np.float64)


def jl_project(seed, wflat):
    p = np.zeros(JL, dtype=np.float64)
    for start, Pi in jl_chunks(seed, wflat.shape[0]):
        p += Pi @ wflat[start:start + Pi.shape[1]]
    return np.rint(p).astype(np.int64)


def jl_aggregate(seed, omega, mpad):
    """Psi_k = sum_j omega_kj sigma_{-1}(pi_j) as elements of R_q^mpad, k = 1..K2"""
    om = omega.astype(np.float64)
    out = np.zeros((omega.shape[0], mpad * D), dtype=np.uint64)
    for start, Pi in jl_chunks(seed, mpad * D):
        out[:, start:start + Pi.shape[1]] = np.mod(om @ Pi, float(Q)).astype(np.uint64)
    return sigma_m1(out.reshape(omega.shape[0], mpad, D))


def level_keys(idx, lv):
    n, r, kap, k1, t1, tw = lv['n'], lv['r'], lv['kap'], lv.get('k1', 0), lv.get('t1', 0), lv.get('tw', 0)
    tag = b"lvl%d|" % idx
    A = crs(tag + b"A", (kap, n, D))
    mode = lv.get('wmode', 'last')
    Bm = crs(tag + b"B", (k1, r * t1 * kap, D)) if mode in ('commit', 'clear') else None
    Cm = crs(tag + b"C", (k1, r * tw, D)) if mode == 'commit' else None
    return A, Bm, Cm


K2 = 4


def derive_agg(tr, stmt):
    rng = tr.rng(b"agg")
    alpha = R.uniform(rng, (stmt.K, D))
    mu = R.uniform(rng, (stmt.L, 4))
    nu = R.uniform(rng, (K2, 4))
    return alpha, mu, nu


def absorb_f(tr, x):
    tr.absorb(np.asarray(x, dtype="<u8"))


# ------------------------------------------------------------------ prover

STATS = []          # one record per level of the last call of prove(): norms and numbers of attempts


def bind_statement(tr, stmt, sched):
    """Absorb the public input into the transcript before the first message: the parameters of all levels
    and the statement.  Without this the challenges would not depend on the statement, and a prover could
    choose the statement after seeing them."""
    par = []
    for lv in sched:
        par.append(";".join("%s=%r" % (k, lv[k]) for k in sorted(lv) if isinstance(lv[k], (bool, int, float, str))))
    tr.absorb("|".join(par).encode())
    for arr in list(stmt.phis) + list(stmt.bs) + list(stmt.fphis) + list(stmt.fbs):
        a = np.ascontiguousarray(arr)
        tr.absorb(repr(a.shape).encode())
        tr.absorb(a)


def prove(stmt, w, sched, wr, log=print, tr=None):
    """tr: a transcript to which the caller has already bound the statement (used by pcs.py, where the
    statement of the first level is determined by the commitment, the evaluation claim and the messages
    of the first fold); if None, a new transcript is created and bind_statement is applied."""
    del STATS[:]
    if tr is None:
        tr = R.Transcript()
        bind_statement(tr, stmt, sched)
    for idx, lv in enumerate(sched):
        t0 = time.time()
        stmt, w = prove_level(idx, lv, stmt, w, tr, wr, log)
        log("  prover level %d done in %.1fs" % (idx + 1, time.time() - t0))
    return wr


def prove_level(idx, lv, stmt, w, tr, wr, log):
    last = (lv['kind'] == 'last')
    mode = 'last' if last else lv['wmode']
    n, r, kap, di = lv['n'], lv['r'], lv['kap'], lv['di']
    beta = lv['beta']
    mpad = r * n
    A, Bm, Cm = level_keys(idx, lv)
    S = np.zeros((mpad, D), dtype=np.int64)
    S[:w.shape[0]] = w
    S = S.reshape(r, n, D)
    nw = float((S.astype(np.float64) ** 2).sum())
    log("level %d (%s): m=%d r=%d n=%d  ||w||/beta = %.4f" % (idx + 1, mode, w.shape[0], r, n, sqrt(nw) / beta))

    # 1. inner commitments, rounded
    T = F.Prepared(A).apply(S)                                              # (r, kap, 64)
    H, Lo = round_hi(T, di)
    if mode in ('last', 'open'):
        put_hi(wr, 'inner', H, di)
        tr.absorb(H)
    else:
        t1, b1, k1, do = lv['t1'], lv['b1'], lv['k1'], lv['do']
        That = digits(H, b1, t1).transpose(1, 0, 2, 3)                      # (r, t1, kap, 64)
        U1 = F.Prepared(Bm).apply(That.reshape(1, -1, D))[0]
        H1, L1 = round_hi(U1, do)
        put_hi(wr, 'outer', H1, do)
        tr.absorb(H1)

    # 2. projection (bjl: the norm bound the verifier accepts; equal to beta unless the schedule says otherwise)
    bjl = lv.get('bjl', beta)
    wflat = S.reshape(-1).astype(np.float64)
    nonce = 0
    while True:
        tj = tr.h.copy()
        tj.update(b"jl" + bytes([nonce]))
        seed = tj.digest(32)
        p = jl_project(seed, wflat)
        if float((p.astype(np.float64) ** 2).sum()) <= 128 * bjl * bjl:
            break
        nonce += 1
        if nonce > 250:
            raise RuntimeError("projection keeps failing; witness too long")
    stat = dict(level=idx + 1, kind=mode, m=int(w.shape[0]), fill=sqrt(nw) / beta, jl_tries=nonce + 1,
                p_ratio=sqrt(float((p.astype(np.float64) ** 2).sum()) / (128 * bjl * bjl)))
    STATS.append(stat)
    wr.put('nonce', bytes([nonce]))
    wr.put_var('jl', C.encode(p, bjl / sqrt(2)))
    tr.absorb(bytes([nonce]))
    tr.absorb(p)

    # 3. lifting the aggregated projections to F
    omega = R.uniform(tr.rng(b"omega"), (K2, JL))
    Psi = jl_aggregate(seed, omega, mpad)
    Sq = R.to_q(S.reshape(mpad, D))
    Bk = np.stack([theta_dot(Psi[k], Sq) for k in range(K2)])
    for k in range(K2):
        assert int(Bk[k, 0]) == int((omega[k].astype(object) * p.astype(object)).sum() % Q)
    wr.put('lift', np.ascontiguousarray(Bk[:, 1:].astype("<u4")).tobytes())
    absorb_f(tr, Bk[:, 1:])

    # 4. aggregation into one constraint over F
    alpha, mu, nu = derive_agg(tr, stmt)
    Phi, bF = stmt.agg(alpha, mu)
    Phi = np.concatenate([Phi, np.zeros((mpad - Phi.shape[0], D), dtype=np.uint64)])
    for k in range(K2):
        Phi = (Phi + R.fmul_ring(nu[k][None, :], Psi[k])) % QU
    bF = fsum([bF] + [R.fmul(nu[k], Bk[k]) for k in range(K2)])
    del Psi
    assert (theta_dot(Phi, Sq) == bF).all(), "aggregated constraint does not hold"

    # 5. sum-check over the position inside a block
    mu_k = lv['mu']
    npad = 1 << mu_k
    Pt = np.zeros((r, npad, D), dtype=np.uint64)
    St = np.zeros((r, npad, D), dtype=np.uint64)
    Pt[:, :n] = Phi.reshape(r, n, D)
    St[:, :n] = Sq.reshape(r, n, D)
    claim = bF
    rho = []
    two = R.fconst(2)
    for j in range(mu_k):
        P0, P1 = Pt[:, 0::2].reshape(-1, D), Pt[:, 1::2].reshape(-1, D)
        S0, S1 = St[:, 0::2].reshape(-1, D), St[:, 1::2].reshape(-1, D)
        g_at0 = theta_dot(P0, S0)
        g_at1 = R.fsub(claim, g_at0)
        PX = (2 * P1 + QU - P0) % QU
        SX = (2 * S1 + QU - S0) % QU
        g_at2 = theta_dot(PX, SX)
        # g(t) = g0 + g1 t + g2 t^2
        g0 = g_at0
        g2 = (((g_at2.astype(object) - 2 * g_at1.astype(object) + g_at0.astype(object)) * R.HALF) % Q).astype(np.uint64)
        g1 = R.fsub(R.fsub(g_at1, g0), g2)
        wr.put('sumcheck', np.ascontiguousarray(np.stack([g0, g2]).astype("<u4")).tobytes())
        absorb_f(tr, np.stack([g0, g2]))
        rj = R.uniform(tr.rng(b"rho"), (4,))
        rho.append(rj)
        claim = fsum([g0, R.fmul(g1, rj), R.fmul(g2, R.fmul(rj, rj))])
        Pt = ((P0 + R.fmul_ring(rj[None, :], (P1 + QU - P0) % QU)) % QU).reshape(r, -1, D)
        St = ((S0 + R.fmul_ring(rj[None, :], (S1 + QU - S0) % QU)) % QU).reshape(r, -1, D)
    PhiP = Pt[:, 0]
    Wv = St[:, 0]                                                           # partial evaluations
    assert (theta_dot(PhiP, Wv) == claim).all()
    E = R.eq_table(np.stack(rho))[:n] if mu_k else R.FONE[None, :].copy()

    # 6. partial evaluations: committed or in the clear
    H2 = L2 = What = None
    if mode == 'commit':
        tw, bw = lv['tw'], lv['bw']
        What = digits(R.center(Wv), bw, tw).transpose(1, 0, 2)              # (r, tw, 64)
        U2 = F.Prepared(Cm).apply(What.reshape(1, -1, D))[0]
        H2, L2 = round_hi(U2, lv['do'])
        put_hi(wr, 'outer', H2, lv['do'])
        tr.absorb(H2)
    else:
        put_ring(wr, 'evals', Wv)
        tr.absorb(Wv)

    # 7. amortisation
    gam2, game2 = lv['gam2'], lv['game2']
    nonce = 0
    while True:
        tc = tr.h.copy()
        tc.update(b"c" + bytes([nonce]))
        c = R.challenges(np.random.Generator(np.random.PCG64(int.from_bytes(tc.digest(32), "little"))), r)
        z = F.lincomb_small(c, S)                                           # (n, 64)
        e = F.lincomb_small(c, Lo)                                          # (kap, 64)
        zz = float((z.astype(np.float64) ** 2).sum())
        ee = float((e.astype(np.float64) ** 2).sum())
        if 'qmax' in lv:                      # joint bound on (z, e), at every level
            stat['q_ratio'] = sqrt((zz + ee) / lv['qmax'])
            if zz + ee <= lv['qmax'] * (1 + 1e-12):
                break
        elif not last:
            stat['z_ratio'] = sqrt(zz / gam2)
            break
        else:
            stat['z_ratio'] = sqrt(zz / gam2)
            if zz <= gam2 and ee <= game2 + 1e-9:
                break
        nonce += 1
        if nonce > 250:
            raise RuntimeError("amortised opening keeps failing")
    stat['c_tries'] = nonce + 1
    wr.put('nonce', bytes([nonce]))
    tr.absorb(bytes([nonce]))
    tr.absorb(c)
    if last:
        wr.put_var('z', C.encode(z, sqrt(gam2 / (n * D))))
        return None, None

    # 8. next witness and next statement
    parts = []
    b = lv['b']
    if lv['dec']:
        zd = digits(z, b, 2)
        parts += [zd[0], zd[1]]
    else:
        parts.append(z)
    if mode != 'open':
        parts.append(That.reshape(-1, D))
    if mode == 'commit':
        parts.append(What.reshape(-1, D))
    if di > 0:
        if lv['dec']:
            ed = digits(e, b, 2)
            parts += [ed[0], ed[1]]
        else:
            parts.append(e)
    if lv['do'] > 0 and mode != 'open':
        parts.append(L1)
        if mode == 'commit':
            parts.append(L2)
    w2 = np.concatenate(parts, axis=0)
    st2 = LevelStatement(lv, A, Bm, Cm, c, E, H, None if mode == 'open' else H1, H2,
                         None if mode == 'commit' else Wv, PhiP, claim)
    return st2, w2


# ------------------------------------------------------------------ verifier

def verify(stmt, sched, proof, log=print, tr=None):
    """tr: see prove()."""
    rd = Reader(proof)
    if tr is None:
        tr = R.Transcript()
        bind_statement(tr, stmt, sched)
    try:
        for idx, lv in enumerate(sched):
            t0 = time.time()
            stmt = verify_level(idx, lv, stmt, tr, rd)
            log("  verifier level %d done in %.1fs" % (idx + 1, time.time() - t0))
    except (ValueError, AssertionError) as ex:
        log("  reject: %s" % ex)
        return False
    return rd.pos == len(proof)


def need(cond, msg):
    if not cond:
        raise ValueError(msg)


def verify_level(idx, lv, stmt, tr, rd):
    last = (lv['kind'] == 'last')
    mode = 'last' if last else lv['wmode']
    n, r, kap, di = lv['n'], lv['r'], lv['kap'], lv['di']
    beta = lv['beta']
    mpad = r * n
    A, Bm, Cm = level_keys(idx, lv)
    H = H1 = H2 = None
    if mode in ('last', 'open'):
        H = get_hi(rd, (r, kap, D), di)
        tr.absorb(H)
    else:
        H1 = get_hi(rd, (lv['k1'], D), lv['do'])
        tr.absorb(H1)
    nonce = rd.get(1)[0]
    tj = tr.h.copy()
    tj.update(b"jl" + bytes([nonce]))
    seed = tj.digest(32)
    pdata = rd.get_var()
    bjl = lv.get('bjl', beta)
    p = C.decode(pdata, JL, bjl / sqrt(2))
    need(C.encode(p, bjl / sqrt(2)) == pdata, "non-canonical projection")
    need(float((p.astype(np.float64) ** 2).sum()) <= 128 * bjl * bjl, "projection too long")
    tr.absorb(bytes([nonce]))
    tr.absorb(p)
    omega = R.uniform(tr.rng(b"omega"), (K2, JL))
    Bk = np.zeros((K2, 4), dtype=np.uint64)
    Bk[:, 1:] = get_ring(rd, (K2, 3))
    for k in range(K2):
        Bk[k, 0] = int((omega[k].astype(object) * p.astype(object)).sum() % Q)
    absorb_f(tr, Bk[:, 1:])
    alpha, mu, nu = derive_agg(tr, stmt)
    Phi, bF = stmt.agg(alpha, mu)
    Phi = np.concatenate([Phi, np.zeros((mpad - Phi.shape[0], D), dtype=np.uint64)])
    Psi = jl_aggregate(seed, omega, mpad)
    for k in range(K2):
        Phi = (Phi + R.fmul_ring(nu[k][None, :], Psi[k])) % QU
    del Psi
    claim = fsum([bF] + [R.fmul(nu[k], Bk[k]) for k in range(K2)])
    mu_k = lv['mu']
    rho = []
    for j in range(mu_k):
        g = get_ring(rd, (2, 4))
        absorb_f(tr, g)
        g0, g2 = g[0], g[1]
        g1 = R.fsub(R.fsub(claim, R.fmul(R.fconst(2), g0)), g2)             # g(0) + g(1) = claim
        rj = R.uniform(tr.rng(b"rho"), (4,))
        rho.append(rj)
        claim = fsum([g0, R.fmul(g1, rj), R.fmul(g2, R.fmul(rj, rj))])
    npad = 1 << mu_k
    eq = R.eq_table(np.stack(rho)) if mu_k else R.FONE[None, :].copy()
    E = eq[:n]
    PhiP = R.fmul_ring(E[None, :, :], Phi.reshape(r, n, D)).sum(axis=1) % QU
    Wv = None
    if mode == 'commit':
        H2 = get_hi(rd, (lv['k1'], D), lv['do'])
        tr.absorb(H2)
    else:
        Wv = get_ring(rd, (r, D))
        tr.absorb(Wv)
        need((theta_dot(PhiP, Wv) == claim).all(), "final sum-check claim")
    nonce = rd.get(1)[0]
    tc = tr.h.copy()
    tc.update(b"c" + bytes([nonce]))
    c = R.challenges(np.random.Generator(np.random.PCG64(int.from_bytes(tc.digest(32), "little"))), r)
    tr.absorb(bytes([nonce]))
    tr.absorb(c)
    if last:
        gam2, game2 = lv['gam2'], lv['game2']
        zdata = rd.get_var()
        z = C.decode(zdata, n * D, sqrt(gam2 / (n * D))).reshape(n, D)
        need(C.encode(z, sqrt(gam2 / (n * D))) == zdata, "non-canonical opening")
        zz = float((z.astype(np.float64) ** 2).sum())
        joint = bool(lv.get('joint'))
        need(zz <= gam2 * (1 + 1e-12 if joint else 1), "amortised opening too long")
        cq = R.to_q(c)
        rhs = R.mul(cq[:, None, :], R.to_q(H << di)).sum(axis=0) % QU       # (kap, 64)
        e = R.center((R.matvec_small(A, z) + QU - rhs) % QU)
        if di == 0:
            need(not e.any(), "opening equation")
        elif joint:
            need(zz + float((e.astype(np.float64) ** 2).sum()) <= gam2 * (1 + 1e-12), "opening too long")
        else:
            need(float((e.astype(np.float64) ** 2).sum()) <= game2 + 1e-9, "low part too long")
        lhs = R.mul(R.emb(E), R.to_q(z)).sum(axis=0) % QU
        need((lhs == R.mul(cq, Wv).sum(axis=0) % QU).all(), "folded evaluation claim")
        return None
    return LevelStatement(lv, A, Bm, Cm, c, E, H, H1, H2, Wv, PhiP, claim)


# ------------------------------------------------------------------ test instances

def random_instance(m, beta, K0=2, L0=1, seed=1, fill=0.985):
    """witness with i.i.d. uniform coefficients of norm about fill*beta and random constraints"""
    rng = np.random.default_rng(seed)
    s = fill * beta / sqrt(m * D)
    half = max(1, int(round(sqrt(3) * s)))
    w = rng.integers(-half, half + 1, size=(m, D), dtype=np.int64)
    wq = R.to_q(w)
    phis = [R.uniform(rng, (m, D)) for _ in range(K0)]
    bs = [R.mul(ph, wq).sum(axis=0) % QU for ph in phis]
    fphis = [R.uniform(rng, (m, D)) for _ in range(L0)]
    fbs = [theta_dot(ph, wq) for ph in fphis]
    return TopStatement(phis, bs, fphis, fbs), w


def model_bits(sched):
    out = {}
    for lv in sched:
        if lv['kind'] == 'rec':
            if lv['wmode'] == 'open':
                out['inner'] = out.get('inner', 0) + lv['outer']
            else:
                out['outer'] = out.get('outer', 0) + lv['outer'] + (lv['garb'] if lv['wmode'] == 'commit' else 0)
            if lv['wmode'] != 'commit':
                out['evals'] = out.get('evals', 0) + lv['garb']
        else:
            out['inner'] = out.get('inner', 0) + lv['t']
            out['evals'] = out.get('evals', 0) + lv['garb']
            out['z'] = lv['z']
        for a, bkey in (('jl', 'jl'), ('lift', 'lift'), ('sumcheck', 'sc'), ('chal', 'chal')):
            out[a] = out.get(a, 0) + lv[bkey]
    return {k: v / 8 for k, v in out.items()}
