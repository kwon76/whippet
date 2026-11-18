"""Polynomial commitment scheme: the commitment and the first fold of Section 5 (Figure 4, equations (10)-(12)),
followed by the linear tail of tail.py.

    pp = setup(Params(logN, var))          var: 'full' (proof-based ranks, expected norms), 'gc' (guaranteed
                                                completeness), 'P' (intermediate option of Section 6.6)
    u, st = commit(pp, f)                  f: N coefficients in [0, q)
    y, proof = open_(pp, st, x)
    ok = verify(pp, u, x, y, proof)

The relation handed to the tail, (12), is the relation of tail.LevelStatement with
    w_hat <-> h_hat,  e_u <-> d_1,  e_v <-> d_2,  a_bar <-> E,  psi_x * b_i <-> Phi'_i,  y_hat <-> claim,
and the witness order [z0 | z1 | t_hat | w_hat | e0 | e1 | e_u | e_v] of LevelStatement.agg.
All parameters (shape, ranks, bases, dropped bits, norm bounds, schedule of the tail) are taken from the
size model: model/final.json (N = 2^26 ... 2^30) or model/small_sched.json (N = 2^18 ... 2^24).
"""
import json
import os
import sys
import time
from math import sqrt, floor, ceil

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "model"))
import params as P          # noqa: E402
import wmodel as W          # noqa: E402
import gh                   # noqa: E402
import ring as R            # noqa: E402
import tail as T            # noqa: E402
import fastring as F        # noqa: E402

Q, QU, D = R.Q, R.QU, R.D
NONCE_MAX = 250
STORE_LIMIT = 6 << 30          # digits of the polynomial are kept in memory up to this size, else recomputed
CHUNK_BYTES = 1 << 28


def load_entry(logN, var):
    key = "whip_gh%d_A_483_%s" % (logN, var)
    for fn in ("final.json", "small_sched.json"):
        path = os.path.join(HERE, "..", "model", fn)
        if os.path.exists(path):
            d = json.load(open(path))
            if key in d and 'trace' in d[key]:
                return key, d[key]
    raise KeyError("no schedule for %s (run model/small_sched.py)" % key)


def _widths(bits, lb, table_rule):
    """digit widths of a value with `bits` bits in base 2^lb, as in gh.first_fold / gh.first_fold_gc"""
    if table_rule:
        k = max(1, floor(bits / lb))
        if bits - (k - 1) * lb > lb + 2:
            k += 1
    else:
        k = max(1, ceil(bits / lb))
    return [lb] * (k - 1) + [bits - (k - 1) * lb]


class Params:
    """public parameters of one size and setting, from the size model"""

    def __init__(self, logN, var):
        self.logN, self.var = logN, var
        self.key, ent = load_entry(logN, var)
        W.set_bkz(483)
        if var == 'gc':
            g = gh.instance_gc(logN, 2.0, 128 ** -0.5, lb0=ent['lb0'], lb=ent['lb'])
        else:
            g = gh.instance(logN, True, 'A')
        for k in ('n', 'n1', 'Dt', 'Du', 'W'):
            assert g[k] == ent[k], (k, g[k], ent[k])
        self.g, self.ent = g, ent
        self.sched = ent['trace']
        self.N = 1 << logN
        self.m, self.r, self.lb0, self.d0, self.lb = g['m'], g['r'], g['lb0'], g['d0'], g['lb']
        self.b0, self.b = 1 << self.lb0, 1 << self.lb
        self.n, self.n1, self.Dt, self.Du = g['n'], g['n1'], g['Dt'], g['Du']
        assert self.Dt > 0 and self.Du > 0
        self.kt = g['kt']
        self.kw = len(_widths(int(P.LOGQ), self.lb, var != 'gc'))
        self.nb = self.d0 * self.m                       # ring elements of s_i
        self.W0 = 2 * self.nb + self.r * self.kt * self.n + self.r * self.kw + 2 * self.n + 2 * self.n1
        assert self.W0 == g['W'] == int(ent['W0']) == self.sched[0]['m'], (self.W0, g['W'])
        assert self.m * self.r * D >= self.N
        self.beta1 = float(ent['beta0'])                 # norm bound of the relation handed to the tail
        self.qmax = g.get('qmax')                        # joint bound on (z, e) (guaranteed setting only)
        self.lv = dict(n=self.nb, r=self.r, kap=self.n, k1=self.n1, t1=self.kt, tw=self.kw, b=self.b,
                       di=self.Dt, do=self.Du, b1=self.b, bw=self.b, dec=True, wmode='commit', m2=self.W0)
        top = (Q // 2) // self.b0 ** (self.d0 - 1) + 1
        self.smax = max(top, self.b0 // 2) + 1           # bound (with margin) on the absolute values of the digits
        self.dtype = np.int8 if self.smax < 127 else np.int16

    def describe(self):
        return ("pcs;N=%d;m=%d;r=%d;b0=%d;d0=%d;b=%d;n=%d;n1=%d;Dt=%d;Du=%d;kt=%d;kw=%d;var=%s"
                % (self.N, self.m, self.r, self.b0, self.d0, self.b, self.n, self.n1, self.Dt, self.Du,
                   self.kt, self.kw, self.var))

    def model_bytes(self):
        """model sizes in bytes without the challenge seeds: first fold, tail, total, commitment"""
        ff = (self.n1 * D * P.hibits(self.Du) + 3 * P.LOGQ) / 8
        tail = self.ent['tail'] * 1024 - 64 * len(self.sched)
        return dict(first_fold=ff, tail=tail, total=ff + tail, commitment=self.g['com'],
                    tail_parts=T.model_bits(self.sched))

    def model_parts(self):
        """expected squared norms of the parts of the first-fold witness in the model (optimised setting),
        resp. their guaranteed bounds (guaranteed setting), for the diagnostics of the experiments"""
        g, Dd = self.g, D
        b, b0, r, m, n, n1 = self.b, self.b0, self.r, self.m, self.n, self.n1
        if self.var != 'gc':
            var = lambda ws: sum(4.0 ** w / 12 for w in ws)
            vz = b0 * b0 / 12 * r * P.TAU
            return dict(z=(b * b / 12 + vz / (b * b)) * self.nb * Dd,
                        t_hat=var(_widths(int(P.LOGQ) - self.Dt, self.lb, False)) * n * r * Dd,
                        w_hat=var(_widths(int(P.LOGQ), self.lb, True)) * r * Dd,
                        e=n * Dd * (b * b / 12 + 4.0 ** self.Dt * P.TAU * r / (P.LOWDIV * b * b)),
                        d=2 * n1 * Dd * 4.0 ** self.Du / P.LOWDIV,
                        ze=g['gbar'] ** 2)
        return dict(ze=self.qmax)


def setup(pp):
    tag = ("pcs|m=%d|r=%d|d0=%d|n=%d|n1=%d|kt=%d|kw=%d|" % (pp.m, pp.r, pp.d0, pp.n, pp.n1, pp.kt, pp.kw)).encode()
    pp.A = T.crs(tag + b"A", (pp.n, pp.nb, D))
    pp.B = T.crs(tag + b"B", (pp.n1, pp.r * pp.kt * pp.n, D))
    pp.Dm = T.crs(tag + b"D", (pp.n1, pp.r * pp.kw, D))
    return pp


def _prep(pp):
    """limb-split copies of the keys for the prover (the verifier does not need them)"""
    if not hasattr(pp, 'Aprep'):
        pp.Aprep, pp.Bprep, pp.Dprep = F.Prepared(pp.A), F.Prepared(pp.B), F.Prepared(pp.Dm)


# ------------------------------------------------------------------ polynomial and digits

def _digits_blocks(pp, f32, i0, i1, out=None):
    """digits s_i (i0 <= i < i1) of the blocks of f: (i1-i0, d0*m, 64); s_i[k*m + j] = k-th digit of f_ij.
    Computed in sub-chunks so that the int64 temporaries stay below CHUNK_BYTES."""
    if out is None:
        out = np.empty((i1 - i0, pp.nb, D), dtype=pp.dtype)
    sub = max(1, CHUNK_BYTES // (pp.m * D * 8 * (pp.d0 + 2)))
    for a in range(i0, i1, sub):
        b = min(i1, a + sub)
        fc = R.center(f32[a * pp.m * D:b * pp.m * D].astype(np.int64)).reshape(b - a, pp.m, D)
        dg = T.digits(fc, pp.b0, pp.d0)                   # (d0, b-a, m, 64)
        out[a - i0:b - i0] = dg.transpose(1, 0, 2, 3).reshape(b - a, pp.nb, D)
    return out


class State:
    def __init__(self, pp, f32):
        self.pp, self.f32 = pp, f32
        self.S = None
        # chunks of blocks of at most 1 GB of digits; all digits are kept if they need at most STORE_LIMIT
        self.rc = max(1, (1 << 30) // (pp.nb * D * np.dtype(pp.dtype).itemsize))
        if pp.r * pp.nb * D * np.dtype(pp.dtype).itemsize <= STORE_LIMIT:
            self.S = np.empty((pp.r, pp.nb, D), dtype=pp.dtype)
            _digits_blocks(pp, f32, 0, pp.r, out=self.S)

    def chunks(self):
        for i0 in range(0, self.pp.r, self.rc):
            i1 = min(self.pp.r, i0 + self.rc)
            yield i0, i1, (self.S[i0:i1] if self.S is not None else _digits_blocks(self.pp, self.f32, i0, i1))


def sqnorm(S):
    """squared Euclidean norm of a (large) integer array, in sub-chunks"""
    S = np.asarray(S).reshape(-1)
    step = 1 << 24
    return float(sum(int((S[a:a + step].astype(np.int64) ** 2).sum()) for a in range(0, S.shape[0], step)))


def commit(pp, f):
    """f: coefficients in [0, q) (at most N).  Returns (u as bytes, state)."""
    f = np.asarray(f)
    assert f.ndim == 1 and f.shape[0] <= pp.N
    _prep(pp)
    f32 = np.zeros(pp.m * pp.r * D, dtype=np.uint32)
    step = 1 << 24
    for a in range(0, f.shape[0], step):
        b = min(f.shape[0], a + step)
        f32[a:b] = (f[a:b].astype(np.int64) % Q).astype(np.uint32)
    st = State(pp, f32)
    Tm = np.empty((pp.r, pp.n, D), dtype=np.uint64)
    for i0, i1, S in st.chunks():
        Tm[i0:i1] = pp.Aprep.apply(S, smax=pp.smax)
    H, L = T.round_hi(Tm, pp.Dt)
    del Tm
    That = T.digits(H, pp.b, pp.kt).transpose(1, 0, 2, 3)                # (r, kt, n, 64)
    U = pp.Bprep.apply(That.reshape(1, -1, D))[0]
    UH, Eu = T.round_hi(U, pp.Du)
    wr = T.Writer()
    T.put_hi(wr, 'u', UH, pp.Du)
    u = wr.bytes()
    st.L, st.That, st.UH, st.Eu, st.u = L, That, UH, Eu, u
    return u, st


def eval_poly(f32, x):
    """f(x) mod q, computed directly from the coefficients (independent check of y)"""
    x = int(x) % Q
    M = f32.shape[0] // D
    xp = np.array([pow(x, k, Q) for k in range(D)], dtype=np.int64)
    xh, xl = xp >> 16, xp & 0xFFFF
    out = 0
    xd = pow(x, D, Q)
    step = 1 << 12
    small = np.array([pow(xd, t, Q) for t in range(step)], dtype=np.uint64)
    big_base = pow(xd, step, Q)
    for c0 in range(0, M, 1 << 18):
        c1 = min(M, c0 + (1 << 18))
        fc = f32[c0 * D:c1 * D].astype(np.int64).reshape(-1, D)
        g = ((fc @ xh) % Q * 65536 + (fc @ xl)) % Q                       # f_t(x), t = c0 .. c1-1
        nb = (c1 - c0 + step - 1) // step
        big = np.array([pow(big_base, (c0 // step) + s, Q) for s in range(nb)], dtype=np.uint64)
        pw = ((big[:, None] * small[None, :]) % QU).reshape(-1)[:c1 - c0]
        assert c0 % step == 0
        out = (out + int((((g.astype(np.uint64) * pw) % QU)).sum(dtype=np.uint64) % QU)) % Q
    return out


def _eval_vectors(pp, x):
    x = int(x) % Q
    xd = pow(x, D, Q)
    a = np.array([pow(xd, j, Q) for j in range(pp.m)], dtype=np.uint64)
    xm = pow(xd, pp.m, Q)
    bvec = np.array([pow(xm, i, Q) for i in range(pp.r)], dtype=np.uint64)
    abar = np.concatenate([(a * np.uint64(pow(pp.b0, k, Q))) % QU for k in range(pp.d0)])
    psi = T.sigma_m1(np.array([pow(x, k, Q) for k in range(D)], dtype=np.uint64)[None, :])[0]
    PhiP = (bvec[:, None] * psi[None, :]) % QU                                # psi_x * b_i, (r, 64)
    E = np.zeros((pp.nb, 4), dtype=np.uint64)
    E[:, 0] = abar
    return abar, bvec, psi, PhiP, E


def bind_public(tr, pp, u, x, y):
    """public parameters, the schedule of all levels, the commitment and the evaluation claim"""
    tr.absorb(pp.describe().encode())
    par = []
    for lv in pp.sched:
        par.append(";".join("%s=%r" % (k, lv[k]) for k in sorted(lv) if isinstance(lv[k], (bool, int, float, str))))
    tr.absorb("|".join(par).encode())
    tr.absorb(bytes(u))
    tr.absorb(int(x).to_bytes(8, "little"))
    tr.absorb(int(y).to_bytes(8, "little"))


def _challenges(tr, nonce, r):
    tc = tr.h.copy()
    tc.update(b"c" + bytes([nonce]))
    return R.challenges(np.random.Generator(np.random.PCG64(int.from_bytes(tc.digest(32), "little"))), r)


FF_STATS = {}           # record of the last call of open_(): first fold norms and attempts, timings


def open_(pp, st, x, log=print):
    """Returns (y, proof bytes, Writer).  Raises RuntimeError if the nonces are exhausted."""
    FF_STATS.clear()
    _prep(pp)
    x = int(x) % Q
    abar, bvec, psi, PhiP, E = _eval_vectors(pp, x)
    t0 = time.time()
    Wv = np.empty((pp.r, D), dtype=np.uint64)
    for i0, i1, S in st.chunks():
        Wv[i0:i1] = F.scalar_dot(abar, S)                                    # w_i = <a_bar, s_i>
    yp = np.zeros(D, dtype=object)
    for i in range(pp.r):
        yp = yp + int(bvec[i]) * Wv[i].astype(object)
    yp = (yp % Q).astype(np.uint64)
    yhat = R.theta(R.mul(psi, yp))
    y = int(yhat[0])
    What = T.digits(R.center(Wv), pp.b, pp.kw).transpose(1, 0, 2)          # (r, kw, 64)
    V = pp.Dprep.apply(What.reshape(1, -1, D))[0]
    VH, Ev = T.round_hi(V, pp.Du)

    tr = R.Transcript(b"whippet-pcs")
    bind_public(tr, pp, st.u, x, y)
    wr = T.Writer()
    T.put_hi(wr, 'first_fold', VH, pp.Du)
    tr.absorb(VH)
    wr.put('first_fold', np.ascontiguousarray(yhat[1:].astype("<u4")).tobytes())
    T.absorb_f(tr, yhat[1:])

    ee_l = st.L
    nonce = 0
    tries = []
    while True:
        c = _challenges(tr, nonce, pp.r)
        z = np.zeros((pp.nb, D), dtype=np.int64)
        for i0, i1, S in st.chunks():
            z += F.lincomb_small(c[i0:i1], S)
        e = F.lincomb_small(c, ee_l)
        zz = float((z.astype(np.float64) ** 2).sum())
        ee = float((e.astype(np.float64) ** 2).sum())
        tries.append(zz + ee)
        if pp.qmax is None or zz + ee <= pp.qmax * (1 + 1e-12):
            break
        nonce += 1
        if nonce > NONCE_MAX:
            FF_STATS.update(c_tries=nonce, failed='first fold')
            raise RuntimeError("first fold: joint bound on (z, e) keeps failing")
    wr.put('first_fold', bytes([nonce]))
    tr.absorb(bytes([nonce]))
    tr.absorb(c)

    zd = T.digits(z, pp.b, 2)
    ed = T.digits(e, pp.b, 2)
    parts = dict(z=[zd[0], zd[1]], t_hat=[st.That.reshape(-1, D)], w_hat=[What.reshape(-1, D)],
                 e=[ed[0], ed[1]], d=[st.Eu, Ev])
    w2 = np.concatenate([np.asarray(a, dtype=np.int64) for k in ('z', 't_hat', 'w_hat', 'e', 'd') for a in parts[k]])
    assert w2.shape[0] == pp.W0
    pn = {k: float(sum((np.asarray(a, dtype=np.float64) ** 2).sum() for a in v)) for k, v in parts.items()}
    pn['ze'] = zz + ee
    nw = sqrt(float((w2.astype(np.float64) ** 2).sum()))
    s2 = sum(sqnorm(S) for _, _, S in st.chunks())
    FF_STATS.update(c_tries=nonce + 1, fill=nw / pp.beta1, z_norm2=zz, e_norm2=ee, s_norm2=s2,
                    parts=pn, model_parts=pp.model_parts(), ze_tries=tries, t_first_fold=time.time() - t0)
    if pp.qmax is not None:
        FF_STATS['q_ratio'] = sqrt((zz + ee) / pp.qmax)
        FF_STATS['bs2'] = pp.g['bs2']
    log("first fold: ||w1||/beta1 = %.4f, ||(z,e)||^2 / model = %.4f, attempts %d"
        % (nw / pp.beta1, (zz + ee) / pp.model_parts()['ze'], nonce + 1))

    stmt = T.LevelStatement(pp.lv, pp.A, pp.B, pp.Dm, c, E, None, st.UH, VH, None, PhiP, yhat)
    T.prove(stmt, w2, pp.sched, wr, log=log, tr=tr)
    return y, wr.bytes(), wr


def verify(pp, u, x, y, proof, log=print):
    try:
        x, y = int(x), int(y)
        T.need(0 <= x < Q and 0 <= y < Q, "non-canonical evaluation claim")
        ru = T.Reader(bytes(u))
        UH = T.get_hi(ru, (pp.n1, D), pp.Du)
        T.need(ru.pos == len(u), "commitment has the wrong length")
        tr = R.Transcript(b"whippet-pcs")
        bind_public(tr, pp, u, x, y)
        rd = T.Reader(proof)
        VH = T.get_hi(rd, (pp.n1, D), pp.Du)
        tr.absorb(VH)
        yhat = np.zeros(4, dtype=np.uint64)
        yhat[0] = y
        yhat[1:] = T.get_ring(rd, (3,))
        T.absorb_f(tr, yhat[1:])
        nonce = rd.get(1)[0]
        c = _challenges(tr, nonce, pp.r)
        tr.absorb(bytes([nonce]))
        tr.absorb(c)
        abar, bvec, psi, PhiP, E = _eval_vectors(pp, x)
        stmt = T.LevelStatement(pp.lv, pp.A, pp.B, pp.Dm, c, E, None, UH, VH, None, PhiP, yhat)
        return T.verify(stmt, pp.sched, proof[rd.pos:], log=log, tr=tr)
    except (ValueError, AssertionError, IndexError, OverflowError) as ex:
        log("  reject: %s" % ex)
        return False
