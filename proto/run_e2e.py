#!/usr/bin/env python3
"""End-to-end runs of the polynomial commitment scheme (pcs.py): commit -> open -> verify on actual
polynomials, for several sizes, settings and kinds of polynomials, and negative tests.  -> e2e.json

usage:
  run_e2e.py runs <logN> <var,var,..> <kind,kind,..> <reps> [--workers W] [--threads T] [--out e2e.json]
  run_e2e.py neg  <logN> <var,var,..> [--flips 60] [--workers W] [--threads T] [--out e2e.json]
var: full (proof-based, expected norms), gc (guaranteed completeness), P (intermediate option).
kind: uniform sparse const period7 blockrep small8 maxdig maxdig_eq.
Every run is done in a fresh process: setup and commit and open in one process (BLAS with T threads),
verify in a second process; wall times and the maximal resident set size of both are recorded.
Results are written after every run; runs that are already in the output file are skipped.
"""
import os
import sys

if __name__ == "__main__" or os.environ.get("WHIPPET_E2E_CHILD"):
    th = "1"
    for i, a in enumerate(sys.argv):
        if a == "--threads":
            th = sys.argv[i + 1]
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, th)

import json
import platform
import resource
import subprocess
import tempfile
import time
from math import sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pcs                  # noqa: E402
import ring as R            # noqa: E402
import tail as T            # noqa: E402

Q = R.Q
KINDS = ("uniform", "sparse", "const", "period7", "blockrep", "small8", "maxdig", "maxdig_eq")
SCRATCH = os.environ.get("WHIPPET_SCRATCH", tempfile.gettempdir())


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def extremal_values(pp):
    """the two coefficients whose centred base-b0 digits have the largest squared norm: all lower digits
    equal to -b0/2 and the top digit as large (resp. as small) as the range (-q/2, q/2] allows"""
    b0, d0 = pp.b0, pp.d0
    low = -(b0 // 2) * (b0 ** (d0 - 1) - 1) // (b0 - 1)
    top = b0 ** (d0 - 1)
    vals = []
    for t in (((Q // 2) - low) // top, -((Q // 2 + low) // top)):
        v = low + t * top
        assert -(Q // 2) < v <= Q // 2
        dg = T.digits(np.array([v]), b0, d0)[:, 0]
        assert T.recompose(dg[:, None], b0)[0] == v and all(dg[:-1] == -(b0 // 2))
        vals.append((v, int((dg.astype(np.int64) ** 2).sum())))
    vals.sort(key=lambda a: -a[1])
    return vals


def polynomial(kind, pp, rng):
    """N coefficients in [0, q) as uint32"""
    N = pp.N
    u32 = np.uint32
    if kind == "uniform":
        return rng.integers(0, Q, size=N, dtype=u32)
    if kind == "sparse":                    # 1% of the coefficients non-zero
        f = np.zeros(N, dtype=u32)
        idx = rng.choice(N, size=N // 100, replace=False)
        f[idx] = rng.integers(1, Q, size=idx.shape[0], dtype=u32)
        return f
    if kind == "const":                     # all coefficients equal
        return np.full(N, int(rng.integers(1, Q)), dtype=u32)
    if kind == "period7":                   # a random pattern of length 7, repeated
        return np.resize(rng.integers(0, Q, size=7, dtype=u32), N)
    if kind == "blockrep":                  # a random block of m*d coefficients, repeated: all blocks f_i equal
        return np.resize(rng.integers(0, Q, size=pp.m * R.D, dtype=u32), N)
    if kind == "small8":                    # coefficients of 8 bits
        return rng.integers(0, 256, size=N, dtype=u32)
    ev = extremal_values(pp)
    if kind == "maxdig":                    # every coefficient one of the two extremal values
        return np.where(rng.integers(0, 2, size=N, dtype=np.uint8) == 0, u32(ev[0][0] % Q), u32(ev[1][0] % Q))
    if kind == "maxdig_eq":                 # all coefficients equal to the extremal value
        return np.full(N, ev[0][0] % Q, dtype=u32)
    raise ValueError(kind)


# ------------------------------------------------------------------ one run (in a fresh process)

def _verify_child(logN, var, x, y, ufile, pfile):
    env = dict(os.environ, WHIPPET_E2E_CHILD="1")
    cmd = [sys.executable, os.path.abspath(__file__), "--verify", str(logN), var, str(x), str(y), ufile, pfile]
    out = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if out.returncode != 0:
        return dict(ok=False, error=out.stderr[-2000:])
    return json.loads(out.stdout.strip().splitlines()[-1])


def verify_main(args):
    logN, var, x, y, ufile, pfile = int(args[0]), args[1], int(args[2]), int(args[3]), args[4], args[5]
    t0 = time.time()
    pp = pcs.setup(pcs.Params(logN, var))
    ts = time.time() - t0
    u, proof = open(ufile, "rb").read(), open(pfile, "rb").read()
    t0 = time.time()
    ok = pcs.verify(pp, u, x, y, proof, log=lambda *a: None)
    tv = time.time() - t0
    print(json.dumps(dict(ok=bool(ok), t_setup=ts, t_verify=tv, rss_mb=rss_mb())))


def one_run(job):
    logN, var, kind, rep = job
    seed = 1000003 * logN + 7919 * KINDS.index(kind) + rep          # the same polynomials in every setting
    rng = np.random.default_rng(seed)
    rec = dict(logN=logN, var=var, kind=kind, rep=rep, seed=seed)
    t0 = time.time()
    pp = pcs.setup(pcs.Params(logN, var))
    rec['t_setup'] = time.time() - t0
    rec['rss_setup_mb'] = rss_mb()
    f = polynomial(kind, pp, rng)
    x = int(rng.integers(0, Q))
    rec['x'] = x
    t0 = time.time()
    u, st = pcs.commit(pp, f)
    rec['t_commit'] = time.time() - t0
    rec['rss_commit_mb'] = rss_mb()
    rec['commitment_bytes'] = len(u)
    t0 = time.time()
    try:
        y, proof, wr = pcs.open_(pp, st, x, log=lambda *a: None)
    except RuntimeError as ex:
        rec['t_open'] = time.time() - t0
        rec['rss_open_mb'] = rss_mb()
        msg = str(ex)
        if 'first fold' in msg:
            lvl, what = 0, 'first fold challenges'
        elif 'projection' in msg:
            lvl, what = len(T.STATS) + 1, 'projection'
        else:
            lvl, what = len(T.STATS), 'challenges'
        rec.update(ok=False, failed=msg, fail_level=lvl, fail_what=what,
                   ff=dict(pcs.FF_STATS), stats=[dict(s) for s in T.STATS])
        return rec
    rec['t_open'] = time.time() - t0
    rec['rss_open_mb'] = rss_mb()
    rec['y_correct'] = bool(y == pcs.eval_poly(st.f32, x))
    rec['proof_bytes'] = len(proof)
    rec['sizes'] = dict(wr.sizes)
    rec['ff'] = {k: v for k, v in pcs.FF_STATS.items()}
    rec['s_norm2'] = pcs.FF_STATS.get('s_norm2')
    rec['stats'] = [dict(s) for s in T.STATS]
    tag = "%d_%s_%s_%d_%d" % (logN, var, kind, rep, os.getpid())
    ufile, pfile = os.path.join(SCRATCH, "u_" + tag), os.path.join(SCRATCH, "p_" + tag)
    open(ufile, "wb").write(u)
    open(pfile, "wb").write(proof)
    v = _verify_child(logN, var, x, y, ufile, pfile)
    os.remove(ufile)
    os.remove(pfile)
    rec['ok'] = bool(v.get('ok')) and rec['y_correct']
    rec['verify'] = v
    return rec


# ------------------------------------------------------------------ negative tests (in a fresh process)

def neg_run(job):
    logN, var, part, flips, nparts = job
    rng = np.random.default_rng(424242 + 31 * logN + part)
    pp = pcs.setup(pcs.Params(logN, var))
    quiet = lambda *a: None
    f = rng.integers(0, Q, size=pp.N, dtype=np.int64)
    u, st = pcs.commit(pp, f)
    x = int(rng.integers(0, Q))
    y, proof, wr = pcs.open_(pp, st, x, log=quiet)
    res = dict(logN=logN, var=var, part=part, honest=bool(pcs.verify(pp, u, x, y, proof, log=quiet)))
    ffb = wr.sizes['first_fold']
    nf = flips // nparts + (1 if part < flips % nparts else 0)
    rej = 0
    rej_ff = 0
    nff = max(1, nf // 6)
    for t in range(nf):
        pos = int(rng.integers(0, ffb)) if t < nff else int(rng.integers(0, len(proof)))
        bad = bytearray(proof)
        bad[pos] ^= 1 << int(rng.integers(0, 8))
        r = not pcs.verify(pp, u, x, y, bytes(bad), log=quiet)
        rej += r
        rej_ff += r if t < nff else 0
    res.update(flips=nf, flips_rejected=rej, flips_in_first_fold=nff, flips_in_first_fold_rejected=rej_ff)
    if part == 0:
        f2 = rng.integers(0, Q, size=pp.N, dtype=np.int64)
        u2, _ = pcs.commit(pp, f2)
        ubad = bytearray(u)
        ubad[len(u) // 2] ^= 8
        # an honest proof for another point, presented for this point
        x3 = (x + 12345) % Q
        y3, proof3, _ = pcs.open_(pp, st, x3, log=quiet)
        res.update(
            wrong_y=not pcs.verify(pp, u, x, (y + 1) % Q, proof, log=quiet),
            noncanonical_y=not pcs.verify(pp, u, x, y + Q, proof, log=quiet),
            wrong_x=not pcs.verify(pp, u, (x + 1) % Q, y, proof, log=quiet),
            other_commitment=not pcs.verify(pp, u2, x, y, proof, log=quiet),
            flipped_commitment=not pcs.verify(pp, bytes(ubad), x, y, proof, log=quiet),
            proof_for_other_point=not pcs.verify(pp, u, x, y, proof3, log=quiet),
            other_point_with_its_value=not pcs.verify(pp, u, x, y3, proof3, log=quiet) if y3 != y else None,
            truncated=not pcs.verify(pp, u, x, y, proof[:-1], log=quiet),
            appended=not pcs.verify(pp, u, x, y, proof + b"\x00", log=quiet),
            empty=not pcs.verify(pp, u, x, y, b"", log=quiet))
        # honest proofs at special points
        edge = {}
        for xe in (0, 1, Q - 1):
            ye, pe, _ = pcs.open_(pp, st, xe, log=quiet)
            edge[str(xe)] = bool(pcs.verify(pp, u, xe, ye, pe, log=quiet)) and ye == pcs.eval_poly(st.f32, xe)
        res['edge_points_ok'] = edge
    return res


# ------------------------------------------------------------------ driver

def meta():
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE).stdout.strip()
    except Exception:
        commit = None
    cpu = ""
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except Exception:
        pass
    mem = None
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem = int(line.split()[1]) / 2 ** 20
    except Exception:
        pass
    return dict(git=commit, cpu=cpu, cores=os.cpu_count(), mem_gb=mem, python=platform.python_version(),
                numpy=np.__version__, blas_threads=os.environ.get("OPENBLAS_NUM_THREADS"), date=time.strftime("%Y-%m-%d %H:%M"))


def params_summary(logN, var):
    pp = pcs.Params(logN, var)
    mb = pp.model_bytes()
    return dict(key=pp.key, N=pp.N, m=pp.m, r=pp.r, b0=pp.b0, d0=pp.d0, b=pp.b, n=pp.n, n1=pp.n1, Dt=pp.Dt,
                Du=pp.Du, kt=pp.kt, kw=pp.kw, W0=pp.W0, beta1=pp.beta1, qmax=pp.qmax, levels=len(pp.sched),
                model_total=mb['total'], model_first_fold=mb['first_fold'], model_tail=mb['tail'],
                model_commitment=mb['commitment'], model_kb=pp.ent['total'], bs2=pp.g.get('bs2'),
                gbar2=pp.g['gbar'] ** 2, extremal=extremal_values(pp))


def load(outp):
    try:
        return json.load(open(outp))
    except Exception:
        return dict(meta=[], params={}, runs={}, neg={})


def main():
    args = [a for a in sys.argv[1:]]
    if args[0] == "--verify":
        return verify_main(args[1:])
    opt = {}
    pos = []
    i = 0
    while i < len(args):
        if args[i].startswith("--"):
            opt[args[i][2:]] = args[i + 1]
            i += 2
        else:
            pos.append(args[i])
            i += 1
    outp = opt.get("out", os.path.join(HERE, "e2e.json"))
    workers = int(opt.get("workers", "4"))
    db = load(outp)
    db['meta'].append(meta())
    mode, logN, vars_ = pos[0], int(pos[1]), pos[2].split(",")
    for var in vars_:
        db['params']["%d_%s" % (logN, var)] = params_summary(logN, var)
    json.dump(db, open(outp, "w"), indent=1)
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    if mode == "runs":
        kinds = KINDS if pos[3] == "all" else pos[3].split(",")
        reps = int(pos[4])
        jobs = []
        for var in vars_:
            for kind in kinds:
                key = "%d_%s_%s" % (logN, var, kind)
                done = {r['rep'] for r in db['runs'].get(key, [])}
                jobs += [(logN, var, kind, rep) for rep in range(reps) if rep not in done]
        fn = one_run
    else:
        flips = int(opt.get("flips", "60"))
        nparts = workers
        jobs = []
        for var in vars_:
            key = "%d_%s" % (logN, var)
            done = {r['part'] for r in db['neg'].get(key, [])}
            jobs += [(logN, var, p, flips, nparts) for p in range(nparts) if p not in done]
        fn = neg_run
    print("%d jobs, %d workers" % (len(jobs), workers), flush=True)
    os.environ["WHIPPET_E2E_CHILD"] = "1"
    with ctx.Pool(workers, maxtasksperchild=1) as pool:
        for rec in pool.imap_unordered(fn, jobs):
            if mode == "runs":
                key = "%d_%s_%s" % (rec['logN'], rec['var'], rec['kind'])
                db['runs'].setdefault(key, []).append(rec)
                st = rec.get('stats', [])
                print("%-26s rep %2d ok=%s bytes=%s ff_fill=%.4f tries ff=%s jl=%s c=%s  t=%.0f/%.0f/%.0fs %s" % (
                    key, rec['rep'], rec['ok'], rec.get('proof_bytes'), rec.get('ff', {}).get('fill', float('nan')),
                    rec.get('ff', {}).get('c_tries'), [s['jl_tries'] for s in st], [s.get('c_tries') for s in st],
                    rec['t_commit'], rec['t_open'], rec.get('verify', {}).get('t_verify', float('nan')),
                    rec.get('failed', '')), flush=True)
            else:
                key = "%d_%s" % (rec['logN'], rec['var'])
                db['neg'].setdefault(key, []).append(rec)
                print(key, json.dumps(rec), flush=True)
            json.dump(db, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
