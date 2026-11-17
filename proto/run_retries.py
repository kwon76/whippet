#!/usr/bin/env python3
"""Numbers of attempts of the honest prover on small instances, for the optimised schedule (expected
norms) and for the schedule with guaranteed completeness, and for three kinds of witnesses:
  uniform  independent uniform coefficients,
  const    all coefficients equal (all blocks equal, every ring element a multiple of 1 + X + ... + X^63),
  sparse   one coefficient in 64 is non-zero.
All witnesses have norm between 0.96 and 1.0 times the bound of the instance.
usage: run_retries.py <output json> [runs per case]"""
import json, os, sys, time
from math import sqrt, log2
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import params as P, wmodel as W
import tail as T, ring as R

outp = sys.argv[1]
RUNS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
M, LOGB = 3000, 13.0
beta = 2.0 ** LOGB


def witness(kind, rng):
    if kind == "uniform":
        half = int(sqrt(3) * beta / sqrt(M * R.D)) + 2
        while True:
            w = rng.integers(-half, half + 1, size=(M, R.D), dtype=np.int64)
            if sqrt(float((w.astype(np.float64) ** 2).sum())) <= beta:
                return w
            half -= 1
    if kind == "const":
        return np.full((M, R.D), int(beta / sqrt(M * R.D)), dtype=np.int64)
    w = np.zeros((M, R.D), dtype=np.int64)
    idx = rng.integers(0, R.D, size=M)
    w[np.arange(M), idx] = int(beta / sqrt(M)) * rng.choice([-1, 1], size=M)
    return w


def statement(w, rng):
    wq = R.to_q(w)
    phis = [R.uniform(rng, (M, R.D)) for _ in range(2)]
    bs = [R.mul(ph, wq).sum(axis=0) % R.QU for ph in phis]
    fphis = [R.uniform(rng, (M, R.D))]
    fbs = [T.theta_dot(fphis[0], wq)]
    return T.TopStatement(phis, bs, fphis, fbs)


W.set_bkz(483)
out = {}
for mode in ("opt", "gc"):
    o = W.WOpt(acc='A', gc=(mode == "gc"))
    tot, sched = W.search(M, beta, o, beam=60)
    for kind in ("uniform", "const", "sparse"):
        rec = dict(mode=mode, kind=kind, levels=len(sched), model=tot / 8, runs=[])
        for seed in range(RUNS):
            rng = np.random.default_rng(1000 + seed)
            w = witness(kind, rng)
            stmt = statement(w, rng)
            run = dict(seed=seed, fill=sqrt(float((w.astype(np.float64) ** 2).sum())) / beta)
            try:
                wr = T.prove(stmt, w, sched, T.Writer(), log=lambda *a: None)
                proof = wr.bytes()
                run.update(ok=bool(T.verify(stmt, sched, proof, log=lambda *a: None)), size=len(proof),
                           jl=[s['jl_tries'] for s in T.STATS], c=[s['c_tries'] for s in T.STATS],
                           fills=[s['fill'] for s in T.STATS])
            except RuntimeError as ex:
                run.update(ok=False, failed=str(ex), jl=[s['jl_tries'] for s in T.STATS],
                           c=[s.get('c_tries', 0) for s in T.STATS], fills=[s['fill'] for s in T.STATS])
            rec['runs'].append(run)
        out["%s_%s" % (mode, kind)] = rec
        ok = sum(1 for r in rec['runs'] if r['ok'])
        att = [sum(r['jl']) + sum(r['c']) for r in rec['runs']]
        print("%s %-7s: %d/%d proofs, attempts per proof mean %.1f max %d (levels %d)" % (
            mode, kind, ok, RUNS, sum(att) / len(att), max(att), len(sched)))
        sys.stdout.flush()
        json.dump(out, open(outp, "w"), indent=1)
