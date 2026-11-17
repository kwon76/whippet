#!/usr/bin/env python3
"""Full-scale run of the reference implementation on an instance with the witness length and norm of
the relation that the (rounded) first fold of Greyhound produces for N = 2^30.

usage: run_full.py <key in ../model/final.json> <output json> [uniform|const]
The schedule (ranks, bases, dropped bits, norm bounds) is the one selected by the size model."""
import json, os, sys, time
from math import sqrt, log2
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import params as P, wmodel as W
import tail as T, ring as R

key, outp = sys.argv[1], sys.argv[2]
ent = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model", "final.json")))[key]
W.set_bkz(int(key.split("_")[3]))
sched = ent['trace']
m, beta = int(ent['W0']), float(ent['beta0'])
log = lambda *a: (print(*a), sys.stdout.flush())
log("instance %s: m = %d ring elements, beta = 2^%.3f, %d levels, model %.2f KB" % (key, m, log2(beta), len(sched), ent['tail']))

rng = np.random.default_rng(2026)
kind = sys.argv[3] if len(sys.argv) > 3 else "uniform"
if kind == "uniform":
    half = int(sqrt(3) * beta / sqrt(m * R.D)) + 2
    while True:
        w = rng.integers(-half, half + 1, size=(m, R.D), dtype=np.int64)
        nw = sqrt(float((w.astype(np.float64) ** 2).sum()))
        if nw <= beta:
            break
        half -= 1
    log("witness: uniform coefficients in [-%d, %d], ||w|| / beta = %.4f" % (half, half, nw / beta))
else:                                   # all coefficients equal: the most structured witness
    half = int(beta / sqrt(m * R.D))
    w = np.full((m, R.D), half, dtype=np.int64)
    nw = sqrt(float((w.astype(np.float64) ** 2).sum()))
    log("witness: all coefficients equal to %d, ||w|| / beta = %.4f" % (half, nw / beta))
wq = R.to_q(w)
phis = [R.uniform(rng, (m, R.D)) for _ in range(2)]
bs = [R.mul(ph, wq).sum(axis=0) % R.QU for ph in phis]
fphis = [R.uniform(rng, (m, R.D))]
fbs = [T.theta_dot(fphis[0], wq)]
del wq
stmt = T.TopStatement(phis, bs, fphis, fbs)

t0 = time.time()
wr = T.prove(stmt, w, sched, T.Writer(), log=log)
tp = time.time() - t0
proof = wr.bytes()
stats = [dict(s) for s in T.STATS]
log("proof: %d bytes; prover %.0f s" % (len(proof), tp))
t0 = time.time()
ok = T.verify(stmt, sched, proof, log=log)
tv = time.time() - t0
log("verified: %s; verifier %.0f s" % (ok, tv))
# a single flipped bit must be rejected
bad = bytearray(proof); bad[len(bad) // 3] ^= 4
ok_bad = T.verify(stmt, sched, bytes(bad), log=lambda *a: None)
mb = T.model_bits(sched)
res = dict(key=key, m=m, logbeta=log2(beta), fill=nw / beta, half=half, levels=len(sched), ok=bool(ok),
           tampered_rejected=(not ok_bad), size=len(proof), sizes=wr.sizes, model=mb,
           model_total=ent['tail'] * 1024, seeds=64 * len(sched), tp=tp, tv=tv, witness=kind,
           stats=[dict(s) for s in stats])
json.dump(res, open(outp, "w"), indent=1)
log(json.dumps(res))
