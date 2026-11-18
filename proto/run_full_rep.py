#!/usr/bin/env python3
"""Repeated full-scale runs of the tail (as run_full.py, uniform witnesses) with different seeds, for the
distribution of the numbers of attempts of the honest prover.
usage: run_full_rep.py <key in ../model/final.json> <output json> <runs>      (resumable)"""
import json, os, sys, time
from math import sqrt, log2
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import params as P, wmodel as W
import tail as T, ring as R

key, outp, runs = sys.argv[1], sys.argv[2], int(sys.argv[3])
ent = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model", "final.json")))[key]
W.set_bkz(int(key.split("_")[3]))
sched = ent['trace']
m, beta = int(ent['W0']), float(ent['beta0'])
try:
    out = json.load(open(outp))
except Exception:
    out = dict(key=key, m=m, logbeta=log2(beta), levels=len(sched), model_total=ent['tail'] * 1024,
               seeds=64 * len(sched), runs=[])
done = {r['seed'] for r in out['runs']}
for seed in range(runs):
    if seed in done:
        continue
    rng = np.random.default_rng(3000 + seed)
    half = int(sqrt(3) * beta / sqrt(m * R.D)) + 2
    while True:
        w = rng.integers(-half, half + 1, size=(m, R.D), dtype=np.int64)
        nw = sqrt(float((w.astype(np.float64) ** 2).sum()))
        if nw <= beta:
            break
        half -= 1
    wq = R.to_q(w)
    phis = [R.uniform(rng, (m, R.D)) for _ in range(2)]
    bs = [R.mul(ph, wq).sum(axis=0) % R.QU for ph in phis]
    fphis = [R.uniform(rng, (m, R.D))]
    fbs = [T.theta_dot(fphis[0], wq)]
    del wq
    stmt = T.TopStatement(phis, bs, fphis, fbs)
    rec = dict(seed=seed, fill=nw / beta)
    t0 = time.time()
    try:
        wr = T.prove(stmt, w, sched, T.Writer(), log=lambda *a: None)
        proof = wr.bytes()
        rec.update(tp=time.time() - t0, size=len(proof))
        t0 = time.time()
        rec.update(ok=bool(T.verify(stmt, sched, proof, log=lambda *a: None)), tv=time.time() - t0)
    except RuntimeError as ex:
        rec.update(ok=False, failed=str(ex))
    rec['stats'] = [dict(s) for s in T.STATS]
    out['runs'].append(rec)
    json.dump(out, open(outp, "w"), indent=1)
    print("seed %d ok=%s size=%s jl=%s c=%s fills=%s" % (seed, rec['ok'], rec.get('size'),
          [s['jl_tries'] for s in rec['stats']], [s.get('c_tries') for s in rec['stats']],
          [round(s['fill'], 3) for s in rec['stats']]), flush=True)
