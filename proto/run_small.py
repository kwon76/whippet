import sys, time, json
from math import sqrt, log2
import numpy as np
import tail as T
import wmodel as W, params as P

def run(m, logbeta, acc='A', seed=1, quiet=False, opts=None):
    o = W.WOpt(acc=acc, **(opts or {}))
    beta = 2.0**logbeta
    tot, sched = W.search(m, beta, o, beam=60)
    log = (lambda *a: None) if quiet else print
    log("schedule: %d levels, model %.2f KB" % (len(sched), tot / P.KB))
    if not quiet: W.show(sched)
    stmt, w = T.random_instance(m, beta, seed=seed)
    t0 = time.time()
    wr = T.prove(stmt, w, sched, T.Writer(), log=log)
    tp = time.time() - t0
    proof = wr.bytes()
    t0 = time.time()
    ok = T.verify(stmt, sched, proof, log=log)
    tv = time.time() - t0
    mb = T.model_bits(sched)
    return dict(ok=ok, size=len(proof), sizes=wr.sizes, model=mb, model_total=tot / 8, tp=tp, tv=tv), stmt, sched, proof

if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    lb = float(sys.argv[2]) if len(sys.argv) > 2 else 12.5
    acc = sys.argv[3] if len(sys.argv) > 3 else 'A'
    res, stmt, sched, proof = run(m, lb, acc)
    print("verified:", res['ok'], " proof bytes:", res['size'], " model bytes (incl. 64 B/level seeds):", round(res['model_total']))
    for k in sorted(set(res['sizes']) | set(res['model'])):
        print("   %-9s measured %7d   model %9.1f" % (k, res['sizes'].get(k, 0), res['model'].get(k, 0)))
    print("prover %.1fs verifier %.1fs" % (res['tp'], res['tv']))
