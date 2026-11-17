#!/usr/bin/env python3
"""All model numbers of the paper; writes final.json."""
import json, sys
from math import sqrt, log2
from multiprocessing import Pool
import params as P
import wmodel as W
import gh

WVAR = {"full": {}, "plain": dict(E=False, O=False, I=False, inner=False),
        "code": dict(O=False, I=False, inner=False), "round": dict(E=False),
        "noinner": dict(inner=False), "f8": dict(fF=8), "P": dict(Pm=True), "gc": dict(gc=True)}
BVAR = {"base": {}, "TL": dict(TL=True), "all": dict(E=True, O=True, I=True, G=True),
        "lin": dict(LIN=True), "linall": dict(LIN=True, E=True, O=True, I=True, G=True)}


def lin_instance(lw):
    return 2**lw, sqrt(2**lw * P.D) * 2.0**5 / sqrt(12)


def job(arg):
    fam, kind, size, acc, bkz, var = arg
    W.set_bkz(bkz)
    P._LC.clear(); P._LA.clear(); P._SP.clear()
    key = f"{fam}_{kind}{size}_{acc}_{bkz}_{var}"
    extra = {}
    if fam == "whip":
        o = W.WOpt(acc=acc, **WVAR[var])
        if kind == "gh" and o.gc:
            cands = [gh.instance_gc(size, o.theta, o.tmar, lb0=l0, lb=l) for l0, l in ((4, 6), (3, 6))]
            cands = [g for g in cands if g is not None]
        elif kind == "gh":
            cands = [gh.instance(size, o.O, acc)]
        else:
            cands = [None]
        best = None
        for g in cands:
            if g is not None:
                W0, b0, own = g['W'], g['beta'], gh.own_bits(g, True) / P.KB
            else:
                (W0, b0), own = lin_instance(size), 0.0
                if o.gc:                      # the stated norm is the guaranteed bound of the witness
                    pass
            for eta in ((0.1, 0.3) if o.O else (0.1,)):
                o.eta = eta
                res = W.search(W0, b0, o, beam=160)
                if res is not None and (best is None or res[0] + own * P.KB < best[0] + best[5] * P.KB):
                    best = (res[0], res[1], eta, W0, b0, own, g)
        tot, tr, eta, W0, b0, own, g = best
        if g is not None:
            extra = {k: g[k] for k in ('n', 'n1', 'Dt', 'Du', 'com', 'W', 'lb0', 'lb', 'd0')}
        comp = W.comps(tr)
    else:
        fl = BVAR[var]
        rounded = bool(fl.get('O'))
        if kind == "gh":
            g = gh.instance(size, rounded, acc)
            W0, b0, own = g['W'], g['beta'], gh.own_bits(g, False) / P.KB
            extra = {k: g[k] for k in ('n', 'n1', 'Dt', 'Du', 'com', 'W')}
        else:
            (W0, b0), own = lin_instance(size), 0.0
        tot, tr, o = P.tune(W0, b0, thorough=True, acc=acc, **fl)
        comp = P.components(tr)
        eta = o.eta
    return {key: dict(tail=tot / P.KB, own=own, total=tot / P.KB + own, levels=len(tr), comp=comp,
                      trace=tr, eta=eta, W0=W0, beta0=b0, **extra)}


if __name__ == "__main__":
    jobs = []
    for bkz in (438, 483):
        for acc in "CA":
            for lN in (30, 26, 28):
                jobs.append(("whip", "gh", lN, acc, bkz, "full"))
                jobs.append(("base", "gh", lN, acc, bkz, "base"))
    for acc in "CA":
        for var in ("plain", "code", "round", "noinner", "f8", "P"):
            jobs.append(("whip", "gh", 30, acc, 438, var))
        for var in ("TL", "all"):
            jobs.append(("base", "gh", 30, acc, 438, var))
        for lw in (12, 14, 16, 18, 20, 22):
            jobs.append(("whip", "lin", lw, acc, 438, "full"))
            jobs.append(("base", "lin", lw, acc, 438, "base"))
    for bkz in (438, 483):
        for acc in "CA":
            for lN in (30, 26, 28):
                jobs.append(("base", "gh", lN, acc, bkz, "lin"))
    for acc in "CA":
        jobs.append(("base", "gh", 30, acc, 438, "linall"))
        for lw in (12, 14, 16, 18, 20, 22):
            jobs.append(("base", "lin", lw, acc, 438, "lin"))
    # revision: guaranteed-completeness setting, and the ablation at block size 483
    for lN in (30, 26, 28):
        jobs.append(("whip", "gh", lN, "A", 483, "gc"))
        jobs.append(("base", "gh", lN, "A", 483, "linall"))
    jobs.append(("whip", "gh", 30, "A", 438, "gc"))
    for var in ("plain", "code", "round", "noinner", "f8", "P"):
        jobs.append(("whip", "gh", 30, "A", 483, var))
    for lw in (12, 16, 20):
        jobs.append(("whip", "lin", lw, "A", 438, "gc"))
    for acc in "CA":
        for lN in (26, 28):
            jobs.append(("base", "gh", lN, acc, 438, "linall"))
    try:
        out = json.load(open("final.json"))
    except Exception:
        out = {}
    jobs = [j for j in jobs if f"{j[0]}_{j[1]}{j[2]}_{j[3]}_{j[4]}_{j[5]}" not in out]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 2, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(job, jobs):
            out.update(res)
            for k, r in res.items():
                print(f"{k:28s} {r['total']:7.2f} KB (tail {r['tail']:.2f}, own {r['own']:.2f}) levels={r['levels']}", flush=True)
            json.dump(out, open("final.json", "w"))
