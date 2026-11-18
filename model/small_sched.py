#!/usr/bin/env python3
"""Schedules of the tail for the small sizes of the end-to-end runs (N = 2^18 ... 2^24), computed exactly
as the entries of final.json (run_final.job), and cached in small_sched.json.
usage: small_sched.py [logN ...]"""
import json, sys, time
import gh, run_final as RF

OUT = "small_sched.json"
try:
    out = json.load(open(OUT))
except Exception:
    out = {}
sizes = [int(a) for a in sys.argv[1:]] or [18, 20, 22, 24]
for lN in sizes:
    for var in ("full", "gc", "P"):
        key = f"whip_gh{lN}_A_483_{var}"
        if key in out:
            continue
        t0 = time.time()
        try:
            res = RF.job(("whip", "gh", lN, "A", 483, var))
        except Exception as ex:                      # e.g. ranks that do not converge
            print(key, "FAILED:", repr(ex), flush=True)
            out[key] = dict(failed=repr(ex))
            continue
        out.update(res)
        r = res[key]
        print(f"{key:24s} {r['total']:7.2f} KB (tail {r['tail']:.2f}, own {r['own']:.2f}) levels={r['levels']} "
              f"W0={r['W0']} n={r.get('n')} n1={r.get('n1')} Dt={r.get('Dt')} Du={r.get('Du')}  {time.time()-t0:.0f}s",
              flush=True)
        json.dump(out, open(OUT, "w"))
