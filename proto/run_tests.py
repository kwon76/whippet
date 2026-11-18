import sys, time, json
import numpy as np
import tail as T, ring as R
from run_small import run

quiet = lambda *a: None
results = []
CONF = [
    (1500, 12.5, 'A', {}),
    (1200, 12.0, 'C', {}),
    (3000, 13.0, 'B', {}),
    (2000, 12.5, 'A', dict(O=False, I=False, inner=False)),      # no rounding anywhere
    (2500, 13.0, 'A', dict(E=False)),
    (6000, 13.5, 'C', {}),
    (900, 12.0, 'A', dict(clear=False)),                          # partial evaluations always committed
]
rng = np.random.default_rng(7)
import fastring
allok = fastring.test(60)          # fast exact arithmetic is bit-identical with the reference routines
for (m, lb, acc, opts) in CONF:
    res, stmt, sched, proof = run(m, lb, acc, seed=int(rng.integers(1, 10**6)), quiet=True, opts=opts)
    modes = [lv.get('wmode', 'last') + ('' if lv.get('dec', True) else '/nodec') for lv in sched]
    comparable = res['model_total'] - 64 * len(sched)
    line = dict(m=m, logbeta=lb, acc=acc, opts=opts, modes=modes, ok=res['ok'], size=res['size'],
                model=round(comparable), dev=round(100 * (res['size'] - comparable) / comparable, 2))
    # every single bit flip must be rejected
    nflip, rejected = 60, 0
    for _ in range(nflip):
        pos = int(rng.integers(0, len(proof))); bit = int(rng.integers(0, 8))
        bad = bytearray(proof); bad[pos] ^= (1 << bit)
        if not T.verify(stmt, sched, bytes(bad), log=quiet):
            rejected += 1
    line['bitflips_rejected'] = "%d/%d" % (rejected, nflip)
    # truncated / extended proofs
    line['trunc_rejected'] = (not T.verify(stmt, sched, proof[:-1], log=quiet)) and (not T.verify(stmt, sched, proof + b"\x00", log=quiet))
    # same proof against a different statement
    st2 = T.TopStatement(stmt.phis, [stmt.bs[0], (stmt.bs[1] + np.uint64(1)) % R.QU], stmt.fphis, stmt.fbs)
    line['wrong_stmt_rejected'] = not T.verify(st2, sched, proof, log=quiet)
    st3 = T.TopStatement(stmt.phis, stmt.bs, stmt.fphis, [(stmt.fbs[0] + np.uint64(5)) % R.QU])
    line['wrong_fstmt_rejected'] = not T.verify(st3, sched, proof, log=quiet)
    allok &= res['ok'] and rejected == nflip and line['trunc_rejected'] and line['wrong_stmt_rejected'] and line['wrong_fstmt_rejected']
    results.append(line)
    print(line, flush=True)
print("ALL OK" if allok else "FAILURES")
json.dump(results, open("tests.json", "w"), indent=1)
