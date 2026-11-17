"""Existence of challenge pairs with difference 3(1-X) and 4(1-X) (used in the remark on weak openings)."""
import numpy as np, json
import ring as R
rng = np.random.default_rng(7)
found = {}
tries = 0
while len(found) < 2 and tries < 400000:
    c = R.challenge(rng); tries += 1
    for k, (a0, a1) in ((3, (2, -1)), (4, (2, -2))):
        if k in found or c[0] != a0 or c[1] != a1:
            continue
        c2 = c.copy(); c2[0], c2[1] = a1, a0
        if R.opnorm(c2) <= 15 and R.opnorm(c) <= 15:
            d = c - c2
            assert d[0] == k and d[1] == -k and not d[2:].any()
            found[k] = dict(c=[int(v) for v in c], op=float(R.opnorm(c)), op2=float(R.opnorm(c2)))
print(tries, {k: (v['op'], v['op2']) for k, v in found.items()})
json.dump(found, open('pairs.json', 'w'))
