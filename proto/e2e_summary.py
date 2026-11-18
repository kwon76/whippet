#!/usr/bin/env python3
"""Summary of proto/e2e.json (end-to-end runs of pcs.py); used by build.py and printed when run directly."""
import json
import os
from math import sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
KINDS = ("uniform", "sparse", "const", "period7", "blockrep", "small8", "maxdig", "maxdig_eq")


def clopper_pearson_upper(k, n, conf=0.95):
    """one-sided upper confidence bound for a binomial proportion with k successes in n trials"""
    if n == 0:
        return 1.0
    if k >= n:
        return 1.0
    lo, hi = k / n, 1.0
    from math import comb
    target = 1 - conf
    for _ in range(100):
        mid = (lo + hi) / 2
        cdf = sum(comb(n, i) * mid ** i * (1 - mid) ** (n - i) for i in range(k + 1))
        if cdf > target:
            lo = mid
        else:
            hi = mid
    return hi


def attempts(r):
    """projections and challenge vectors drawn in one proof: first fold, then every level of the tail"""
    ff = r.get('ff', {}).get('c_tries', 1) or 1
    return ff + sum(s.get('jl_tries', 0) for s in r.get('stats', [])) + sum(s.get('c_tries', 0) or 0 for s in r.get('stats', []))


def load(path=None):
    return json.load(open(path or os.path.join(HERE, "e2e.json")))


def groups(db):
    """{(logN, var): {kind: [runs]}}"""
    out = {}
    for key, runs in db['runs'].items():
        lN, var, kind = key.split("_", 2)
        out.setdefault((int(lN), var), {})[kind] = sorted(runs, key=lambda r: r['rep'])
    return out


def summary(db):
    g = groups(db)
    res = {}
    for (lN, var), kinds in sorted(g.items()):
        par = db['params']["%d_%s" % (lN, var)]
        allr = [r for k in kinds.values() for r in k]
        ok = [r for r in allr if r['ok']]
        fail = [r for r in allr if not r['ok']]
        uni = [r for r in kinds.get('uniform', []) if r['ok']]
        devs = [100.0 * (r['proof_bytes'] / par['model_total'] - 1) for r in ok]
        rec = dict(logN=lN, var=var, runs=len(allr), ok=len(ok), failed=len(fail),
                   fail_levels=sorted({(r.get('fail_level'), r.get('fail_what')) for r in fail}),
                   model_total=par['model_total'], model_ff=par['model_first_fold'], model_com=par['model_commitment'],
                   bytes_min=min(r['proof_bytes'] for r in ok) if ok else None,
                   bytes_max=max(r['proof_bytes'] for r in ok) if ok else None,
                   bytes_uniform=(sum(r['proof_bytes'] for r in uni) / len(uni)) if uni else None,
                   dev_min=min(devs) if devs else None, dev_max=max(devs) if devs else None,
                   ff_bytes=sorted({r['sizes']['first_fold'] for r in ok}),
                   com_bytes=sorted({r['commitment_bytes'] for r in allr}),
                   levels=par['levels'], kinds={})
        for kind, runs in kinds.items():
            okk = [r for r in runs if r['ok']]
            att = [attempts(r) for r in okk]
            fills0 = [r['ff']['fill'] for r in runs if 'fill' in r.get('ff', {})]
            ze = [r['ff']['parts']['ze'] / r['ff']['model_parts']['ze'] for r in runs if 'parts' in r.get('ff', {})]
            lvl = [[s['fill'] for s in r.get('stats', [])] for r in runs]
            rec['kinds'][kind] = dict(
                runs=len(runs), ok=len(okk), fail_upper=clopper_pearson_upper(len(runs) - len(okk), len(runs)),
                att_mean=(sum(att) / len(att)) if att else None, att_max=max(att) if att else None,
                ff_tries_max=max((r.get('ff', {}).get('c_tries') or 0) for r in runs),
                jl_max=max((max((s['jl_tries'] for s in r.get('stats', [])), default=0)) for r in runs),
                c_max=max((max(((s.get('c_tries') or 0) for s in r.get('stats', [])), default=0)) for r in runs),
                fill1_min=min(fills0) if fills0 else None, fill1_max=max(fills0) if fills0 else None,
                fill_max=max((max(l) for l in lvl if l), default=None),
                ze_min=min(ze) if ze else None, ze_max=max(ze) if ze else None,
                bytes_mean=(sum(r['proof_bytes'] for r in okk) / len(okk)) if okk else None,
                t_commit=_median([r['t_commit'] for r in runs]), t_open=_median([r['t_open'] for r in okk]),
                t_verify=_median([r['verify']['t_verify'] for r in okk]),
                rss_prover=max(r.get('rss_open_mb', 0) for r in runs),
                rss_verifier=max((r['verify'].get('rss_mb', 0) for r in okk), default=None))
        # per level of the tail, over all runs of this size and setting
        nl = par['levels']
        lv = []
        for i in range(nl):
            fl = [r['stats'][i]['fill'] for r in ok if len(r['stats']) > i]
            jl = [r['stats'][i]['jl_tries'] for r in ok if len(r['stats']) > i]
            cc = [r['stats'][i]['c_tries'] for r in ok if len(r['stats']) > i]
            lv.append(dict(fill_min=min(fl), fill_max=max(fl), jl_mean=sum(jl) / len(jl), jl_max=max(jl),
                           c_mean=sum(cc) / len(cc), c_max=max(cc)) if fl else None)
        ffc = [r['ff']['c_tries'] for r in ok]
        rec['first_fold_tries_mean'] = sum(ffc) / len(ffc) if ffc else None
        rec['first_fold_tries_max'] = max(ffc) if ffc else None
        rec['level_stats'] = lv
        rec['t_commit'] = _median([r['t_commit'] for r in uni]) if uni else None
        rec['t_open'] = _median([r['t_open'] for r in uni]) if uni else None
        rec['t_verify'] = _median([r['verify']['t_verify'] for r in uni]) if uni else None
        rec['t_setup_verifier'] = _median([r['verify']['t_setup'] for r in uni]) if uni else None
        rec['rss_prover'] = max(r.get('rss_open_mb', 0) for r in uni) if uni else None
        rec['rss_verifier'] = max(r['verify'].get('rss_mb', 0) for r in uni) if uni else None
        res[(lN, var)] = rec
    return res


def _median(a):
    a = sorted(x for x in a if x is not None)
    if not a:
        return None
    k = len(a)
    return a[k // 2] if k % 2 else 0.5 * (a[k // 2 - 1] + a[k // 2])


def neg_summary(db):
    out = {}
    for key, parts in db.get('neg', {}).items():
        p0 = [p for p in parts if p['part'] == 0]
        rec = dict(honest=all(p['honest'] for p in parts), flips=sum(p['flips'] for p in parts),
                   flips_rejected=sum(p['flips_rejected'] for p in parts),
                   flips_ff=sum(p['flips_in_first_fold'] for p in parts),
                   flips_ff_rejected=sum(p['flips_in_first_fold_rejected'] for p in parts))
        if p0:
            for k in ('wrong_y', 'noncanonical_y', 'wrong_x', 'other_commitment', 'flipped_commitment',
                      'proof_for_other_point', 'other_point_with_its_value', 'truncated', 'appended', 'empty'):
                rec[k] = p0[0].get(k)
            rec['edge_points_ok'] = p0[0].get('edge_points_ok')
        out[key] = rec
    return out


if __name__ == "__main__":
    import sys
    db = load(sys.argv[1] if len(sys.argv) > 1 else None)
    for (lN, var), r in summary(db).items():
        print("N=2^%d %-4s runs %3d ok %3d failed %d %s  bytes %s..%s (model %.0f, dev %+.3f..%+.3f%%)  ff %s (model %.1f)  com %s (model %.1f)" % (
            lN, var, r['runs'], r['ok'], r['failed'], r['fail_levels'] or '', r['bytes_min'], r['bytes_max'],
            r['model_total'], r['dev_min'] or 0, r['dev_max'] or 0, r['ff_bytes'], r['model_ff'], r['com_bytes'], r['model_com']))
        print("   t commit/open/verify %s/%s/%s s, rss prover/verifier %s/%s MB; first fold tries mean %.2f max %s" % (
            *("%.1f" % x if x is not None else "-" for x in (r['t_commit'], r['t_open'], r['t_verify'], r['rss_prover'], r['rss_verifier'])),
            r['first_fold_tries_mean'] or 0, r['first_fold_tries_max']))
        for k in KINDS:
            if k in r['kinds']:
                s = r['kinds'][k]
                print("   %-10s %2d/%2d ok (fail<=%.2f) att mean %s max %s  ff_tries<=%d jl<=%d c<=%d  fill1 %.3f..%.3f  fillmax %s  ze %.2f..%.2f" % (
                    k, s['ok'], s['runs'], s['fail_upper'], "%.1f" % s['att_mean'] if s['att_mean'] else '-', s['att_max'],
                    s['ff_tries_max'], s['jl_max'], s['c_max'], s['fill1_min'] or 0, s['fill1_max'] or 0,
                    "%.3f" % s['fill_max'] if s['fill_max'] else '-', s['ze_min'] or 0, s['ze_max'] or 0))
        for i, l in enumerate(r['level_stats']):
            if l:
                print("   level %d: fill %.3f..%.3f  jl mean %.2f max %d  c mean %.2f max %d" % (
                    i + 1, l['fill_min'], l['fill_max'], l['jl_mean'], l['jl_max'], l['c_mean'], l['c_max']))
    for k, v in neg_summary(db).items():
        print("neg", k, v)
