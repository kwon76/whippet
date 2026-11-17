"""The first fold of Greyhound [NS24, Fig. 1 and Table 4] as the source of the tail instance.

shape parameters (m, r, b0, delta0, b, delta) are those of [NS24, Table 4].  Ranks:
  accounting C: the published ranks (n, n1) = (18, 7) and the challenge space of [NS24] (tau = 64);
  accountings A, B: the challenge space of [BS23] (tau = 71, T = 15) and ranks recomputed from the norms
     required by the security proof: the outer commitments u, v for 2*SL*beta', the inner commitments
     for 8*T*(b+1)*SL*beta' (A), resp. 2*beta' and 8*T*gamma (B), where beta' is the norm bound of the
     relation handed to the tail.
"""
from math import log2, sqrt, ceil, floor
import params as P

SHAPE = {26: (3156, 333, 6, 5, 7, 5), 28: (6312, 665, 5, 6, 6, 5), 30: (12625, 1329, 4, 8, 6, 5)}


def first_fold(logN, rounded, n=18, n1=7, TG=64, eta=0.10):
    """Witness handed to the tail.  Rounded variant: the inner commitments lose Dt low bits and the two
    outer commitments Du low bits (rules as in the docstring of run_eval.greyhound_instance)."""
    m, r, lb0, d0, lb, dl = SHAPE[logN]
    D, LQ = P.D, P.LOGQ
    b0, b = 2.0**lb0, 2.0**lb
    vz = b0*b0/12*r*TG

    def widths(bits, table_rule):
        if table_rule:
            k = max(1, floor(bits/lb))
            if bits - (k - 1)*lb > lb + 2:
                k += 1
        else:
            k = max(1, ceil(bits/lb))
        return [lb]*(k - 1) + [bits - (k - 1)*lb]
    var = lambda ws: sum(4.0**w/12 for w in ws)
    wc = lambda ws: sum(4.0**(w - 1) for w in ws)
    w_full = widths(LQ, True)
    zvar = b*b/12 + vz/(b*b)
    norm2 = zvar*d0*m*D + var(w_full)*r*D
    ref = norm2 + var(w_full)*n*r*D
    Dt = Du = 0
    if rounded:
        cap = lambda x: max(0, min(20, floor(0.5*log2(max(1.0, x)))))
        Dt = min(cap(0.1*b0*b0*d0*m*P.LOWDIV/(12*n)), cap(eta*ref*P.LOWDIV*b*b/(n*D*TG*r)))
        w_t = widths(LQ - Dt, False)
        freed = (wc(w_full) - wc(w_t))*n*r*D
        Du = min(1 + cap(freed/(n1*D)) if freed > n1*D else 0, cap(eta*ref*P.LOWDIV/(2*n1*D)))
    else:
        w_t = w_full
    kw, kt = len(w_full), len(w_t)
    W = kw*r + kt*n*r + 2*d0*m
    norm2 += var(w_t)*n*r*D
    gz2 = vz*d0*m*D                                       # expected squared norm of z
    ge2 = 0.0
    if rounded:
        W += 2*n1 + 2*n
        ge2 = n*D*4.0**Dt*TG*r/P.LOWDIV
        norm2 += 2*n1*D*4.0**Du/P.LOWDIV + n*D*(b*b/12 + 4.0**Dt*TG*r/(P.LOWDIV*b*b))
    return dict(W=W, beta=sqrt(norm2), Dt=Dt, Du=Du, n=n, n1=n1, m=m, r=r, lb=lb, lb0=lb0, d0=d0,
                kt=kt, TG=TG, gbar=sqrt(gz2 + ge2), com=n1*D*P.hibits(Du)/8)


def instance(logN, rounded, acc):
    if acc == 'C':
        return first_fold(logN, rounded)
    n, n1 = 18, 7
    for _ in range(50):
        g = first_fold(logN, rounded, n=n, n1=n1, TG=P.TAU)
        if acc == 'A':
            nn1 = P.rank(2*P.SL*g['beta'])
            nn = P.rank(8*P.T*(2.0**g['lb'] + 1)*P.SL*g['beta'])
        else:
            nn1 = P.rank(2*g['beta'])
            nn = P.rank(8*P.T*g['gbar'])
        if (nn, nn1) == (n, n1):
            return g
        n, n1 = nn, nn1
    raise RuntimeError("ranks do not converge")


def first_fold_gc(logN, n, n1, theta, tmar, eta=0.10, lb0=None, lb=None):
    """First fold with norm bounds that hold for every polynomial: all digits and low parts are counted
    with their maximal absolute values, and the prover repeats the challenge until
    ||z||^2 + ||e||^2 <= theta*tau*(||s||_max^2 + Lambda^2)."""
    m, r, lb0_t, d0_t, lb_t, dl = SHAPE[logN]
    D, LQ, TAU = P.D, P.LOGQ, P.TAU
    # digit base of the polynomial and base of the commitment digits: free in this setting (the bases of
    # [NS24, Table 4] leave a top digit of up to eight bits, whose maximal value dominates the norm)
    lb0 = lb0_t if lb0 is None else lb0
    lb = lb_t if lb is None else lb
    d0 = ceil(LQ/lb0)
    b0, b = 2.0**lb0, 2.0**lb

    def widths(bits, table_rule):
        k = max(1, ceil(bits/lb))
        return [lb]*(k - 1) + [bits - (k - 1)*lb]
    wcb = lambda ws: sum(4.0**(w - 1) for w in ws[:-1]) + (2.0**(ws[-1] - 1) + 1)**2
    cap = lambda x: max(0, min(20, floor(0.5*log2(max(1.0, x)))))
    top0 = 2.0**(LQ - 1)/b0**(d0 - 1) + 1
    bs2 = m*r*D*((d0 - 1)*b0*b0/4 + top0*top0)            # ||s||^2 for the digits of any polynomial
    Dt = 1 + cap(0.1*bs2/(r*n*D)) if 0.1*bs2 >= r*n*D else 0
    w_full, w_t = widths(LQ, True), widths(LQ - Dt, False)
    lam2 = r*n*D*4.0**(Dt - 1) if Dt > 0 else 0.0
    qmax = theta*TAU*(bs2 + lam2)
    nz = d0*m + (n if Dt > 0 else 0)
    det0 = nz*D*b*b/4
    norm2 = det0 + ((sqrt(qmax) + sqrt(det0))/b)**2 + wcb(w_full)*r*D + wcb(w_t)*n*r*D
    Du = 1 + cap(eta*norm2/(2*n1*D)) if eta*norm2 >= 2*n1*D else 0
    if Du > 0:
        norm2 += 2*n1*D*4.0**(Du - 1)
    kw, kt = len(w_full), len(w_t)
    W = kw*r + kt*n*r + 2*d0*m + (2*n1 if Du > 0 else 0) + (2*n if Dt > 0 else 0)
    return dict(W=W, beta=sqrt(norm2), Dt=Dt, Du=Du, n=n, n1=n1, m=m, r=r, lb=lb, lb0=lb0, d0=d0,
                kt=kt, TG=TAU, gbar=sqrt(qmax), qmax=qmax, bs2=bs2, com=n1*D*P.hibits(Du)/8)


def instance_gc(logN, theta, tmar, lb0=4, lb=6):
    n, n1 = 18, 7
    tm = sqrt(1.0 + tmar)
    for _ in range(50):
        g = first_fold_gc(logN, n, n1, theta, tmar, lb0=lb0, lb=lb)
        nn1 = P.rank(2*P.SL*tm*g['beta'])
        nn = P.rank(8*P.T*(2.0**g['lb'] + 1)*P.SL*tm*g['beta'])
        if nn is None or nn1 is None:
            return None
        if (nn, nn1) == (n, n1):
            return g
        n, n1 = nn, nn1
    raise RuntimeError("ranks do not converge")


def own_bits(g, lifted):
    """messages of the first fold: v (or its high part), the lifted evaluation, two seeds"""
    y = 3*P.LOGQ if lifted else P.D*P.LOGQ
    return g['n1']*P.D*P.hibits(g['Du']) + y + 2*128
