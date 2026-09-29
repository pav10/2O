"""
refined_bound.py -- integrality-refined version of Corollary 4.2 (nu bound).

For the shortest vector x_1 (Tr x_1 = t, tau = t/d, normalised t in {0..floor(d/2)}), the
downspread conditions (Lemma 4.1 with Proposition 3.3) are
    tau - ceil(theta_1) + 1 >= V  and  floor(theta_d) + 1 - tau >= V,   V = Var(x_1) = S/d,
i.e.  theta_1 <= L := floor(tau - V + 1)  and  theta_d >= U := ceil(tau + V - 1).
With a = tau - L, b = U - tau the remaining d-2 deviations sum to a' - b', so
    S = sum (theta_j - tau)^2 >= min_{a'>=a, b'>=b} a'^2 + b'^2 + (a'-b')^2/(d-2).
S takes the values (d-1) t^2/d - 2 a_2 (a_2 in Z).  We list, per t, the largest feasible S.
"""
import math, sys

def minquad(a, b, d):
    # min of a'^2 + b'^2 + (a'-b')^2/(d-2) over a' >= a, b' >= b (convex; check KKT candidates)
    best = float('inf')
    k = 1.0 / (d - 2)
    cands = [(a, b)]
    # a' free (b' = b): derivative 2a' + 2k(a'-b) = 0 -> a' = k b/(1+k)
    ap = k * b / (1 + k); cands.append((max(a, ap), b))
    bp = k * a / (1 + k); cands.append((a, max(b, bp)))
    for (x, y) in cands:
        if x >= a - 1e-12 and y >= b - 1e-12:
            best = min(best, x * x + y * y + k * (x - y) ** 2)
    return best

def feasible(S, t, d):
    tau = t / d; V = S / d
    L = math.floor(tau - V + 1 + 1e-12); U = math.ceil(tau + V - 1 - 1e-12)
    a = tau - L; b = U - tau
    return S >= minquad(a, b, d) - 1e-9

for d in map(int, sys.argv[1:] or [6, 7, 8]):
    smax = d * ((d + 4) + math.sqrt((d + 4) ** 2 - 16)) / 4
    print("d=%d  continuous bound: nu <= %.4f, S <= %.3f" % (d, smax / d, smax))
    worst = 0
    for t in range(0, d // 2 + 1):
        best = None
        for a2 in range(0, -200, -1):
            S = (d - 1) * t * t / d - 2 * a2
            if S <= 0: continue
            if S > smax + 1e-9: break
            if feasible(S, t, d): best = (a2, S)
        worst = max(worst, best[1])
        print("   Tr x_1 = %d: largest feasible S = %.3f (a2 = %d), nu <= %.4f" % (t, best[1], best[0], best[1] / d))
    print("   => refined: nu <= %.4f  (S <= %.3f)" % (worst / d, worst))
