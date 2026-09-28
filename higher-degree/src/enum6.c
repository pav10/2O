/*
 * enum6.c -- enumerate the possible shortest vectors x_1 of the variance
 * lattice of a totally real sextic field K admitting a universal Z-form,
 * in the case K = Q(x_1).
 *
 * Output: every monic f = x^6 + a1 x^5 + ... + a6 in Z[x] such that
 *   (E1) f has only real roots,
 *   (E2) a1 in {0,-1,-2,-3}              (normalisation x -> +-x + k),
 *   (E3) S := sum (theta_i - tau)^2 <= SMAX = 15 + 3 sqrt(21)
 *        (i.e. Var = S/6 <= (5+sqrt21)/2, Lemma "nu bound"),
 * and whose roots pass the filters (with V = S/6, tau = -a1/6):
 *   (F1) downspread  tau - (ceil(theta_min) - 1) >= V,
 *   (F2) upspread    (floor(theta_max) + 1) - tau >= V,
 *   (F3) no monic P in Z[t] of degree 2 with P(theta_i) > 0 for all i and
 *        P(tau) < 0.
 * Each filter expresses tau(alpha) >= nu(K) = V for a specific totally
 * positive non-rational alpha in Z[x_1] (x_1 - n, -x_1 - n', P(x_1)).
 * All comparisons are made with safety margins in the direction that keeps
 * polynomials; exact re-verification of everything is done in PARI/GP.
 *
 * Enumeration: Rolle/interlacing ("Robinson") method.  With
 *   g_m(x) = sum_{i=0}^m a_i C(6-i, m-i) x^{m-i}   (so g_6 = f, g_m' = (7-m) g_{m-1}),
 * f real-rooted implies every g_m real-rooted with roots interlacing those of
 * g_{m-1}; given a_0..a_{m-1} this confines a_m (the constant term of g_m)
 * to an explicit interval.  The enumeration is exhaustive.
 *
 * usage: enum6 a1 a2 [smax [all]]   (the subtree with the given a1,a2)
 * output lines: a1 a2 a3 a4 a5 a6
 * stderr: statistics.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef long double ld;
#ifndef DEG
#define DEG 6
#endif
#define N DEG
#define EPS 1e-9L

static long a[N + 1];
static ld Cb[N + 1][N + 1];
static ld Rt[N + 1][N + 1];      /* Rt[m][0..m-1]: roots of g_m, ascending */
static ld LB, UB, SMAX, TAU, VAR;
static long long n_leaf = 0, n_out = 0, n_f1 = 0, n_f2 = 0, n_f3 = 0;
static int NOFILTER = 0;

static ld g(int m, ld x)
{
    ld v = 0;
    for (int i = 0; i <= m; i++) v = v * x + (ld)a[i] * Cb[N - i][m - i];
    return v;
}

static ld findroot(int m, ld lo, ld hi)
{
    ld flo = g(m, lo), fhi = g(m, hi);
    if (flo == 0) return lo;
    if (fhi == 0) return hi;
    if ((flo > 0) == (fhi > 0)) return fabsl(flo) < fabsl(fhi) ? lo : hi; /* tangency */
    ld x = 0.5L * (lo + hi);
    for (int it = 0; it < 300; it++) {
        ld fx = g(m, x);
        if (fx == 0) return x;
        if ((fx > 0) == (flo > 0)) { lo = x; flo = fx; } else { hi = x; fhi = fx; }
        ld d = (ld)(N - m + 1) * g(m - 1, x);
        ld xn = (d != 0) ? x - fx / d : 0.5L * (lo + hi);
        if (!(xn > lo && xn < hi)) xn = 0.5L * (lo + hi);
        if (fabsl(xn - x) <= 1e-18L * (1 + fabsl(x))) return xn;
        x = xn;
        if (hi - lo <= 1e-18L * (1 + fabsl(x))) return x;
    }
    return x;
}

/* F3: is there a monic integer quadratic with both roots in the gap
 * (th[k], th[k+1]) that contains tau, and with tau strictly between them? */
static int quad_violation(const ld *th)
{
    int k = -1;
    for (int i = 0; i + 1 < N; i++)
        if (th[i] < TAU && TAU < th[i + 1]) { k = i; break; }
    if (k < 0) return 0;               /* tau (numerically) at a root: keep */
    ld l = th[k], r = th[k + 1];
    long bmin = (long)ceill(-2 * r), bmax = (long)floorl(-2 * l);
    for (long b = bmin; b <= bmax; b++) {
        ld lo1 = -l * l - b * l, lo2 = -r * r - b * r;
        ld lo = lo1 > lo2 ? lo1 : lo2;
        ld up = -TAU * TAU - b * TAU;
        /* need integer c with lo < c < up; declare violation only with margin */
        long c = (long)floorl(lo + 1e-7L) + 1;
        if ((ld)c > lo + 1e-7L && (ld)c < up - 1e-7L) return 1;
    }
    return 0;
}

static void leaf(void)
{
    ld br[N + 2], th[N];
    n_leaf++;
    br[0] = LB;
    for (int i = 0; i < N - 1; i++) br[i + 1] = Rt[N - 1][i];
    br[N] = UB;
    for (int i = 0; i < N; i++) th[i] = findroot(N, br[i], br[i + 1]);
    if (NOFILTER) { n_out++; for (int i = 1; i <= N; i++) printf("%ld%c", a[i], i == N ? '\n' : ' '); return; }
    /* F1: downspread.  theta_min is never an integer for irreducible f. */
    ld dm = TAU - (ceill(th[0]) - 1);
    if (dm < VAR - 1e-7L) { n_f1++; return; }
    ld dp = (floorl(th[N - 1]) + 1) - TAU;
    if (dp < VAR - 1e-7L) { n_f2++; return; }
    if (N >= 3 && quad_violation(th)) { n_f3++; return; }   /* P(x_1) non-rational needs deg x_1 >= 3 */
    n_out++;
    for (int i = 1; i <= N; i++) printf("%ld%c", a[i], i == N ? '\n' : ' ');
}

static void rec(int m)
{
    ld lo = -1e30L, hi = 1e30L;
    a[m] = 0;
    for (int i = 0; i < m - 1; i++) {
        ld hv = g(m, Rt[m - 1][i]);
        int k = m - (i + 1);
        if (k % 2 == 0) { if (-hv > lo) lo = -hv; }
        else            { if (-hv < hi) hi = -hv; }
    }
    long amin = (long)ceill(lo - EPS * (1 + fabsl(lo)));
    long amax = (long)floorl(hi + EPS * (1 + fabsl(hi)));
    for (long am = amin; am <= amax; am++) {
        a[m] = am;
        if (m == N) { leaf(); continue; }
        ld br[N + 2];
        br[0] = LB;
        for (int i = 0; i < m - 1; i++) br[i + 1] = Rt[m - 1][i];
        br[m] = UB;
        for (int i = 0; i < m; i++) Rt[m][i] = findroot(m, br[i], br[i + 1]);
        rec(m + 1);
    }
    a[m] = 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s a1 a2\n", argv[0]); return 1; }
    for (int n = 0; n <= N; n++)
        for (int k = 0; k <= n; k++)
            Cb[n][k] = (k == 0 || k == n) ? 1 : Cb[n - 1][k - 1] + Cb[n - 1][k];
    /* d*nu_max with d*nu >= 2(nu-1)^2:  nu_max = ((d+4)+sqrt((d+4)^2-16))/4 */
    SMAX = (ld)N * ((N + 4) + sqrtl((ld)(N + 4) * (N + 4) - 16)) / 4 + 1e-9L;
    a[0] = 1;
    a[1] = atol(argv[1]);
    a[2] = atol(argv[2]);
    /* optional: argv[3] = SMAX override (0 = default), argv[4] = "all" (no filters) */
    if (argc > 3 && strtold(argv[3], NULL) > 0) SMAX = strtold(argv[3], NULL);
    if (argc > 4) NOFILTER = 1;
    ld S = (ld)(N - 1) * a[1] * a[1] / N - 2 * (ld)a[2];
    if (S <= 0 || S > SMAX) return 0;
    TAU = -(ld)a[1] / N;
    VAR = S / N;
    ld rad = sqrtl((N - 1) * S / N);
    LB = TAU - rad - 0.01L;
    UB = TAU + rad + 0.01L;
    Rt[1][0] = TAU;
    if (N == 2) {
        leaf();
    } else {
        /* roots of g_2 */
        ld br[3] = {LB, TAU, UB};
        for (int i = 0; i < 2; i++) Rt[2][i] = findroot(2, br[i], br[i + 1]);
        rec(3);
    }
    fprintf(stderr, "a1=%ld a2=%ld S=%.4Lf leaves=%lld F1=%lld F2=%lld F3=%lld out=%lld\n",
            a[1], a[2], S, n_leaf, n_f1, n_f2, n_f3, n_out);
    return 0;
}
