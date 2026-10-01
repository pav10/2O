"""
deltaF.py -- rigorous upper bound for the positive covering constant

    delta_F = sup_{u in R^f}  min { (1/f) sum_k (p_k - u_k) : p in L, p > u componentwise }

of the Minkowski lattice L = sigma(O_F) of a totally real field F of degree f.
(For F = Q, delta = 1.)  It enters the relative downspread lemma:
for y in O_K \\ F there is gamma in O_F with y + gamma totally positive and
tau(y + gamma) <= (mean over blocks of the within-block downspread) + delta_F.

Method: h(u) = min_{p>u} sum(p-u) is L-periodic and satisfies
h(u) <= h(u') + sum(u'-u) whenever u <= u' componentwise.  Hence for a grid of
step eps covering a fundamental parallelepiped (every u has a grid point u' >= u
with u'-u <= eps coordinatewise),  delta_F <= max_grid h / f + eps.

Rigour in floating point: h(u') is only ever bounded from ABOVE, by sum(p - u') for a lattice
point p that is accepted only if p > u' + MARGIN in double precision (MARGIN = 1e-9, far above
the rounding error of the embeddings, ~1e-14).  Every accepted p therefore satisfies p > u'
exactly, and the reported bound adds MARGIN for the rounding of the sums.  The grid is built
from integer multiples of eps, so it covers [lo, hi] with an upper neighbour within eps.
"""
import numpy as np, itertools, sys, os

MARGIN = 1e-9

def delta_bound(B, eps, H=6.0, chunk=200000):
    B = np.array(B, dtype=float)          # f x f, columns = embeddings of a Z-basis
    f = B.shape[0]
    # bounding box of the fundamental parallelepiped
    corners = np.array([B @ np.array(c) for c in itertools.product([0, 1], repeat=f)])
    lo, hi = corners.min(0), corners.max(0)
    # lattice points in a box [lo, hi + H]
    Binv = np.linalg.inv(B)
    box = np.array(list(itertools.product(*[[lo[k], hi[k] + H] for k in range(f)])))
    coef = box @ Binv.T
    cmin, cmax = np.floor(coef.min(0)) - 1, np.ceil(coef.max(0)) + 1
    grids = np.meshgrid(*[np.arange(cmin[i], cmax[i] + 1) for i in range(f)], indexing='ij')
    C = np.stack([g.ravel() for g in grids], 1)
    P = C @ B.T
    keep = np.all(P > lo - 1e-9, 1) & np.all(P < hi + H + 1e-9, 1)
    P = P[keep]
    Ps = P.sum(1)
    order = np.argsort(Ps); P, Ps = P[order], Ps[order]
    axes = [lo[k] + eps * np.arange(0, int(np.ceil((hi[k] - lo[k]) / eps)) + 2) for k in range(f)]
    G = np.stack([g.ravel() for g in np.meshgrid(*axes, indexing='ij')], 1)
    # restrict to grid points u' that can be the upper neighbour of some u in the parallelepiped
    best = -np.inf; arg = None
    for s in range(0, len(G), chunk):
        U = G[s:s + chunk]
        mask = np.all(P[None, :, :] > U[:, None, :] + MARGIN, axis=2)
        first = np.where(mask.any(1), mask.argmax(1), -1)
        if (first < 0).any():
            raise RuntimeError("H too small")
        h = Ps[first] - U.sum(1)
        i = h.argmax()
        if h[i] > best: best, arg = h[i], U[i]
    return best / f + eps + MARGIN, best / f, arg

FIELDS = {
  "Q(sqrt5)":  [[1, -0.61803398874989484820], [1, 1.61803398874989484820]],
  "Q(sqrt2)":  [[1, -1.41421356237309504880], [1, 1.41421356237309504880]],
  "Q(z7)+":    [[1, -1.24697960371746706105, 1.80193773580483825247],
                [1, 0.44504186791262880858, -1.24697960371746706105],
                [1, 1.80193773580483825247, 0.44504186791262880858]],
  "Q(z9)+":    [[1, -1.53208888623795607040, 1.87938524157181676811],
                [1, -0.34729635533386069770, -1.53208888623795607040],
                [1, 1.87938524157181676811, -0.34729635533386069770]],
  # degree 8, Case B: the intermediate quartic fields F' containing Q(sqrt5)
  "F725":      [[1, -1.61803398874989484820, -1.35567429397808222658, 2.19352708533105393856],
                [1, 0.61803398874989484820, -0.47725999647401964454, -0.29496289929159911192],
                [1, -1.61803398874989484820, 0.73764030522818737837, -1.19352708533105393856],
                [1, 0.61803398874989484820, 2.09529398522391449275, 1.29496289929159911192]],
  "F1600":     [[1, 1.61803398874989484820, -1.41421356237309504880, -2.28824561127073719040],
                [1, -0.61803398874989484820, 1.41421356237309504880, -0.87403204889764214160],
                [1, -0.61803398874989484820, -1.41421356237309504880, 0.87403204889764214160],
                [1, 1.61803398874989484820, 1.41421356237309504880, 2.28824561127073719040]],
}
EPS = {2: 0.0005, 3: 0.01, 4: 0.1}
CHUNK = {2: 200000, 3: 20000, 4: 500}

if __name__ == "__main__":
    for name, E in FIELDS.items():
        f = len(E)
        if len(sys.argv) > 1 and name not in sys.argv[1:]: continue
        eps = float(os.environ.get("DELTA_EPS", EPS[f]))
        ub, grid_max, arg = delta_bound(E, eps, H=4.0 if f == 4 else 6.0, chunk=CHUNK[f])
        print("%-10s f=%d  max_grid h/f = %.10f   =>  delta_F <= %.10f   (worst u' = %s)" % (name, f, grid_max, ub, arg))
