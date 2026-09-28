\\ lift6.gp -- variance geometry of orders in totally real fields, for the
\\ lifting problem (universal Z-forms) in degree 6.
\\
\\ An "order structure" Ob (made by ordinit) stores, for a Z-basis b_1=1,b_2,...,b_d
\\ of an order in K = Q[x]/(f):
\\   Ob.f      the defining polynomial (totally real, irreducible)
\\   Ob.B      the basis as polynomials in x
\\   Ob.th     the real roots of f (all embeddings), high precision
\\   Ob.E      d x d real matrix, E[j,i] = sigma_j(b_i)
\\   Ob.T      exact trace Gram  Tr(b_i b_k)
\\   Ob.tr     exact traces      Tr(b_i)
\\   Ob.G      exact integral Gram of d^2 * Var on Lambda = Ob/Z, in the basis b_2..b_d
\\ All decisions are made with exact arithmetic; floating point is only used to
\\ locate candidates, and every candidate is re-checked exactly.

default(realprecision, 60);

\\ ---------- construction ----------------------------------------------------
\\ record as a vector: [f, d, B, th, E, T, tr, G, den]
ordinit(f, B = 0) =
{
  my(d = poldegree(f), th, E, T, tr, G, den);
  if (B === 0, B = vector(d, i, x^(i-1)));
  if (B[1] != 1, error("ordinit: first basis element must be 1"));
  th = polroots(f); th = vecsort(real(th));
  E = matrix(d, d, j, i, subst(B[i], x, th[j]));
  T = matrix(d, d, i, k, trace(Mod(B[i]*B[k], f)));
  tr = vector(d, i, trace(Mod(B[i], f)));
  G = matrix(d-1, d-1, i, k, d*T[i+1,k+1] - tr[i+1]*tr[k+1]);
  den = denominator(G); G = G*den;   \\ G = den * d^2 * Var  (den=1 for nice bases)
  [f, d, B, th, E, T, tr, G, den];
}
ordZx(f) = ordinit(f);
ordmax(f) = my(nf = nfinit(f)); ordinit(f, nf.zk);

Ob_f(Ob) = Ob[1];
Ob_d(Ob) = Ob[2];
Ob_B(Ob) = Ob[3];
Ob_th(Ob) = Ob[4];
Ob_E(Ob) = Ob[5];
Ob_T(Ob) = Ob[6];
Ob_tr(Ob) = Ob[7];
Ob_G(Ob) = Ob[8];
Ob_den(Ob) = Ob[9];

\\ exact variance of the element with class coordinates c (in basis b_2..b_d)
ovar(Ob, c) = (c~ * Ob_G(Ob) * c) / (Ob_den(Ob) * Ob_d(Ob)^2);
\\ element (as polynomial) from class coords + integer shift n
oelt(Ob, c, n = 0) = n + sum(i = 1, #c, c[i]*Ob_B(Ob)[i+1]);
\\ embeddings of class coords (without shift)
oemb(Ob, c) = Ob_E(Ob)[, 2..Ob_d(Ob)] * c;
\\ exact tau = Tr/d of the element with class coords c and shift n
otau(Ob, c, n = 0) = n + sum(i = 1, #c, c[i]*Ob_tr(Ob)[i+1]) / Ob_d(Ob);

\\ nu(Ob) = least variance of a non-rational element of Ob
onu(Ob) = my(q = qfminim(Ob_G(Ob))); q[2] / (Ob_den(Ob) * Ob_d(Ob)^2);
\\ all shortest class vectors (up to sign)
oshortest(Ob) = my(q = qfminim(Ob_G(Ob))); q[3];

\\ exact total positivity test of a polynomial element a (in x) of K=Q[x]/f
istotpos(f, a) =
{
  my(cp = charpoly(Mod(a, f)));
  \\ all roots of the char poly positive  <=>  no sign changes issue: use polsturm on (0,+oo)
  polsturm(cp, [-oo, 0]) == 0 && subst(cp, x, 0) != 0;
}

\\ ---------- totally positive elements of small trace ------------------------
\\ Returns all [tau, c, n] with alpha = n + sum c_i b_{i+1} totally positive,
\\ non-rational, tau(alpha) < Tmax (strict).  Uses Var(alpha) <= (d-1) tau^2.
otpos(Ob, Tmax) =
{
  my(d = Ob_d(Ob), G = Ob_G(Ob), res = List(), q, bound, V, emb, n, t0, a);
  bound = floor((d-1) * Tmax^2 * d^2 * Ob_den(Ob)) + 1;
  q = qfminim(G, bound)[3];
  for (k = 1, #q,
    for (s = 0, 1,
      my(c = if (s, -q[,k], q[,k]));
      emb = oemb(Ob, c);
      t0 = otau(Ob, c);
      n = floor(-vecmin(emb)) + 1;          \\ least n with all emb + n > 0
      if (vecmin(emb) + n < 1e-40 || vecmin(emb) + n - 1 > -1e-40,
        \\ numerically ambiguous: decide exactly
        while (istotpos(Ob_f(Ob), oelt(Ob, c, n - 1)), n--);
        while (!istotpos(Ob_f(Ob), oelt(Ob, c, n)), n++));
      while (t0 + n < Tmax,
        listput(res, [t0 + n, c, n]);
        n++)));
  res = Vec(res);
  vecsort(res, 1);
}

\\ least trace (tau) of a non-rational totally positive element: search below Tmax
ot(Ob, Tmax) = my(L = otpos(Ob, Tmax)); if (#L, L[1][1], oo);

\\ ---------- downspread quantities -------------------------------------------
\\ for class coords c: tau of the least totally positive shift of +x and of -x
odown(Ob, c) =
{
  my(e = oemb(Ob, c), t0 = otau(Ob, c), nm, np);
  nm = floor(-vecmin(e)) + 1;  np = floor(vecmax(e)) + 1;
  [t0 + nm, -t0 + np];
}

\\ ---------- rigidity (shape) test ------------------------------------------
\\ Is alpha = y^2 + c y + b with y in Ob \ Z, b,c in Z and Var(y) <= tau(alpha)?
\\ Here the search for y is over the order Ob itself.
oshape(Ob, c, n) =
{
  my(d = Ob_d(Ob), f = Ob_f(Ob), a = oelt(Ob, c, n), ta = otau(Ob, c, n), q, G = Ob_G(Ob), cls);
  q = qfminim(G, floor(ta * d^2 * Ob_den(Ob)))[3];
  for (k = 1, #q,
    my(yc = q[,k], y = oelt(Ob, yc), r, rc, lam);
    if (ovar(Ob, yc) > ta, next);
    r = lift(Mod(a - y^2, f));           \\ must be c*y + b
    rc = oclass(Ob, r);
    if (rc === 0, next);
    \\ rc must be an integer multiple of yc
    lam = 0; my(ok = 1);
    for (i = 1, #yc,
      if (yc[i] != 0, lam = rc[i] / yc[i]; break));
    if (denominator(lam) != 1, next);
    if (rc == lam * yc, return([y, lam])));
  0;
}

\\ class coordinates (in b_2..b_d) of a polynomial element r, or 0 if r not in Ob
oclass(Ob, r) =
{
  my(d = Ob_d(Ob), B = Ob_B(Ob), M, v);
  M = matrix(d, d, i, k, polcoeff(lift(Mod(B[k], Ob_f(Ob))), i-1));
  v = matsolve(M, vector(d, i, polcoeff(r, i-1))~);
  if (denominator(v) != 1, return(0));
  v[2..d];
}

\\ ---------- rigidity via squares (no search in Ob needed) --------------------
\\ alpha (polynomial in x) totally positive, generating K.  Rigidity requires
\\   4 alpha - kappa = z^2,  z in Ob_K,  kappa = 4b - c^2 in Z (kappa = 0,3 mod 4),
\\   4 nu <= Var(z) <= 4 tau(alpha).
\\ Returns the list of admissible [kappa, minpoly(z), Var(z)] (empty => excluded).
sqshape(f, a, nu) =
{
  my(m = charpoly(Mod(a, f)), d = poldegree(f), ev, mn, ta, va, A, umax, res = List());
  if (!issquarefree(m), return(-1));      \\ alpha in a proper subfield: not handled
  ev = real(polroots(m)); mn = vecmin(ev);
  ta = -polcoeff(m, d-1) / d;
  va = (polcoeff(m, d-1)^2 - 2*polcoeff(m, d-2)) / d - ta^2;
  A = 4 * ta;
  umax = max(1.8 * A, va / nu) + 1;
  forstep (kap = floor(4*mn - umax) - 1, ceil(4*mn), 1,
    if (kap >= 4*mn, break);
    if (kap % 4 == 1 || kap % 4 == 2, next);
    my(h = subst(m, x, (x^2 + kap)/4) * 4^d, F);
    if (polisirreducible(h), next);
    F = factor(h)[,1];
    for (i = 1, #F,
      my(g = F[i]);
      if (poldegree(g) != d, next);
      my(s1 = -polcoeff(g, d-1), s2 = polcoeff(g, d-1)^2 - 2*polcoeff(g, d-2), vz);
      vz = s2/d - (s1/d)^2;
      if (vz >= 4*nu && vz <= A, listput(res, [kap, g, vz]))));
  Vec(res);
}
