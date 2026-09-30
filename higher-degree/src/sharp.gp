\\ sharp.gp -- sharpened rigidity (requires lift6.gp, field6.gp).
\\ In Proposition "rigidity", alpha = y^2 + c y + b with c = B(e,r), b = Q(r)/2, hence
\\   kappa := 4b - c^2 = Q(e)Q(r) - B(e,r)^2 >= 0          (Cauchy-Schwarz in L')
\\ and tau(alpha) = Var(y) + (tau(y) + c/2)^2 + kappa/4, so kappa + tau(z)^2 <= 4(tau(alpha) - nu)
\\ for z = 2y + c.  Relative version: alpha = y^2 + g y + b with g, b in O_F and
\\ 4b - g^2 totally >= 0 in O_F (Cauchy-Schwarz at every real place of F).
\\ Square dichotomy: tau(alpha) < min(3/2 nu, nu + 3/4)  =>  alpha is a square in O_K.

\\ sharp version of sqshape (alpha generating K): kappa in [0, min(4 min alpha, 4(tau - nu))]
sqshape2(f, a, nu) =
{
  my(m = charpoly(Mod(a, f)), d = poldegree(f), ev, mn, ta, res = List());
  if (!issquarefree(m), return(-1));
  ev = real(polroots(m)); mn = vecmin(ev);
  ta = -polcoeff(m, d-1) / d;
  for (kap = 0, min(ceil(4*mn) - 1, floor(4*(ta - nu))),
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
      if (vz >= 4*nu && kap + (s1/d)^2 <= 4*(ta - nu), listput(res, [kap, g, vz]))));
  Vec(res);
}

\\ sharp version of oshape (search y in the order): kappa = 4b - c^2 >= 0 required
oshape2(Ob, c, n, nu) =
{
  my(d = Ob_d(Ob), f = Ob_f(Ob), a = oelt(Ob, c, n), ta = otau(Ob, c, n), q);
  q = qfminim(Ob_G(Ob), floor(ta * d^2 * Ob_den(Ob)))[3];
  for (k = 1, #q,
    for (sg = 0, 1,
      my(yc = if (sg, -q[,k], q[,k]), y = oelt(Ob, yc), r, rc, lam, b);
      if (ovar(Ob, yc) > ta, next);
      r = lift(Mod(a - y^2, f));
      rc = oclass(Ob, r);
      if (rc === 0, next);
      lam = 0;
      for (i = 1, #yc, if (yc[i] != 0, lam = rc[i] / yc[i]; break));
      if (denominator(lam) != 1 || rc != lam * yc, next);
      b = lift(Mod(a - y^2 - lam*y, f));
      if (poldegree(b) > 0, next);
      b = polcoeff(b, 0);
      if (denominator(b) != 1, next);
      if (4*b - lam^2 >= 0, return([y, lam, b]))));
  0;
}

\\ sharp maximal-order test: absolute rigidity with kappa >= 0 and relative rigidity with
\\ 4 beta - gamma^2 totally >= 0; also reports window elements that are not squares.
fieldtest2(P) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), d = poldegree(P), nu = onu(Ob), L, S, out = List());
  L = otpos(Ob, nu);
  if (#L, return(["BC", nu, L[1][1], oelt(Ob, L[1][2], L[1][3])]));
  L = otpos(Ob, 3/2*nu);
  for (i = 1, #L,
    if (oshape2(Ob, L[i][2], L[i][3], nu) === 0,
      return(["RG2", nu, L[i][1], oelt(Ob, L[i][2], L[i][3])])));
  S = select(s -> my(k = poldegree(s[1])); k > 1 && k < d, nfsubfields(nf));
  for (si = 1, #S,
    my(R = relstruct(Ob, S[si]), Gr = R[1], U = R[2], Qm = R[3], OF = R[4], f = R[5], e = R[6], nur, LL);
    nur = qfminim(Qm)[2] / (d*e);
    LL = otpos(Ob, nur);
    for (i = 1, #LL,
      my(c = concat([0], LL[i][2]~)~);
      if (c~ * Gr * c != 0, return(["BCR", S[si][1], nur, LL[i][1], oelt(Ob, LL[i][2], LL[i][3])])));
    LL = otpos(Ob, 3/2*nur);
    for (i = 1, #LL,
      my(c = concat([0], LL[i][2]~)~, al, ta = LL[i][1], ok = 0, qv);
      if (c~ * Gr * c == 0, next);
      al = oelt(Ob, LL[i][2], LL[i][3]);
      qv = qfminim(Qm, floor(ta * d * e))[3];
      for (k = 1, #qv,
        for (sg = 0, 1,
          my(yc = U[, f+1..d] * (if (sg, -1, 1) * qv[,k]), yv, r, A, sol);
          yv = sum(i = 1, d, yc[i] * Ob_B(Ob)[i]);
          r = kcoords(Ob, al - yv^2);
          A = concat(OF, matrix(d, f, i, j, kcoords(Ob, sum(t = 1, d, OF[t,j]*Ob_B(Ob)[t]) * yv)[i]));
          sol = iferr(matinverseimage(A, r), E, 0);
          if (sol === 0 || #sol == 0 || denominator(sol) != 1, next);
          \\ beta = sum sol[1..f] * OF-basis, gamma = sum sol[f+1..2f] * OF-basis
          my(be = sum(j = 1, f, sol[j] * sum(t = 1, d, OF[t,j]*Ob_B(Ob)[t])),
             ga = sum(j = 1, f, sol[f+j] * sum(t = 1, d, OF[t,j]*Ob_B(Ob)[t])), kk);
          kk = lift(Mod(4*be - ga^2, P));
          if (kk == 0 || polsturm(charpoly(Mod(kk, P)), [-oo, 0]) == 0,
            ok = 1; break));
        if (ok, break));
      if (!ok, return(["RGR2", S[si][1], nur, ta, al])));
    listput(out, [S[si][1], nur]));
  ["OK", nu, Vec(out)];
}
