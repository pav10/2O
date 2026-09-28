\\ field6.gp -- maximal-order tests for a candidate field K (requires lift6.gp).
\\
\\ For K admitting a universal Z-form (Kitaoka / tensor-rank argument, see notes):
\\  (BC)  every alpha in O_K^+ \ Z has tau(alpha) >= nu(K)
\\  (RG)  every alpha in O_K^+ \ Z with tau(alpha) < 3/2 nu(K) is y^2 + c y + b,
\\        y in O_K \ Z, b,c in Z, Var(y) <= tau(alpha)
\\ and for every proper subfield F of K (e = [K:F]):
\\  (BCR) every alpha in O_K^+ \ F has tau(alpha) >= nu_{K/F}
\\  (RGR) every alpha in O_K^+ \ F with tau(alpha) < 3/2 nu_{K/F} is
\\        y^2 + g y + b with y in O_K \ F, g,b in O_F, Var_{K/F}(y) <= tau(alpha).

\\ ---- relative structure w.r.t. a subfield given by [g, h] from nfsubfields ----
\\ returns [Grel (psd, integral = d*e*Var_{K/F} on Z^d coords of O_K), U, Q, OFcoords]
\\ where U is unimodular with first f columns spanning O_F, Q = positive definite
\\ Gram on O_K/O_F in the basis given by the remaining columns of U.
relstruct(Ob, sub) =
{
  my(P = Ob_f(Ob), d = Ob_d(Ob), g = sub[1], h = sub[2], f = poldegree(g), e = d / f,
     th = Ob_th(Ob), E = Ob_E(Ob), hv, cl = List(), blocks, G, Gr, err, Kr, U, Qm, OF);
  hv = vector(d, j, subst(lift(h), x, th[j]));
  my(vals = vecsort(hv, , 8)); for (i = 1, #vals, if (#cl == 0 || abs(vals[i] - cl[#cl]) > 1e-20, listput(cl, vals[i])));
  if (#cl != f, error("relstruct: blocks"));
  blocks = vector(f, b, select(j -> abs(hv[j] - cl[b]) < 1e-20, [1..d], 1));
  G = matrix(d, d, s, t, my(u = E[,s], v = E[,t], acc = 0);
        for (b = 1, f, my(J = blocks[b], mu = sum(i = 1, #J, u[J[i]])/e, mv = sum(i = 1, #J, v[J[i]])/e);
             acc += sum(i = 1, #J, (u[J[i]] - mu)*(v[J[i]] - mv)));
        acc * e);
  Gr = round(G, &err); if (err > -40, error("relstruct: Gram not integral"));
  Kr = matkerint(Gr);                       \\ coords of O_F (saturated, rank f)
  if (#Kr != f, error("relstruct: kernel rank"));
  U = mycompletebasis(Kr);
  Qm = (U~ * Gr * U)[f+1..d, f+1..d];
  [Gr, U, Qm, Kr, f, e];
}

\\ unimodular U whose first columns are the (saturated) columns of Kr
mycompletebasis(Kr) =
{
  my(d = #Kr~, f = #Kr, R = mathnf(Kr~, 1), V = R[2], W = (V~)^-1, U);
  U = matconcat([Kr, W[, 1..d-f]]);
  if (abs(matdet(U)) != 1, error("mycompletebasis: not unimodular"));
  U;
}

\\ coordinates (in the O_K basis) of a polynomial element
kcoords(Ob, r) =
{
  my(d = Ob_d(Ob), B = Ob_B(Ob), M);
  M = matrix(d, d, i, k, polcoeff(B[k], i-1));
  matsolve(M, vector(d, i, polcoeff(lift(Mod(r, Ob_f(Ob))), i-1))~);
}

\\ is the integral vector v in the Z-span of the columns of A ?
inspan(A, v) = my(H = mathnf(A), s); iferr(s = matinverseimage(H, v), E, return(0)); if (#s == 0, return(0)); denominator(s) == 1;

fieldtest(P) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), d = poldegree(P), nu = onu(Ob), L, S, out = List());
  \\ (BC)
  L = otpos(Ob, nu);
  if (#L, return(["BC", nu, L[1][1], oelt(Ob, L[1][2], L[1][3])]));
  \\ (RG)
  L = otpos(Ob, 3/2*nu);
  for (i = 1, #L,
    if (oshape(Ob, L[i][2], L[i][3]) === 0,
      return(["RG", nu, L[i][1], oelt(Ob, L[i][2], L[i][3])])));
  \\ relative tests for every proper subfield
  S = select(s -> my(k = poldegree(s[1])); k > 1 && k < d, nfsubfields(nf));
  for (si = 1, #S,
    my(R = relstruct(Ob, S[si]), Gr = R[1], U = R[2], Qm = R[3], OF = R[4], f = R[5], e = R[6], nur, LL);
    nur = qfminim(Qm)[2] / (d*e);
    \\ (BCR)
    LL = otpos(Ob, nur);
    for (i = 1, #LL,
      my(c = concat([0], LL[i][2]~)~);          \\ full coords of the class vector
      if (c~ * Gr * c != 0, return(["BCR", S[si][1], nur, LL[i][1], oelt(Ob, LL[i][2], LL[i][3])])));
    \\ (RGR)
    LL = otpos(Ob, 3/2*nur);
    for (i = 1, #LL,
      my(c = concat([0], LL[i][2]~)~, al, ta = LL[i][1], ok = 0, qv);
      if (c~ * Gr * c == 0, next);                \\ alpha in F
      al = oelt(Ob, LL[i][2], LL[i][3]);
      qv = qfminim(Qm, floor(ta * d * e))[3];
      for (k = 1, #qv,
        for (sg = 0, 1,
          my(yc = U[, f+1..d] * (if (sg, -1, 1) * qv[,k]), yv, r, A);
          yv = sum(i = 1, d, yc[i] * Ob_B(Ob)[i]);
          r = kcoords(Ob, al - yv^2);
          A = concat(OF, matrix(d, f, i, j, kcoords(Ob, sum(t = 1, d, OF[t,j]*Ob_B(Ob)[t]) * yv)[i]));
          if (inspan(A, r), ok = 1; break));
        if (ok, break));
      if (!ok, return(["RGR", S[si][1], nur, ta, al])));
    listput(out, [S[si][1], nur]));
  ["OK", nu, Vec(out)];
}
