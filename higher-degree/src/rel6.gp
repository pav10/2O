\\ rel6.gp -- the imprimitive case: the shortest vector x_1 of Lambda_K lies in a
\\ proper subfield F (quadratic or cubic).  Then nu(K) = nu(F), and we enumerate
\\ K = F(y) with y a relative shortest vector: Var_{K/F}(y) = nu_{K/F} <= NUREL(F).
\\ Requires lift6.gp.
\\
\\ Tests are made in the order O_F[y] (a suborder of O_K; no maximal order needed):
\\   NSA  some element of O_F[y] \ Z has Var < nu(F)          (nu(K) = nu(F) violated)
\\   BCA  some alpha in O_F[y]^+ \ Z has tau(alpha) < nu(F)      (absolute budget-cost)
\\   NSR  some element of O_F[y] \ F has Var_{K/F} < Var_{K/F}(y) (y not relative shortest)
\\   BCR  some alpha in O_F[y]^+ \ F has tau(alpha) < Var_{K/F}(y) (relative budget-cost)
\\   RGA  some alpha in O_F[y]^+, generating K, tau < 3/2 nu(F), has no shape y'^2+cy'+b
\\   OK   survives

\\ ---- data for the four possible subfields --------------------------------
\\ [name, polynomial, nu(F), delta_F upper bound (deltaF.py), e]
\\ F is written in the variable y; K = F(beta) in the variable x.
{SUBF = [["Q(sqrt5)", y^2-y-1, 5/4, 1.3110, 3],
        ["Q(sqrt2)", y^2-2,   2,   1.7091, 3],
        ["Q(z7)+",  y^3-y^2-2*y+1, 14/9, 1.5513, 2],
        ["Q(z9)+",  y^3-3*y-1,     2,    2.1529, 2]];}

\\ largest nu with e*nu >= 2*(nu-delta)^2
nurelmax(e, dl) = my(s = (sqrt(e/2) + sqrt(e/2 + 4*dl)) / 2); s^2;

\\ ---- build the order O_F[y] inside K = F(y) -----------------------------
\\ nfF = nfinit(F) (variable y), rp = relative polynomial in x with coefficients in Z[y]
\\ returns [Ob, P, e, blocks, Grel, relbasisidx] or 0 if K not totally real / not a field
relorder(nfF, rp) =
{
  my(f = poldegree(nfF.pol), e = poldegree(rp, x), R, P, a, k, th, B, Ob, al, blocks, Grel, E, n);
  R = rnfequation(nfF, rp, 1);          \\ [P(x), a(x) = root of F-pol, k], root theta = beta + k*alpha
  P = R[1]; a = lift(R[2]); k = R[3];
  if (!polisirreducible(P), return(0));
  if (polsturm(P) != poldegree(P), return(0));
  \\ basis: omega_i * beta^j, with beta = x - k*a(x)
  my(bb = lift(Mod(x - k*a, P)));
  B = vector(f*e, t, my(i = (t-1) % f + 1, j = (t-1) \ f);
         lift(Mod(subst(nfF.zk[i], y, a) * bb^j, P)));
  Ob = ordinit(P, B);
  \\ blocks: group embeddings of K by the value of the F-generator a(theta_j)
  th = Ob_th(Ob);
  al = vector(#th, j, subst(a, x, th[j]));
  my(vals = vecsort(al, , 8), eps = 1e-20);
  \\ merge numerically equal values
  my(cl = List()); for (i = 1, #vals, if (#cl == 0 || abs(vals[i] - cl[#cl]) > eps, listput(cl, vals[i])));
  if (#cl != f, error("relorder: block detection failed"));
  blocks = vector(f, b, select(j -> abs(al[j] - cl[b]) < eps, [1..#th], 1));
  \\ Gram of  d * e * Var_{K/F}  on the classes of basis elements omega_i y^j, j >= 1 (exact after rounding)
  E = Ob_E(Ob); n = f*e;
  my(idx = [f+1..n]);
  Grel = matrix(#idx, #idx, s, t,
          my(u = E[, idx[s]], v = E[, idx[t]], acc = 0);
          for (b = 1, f, my(J = blocks[b], mu = sum(i = 1, #J, u[J[i]])/e, mv = sum(i = 1, #J, v[J[i]])/e);
               acc += sum(i = 1, #J, (u[J[i]] - mu) * (v[J[i]] - mv)));
          acc * e);
  my(Gr = round(Grel, &err));
  if (err > -40, error("relorder: relative Gram not integral"));
  [Ob, P, e, blocks, Gr, idx];
}

\\ relative variance of an element given by FULL coordinates c (length n, basis B, B[1] = 1)
relvar(RO, c) =
{
  my(Ob = RO[1], e = RO[3], blocks = RO[4], E = Ob_E(Ob), v = E*c, acc = 0, d = Ob_d(Ob));
  for (b = 1, #blocks, my(J = blocks[b], m = sum(i = 1, #J, v[J[i]])/e);
       acc += sum(i = 1, #J, (v[J[i]] - m)^2));
  acc / d;
}

\\ is the element with class coords c (basis B[2..n]) in F?  (all y-coordinates zero)
relinF(RO, c) = my(f = RO[6][1] - 1); c[f..#c] == vector(#c - f + 1)~;

\\ ---- the tests ------------------------------------------------------------
reltest(RO, nuF, vrel) =
{
  my(Ob = RO[1], f = RO[6][1] - 1, L, nr);
  \\ NSA: absolute nu of the order must be >= nu(F)
  if (onu(Ob) < nuF, return(["NSA", onu(Ob)]));
  \\ NSR: relative minimum of the order must equal vrel
  nr = qfminim(RO[5])[2] / (Ob_d(Ob) * RO[3]);
  if (nr < vrel, return(["NSR", nr]));
  \\ BCA / BCR
  L = otpos(Ob, max(nuF, vrel));
  for (i = 1, #L,
    my(t = L[i][1], c = L[i][2]);
    if (t < nuF, return(["BCA", t, oelt(Ob, c, L[i][3])]));
    if (t < vrel && !relinF(RO, c), return(["BCR", t, oelt(Ob, c, L[i][3])])));
  \\ RGA: absolute rigidity for alpha generating K
  L = otpos(Ob, 3/2*nuF);
  for (i = 1, #L,
    my(al = oelt(Ob, L[i][2], L[i][3]), s = sqshape(Ob_f(Ob), al, nuF));
    if (s === -1, next);
    if (#s == 0, return(["RGA", L[i][1], al])));
  ["OK", vrel];
}

\\ ---- enumeration, e = 2 (F cubic):  y^2 + b1 y + (b1^2 - D)/4,  D >> 0 -------
enum_e2(Fi, outfile) =
{
  my(nm = SUBF[Fi][1], nfF = nfinit(SUBF[Fi][2]), nuF = SUBF[Fi][3], dl = SUBF[Fi][4],
     numax = nurelmax(2, dl), Tmax = 12*numax, zk = nfF.zk, f = 3, T2, q, reps, cnt = Map(), pF = nfF.pol);
  \\ residues b1 mod 2 O_F
  reps = vector(8, t, my(v = binary(t - 1 + 8)[2..4]); sum(i = 1, 3, v[i]*zk[i]));
  \\ totally positive D with Tr D <= Tmax  (then T2(D) <= Tr(D)^2)
  T2 = matrix(3, 3, i, k, trace(Mod(zk[i]*zk[k], nfF.pol)));
  q = qfminim(T2, floor(Tmax^2))[3];
  my(Ds = List());
  for (k = 1, #q, for (s = 0, 1, my(D = (1 - 2*s) * sum(i = 1, 3, q[i,k]*zk[i]));
      if (istotpos(pF, D) && trace(Mod(D, pF)) <= Tmax, listput(Ds, D))));
  print(nm, ": nu_rel <= ", numax, ", Tr(Delta) <= ", Tmax, ", ", #Ds, " totally positive Delta");
  for (i = 1, #Ds,
    my(D = Ds[i]);
    if (#nfroots(nfF, x^2 - D), next);     \\ D a square: K would not be a field
    for (r = 1, #reps,
      my(b1 = reps[r], c0 = lift(Mod(b1^2 - D, pF) / 4));
      if (!nfeltisintegral(nfF, c0), next);
      my(rp = x^2 + b1*x + c0, RO = relorder(nfF, rp), vrel, res);
      if (RO === 0, next);
      vrel = trace(Mod(D, pF)) / 12;
      res = reltest(RO, nuF, vrel);
      mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
      write(outfile, [nm, D, b1, RO[2], res])));
  print(nm, ": ", Mat(cnt));
}

\\ helper: element of O_F integral?
nfeltisintegral(nf, a) = my(v = nfalgtobasis(nf, a)); denominator(v) == 1;

\\ ---- enumeration, e = 3 (F quadratic): y^3 + b1 y^2 + b2 y + b3 -------------
\\ within-block centred T2 in block k:  s_k = rho_k( (2/3) b1^2 - 2 b2 ) >= 0,
\\ Var_{K/F}(y) = (s_1 + s_2)/6 <= numax.
enum_e3(Fi, outfile) =
{
  my(nm = SUBF[Fi][1], nfF = nfinit(SUBF[Fi][2]), nuF = SUBF[Fi][3], dl = SUBF[Fi][4],
     numax = nurelmax(3, dl), zk = nfF.zk, pF = nfF.pol, th = vecsort(real(polroots(nfF.pol))), cnt = Map(),
     emb, reps, ncand = 0);
  emb = (a -> vector(2, k, subst(lift(Mod(a, pF)), y, th[k])));
  reps = vector(9, t, ((t-1) % 3) * zk[1] + ((t-1) \ 3) * zk[2]);
  print(nm, ": nu_rel <= ", numax);
  for (r = 1, #reps,
    my(b1 = reps[r], e1 = emb(b1), S = 6*numax);
    \\ b2 with p_k := rho_k(b2) - rho_k(b1)^2/3 <= 0 and  -2(p_1 + p_2) <= S
    \\ i.e. rho_k(b2) in [rho_k(b1)^2/3 - S/2, rho_k(b1)^2/3]
    my(lo = vector(2, k, e1[k]^2/3 - S/2), hi = vector(2, k, e1[k]^2/3), M, Mi, box);
    M = matrix(2, 2, k, i, subst(zk[i], y, th[k])); Mi = M^-1;
    \\ coefficient box from the embedding box
    box = vector(2, i, my(vs = [Mi[i,1]*lo[1] + Mi[i,2]*lo[2], Mi[i,1]*lo[1] + Mi[i,2]*hi[2],
                                Mi[i,1]*hi[1] + Mi[i,2]*lo[2], Mi[i,1]*hi[1] + Mi[i,2]*hi[2]]);
                        [floor(vecmin(vs)) - 1, ceil(vecmax(vs)) + 1]);
    forvec (u = box,
      my(b2 = u[1]*zk[1] + u[2]*zk[2], e2 = emb(b2), s);
      if (e2[1] > hi[1] + 1e-20 || e2[2] > hi[2] + 1e-20, next);
      s = vector(2, k, (2/3)*e1[k]^2 - 2*e2[k]);
      if (s[1] + s[2] > S + 1e-20, next);
      \\ b3: the cubic t^3 + b1 t^2 + b2 t + b3 must be real-rooted in both embeddings:
      \\ -b3 lies between the local max and min of t^3 + b1 t^2 + b2 t
      my(lo3 = vector(2), hi3 = vector(2), ok = 1);
      for (k = 1, 2,
        my(dq = 3*'t^2 + 2*e1[k]*'t + e2[k], rr, g = 't^3 + e1[k]*'t^2 + e2[k]*'t, v);
        if (poldisc(dq) < 0, ok = 0; break);
        rr = polroots(dq);
        rr = vecsort(real(rr));
        v = [subst(g, 't, rr[1]), subst(g, 't, rr[2])];     \\ local max, local min
        \\ need g(rr1) + b3 >= 0 and g(rr2) + b3 <= 0
        lo3[k] = -v[1]; hi3[k] = -v[2]);
      if (!ok, next);
      my(box3 = vector(2, i, my(vs = [Mi[i,1]*lo3[1] + Mi[i,2]*lo3[2], Mi[i,1]*lo3[1] + Mi[i,2]*hi3[2],
                                      Mi[i,1]*hi3[1] + Mi[i,2]*lo3[2], Mi[i,1]*hi3[1] + Mi[i,2]*hi3[2]]);
                              [floor(vecmin(vs)) - 1, ceil(vecmax(vs)) + 1]));
      forvec (w = box3,
        my(b3 = w[1]*zk[1] + w[2]*zk[2], e3 = emb(b3));
        if (e3[1] < lo3[1] - 1e-20 || e3[1] > hi3[1] + 1e-20 || e3[2] < lo3[2] - 1e-20 || e3[2] > hi3[2] + 1e-20, next);
        ncand++;
        my(rp = x^3 + b1*x^2 + b2*x + b3, RO, vrel, res);
        if (#nfroots(nfF, rp), next);          \\ reducible over F
        RO = relorder(nfF, rp);
        if (RO === 0, next);
        vrel = trace(Mod((2/3)*b1^2 - 2*b2, pF)) / 6;   \\ exact Var_{K/F}(beta)
        res = reltest(RO, nuF, vrel);
        mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
        write(outfile, [nm, rp, RO[2], res]))));
  print(nm, ": ", ncand, " relative cubics in the box; ", Mat(cnt));
}

\\ ---- general relative quadratic K = F(sqrt D), any base degree f (d = 2f) -------
\\ Var_{K/F}(beta) = Tr_F(D) / (2d) for beta = (-b1 + sqrt D)/2.
enum_e2g(nm, Fpol, nuF, dl, outfile) =
{
  my(nfF = nfinit(Fpol), f = poldegree(Fpol), d = 2*f, numax = nurelmax(2, dl), Tmax = 2*d*numax,
     zk = nfF.zk, pF = nfF.pol, T2, q, reps, cnt = Map(), Ds = List());
  reps = vector(2^f, t, my(v = binary(t - 1 + 2^f)[2..f+1]); sum(i = 1, f, v[i]*zk[i]));
  T2 = matrix(f, f, i, k, trace(Mod(zk[i]*zk[k], pF)));
  q = qfminim(T2, floor(Tmax^2))[3];
  for (k = 1, #q, for (s = 0, 1, my(D = (1 - 2*s) * sum(i = 1, f, q[i,k]*zk[i]));
      if (istotpos(pF, D) && trace(Mod(D, pF)) <= Tmax, listput(Ds, D))));
  print(nm, " (d=", d, "): nu_rel <= ", numax, ", Tr(Delta) <= ", Tmax, ", ", #Ds, " totally positive Delta");
  for (i = 1, #Ds,
    my(D = Ds[i]);
    if (#nfroots(nfF, x^2 - D), next);
    for (r = 1, #reps,
      my(b1 = reps[r], c0 = lift(Mod(b1^2 - D, pF) / 4));
      if (!nfeltisintegral(nfF, c0), next);
      my(rp = x^2 + b1*x + c0, RO = relorder(nfF, rp), vrel, res);
      if (RO === 0, next);
      vrel = trace(Mod(D, pF)) / (2*d);
      res = reltest(RO, nuF, vrel);
      mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
      write(outfile, [nm, D, b1, RO[2], res])));
  print(nm, ": ", Mat(cnt));
}
