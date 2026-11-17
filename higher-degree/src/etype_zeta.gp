\\ etype_zeta.gp -- Siegel's formula by ideal classes and the trace-minimal codifferent condition.
\\ Companion to notes/etype-zeta.tex.  PARI/GP >= 2.15.  Run scripts with stdin closed:  gp -q file.gp < /dev/null
\\
\\ siegelb(w)        Siegel's coefficients b_1(w),...,b_r(w), r = dim M_w(SL_2(Z)):  a_0(f) = sum_l b_l a_l(f)
\\ traceL(nf, l)     totally positive elements of the codifferent of trace l (coordinates on nf.zk); floating-point
\\                   enumeration (qfminim, margin 1e-9), each element then checked exactly
\\ partialzeta2(bnf) Map: wide class (vector mod cyc) -> zeta(2; C), via Hecke L-values of class characters
\\ csiegel(P)        both sides of (2.1) of the notes for every wide class C
\\ condP(P)          the condition (1.2) of the notes for m = 1 (trace-one elements)
\\ certify(P)        an exactly checked certificate [nu, alpha, N(alpha)] that (1.2) fails, or 0
\\ smallideal(P)     the norm of xi*diff for a shortest nonzero trace-zero xi in the codifferent (Lemma 4.1)
\\ etad, gam, Gd, eps, cnt   the constants eta_d, gamma_n, G_d, epsilon_{d,l}, c_{d,l} of Lemmas 4.1 and 4.2
default(parisizemax, 2*10^9);
default(threadsizemax, 2*10^9);
default(realprecision, 38);

siegelb(w) =
{
  my(mf = mfinit([1, w], 4), B = mfbasis(mf), r = #B, M, y);
  M = matrix(#B, r, j, l, mfcoefs(B[j], r)[l + 1]);
  y = vector(#B, j, mfcoefs(B[j], 0)[1])~;
  matsolve(M, y);
}

traceL(nf, l) =
{
  my(d = poldegree(nf.pol), zk = nf.zk, Tm = matrix(d, d, i, j, nfelttrace(nf, nfeltmul(nf, zk[i], zk[j]))),
     Bd = Tm^-1, r = real(polroots(nf.pol)), E, G, q, res = List());
  E = matrix(d, d, k, i, subst(lift(zk[i]), variable(nf.pol), r[k]));
  G = Bd~ * (E~ * E) * Bd;                         \\ T2 Gram matrix on the dual basis
  q = qfminim(G, l^2 + 1e-9, , 2)[3];              \\ Tr nu = l and nu >> 0 force T2(nu) <= l^2
  for (k = 1, #q, for (s = 0, 1,
    my(c = if (s, -q[, k], q[, k]), g = Bd * c);
    if (nfelttrace(nf, g) == l && nfeltsign(nf, g) == vector(d, i, 1), listput(res, g))));
  Vec(res);
}

cls(bnf, J) = my(e = bnfisprincipal(bnf, J, 0)); vector(#e, i, e[i] % bnf.cyc[i]);

partialzeta2(bnf) =
{
  my(bnr = bnrinit(bnf, 1), cyc = bnr.cyc, h = vecprod(cyc), chars = List(), Lv = List(), classes = List(), res = Map());
  forvec(e = vector(#cyc, i, [0, cyc[i] - 1]), listput(classes, e));
  forvec(c = vector(#cyc, i, [0, cyc[i] - 1]), listput(chars, c); listput(Lv, lfun(lfuncreate([bnr, c]), 2)));
  for (k = 1, #classes,
    my(C = classes[k], z = sum(j = 1, #chars, conj(exp(2*Pi*I*sum(i = 1, #cyc, chars[j][i] * C[i] / cyc[i]))) * Lv[j]) / h);
    mapput(res, C, real(z)));
  res;
}

csiegel(P) =
{
  my(bnf = bnfinit(P, 1), nf = bnf.nf, d = poldegree(P), b = siegelb(2 * d), r = #b, cyc = bnf.cyc,
     dif = nf.diff, cD = cls(bnf, dif), Z2 = partialzeta2(bnf), A = Map(), classes = List(), Dk = abs(nf.disc));
  forvec(e = vector(#cyc, i, [0, cyc[i] - 1]), listput(classes, e); mapput(A, e, 0));
  for (l = 1, r,
    my(T = traceL(nf, l));
    for (t = 1, #T,
      my(Ai = idealmul(nf, T[t], dif), F = idealfactor(nf, Ai), np = matsize(F)[1]);
      if (np == 0, my(e = cls(bnf, idealhnf(nf, 1))); mapput(A, e, mapget(A, e) + b[l]); next);
      forvec(ex = vector(np, i, [0, F[i, 2]]),
        my(J = idealfactorback(nf, F[, 1], ex~), e = cls(bnf, J));
        mapput(A, e, mapget(A, e) + b[l] * idealnorm(nf, J)))));
  print(P, "  disc=", nf.disc, "  Cl=", cyc, "  [diff]=", cD, "  b=", b~);
  for (k = 1, #classes,
    my(C = classes[k], Ct = vector(#cyc, i, (-C[i] + cD[i]) % cyc[i]),
       rhs = (-1)^d * Dk^(3/2) * mapget(Z2, Ct) / (4 * Pi^2)^d);
    printf("  C=%s  LHS=%s  RHS=%.30g  |LHS-RHS|=%.3g\n", C, mapget(A, C), rhs, abs(mapget(A, C) - rhs)));
}

condP(P) =
{
  my(nf = nfinit(P), T = traceL(nf, 1), D = abs(nf.disc), Ns, bnf, dif);
  if (#T == 0, return([0, "EMPTY"]));
  Ns = vector(#T, t, abs(nfeltnorm(nf, T[t])) * D);
  if (vecmin(Ns) == 1, return([#T, if (vecmax(Ns) == 1, "PASS-principal", "FAIL-principal")]));
  bnf = bnfinit(nf); dif = nf.diff;
  if (bnfisprincipal(bnf, dif, 0) == 0 * bnf.cyc~, return([#T, "FAIL-principal"]));
  for (t = 1, #T,
    my(A = idealmul(nf, T[t], dif), F = idealfactor(nf, A), np = matsize(F)[1]);
    forvec(ex = vector(np, i, [0, F[i, 2]]),
      if (ex == vector(np), next);
      if (bnfisprincipal(bnf, idealfactorback(nf, F[, 1], ex~), 0) == 0 * bnf.cyc~, return([#T, "FAIL-nonprincipal", t]))));
  [#T, "PASS-nonprincipal"];
}

\\ fraction of trace-one elements that individually satisfy (1.2); used for non-principal different
condPcount(bnf) =
{
  my(nf = bnf.nf, T = traceL(nf, 1), dif = nf.diff, nfail = 0);
  for (t = 1, #T,
    my(A = idealmul(nf, T[t], dif), F = idealfactor(nf, A), np = matsize(F)[1], hit = 0);
    if (np == 0, next);
    forvec(ex = vector(np, i, [0, F[i, 2]]),
      if (ex == vector(np), next);
      if (bnfisprincipal(bnf, idealfactorback(nf, F[, 1], ex~), 0) == 0 * bnf.cyc~, hit = 1; break));
    if (hit, nfail++));
  [#T, nfail];
}

isindual(nf, x) = for (i = 1, poldegree(nf.pol), if (denominator(nfelttrace(nf, nfeltmul(nf, x, nf.zk[i]))) != 1, return(0))); 1;

certify(P) =
{
  my(bnf = bnfinit(P, 1), nf = bnf.nf, T = traceL(nf, 1), dif = nf.diff, d = poldegree(P));
  for (t = 1, #T,
    my(g = T[t], A = idealmul(nf, g, dif), F = idealfactor(nf, A), np = matsize(F)[1]);
    if (np == 0, next);
    forvec(ex = vector(np, i, [0, F[i, 2]]),
      if (ex == vector(np), next);
      my(r = bnfisprincipal(bnf, idealfactorback(nf, F[, 1], ex~)));
      if (r[1] == 0 * r[1],
        my(al = r[2], de = nfeltdiv(nf, g, al));
        if (nfelttrace(nf, g) == 1 && nfeltsign(nf, g) == vector(d, i, 1) && isindual(nf, g)
            && denominator(content(al)) == 1 && abs(nfeltnorm(nf, al)) > 1 && isindual(nf, de),
          return([lift(nf.zk * g), lift(nf.zk * al), nfeltnorm(nf, al)])))));
  0;
}

smallideal(P) =
{
  my(nf = nfinit(P), d = poldegree(P), D = abs(nf.disc), zk = nf.zk, r = real(polroots(P)),
     E = matrix(d, d, k, i, subst(lift(zk[i]), variable(P), r[k])),
     Tm = matrix(d, d, i, j, nfelttrace(nf, nfeltmul(nf, zk[i], zk[j]))), Bd = Tm^-1, K, Gm, sv, xi);
  K = matkerint(Mat(vector(d, i, nfelttrace(nf, Bd[, i]))));
  Gm = K~ * Bd~ * (E~ * E) * Bd * K;
  sv = qfminim(Gm, , 1, 2); xi = Bd * K * sv[3][, 1];
  abs(nfeltnorm(nf, xi)) * D;
}

\\ constants of Lemmas 4.1 and 4.2
etad(d) = { my(best = 0.); for (k = 1, d - 1, my(a = d - k, b = k, s = sqrt(k * a^2 + (d - k) * b^2)); best = max(best, (a / s)^k * (b / s)^(d - k))); best; }
gam(n) = { my(ex = [1, 2/sqrt(3), 2^(1/3), sqrt(2), 8^(1/5), (64/3)^(1/6), 64^(1/7), 2]); if (n <= 8, ex[n], (2 / Pi) * gamma(2 + n / 2)^(2 / n)); }
Gd(d) = my(n = d - 1); etad(d) * (gam(n) * d^(1 / n))^(d / 2);
eps(d, l) = my(a = 2 * Pi^2 * l^2); a * intnum(V = 1/3, [oo, a], exp(-a * V) * ((1 + 2 * sqrt(3 * V))^(d - 1) - 1));
cnt(d, l) = (2 * Pi * exp(1) * l^2 / d)^((d - 1) / 2) * (1 + eps(d, l)) / sqrt(d);
