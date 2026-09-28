\\ condN.gp -- exact decision procedure for condition (N) (requires lift6.gp, field6.gp).
\\
\\ (N)(alpha):  there is a positive semidefinite symmetric M with M_ii in Z, 2 M_ij in Z and
\\              sum_{i,j} M_ij w_i w_j = alpha,  (w_1..w_d) the integral basis of O_K.
\\ If K has a universal Z-form, (N)(alpha) holds for every alpha in O_K^+
\\ ([KY1, Prop. 3.1]: M = V^t M_Q V).  (N) is invariant under alpha -> eps^2 alpha
\\ and under automorphisms of K.
\\
\\ Bounds used (all rigorous):
\\  * <M, T> = Tr(alpha)  with T the trace Gram, so tr(M) <= Tr(alpha) / lambda_min(T);
\\  * M = E^{-1} H E^{-t} with H psd, diag(H) = sigma(alpha)  (E = embedding matrix), hence
\\    M_ii <= ( sum_j |sigma_j(w_i^*)| sqrt(sigma_j(alpha)) )^2  (w^* = trace-dual basis);
\\  * a zero diagonal entry forces its row to vanish;
\\  * for psd M with diagonal a (support of size s):  sum_{i<j} (2M_ij)^2/(4 a_i a_j) <= s(s-1)/2.
\\ For each support S and diagonal a, the equation is linear in b_ij = 2 M_ij (i<j in S);
\\ its integer solutions form a coset b0 + Ker, enumerated by Fincke-Pohst, and every
\\ candidate is tested for positive semidefiniteness exactly (char. polynomial signs).

\\ exact psd test for a rational symmetric matrix
ispsd(M) =
{
  my(n = #M, p = charpoly(M), c);
  \\ eigenvalues all >= 0  <=>  coefficients of det(tI - M) alternate in sign (weakly)
  for (k = 0, n, c = polcoeff(p, n - k); if (c != 0 && sign(c) != (-1)^k, return(0)));
  1;
}

\\ precomputation for a field: products w_i w_j in coordinates, dual-basis embeddings
Ninit(Ob) =
{
  my(d = Ob_d(Ob), B = Ob_B(Ob), P, T = Ob_T(Ob), U, lam);
  P = matrix(d, d, i, j, kcoords(Ob, B[i]*B[j]));
  U = Ob_E(Ob) * T^-1;                       \\ U[j,i] = sigma_j(w_i^*)
  lam = vecmin(real(mateigen(T * 1.0, 1)[1]));
  [Ob, P, U, lam];
}

\\ all integer vectors z with (b0 + K z)^t W (b0 + K z) <= R, W = diag(w) (w > 0 rational)
cosetFP(b0, K, w, R) =
{
  my(p = #b0, k = #K, L, Wd, G, res = List(), c, q);
  if (k == 0, if (sum(i = 1, p, w[i]*b0[i]^2) <= R, return([b0]), return([])));
  L = denominator(w); Wd = matdiagonal(w * L);
  c = L * R + 1;
  \\ homogenised form on (z, t):  |t b0 + K z|_W^2 + c t^2
  my(A = matconcat([K, b0]));
  G = A~ * Wd * A; G[k+1, k+1] += c;
  q = qfminim(G, floor(L*R + c))[3];
  for (i = 1, #q,
    my(t = q[k+1, i]);
    if (abs(t) != 1, next);
    my(z = t * q[1..k, i], b = b0 + K*z);
    if (sum(j = 1, p, w[j]*b[j]^2) <= R, listput(res, b)));
  Vec(res);
}

\\ decide (N) for alpha (polynomial in x).  Returns a psd matrix M or 0.
condN(NI, al) =
{
  my(Ob = NI[1], P = NI[2], U = NI[3], lam = NI[4], d = Ob_d(Ob), f = Ob_f(Ob),
     av = kcoords(Ob, al), ta, Bt, sa, bnd, emb);
  ta = trace(Mod(al, f));
  Bt = floor(ta / (lam * (1 - 1e-30)) + 1e-30);
  emb = vector(d, j, subst(lift(Mod(al, f)), x, Ob_th(Ob)[j]));
  if (vecmin(emb) <= 0, error("condN: alpha not totally positive"));
  sa = vector(d, j, sqrt(emb[j]));
  bnd = vector(d, i, min(Bt, floor(sum(j = 1, d, abs(U[j,i]) * sa[j])^2 * (1 + 1e-30) + 1e-30)));
  forsubset (d,
    S0,
    my(S = Vec(S0), s = #S);
    if (s == 0, next);
    if (vecsum(vector(s, t, 1)) > Bt, next);
    \\ diagonal entries a_t >= 1, a_t <= bnd[S[t]], sum <= Bt
    forvec (a = vector(s, t, [1, bnd[S[t]]]),
      if (vecsum(a) > Bt, next);
      my(rhs = av - sum(t = 1, s, a[t] * P[S[t], S[t]]), pairs = List());
      for (u = 1, s, for (v = u+1, s, listput(pairs, [u, v])));
      pairs = Vec(pairs);
      if (#pairs == 0,
        if (rhs == 0, my(Mf = matrix(d, d)); for (u = 1, s, Mf[S[u], S[u]] = a[u]); return(Mf));
        next);
      my(A = matrix(d, #pairs, r, c, P[S[pairs[c][1]], S[pairs[c][2]]][r]), sol, K, w, cand);
      sol = matsolvemod(A, 0, rhs);
      if (sol === 0, next);
      K = matkerint(A);
      w = vector(#pairs, c, 1 / (4 * a[pairs[c][1]] * a[pairs[c][2]]));
      cand = cosetFP(sol, K, w, #pairs);
      for (ci = 1, #cand,
        my(b = cand[ci], M = matdiagonal(a) * 1);
        if (vecmax(vector(#b, c, b[c]^2 - 4*a[pairs[c][1]]*a[pairs[c][2]])) > 0, next);
        for (c = 1, #pairs, M[pairs[c][1], pairs[c][2]] = b[c]/2; M[pairs[c][2], pairs[c][1]] = b[c]/2);
        if (ispsd(M),
          \\ return the full d x d matrix
          my(Mf = matrix(d, d));
          for (u = 1, s, for (v = 1, s, Mf[S[u], S[v]] = M[u, v]));
          return(Mf)))));
  0;
}

\\ verify a returned certificate exactly
checkN(NI, al, Mf) =
{
  my(Ob = NI[1], d = Ob_d(Ob), B = Ob_B(Ob));
  if (!ispsd(Mf), return(0));
  for (i = 1, d, if (denominator(Mf[i,i]) != 1, return(0)); for (j = 1, d, if (denominator(2*Mf[i,j]) != 1, return(0))));
  lift(Mod(sum(i = 1, d, sum(j = 1, d, Mf[i,j]*B[i]*B[j])) - al, Ob_f(Ob))) == 0;
}
