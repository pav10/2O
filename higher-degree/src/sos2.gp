\\ sos2.gp -- independent cross-check of (N) for small-trace alpha, via Mordell-Ko + E6.
\\ (requires lift6.gp, field6.gp)
\\
\\ If (N)(alpha) holds with Gram 2M of rank r, then 2M is the Gram of an even lattice L'
\\ of rank r <= 6.  Mordell-Ko: every positive definite integral lattice of rank <= 5, and
\\ every one of rank 6 other than E6, is a sum of squares of integral linear forms.
\\ Hence (N)(alpha) <=> (SoS2)  2 alpha = sum_k x_k^2,  x_k in O_K,  sum_k x_k in 2 O_K
\\                     or (E6)   2 alpha = Q_{E6}(v), v in E6 (x) O_K with coefficients a Z-basis
\\                               of O_K; then 2 Tr(alpha) >= 6 (3 d_K)^{1/6}  (AM-GM).
\\ (The parity condition: (U^t U)_ii = sum_k U_ki^2 == sum_k U_ki mod 2.)
\\ So for Tr(alpha) < 3 (3 d_K)^{1/6},  (N) <=> (SoS2).
\\
\\ This routine decides (SoS2) by a completely different search from condN.gp:
\\ all x with sigma_j(x)^2 <= 2 sigma_j(alpha) (up to sign), then a depth-first
\\ multiset search on the remainder 2 alpha - sum x_k^2 (must stay totally >= 0),
\\ and finally the parity condition.

sos2(P, al) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), d = poldegree(P), B = Ob_B(Ob), th = Ob_th(Ob),
     a2 = 2*al, ea, T = Ob_T(Ob), tr2, q, X = List(), found = 0, sols = List());
  ea = vector(d, j, subst(lift(Mod(a2, P)), x, th[j]));
  tr2 = trace(Mod(a2, P));
  \\ candidates: T2(x) = Tr(x^2) <= Tr(2 alpha)
  q = qfminim(T, floor(tr2))[3];
  for (k = 1, #q,
    my(c = q[,k], xe = Ob_E(Ob) * c, ok = 1);
    for (j = 1, d, if (xe[j]^2 > ea[j] + 1e-30, ok = 0; break));
    if (ok, listput(X, [c, sum(i = 1, d, c[i]*B[i])])));
  X = Vec(X);
  \\ squares as coordinate vectors, sorted by trace (descending) for the search
  my(Sq = vector(#X, i, kcoords(Ob, X[i][2]^2)), trs = vector(#X, i, trace(Mod(X[i][2]^2, P))),
     perm = vecsort(trs, , 5), target = kcoords(Ob, a2), E = Ob_E(Ob));
  X = vector(#X, i, X[perm[i]]); Sq = vector(#Sq, i, Sq[perm[i]]); trs = vector(#trs, i, trs[perm[i]]);
  \\ depth-first search over multisets (non-increasing index order); track sum of x's mod 2
  my(D = [X, Sq, trs, E, Ob_tr(Ob), d], r);
  r = sos2rec(D, 1, target, vector(d)~, []);
  found = (type(r) == "t_VEC");
  [found, if (found, apply(i -> X[i][2], r), 0), #X];
}

\\ recursive multiset search; D = [X, Sq, trs, E, tr, d].  Returns index list or 0.
sos2rec(D, start, rem, par, used) =
{
  my(X = D[1], Sq = D[2], trs = D[3], E = D[4], tr = D[5], d = D[6], remtr);
  if (rem == 0, return(if (par % 2 == 0, used, 0)));
  remtr = sum(i = 1, d, rem[i]*tr[i]);
  for (i = start, #X,
    if (trs[i] > remtr, next);
    my(r2 = rem - Sq[i], e = E * r2, ok = 1, res);
    for (j = 1, d, if (e[j] < -1e-30, ok = 0; break));
    if (!ok, next);
    res = sos2rec(D, i, r2, par + X[i][1], concat(used, [i]));
    if (type(res) == "t_VEC", return(res)));
  0;
}
\\ E6 trace threshold  3 (3 d_K)^(1/6)
e6bound(P) = 3 * (3 * abs(nfdisc(P)))^(1/6);
