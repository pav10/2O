\\ enum_check.gp -- independent exact re-implementation of the Rolle/interlacing
\\ enumeration (no filters), for validating enum6.c.  Candidate ranges for a_m are
\\ obtained from floating roots of g_{m-1} and then widened by 2 on each side;
\\ membership is decided EXACTLY by counting real roots of g_m (with multiplicity).
nreal(g) = my(h = g / gcd(g, g')); if (poldegree(h) == 0, 0, polsturm(h)) + (poldegree(g) - poldegree(h)) * 0;
allreal(g) = {my(f = g, n = poldegree(g), c = 0);
  \\ count real roots with multiplicity via squarefree decomposition
  my(F = factor(f)); for (i = 1, #F~, if (poldegree(F[i,1]) > 0, c += F[i,2] * polsturm(F[i,1]))); c == n;}
gm(a, m, N) = sum(i = 0, m, a[i+1] * binomial(N - i, m - i) * 'x^(m - i));
encheck(N, a1, a2) =
{
  my(res = List(), S = (N-1)*a1^2/N - 2*a2);
  if (S <= 0, return([]));
  encheck_rec(N, 3, [1, a1, a2], ~res);
  Vec(res);
}
encheck_rec(N, m, a, ~res) =
{
  if (m > N, listput(res, a[2..#a]); return);
  my(gprev = gm(a, m-1, N), r = vecsort(real(polroots(gprev))), h = gm(concat(a, [0]), m, N), vals, lo, hi);
  vals = vector(#r, i, -subst(h, 'x, r[i]));
  lo = -oo; hi = oo;
  for (i = 1, #r, my(k = m - i); if (k % 2 == 0, lo = max(lo, vals[i]), hi = min(hi, vals[i])));
  if (#r == 0, error("m>=2 expected"));
  for (am = floor(lo) - 2, ceil(hi) + 2,
    my(b = concat(a, [am]));
    if (allreal(gm(b, m, N)), encheck_rec(N, m + 1, b, ~res)));
}
