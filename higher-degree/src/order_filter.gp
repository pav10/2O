\\ order_filter.gp -- order-level (Z[x_1]) necessary conditions for x_1 to be the
\\ shortest vector of Lambda_K for a sextic K = Q(x_1) with a universal Z-form.
\\ usage (from shell):  gp -q lift6.gp order_filter.gp  with INFILE, OUTFILE, PART, NPARTS set
\\ status codes:
\\   NS  a non-rational y in Z[x_1] has Var(y) < Var(x_1): x_1 is not a shortest vector
\\   BC  a non-rational totally positive alpha in Z[x_1] has tau(alpha) < Var(x_1)  (budget-cost)
\\   RG  a totally positive alpha in Z[x_1] with tau < 3/2 Var(x_1) has no shape y^2+cy+b (rigidity)
\\   OK  survives

filt(f) =
{
  my(Ob = ordZx(f), a = Vec(f), d = poldegree(f), V, nu, L, r);
  V = ((d-1)*a[2]^2/d - 2*a[3]) / d;
  nu = onu(Ob);
  if (nu < V, return(["NS", nu]));
  L = otpos(Ob, V);
  if (#L, return(["BC", L[1][1], oelt(Ob, L[1][2], L[1][3])]));
  L = otpos(Ob, 3/2*V);
  for (i = 1, #L,
    my(al = oelt(Ob, L[i][2], L[i][3]), s = sqshape(f, al, V));
    if (s === -1, next);
    if (#s == 0, return(["RG", L[i][1], al])));
  ["OK", V];
}

run(infile, outfile, part, nparts) =
{
  my(P = readvec(infile), cnt = Map());
  for (i = 1, #P,
    if ((i - 1) % nparts != part, next);
    my(f = P[i], r = filt(f));
    write(outfile, [f, r]));
}
