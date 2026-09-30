\\ verify_split.gp -- for an irreducible polynomial kept by the plain enumerator but not by
\\ the split enumeration (modes X/YL/YU), certify the exclusion exactly:
\\ some forced side (tau_+- < min(3/2 V, V+3/4)) has NO square root in O_K,
\\ or its square root y has a downspread < V (then y - m or m' - y violates budget-cost).
explain(f) =
{
  my(d = poldegree(f), a = Vec(f), tau = -a[2]/d, V, th, m, nm, np, nz = nfinit(subst(f, x, 'w)), res = List());
  V = ((d-1)*a[2]^2/d - 2*a[3]) / d;
  th = vecsort(real(polroots(f))); m = min(3/2*V, V + 3/4);
  nm = ceil(th[1]) - 1; np = floor(th[d]) + 1;
  foreach ([[tau - nm, 'w - nm], [np - tau, np - 'w]], side,
    if (side[1] >= m, next);
    my(r = nfroots(nz, x^2 - side[2]));
    if (#r == 0, return(["no-sqrt"]));
    my(g = charpoly(Mod(lift(r[1]), nz.pol), x), ys = vecsort(real(polroots(g))), ty = vecsum(ys)/d);
    if (ty - (ceil(ys[1]) - 1) < V || floor(ys[d]) + 1 - ty < V, return(["sqrt-fails-budget"])));
  ["UNEXPLAINED"];
}
checkfile(file) =
{
  my(L = externstr(concat("cat ", file)), cnt = Map());
  for (i = 1, #L,
    my(v = eval(concat(["[", strjoin(strsplit(L[i], " "), ","), "]"])), f = Pol(concat([1], v)), r);
    if (!polisirreducible(f), r = "reducible", r = explain(f)[1]);
    mapput(cnt, r, if (mapisdefined(cnt, r), mapget(cnt, r), 0) + 1));
  Mat(cnt);
}
