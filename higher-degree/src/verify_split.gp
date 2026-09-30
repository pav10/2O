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

\\ vsharp: exact re-check of every exclusion printed by the split enumerator in mode whysharp
\\ (modes X, YL, YU).  Lines: "NS 0 a..", "RIGM n a..", "RIGP n a.." (x-polynomial) and
\\ "YBL n g..", "YBU n g.." (y-polynomial; x = n + y^2 resp. n - y^2).  For irreducible f_x:
\\   NS   : the exact minimum of Var on Z[x_1] \ Z is < V = Var(x_1);
\\   RIG* : alpha = x - n (n - x) has tau < 3/2 V and sqshape2(f, alpha, V) = [] (sharp rigidity);
\\   YB*  : y is in O_K (K = Q(x) = Q(y)) and a downspread of y is < V (budget-cost fails).
\\ Returns [checked, reducible skipped, disagreements, counts by tag].
vsharp(infile, part, nparts) =
{
  my(L = externstr(concat("cat ", infile)), chk = 0, red = 0, bad = 0, cnt = Map());
  for (i = 1, #L,
    if ((i - 1) % nparts != part, next);
    my(w = strsplit(L[i], " "), tag = w[1], n = eval(w[2]), c = apply(eval, w[3..#w]), f, d, V, ok);
    if (tag == "OUT", next);
    if (tag == "YBL" || tag == "YBU",
      my(g = Pol(concat([1], c)), ys, ty);
      f = charpoly(Mod(if (tag == "YBL", n + x^2, n - x^2), g));
      if (!polisirreducible(f), red++; next);
      d = poldegree(f); V = ((d-1)*polcoeff(f, d-1)^2/d - 2*polcoeff(f, d-2)) / d;
      ys = vecsort(real(polroots(g))); ty = -polcoeff(g, d-1)/d;
      ok = (ty - (ceil(ys[1]) - 1) < V) || (floor(ys[d]) + 1 - ty < V),
    \\ else: x-polynomial
      f = Pol(concat([1], c)); d = poldegree(f);
      if (!polisirreducible(f), red++; next);
      V = ((d-1)*c[1]^2/d - 2*c[2]) / d;
      if (tag == "NS", ok = (onu(ordZx(f)) < V),
        my(al = if (tag == "RIGM", x - n, n - x), ta = trace(Mod(al, f)) / d, s);
        s = sqshape2(f, al, V);
        ok = (ta < 3/2*V) && !(s === -1) && (#s == 0)));
    chk++;
    mapput(cnt, tag, if (mapisdefined(cnt, tag), mapget(cnt, tag), 0) + 1);
    if (!ok, bad++; print("DISAGREE: ", L[i])));
  [chk, red, bad, Mat(cnt)];
}
