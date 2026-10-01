\\ sqforce8.gp -- degree 8, Case B: square forcing (Proposition "square dichotomy").
\\ Requires lift6.gp, order_filter.gp, field6.gp, sharp.gp, pipeline.gp.
\\
\\ If x_1 lies in a subfield L of K, then nu(K) = nu(L) = nu, and every totally positive
\\ alpha in O_L \ Z with tau(alpha) < min(3/2 nu, nu + 3/4) is a square in O_K.  So a window
\\ non-square alpha of L forces L(sqrt alpha) in K.  forceup iterates this up to degree dK.

issq(f, a) = my(nz = nfinit(subst(f, x, 'w))); #nfroots(nz, x^2 - subst(lift(Mod(a, f)), x, 'w)) > 0;

\\ returns ["free", L] (no window non-square in L; L needs an enumeration),
\\ ["impossible", reason] or ["forced", K, fieldtest2(K)]
forceup(f, nu, dK) =
{
  my(nf = nfinit(f), Ob = ordinit(f, nf.zk), L = otpos(Ob, min(3/2*nu, nu + 3/4)), ns = List(), Kp);
  if (onu(Ob) < nu, return(["impossible", "nu(L) < nu"]));
  for (j = 1, #L, my(al = oelt(Ob, L[j][2], L[j][3])); if (!issq(f, al), listput(ns, al)));
  if (#ns == 0, return(if (poldegree(f) == dK, ["forced", f, fieldtest2(f)], ["free", f])));
  \\ L(sqrt alpha) has degree 2 [L:Q], which must divide [K:Q]
  if (dK % (2*poldegree(f)) != 0, return(["impossible", if (poldegree(f) == dK, "window non-square in K", "window non-square, [K:L] odd")]));
  Kp = polredabs(rnfequation(nfinit(subst(f, x, y)), x^2 - subst(lift(Mod(ns[1], f)), x, y)));
  if (polsturm(Kp) != poldegree(Kp), return(["impossible", "forced extension not totally real"]));
  forceup(Kp, nu, dK);
}

\\ subfields F = Q(x_1) (degree 2 or 4) passing the tests valid in a subfield (NS, BC); the
\\ downspread bound only involves the conjugates of x_1, so the degree-[F:Q] enumeration applies
subF(infile, dK = 8) =
{
  my(P = select(p -> polisirreducible(p), readenum(infile)), M = Map(), C);
  \\ x_1 realises nu(F), so keep the least Var over all generators of the same field
  for (i = 1, #P, my(r = filtNoRG(P[i]), g);
    if (r[1] == "OK", g = polredabs(P[i]);
      mapput(M, g, if (mapisdefined(M, g), min(mapget(M, g), r[2]), r[2]))));
  C = Mat(M);
  print(#C~, " fields F with x_1 in F possible (NS/BC):");
  for (i = 1, #C~, my(r = forceup(C[i,1], C[i,2], dK));
    print("  ", C[i,1], "  d=", nfdisc(C[i,1]), "  nu=", C[i,2], "  -> ", r[1], "  ",
          if (r[1] == "forced", [r[2], nfdisc(r[2]), r[3][1]], r[2])));
}
