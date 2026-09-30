\\ certN8.gp -- independent certificates for the degree-8 Case B exclusions (requires lift6.gp,
\\ field6.gp, condN.gp).  Every exclusion by fieldtest2 (RG2, BCR, RGR2) names an element alpha;
\\ since budget-cost and (relative, sharpened) rigidity are consequences of condition (N),
\\ (N) must fail at alpha.  condN decides (N) at alpha directly.
\\ Output per field: [P, status, "N-FAILS" | "N-HOLDS" (= contradiction) | "UNDECIDED", ms].
certN(infile, outfile, part = 0, nparts = 1) =
{
  my(L = readvec(infile), cnt = Map());
  for (i = 1, #L,
    if ((i - 1) % nparts != part, next);
    my(P = L[i][1], r = L[i][#L[i]], al, NI, M, t = getabstime(), v);
    if (r[1] == "OK", next);
    al = r[#r];
    NI = Ninit(ordinit(P, nfinit(P).zk));
    M = iferr(condN(NI, al), E, "ERR");
    v = if (M === 0, "N-FAILS", if (M === "ERR", "UNDECIDED", "N-HOLDS"));
    mapput(cnt, v, if (mapisdefined(cnt, v), mapget(cnt, v), 0) + 1);
    write(outfile, [P, r[1], v, getabstime() - t]));
  Mat(cnt);
}
