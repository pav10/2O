\\ search, by increasing trace, for alpha in O_K^+ failing (N)
findwitness(P, Tmax) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), NI = Ninit(Ob), L, t0, n = 0);
  L = otpos(Ob, Tmax);
  for (i = 1, #L,
    my(al = oelt(Ob, L[i][2], L[i][3]), M);
    t0 = getabstime();
    M = condN(NI, al);
    n++;
    if (M === 0,
      return([P, "WITNESS", al, "Tr", trace(Mod(al, P)), "N", norm(Mod(al, P)), "minpoly", minpoly(Mod(al, P)), "tested", n, "ms", getabstime() - t0])));
  [P, "none below tau", Tmax, "tested", n];
}

\\ as findwitness, but an alpha on which condN cannot finish (stack) is recorded and skipped
findwitness2(P, Tmax) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), NI = Ninit(Ob), L, n = 0, skipped = List());
  L = otpos(Ob, Tmax);
  for (i = 1, #L,
    my(al = oelt(Ob, L[i][2], L[i][3]), M);
    M = iferr(condN(NI, al), E, listput(skipped, [trace(Mod(al, P)), minpoly(Mod(al, P))]); 1);
    n++;
    if (M === 0,
      return([P, "WITNESS", al, "Tr", trace(Mod(al, P)), "minpoly", minpoly(Mod(al, P)), "tested", n, "skipped", Vec(skipped)])));
  [P, "none below tau", Tmax, "tested", n, "skipped", Vec(skipped)];
}
