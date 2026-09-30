\\ collect8.gp -- degree 8, Case B over Q(sqrt5): gather the order-level survivors ("OK")
\\ of enum_q4 / enum_e2g, pass to the maximal order, dedupe, and run fieldtest2.
\\ Requires lift6.gp, rel6.gp, field6.gp, sharp.gp.
collect8(files, outfile, part = 0, nparts = 1) =
{
  my(F = Map(), M, cnt = Map());
  for (i = 1, #files,
    my(L = readvec(files[i]));
    for (j = 1, #L, my(e = L[j], res = e[#e], P = e[#e - 1]);
      if (res[1] == "OK", mapput(F, polredabs(P), 1))));
  M = if (#F, Mat(F)[,1], []);
  print("distinct fields among order-level survivors: ", #M);
  for (i = 1, #M,
    if ((i - 1) % nparts != part, next);
    my(r = fieldtest2(M[i]));
    mapput(cnt, r[1], if (mapisdefined(cnt, r[1]), mapget(cnt, r[1]), 0) + 1);
    write(outfile, [M[i], nfdisc(M[i]), r]));
  print(Mat(cnt));
}
