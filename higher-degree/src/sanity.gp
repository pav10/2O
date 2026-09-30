\\ run the order-level pipeline on the enumerator output for a given degree
sanity(infile) =
{
  my(L = externstr(concat("cat ", infile)), P, irr, res = Map(), surv = List());
  P = apply(s -> Pol(concat([1], eval(concat(["[", strjoin(strsplit(s, " "), ","), "]"])))), L);
  irr = select(p -> polisirreducible(p), P);
  for (i = 1, #irr,
    my(r = filt(irr[i]));
    mapput(res, r[1], if (mapisdefined(res, r[1]), mapget(res, r[1]), 0) + 1);
    if (r[1] == "OK", listput(surv, [irr[i], polredabs(irr[i]), nfdisc(irr[i])])));
  print(infile, ": ", #P, " polys, ", #irr, " irreducible; ", Mat(res));
  print("  survivors: ", Vec(surv));
}
