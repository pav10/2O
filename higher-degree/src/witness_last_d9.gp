default(parisizemax, 6*10^9);
P = x^9 - 4*x^8 - 7*x^7 + 34*x^6 + 15*x^5 - 85*x^4 - 13*x^3 + 56*x^2 + 4*x - 8;
print(nfdisc(P), "  ", factor(nfdisc(P)), "  subfields: ", apply(s -> s[1], select(s -> poldegree(s[1]) == 3, nfsubfields(P))));
nf = nfinit(P); Ob = ordinit(P, nf.zk); NI = Ninit(Ob); L = otpos(Ob, 11/3);
print(#L, " elements with tau < 11/3");
{for (i = 1, #L, my(al = oelt(Ob, L[i][2], L[i][3]), t = getabstime(), M);
  M = iferr(condN(NI, al), E, "ERR");
  print(i, " Tr=", trace(Mod(al, P)), " -> ", if (M === 0, "FAILS (N)", if (M === "ERR", "undecided", "holds")), "  ", getabstime() - t, "ms");
  if (M === 0, print("WITNESS ", al, "  minpoly ", minpoly(Mod(al, P))); break))}
