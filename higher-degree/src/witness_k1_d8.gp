default(parisizemax, 2*10^9);
P = x^8 - 2*x^7 - 12*x^6 + 26*x^5 + 17*x^4 - 36*x^3 - 5*x^2 + 11*x - 1;
nf = nfinit(P); Ob = ordinit(P, nf.zk); NI = Ninit(Ob); L = otpos(Ob, 5/2);
print(#L, " elements with tau < 5/2");
{for (i = 1, #L, my(al = oelt(Ob, L[i][2], L[i][3]), t = getabstime(), M);
  M = iferr(condN(NI, al), E, "ERR");
  print(i, " Tr=", trace(Mod(al, P)), " ", minpoly(Mod(al, P)), " -> ", if (M === 0, "FAILS (N)", if (M === "ERR", "undecided (stack)", "holds")), "  ", getabstime() - t, "ms");
  if (M === 0, print("WITNESS ", al); break))}
