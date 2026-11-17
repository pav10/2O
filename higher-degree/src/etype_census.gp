\\ census drivers for notes/etype-zeta.tex, Section 6.  Load after etype_zeta.gp.
\\ Field lists: one coefficient vector per line, ascending, leading coefficient omitted
\\ (the format of data/fields/degN.txt in pav10/variance, a mirror of the LMFDB lists).
\\ runP:  fields satisfying (1.2) (with m = 1); prints only the passing fields
\\ runNP: fields with non-principal different; prints [#T_1, #(nu failing (1.2))]
\\ runCert: exact certificates against (1.2) for all fields with |disc| <= dmax
runP(file, shard, nshard, out, maxn) =
{
  my(L = readvec(file));
  for (i = 1, min(#L, maxn), if (i % nshard != shard, next);
    my(P = Polrev(concat(L[i], [1])), r = condP(P));
    if (r[2] != "FAIL-principal" && r[2] != "FAIL-nonprincipal", write(out, i, " ", P, " ", nfdisc(P), " ", r)));
  write(out, "#done");
}
runNP(file, shard, nshard, out) =
{
  my(L = readvec(file));
  for (i = 1, #L, if (i % nshard != shard, next);
    my(P = Polrev(concat(L[i], [1])), bnf = bnfinit(P), cD);
    if (bnf.no == 1, next);
    cD = bnfisprincipal(bnf, bnf.nf.diff, 0);
    if (cD == 0 * cD, next);
    write(out, i, " ", P, " ", bnf.disc, " ", bnf.cyc, " ", cD~, " ", condPcount(bnf)));
  write(out, "#done");
}
runCert(file, dmax, out) =
{
  my(L = readvec(file), n = 0, bad = List());
  for (i = 1, #L, my(P = Polrev(concat(L[i], [1])), D = nfdisc(P));
    if (D > dmax, next); n++;
    my(c = certify(P)); write(out, D, " ", P, " ", c); if (c == 0, listput(bad, P)));
  write(out, "# fields with disc <= ", dmax, ": ", n, "; without certificate: ", Vec(bad));
}
