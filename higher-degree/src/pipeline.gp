\\ pipeline.gp -- drivers for the degree-6 proof (requires lift6.gp order_filter.gp
\\ rel6.gp field6.gp condN.gp witness.gp sos2.gp).  See ../README.md.

\\ read the enumerator output (lines "a1 a2 ... ad") into polynomials
readenum(file) =
{
  my(L = externstr(concat("cat ", file)));
  apply(s -> Pol(concat([1], eval(concat(["[", strjoin(strsplit(s, " "), ","), "]"])))), L);
}

\\ ---- step A: primitive case, order-level filters -------------------------
\\ writes one line [f, status] per irreducible polynomial; returns the counts.
\\ (part, nparts) splits the work: polynomial number i is handled iff i = part mod nparts.
primitive(infile, outfile, part = 0, nparts = 1) =
{
  my(P = readenum(infile), irr, cnt = Map(), surv = List());
  irr = select(p -> polisirreducible(p), P);
  write1(outfile, "");
  for (i = 1, #irr,
    if ((i - 1) % nparts != part, next);
    my(r = filt(irr[i]));
    write(outfile, [irr[i], r]);
    mapput(cnt, r[1], if (mapisdefined(cnt, r[1]), mapget(cnt, r[1]), 0) + 1);
    if (r[1] == "OK", listput(surv, irr[i])));
  [#P, #irr, Mat(cnt), Vec(surv)];
}

\\ ---- step B: which subfields F = Q(x_1) are possible ----------------------
\\ Only tests valid for x_1 in a PROPER subfield of K are used: NS and BC
\\ (elements of Z[x_1] lie in O_K with the same tau and Var).  Rigidity is NOT used,
\\ since the y in y^2 + c y + b may lie outside F.
filtNoRG(f) =
{
  my(Ob = ordZx(f), a = Vec(f), d = poldegree(f), V, nu, L);
  V = ((d-1)*a[2]^2/d - 2*a[3]) / d;
  nu = onu(Ob);
  if (nu < V, return(["NS", nu]));
  L = otpos(Ob, V);
  if (#L, return(["BC", L[1][1], oelt(Ob, L[1][2], L[1][3])]));
  ["OK", V];
}
subfieldcands(infile) =
{
  my(P = select(p -> polisirreducible(p), readenum(infile)), res = List());
  for (i = 1, #P, my(r = filtNoRG(P[i])); if (r[1] == "OK", listput(res, [polredabs(P[i]), r[2]])));
  \\ distinct fields with their nu
  my(M = Map()); for (i = 1, #res, mapput(M, res[i][1], res[i][2]));
  Mat(M);
}

\\ ---- step C: collect and deduplicate the imprimitive survivors -------------
impfields(files) =
{
  my(L = concat(apply(f -> readvec(f), files)), S, M = Map());
  S = select(v -> v[#v][1] == "OK", L);
  for (i = 1, #S,
    my(R = polredabs(S[i][#S[i]-1]), l);
    if (!mapisdefined(M, R), mapput(M, R, List()));
    l = mapget(M, R); listput(l, [S[i][1], S[i][#S[i]][2]]); mapput(M, R, l));
  my(K = Mat(M), out);
  out = vector(#K~, i, [K[i,1], nfdisc(K[i,1]), Vec(K[i,2])]);
  [#S, vecsort(out, 2)];
}

\\ ---- step D/E: field tests, then (N)-witnesses -----------------------------
finalfields(fields, outfile) =
{
  my(res = List());
  write1(outfile, "");
  for (i = 1, #fields,
    my(P = fields[i][1], r = fieldtest(P), w = 0, s = 0);
    if (r[1] == "OK",
      w = findwitness(P, 4);
      if (w[2] == "WITNESS",
        \\ independent cross-check via Mordell-Ko + even rank-6 determinant bound
        s = [sos2(P, w[3])[1], trace(Mod(w[3], P)), e6bound(P)]));
    write(outfile, [P, nfdisc(P), r, w, s]);
    listput(res, [P, r[1], if (w === 0, "-", w[2])]));
  Vec(res);
}
