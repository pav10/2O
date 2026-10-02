\\ rel9.gp -- degree 9, Case B: the shortest vector x_1 lies in a cubic subfield F (the only proper
\\ subfields of a nonic field are cubic).  Square forcing (sqforce8.gp, subF(..., 9)) leaves only
\\ F = Q(zeta_7)^+ (nu = 14/9).  Then [K:F] = 3, so K = F(beta) for every beta in O_K \ F, and we
\\ enumerate the relative cubic polynomials of a beta realising nu_{K/F}.
\\ Requires lift6.gp rel6.gp field6.gp sharp.gp rel8.gp (nurelmax, relorder, reltest, shift_exists).
\\ F is written in the variable y, K in x.

\\ all u in Z^f with E*u in the box prod [lo_k, hi_k]  (E: embeddings of the integral basis)
latbox(lo, hi, E, Ei) =
{
  my(f = #lo, cs = List(), res = List(), rng);
  forvec (t = vector(f, k, [0, 1]), listput(cs, vector(f, k, if (t[k], hi[k], lo[k]))~));
  rng = vector(f, i, my(v = apply(c -> (Ei * c)[i], Vec(cs))); [floor(vecmin(v)) - 1, ceil(vecmax(v)) + 1]);
  forvec (u = rng,
    my(r = E * u~, ok = 1);
    for (k = 1, f, if (r[k] < lo[k] - 1e-20 || r[k] > hi[k] + 1e-20, ok = 0; break));
    if (ok, listput(res, u~)));
  Vec(res);
}

\\ relative cubics t^3 + b1 t^2 + b2 t + b3 over O_F (F cubic, e = 3, d = 9).
\\ beta realises nu_{K/F}; Var_{K/F}(beta) = (1/9) sum_k s_k with s_k = rho_k(2/3 b1^2 - 2 b2) >= 0.
\\ b1 runs over O_F / 3 O_F (translation beta -> beta + gamma), b2 over the box s_k >= 0,
\\ sum_k s_k <= 9 nu_max, b3 over the real-rootedness interval at each real place (critical values).
\\ Prefilter: exact relative downspread (shift_exists of rel8.gp) for beta and -beta.
enum_c3(nm, Fpol, nuF, dl, outfile, part = 0, nparts = 1) =
{
  my(nfF = nfinit(Fpol), f = 3, d = 9, numax = nurelmax(3, dl), S = d*numax, zk = nfF.zk, pF = nfF.pol,
     th = vecsort(real(polroots(pF))), E, Ei, emb, cnt = Map(), ncand = 0, B1 = List(), j = 0);
  emb = (a -> vector(f, k, subst(lift(Mod(a, pF)), variable(pF), th[k])));
  E = matrix(f, f, k, i, emb(zk[i])[k]); Ei = E^-1;
  forvec (t = vector(f, k, [0, 2]), listput(B1, sum(i = 1, f, t[i]*zk[i])));
  print(nm, ": nu_rel <= ", numax, ", sum s_k <= ", S);
  foreach (B1, b1,
    j++; if ((j - 1) % nparts != part, next);
    my(e1 = emb(b1), lo2 = vector(f, k, e1[k]^2/3 - S/2), hi2 = vector(f, k, e1[k]^2/3));
    foreach (latbox(lo2, hi2, E, Ei), u2,
      my(b2 = sum(i = 1, f, u2[i]*zk[i]), e2 = emb(b2), s = vector(f, k, 2/3*e1[k]^2 - 2*e2[k]),
         lo3 = vector(f), hi3 = vector(f), ok = 1);
      if (vecmin(s) < -1e-20 || vecsum(s) > S + 1e-20, next);
      for (k = 1, f,
        my(dq = 3*'t^2 + 2*e1[k]*'t + e2[k], g = 't^3 + e1[k]*'t^2 + e2[k]*'t, rr);
        if (poldisc(dq) < 0, ok = 0; break);
        rr = vecsort(real(polroots(dq)));
        lo3[k] = -subst(g, 't, rr[1]) - 1e-15; hi3[k] = -subst(g, 't, rr[2]) + 1e-15);
      if (!ok, next);
      foreach (latbox(lo3, hi3, E, Ei), u3,
        my(b3 = sum(i = 1, f, u3[i]*zk[i]), rp = x^3 + b1*x^2 + b2*x + b3, vrel, mu, rts, sm, sp, RO, res);
        ncand++;
        vrel = trace(Mod(2/3*b1^2 - 2*b2, pF)) / d;
        if (vrel <= 0, next);
        \\ conjugates of beta over the k-th place of F: mean mu_k, spreads to the minimum / maximum
        rts = vector(f, k, vecsort(real(polroots(x^3 + e1[k]*x^2 + e2[k]*x + emb(b3)[k]))));
        mu = vector(f, k, -e1[k]/3);
        sm = vector(f, k, mu[k] - rts[k][1]); sp = vector(f, k, rts[k][3] - mu[k]);
        if (shift_exists(E, Ei, mu, sm, vrel) || shift_exists(E, Ei, -mu, sp, vrel),
          mapput(cnt, "RD", if (mapisdefined(cnt, "RD"), mapget(cnt, "RD"), 0) + 1); next);
        if (#nfroots(nfF, rp), next);                    \\ reducible over F
        RO = relorder(nfF, rp);
        if (RO === 0, next);
        res = reltest(RO, nuF, vrel);
        mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
        write(outfile, [nm, rp, RO[2], res]))));
  print(nm, " part ", part, "/", nparts, ": ", ncand, " relative cubics in the boxes; ", Mat(cnt));
}

\\ absolute part of fieldtest2 (sharp.gp): BC and sharpened rigidity in the maximal order.
\\ The relative part of fieldtest2 enumerates elements of trace < 3/2 nu_{K/F}, which is too
\\ large in degree 9; it is applied separately to the fields surviving this test.
fieldtestA(P) =
{
  my(nf = nfinit(P), Ob = ordinit(P, nf.zk), nu = onu(Ob), L);
  L = otpos(Ob, nu);
  if (#L, return(["BC", nu, L[1][1], oelt(Ob, L[1][2], L[1][3])]));
  L = otpos(Ob, 3/2*nu);
  for (i = 1, #L,
    if (oshape2(Ob, L[i][2], L[i][3], nu) === 0,
      return(["RG2", nu, L[i][1], oelt(Ob, L[i][2], L[i][3])])));
  ["OK", nu];
}

\\ distinct fields (polredabs) of the order-level survivors in the given files
c3fields(files) =
{
  my(F = Map());
  for (i = 1, #files, my(L = readvec(files[i]));
    for (j = 1, #L, my(e = L[j], res = e[#e]); if (res[1] == "OK", mapput(F, polredabs(e[#e - 1]), 1))));
  if (#F, Mat(F)[,1], []);
}
