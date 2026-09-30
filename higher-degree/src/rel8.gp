\\ rel8.gp -- degree 8, Case B with F = Q(sqrt5) (the shortest vector x_1 lies in F, so
\\ nu(K) = 5/4), [K:F] = 4.  Requires lift6.gp, rel6.gp, field6.gp, sharp.gp.
\\ F is written in the variable y, K in x (as in rel6.gp).
\\
\\ Let beta in O_K \ F realise nu_{K/F}.  By Corollary 4.6 (e = 4, delta_F <= 1.3110)
\\ nu_{K/F} <= NUMAX8 = 4.2105...  Either F(beta) = K (relative quartic, enum_q4) or
\\ F(beta) = F' is a quartic field (relative quadratic, enum_q2), and then K/F' is quadratic.

FQ = y^2 - y - 1;
NUMAX8 = nurelmax(4, 1.3110);

\\ ---- exact relative downspread for given block minima/maxima ------------------
\\ min over gamma in O_F with rho_k(gamma) > -m_k of (1/2) sum_k (rho_k(gamma) + m_k)
\\ (the value delta_F(m) of Lemma 4.5 for this particular target), by a small search.
covdelta(mk, th) =
{
  my(best = oo);
  \\ gamma = u + v*phi ; rho_k(gamma) = u + v*th[k]
  for (v = floor(-mk[1] - mk[2]) - 8, ceil(-mk[1] - mk[2]) + 8,
    my(u0 = max(floor(-mk[1] - v*th[1]), floor(-mk[2] - v*th[2])) + 1);
    for (u = u0, u0 + 2,
      my(g1 = u + v*th[1] + mk[1], g2 = u + v*th[2] + mk[2]);
      if (g1 > 0 && g2 > 0, best = min(best, (g1 + g2)/2))));
  best;
}

\\ block data of a relative polynomial rp (in x, coefficients in Z[y]) at both embeddings
blockdata(rp, th) =
{
  vector(2, k, my(p = substpol(lift(Mod(rp, FQ)), y, th[k]), r = vecsort(real(polroots(p))));
    [r[1], r[#r], vecsum(r)/#r, r]);
}

\\ relative downspread test: both D_- = mean_k(mu_k - m_k) + delta(m) and
\\ D_+ = mean_k(M_k - mu_k) + delta(-M) must be >= vrel
reldown_ok(B, th, vrel) =
{
  my(Dm = sum(k = 1, 2, B[k][3] - B[k][1])/2 + covdelta([B[1][1], B[2][1]], th),
     Dp = sum(k = 1, 2, B[k][2] - B[k][3])/2 + covdelta([-B[1][2], -B[2][2]], th));
  Dm >= vrel - 1e-20 && Dp >= vrel - 1e-20;
}

\\ lattice points c = u + v*phi of O_F with rho_k(c) in [lo[k], hi[k]]
boxpts(lo, hi, th) =
{
  my(res = List(), M = [1, th[1]; 1, th[2]], Mi = M^-1, cs, bu, bv);
  if (lo[1] > hi[1] || lo[2] > hi[2], return([]));
  cs = [[lo[1], lo[2]], [lo[1], hi[2]], [hi[1], lo[2]], [hi[1], hi[2]]];
  bu = apply(c -> (Mi * c~)[1], cs); bv = apply(c -> (Mi * c~)[2], cs);
  for (u = floor(vecmin(bu)) - 1, ceil(vecmax(bu)) + 1,
    for (v = floor(vecmin(bv)) - 1, ceil(vecmax(bv)) + 1,
      my(r1 = u + v*th[1], r2 = u + v*th[2]);
      if (r1 >= lo[1] - 1e-20 && r1 <= hi[1] + 1e-20 && r2 >= lo[2] - 1e-20 && r2 <= hi[2] + 1e-20,
        listput(res, u + v*y))));
  Vec(res);
}

\\ Rolle interval for the constant term of g_m at embedding k (real coefficients e[0..m-1]),
\\ g_m(t) = sum_{i<=m} c_i C(4-i, m-i) t^(m-i), n = 4
rolle_interval(ce, m) =
{
  my(n = 4, gprev = sum(i = 0, m-1, ce[i+1] * binomial(n-i, m-1-i) * 't^(m-1-i)),
     h = sum(i = 0, m-1, ce[i+1] * binomial(n-i, m-i) * 't^(m-i)), r, lo = -oo, hi = oo);
  r = vecsort(real(polroots(gprev)));
  for (i = 1, #r, my(v = -subst(h, 't, r[i]), kk = m - i);
    if (kk % 2 == 0, lo = max(lo, v), hi = min(hi, v)));
  [lo - 1e-15, hi + 1e-15];
}

\\ ---- relative quartics ------------------------------------------------------
enum_q4(outfile, part = -1, nparts = 1) =
{
  my(nfF = nfinit(FQ), th = vecsort(real(polroots(FQ))), S = 8*NUMAX8, cnt = Map(), ncand = 0,
     emb = (a -> vector(2, k, subst(lift(Mod(a, FQ)), y, th[k]))));
  for (i = 0, 3, for (j = 0, 3,
    if (part >= 0 && (4*i + j) % nparts != part, next);
    my(c1 = i + j*y, e1 = emb(c1));
    \\ c2: s_k = rho_k(3/4 c1^2 - 2 c2) >= 0, s_1 + s_2 <= S
    my(C2 = boxpts(vector(2, k, 3*e1[k]^2/8 - S/2), vector(2, k, 3*e1[k]^2/8), th));
    foreach (C2, c2,
      my(e2 = emb(c2), s = vector(2, k, 3*e1[k]^2/4 - 2*e2[k]));
      if (s[1] + s[2] > S + 1e-20 || s[1] < -1e-20 || s[2] < -1e-20, next);
      my(I3 = vector(2, k, rolle_interval([1, e1[k], e2[k]], 3)),
         C3 = boxpts([I3[1][1], I3[2][1]], [I3[1][2], I3[2][2]], th));
      foreach (C3, c3,
        my(e3 = emb(c3), I4 = vector(2, k, rolle_interval([1, e1[k], e2[k], e3[k]], 4)),
           C4 = boxpts([I4[1][1], I4[2][1]], [I4[1][2], I4[2][2]], th));
        foreach (C4, c4,
          ncand++;
          my(rp = x^4 + c1*x^3 + c2*x^2 + c3*x + c4, vrel = trace(Mod(3/4*c1^2 - 2*c2, FQ))/8, B, res);
          if (vrel <= 0, next);
          B = blockdata(rp, th);
          if (!reldown_ok(B, th, vrel), mapput(cnt, "RD", if (mapisdefined(cnt, "RD"), mapget(cnt, "RD"), 0) + 1); next);
          my(fa = nffactor(nfF, rp));
          if (#fa[,1] != 1 || fa[1,2] != 1, next);           \\ reducible over F (or a power)
          my(RO = relorder(nfF, rp));
          if (RO === 0, next);
          res = reltest(RO, 5/4, vrel);
          mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
          write(outfile, [rp, RO[2], res]))))));
  print("relative quartics over Q(sqrt5): ", ncand, " candidates in the boxes; ", Mat(cnt));
}

\\ ---- relative quadratics: intermediate quartic fields F' = F(sqrt D) ----------
\\ beta = (-b1 + sqrt D)/2, Var_{K/F}(beta) = Tr_F(D)/8 <= NUMAX8
enum_q2() =
{
  my(nfF = nfinit(FQ), th = vecsort(real(polroots(FQ))), zk = nfF.zk, Tmax = 8*NUMAX8, T2, q, res = Map());
  T2 = matrix(2, 2, i, k, trace(Mod(zk[i]*zk[k], FQ)));
  q = qfminim(T2, floor(Tmax^2))[3];
  for (k = 1, #q, for (s = 0, 1,
    my(D = (1 - 2*s) * (q[1,k]*zk[1] + q[2,k]*zk[2]));
    if (!istotpos(FQ, D) || trace(Mod(D, FQ)) > Tmax, next);
    if (#nfroots(nfF, x^2 - D), next);
    for (r = 0, 3, my(b1 = (r % 2) + (r \ 2)*y, c0 = lift(Mod(b1^2 - D, FQ)/4));
      if (!nfeltisintegral(nfF, c0), next);
      my(rp = x^2 + b1*x + c0, vrel = trace(Mod(D, FQ))/8, B = blockdata(rp, th));
      if (!reldown_ok(B, th, vrel), next);
      my(P = polredabs(rnfequation(nfF, rp)));
      mapput(res, P, 1))));
  Mat(res)[,1];
}

\\ ---- relative quadratics over a general base F (used for the intermediate quartic F') ----
\\ Exact relative downspread prefilter.  beta has conjugates mu_k +- s_k over the k-th real place
\\ of F.  If beta realises nu_{K/F} then, by relative budget-cost, tau(beta + gamma) >= nu_{K/F}
\\ for every gamma in O_F with beta + gamma >> 0 (and likewise for gamma - beta).
\\ shift_exists(E, Ei, mu, s, vrel): is there gamma = E u (u integral) with
\\   z_k = rho_k(gamma) + mu_k - s_k > 0 for all k and  sum_k (mu_k + rho_k(gamma)) < f * vrel ?
\\ z lies in the simplex {z > 0, sum z < Bp}; u = Ei (z - m) ranges over the image of its vertices.
shift_exists(E, Ei, mu, s, vrel) =
{
  my(f = #mu, m = vector(f, k, mu[k] - s[k]), Bp = f*vrel - vecsum(mu) + vecsum(m), V, lo, hi, cs);
  if (Bp <= 0, return(0));
  V = concat([vector(f)~], vector(f, k, Bp * matid(f)[,k]));
  V = apply(v -> Ei * (v - m~), V);
  lo = vector(f, i, floor(vecmin(vector(#V, j, V[j][i])) - 1e-9));
  hi = vector(f, i, ceil(vecmax(vector(#V, j, V[j][i])) + 1e-9));
  cs = sum(k = 1, f, E[k, f]);
  forvec (w = vector(f - 1, i, [lo[i], hi[i]]),
    my(a = vector(f, k, m[k] + sum(i = 1, f - 1, E[k, i]*w[i])), l = -oo, h = oo, A = vecsum(a));
    \\ z_k = a_k + E[k,f] u_f > 0
    for (k = 1, f, my(c = E[k, f]); if (c > 0, l = max(l, -a[k]/c), if (c < 0, h = min(h, -a[k]/c), if (a[k] <= 0, next(2)))));
    \\ sum z = A + cs u_f < Bp   (strict, with margin: only certain shifts count)
    if (cs > 0, h = min(h, (Bp - A - 1e-20)/cs), if (cs < 0, l = max(l, (Bp - A - 1e-20)/cs), if (A >= Bp - 1e-20, next)));
    my(ul = floor(l) + 1, uh = ceil(h) - 1);
    if (ul <= uh, return(1)));
  0;
}

\\ as enum_e2g (rel6.gp) with the exact relative downspread prefilter before relorder/reltest
enum_e2g8(nm, Fpol, nuF, dl, outfile) =
{
  my(nfF = nfinit(Fpol), f = poldegree(Fpol), d = 2*f, numax = nurelmax(2, dl), Tmax = 2*d*numax,
     zk = nfF.zk, pF = nfF.pol, T2, q, reps, cnt = Map(), Ds = List(), th = real(polroots(pF)), E, Ei,
     emb = (a -> vector(f, k, subst(lift(Mod(a, pF)), variable(pF), th[k]))));
  E = matrix(f, f, k, i, emb(zk[i])[k]); Ei = E^-1;
  reps = vector(2^f, t, my(v = binary(t - 1 + 2^f)[2..f+1]); sum(i = 1, f, v[i]*zk[i]));
  T2 = matrix(f, f, i, k, trace(Mod(zk[i]*zk[k], pF)));
  q = qfminim(T2, floor(Tmax^2))[3];
  for (k = 1, #q, for (s = 0, 1, my(D = (1 - 2*s) * sum(i = 1, f, q[i,k]*zk[i]));
      if (istotpos(pF, D) && trace(Mod(D, pF)) <= Tmax, listput(Ds, D))));
  print(nm, " (d=", d, "): nu_rel <= ", numax, ", Tr(Delta) <= ", Tmax, ", ", #Ds, " totally positive Delta");
  for (i = 1, #Ds,
    my(D = Ds[i], sD);
    if (#nfroots(nfF, x^2 - D), next);
    sD = apply(t -> sqrt(t)/2, emb(D));
    for (r = 1, #reps,
      my(b1 = reps[r], c0 = lift(Mod(b1^2 - D, pF) / 4), mu, vrel);
      if (!nfeltisintegral(nfF, c0), next);
      vrel = trace(Mod(D, pF)) / (2*d);
      mu = apply(t -> -t/2, emb(b1));
      if (shift_exists(E, Ei, mu, sD, vrel) || shift_exists(E, Ei, -mu, sD, vrel),
        mapput(cnt, "RD", if (mapisdefined(cnt, "RD"), mapget(cnt, "RD"), 0) + 1); next);
      my(rp = x^2 + b1*x + c0, RO = relorder(nfF, rp), res);
      if (RO === 0, next);
      res = reltest(RO, nuF, vrel);
      mapput(cnt, res[1], if (mapisdefined(cnt, res[1]), mapget(cnt, res[1]), 0) + 1);
      write(outfile, [nm, D, b1, RO[2], res])));
  print(nm, ": ", Mat(cnt));
}
