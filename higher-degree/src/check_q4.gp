\\ check_q4.gp -- validation of the relative-quartic candidate generation of enum_q4 (rel8.gp).
\\ Requires lift6.gp rel6.gp field6.gp sharp.gp rel8.gp.
\\
\\ q4cands(c1): the (c2, c3, c4) that enum_q4 visits for this c1 (same code as enum_q4).
\\ q4brute(c1, smax): all (c2, c3, c4) with 0 <= s_k, s_1 + s_2 <= smax, such that
\\   x^4 + c1 x^3 + c2 x^2 + c3 x + c4 is real-rooted at both real places of F = Q(sqrt5),
\\   found by an independent search: coefficient boxes from power-sum bounds and an exact
\\   real-rootedness test (all roots of the norm polynomial real, by Sturm).
\\ Every brute-force polynomial must be among the enum_q4 candidates.

q4cands(c1) =
{
  my(th = vecsort(real(polroots(FQ))), S = 8*NUMAX8, res = List(),
     emb = (a -> vector(2, k, subst(lift(Mod(a, FQ)), y, th[k]))), e1 = emb(c1));
  my(C2 = boxpts(vector(2, k, 3*e1[k]^2/8 - S/2), vector(2, k, 3*e1[k]^2/8), th));
  foreach (C2, c2,
    my(e2 = emb(c2), s = vector(2, k, 3*e1[k]^2/4 - 2*e2[k]));
    if (s[1] + s[2] > S + 1e-20 || s[1] < -1e-20 || s[2] < -1e-20, next);
    my(I3 = vector(2, k, rolle_interval([1, e1[k], e2[k]], 3)),
       C3 = boxpts([I3[1][1], I3[2][1]], [I3[1][2], I3[2][2]], th));
    foreach (C3, c3,
      my(e3 = emb(c3), I4 = vector(2, k, rolle_interval([1, e1[k], e2[k], e3[k]], 4)),
         C4 = boxpts([I4[1][1], I4[2][1]], [I4[1][2], I4[2][2]], th));
      foreach (C4, c4, listput(res, [c2, c3, c4]))));
  Vec(res);
}

\\ all u + v*y with rho_k in [lo_k, hi_k], by a plain double loop (independent of boxpts)
latpts(lo, hi, th) =
{
  my(res = List(), sq5 = th[2] - th[1], vlo = floor((lo[2] - hi[1]) / sq5) - 1, vhi = ceil((hi[2] - lo[1]) / sq5) + 1);
  for (v = vlo, vhi,
    for (u = floor(lo[1] - v*th[1]) - 1, ceil(hi[1] - v*th[1]) + 1,
      my(r1 = u + v*th[1], r2 = u + v*th[2]);
      if (r1 >= lo[1] && r1 <= hi[1] && r2 >= lo[2] && r2 <= hi[2], listput(res, u + v*y))));
  Vec(res);
}

realrooted2(rp) =
{
  my(Q = polresultant(rp, FQ, y), g);
  g = Q / gcd(Q, Q');
  polsturm(g) == poldegree(g);
}

q4brute(c1, smax) =
{
  my(th = vecsort(real(polroots(FQ))), res = List(), emb = (a -> vector(2, k, subst(lift(Mod(a, FQ)), y, th[k]))), e1 = emb(c1));
  \\ c2: s_k = rho_k(3/4 c1^2 - 2 c2) in [0, smax]
  my(C2 = latpts(vector(2, k, (3/4*e1[k]^2 - smax)/2 - 1e-9), vector(2, k, 3/8*e1[k]^2 + 1e-9), th));
  foreach (C2, c2,
    my(e2 = emb(c2), s = vector(2, k, 3/4*e1[k]^2 - 2*e2[k]));
    if (s[1] < -1e-20 || s[2] < -1e-20 || s[1] + s[2] > smax + 1e-20, next);
    \\ centred deviations d (sum 0, sum d^2 = s): |p3(d)| <= s^{3/2}, s^2/4 <= p4(d) <= s^2,
    \\ e2(d) = -s/2, e3(d) = p3/3, e4(d) = (s^2/2 - p4)/4.  With mu = -c1/4:
    \\ c3 = -(e3(d) + 2 mu e2(d) + 4 mu^3),  c4 = e4(d) + mu e3(d) + mu^2 e2(d) + mu^4.
    my(lo3 = vector(2), hi3 = vector(2), lo4 = vector(2), hi4 = vector(2));
    for (k = 1, 2,
      my(mu = -e1[k]/4, sk = max(s[k], 0), B3 = sk^(3/2)/3, e2d = -sk/2, base3 = 2*mu*e2d + 4*mu^3,
         base4 = mu^2*e2d + mu^4);
      lo3[k] = -(B3 + base3) - 1e-9; hi3[k] = -(-B3 + base3) + 1e-9;
      lo4[k] = base4 - sk^2/8 - abs(mu)*B3 - 1e-9; hi4[k] = base4 + sk^2/16 + abs(mu)*B3 + 1e-9);
    my(C3 = latpts(lo3, hi3, th), C4 = latpts(lo4, hi4, th));
    foreach (C3, c3, foreach (C4, c4,
      my(rp = x^4 + c1*x^3 + c2*x^2 + c3*x + c4);
      if (realrooted2(rp), listput(res, [c2, c3, c4])))));
  Vec(res);
}

\\ compare for c1 = i + j*y: brute force (s_1 + s_2 <= smax) must be contained in q4cands
cmp_q4(i, j, smax) =
{
  my(c1 = i + j*y, A = q4cands(c1), B = q4brute(c1, smax), SA = Set(apply(t -> apply(z -> lift(Mod(z, FQ)), t), A)), miss = 0);
  foreach (B, t, if (!setsearch(SA, apply(z -> lift(Mod(z, FQ)), t)), miss++; print("MISSING: c1=", c1, " ", t)));
  [c1, "enum_q4 candidates", #A, "real-rooted (s1+s2 <= ", smax, ")", #B, "missing", miss];
}
