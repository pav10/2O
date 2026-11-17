\\ gp -q etype_zeta.gp etype_constants.gp < /dev/null : the constants of Theorem 1.1 and Table 1 of notes/etype-zeta.tex
\\ (assuming m(K) = r_d; for d in {2,3,4,5,7} this is automatic)
{
for (d = 2, 43,
  my(b = siegelb(2 * d), r = #b, A = (4 * Pi^2)^d * abs(b[r]), c = cnt(d, r), G = Gd(d), base = A * c * G^2);
  printf("d=%2d r_d=%d  A_d=%.4g  eps=%.3g  c_d=%.4g  G_d=%.4g  A_d*c_d*G_d^2=%.6g  log10(disc bound)=%.3f  regulator bound=%.4g\n",
    d, r, A, eps(d, r), c, G, base, (d - 1) * log(base) / log(10), gam(d)^d * A * c / 2^(d - 2)));
}
