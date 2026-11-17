\\ gp -q etype_zeta.gp etype_zeta_checks6.gp < /dev/null : formula (2.1) in degrees 4 and 6 (r_6 = 2, trace-2 terms enter)
foreach([x^4 - x^3 - 19*x^2 + 4*x + 76, x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1], P, csiegel(P));
\\ the first sextic field with class number > 1 in the LMFDB list (entry 3427)
csiegel(Polrev([-1, -16, 58, 17, -18, -1, 1]));
