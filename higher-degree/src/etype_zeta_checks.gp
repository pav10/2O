\\ gp -q etype_zeta.gp etype_zeta_checks.gp < /dev/null
\\ (a) formula (2.1) for three cubic fields with class number 3 (two with non-principal different)
foreach([x^3 - 21*x - 28, x^3 - x^2 - 16*x + 8, x^3 - x^2 - 9*x + 8], P, csiegel(P));
\\ (b) signs of (-1)^d b_l(2d) for 2 <= d <= 43
{ for (d = 2, 43, my(b = siegelb(2 * d), pat = "");
    for (l = 1, #b, pat = concat(pat, if ((-1)^d * b[l] > 0, "+", "-")));
    print("d=", d, " r_d=", #b, " signs of (-1)^d b_l(2d): ", pat)); }
\\ (c) exact certificates that (1.2) fails: quadratic fields of discriminant 8, 12, 13
foreach([x^2 - 2, x^2 - 3, x^2 - x - 3], P, print(nfdisc(P), " ", P, " ", certify(P)));
