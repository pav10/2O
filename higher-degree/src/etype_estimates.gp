\\ gp -q etype_zeta.gp etype_estimates.gp < /dev/null : Lemmas 4.1-4.3 of notes/etype-zeta.tex on sample fields
est(P) =
{
  my(bnf = bnfinit(P, 1), d = poldegree(P), D = abs(bnf.disc), T1 = #traceL(bnf.nf, 1), Z = Mat(partialzeta2(bnf)),
     zmin = vecmin(Z[, 2]), lb = 2^(d - 2) * bnf.reg / (gam(d)^d * D));
  printf("%s  disc=%d h=%d | N(xi*D)=%d vs G_d*D^((d-2)/(2d-2))=%.3f | #T1=%d vs c_d1*sqrt(D)=%.1f, sqrt(D)/(d-1)!=%.1f | min_C zeta(2;C)=%.4g vs 2^(d-2)R/(gamma_d^d D)=%.4g (ratio %.1f)\n",
    P, D, bnf.no, smallideal(P), Gd(d) * D^((d - 2) / (2 * d - 2)), T1, cnt(d, 1) * sqrt(D), sqrt(D) / (d - 1)!, zmin, lb, zmin / lb);
}
foreach([x^3-x^2-2*x+1, x^3-21*x-28, x^3-x^2-16*x+8, x^3-435*x-2774, x^3-210*x-497, x^4-x^3-3*x^2+x+1, x^4-6*x^2+4, x^4-x^3-19*x^2+4*x+76, x^5-x^4-4*x^3+3*x^2+3*x-1, x^2-79], P, est(P));
