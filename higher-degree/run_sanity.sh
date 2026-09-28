#!/bin/sh
# run_sanity.sh -- run the same machinery in degrees 2..5 (answer key: Kala-Yatsyna,
# arXiv:2402.03850, Thm 1.2).  Expected: Q(sqrt5), Q(zeta7)+ survive everything;
# Q(sqrt2), Q(zeta9)+ and three quartic fields (d_K = 725, 1600, 7168) fail (N).
set -e
cd "$(dirname "$0")/src"
R=../results
GP="gp -q lift6.gp order_filter.gp rel6.gp field6.gp condN.gp witness.gp sos2.gp pipeline.gp sanity.gp"
for D in 2 3 4 5; do gcc -O2 -DDEG=$D -o enum_d$D enum6.c -lm; ./run_enum.sh $D $R/sanity_d$D.txt 2>/dev/null; done
echo "== primitive case (x_1 generates K), degrees 2..5"
for D in 2 3 4 5; do echo "sanity(\"$R/sanity_d$D.txt\")" | $GP; done
echo "== (N) for the degree 2 and 3 survivors"
echo "print(findwitness(x^2-2, 3)); print(findwitness(x^3-3*x-1, 3)); \
      print(findwitness(x^2-x-1, 3)); print(findwitness(x^3-x^2-2*x+1, 3))" | $GP
echo "== degree 4, imprimitive case (F = Q(sqrt5), Q(sqrt2))"
rm -f $R/sanity_d4_imp.txt
echo "enum_e2g(\"Q(sqrt5)\", y^2-y-1, 5/4, 1.3110, \"$R/sanity_d4_imp.txt\"); \
      enum_e2g(\"Q(sqrt2)\", y^2-2, 2, 1.7091, \"$R/sanity_d4_imp.txt\"); \
      F = impfields([\"$R/sanity_d4_imp.txt\"]); print(F[1], \" survivors, fields: \", F[2]); \
      for (i = 1, #F[2], my(P = F[2][i][1], r = fieldtest(P)); print(P, \" \", r); \
        if (r[1] == \"OK\", print(\"   \", findwitness(P, 4))))" | $GP
