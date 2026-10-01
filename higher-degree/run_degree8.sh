#!/bin/sh
# run_degree8.sh -- reproduce the degree-8 computation end to end.
# Requirements: gcc, PARI/GP >= 2.15, python3 (+numpy for deltaF.py).
# Case A (x_1 generates K): rung-split enumeration, < 1 minute on 4 cores.
# Case B (x_1 in a proper subfield): square forcing; for F = Q(sqrt5) relative quartics
# (about 1.5 hours on 4 cores) and relative quadratics over the two intermediate quartic fields.
set -e
cd "$(dirname "$0")/src"
R=../results/d8
GP="gp -q lift6.gp order_filter.gp rel6.gp field6.gp sharp.gp pipeline.gp"
mkdir -p $R

echo "== Case A: build and split enumeration (modes X / YL / YU)"
gcc -O2 -DDEG=8 -o enum_d8 enum6.c -lm
rm -rf $R/caseA
./run_split_par.sh 8 $R/caseA sharp 4
cat $R/caseA/out_*.txt | sort -u > $R/d8_split_raw.txt
echo "distinct polynomials: $(wc -l < $R/d8_split_raw.txt)"

echo "== Case A: exact order-level tests (NS / BC / RG), then maximal order"
for p in 0 1 2 3; do
  echo "primitive(\"$R/d8_split_raw.txt\", \"$R/prim8_status_$p.tmp\", $p, 4)" | $GP > /dev/null &
done
wait
cat $R/prim8_status_*.tmp > $R/prim8_status.txt; rm -f $R/prim8_status_*.tmp
echo "irreducible: $(wc -l < $R/prim8_status.txt)"
for s in NS BC RG OK; do echo "  $s: $(grep -c "\"$s\"" $R/prim8_status.txt)"; done
echo "L = readvec(\"$R/prim8_status.txt\"); for (i = 1, #L, if (L[i][2][1] == \"OK\", \
  my(P = polredabs(L[i][1])); print(P, \"  d=\", nfdisc(P), \"  \", fieldtest2(P)[1])))" | $GP

echo "== Case B: square forcing (F = Q(x_1) of degree 4 and 2, from the degree-4 and degree-2 enumerations)"
[ -s ../results/sanity_d4.txt ] || { gcc -O2 -DDEG=4 -o enum_d4 enum6.c -lm; ./run_enum.sh 4 ../results/sanity_d4.txt 2>/dev/null; }
printf 'subF("../results/sanity_d4.txt")\nsubF(\"../results/sanity_d2.txt\")\n' \
  | $GP sqforce8.gp

echo "== Case B, F = Q(sqrt5): intermediate quartic fields F' (relative quadratics beta)"
echo 'L = enum_q2(); for (i = 1, #L, print(L[i], "  d=", nfdisc(L[i]), "  ", forceup(L[i], 5/4, 8)))' \
  | $GP sqforce8.gp rel8.gp
echo "== Case B, F = Q(sqrt5): K/F' quadratic (delta_F' bounds from deltaF.py: F725 grid 0.05, F1600 grid 0.1)"
mkdir -p $R/caseB_sqrt5
rm -f $R/caseB_sqrt5/q2_rel_F725.txt $R/caseB_sqrt5/q2_rel_F1600.txt
echo "enum_e2g8(\"F725\", y^4-y^3-3*y^2+y+1, 5/4, 1.9144, \"$R/caseB_sqrt5/q2_rel_F725.txt\")" | $GP rel8.gp &
echo "enum_e2g8(\"F1600\", y^4-6*y^2+4, 5/4, 2.1567, \"$R/caseB_sqrt5/q2_rel_F1600.txt\")" | $GP rel8.gp &
wait
echo "== Case B, F = Q(sqrt5): relative quartics (10 residues of c1 up to conjugation)"
rm -rf $R/caseB_sqrt5/q4
./run_q4_par.sh $R/caseB_sqrt5/q4 4
cat $R/caseB_sqrt5/q4/log_*.txt

echo "== Case B, F = Q(sqrt5): maximal-order tests on the order-level survivors"
FILES="[\"$R/caseB_sqrt5/q2_rel_F725.txt\", \"$R/caseB_sqrt5/q2_rel_F1600.txt\""
for j in 0 1 2 3 4 5 6 8 9 12; do FILES="$FILES, \"$R/caseB_sqrt5/q4/q4_$j.txt\""; done
FILES="$FILES]"
for p in 0 1 2 3; do rm -f $R/caseB_sqrt5/fields_$p.txt
  echo "collect8($FILES, \"$R/caseB_sqrt5/fields_$p.txt\", $p, 4)" | $GP collect8.gp > /dev/null &
done
wait
cat $R/caseB_sqrt5/fields_[0-3].txt | sort -u -t, -k1,1 > $R/caseB_sqrt5/fields_all.txt
rm -f $R/caseB_sqrt5/fields_[0-3].txt
grep -o '\["[A-Z0-9]*"' $R/caseB_sqrt5/fields_all.txt | sort | uniq -c
grep '"OK"' $R/caseB_sqrt5/fields_all.txt | cut -c1-100

echo "== condition (N) for the two surviving fields"
echo 'print(findwitness(x^8 - 3*x^7 - 4*x^6 + 13*x^5 + 5*x^4 - 13*x^3 - 4*x^2 + 3*x + 1, 5/2))' \
  | gp -q lift6.gp field6.gp condN.gp witness.gp
gp -q lift6.gp field6.gp condN.gp witness.gp witness_k1_d8.gp < /dev/null
