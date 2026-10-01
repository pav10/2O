#!/bin/sh
# validate8.sh -- independent checks of the degree-8 computation (run after run_degree8.sh).
#  (1) every exclusion made inside the split enumerator (NS, sharpened RIG, y-budget) on an
#      irreducible polynomial is re-checked in exact arithmetic (verify_split.gp, vsharp);
#  (2) condition (N) is decided directly (condN.gp) at the element named by every
#      maximal-order exclusion (Case A: RG2; Case B: RGR2, BCR), which must fail (N);
#  (3) the relative-quartic candidate generation of enum_q4 is compared with an independent
#      brute-force search (check_q4.gp), for all 16 residues of c1;
#  (4) the rung split is compared with the plain enumerator on all 115 subtrees of Case A
#      (compare_split.sh), every difference certified by verify_split.gp.
# Running time: (1) 5 min, (2) 1 min, (3) about 1 hour, (4) several hours, on 4 cores.
set -e
cd "$(dirname "$0")/src"
V=../results/d8/validate; mkdir -p $V
gcc -O2 -DDEG=8 -o enum_d8 enum6.c -lm

echo "== (1) exact re-check of the in-enumerator exclusions (Case A)"
./run_split_par.sh 8 $V/why whysharp 4
cat $V/why/out_*.txt | grep -v '^OUT' > $V/why.txt
for p in 0 1 2; do
  echo "default(realprecision, 80); print(vsharp(\"$V/why.txt\", $p, 3))" \
    | gp -q lift6.gp field6.gp sharp.gp verify_split.gp > $V/vsharp_$p.txt 2>&1 &
done
wait
cat $V/vsharp_*.txt
rm -rf $V/why $V/why.txt

echo "== (2) condition (N) at every certificate element"
rm -f $V/caseA_fields.txt $V/certN_*.txt
echo "foreach([x^8-8*x^6+16*x^4-8*x^2+1, x^8-8*x^6+18*x^4-12*x^2+2, x^8-8*x^6+17*x^4-10*x^2+1, \
  x^8-8*x^6+19*x^4-14*x^2+1], f, my(P = polredabs(f)); write(\"$V/caseA_fields.txt\", [P, nfdisc(P), fieldtest2(P)]))" \
  | gp -q lift6.gp rel6.gp field6.gp sharp.gp
for p in 0 1 2; do
  printf "default(parisizemax, 2*10^9)\nprint(certN(\"../results/d8/caseB_sqrt5/fields_all.txt\", \"$V/certN_B_$p.txt\", $p, 3))\n" \
    | gp -q lift6.gp field6.gp condN.gp certN8.gp 2>/dev/null &
done
wait
echo "print(certN(\"$V/caseA_fields.txt\", \"$V/certN_A.txt\"))" | gp -q lift6.gp field6.gp condN.gp certN8.gp

echo "== (3) relative quartics: enum_q4 candidates versus brute force"
echo 'for (i = 0, 3, for (j = 0, 3, print(cmp_q4(i, j, 8*NUMAX8))))' \
  | gp -q lift6.gp rel6.gp field6.gp sharp.gp rel8.gp check_q4.gp 2>/dev/null > $V/check_q4.txt
cat $V/check_q4.txt

echo "== (4) rung split versus plain enumeration"
./compare_split.sh 8 0 24 $V/split_0_24 4 3600
./compare_split.sh 8 24 30 $V/split_24_30 4 3600
./compare_split.sh 8 30 36 $V/split_30_36 3 5400
./compare_split.sh 8 36 47 $V/split_36_47 4 86400
