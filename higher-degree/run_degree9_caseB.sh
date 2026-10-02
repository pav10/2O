#!/bin/sh
# run_degree9_caseB.sh -- degree 9, Case B (the shortest vector lies in a cubic subfield).
# Case A of degree 9 is NOT done (see the notes); this script only handles Case B.
# Running time: about 4 hours on 4 cores (relative cubics), then about 1 hour.
set -e
cd "$(dirname "$0")/src"
D=../results/d9/caseB_z7
GP="gp -q lift6.gp order_filter.gp rel6.gp field6.gp sharp.gp pipeline.gp"

echo "== square forcing over the cubic subfields (degree-3 enumeration, [K:F] = 3 odd)"
printf 'subF("../results/sanity_d3.txt", 9)\n' | $GP sqforce8.gp

echo "== F = Q(zeta_7)^+: relative cubics (27 residues of b1 mod 3 O_F; resumable)"
mkdir -p $D
./run_c3_par.sh $D 4
cat $D/log_*.txt

echo "== maximal-order tests (absolute BC / sharpened rigidity)"
rm -f $D/fieldsA_26.txt
{ echo 'default(parisizemax, 4*10^9);'
  printf 'M = c3fields(['
  for j in $(seq 0 26); do printf '"%s/c3_%d.txt"' $D $j; [ $j -lt 26 ] && printf ', '; done
  echo ']);'
  echo "{for (i = 1, #M, write(\"$D/fieldsA_26.txt\", [M[i], nfdisc(M[i]), fieldtestA(M[i])]))}"
} > $D/ftA.gp
gp -q lift6.gp rel6.gp field6.gp sharp.gp rel8.gp rel9.gp $D/ftA.gp < /dev/null > /dev/null
grep -o '\["[A-Z0-9]*"' $D/fieldsA_26.txt | sort | uniq -c

echo "== condition (N): witness search (tau < 3) for the fields passing the absolute tests"
rm -f $D/witness_[0-2].txt
for p in 0 1 2; do
  printf 'default(parisizemax, 3*10^9)\nL = apply(e -> e[1], select(e -> e[3][1] == "OK", readvec("%s/fieldsA_26.txt")))\nwitnesslist(L, 3, "%s/witness_%d.txt", %d, 3)\n' $D $D $p $p \
    | gp -q lift6.gp field6.gp condN.gp witness.gp > /dev/null 2>&1 &
done
wait
cat $D/witness_[0-2].txt | sort > $D/witness_all.txt; rm -f $D/witness_[0-2].txt
echo "witnesses: $(grep -c WITNESS $D/witness_all.txt) of $(wc -l < $D/witness_all.txt)"
echo "== the one field without a witness below tau = 3 (witness at trace 27)"
gp -q lift6.gp field6.gp condN.gp witness.gp witness_last_d9.gp < /dev/null | tail -2
