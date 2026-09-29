#!/bin/sh
# run_degree6.sh -- reproduce the degree-6 computation end to end.
# Requirements: gcc, PARI/GP >= 2.15, python3 + numpy (only for deltaF.py).
# Total running time on 4 cores: roughly 10 minutes.
set -e
cd "$(dirname "$0")/src"
R=../results
GP="gp -q lift6.gp order_filter.gp rel6.gp field6.gp condN.gp witness.gp sos2.gp pipeline.gp"

echo "== build enumerators"
for D in 2 3 4 5 6; do gcc -O2 -DDEG=$D -o enum_d$D enum6.c -lm; done

echo "== step A: primitive case (x_1 generates K)"
./run_enum.sh 6 $R/prim_raw.txt 2> $R/prim_enum.log
for p in 0 1 2 3; do
  echo "primitive(\"$R/prim_raw.txt\", \"$R/prim_status_$p.tmp\", $p, 4)" | $GP &
done
wait
cat $R/prim_status_*.tmp > $R/prim_status.txt; rm -f $R/prim_status_*.tmp
echo "irreducible: $(wc -l < $R/prim_status.txt)"
for s in NS BC RG OK; do echo "  $s: $(grep -c "\"$s\"" $R/prim_status.txt)"; done

echo "== step B: possible subfields F = Q(x_1) (degrees 2 and 3; NS/BC only)"
./run_enum.sh 2 $R/sub_d2_raw.txt 2>/dev/null
./run_enum.sh 3 $R/sub_d3_raw.txt 2>/dev/null
echo "print(subfieldcands(\"$R/sub_d2_raw.txt\")); print(subfieldcands(\"$R/sub_d3_raw.txt\"))" | $GP

echo "== step C: positive covering constants delta_F"
python3 deltaF.py

echo "== step D: relative enumeration K = F(y)"
rm -f $R/imp_F*.txt
for spec in "e3 1" "e3 2" "e2 3" "e2 4"; do
  set -- $spec
  echo "enum_$1($2, \"$R/imp_F$2.txt\")" | $GP > $R/imp_F$2.log 2>&1 &
done
wait
cat $R/imp_F*.log

echo "== step E: maximal-order tests and (N)-witnesses"
rm -f $R/final_fields.txt
echo "F = impfields([\"$R/imp_F1.txt\",\"$R/imp_F2.txt\",\"$R/imp_F3.txt\",\"$R/imp_F4.txt\"]); \
      print(F[1], \" order-level survivors, \", #F[2], \" fields\"); \
      print(finalfields(F[2], \"$R/final_fields.txt\"))" | $GP
