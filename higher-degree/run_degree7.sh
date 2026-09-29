#!/bin/sh
# run_degree7.sh -- reproduce the degree-7 computation end to end.
# A septic field has no proper subfields, so only Case A (x_1 generates K) occurs.
# Requirements: gcc, PARI/GP >= 2.15, python3.  Running time: about 1 minute on 4 cores for the
# enumeration (7.6 * 10^6 leaves after tree pruning; 5.3 * 10^8 without), a few minutes for the rest.
set -e
cd "$(dirname "$0")/src"
R=../results/d7
GP="gp -q lift6.gp order_filter.gp rel6.gp field6.gp condN.gp witness.gp pipeline.gp"
mkdir -p $R

echo "== build"
gcc -O2 -DDEG=7 -o enum_d7 enum6.c -lm

echo "== enumeration with the strong filters (downspread, quadratic, NS, RIG in Z[x_1])"
rm -rf $R/enum
./run_enum_par.sh 7 $R/enum strong 4
cat $R/enum/out_*.txt > $R/d7_strong_raw.txt
cat $R/enum/log_* | python3 -c "
import re, sys
tot = {}
for line in sys.stdin:
    for k, v in re.findall(r'(\w+)=(\d+)', line):
        if k not in ('a1', 'a2', 'S'): tot[k] = tot.get(k, 0) + int(v)
print(tot)"

echo "== exact order-level tests (NS / BC / RG) on the survivors"
for p in 0 1 2 3; do
  echo "primitive(\"$R/d7_strong_raw.txt\", \"$R/prim7_status_$p.tmp\", $p, 4)" | $GP > /dev/null &
done
wait
cat $R/prim7_status_*.tmp > $R/prim7_status.txt; rm -f $R/prim7_status_*.tmp
echo "irreducible: $(wc -l < $R/prim7_status.txt)"
for s in NS BC RG OK; do echo "  $s: $(grep -c "\"$s\"" $R/prim7_status.txt)"; done
