#!/bin/sh
# run_enum.sh D OUT : run the degree-D enumerator over all (a1,a2), a1 in [-floor(D/2),0]
D=$1; OUT=$2; : > $OUT
for a1 in $(seq 0 -1 -$((D/2))); do
  for a2 in $(seq -40 20); do ./enum_d$D $a1 $a2 >> $OUT; done
done
