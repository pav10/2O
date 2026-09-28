#!/bin/sh
# run_enum.sh D OUT [MODE] : run the degree-D enumerator over all (a1,a2), a1 in [-floor(D/2),0]
# MODE: empty (downspread + quadratic filters), "strong" (+ NS/RIG tests in Z[x_1]).
D=$1; OUT=$2; MODE=$3; : > $OUT
for a1 in $(seq 0 -1 -$((D/2))); do
  for a2 in $(seq -40 20); do
    if [ -n "$MODE" ]; then ./enum_d$D $a1 $a2 0 $MODE >> $OUT; else ./enum_d$D $a1 $a2 >> $OUT; fi
  done
done
