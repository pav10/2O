#!/bin/sh
# run_q4_par.sh OUTDIR [NPROC] [JOBS] -- relative quartics over Q(sqrt5) (degree 8, Case B),
# one job per residue c1 = i + j*y mod 4, job number 4*i + j.
# The automorphism sqrt5 -> -sqrt5 (y -> 1 - y) maps c1 = i + j*y to (i+j) - j*y and gives the
# same field K, so only the 10 orbit representatives are needed:
#   j = 0: 0 4 8 12;  (0,1) 1, (0,2) 2, (0,3) 3, (1,1) 5, (1,2) 6, (2,1) 9.
OUT=$1; NP=${2:-4}; JOBS=${3:-"0 1 2 3 4 5 6 8 9 12"}
mkdir -p "$OUT"
# a residue whose log already reports its totals is skipped (resumable after an interruption)
for j in $JOBS; do grep -qs "candidates" "$OUT/log_$j.txt" || echo $j; done | xargs -P "$NP" -I{} sh -c "rm -f $OUT/q4_{}.txt; echo 'enum_q4(\"$OUT/q4_{}.txt\", {}, 16)' | gp -q lift6.gp rel6.gp field6.gp sharp.gp rel8.gp 2>&1 | grep -v Warning > $OUT/log_{}.txt"
