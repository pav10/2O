#!/bin/sh
# run_c3_par.sh OUTDIR [NPROC] -- degree 9, Case B: relative cubics over F = Q(zeta_7)^+
# (rel9.gp), one job per residue b1 mod 3 O_F (27 jobs).  Resumable: a job whose log already
# reports its totals is skipped, and an unfinished job is restarted from scratch.
OUT=$1; NP=${2:-4}
mkdir -p "$OUT"
for j in $(seq 0 26); do grep -qs "relative cubics" "$OUT/log_$j.txt" || echo $j; done | xargs -P "$NP" -I{} sh -c "rm -f $OUT/c3_{}.txt; echo 'enum_c3(\"Q(z7)+\", y^3-y^2-2*y+1, 14/9, 1.5513, \"$OUT/c3_{}.txt\", {}, 27)' | gp -q lift6.gp rel6.gp field6.gp sharp.gp rel8.gp rel9.gp 2>&1 | grep -v Warning > $OUT/log_{}.txt"
