#!/bin/sh
# compare_split.sh D SMIN SMAX OUTDIR [NPROC] [TIMEOUT]
# Runs the plain enumerator (mode sharp) and the rung split (modes X, YL, YU) on every subtree
# (a1, a2) with SMIN < S <= SMAX (largest S first), and certifies each polynomial found by only
# one of the two runs with verify_split.gp (reducible / no square root / square root fails
# budget-cost).  Each subtree is compared as soon as it finishes (op_*.txt, os_*.txt), so
# partial results survive an interruption.  Subtrees whose plain run exceeds TIMEOUT seconds
# (default 3600) are reported and skipped.
D=$1; SMIN=$2; SMAX=$3; OUT=$4; NP=${5:-4}; TO=${6:-3600}
cd "$(dirname "$0")"
mkdir -p "$OUT"
python3 -c "
d = $D
jobs = []
for a1 in range(0, -(d // 2) - 1, -1):
    for a2 in range(-80, 40):
        S = (d - 1) * a1 * a1 / d - 2 * a2
        if $SMIN + 1e-9 < S <= $SMAX + 1e-9: jobs.append((S, a1, a2))
for S, a1, a2 in sorted(jobs, reverse=True): print(a1, a2)" > "$OUT/jobs.txt"
cat "$OUT/jobs.txt" | xargs -P "$NP" -n 2 sh -c "
  t0=\$(date +%s)
  if timeout $TO ./enum_d$D \$0 \$1 0 sharp > $OUT/b_\$0_\$1.txt 2>/dev/null; then r=done; else r=TIMEOUT; fi
  for m in X YL YU; do ./enum_d$D \$0 \$1 0 sharp \$m; done > $OUT/s_\$0_\$1.txt 2>/dev/null
  if [ \$r = done ]; then
    sort -u $OUT/b_\$0_\$1.txt -o $OUT/b_\$0_\$1.txt; sort -u $OUT/s_\$0_\$1.txt -o $OUT/s_\$0_\$1.txt
    comm -23 $OUT/b_\$0_\$1.txt $OUT/s_\$0_\$1.txt > $OUT/op_\$0_\$1.txt
    comm -13 $OUT/b_\$0_\$1.txt $OUT/s_\$0_\$1.txt > $OUT/os_\$0_\$1.txt
    np=\$(wc -l < $OUT/b_\$0_\$1.txt); ns=\$(wc -l < $OUT/s_\$0_\$1.txt)
  else np=-; ns=-; fi
  rm -f $OUT/b_\$0_\$1.txt $OUT/s_\$0_\$1.txt
  echo \"\$0 \$1 plain=\$r \$((\$(date +%s) - t0))s plain=\$np split=\$ns\" >> $OUT/times.txt"
cat "$OUT"/op_*.txt 2>/dev/null | sort -u > "$OUT/only_plain.txt"
cat "$OUT"/os_*.txt 2>/dev/null | sort -u > "$OUT/only_split.txt"
rm -f "$OUT"/op_*.txt "$OUT"/os_*.txt
{
  echo "subtrees $(wc -l < "$OUT/jobs.txt"), compared $(grep -c plain=done "$OUT/times.txt"), timed out $(grep -c TIMEOUT "$OUT/times.txt")"
  echo "only-plain $(wc -l < "$OUT/only_plain.txt")  only-split $(wc -l < "$OUT/only_split.txt")"
  echo "print(\"only-plain: \", checkfile(\"$OUT/only_plain.txt\")); print(\"only-split: \", checkfile(\"$OUT/only_split.txt\"))" | gp -q verify_split.gp
} > "$OUT/summary.txt"
cat "$OUT/summary.txt"
