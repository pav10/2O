#!/bin/sh
# run_split_par.sh D OUTDIR [MODE] [NPROC] : split enumeration (modes X, YL, YU per (a1,a2)),
# one job per (a1,a2), largest first.  MODE = sharp (default) or strong.
D=$1; OUT=$2; MODE=${3:-sharp}; NP=${4:-4}
mkdir -p $OUT
python3 - "$D" > $OUT/jobs.txt <<'PY'
import math, sys
d = int(sys.argv[1])
smax = d * ((d + 4) + math.sqrt((d + 4) ** 2 - 16)) / 4
jobs = []
for a1 in range(0, -(d // 2) - 1, -1):
    for a2 in range(-80, 40):
        S = (d - 1) * a1 * a1 / d - 2 * a2
        if 0 < S <= smax + 1e-9:
            jobs.append((S, a1, a2))
jobs.sort(reverse=True)
for S, a1, a2 in jobs:
    print(a1, a2)
PY
cat $OUT/jobs.txt | xargs -P $NP -n 2 sh -c "for m in X YL YU; do ./enum_d$D \$0 \$1 0 $MODE \$m; done > $OUT/out_\$0_\$1.txt 2> $OUT/log_\$0_\$1.txt"
