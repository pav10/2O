#!/bin/sh
# validate.sh D -- independent checks of the degree-D computation (D = 3..7):
#  (1) every exclusion made by the in-enumerator tests (NS / RIG, mode "strong") on an
#      irreducible polynomial is re-checked in exact arithmetic (verify_strong.gp);
#  (2) enum6.c is compared leaf by leaf with the exact re-implementation enum_check.gp
#      (all subtrees for D = 6; the subtrees with the fewest leaves, 2.2 million leaves
#      in total, for D = 7).
# Running time: D = 6 about 45 minutes, D = 7 about 1.5 hours on 4 cores.
set -e
D=$1
cd "$(dirname "$0")/src"
V=../results/validate_d$D; rm -rf $V; mkdir -p $V
gcc -O2 -DDEG=$D -o enum_d$D enum6.c -lm

echo "== (1) exact re-check of the strong-mode exclusions"
./run_enum_par.sh $D $V/why why 4
cat $V/why/out_*.txt > $V/why.txt
for p in 0 1 2 3; do
  echo "print(vstrong(\"$V/why.txt\", $p, 4))" | gp -q lift6.gp verify_strong.gp > $V/vstrong_$p.txt &
done
wait
cat $V/vstrong_*.txt | tr -d '[]' | awk -F, '{c+=$1; r+=$2; b+=$3} END{print "checked (irreducible):", c, " reducible:", r, " disagreements:", b}'

echo "== (2) enumerator versus exact re-implementation"
# leaf counts per subtree are in the logs of step (1)
cat $V/why/log_* | awk '{split($1,a,"="); split($2,b,"="); split($4,c,"="); print c[2], a[2], b[2]}' \
  | sort -n | awk -v D=$D '{ s += $1; if (D == 6 || s < 2600000) print $2, $3 }' > $V/jobs.txt
rm -rf $V/why $V/why.txt
cat $V/jobs.txt | xargs -P 4 -n 2 sh -c "
  ./enum_d$D \$0 \$1 0 all 2>/dev/null | sort > $V/c_\$0_\$1.txt
  echo \"v=encheck($D,\$0,\$1); for(i=1,#v, print(strjoin(apply(x->Str(x),v[i]),\\\" \\\")))\" | gp -q enum_check.gp | sort > $V/g_\$0_\$1.txt
  if cmp -s $V/c_\$0_\$1.txt $V/g_\$0_\$1.txt; then r=SAME; else r=DIFF; fi
  echo \"\$0 \$1 \$(wc -l < $V/c_\$0_\$1.txt) \$r\" >> $V/enum_compare.txt
  rm -f $V/c_\$0_\$1.txt $V/g_\$0_\$1.txt"
awk '{n+=$3; s[$4]++} END{print "subtrees:", NR, " leaves:", n; for (k in s) print "  ", k, s[k]}' $V/enum_compare.txt
