#!/bin/sh
# Census of condition (1.2) of notes/etype-zeta.tex.  Usage: sh src/etype_census.sh FIELDDIR [JOBS]
# FIELDDIR holds deg3.txt, deg4.txt, deg5.txt in the format of pav10/variance data/fields (header line + label|coeffs|...).
# About 15 minutes on 4 cores.  Outputs go to results/etype/.
set -e
F=${1:?field directory}; J=${2:-4}; OUT=results/etype; TMP=$(mktemp -d)
for d in 3 4 5; do sed -n '2,$p' "$F/deg$d.txt" | awk -F'|' '{print "["$2"]"}' > "$TMP/d$d.txt"; done
s=0
while [ $s -lt "$J" ]; do
  echo "runP(\"$TMP/d3.txt\",$s,$J,\"$TMP/P3_$s\",40000); runP(\"$TMP/d4.txt\",$s,$J,\"$TMP/P4_$s\",10^9); runP(\"$TMP/d5.txt\",$s,$J,\"$TMP/P5_$s\",10^9); runNP(\"$TMP/d3.txt\",$s,$J,\"$TMP/N3_$s\"); runNP(\"$TMP/d4.txt\",$s,$J,\"$TMP/N4_$s\"); runNP(\"$TMP/d5.txt\",$s,$J,\"$TMP/N5_$s\");" \
    | gp -q src/etype_zeta.gp src/etype_census.gp > /dev/null 2>&1 &
  s=$((s + 1))
done
echo "runCert(\"$TMP/d3.txt\", 3621, \"$OUT/cert_cubic_disc_le_3621.txt\");" | gp -q src/etype_zeta.gp src/etype_census.gp > /dev/null 2>&1
wait
for d in 3 4 5; do cat "$TMP"/P${d}_* | grep -v '#done' | sort -n > "$OUT/condP_deg$d.txt"; cat "$TMP"/N${d}_* | grep -v '#done' | sort -n > "$OUT/nonprincipal_deg$d.txt"; done
rm -rf "$TMP"
