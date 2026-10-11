#!/bin/bash
# Local answer-check matrix for the perf driver (plan Stage 6.8, gates B2/B3).
#
#   run_perf_local.sh <case run dir: lnd_in, drv_in, lnd_modelio.nml> <elmxx_perf exe> [nsteps]
#
# 1. Reference: 1 cell on 1 rank writes ref.bin and a restart.
# 2. N = 1, 10, 100 cells on 1, 2, 4 ranks (ranks <= N): every copy must be
#    bit-identical to its first copy (in-run check) AND to ref.bin.
# 3. For each N, the restarts from every rank count must be identical
#    (ncdump at 17 digits; restarts are layout-independent, elmxxRestMod).
# 4. Negative control: 10 cells with cell 7 perturbed by 1e-9 K must FAIL.
#
# Runs in <run dir>/perf_check/, one subdirectory per run. Exit status 0
# only if every expected pass passes and the negative control fails.
set -u

RUNDIR=$(cd "$1" && pwd)
EXE=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
NSTEPS=${3:-48}
# Shared-memory transport only: on this laptop Open MPI's default transport
# selection hangs in the first collective (it does not in the coupled runs'
# environment). Override with MPIRUN=... elsewhere.
MPIRUN=${MPIRUN:-"mpirun --mca pml ob1 --mca btl self,sm"}
OUT=$RUNDIR/perf_check
rm -rf "$OUT" && mkdir -p "$OUT"
nfail=0

run() {   # run <name> <np> <args...>
  local name=$1 np=$2; shift 2
  mkdir -p "$OUT/$name"
  cp "$RUNDIR/lnd_in" "$RUNDIR/drv_in" "$RUNDIR/lnd_modelio.nml" "$OUT/$name/"
  (cd "$OUT/$name" && $MPIRUN -np "$np" "$EXE" --nsteps "$NSTEPS" "$@" > log 2>&1)
  echo $?
}

restart_dump() {  # canonical text of a restart: data section only, 17 digits
  ncdump -p 9,17 "$1" | sed -n '/^data:/,$p'
}

echo "== reference: 1 cell, 1 rank, $NSTEPS steps"
rc=$(run ref 1 --ncells 1 --write-reference "$OUT/ref.bin" --write-restart)
grep -E "PERF |PERF CHECK|write_reference" "$OUT/ref/log"
[ "$rc" = 0 ] || { echo "reference run FAILED (rc=$rc)"; tail -20 "$OUT/ref/log"; exit 1; }

for n in 1 10 100; do
  for np in 1 2 4; do
    [ "$np" -le "$n" ] || continue
    name=n${n}_np${np}
    rc=$(run $name $np --ncells $n --reference "$OUT/ref.bin" --write-restart)
    res=$(grep -E "PERF CHECK" "$OUT/$name/log")
    t=$(grep -E "^ PERF ncells" "$OUT/$name/log" | sed 's/.*t_step_max= *\([^ ]*\).*/\1/')
    printf "  %-10s rc=%s  %-18s t_step_max=%s\n" "$name" "$rc" "$res" "$t"
    if [ "$rc" != 0 ] || ! echo "$res" | grep -q PASS; then
      nfail=$((nfail+1)); grep -E "MISMATCH|MISSING|SHAPE|ERROR" "$OUT/$name/log" | head -5
    fi
  done
  # cross-layout: restarts identical for every rank count
  base=""
  for np in 1 2 4; do
    [ "$np" -le "$n" ] || continue
    f=$(ls "$OUT/n${n}_np${np}"/*.elmxx.r.*.nc 2>/dev/null | head -1)
    [ -n "$f" ] || { echo "  n=$n np=$np: no restart written"; nfail=$((nfail+1)); continue; }
    restart_dump "$f" > "$OUT/n${n}_np${np}.cdl"
    if [ -z "$base" ]; then base=$OUT/n${n}_np${np}.cdl; continue; fi
    if cmp -s "$base" "$OUT/n${n}_np${np}.cdl"; then
      echo "  restart n=$n: np=$np identical to np=1"
    else
      echo "  restart n=$n: np=$np DIFFERS from np=1 ($(diff "$base" "$OUT/n${n}_np${np}.cdl" | grep -c '^<') lines)"
      nfail=$((nfail+1))
    fi
  done
done

echo "== negative control: 10 cells, 2 ranks, cell 7 +1e-9 K (must FAIL)"
rc=$(run negctl 2 --ncells 10 --reference "$OUT/ref.bin" --perturb-cell 7)
grep -E "PERF CHECK|fields checked|fields against" "$OUT/negctl/log"
grep -E "MISMATCH" "$OUT/negctl/log" | head -3
if grep -q "PERF CHECK FAIL" "$OUT/negctl/log"; then
  echo "  negative control failed as it must"
else
  echo "  NEGATIVE CONTROL PASSED -- the checks are blind"; nfail=$((nfail+1))
fi

echo "== $nfail problem(s)"
[ "$nfail" = 0 ]
