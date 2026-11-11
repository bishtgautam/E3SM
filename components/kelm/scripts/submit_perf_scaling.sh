#!/bin/bash

usage() {
  cat << USAGE_EOF
Usage: $0 [-h] [-n] <case dir>

Submit a strong scaling study of kelm_perf on Perlmutter. The machine comes
from the case, so run this once for a pm-cpu case and once for a pm-gpu case.
Create and build the cases with perf_brazil_kelm.sh.

Submits these jobs:
  reference  1 node: 1 cell on 1 rank writes ref.bin and a restart, then a
             negative control with one perturbed cell must fail the checks.
  scaling    One job for each node count, after the reference succeeds. Every
             cell count runs with every copy checked against ref.bin.

Each run fills a node: 128 ranks per pm-cpu node, 4 ranks (one per GPU) per
pm-gpu node. Results are written to <RUNDIR>/perf_scaling/<date>_<time>/, with
one results.nodes<K>.csv for each scaling job. Compare CPU and GPU answers
afterwards from the restarts in ref/ of each machine's results directory.

Options:
  -n, --dry-run  Print the sbatch commands without submitting them
  -h, --help     Show this help and exit

Settings, from the environment:
  CELLS="100000"  Cell counts. The same counts run on every node count
  NODES="1 2 4"   Node counts
  NDAYS=5         Simulated days, at 48 steps of 1800 s per day
  QOS=debug       Queue. debug allows at most 8 nodes
USAGE_EOF
}

# Set the per-machine values used by both the submission and the batch jobs.
set_machine() {
  case $1 in
    pm-cpu)
      ACCOUNT=m5263
      CONSTRAINT=cpu
      RANKS_PER_NODE=128
      CPUS_PER_RANK=2
      SBATCH_GPU_FLAGS=""
      SRUN_GPU_FLAGS=""
      ;;
    pm-gpu)
      ACCOUNT=m5263_g
      # gpu alone can land on 40 GB or 80 GB A100s. Pin one, so that a problem
      # size fits every time and the timings stay comparable.
      CONSTRAINT='gpu&hbm80g'
      RANKS_PER_NODE=4
      CPUS_PER_RANK=32
      SBATCH_GPU_FLAGS="--gpus-per-node=4"
      SRUN_GPU_FLAGS="--gpus-per-task=1"
      ;;
    *)
      echo "ERROR: unsupported machine '$1': expected pm-cpu or pm-gpu" >&2
      exit 1
      ;;
  esac
}

# run_perf <run name> <number of ranks> <kelm_perf arguments...>
# Each run has its own directory, so concurrent jobs cannot collide.
run_perf() {
  local name=$1 ntasks=$2
  shift 2
  mkdir -p "${OUT}/${name}"
  cp "${RUNDIR}/lnd_in" "${RUNDIR}/drv_in" "${RUNDIR}/lnd_modelio.nml" "${OUT}/${name}/"
  (
    cd "${OUT}/${name}" &&
      srun -n "${ntasks}" -c "${CPUS_PER_RANK}" --cpu-bind=cores ${SRUN_GPU_FLAGS} \
        "${EXE}" --nsteps "${NSTEPS}" "$@" > log 2>&1
  )
}

# perf_value <log> <field>: a value from the driver's one-line PERF summary.
perf_value() {
  grep "PERF ncells=" "$1" | sed -n "s/.* $2= *\([^ ]*\).*/\1/p"
}

job_reference() {
  local rc
  echo "== reference: 1 cell, 1 rank, ${NSTEPS} steps"
  run_perf ref 1 --ncells 1 --write-reference "${OUT}/ref.bin" --write-restart
  rc=$?
  grep -E "PERF ncells=|PERF CHECK" "${OUT}/ref/log"
  if [ "${rc}" != 0 ] || [ ! -s "${OUT}/ref.bin" ]; then
    echo "ERROR: reference run failed (rc=${rc}), see ${OUT}/ref/log"
    return 1
  fi

  echo "== negative control: 10 cells, 2 ranks, cell 7 perturbed (must FAIL)"
  run_perf negctl 2 --ncells 10 --reference "${OUT}/ref.bin" --perturb-cell 7
  if ! grep -q "PERF CHECK FAIL" "${OUT}/negctl/log"; then
    echo "ERROR: negative control passed, so the checks are blind"
    return 1
  fi
  echo "negative control failed as it must"
}

job_scaling() {
  local nodes=$1 cells=$2 ranks n name log rc check nfail=0
  local csv=${OUT}/results.nodes${nodes}.csv
  ranks=$((nodes * RANKS_PER_NODE))
  echo "machine,nodes,ranks,ncells,ntimed,t_step_max,t_step_mean,t_init,check,rc,jobid" > "${csv}"

  # Results are written after each run, so a job that reaches its time limit
  # keeps what it finished.
  for n in ${cells}; do
    name=n${n}_nodes${nodes}
    log=${OUT}/${name}/log
    echo "== ${name}: ${ranks} ranks, ${NSTEPS} steps"
    run_perf "${name}" "${ranks}" --ncells "${n}" --reference "${OUT}/ref.bin"
    rc=$?
    check=FAIL
    grep -q "PERF CHECK PASS" "${log}" && check=PASS
    echo "${MACH},${nodes},${ranks},${n},$(perf_value "${log}" ntimed)," \
      "$(perf_value "${log}" t_step_max),$(perf_value "${log}" t_step_mean)," \
      "$(perf_value "${log}" t_init),${check},${rc},${SLURM_JOB_ID}" |
      tr -d ' ' >> "${csv}"
    tail -n 1 "${csv}"
    if [ "${rc}" != 0 ] || [ "${check}" != PASS ]; then
      nfail=$((nfail + 1))
      grep -E "MISMATCH|MISSING|SHAPE|ERROR" "${log}" | head -5
    fi
  done

  echo "== ${nfail} problem(s)"
  [ "${nfail}" = 0 ]
}

# Batch entry point: <reference|scaling> <case dir> <output dir> <machine>
#                    <exe> <rundir> <nsteps> [<nodes> <cells>]
run_job() {
  local kind=$1 case_dir=$2
  OUT=$3
  MACH=$4
  EXE=$5
  RUNDIR=$6
  NSTEPS=$7
  set_machine "${MACH}"
  cd "${case_dir}" || exit 1
  . ./.env_mach_specific.sh
  case ${kind} in
    reference) job_reference ;;
    scaling) job_scaling "$8" "$9" ;;
    *)
      echo "ERROR: unknown job '${kind}'" >&2
      exit 1
      ;;
  esac
}

# submit <job name> <nodes> <time limit> <sbatch flags...> -- <job arguments...>
# Prints the job id. In a dry run, prints the command and a placeholder id.
submit() {
  local name=$1 nodes=$2 walltime=$3
  shift 3
  local flags=()
  while [ "$1" != "--" ]; do
    flags+=("$1")
    shift
  done
  shift
  local cmd=(sbatch --parsable --job-name="${name}" --account="${ACCOUNT}" --qos="${QOS}"
    --constraint="${CONSTRAINT}" --nodes="${nodes}" --time="${walltime}"
    --output="${OUT}/${name}.%j.out" ${SBATCH_GPU_FLAGS} "${flags[@]}" "${SELF}" --job "$@")
  if [ "${DRY_RUN}" = 1 ]; then
    echo "${cmd[*]}" >&2
    echo "<${name}>"
  else
    "${cmd[@]}" | cut -d';' -f1
  fi
}

if [ "${1:-}" = "--job" ]; then
  shift
  run_job "$@"
  exit $?
fi

DRY_RUN=0
CASE_DIR=""
while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -n|--dry-run) DRY_RUN=1 ;;
    -*)
      usage >&2
      exit 1
      ;;
    *)
      if [ -n "${CASE_DIR}" ]; then
        usage >&2
        exit 1
      fi
      CASE_DIR=$1
      ;;
  esac
  shift
done
if [ -z "${CASE_DIR}" ]; then
  usage >&2
  exit 1
fi

CELLS=${CELLS:-100000}
NODES=${NODES:-1 2 4}
NDAYS=${NDAYS:-5}
QOS=${QOS:-debug}
NSTEPS=$((NDAYS * 48))

SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
CASE_DIR=$(cd "${CASE_DIR}" && pwd) || exit 1
cd "${CASE_DIR}" || exit 1
MACH=$(./xmlquery --value MACH) || exit 1
EXEROOT=$(./xmlquery --value EXEROOT) || exit 1
RUNDIR=$(./xmlquery --value RUNDIR) || exit 1
set_machine "${MACH}"

EXE=${EXEROOT}/kelm_perf
if [ ! -x "${EXE}" ]; then
  echo "ERROR: ${EXE} not found: build it with perf_brazil_kelm.sh" >&2
  exit 1
fi
for f in lnd_in drv_in lnd_modelio.nml; do
  if [ ! -f "${RUNDIR}/${f}" ]; then
    echo "ERROR: ${RUNDIR}/${f} not found: run the case once to create it" >&2
    exit 1
  fi
done
for k in ${NODES}; do
  for n in ${CELLS}; do
    if [ "${n}" -lt $((k * RANKS_PER_NODE)) ]; then
      echo "ERROR: ${n} cells is fewer than the ranks on ${k} node(s)" >&2
      exit 1
    fi
  done
done

OUT=${RUNDIR}/perf_scaling/$(date +%Y%m%d_%H%M%S)
[ "${DRY_RUN}" = 1 ] || mkdir -p "${OUT}" || exit 1

JOB_ARGS=("${CASE_DIR}" "${OUT}" "${MACH}" "${EXE}" "${RUNDIR}" "${NSTEPS}")
REF_JOB=$(submit perf_ref 1 00:10:00 -- reference "${JOB_ARGS[@]}") || exit 1
[ -n "${REF_JOB}" ] || exit 1
echo "reference job: ${REF_JOB}"
for k in ${NODES}; do
  JOB=$(submit "perf_nodes${k}" "${k}" 00:30:00 \
    --dependency="afterok:${REF_JOB}" --kill-on-invalid-dep=yes \
    -- scaling "${JOB_ARGS[@]}" "${k}" "${CELLS}") || exit 1
  echo "scaling job on ${k} node(s): ${JOB}"
done
echo "results: ${OUT}"
