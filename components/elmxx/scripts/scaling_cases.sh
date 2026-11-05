#!/bin/bash

usage() {
  cat << USAGE_EOF
Usage: $0 [-h] [-m elmxx|elm] build <pm-cpu|pm-gpu>
       $0 [-h] [-m elmxx|elm] case  <pm-cpu|pm-gpu> <ncells> <nodes>

Strong scaling cases for ELMxx inside E3SM, on replicated Brazil cells. Make
the inputs first with make_replicated_inputs.py.

Commands:
  build  Create and build one case for the machine. Every scaling case of that
         machine and model runs its e3sm.exe, so build once.
  case   Create a case of <ncells> on <nodes> full nodes, using that build.
         Submit it with case.submit, adding --constraint for 80 GB GPUs:
           ./case.submit -a="--constraint='gpu&hbm80g'"

Options:
  -m MODEL    elmxx (default) or the ELM twin, I1850ELM on pm-cpu
  -h, --help  Show this help and exit

Settings, from the environment:
  INPUTS=\$PSCRATCH/e3sm_scratch/elmxx_scaling/inputs  Replicated inputs
  NDAYS=5             Simulated days
  QUEUE=debug         Queue
  WALLCLOCK=00:30:00  Time limit
  REST_OPTION=never   Restarts. Set to ndays to write one at the end of the run
  RANKS=N             Run on N ranks instead of full nodes, for small tests
USAGE_EOF
}

MODEL=elmxx
while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
      ;;
    -m)
      MODEL=$2
      shift 2
      ;;
    *) break ;;
  esac
done
COMMAND=$1
MACH=$2
NCELLS=$3
NODES=$4

case ${MODEL} in
  elmxx) COMPSET=I1850ELMXX ;;
  elm) COMPSET=I1850ELM ;;
  *)
    usage >&2
    exit 1
    ;;
esac
case ${MACH} in
  pm-cpu)
    COMPILER=gnu
    PROJECT=m5263
    RANKS_PER_NODE=128
    ;;
  pm-gpu)
    COMPILER=gnugpu
    PROJECT=m5263_g
    RANKS_PER_NODE=4
    ;;
  *)
    usage >&2
    exit 1
    ;;
esac
if [ "${COMMAND}" != build ] && { [ "${COMMAND}" != case ] || [ -z "${NODES}" ]; }; then
  usage >&2
  exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SRC_DIR=$(cd "${SCRIPT_DIR}/../../.." && pwd)
GIT_HASH=$(git -C "${SRC_DIR}" rev-parse --short HEAD)
CASE_ROOT=${SRC_DIR}/cime/scripts
INPUTS=${INPUTS:-${PSCRATCH}/e3sm_scratch/elmxx_scaling/inputs}
NDAYS=${NDAYS:-5}
QUEUE=${QUEUE:-debug}
WALLCLOCK=${WALLCLOCK:-00:30:00}
REST_OPTION=${REST_OPTION:-never}
BUILD_CASE=${CASE_ROOT}/scaling_build.${MODEL}.${MACH}.${COMPILER}.${GIT_HASH}

new_case() {
  "${CASE_ROOT}/create_newcase" --case "$1" --res 1x1_brazil --mach "${MACH}" \
    --compiler "${COMPILER}" --compset "${COMPSET}" --driver mct --project "${PROJECT}" || exit 1
  cd "$1" || exit 1
  ./xmlchange NTHRDS=1
}

# The parameter files come from the case's defaults for ELM, as its twin needs.
write_namelist() {
  local surfdata=$1
  if [ "${MODEL}" = elm ]; then
    echo "fsurdat = '${surfdata}'" >> user_nl_elm
    return
  fi
  cat >> user_nl_elmxx << NL_EOF
fsurdat = '${surfdata}'
fparamfile = '${CASE_ROOT}/clm_params_c180301.nc'
fsnowoptics = '${CASE_ROOT}/snicar_optics_5bnd_mam_c160322.nc'
fsnowaging = '${CASE_ROOT}/snicar_drdt_bst_fit_60_c070416.nc'
elmxx_kernels = 'begwaterbal,canhydro,cansunshade,surfrad,cantemp,baregrnd,canflux,soiltemp,soilflux,snowwater,surfrunoff,rootwater,soilwater,watertable,hydrodrain,snowlayers,snowage,balancecheck'
elmxx_do_albedo = .true.
elmxx_do_photosynthesis = .true.
elmxx_check_boundary = .false.
NL_EOF
}

if [ "${COMMAND}" = build ]; then
  new_case "${BUILD_CASE}"
  ./xmlchange NTASKS=1
  ./case.setup || exit 1
  ./case.build
  exit $?
fi

DOMAIN=domain.lnd.brazil_${NCELLS}.nc
SURFDATA=${INPUTS}/surfdata_brazil_${NCELLS}.nc
if [ ! -f "${INPUTS}/${DOMAIN}" ] || [ ! -f "${SURFDATA}" ]; then
  echo "ERROR: inputs for ${NCELLS} cells not in ${INPUTS}: run make_replicated_inputs.py" >&2
  exit 1
fi
EXEROOT=$(cd "${BUILD_CASE}" 2>/dev/null && ./xmlquery --value EXEROOT)
if [ ! -x "${EXEROOT}/e3sm.exe" ]; then
  echo "ERROR: no e3sm.exe from ${BUILD_CASE}: run '$0 -m ${MODEL} build ${MACH}'" >&2
  exit 1
fi

RANKS=${RANKS:-$((NODES * RANKS_PER_NODE))}
if [ "${RANKS}" -gt "${NCELLS}" ]; then
  echo "ERROR: ${RANKS} ranks is more than ${NCELLS} cells" >&2
  exit 1
fi

CASE=${CASE_ROOT}/scaling.${MODEL}.${MACH}.n${NCELLS}.nodes${NODES}.r${RANKS}.${GIT_HASH}
new_case "${CASE}"
./xmlchange NTASKS="${RANKS}"
./xmlchange MAX_MPITASKS_PER_NODE="${RANKS_PER_NODE}"
./xmlchange LND_DOMAIN_PATH="${INPUTS}",LND_DOMAIN_FILE="${DOMAIN}"
./xmlchange ATM_DOMAIN_PATH="${INPUTS}",ATM_DOMAIN_FILE="${DOMAIN}"
./xmlchange LND_NX="${NCELLS}",LND_NY=1,ATM_NX="${NCELLS}",ATM_NY=1
./xmlchange STOP_OPTION=ndays,STOP_N="${NDAYS}"
./xmlchange REST_OPTION="${REST_OPTION}",REST_N="${NDAYS}",DOUT_S=FALSE
./xmlchange JOB_QUEUE="${QUEUE}",JOB_WALLCLOCK_TIME="${WALLCLOCK}"
./xmlchange TIMER_LEVEL=10
./xmlchange EXEROOT="${EXEROOT}",BUILD_COMPLETE=TRUE

# GPU-aware MPI slows PIO's rearranger about 100 times on host buffers, and
# ELMxx passes no device buffers to MPI.
if [ "${MACH}" = pm-gpu ]; then
  sed -i 's|\(<env name="MPICH_GPU_SUPPORT_ENABLED">\)1<|\10<|' env_mach_specific.xml
fi
write_namelist "${SURFDATA}"
./case.setup || exit 1
echo "case: ${CASE}"
