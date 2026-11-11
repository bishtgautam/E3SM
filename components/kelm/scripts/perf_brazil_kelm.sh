#!/bin/sh

usage() {
  cat << USAGE_EOF
Usage: $0 [-h] [pm-gpu|pm-cpu]

Create and build a 1x1_brazil kELM case on Perlmutter, then build the
kelm_perf performance driver, which is written to \${EXEROOT}/kelm_perf.
The case is created in cime/scripts.

Arguments:
  pm-gpu      Build with gnugpu so kelm_perf runs on GPUs (default)
  pm-cpu      Build with gnu so kelm_perf runs on CPUs
  -h, --help  Show this help and exit
USAGE_EOF
}

MACH=${1:-pm-gpu}
case ${MACH} in
  -h|--help)
    usage
    exit 0
    ;;
  pm-gpu)
    COMPILER=gnugpu
    PROJECT=e3sm_g
    ;;
  pm-cpu)
    COMPILER=gnu
    PROJECT=e3sm
    ;;
  *)
    usage >&2
    exit 1
    ;;
esac

RES=1x1_brazil
COMPSET=I1850KELM
DRIVER=mct

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SRC_DIR=$(cd "${SCRIPT_DIR}/../../.." && pwd)
GIT_HASH=$(git -C ${SRC_DIR} rev-parse --short HEAD)

CASE_DIR=${SRC_DIR}/cime/scripts
CASE_NAME=perf_brazil.${RES}.${COMPSET}.${MACH}.${COMPILER}.${GIT_HASH}

cd ${SRC_DIR}/cime/scripts
./create_newcase -case ${CASE_DIR}/${CASE_NAME} \
-res ${RES} -mach ${MACH} -compiler ${COMPILER} -compset ${COMPSET} --driver ${DRIVER} \
--project ${PROJECT}

cd ${CASE_DIR}/${CASE_NAME}
./xmlchange NTASKS=1
./xmlchange JOB_QUEUE=debug,JOB_WALLCLOCK_TIME=00:10:00
./xmlchange STOP_OPTION=ndays
./xmlchange STOP_N=2

INP=/global/cfs/projectdirs/m5263/gbisht/e3sm-elmxx/cime/scripts
cat >> user_nl_kelm << NL_EOF
fsurdat = '${INP}/surfdata_1x1_brazil_simyr1850_c140610.nc'
fparamfile = '${INP}/clm_params_c180301.nc'
fsnowoptics = '${INP}/snicar_optics_5bnd_mam_c160322.nc'
fsnowaging = '${INP}/snicar_drdt_bst_fit_60_c070416.nc'
kelm_kernels = 'begwaterbal,canhydro,cansunshade,surfrad,cantemp,baregrnd,canflux,soiltemp,soilflux,snowwater,surfrunoff,rootwater,soilwater,watertable,hydrodrain,snowlayers,snowage,balancecheck'
kelm_do_albedo = .true.
kelm_do_photosynthesis = .true.
kelm_check_boundary = .false.
NL_EOF

./case.setup
./case.build || exit 1

# kelm_perf is an extra target in the case's CMake build, not part of
# case.build. Use cmake --build so this works with either the make or the
# ninja backend. It writes ${EXEROOT}/kelm_perf.
. ./.env_mach_specific.sh
EXEROOT=$(./xmlquery --value EXEROOT)
GMAKE_J=$(./xmlquery --value GMAKE_J)
cd ${EXEROOT}/cmake-bld
cmake --build . --target kelm_perf --parallel ${GMAKE_J}
