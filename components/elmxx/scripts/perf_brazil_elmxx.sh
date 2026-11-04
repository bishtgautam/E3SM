#!/bin/sh

usage() {
  cat << USAGE_EOF
Usage: $0 [-h] [pm-gpu|pm-cpu]

Create and build a 1x1_brazil ELMxx case on Perlmutter, then build the
elmxx_perf performance driver, which is written to \${EXEROOT}/elmxx_perf.
The case is created in cime/scripts.

Arguments:
  pm-gpu      Build with gnugpu so elmxx_perf runs on GPUs (default)
  pm-cpu      Build with gnu so elmxx_perf runs on CPUs
  -h, --help  Show this help and exit
USAGE_EOF
}

MACH=${1:-pm-gpu}
case ${MACH} in
  -h|--help)
    usage
    exit 0
    ;;
  pm-gpu) COMPILER=gnugpu ;;
  pm-cpu) COMPILER=gnu ;;
  *)
    usage >&2
    exit 1
    ;;
esac

RES=1x1_brazil
COMPSET=I1850ELMXX
DRIVER=mct

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
SRC_DIR=$(cd "${SCRIPT_DIR}/../../.." && pwd)
GIT_HASH=$(git -C ${SRC_DIR} rev-parse --short HEAD)

CASE_DIR=${SRC_DIR}/cime/scripts
CASE_NAME=perf_brazil.${RES}.${COMPSET}.${MACH}.${COMPILER}.${GIT_HASH}

cd ${SRC_DIR}/cime/scripts
./create_newcase -case ${CASE_DIR}/${CASE_NAME} \
-res ${RES} -mach ${MACH} -compiler ${COMPILER} -compset ${COMPSET} --driver ${DRIVER}

cd ${CASE_DIR}/${CASE_NAME}
./xmlchange NTASKS=1
./xmlchange STOP_OPTION=ndays
./xmlchange STOP_N=2

INP=/global/cfs/projectdirs/m5263/gbisht/e3sm-elmxx/cime/scripts
cat >> user_nl_elmxx << NL_EOF
fsurdat = '${INP}/surfdata_1x1_brazil_simyr1850_c140610.nc'
fparamfile = '${INP}/clm_params_c180301.nc'
fsnowoptics = '${INP}/snicar_optics_5bnd_mam_c160322.nc'
fsnowaging = '${INP}/snicar_drdt_bst_fit_60_c070416.nc'
elmxx_kernels = 'begwaterbal,canhydro,cansunshade,surfrad,cantemp,baregrnd,canflux,soiltemp,soilflux,snowwater,surfrunoff,rootwater,soilwater,watertable,hydrodrain,snowlayers,snowage,balancecheck'
elmxx_do_albedo = .true.
elmxx_do_photosynthesis = .true.
elmxx_check_boundary = .false.
NL_EOF

./case.setup
./case.build || exit 1

# elmxx_perf is an extra target in the case's CMake build, not part of
# case.build. Use cmake --build so this works with either the make or the
# ninja backend. It writes ${EXEROOT}/elmxx_perf.
. ./.env_mach_specific.sh
EXEROOT=$(./xmlquery --value EXEROOT)
GMAKE_J=$(./xmlquery --value GMAKE_J)
cd ${EXEROOT}/cmake-bld
cmake --build . --target elmxx_perf --parallel ${GMAKE_J}
