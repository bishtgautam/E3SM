#!/bin/sh

cd cime/scripts
./create_newcase --case f19_g16.ICLM45 --res f19_g16 --compset ICLM45 --mach gitlab-ci-linux --compiler gnu
cd f19_g16.ICLM45

./xmlchange DATM_CLMNCEP_YR_END=1972
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUNDIR=${PWD}/run
./xmlchange EXEROOT=${PWD}/bld
./xmlchange NTASKS=1
./xmlchange DIN_LOC_ROOT=$PWD
./xmlchange MPILIB=mpich

./case.setup

./case.build

