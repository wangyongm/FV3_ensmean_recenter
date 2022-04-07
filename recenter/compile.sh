#!/bin/sh

F90=gen_be_ensmeanrecenter.f90

# ==========================================
PPREFIX=`basename ${F90} .f90`

EXEC=${PPREFIX}.x

mpif90 -f90=ifort -O3 -r8 -warn all -implicitnone -g -traceback -fp-model strict -qopenmp -module ${NETCDF}/include  ${F90} -o ${EXEC} -L${NETCDF}/lib -lnetcdf -lnetcdff # -warn all,noexternal

#mpif90 -f90=ifort -mkl  -o ${EXEC}  -L${NETCDF}/lib -lnetcdf

