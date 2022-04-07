#!/bin/sh

F90=gen_be_ensmean.f90
F901=gen_be_ensmean_ref.f90


# ==========================================
PPREFIX=`basename ${F90} .f90`
PPREFIX1=`basename ${F901} .f90`

EXEC=${PPREFIX}.x
EXEC1=${PPREFIX1}.x

mpif90 -f90=ifort -O3 -assume byterecl -warn all -implicitnone -g -traceback -fp-model strict -qopenmp -module ${NETCDF}/include  ${F90} -o ${EXEC} -L${NETCDF}/lib -lnetcdf -lnetcdff -warn all,noexternal

mpif90 -f90=ifort -O3 -assume byterecl -warn all -implicitnone -g -traceback -fp-model strict -qopenmp -module ${NETCDF}/include  ${F901} -o ${EXEC1} -L${NETCDF}/lib -lnetcdf -lnetcdff -warn all,noexternal

#mpif90 -f90=ifort -mkl  -o ${EXEC}  -L${NETCDF}/lib -lnetcdf

