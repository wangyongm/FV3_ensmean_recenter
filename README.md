# FV3_ensmean_recenter
 Run compile.sh to obtain the corresponding executables

 Usage:

 To calculate ensemble mean:

  1) Link all ensemble members
  imem=1
  while [ ${imem} -le ${ens_size} ]; do
    memstr4=`printf %04i ${imem}`
    memstr3=`printf %03i ${imem}`
    FG_MEM=${WORK_ROOT}/fv3prd_mem${memstr4}/ANA/fv_core.res.tile1_new.nc
    ${LN} -sf ${FG_MEM} fv3sar_tile1_mem${memstr3}_dynvar
    FG_MEM=${WORK_ROOT}/fv3prd_mem${memstr4}/ANA/fv_tracer.res.tile1_new.nc
    ${LN} -sf ${FG_MEM} fv3sar_tile1_mem${memstr3}_tracer
    FG_MEM=${WORK_ROOT}/fv3prd_mem${memstr4}/ANA/phy_data.nc
    ${LN} -sf ${FG_MEM} fv3sar_tile1_mem${memstr3}_phyvar

    (( imem = imem + 1 ))
  done

  2) Prepare the data structure for ensemble mean
  cp -f ./fv3sar_tile1_mem001_dynvar fv3sar_tile1_dynvar
  cp -f ./fv3sar_tile1_mem001_tracer fv3sar_tile1_tracer
  cp -f ./fv3sar_tile1_mem001_phyvar fv3sar_tile1_phyvar
  rm -f ./.hostfile_ensmean_*
  rm -f ./ensmean.output_*

  3) Run exectuable
  mpirun -n ${PROC} ${ENSMEAN_EXE} ./ fv3sar_tile1 ${ens_size} ${varname} ${ftail}

  For example, 
  mpirun -n 40 ${ENSMEAN_EXE} ./ fv3sar_tile1 40 u dynvar
 
  The above command provides the ensemble mean for u in fv3sar_tile1_dynvar

###########################################
  To calculate recentering after obtaining the ensemble mean

  1) Link control member
   ln -sf ${WORK_ROOT}/fv3prd_mem0000/ANA/fv_core.res.tile1_new.nc control_dynvar
   ln -sf ${WORK_ROOT}/fv3prd_mem0000/ANA/fv_tracer.res.tile1_new.nc control_tracer
   ln -sf ${WORK_ROOT}/fv3prd_mem0000/ANA/phy_data.nc control_phyvar

  2) Link all members to be recentered
   iens=1
   while [ ${iens} -le ${ens_size} ]
   do
     memberstring=`printf %03i ${iens}`
     ln -sf ./fv3sar_tile1_mem${memberstring}_dynvar  ./rec_fv3sar_tile1_mem${memberstring}_dynvar
     ln -sf ./fv3sar_tile1_mem${memberstring}_tracer  ./rec_fv3sar_tile1_mem${memberstring}_tracer
     ln -sf ./fv3sar_tile1_mem${memberstring}_phyvar  ./rec_fv3sar_tile1_mem${memberstring}_phyvar
     iens=`expr ${iens} + 1`
   done

  3) Run executable
    ${MPIRUN} -n ${PROC} ${RECENTER_EXE} ./ rec_fv3sar_tile1 fv3sar_tile1 control ${ens_size} ${varname} ${ftail}
