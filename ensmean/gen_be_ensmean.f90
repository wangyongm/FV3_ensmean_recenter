program gen_be_ensmean
!
!---------------------------------------------------------------------- 
!  Purpose: Calculate ensemble mean file from input FV3LAM NETCDF input
!  ensemble members.
!
!  2021-04 Yongming Wang and X. Wang - Initial codes from WRFDA were changed for FV3LAM
!                              - Enable parallel ensemble IO
!                                poc: xuguang.wang@ou.edu
!
!----------------------------------------------------------------------

   use netcdf 
   implicit none

   integer, parameter    :: max_num_dims = 4          ! Maximum number of dimensions.
   integer, parameter    :: max_num_vars = 50         ! Maximum number of variables.
   integer, parameter    :: unit = 100                ! Unit number.
   integer, parameter    :: filename_len=500

   character (len=filename_len)   :: directory                 ! General filename stub.
   character (len=filename_len)   :: filename                  ! General filename stub.
   character (len=filename_len)   :: input_file                  ! General filename stub.
   character (len=filename_len)   :: filetail
   character (len=10)    :: varname                       ! Variable to search for.
   character (len=10)    :: charnanal                     ! ensmeble size
   character (len=3)     :: ce                        ! Member index -> character.

   integer               :: nanals                    ! Ensemble size.
   integer               :: nsize                     ! size of array
   integer               :: k,i                       ! Loop counters.
   integer               :: length                    ! Filename length.
   integer               :: rcode                     ! NETCDF return code.
   integer               :: cdfid                     ! NETCDF file IDs.
   integer               :: cdfid_mean                ! NETCDF file IDs.
   integer               :: id_var                    ! NETCDF variable ID.
   integer               :: ivtype                    ! 4=integer, 5=real, 6=d.p.
   integer               :: ndims                     ! Number of field dimensions.
   integer               :: natts                     ! Number of field attributes.
   real                  :: rnanals                   ! 1 / ensemble size.

   integer               :: dimids(1:10)              ! Array of dimension IDs.
   integer               :: dims(1:max_num_dims)      ! Array of dimensions.
   integer               :: istart(1:max_num_dims)    ! Array of dimension starts.
   integer               :: iend(1:max_num_dims)      ! Array of dimension ends.

   real (kind=4), allocatable     :: data_r(:,:,:)             ! Data array.
   real (kind=4), allocatable     :: data_r_mean(:,:,:)        ! Data array mean.

   ! === variables for mpi 
   integer                :: iret, mype, npe, mype1, orig_group, new_group, new_comm
   integer,allocatable,dimension(:) :: new_group_members
   !
  
  include 'mpif.h'

! Initialize mpi, mype is process number, npe is total number of processes.
  call mpi_init(iret)
  call mpi_comm_rank(mpi_comm_world,mype,iret)
  call mpi_comm_size(mpi_comm_world,npe,iret)

  mype1 = mype + 1

  if (mype == 0) print*, "Calculate ensemble mean for FV3LAM"

  ! Get user input from command line
  call getarg(1,directory)
  call getarg(2,filename)
  call getarg(3,charnanal)
  call getarg(4,varname)
  call getarg(5,filetail)
  read(charnanal,'(i2)') nanals

  rnanals = 1.0_8/nanals
  filename = trim(adjustl(directory)) // trim(adjustl(filename))

  if ( mype == 0 ) then
     write(6,'(a)')  'Command line input'
     write(6,'(a,a)')' directory    = ',trim(directory)
     write(6,'(a,a)')' filename    = ',trim(filename)//'_'//trim(filetail)
     write(6,'(a,a)')' nanals      = ',trim(charnanal)
     write(6,'(a,a)')' Processing ...', trim(varname)
     write(6,'(a)')  ' '
  endif

  if ( npe < nanals ) then
     write(6,'(2(a,i4))')'***ERROR***  npe too small.  npe = ',npe,' < nanals = ',nanals
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  end if

! Create sub-communicator to handle number of cases (nanals)
  call mpi_comm_group(mpi_comm_world,orig_group,iret)

  allocate(new_group_members(nanals))
  do k=1,nanals
     new_group_members(k)=k-1
  end do

  call mpi_group_incl(orig_group,nanals,new_group_members,new_group,iret)
  call mpi_comm_create(mpi_comm_world,new_group,new_comm,iret)
  if ( iret /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',iret
     call mpi_abort(mpi_comm_world,101,iret)
  endif

  ! Process input files (one file per task)
  if ( mype1 <= nanals ) then

!  Open template ensemble mean with write access:

   if ( mype == 0  )then
     input_file = trim(directory)//'/'//trim(filename)//'_'//trim(filetail)
     print*,'mean file',input_file
     length = len_trim(input_file)
     rcode = nf90_open(input_file(1:length), NF90_WRITE, cdfid_mean )
     if ( rcode /= 0) then
        write(6,FMT='(A,A)') &
          ' Error opening netcdf file ', input_file(1:length)
     end if
   end if


   write(6,'(2a)')' Computing ensemble mean for variable ', varname

   write(UNIT=ce,FMT='(i3.3)')mype1

!  Open file:
   input_file = trim(directory)//'/'//trim(filename)//'_mem'//trim(ce)//'_'//trim(filetail)
   print *, 'APM open ',input_file
   length = len_trim(input_file)
   rcode = nf90_open( input_file(1:length), NF90_NOWRITE, cdfid )

   rcode = nf90_inq_varid ( cdfid, trim(varname), id_var )

!  Check variable is in file:
   if ( rcode /= 0 ) then
      write(6,FMT='(A,A)') &
            varname, ' variable is not in input file'
   end if

!  Get dimension information for this field:
   dimids = 0
   dims = 0
   rcode = nf90_inquire_variable( cdfid, id_var, varname, ivtype, ndims, dimids, natts )
   do i = 1, ndims
      rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
   end do
   istart = 1
   iend = 1
   do i = 1, ndims-1
      iend(i) = dims(i)
   end do

print*,iend,cdfid,ndims

!  Allocate and initialize data:
   !if ( ivtype == 5 ) then
      allocate( data_r(iend(1),iend(2),iend(3)))
      allocate( data_r_mean(iend(1),iend(2),iend(3)))
      data_r_mean = 0.0
   !else
   !   write(6,FMT='(A,A)') &
   !         varname, ' variable is not real type'
   !end if

!  Calculate accumulating mean and variance:
   call ncvgt( cdfid, id_var, istart, iend, data_r, rcode)
   if( varname == "ref_f3d" ) then
       where( data_r .lt.  0.0 )
            data_r = 0.0
       end where
   end if
   nsize=iend(1)*iend(2)*iend(3)

   call mpi_allreduce(data_r,data_r_mean,nsize,mpi_real,mpi_sum,new_comm,iret)

   print*,maxval(data_r_mean),minval(data_r_mean)
   data_r_mean = data_r_mean * rnanals
   print*,maxval(data_r_mean),minval(data_r_mean)

   rcode = nf90_close( cdfid )

   if ( mype == 0 )then

!   Open file: dif dim/varids for dif members so use member1
    input_file = trim(directory)//'/'//trim(filename)//'_mem001'//'_'//trim(filetail)
    length = len_trim(input_file)
    rcode = nf90_open( input_file(1:length), NF90_NOWRITE, cdfid )
!   Get dimension information for this field:
    dimids = 0
    dims = 0
    rcode = nf90_inquire_variable( cdfid, id_var, varname, ivtype, ndims, dimids, natts )
    do i = 1, ndims
       rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
    end do
    istart = 1
    iend = 1
    do i = 1, ndims-1
       iend(i) = dims(i)
    end do

    call ncvpt( cdfid_mean, id_var, istart, iend, data_r_mean, rcode)

    rcode = nf90_close( cdfid_mean )
    rcode = nf90_close( cdfid )
   end if
   deallocate( data_r )
   deallocate( data_r_mean )
  else
    write(6,'(a,i5)') 'No files to process for mpi task = ',mype
  end if
  call mpi_barrier(mpi_comm_world,iret)

  deallocate(new_group_members)

  if ( mype == 0 ) write(6,'(a)')'ensmean_done'

  call mpi_finalize(iret)

end program gen_be_ensmean

