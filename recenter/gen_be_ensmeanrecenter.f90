program gen_be_ensmean_recenter
!
!---------------------------------------------------------------------- 
!  Purpose: Recenter ensemble members around EnVar control
!
!  2021-04 Yongming Wang and X. Wang - Changed from gen_be_ensmean
!                              - Enable parallel ensemble IO
!                                poc: xuguang.wang@ou.edu
!
!----------------------------------------------------------------------


   use netcdf
   implicit none

   integer, parameter    :: max_num_dims = 4          ! Maximum number of dimensions.
   integer, parameter    :: max_num_vars = 50         ! Maximum number of variables.
   integer, parameter    :: unit = 100                ! Unit number.
   integer, parameter    :: filename_len = 500

   character (len=filename_len)   :: directory                 ! General filename stub.
   character (len=filename_len)   :: filename                  ! General filename stub.
   character (len=filename_len)   :: filename_mean          ! name of the file of the new mean ( to be recentered at) 
   character (len=filename_len)   :: filename_meannew          ! name of the file of the new mean ( to be recentered at) 
   character (len=filename_len)   :: filetail
   character (len=filename_len)   :: input_file                ! Input file. 
   character (len=filename_len)   :: input_file_meanew                ! Input file. 
   character (len=10)    :: var                       ! Variable to search for.
   character (len=10)    :: charnanal                 ! ensmeble size
   character (len=3)     :: ce                        ! Member index -> character.

   integer               :: num_members               ! Ensemble size.
   integer               :: i              ! Loop counters.
   integer               :: length                    ! Filename length.
   integer               :: rcode                     ! NETCDF return code.
   integer               :: cdfid                     ! NETCDF file IDs.
   integer               :: cdfid_mean                ! NETCDF file IDs.
   integer               :: cdfid_meannew                ! NETCDF file IDs.
   integer               :: id_var                    ! NETCDF variable ID.
   integer               :: id_varmean                    ! NETCDF variable ID.
   integer               :: id_varmeannew                    ! NETCDF variable ID.
   integer               :: ivtype                    ! 4=integer, 5=real, 6=d.p.
   integer               :: ndims                     ! Number of field dimensions.
   integer               :: natts                     ! Number of field attributes.
   real                  :: rnanals                ! 1 / ensemble size.

   integer               :: dimids(1:10)              ! Array of dimension IDs.
   integer               :: dims(1:max_num_dims)      ! Array of dimensions.
   integer               :: istart(1:max_num_dims)    ! Array of dimension starts.
   integer               :: iend(1:max_num_dims)      ! Array of dimension ends.

   real (kind=4), allocatable     :: data_r(:,:,:)             ! Data array.
   real (kind=4), allocatable     :: data_r_mean(:,:,:)        ! Data array mean.
   real (kind=4), allocatable     :: data_r_meannew(:,:,:)        ! Data array mean.
   real (kind=4), allocatable     :: diff_data_r_mean(:,:,:)        ! Data array mean.
 
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

  if (mype == 0) print*, "Do recenter for FV3LAM"


  directory = './'
  filename = 'test'
  filename_mean = 'test'
  filename_meannew = 'test'
  
  var = 'U'
  charnanal = '40'

  ! Get user input from command line
  call getarg(1,directory)
  call getarg(2,filename)
  call getarg(3,filename_mean)
  call getarg(4,filename_meannew)
  call getarg(5,charnanal)
  call getarg(6,var)
  call getarg(7,filetail)

  read(charnanal,'(i2)') num_members
  rnanals = 1.0_8/num_members
  
  if ( mype == 0 )then
   write(6,'(a,a)')'   Directory = ', trim(directory)
   write(6,'(a,a)')'   filename = ', trim(filename)
   write(6,'(a,a)')'   filename_mean = ', trim(filename_mean)
   write(6,'(a,a)')'   filename_meannew = ', trim(filename_meannew)
   write(6,'(a,i4)')'   Number of ensemble members = ', num_members
   write(6,'(50a)')'   List of variables to average = ', var
  end if

  if ( npe < num_members ) then
     write(6,'(2(a,i4))')'***ERROR***  npe too small.  npe = ',npe,' < num_members = ',num_members
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  end if

! Create sub-communicator to handle number of cases (nanals)
  call mpi_comm_group(mpi_comm_world,orig_group,iret)

  allocate(new_group_members(num_members))
  do i=1,num_members
     new_group_members(i)=i-1
  end do

  call mpi_group_incl(orig_group,num_members,new_group_members,new_group,iret)
  call mpi_comm_create(mpi_comm_world,new_group,new_comm,iret)
  if ( iret /= 0 ) then
     write(6,'(a,i5)')'***ERROR*** after mpi_comm_create with iret = ',iret
     call mpi_abort(mpi_comm_world,101,iret)
  endif

  ! Process input files (one file per task)
  if ( mype1 <= num_members ) then

!  Open template ensemble mean with write access:
   input_file = trim(directory)//trim(filename_mean)//'_'//trim(filetail)
   length = len_trim(input_file)
   rcode = nf90_open(input_file(1:length), NF90_NOWRITE, cdfid_mean )
   if ( rcode /= 0) then
      write(UNIT=6,FMT='(A,A)') &
         ' Error opening netcdf file ', input_file(1:length)
   end if
!clt open the new mean files
   input_file_meanew = trim(directory)//trim(filename_meannew)//'_'//trim(filetail)
   length = len_trim(input_file_meanew)
   rcode = nf90_open(input_file_meanew(1:length), NF90_NOWRITE, cdfid_meannew )
   if ( rcode /= 0) then
      write(UNIT=6,FMT='(A,A)') &
         ' Error opening new mean netcdf file ', input_file(1:length)
   end if

   write(6,'(2a)')' Computing recenter for variable ', var
!  Open file:
!  Get variable ID:
   rcode = nf90_inq_varid ( cdfid_mean, var, id_varmean )

!  Check variable is in file:
   if ( rcode /= 0 ) then
        write(UNIT=6,FMT='(A,A)') &
        var, ' variable is not in input file'
   end if

!  Get dimension information for this field:
   dimids = 0
   dims = 0
   rcode = nf90_inquire_variable( cdfid_mean, id_varmean, var, ivtype, ndims, dimids, natts )
   do i = 1, ndims
      rcode = nf90_inquire_dimension( cdfid_mean, dimids(i), len=dims(i) )
   end do
   istart = 1
   iend = 1
   do i = 1, ndims-1
      iend(i) = dims(i)
   end do

!  Allocate and initialize data:
   allocate( data_r_mean(iend(1),iend(2),iend(3)))

!  get the old mean:
   call ncvgt( cdfid_mean, id_varmean, istart, iend, data_r_mean, rcode)

   if( var .eq. "ref_f3d" .or. var .eq. "sphum" .or. var .eq. "liq_wat" .or.   &
       var .eq. "rainwat" .or. var .eq. "ice_wat" .or. var .eq. "snowwat" .or. &
       var .eq. "graupel" )then
       where( data_r_mean .lt. 0.0 )
            data_r_mean = 0.0
       end where
   end if

!  Get variable ID:
   rcode = nf90_inq_varid ( cdfid_meannew, var, id_varmeannew )

!  Check variable is in file:
   if ( rcode /= 0 ) then
       write(UNIT=6,FMT='(A,A)') &
             var, ' variable is not in input file'
   end if

!  Get dimension information for this field:
   dimids = 0
   dims = 0
   rcode = nf90_inquire_variable( cdfid_meannew, id_varmeannew, var, ivtype, ndims, dimids, natts )
   do i = 1, ndims
      rcode = nf90_inquire_dimension( cdfid_mean, dimids(i), len=dims(i) )
   end do
   istart = 1
   iend = 1
   do i = 1, ndims-1
      iend(i) = dims(i)
   end do

!  Allocate and initialize data:
   allocate( data_r_meannew(iend(1),iend(2),iend(3)))

!  get the new mean:
   call ncvgt( cdfid_meannew, id_varmeannew, istart, iend, data_r_meannew, rcode)

   if( var .eq. "ref_f3d" .or. var .eq. "sphum" .or. var .eq. "liq_wat" .or.   &
       var .eq. "rainwat" .or. var .eq. "ice_wat" .or. var .eq. "snowwat" .or. & 
       var .eq. "graupel" )then
       where( data_r_meannew .lt. 0.0 )
             data_r_meannew = 0.0
       end where
   end if

   allocate( diff_data_r_mean(iend(1),iend(2),iend(3)))
   diff_data_r_mean=data_r_meannew-data_r_mean
   print*,maxval(diff_data_r_mean),minval(diff_data_r_mean)

   write(UNIT=ce,FMT='(i3.3)') mype1

!  Open file:
   input_file = trim(directory)//'/'//trim(filename)//'_mem'//trim(ce)//'_'//trim(filetail)
   print *, 'APM open ',input_file
   length = len_trim(input_file)
   rcode = nf90_open( input_file(1:length), NF90_WRITE, cdfid )

!  Get variable ID:
   rcode = nf90_inq_varid ( cdfid, var, id_var )

!  Check variable is in file:
   if ( rcode /= 0 ) then
        write(UNIT=6,FMT='(A,A)') &
              var, ' variable is not in input file'
   end if

!  Get dimension information for this field:
   dimids = 0
   dims = 0
   rcode = nf90_inquire_variable( cdfid, id_var, var, ivtype, ndims, dimids, natts )
   do i = 1, ndims
       rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
   end do
   istart = 1
   iend = 1
   do i = 1, ndims-1
      iend(i) = dims(i)
   end do

!  Allocate and initialize data:
   allocate( data_r(iend(1),iend(2),iend(3)))

!  Add diff to each member
   call ncvgt( cdfid, id_var, istart, iend, data_r, rcode)
   data_r=data_r+diff_data_r_mean

   if( var .eq. "ref_f3d" .or. var .eq. "sphum" .or. var .eq. "liq_wat" .or.   &
       var .eq. "rainwat" .or. var .eq. "ice_wat" .or. var .eq. "snowwat" .or. & 
       var .eq. "graupel" )then
        where( data_r .lt. 0.0 )
             data_r = 0.0
        end where
   end if

   call ncvpt( cdfid, id_var, istart, iend, data_r, rcode)

   rcode = nf90_close( cdfid )


   deallocate( data_r )
   deallocate( data_r_mean )
   deallocate( data_r_meannew )
   deallocate( diff_data_r_mean)


   rcode = nf90_close( cdfid_mean )
   rcode = nf90_close( cdfid_meannew )

  else
    write(6,'(a,i5)') 'No files to process for mpi task = ',mype
  end if
  call mpi_barrier(mpi_comm_world,iret)

  deallocate(new_group_members)

  if ( mype == 0 ) write(6,'(a)')'ensmean_done'

  call mpi_finalize(iret)

end program gen_be_ensmean_recenter

