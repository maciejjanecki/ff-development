!===============================================================================
! SVN $Id: seq_io_mod.F90 24486 2010-08-25 21:28:04Z mvr $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/branch_tags/cesm1_0_rel_tags/cesm1_0_rel03_drvseq3_1_35/shr/seq_io_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_io_mod -- reads and writes driver files
!
! !DESCRIPTION:
!    Writes attribute vectors to netcdf
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig first version
!    2007-Dec-06 - T. Craig update and improve
!    2008-Feb-16 - J. Edwards convert to PIO
!
! !INTERFACE: ------------------------------------------------------------------

module seq_io_mod

  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_sys_mod       ! system calls
  use seq_comm_mct, only : logunit, seq_comm_setptrs, CPLID
  use seq_flds_mod, only : seq_flds_lookup
  use mct_mod           ! mct wrappers
  use pio

  implicit none
  private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

  public seq_io_init
  public seq_io_wopen
  public seq_io_close
  public seq_io_redef
  public seq_io_enddef
  public seq_io_date2yyyymmdd
  public seq_io_sec2hms
  public seq_io_read
  public seq_io_write

! !PUBLIC DATA MEMBERS

  ! none

!EOP

  interface seq_io_read
     module procedure seq_io_read_av
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
  end interface
  interface seq_io_write
     module procedure seq_io_write_av
     module procedure seq_io_write_int
     module procedure seq_io_write_int1d
     module procedure seq_io_write_r8
     module procedure seq_io_write_r81d
     module procedure seq_io_write_char
     module procedure seq_io_write_time
  end interface

!-------------------------------------------------------------------------------
! Local data
!-------------------------------------------------------------------------------

   character(*),parameter :: prefix = "seq_io_"
   character(CL)          :: wfilename = ''
   real(r8)    ,parameter :: fillvalue = SHR_CONST_SPVAL

   character(*),parameter :: modName = "(seq_io_mod) "
   integer(in) ,parameter :: debug = 0 ! internal debug level

   integer(IN)          ,save :: cpl_io_type
   type(file_desc_t)    ,save :: cpl_io_file
   type(iosystem_desc_t),save :: cpl_io_subsystem

   character(*),parameter :: version ='cpl7v10'
   character(*),parameter :: version0='cpl7v00'

   character(CL) :: charvar   ! buffer for string read/write

!=================================================================================
contains
!=================================================================================

!=================================================================================
!BOP =============================================================================
!
! !IROUTINE: seq_io_init - initialize io for coupler
!
! !DESCRIPTION:
!    Read the pio_inparm namelist and initialize the pio subsystem
!
! !REVISION HISTORY:
!    2009-Sep-30 - B. Kauffman - consolidation, clean up
!    2009-Feb-17 - J. Edwards - initial version
!
! !INTERFACE: --------------------------------------------------------------------

subroutine seq_io_init(nlfilename)

   use shr_string_mod, only : shr_string_toupper
   use shr_file_mod,   only : shr_file_getunit, shr_file_freeunit
   use shr_mpi_mod,    only : shr_mpi_bcast
   use pio, only : pio_iotype_netcdf, pio_iotype_pnetcdf, pio_iotype_netcdf4c, pio_iotype_netcdf4p

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*), intent(in) :: nlfilename

! LOCAL

   integer(IN)   :: iam, mpicom, npes, ierr, unitn
   logical       :: iamroot
   character(CS) :: cpl_io_typename
   integer(IN)   :: cpl_io_stride
   integer(IN)   :: cpl_io_numtasks
   integer(IN)   :: cpl_io_root
   integer(IN),parameter :: cpl_io_root_default = 0

   namelist /pio_inparm/ cpl_io_stride, cpl_io_root, cpl_io_numtasks, cpl_io_typename

   character(*),parameter :: subName =   '(seq_io_init) '
   character(*),parameter :: F00     = "('(seq_io_init) ',4a)"
   character(*),parameter :: F01     = "('(seq_io_init) ',a,i6)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID, iam=iam, mpicom=mpicom,npes=npes, iamroot=iamroot)

    !--------------------------------------------------------------------------
    ! read io nml parameters 
    !--------------------------------------------------------------------------
    cpl_io_stride   = -1 ! set based on cpl_io_numtasks value when initialized < 0
    cpl_io_numtasks = -1 ! set based on cpl_io_stride   value when initialized < 0
    cpl_io_root     =  cpl_io_root_default
    cpl_io_typename = 'netcdf'

    if (iamroot) then
       if (debug > 0) then
          write(logunit,F00) 'pio init parameters: before nml read'
          write(logunit,F01) '   cpl_io_stride   = ',cpl_io_stride
          write(logunit,F01) '   cpl_io_root     = ',cpl_io_root
          write(logunit,F00) '   cpl_io_typename = ',cpl_io_typename
          write(logunit,F01) '   cpl_io_numtasks = ',cpl_io_numtasks
       end if

       unitn=shr_file_getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ierr = 1
       do while( ierr /= 0 )
          read(unitn,nml=pio_inparm,iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( subname//':: namelist read returns an'// &
                                 ' end of file or end of record condition' )
          end if
       end do
       close(unitn)
       call shr_file_freeUnit( unitn )

       if (debug > 0) then
          write(logunit,F00) 'pio init parameters: after nml read'
          write(logunit,F01) '   cpl_io_stride   = ',cpl_io_stride
          write(logunit,F01) '   cpl_io_root     = ',cpl_io_root
          write(logunit,F00) '   cpl_io_typename = ',cpl_io_typename
          write(logunit,F01) '   cpl_io_numtasks = ',cpl_io_numtasks
       end if

       if      ( shr_string_toupper(cpl_io_typename) .eq. 'NETCDF' ) then
          cpl_io_type = pio_iotype_netcdf
       else if ( shr_string_toupper(cpl_io_typename) .eq. 'PNETCDF') then
          cpl_io_type = pio_iotype_pnetcdf
       else if ( shr_string_toupper(cpl_io_typename) .eq. 'NETCDF4P') then
          cpl_io_type = pio_iotype_netcdf4p
       else if ( shr_string_toupper(cpl_io_typename) .eq. 'NETCDF4C') then
          cpl_io_type = pio_iotype_netcdf4c
       else
          write(logunit,*) subName,'Bad io_type argument - using iotype_netcdf'
          cpl_io_type=pio_iotype_netcdf
       end if
    end if
    call shr_mpi_bcast(cpl_io_type , mpicom)
    call shr_mpi_bcast(cpl_io_stride  , mpicom)
    call shr_mpi_bcast(cpl_io_root    , mpicom)
    call shr_mpi_bcast(cpl_io_numtasks, mpicom)

    !--------------------------------------------------------------------------
    ! check/set/correct io pio parameters
    !--------------------------------------------------------------------------


    if (cpl_io_stride>0.and.cpl_io_numtasks<0) then
       cpl_io_numtasks = npes/cpl_io_stride
    else if(cpl_io_numtasks>0 .and. cpl_io_stride<0) then
       cpl_io_stride = npes/cpl_io_numtasks
    else if(cpl_io_numtasks<0 .and. cpl_io_stride<0) then
       cpl_io_stride = 4
       cpl_io_numtasks = npes/cpl_io_stride
       cpl_io_numtasks = max(1, cpl_io_numtasks)
    end if

    if (cpl_io_root<0) then
       cpl_io_root = cpl_io_root_default
    endif
    cpl_io_root = min(cpl_io_root,npes-1)

    if (cpl_io_root + (cpl_io_stride)*(cpl_io_numtasks-1) >= npes .or. &
       cpl_io_stride<=0 .or. cpl_io_numtasks<=0 .or. cpl_io_root < 0 .or. &
       cpl_io_root > npes-1) then
       write(logunit,*) subName,'cpl_io_stride, iotasks or root out of bounds - resetting to defaults ', &
		cpl_io_stride, cpl_io_numtasks, cpl_io_root
       cpl_io_stride = max(1,npes/4)
       cpl_io_numtasks = npes/cpl_io_stride
       cpl_io_root = min(1,npes-1)
    end if

    !--------------------------------------------------------------------------
    ! init pio library
    !--------------------------------------------------------------------------
    if (iamroot) then
       write(logunit,F00) 'pio init parameters: '
       write(logunit,F01) '   cpl_io_stride   = ',cpl_io_stride
       write(logunit,F01) '   cpl_io_root     = ',cpl_io_root
       write(logunit,F00) '   cpl_io_typename = ',cpl_io_typename
       write(logunit,F01) '   cpl_io_numtasks = ',cpl_io_numtasks
    end if

!   call pio_init(iam, mpicom, cpl_io_numtasks, cpl_io_root, cpl_io_stride, PIO_REARR_BOX, cpl_io_subsystem)
    call pio_init(iam, mpicom, cpl_io_numtasks, 0, cpl_io_stride, PIO_REARR_BOX, &
                  cpl_io_subsystem, base=cpl_io_root)
    
    cpl_io_file%fh=-1

  end subroutine seq_io_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_wopen - open netcdf file
!
! !DESCRIPTION:
!    open netcdf file
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig - initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_io_wopen(filename,clobber,cdf64)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(*),intent(in) :: filename
    logical,optional,intent(in):: clobber
    logical,optional,intent(in):: cdf64

    !EOP

    logical :: exists
    logical :: lclobber
    logical :: lcdf64
    integer :: iam,mpicom
    integer :: rcode
    integer :: nmode
    character(CL)  :: lversion
    character(*),parameter :: subName = '(seq_io_wopen) '
    
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    lversion=trim(version0)

    lclobber = .false.
    if (present(clobber)) lclobber=clobber

    lcdf64 = .false.
    if (present(cdf64)) lcdf64=cdf64

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    if (cpl_io_file%fh<0) then
       ! filename not open
       if (iam==0) inquire(file=trim(filename),exist=exists)
       call shr_mpi_bcast(exists,mpicom,'seq_io_wopen exists')
       if (exists) then
          if (lclobber) then
             nmode = pio_clobber
             if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
             rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_io_type, trim(filename), nmode)
             if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
             rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
          else
             rcode = pio_openfile(cpl_io_subsystem, cpl_io_file, cpl_io_type, trim(filename), pio_write)
             if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
             call pio_seterrorhandling(cpl_io_file,PIO_BCAST_ERROR)
             rcode = pio_get_att(cpl_io_file,pio_global,"file_version",lversion)
             call pio_seterrorhandling(cpl_io_file,PIO_INTERNAL_ERROR)
             if (trim(lversion) /= trim(version)) then
                rcode = pio_redef(cpl_io_file)
                rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
                rcode = pio_enddef(cpl_io_file)
             endif
          endif
       else
          nmode = pio_noclobber
          if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
          rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_io_type, trim(filename), nmode)
          if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
          rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
       endif
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename is open, better match open filename
       if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
       call shr_sys_abort()
    else
       ! filename is already open, just return
    endif

end subroutine seq_io_wopen

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_close - close netcdf file
!
! !DESCRIPTION:
!    close netcdf file
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig - initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_io_close(filename)

    use pio, only : pio_closefile

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(*),intent(in) :: filename

    !EOP

    integer :: iam
    integer :: rcode
    character(*),parameter :: subName = '(seq_io_close) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam)


    if (cpl_io_file%fh<0) then
       ! filename not open, just return
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename matches, close it
       call pio_closefile(cpl_io_file)
       cpl_io_file%fh=-1
    else
       ! different filename is open, abort
       if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
       call shr_sys_abort()
    endif

    wfilename = ''

end subroutine seq_io_close

!===============================================================================

subroutine seq_io_redef(filename)
    character(len=*), intent(in) :: filename
    integer :: rcode

    rcode = pio_redef(cpl_io_file)
end subroutine seq_io_redef

!===============================================================================

subroutine seq_io_enddef(filename)
    character(len=*), intent(in) :: filename
    integer :: rcode

    rcode = pio_enddef(cpl_io_file)
end subroutine seq_io_enddef

!===============================================================================

character(len=10) function seq_io_date2yyyymmdd (date)

! Input arguments

   integer, intent(in) :: date

! Local workspace

   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

!-------------------------------------------------------------------------------

   if (date < 0) then
      call shr_sys_abort ('seq_io_date2yyyymmdd: negative date not allowed')
   end if

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   write(seq_io_date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)

end function seq_io_date2yyyymmdd

!===============================================================================

character(len=8) function seq_io_sec2hms (seconds)

! Input arguments

   integer, intent(in) :: seconds

! Local workspace

   integer :: hours     ! hours of hh:mm:ss
   integer :: minutes   ! minutes of hh:mm:ss
   integer :: secs      ! seconds of hh:mm:ss

!-------------------------------------------------------------------------------

   if (seconds < 0 .or. seconds > 86400) then
      write(logunit,*)'seq_io_sec2hms: bad input seconds:', seconds
      call shr_sys_abort()
   end if

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   if (minutes < 0 .or. minutes > 60) then
      write(logunit,*)'seq_io_sec2hms: bad minutes = ',minutes
      call shr_sys_abort()
   end if

   if (secs < 0 .or. secs > 60) then
      write(logunit,*)'seq_io_sec2hms: bad secs = ',secs
      call shr_sys_abort()
   end if

   write(seq_io_sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)

end function seq_io_sec2hms

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_av - write AV to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_av(filename,gsmap,AV,dname,whead,wdata,nx,ny,nt,fillval,pre,tavg,&
	                     use_float)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(in) :: AV       ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer(in),optional,intent(in) :: nx   ! 2d grid size if available
    integer(in),optional,intent(in) :: ny   ! 2d grid size if available
    integer(in),optional,intent(in) :: nt   ! time sample
    real(r8),optional,intent(in) :: fillval ! fill value
    character(len=*),optional,intent(in) :: pre      ! prefix to variable name
    logical,optional,intent(in) :: tavg     ! is this a tavg
    logical,optional,intent(in) :: use_float ! write output as float rather than double

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: nf,ns,ng
    integer(in) :: i,j,k,n
    integer(in),target  :: dimid2(2)
    integer(in),target  :: dimid3(3)
    integer(in),pointer :: dimid(:)
    type(var_desc_t) :: varid
    type(io_desc_t) :: iodesc
    integer(kind=PIO_OffSet) :: frame
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    character(CL)    :: name1       ! var name
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    character(CL)    :: lpre        ! local prefix
    logical :: exists
    logical :: lwhead, lwdata
    integer(in) :: lnx,lny
    real(r8) :: lfillvalue
    type(mct_aVect) :: AVroot
    real(r8),pointer :: fld1(:,:)  ! needed to convert AVroot ng rAttr to 2d nx,ny
    character(*),parameter :: subName = '(seq_io_write_av) '

    integer :: lbnum
    integer, pointer :: Dof(:)

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lfillvalue = fillvalue
    if (present(fillval)) then
       lfillvalue = fillval
    endif

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    call seq_comm_setptrs(CPLID,iam=iam)	

    ng = mct_gsmap_gsize(gsmap)
    lnx = ng
    lny = 1
	
    nf = mct_aVect_nRattr(AV)
    if (nf < 1) then
       write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
       call shr_sys_abort()
    endif

    if (present(nx)) then
       if (nx /= 0) lnx = nx
    endif
    if (present(ny)) then
       if (ny /= 0) lny = ny
    endif
    if (lnx*lny /= ng) then
       if(iam==0) write(logunit,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny,trim(dname)
       call shr_sys_abort()
    endif

    if (lwhead) then
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))

       if (present(nt)) then
          dimid3(1:2) = dimid2
          rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
          dimid => dimid3
       else
          dimid => dimid2
       endif

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AV)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
! "v0"    name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
          name1 = trim(lpre)//'_'//trim(itemc)
          call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
	  if (present(use_float)) then 
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
          else
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
          end if
          rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
          rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
          rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
          if (present(tavg)) then
             if (tavg) then
                rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
             endif
          endif
       enddo
       if (lwdata) call seq_io_enddef(filename)
    end if

    if (lwdata) then
       call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
       call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       deallocate(dof)

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AV)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
! "v0"    name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
          name1 = trim(lpre)//'_'//trim(itemc)
          rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
          if (present(nt)) then
             frame = nt
          else
             frame = 1
          endif
          call pio_setframe(varid,frame)
          call pio_write_darray(cpl_io_file, varid, iodesc, av%rattr(k,:), rcode, fillval=lfillvalue)
       enddo

       call pio_freedecomp(cpl_io_file, iodesc)

    end if
  end subroutine seq_io_write_av

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int - write scalar integer to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar integer to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_int(filename,idata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: exists
    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_int) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
!       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
!       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,idata)

       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
    endif

  end subroutine seq_io_write_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int1d - write 1d integer array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d integer array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_int1d(filename,idata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_int1d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(idata)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,idata)
    endif

       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata

  end subroutine seq_io_write_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r8 - write scalar double to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar double to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_r8(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: exists
    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_r8) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
!       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
!       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,rdata)
    endif


  end subroutine seq_io_write_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r81d - write 1d double array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d double array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_r81d(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_r81d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(rdata)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,rdata)

       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
    endif

  end subroutine seq_io_write_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_char - write char string to netcdf file
  !
  ! !DESCRIPTION:
  !    Write char string to netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_char(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    character(len=*),intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_char) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = len(charvar)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_len',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_CHAR,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       charvar = ''
       charvar = trim(rdata)
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,charvar)
    endif

  end subroutine seq_io_write_char

  !===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_io_write_time - write time variable to netcdf file
!
! !DESCRIPTION:
!    Write time variable to netcdf file
!
! !REVISION HISTORY:
!    2009-Feb-11 - M. Vertenstein - initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_io_write_time(filename,time_units,time_cal,time_val,nt,whead,wdata,tbnds)

! !INPUT/OUTPUT PARAMETERS:
   implicit none
   character(len=*),intent(in) :: filename      ! file
   character(len=*),intent(in) :: time_units    ! units of time
   character(len=*),intent(in) :: time_cal      ! calendar type
   real(r8)        ,intent(in) :: time_val      ! data to be written
   integer(in),optional,intent(in) :: nt
   logical,optional,intent(in) :: whead         ! write header
   logical,optional,intent(in) :: wdata         ! write data
   real(r8),optional,intent(in) :: tbnds(2)     ! time bounds

!EOP

   integer(in) :: rcode
   integer(in) :: iam
   integer(in) :: dimid(1)
   integer(in) :: dimid2(2)
   type(var_desc_t) :: varid
   integer(in) :: lnx
   logical :: exists
   logical :: lwhead, lwdata
   integer :: start(4),count(4)
   real(r8) :: time_val_1d(1)
   character(*),parameter :: subName = '(seq_io_write_time) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lwhead = .true.
   lwdata = .true.
   if (present(whead)) lwhead = whead
   if (present(wdata)) lwdata = wdata

   if (.not.lwhead .and. .not.lwdata) then
      ! should we write a warning?
      return
   endif

   call seq_comm_setptrs(CPLID,iam=iam)

   if (lwhead) then
      rcode = pio_def_dim(cpl_io_file,'time',PIO_UNLIMITED,dimid(1))
      rcode = pio_def_var(cpl_io_file,'time',PIO_DOUBLE,dimid,varid)
      rcode = pio_put_att(cpl_io_file,varid,'units',trim(time_units))
      if (trim(time_cal) == 'NO_LEAP') then
         rcode = pio_put_att(cpl_io_file,varid,'calendar','noleap')
      else if (trim(time_cal) == 'GREGORIAN') then
         rcode = pio_put_att(cpl_io_file,varid,'calendar','365_day')
      else
         rcode = pio_put_att(cpl_io_file,varid,'calendar','time_cal')
      end if
      if (present(tbnds)) then
         rcode = pio_put_att(cpl_io_file,varid,'bounds','time_bnds')
         dimid2(2)=dimid(1)
         rcode = pio_def_dim(cpl_io_file,'ntb',2,dimid2(1))
         rcode = pio_def_var(cpl_io_file,'time_bnds',PIO_DOUBLE,dimid2,varid)
      endif
      if (lwdata) call seq_io_enddef(filename)
   endif

   if (lwdata) then
      start = 1
      count = 1
      if (present(nt)) then
         start(1) = nt
      endif
      time_val_1d(1) = time_val
      rcode = pio_inq_varid(cpl_io_file,'time',varid)
      rcode = pio_put_var(cpl_io_file,varid,start,count,time_val_1d)
      if (present(tbnds)) then
         rcode = pio_inq_varid(cpl_io_file,'time_bnds',varid)
         start = 1
         count = 1
         if (present(nt)) then
            start(2) = nt
         endif
         count(1) = 2
         rcode = pio_put_var(cpl_io_file,varid,start,count,tbnds)
      endif

      !      write(logunit,*) subname,' wrote time ',lwhead,lwdata
   endif

 end subroutine seq_io_write_time

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_av - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_av(filename,gsmap,AV,dname,pre)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(inout):: AV     ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    character(len=*),intent(in),optional :: pre      ! prefix name

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    integer(in) :: nf,ns,ng
    integer(in) :: i,j,k,n, ndims
    type(file_desc_t) :: pioid
    integer(in) :: dimid(2)
    type(var_desc_t) :: varid
    integer(in) :: lnx,lny
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    logical :: exists
    type(io_desc_t) :: iodesc
    integer(in), pointer :: dof(:)
    character(CL)  :: lversion
    character(CL)  :: name1
    character(CL)  :: lpre
    character(*),parameter :: subName = '(seq_io_read_av) '
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lversion = trim(version0)

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
    call mct_gsmap_OrderedPoints(gsmap, iam, Dof)

    ns = mct_aVect_lsize(AV)
    nf = mct_aVect_nRattr(AV)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_av exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_io_type, trim(filename),pio_nowrite)
       if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    do k = 1,nf
       call mct_aVect_getRList(mstring,k,AV)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       if (trim(lversion) == trim(version)) then
          name1 = trim(lpre)//'_'//trim(itemc)
       else
          name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
       endif
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid,trim(name1),varid)
       if (rcode == pio_noerr) then
          if (k==1) then
             rcode = pio_inq_varndims(pioid, varid, ndims)
             rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
             rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
             if (ndims==2) then
                rcode = pio_inq_dimlen(pioid, dimid(2), lny)
             else
                lny = 1
             end if
             ng = lnx * lny
             if (ng /= mct_gsmap_gsize(gsmap)) then
                if (iam==0) write(logunit,*) subname,' ERROR: dimensions do not match',&
                     lnx,lny,mct_gsmap_gsize(gsmap)
                call shr_sys_abort()
             end if
             call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
             deallocate(dof)
          end if
          call pio_read_darray(pioid,varid,iodesc, av%rattr(k,:), rcode)
       else
          write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
          write(logunit,*)'for backwards compatibility will set it to 0'
          av%rattr(k,:) = 0.0_r8
       end if
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
          
       !--- zero out fill value, this is somewhat arbitrary
       do n = 1,ns
          if (AV%rAttr(k,n) == fillvalue) then
              AV%rAttr(k,n) = 0.0_r8
          endif
       enddo
    enddo

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine seq_io_read_av

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int - read scalar integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar integer from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int(filename,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(inout):: idata  ! integer data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer :: i1d(1)
    character(*),parameter :: subName = '(seq_io_read_int) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_int1d(filename,i1d,dname)
    idata = i1d(1)

  end subroutine seq_io_read_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int1d - read 1d integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d integer array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int1d(filename,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(inout):: idata(:)  ! integer data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_t) :: pioid 
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_int1d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_int1d exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_io_type, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)

    call pio_closefile(pioid)

    !      write(logunit,*) subname,' read int ',trim(dname)


  end subroutine seq_io_read_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r8 - read scalar double from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar double from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r8(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(inout):: rdata  ! real data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(seq_io_read_r8) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_r81d(filename,r1d,dname)
    rdata = r1d(1)

  end subroutine seq_io_read_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r81d - read 1d double array from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d double array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r81d(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(inout):: rdata(:)  ! real data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_T) :: pioid 
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_r81d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_r81d exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_io_type, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

    call pio_closefile(pioid)

    !      write(logunit,*) subname,' read int ',trim(dname)

  end subroutine seq_io_read_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_char - read char string from netcdf file
  !
  ! !DESCRIPTION:
  !    Read char string from netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_char(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    character(len=*),intent(inout):: rdata  ! character data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_T) :: pioid 
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_char) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_char exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_io_type, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,charvar)
    rdata = trim(charvar)

    call pio_closefile(pioid)

  end subroutine seq_io_read_char

  !===============================================================================
!===============================================================================
end module seq_io_mod
