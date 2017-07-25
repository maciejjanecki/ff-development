
module mo_flbc
  !---------------------------------------------------------------
  ! 	... lower boundary module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use m_types,      only : time_ramp
  use chem_mods,    only : pcnstm1=>gas_pcnst
  use spmd_utils,   only : masterproc,iam
  use abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use ppgrid,       only : pcols, begchunk, endchunk, pver
  use time_manager, only : get_curr_date, get_curr_calday
  use time_utils,   only : flt_date
  use mo_tracname,  only : tracnam=>solsym
  use wrap_nf
  use cam_logfile,  only : iulog
#ifdef SPMD
  use mpishorthand, only : mpicom, mpiint
#endif

  implicit none

  type :: flbc
     integer           :: spc_ndx
     real(r8), pointer     :: vmr(:,:,:)
     character(len=16)  :: species
  end type flbc

  private
  public  :: flbc_inti, flbc_set, flbc_chk, has_flbc

  save

  integer, parameter :: time_span = 1

  integer :: ntimes
  integer :: flbc_cnt
  integer :: gndx
  integer :: tim_ndx(2)
  integer :: jlim(2)
  integer, allocatable  :: dates(:)
  real(r8), allocatable     :: times(:)
  logical :: has_flbc(pcnstm1)
  character(len=256) :: filename, lpath, mspath

  type(time_ramp)             :: flbc_timing
  type(flbc), allocatable     :: flbcs(:)
  integer ::  ncdate, ncsec

contains

  subroutine flbc_inti( flbc_file, flbc_list, flbc_timing_in )
    !-----------------------------------------------------------------------
    ! 	... initialize the fixed lower bndy cond
    !-----------------------------------------------------------------------

    use chem_mods,     only : adv_mass
    use mo_constants,  only : d2r, pi, rearth
    use string_utils,  only : to_upper
    use mo_chem_utls,  only : get_spc_ndx
    use constituents,  only : pcnst

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: flbc_file
    character(len=*), intent(in) :: flbc_list(:)
    type(time_ramp),  intent(in) :: flbc_timing_in

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n                           ! Indices
    integer :: t1, t2
    integer :: ncid
    integer :: dimid
    integer :: varid
    integer :: yr, mon, day, wrk_date, wrk_sec
    real(r8)    :: seq
    real(r8)    :: wrk_time
    character(len=16)  :: species
    character(len=16)  :: spc_name
    character(len=8)  :: time_type

#include <netcdf.inc>

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    !-----------------------------------------------------------------------
    ! 	... check timing
    !-----------------------------------------------------------------------
    flbc_timing = flbc_timing_in
    time_type = to_upper(flbc_timing%type)
    flbc_timing%type = time_type
    if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' &
         .and. time_type /= 'FIXED' ) then
       write(iulog,*) 'flbc_inti: time type ',trim(time_type),' is not SERIAL,CYCLICAL, or FIXED'
       call endrun
    end if

    wrk_sec  = ncsec
    if( time_type == 'SERIAL' ) then
       wrk_date = ncdate + flbc_timing%yr_offset*10000
    else if( time_type == 'CYCLICAL' ) then

    	! If this is a leap-day, we have to avoid asking for a non-leap-year
    	! on a cyclical dataset. When this happens, just use Feb 28 instead
    	if (( mon .eq. 2 ) .and. ( day.eq.29 )) then
	   ncdate = yr*10000 + mon*100 + (day-1)
           write(iulog,*)'WARNING: flbc_inti using Feb 28 instead of Feb 29 for cyclical dataset'
        endif 	
       wrk_date = (flbc_timing%date/10000)*10000 + mod(ncdate,10000)
    else
       wrk_date = flbc_timing%date
       wrk_sec  = 0
    end if
    wrk_time = flt_date( wrk_date, wrk_sec )
    if (masterproc) write(iulog,*) 'flbc_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time

    !-----------------------------------------------------------------------
    ! 	... species with fixed lbc ?
    !-----------------------------------------------------------------------
    has_flbc(:) = .false.
    flbc_cnt = 0

    do m = 1,pcnst

       if ( len_trim(flbc_list(m))==0 ) exit

       n = get_spc_ndx(flbc_list(m))

       if (n > 0) then
          has_flbc(n) = .true.
          flbc_cnt = flbc_cnt + 1
       else
          write(iulog,*) 'flbc_init: '//flbc_list(m)//' is not included in species set'
          call endrun('flbc_init: invalid fixed lower boundary species')
       endif

    enddo

    if( flbc_cnt == 0 ) then
       return
    end if
    !-----------------------------------------------------------------------
    ! 	... allocate type array
    !-----------------------------------------------------------------------
    allocate( flbcs(flbc_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate flbc array; error = ',astat
       call endrun
    end if

    n = 0
    do m = 1,pcnstm1
       if( has_flbc(m) ) then
          n = n + 1

          flbcs(n)%spc_ndx = m

          if( trim(tracnam(m)) == 'CFC11' ) then
             flbcs(n)%species = 'CFCL3'
          else if( trim(tracnam(m)) == 'CFC12' ) then
             flbcs(n)%species = 'CF2CL2'
          else
             flbcs(n)%species = trim( tracnam(m) )
          end if

       end if
    end do

    if (masterproc) then

       write(iulog,*) ' '
       if( flbc_cnt > 0 ) then
          write(iulog,*) 'Species with specified lower boundary values'
          do m = 1,pcnstm1
             if( has_flbc(m) ) then
                write(iulog,*) trim(tracnam(m))
             end if
          end do
       else
          write(iulog,*) 'There are no species with specified lower boundary values'
       end if
       write(iulog,*) ' '

       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'flbc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'lower bndy timing specs'
       write(iulog,*) 'type = ',flbc_timing%type
       if( time_type == 'SERIAL' ) then
          write(iulog,*) 'year offset = ',flbc_timing%yr_offset
       else if( time_type == 'CYCLICAL' ) then
          write(iulog,*) 'year = ',flbc_timing%date/10000
       else
          write(iulog,*) 'date = ',flbc_timing%date
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',flbc_cnt,' species with specified lower bndy values'
       write(iulog,*) ' '

       !-----------------------------------------------------------------------
       ! 	... get timing information, allocate arrays, and read in dates
       !-----------------------------------------------------------------------
       call getfil ( flbc_file, filename, 0)
       call wrap_open (trim(filename), NF_NOWRITE, ncid)
       call wrap_inq_dimid( ncid, 'time', dimid )
       call wrap_inq_dimlen( ncid, dimid, ntimes )
    endif
#if ( defined SPMD )
    call mpibcast( ntimes, 1, mpiint, 0, mpicom )
#endif

    allocate( dates(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate dates array; error = ',astat
       call endrun
    end if
    allocate( times(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate times array; error = ',astat
       call endrun
    end if

    if (masterproc) then
       call wrap_inq_varid( ncid, 'date', varid )
       call wrap_get_var_int( ncid, varid, dates )
    endif
#if ( defined SPMD )
    call mpibcast( dates, ntimes, mpiint, 0, mpicom )
#endif

    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    if( time_type /= 'CYCLICAL' ) then
       if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       do n = 2,ntimes
          if( wrk_time <= times(n) ) then
             exit
          end if
       end do
       tim_ndx(1) = n - 1
    else
       yr = flbc_timing%date/10000
       do n = 1,ntimes
          if( yr == dates(n)/10000 ) then
             exit
          end if
       end do
       if( n >= ntimes ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       tim_ndx(1) = n
    end if
    select case( time_type )
    case( 'FIXED' )
       tim_ndx(2) = n
    case( 'CYCLICAL' )
       do n = tim_ndx(1),ntimes
          if( yr /= dates(n)/10000 ) then
             exit
          end if
       end do
       tim_ndx(2) = n - 1
       if( (tim_ndx(2) - tim_ndx(1)) < 2 ) then
          write(iulog,*) 'flbc_inti: cyclical lb conds require at least two time points'
          call endrun
       end if
    case( 'SERIAL' )
       tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
    end select
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)

    if( masterproc ) then
       write(iulog,*) ' '
       write(iulog,*) 'flbc time cnt = ',ntimes
       write(iulog,*) 'flbc times'
       write(iulog,'(10i10)') dates(:)
       write(iulog,'(1p,5g15.7)') times(:)
       write(iulog,*) 'flbc time indicies = ',tim_ndx(:)
       write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
       write(iulog,*) ' '
    endif

    do m = 1,flbc_cnt
       !-----------------------------------------------------------------------
       ! 	... allocate array
       !-----------------------------------------------------------------------
       allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
       if( astat/= 0 ) then
          write(iulog,*) 'flbc_inti: failed to allocate lbc vmr; error = ',astat
          call endrun
       end if
       !-----------------------------------------------------------------------
       ! 	... readin the flbc vmr
       !-----------------------------------------------------------------------
       call flbc_get( ncid, flbcs(m), .true. )
    end do

    !-----------------------------------------------------------------------
    ! 	... close the file
    !-----------------------------------------------------------------------
    if( masterproc ) then
       call wrap_close( ncid )
    endif
  end subroutine flbc_inti

  subroutine flbc_chk( )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    integer                     :: ncid
    real(r8)                        :: wrk_time
    integer ::  yr, mon, day

#include <netcdf.inc>

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    if( flbc_cnt > 0 .and. flbc_timing%type == 'SERIAL' ) then
       wrk_time = flt_date( ncdate + flbc_timing%yr_offset*10000, ncsec )
       if( wrk_time > times(tim_ndx(2)) ) then
          tcnt = tim_ndx(2) - tim_ndx(1)
          tim_ndx(1) = tim_ndx(2)
          tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
          t1 = tim_ndx(1)
          t2 = tim_ndx(2)
!!$          if( tcnt /= (t2 - t1) ) then
             !-----------------------------------------------------------------------
             ! 	... allocate array
             !-----------------------------------------------------------------------
             do m = 1,flbc_cnt
                if( associated( flbcs(m)%vmr ) ) then
                   deallocate( flbcs(m)%vmr,stat=astat )
                   if( astat/= 0 ) then
                      write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                      call endrun
                   end if
                end if
                allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                   call endrun
                end if
             end do
!!$          end if

          if( masterproc ) then
             call wrap_open (trim(filename), NF_NOWRITE, ncid)
          endif
          !-----------------------------------------------------------------------
          ! 	... readin the lb concentrations
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
!!$             call flbc_get( ncid, flbcs(m), .false. )
             call flbc_get( ncid, flbcs(m), .true. )
          end do

          !-----------------------------------------------------------------------
          ! 	... close the file
          !-----------------------------------------------------------------------
          if( masterproc ) then
             call wrap_close( ncid )
          endif

       end if
    end if

  end subroutine flbc_chk

  subroutine flbc_get( ncid, lbcs, initial )
    !-----------------------------------------------------------------------
    !       ... read lower bndy values
    !-----------------------------------------------------------------------
    use dyn_grid,      only : get_dyn_grid_parm
    use mo_constants,  only : d2r
    use mo_regrider,   only : regrid_inti, regrid_2d, regrid_lat_limits
    use commap,        only : clat, clon
    use phys_grid,     only: scatter_field_to_chunk

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)           :: ncid
    logical, intent(in)           :: initial
    type(flbc), intent(inout) :: lbcs

#include <netcdf.inc>
    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: j, jl, ju, m               ! Indices
    integer                     :: t1, t2, tcnt
    integer                     :: ierr
    integer                     :: vid, nlat, nlon
    integer                     :: dimid_lat, dimid_lon
    integer                     :: plon, plat
    real(r8), allocatable           :: lat(:)
    real(r8), allocatable           :: lon(:)
    real(r8), allocatable           :: wrk(:,:,:), wrk_zonal(:,:)
    real(r8), allocatable           :: wrk2d(:,:)
    character(len=nf_max_name)  :: varname
    real(r8), allocatable       :: glob_vmr(:,:,:)
    real(r8), allocatable       :: locl_vmr(:,:,:)
    integer :: ndims, t

    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1
    allocate( locl_vmr(pcols,begchunk:endchunk,tcnt), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'srf_emis_get: locl_emis allocation error = ',ierr
       call endrun
    end if

    locl_vmr(:,:,:) = 0._r8

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    allocate( glob_vmr(plon,plat,tcnt), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'srf_emis_get: glob_emis allocation error = ',ierr
       call endrun
    end if

    Masterproc_only :  if( masterproc ) then

       initialization :  if( initial ) then
          !-----------------------------------------------------------------------
          !       ... get grid dimensions from file
          !-----------------------------------------------------------------------
          !           latitudes
          !-----------------------------------------------------------------------
          call wrap_inq_dimid( ncid, 'lat', dimid_lat )
          call wrap_inq_dimlen( ncid, dimid_lat, nlat )
          allocate( lat(nlat),stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'flbc_get: lat allocation error = ',ierr
             call endrun
          end if
          call wrap_inq_varid( ncid, 'lat', vid )
          call wrap_get_var_realx( ncid, vid, lat )
          lat(:nlat) = lat(:nlat) * d2r

          !-----------------------------------------------------------------------
          !           longitudes
          !-----------------------------------------------------------------------
          call wrap_inq_dimid( ncid, 'lon', dimid_lon )
          call wrap_inq_dimlen( ncid, dimid_lon, nlon )
          allocate( lon(nlon),stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'flbc_get: lon allocation error = ',ierr
             call endrun
          end if
          call wrap_inq_varid( ncid, 'lon', vid )
          call wrap_get_var_realx( ncid, vid, lon )
          lon(:nlon) = lon(:nlon) * d2r
          !-----------------------------------------------------------------------
          !       ... set up regridding
          !-----------------------------------------------------------------------

          gndx = regrid_inti( nlat, plat, &
               nlon, plon, &
               lon,  clon(:,1), &
               lat,  clat, &
               0, &
               do_lons=.true.,do_lats=.true. )
          deallocate( lat,lon,stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'flbc_get: Failed to deallocate lat,lon; ierr = ',ierr
             call endrun
          end if
          jl   = 1
          ju   = plat
          if( gndx /= 0 ) then
             jlim = regrid_lat_limits( gndx )
          else
             jlim = (/ jl,ju /)
          end if
          write(iulog,'(1x,''flbc_get: gndx='',i2,'', grid limits = '',2i4,'', jl,ju='',2i4)') gndx,jlim,jl,ju
       end if initialization

       allocate( wrk(nlon,jlim(1):jlim(2),tcnt), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: wrk allocation error = ',ierr
          call endrun
       end if
       !-----------------------------------------------------------------------
       !       ... read data
       !-----------------------------------------------------------------------
       varname = trim(lbcs%species) // '_LBC'
       call wrap_inq_varid( ncid, trim(varname), vid )
       call wrap_inq_varndims (ncid, vid, ndims)
       allocate(wrk2d(plon,plat),stat=ierr)
       if (ndims==2) then
         allocate( wrk_zonal(nlat,tcnt), stat=ierr )
         if( ierr /= 0 ) then
            write(iulog,*) 'flbc_get: wrk_zonal allocation error = ',ierr
            call endrun
         end if
       endif

       if (ndims==2) then
          call wrap_get_vara_realx( ncid, vid, (/ 1, t1/), &
               (/ nlat, tcnt /), wrk_zonal )
          do t = 1,tcnt
             do j = 1,nlat
                wrk(:nlon,j,t) = wrk_zonal(j,t)
             enddo
          enddo
       else
          call wrap_get_vara_realx( ncid, vid, (/ 1, jlim(1), t1/), &
               (/ nlon, jlim(2)-jlim(1)+1, tcnt /), wrk )
       endif

       do m = 1,tcnt
          call regrid_2d( wrk(:,jlim(1):jlim(2),m), wrk2d, gndx, jl, ju, &
               do_poles=.true. )
          glob_vmr(:,:,m) = wrk2d(:,:)
       end do
       deallocate( wrk,stat=ierr )
       deallocate( wrk2d,stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: Failed to deallocate wrk, ierr = ',ierr
          call endrun
       end if
       if (ndims==2) then
          deallocate( wrk_zonal,stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'flbc_get: Failed to deallocate wrk_zonal, ierr = ',ierr
             call endrun
          end if
       end if

    end if Masterproc_only

    call scatter_field_to_chunk(1, 1, tcnt, plon, glob_vmr, locl_vmr)

    deallocate(glob_vmr, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: Failed to deallocate glob_vmr; ierr = ',ierr
       call endrun
    end if

    do m = t1,t2
       lbcs%vmr(:,:,m) = locl_vmr(:,:,m-t1+1)
    enddo

    deallocate(locl_vmr, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: Failed to deallocate locl_vmr; ierr = ',ierr
       call endrun
    end if

  end subroutine flbc_get

  subroutine flbc_set( vmr, ncol, lchnk )
    !--------------------------------------------------------
    !	... set the lower bndy values
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    ::   ncol
    integer,  intent(in)    ::   lchnk
    real(r8), intent(inout) ::   vmr(ncol,pver,pcnstm1)    ! lower bndy concentrations( mol/mol )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    integer  ::  last, next
    integer  ::  wrk_date, wrk_sec
    integer  ::  tcnt
    integer  ::  astat
    real(r8)     ::  dels
    real(r8)     ::  wrk_time

    if( flbc_cnt < 1 ) then
       return
    end if
    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    wrk_sec  = ncsec
    select case( flbc_timing%type )
    case( 'SERIAL' )
       wrk_date = ncdate + flbc_timing%yr_offset*10000
    case( 'CYCLICAL' )
       wrk_date = (flbc_timing%date/10000)*10000 + mod( ncdate,10000 )
    case( 'FIXED' )
       wrk_date = flbc_timing%date
       wrk_sec  = 0
    end select
    wrk_time = flt_date( wrk_date, wrk_sec )

    !--------------------------------------------------------
    !	... set time interpolation factor
    !--------------------------------------------------------
    if( flbc_timing%type /= 'CYCLICAL' ) then
       do n = tim_ndx(1)+1,tim_ndx(2)
          if( wrk_time <= times(n) ) then
             last = n - 1
             next = n
             exit
          end if
       end do
       if( n > ntimes ) then
          write(iulog,*) 'flbc_set: interp time is out of bounds'
          call endrun
       end if
       dels = (wrk_time - times(last))/(times(next) - times(last))
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    else
       tcnt = tim_ndx(2) - tim_ndx(1) + 1
       call findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
       if( n < tcnt ) then
          last = tim_ndx(1) + n - 1
          next = last + 1
          dels = (wrk_time - times(last))/(times(next) - times(last))
       else
          next = tim_ndx(1)
          last = tim_ndx(2)
          dels = wrk_time - times(last)
          if( dels < 0._r8 ) then
             dels = 365._r8 + dels
          end if
          dels = dels/(365._r8 + times(next) - times(last))
       end if
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    end if

    dels = max( min( 1._r8,dels ),0._r8 )

    do m = 1,flbc_cnt
       n = flbcs(m)%spc_ndx
       vmr(:ncol,pver,n) = flbcs(m)%vmr(:ncol,lchnk,last) &
            + dels * (flbcs(m)%vmr(:ncol,lchnk,next) - flbcs(m)%vmr(:ncol,lchnk,last))
    end do

  end subroutine flbc_set

end module mo_flbc
