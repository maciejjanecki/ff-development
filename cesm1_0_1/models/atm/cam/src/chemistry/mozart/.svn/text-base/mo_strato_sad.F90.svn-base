
module mo_strato_sad
!---------------------------------------------------------------
! 	... prescribed strat aero surf area density module
!---------------------------------------------------------------

  use m_types,      only : time_ramp
  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pver, pcols, pverp
  use time_manager, only : get_curr_date
  use time_utils,   only : flt_date, moz_findplb
  use abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use wrap_nf
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog
  use dyn_grid, only : get_dyn_grid_parm

  implicit none

  private
  public  :: strato_sad_inti
  public  :: strato_sad_timestep_init
  public  :: strato_sad_set

  save

  integer, parameter    :: time_span = 1

  integer               :: ntimes
  integer               :: nlon
  integer               :: nlat
  integer               :: nlev
  integer               :: gndx = 0
  integer               :: tim_ndx(2)
  integer               :: jlim(2)
  integer,  allocatable :: dates(:)
  real(r8), allocatable :: times(:)
  real(r8), allocatable :: sad_lats(:)
  real(r8), allocatable :: sad_lons(:)
  real(r8), allocatable :: sad_levs(:)
  real(r8), allocatable :: sage_sad(:,:,:,:)
  character(len=256)     :: filename

  logical :: has_sulfate_rxts
  type(time_ramp)   :: strato_sad_timing

contains

  subroutine strato_sad_inti( sad_file, sad_timing )
!-----------------------------------------------------------------------
! 	... initialize the strato sad module
!-----------------------------------------------------------------------
    
    use mo_constants,  only : d2r
    use commap,        only : clat, clon
    use mo_regrider,   only : regrid_inti, regrid_lat_limits
    use mo_chem_utls,  only : get_rxt_ndx
    use mo_strato_rates,only : has_strato_chem
    use string_utils,  only : to_upper

    implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: sad_file
    type(time_ramp),  intent(in) :: sad_timing

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
    real(r8), parameter :: hPa2Pa = 100.
    integer :: astat
    integer :: j, l, m, n                           ! Indices
    integer :: t1, t2
    integer :: ncid
    integer :: dimid
    integer :: varid
    integer :: yr, wrk_date, wrk_sec
    integer :: plon, plat
    real(r8)    :: seq
    real(r8)    :: wrk_time
    character(len=8)  :: time_type

    integer :: mon, day, tod, ncdate, ncsec

    integer :: usrrxt_ndx(4)

    usrrxt_ndx(1) = get_rxt_ndx( 'usr_N2O5_aer' )
    usrrxt_ndx(2) = get_rxt_ndx( 'usr_NO3_aer' )
    usrrxt_ndx(3) = get_rxt_ndx( 'usr_NO2_aer' )
    usrrxt_ndx(4) = get_rxt_ndx( 'usr_HO2_aer' )

    has_sulfate_rxts = any( usrrxt_ndx > 0 )
    has_sulfate_rxts = has_sulfate_rxts .or. has_strato_chem

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod

!-----------------------------------------------------------------------
! 	... check timing
!-----------------------------------------------------------------------
    strato_sad_timing = sad_timing 
    strato_sad_timing%type = to_upper(strato_sad_timing%type)
    time_type = strato_sad_timing%type 


    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' .and. time_type /= 'FIXED' ) then
       write(iulog,*) 'strato_sad_inti: time type ',trim(time_type),' is not SERIAL, CYCLICAL, or FIXED'
       call endrun
    end if

    wrk_sec  = ncsec
    if( (time_type == 'SERIAL') ) then
       wrk_date = ncdate + strato_sad_timing%yr_offset*10000
    else if( time_type == 'CYCLICAL' ) then
       wrk_date = (strato_sad_timing%date/10000)*10000 + mod(ncdate,10000)
    else
       wrk_date = strato_sad_timing%date
       wrk_sec  = 0
    end if
    wrk_time = flt_date( wrk_date, wrk_sec )
    if (masterproc) then
       write(iulog,*) 'strato_sad_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time

       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'strato sad timing specs'
       write(iulog,*) 'type = ',strato_sad_timing%type
       if( time_type == 'SERIAL' ) then
          write(iulog,*) 'year offset = ',strato_sad_timing%yr_offset
       else if( time_type == 'CYCLICAL' ) then
          write(iulog,*) 'year = ',strato_sad_timing%date/10000
       else
          write(iulog,*) 'date = ',strato_sad_timing%date
       end if
    endif

!-----------------------------------------------------------------------
! 	... open netcdf file
!-----------------------------------------------------------------------
    call getfil ( sad_file, filename, 0 )
    call wrap_open (trim(filename), NF_NOWRITE, ncid)
!-----------------------------------------------------------------------
! 	... get timing information, allocate arrays, and read in dates
!-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'time', dimid )
    call wrap_inq_dimlen( ncid, dimid, ntimes )
    allocate( dates(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'strato_sad_inti: failed to allocate dates array; error = ',astat
       call endrun
    end if
    allocate( times(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'strato_sad_inti: failed to allocate times array; error = ',astat
       call endrun
    end if
    call wrap_inq_varid( ncid, 'date', varid )
    call wrap_get_var_int( ncid, varid, dates )
    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    if( time_type /= 'CYCLICAL' ) then
       if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
          write(iulog,*) 'strato_sad_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       do n = 2,ntimes
          if( wrk_time <= times(n) ) then
             exit
          end if
       end do
       tim_ndx(1) = n - 1
    else
       yr = strato_sad_timing%date/10000
       do n = 1,ntimes
          if( yr == dates(n)/10000 ) then
             exit
          end if
       end do
       if( n >= ntimes ) then
          write(iulog,*) 'strato_sad_inti: time out of bounds for dataset = ',trim(filename)
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
          write(iulog,*) 'strato_sad_inti: cyclical sad require at least two time points'
          call endrun
       end if
    case( 'SERIAL' )
       tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
    end select
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    
    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad time cnt = ',ntimes
       write(iulog,*) 'strato_sad times'
       write(iulog,'(10i10)') dates(:)
       write(iulog,'(1p,5g15.7)') times(:)
       write(iulog,*) 'strato_sad time indicies = ',tim_ndx(:)
       write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
       write(iulog,*) ' '
    endif

!-----------------------------------------------------------------------
!     	... inquire about latitudes
!-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'lat', dimid )
    call wrap_inq_dimlen( ncid, dimid, nlat )
!-----------------------------------------------------------------------
!     	... allocate space for latitude coordinates
!-----------------------------------------------------------------------
    allocate( sad_lats(nlat),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate latitudes array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get latitudes
!-----------------------------------------------------------------------
    call wrap_inq_varid( ncid, 'lat', varid )
    call wrap_get_var_realx( ncid, varid, sad_lats  )

!-----------------------------------------------------------------------
!     	... inquire about longitudes
!-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'lon', dimid )
    call wrap_inq_dimlen( ncid, dimid, nlon )
!-----------------------------------------------------------------------
!     	... allocate space for longitude coordinates
!-----------------------------------------------------------------------
    allocate( sad_lons(nlon),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate longitudes array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get longitudes
!-----------------------------------------------------------------------
    call wrap_inq_varid( ncid, 'lon', varid )
    call wrap_get_var_realx( ncid, varid, sad_lons )

!-----------------------------------------------------------------------
!     	... convert to radians and setup regridding
!-----------------------------------------------------------------------
    sad_lats(:) = d2r * sad_lats(:)
    sad_lons(:) = d2r * sad_lons(:)
    gndx = regrid_inti( nlat, plat, &
                        nlon, plon, &
                        sad_lons,  clon(:,1), &
                        sad_lats,  clat, &
                         0 )
    if( gndx < 0 ) then
       write(iulog,*) 'strato_sad_inti: regrid_inti error = ',gndx
       call endrun
    else
       write(iulog,*) 'strato_sad_inti: regrid_inti index = ',gndx
    end if
    
    if( gndx /= 0 ) then
       jlim = regrid_lat_limits( gndx )
    else
       jlim(:) = (/ 1, plat /)
    end if
    write(iulog,*) 'strato_sad_inti: jlim = ',jlim

!-----------------------------------------------------------------------
!     	... inquire about levels
!-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'lev', dimid )
    call wrap_inq_dimlen( ncid, dimid, nlev )
!-----------------------------------------------------------------------
!     	... allocate space for level coordinates
!-----------------------------------------------------------------------
    allocate( sad_levs(nlev),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate levels array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get levels
!-----------------------------------------------------------------------
    call wrap_inq_varid( ncid, 'lev', varid )
    call wrap_get_var_realx( ncid, varid, sad_levs )
    sad_levs(:) = hPa2Pa * sad_levs(:)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad_inti: sad_levs'
       write(iulog,'(10g12.5)') sad_levs(:)
       write(iulog,*) ' '
    endif

!-----------------------------------------------------------------------
!     	... allocate module sad variable
!-----------------------------------------------------------------------
    allocate( sage_sad(plon,plat,nlev,t1:t2),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate sage_sad array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
! 	... read and interpolate sage sad data
!-----------------------------------------------------------------------
    call strato_sad_get( ncid )
!-----------------------------------------------------------------------
! 	... close netcdf file
!-----------------------------------------------------------------------
    call wrap_close( ncid )

  end subroutine strato_sad_inti

  subroutine strato_sad_timestep_init( )
!-----------------------------------------------------------------------
!       ... check serial case for time span
!-----------------------------------------------------------------------
    implicit none

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    integer                     :: ncid
    real(r8)                    :: wrk_time

    integer ::  yr, mon, day, tod, ncdate, ncsec
    integer :: plon, plat

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod
    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')
	
    if( strato_sad_timing%type == 'SERIAL' ) then
       wrk_time = flt_date( ncdate + strato_sad_timing%yr_offset*10000, ncsec )
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
             if( allocated( sage_sad ) ) then
                deallocate( sage_sad,stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'strato_sad_timestep_init: failed to deallocate strato_sad vmr; error = ',astat
                   call endrun
                end if
             end if
             allocate( sage_sad(plon,plat,nlev,t1:t2),stat=astat )
             if( astat/= 0 ) then
                write(iulog,*) 'strato_sad_timestep_init: failed to allocate sage_sad; error = ',astat
                call endrun
             end if
!!$          end if
!-----------------------------------------------------------------------
! 	... open netcdf file
!-----------------------------------------------------------------------
          call wrap_open (trim(filename), NF_NOWRITE, ncid)
!-----------------------------------------------------------------------
! 	... read and interpolate sage sad data
!-----------------------------------------------------------------------
          call strato_sad_get( ncid )
!-----------------------------------------------------------------------
! 	... close netcdf file
!-----------------------------------------------------------------------
          call wrap_close( ncid )
       end if
    end if

  end subroutine strato_sad_timestep_init

  subroutine strato_sad_get( ncid )
!-----------------------------------------------------------------------
!       ... read sad values
!-----------------------------------------------------------------------
    
    use mo_regrider,   only : regrid_2d

    implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
    integer, intent(in)           :: ncid

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer               :: jl, ju
    integer               :: k, n                            ! Indices
    integer               :: t1, t2, tcnt
    integer               :: ierr
    integer               :: vid
    integer               :: plon, plat
    real(r8), allocatable :: sage_sad_in(:,:,:,:)
    real(r8), allocatable :: wrk3d(:,:,:)
    real(r8), allocatable :: tmp_3d(:,:,:)
    integer :: ndims, i

!-----------------------------------------------------------------------
!       ... read sage data
!-----------------------------------------------------------------------
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1

!-----------------------------------------------------------------------
!       ... allocate local wrk arrays
!-----------------------------------------------------------------------
    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    allocate( sage_sad_in(nlon,jlim(1):jlim(2),nlev,t1:t2), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'strato_sad_get: wrk allocation error = ',ierr
       call endrun
    end if
    allocate( wrk3d(plon,plat,nlev),stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'strato_sad_get: failed to allocate wrk3d; error = ',ierr
       call endrun
    end if
    call wrap_inq_varid(  ncid, 'sad_sage', vid )
    call wrap_inq_varndims (ncid, vid, ndims)

    if ( ndims == 4 ) then
       call wrap_get_vara_realx( ncid, vid, &
            (/ 1, jlim(1), 1, t1/), &                     ! start
            (/ nlon, jlim(2)-jlim(1)+1, nlev, tcnt /), &  ! count
            sage_sad_in   )
    else
       allocate( tmp_3d(jlim(1):jlim(2),nlev,t1:t2), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'strato_sad_get: tmp_3d allocation error = ',ierr
          call endrun
       end if

       call wrap_get_vara_realx( ncid, vid, &
            (/ jlim(1), 1, t1/), &                     ! start
            (/ jlim(2)-jlim(1)+1, nlev, tcnt /), &  ! count
            tmp_3d   )

       do i = 1,nlon
          sage_sad_in( i, :,:,:) = tmp_3d(:,:,:)
       enddo

       deallocate( tmp_3d, stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'strato_sad_get: failed to deallocate tmp_3d, ierr = ',ierr
          call endrun
       end if

    endif

    jl = 1
    ju = plat
    do n = t1,t2
       if( gndx /= 0 ) then
          do k = 1,nlev
             call regrid_2d( sage_sad_in(:,jlim(1):jlim(2),k,n), wrk3d(:,:,k), gndx, jl, ju, do_poles=.true. ) 
          end do
       else
          do k = 1,nlev
             wrk3d(:plon,:plat,k) = sage_sad_in(:plon,jlim(1):jlim(2),k,n)
          end do
       end if
       do k = 1,nlev
          sage_sad(:,:,k,n) = wrk3d(:,:,k)
       end do
    end do
    deallocate( sage_sad_in, wrk3d, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'strato_sad_get: failed to deallocate sage_sad_in,wrk3d, ierr = ',ierr
       call endrun
    end if

  end subroutine strato_sad_get

  subroutine strato_sad_set( pmid, latndx, lonndx, sad, ncol )
!--------------------------------------------------------
!	... set the sad values
!--------------------------------------------------------
    
    implicit none

!--------------------------------------------------------
!	... dummy arguments
!--------------------------------------------------------
    integer, intent(in)     :: ncol
    integer, intent(in)     :: latndx(pcols)    ! latitude index
    integer, intent(in)     :: lonndx(pcols)    ! longitude tile index
    real(r8), intent(in)    :: pmid(pcols,pver) ! midpoint pressure (Pa)
    real(r8), intent(inout) :: sad(pcols, pver) ! stratospheric aerosol surface area density (1/cm)

!--------------------------------------------------------
!	... local variables
!--------------------------------------------------------
    integer  ::  i, k, n
    integer  ::  last, next
    integer  ::  wrk_date, wrk_sec
    integer  ::  tcnt
    integer  ::  astat
    real(r8) ::  dels
    real(r8) ::  wrk_time
    real(r8) ::  sad_int(ncol,pver,2)
    integer  ::  yr, mon, day, tod, ncdate, ncsec

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod

!--------------------------------------------------------
!	... setup the time interpolation
!--------------------------------------------------------
    wrk_sec  = ncsec
    select case( strato_sad_timing%type )
    case( 'SERIAL' )
       wrk_date = ncdate + strato_sad_timing%yr_offset*10000
    case( 'CYCLICAL' )
       wrk_date = (strato_sad_timing%date/10000)*10000 + mod( ncdate,10000 )
    case( 'FIXED' )
       wrk_date = strato_sad_timing%date
       wrk_sec  = 0
    end select
    wrk_time = flt_date( wrk_date, wrk_sec )

!--------------------------------------------------------
!	... set time interpolation factor
!--------------------------------------------------------
    if( strato_sad_timing%type /= 'CYCLICAL' ) then
       do n = tim_ndx(1)+1,tim_ndx(2)
          if( wrk_time <= times(n) ) then
             last = n - 1
             next = n
             exit
          end if
       end do
       if( n > ntimes ) then
          write(iulog,*) 'strato_sad_set: interp time is out of bounds'
          call endrun
       end if
       dels = (wrk_time - times(last))/(times(next) - times(last))
    else
       tcnt = tim_ndx(2) - tim_ndx(1) + 1
       call moz_findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
       if( n < tcnt ) then
          last = tim_ndx(1) + n - 1
          next = last + 1
          dels = (wrk_time - times(last))/(times(next) - times(last))
       else
          next = tim_ndx(1)
          last = tim_ndx(2)
          dels = wrk_time - times(last)
          if( dels < 0. ) then
             dels = 365. + dels
          end if
          dels = dels/(365. + times(next) - times(last))
       end if
    end if
    
    dels = max( min( 1._r8,dels ),0._r8 )

#ifdef DEBUG
    write(iulog,*) ' '
    write(iulog,*) 'strato_sad_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    write(iulog,*) 'strato_sad_set: dates(last),dates(next)     = ',dates(last),dates(next)
#endif

    call vinterp( pmid, sage_sad(:,:,:,last), sad_int(:,:,1), ncol, lonndx, latndx )
    call vinterp( pmid, sage_sad(:,:,:,next), sad_int(:,:,2), ncol, lonndx, latndx )

#ifdef DEBUG
    write(iulog,*) 'strato_sad_set: pmid(1,:)'
    write(iulog,'(1p,5g15.7)') pmid(1,:)
    write(iulog,*) 'strato_sad_set: sad_levs'
    write(iulog,'(1p,5g15.7)') sad_levs(:)
    write(iulog,*) 'strato_sad_set: sage_sad(last)'
    write(iulog,'(1p,5g15.7)') sage_sad(1,:,lat,ip,last)
    write(iulog,*) 'strato_sad_set: sage_sad(next)'
    write(iulog,'(1p,5g15.7)') sage_sad(1,:,lat,ip,next)
#endif

    do i=1,ncol
       sad(i,:) = sad_int(i,:,1) + dels * (sad_int(i,:,2) - sad_int(i,:,1))
    enddo

#ifdef DEBUG
    write(iulog,*) 'strato_sad_set: sad'
    write(iulog,'(1p,5g15.7)') sad(1,:)
    write(iulog,*) ' '
    !call endrun 
#endif

  end subroutine strato_sad_set

  subroutine vinterp( pmid, sad_src, sad_int, ncol, lonndx, latndx )
!-----------------------------------------------------------------------
!   	... vertically interpolate input data
!-----------------------------------------------------------------------
    implicit none

!-----------------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------------
    integer,  intent(in)  :: ncol
    integer,  intent(in), dimension(:) :: lonndx,latndx
    real(r8), intent(in)  :: pmid(:,:)
    real(r8), intent(in)  :: sad_src(:,:,:)
    real(r8), intent(out) :: sad_int(:,:)

!-----------------------------------------------------------------------
!   	... local variables
!-----------------------------------------------------------------------
    integer :: i, lon,lat
    integer :: k, kl, ku
    real(r8)    :: delp, pinterp

    level_loop : do k = 1,pver
       long_loop : do i = 1,ncol
          lon = lonndx(i)
          lat = latndx(i)            
          pinterp = pmid(i,k)
          if( pinterp <= sad_levs(1) ) then
             sad_int(i,k) = sad_src(lon,lat,1)
          else if( pinterp > sad_levs(nlev) ) then
             sad_int(i,k) = sad_src(lon,lat,nlev)
          else
             do ku = 2,nlev
                if( pinterp <= sad_levs(ku) ) then
                   kl = ku - 1
                   delp = log( pinterp/sad_levs(kl) ) &
                        / log( sad_levs(ku)/sad_levs(kl) )
                   sad_int(i,k) = sad_src(lon,lat,kl) + delp * (sad_src(lon,lat,ku) - sad_src(lon,lat,kl))
                   exit
                end if
             end do
          end if
       end do long_loop
    end do level_loop

  end subroutine vinterp

end module mo_strato_sad
