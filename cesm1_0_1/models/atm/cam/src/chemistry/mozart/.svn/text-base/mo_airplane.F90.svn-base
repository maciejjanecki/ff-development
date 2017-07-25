module mo_airplane
  !--------------------------------------------------------------------
  !	... Airplane insitu emission sources
  !--------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use pio,          only : pio_inq_dimid, pio_inq_dimlen, pio_get_var, &
       file_desc_t, var_desc_t, pio_inq_vardimid, pio_inq_varndims, pio_nowrite, &
       pio_inq_varid, pio_closefile
       
  use cam_pio_utils,only : cam_pio_openfile
  use cam_logfile,  only : iulog
  implicit none

  private
  public  :: air_altitude, pno, pco, airpl_src, has_airpl_src


  save 

  real(r8), allocatable :: &
       pno(:,:,:), &
       pco(:,:,:), &
       air_altitude(:)

  logical :: has_airpl_src = .false.

contains

  subroutine airpl_src( airpl_emis_file )
    !-----------------------------------------------------------------------
    ! 	... Initialize airplane emissions
    !	    Note: The emissions are read in in units of molecules/cm**2/s
    !	          on a vertically resolved grid.
    !		  Conversion to units of molec/cm**3/s is done in SETEXT
    !-----------------------------------------------------------------------
    use dyn_grid,      only : get_dyn_grid_parm, get_horiz_grid_d 
    use spmd_utils,    only : masterproc
    use mo_regrider,   only : regrid_inti, regrid_2d, regrid_lat_limits
    use chem_mods,     only : adv_mass
    use physconst,     only : rearth
    use ioFileMod,     only : getfil
    use mo_chem_utls,  only : get_spc_ndx, get_extfrc_ndx

    implicit none

    !-----------------------------------------------------------------------
    ! 	... Dummy args
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: airpl_emis_file

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: msq2cmsq = 1.e4_r8
    integer  :: ios, k, jl, ju, j
    integer  :: jlim(2)
    integer  :: nlat, nlon, nlev, ndims
    integer  :: ierr
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    integer  :: dimid_lat, dimid_lon, dimid_lev
    integer  :: gndx
    integer  :: dimid(3)
    real(r8) :: pi, d2r
    real(r8), allocatable :: lat(:), lon(:), sf(:)
    real(r8), allocatable :: pno_in(:,:,:), pco_in(:,:,:)
    real(r8)  :: total(3)
    real(r8)  :: seq, factor
    character(len=256) :: locfn
    integer :: co_ndx, no_ndx
    integer :: plat, plon
    real(r8), allocatable :: latwts(:), clat(:), clon(:)    

    co_ndx = get_extfrc_ndx('CO')
    no_ndx = get_extfrc_ndx('NO')

    if ( co_ndx < 0 .and. no_ndx < 0 ) then
       if( masterproc ) then
          write(iulog,*) 'airpl_src: NO and CO do not have external source --> no aircraft sources will be applied'      
       endif
       return
    endif

    if ( len_trim(airpl_emis_file) == 0 ) then
       return
    endif

    has_airpl_src = .true.

    co_ndx = get_spc_ndx('CO')
    no_ndx = get_spc_ndx('NO')
    plat = get_dyn_grid_parm('plat')
    plon = get_dyn_grid_parm('plon')

    !-----------------------------------------------------------------------
    !	... Open NetCDF file
    !-----------------------------------------------------------------------
    call getfil (airpl_emis_file, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)

    pi  = 4._r8 * atan( 1._r8 )
    d2r = pi/180._r8

    !-----------------------------------------------------------------------
    !       ... Get grid dimensions from file
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
    ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lat', vid )
    ierr = pio_get_var( piofile, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r

    ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
    ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lon', vid )
    ierr = pio_get_var( piofile, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    ierr = pio_inq_dimid( piofile, 'altitude', dimid_lev )
    ierr = pio_inq_dimlen( piofile, dimid_lev, nlev )
    allocate( air_altitude(nlev+1), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: air_altitude allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'altitude', vid )
    ierr = pio_get_var( piofile, vid, air_altitude(1:nlev) )
    air_altitude(nlev+1) = air_altitude(nlev) + (air_altitude(nlev) - air_altitude(nlev-1))

    !-----------------------------------------------------------------------
    !       ... Set up regridding
    !-----------------------------------------------------------------------
    allocate(latwts(plat), clat(plat), clon(plon))

    call get_horiz_grid_d(plat, clat_d_out=clat, wght_d_out=latwts)
    call get_horiz_grid_d(plon, clon_d_out=clon)

    gndx = regrid_inti( nlat, plat, &
         nlon, plon, &
         lon,  clon, &
         lat,  clat, &
         0, &
         do_lons=.true.,do_lats=.true. )
    deallocate( lat, lon, clat, clon, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: Failed to deallocate lat,lon; ierr = ',ierr
       call endrun
    end if
    jl   = 1
    ju   = plat
    if( gndx /= 0 )then
       jlim = regrid_lat_limits( gndx )
    else
       jlim = (/ jl,ju /)
    end if
    allocate( pno_in(nlon,jlim(1):jlim(2),nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pno_in allocation error = ',ierr
       call endrun
    end if
    allocate( pco_in(nlon,jlim(1):jlim(2),nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pco_in allocation error = ',ierr
       call endrun
    end if
    allocate(pno(plon,plat,nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pno allocation error = ',ierr
       call endrun
    end if
    allocate( pco(plon,plat,nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pco allocation error = ',ierr
       call endrun
    end if

    !-----------------------------------------------------------------------
    !	... Read emissions
    !-----------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'nox', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )
    if( ndims /= 3 ) then
       write(iulog,*) 'airpl_src: variable nox has ndims = ',ndims,', expecting 3'
       call endrun
    end if
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    if( dimid(1) /= dimid_lon  .or. dimid(2) /= dimid_lat .or.  dimid(3) /= dimid_lev ) then
       write(iulog,*) 'airpl_src: Dimensions in wrong order for variable nox'
       write(iulog,*) '...      Expecting (lon, lat, lev)'
       call endrun
    end if
    ierr = pio_get_var( piofile, vid, &
         (/ 1, jlim(1), 1/), &                    ! start
         (/ nlon, jlim(2)-jlim(1)+1, nlev /), &   ! count
         pno_in )

    ierr = pio_inq_varid( piofile, 'co', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )

    if( ndims /= 3 ) then
       write(iulog,*) 'READ_SFLX: variable co has ndims = ',ndims,', expecting 3'
       call endrun
    end if
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. dimid(3) /= dimid_lev ) then
       write(iulog,*) 'airpl_src: Dimensions in wrong order for variable co'
       write(iulog,*) '...      Expecting (lon, lat, lev)'
       call endrun
    end if
    ierr = pio_get_var( piofile, vid, &
         (/ 1, jlim(1), 1/), &                    ! start
         (/ nlon, jlim(2)-jlim(1)+1, nlev /), &   ! count
         pco_in )
    call pio_closefile( piofile )

    !-----------------------------------------------------------------------
    !	... Regrid emissions
    !-----------------------------------------------------------------------
    do k = 1,nlev
       call regrid_2d( pno_in(:,jlim(1):jlim(2),k), pno(:,:,k), gndx, jl, ju, do_poles=.true. )
       call regrid_2d( pco_in(:,jlim(1):jlim(2),k), pco(:,:,k), gndx, jl, ju, do_poles=.true. )
    end do
    deallocate( pno_in, pco_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: Failed to deallocate pno_in,pco_in; ierr = ',ierr
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... Get global emission from this source
    !-----------------------------------------------------------------------
    seq = 2._r8*pi*msq2cmsq*rearth**2/real(plon,r8)
    allocate(sf(plat))
    do j = 1,plat
       sf(j) = seq*latwts(j)
    end do
    deallocate(latwts)

    factor = 86400._r8 * 365._r8 &   ! sec / year
         / 6.022e23_r8 &           ! molec / mole
         * 1.e-12_r8               ! Tg / g

    total(:) = 0._r8
    do j = 1,plat
       total(1) = total(1) + sum( pno(:plon,j,:nlev ) ) * sf(j)
       total(2) = total(2) + sum( pco(:plon,j,:nlev ) ) * sf(j)
    end do
    deallocate(sf)

    if(masterproc) write(iulog,*) 'airpl_src: nlev = ',nlev
    !-----------------------------------------------------------------------
    !       ... Convert totals from molec cm^-2 s^-1 to Tg y^-1
    !-----------------------------------------------------------------------
    if (no_ndx .gt. 0) then
       total(1) = total(1) * adv_mass(no_ndx) * factor
       if(masterproc) write(iulog,'('' airpl_src Aircraft emissions: '',a6,'' = '',f10.3,1X,a6)') 'NO',total(1),'TgN/y'
    endif
    if (co_ndx .gt. 0) then
       total(2) = total(2) * adv_mass(co_ndx) * factor
       if(masterproc) write(iulog,'('' airpl_src Aircraft emissions: '',a6,'' = '',f10.3,1X,a6)') 'CO',total(2),'Tg/y'
    endif


  end subroutine airpl_src

end module mo_airplane
