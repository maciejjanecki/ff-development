module tgcm_forcing

!
! Provides TGCM forcing data for use without interactive chemistry
!
  use ppgrid
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use spmd_utils,         only: masterproc
  use pmgrid,             only: plon, plat, plev
  use physconst,          only: cpair
  use abortutils,         only: endrun
  use interpolate_data,   only: lininterp
  use wrap_nf
  use cam_logfile,        only: iulog

  implicit none
  private
  save

! Public interfaces
  public &
       tgcm_init,          &
       tgcm_timestep_init, &
       tgcm_get_o3,        &   ! return tgcm ozone mixing ratio
       tgcm_get_cnst,      &   ! return tgcm constituents for nlte (except o3)
       tgcm_get_solar          ! return tgcm net solar heating rate

  real(r8), parameter, public :: preftgcm=5.e-5_r8    ! TIME GCM reference pressure (Pa)

! Private module data

  integer :: itgcmcyc                 ! flag for cycling TIME/GCM input dataset:
!  = 0 : one time sample only
!  = 1 : full annual cycle (default)
!  = 2 : two time samples
  character(len=256) :: cftgcm        ! Pathname of time-variant TIME/GCM output

  integer, parameter :: plontgcm=72
  integer, parameter :: plattgcm=36
  integer, parameter :: plevtgcm=45
  integer, parameter :: ptime=12

  integer i,ix,js

  logical donlte                         ! True => non-LTE calculation

  integer ncid_tgcm  ! TIME/GCM file

  real(r8), allocatable, dimension(:,:,:,:) :: o1     ! O   interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: o2     ! O2  interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: o3     ! O3  interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: n2     ! N2  interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: co2    ! CO2 interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: no     ! NO interpolated to CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:,:) :: qsolar ! QSOLAR interpolated to CCM grid (K/s)
  real(r8) caldtgcm(ptime)               ! TIME/GCM calendar day for each sample in input dataset

  real(r8), allocatable, dimension(:,:,:) :: o1c          ! O   instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: o2c          ! O2  instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: o3c          ! O3  instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: n2c          ! N2  instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: co2c         ! CO2 instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: noc          ! NO instantaneous in CCM grid (mmr)
  real(r8), allocatable, dimension(:,:,:) :: qsolarc      ! QSOLAR instantaneous in CCM grid (K/s)

  real(r8) t_top_msis(plon,1,plat,ptime)
  real(r8) t_top_msisc(plon,plat)

  real(r8) z_top_msis(plon,1,plat,ptime)
  real(r8) z_top_msisc(plon,plat)

  real(r8),target, allocatable, dimension(:,:,:) :: o1_to_chunk       ! instant. O1 scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: o2_to_chunk       ! instant. O2 scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: o3_to_chunk       ! instant. O3 scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: n2_to_chunk       ! instant. N2 scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: co2_to_chunk      ! instant. CO2 scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: no_to_chunk       ! instant. NO scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: qsolar_to_chunk   ! instant. QSOLAR scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: ttop_to_chunk     ! instant. TTOP scattered to chunks
  real(r8),target, allocatable, dimension(:,:,:) :: zalt_to_chunk     ! instant. ZALT scattered to chunks

!================================================================================================
contains
!================================================================================================

  subroutine tgcm_init (hypm, nlte_use_mo, use_data_o3, itgcmcyc_in, cftgcm_in)
!
!     This routine reads a netCDF file containing instantaneous fields
!     from TIME/GCM output. Fields are:
!     O,O2,O3,N2,CO2 : used in Fomichev's mods
!     QSOLAR         : SW heating (to disappear)
!
!------------------------------------------------------------------------
    use ioFileMod
    use pspect
    use commap
    use constituents,     only: cnst_mw, cnst_get_ind
    use physconst,        only: mwco2
    use cam_history,      only: add_default, addfld, phys_decomp
    use interpolate_data, only: bilin

    real(r8),         intent(in) :: hypm(plev)
    logical,          intent(in) :: nlte_use_mo
    logical,          intent(in) :: use_data_o3
    integer,          intent(in) :: itgcmcyc_in
    character(len=*), intent(in) :: cftgcm_in

!------------------------------------------------------------------------

!---------------------------Local workspace-------------------------------------------
    integer londimid                        ! netcdf id for longitude dimension
    integer latdimid                        ! netcdf id for latitude dimension
    integer levdimid                        ! netcdf id for level dimension
    integer lonid                           ! netcdf id for longitude variable
    integer latid                           ! netcdf id for latitude variable
    integer levid                           ! netcdf id for level variable
    integer tdimid                          ! netcdf id for time dimension
    integer mtdimid                         ! netcdf id for mtime dimension
    integer lonsiz                          ! netcdf longitude size
    integer latsiz                          ! netcdf latitude size
    integer levsiz                          ! netcdf level size
    integer timesiz                         ! netcdf time size
    integer mtimesiz                        ! netcdf mtime size
    integer noid                            ! netcdf NO id
    integer o1id                            ! netcdf O1 id
    integer o2id                            ! netcdf O2 id
    integer o3id                            ! netcdf O3 id
    integer n2id                            ! netcdf N2 id
    integer co2id                           ! netcdf CO2 id
    integer qsolarid                        ! netcdf QSOLAR id
    integer tnid                            ! netcdf TN id
    integer zaltid                            ! netcdf TN id
    integer mtimeid                         ! netcdf MTIME id
    integer dimids(4)                       ! variables shape
    integer k                               ! index
    integer j                               ! index
    integer i                               ! index
    integer l
    integer cnt4(4)                         ! array of counts for each dimension
    integer strt4(4)                        ! array of starting indices
    integer cnt3(3)                         ! array of counts for each dimension
    integer strt3(3)                        ! array of starting indices
    integer cnt2(2)                         ! array of coutns for each dim (MTIME only)
    integer strt2(2)                        ! array of starting indices (MTIME only)
    integer ierr                            ! error status
    integer nt                              ! no. of time samples in TIME/GCM file
    integer kinv                            ! index
    integer status                          ! returned by netCDF file


    real(r8) lontgcm(plontgcm)              ! TIME/GCM longitude array
    real(r8) lattgcm(plattgcm)              ! TIME/GCM latitude array
    real(r8) xtgcm(plevtgcm)                ! TIME/GCM psh
    real(r8) xtgcm_ccm(plevtgcm)            ! TIME/GCM psh mapped using CCM Pref
    real(r8) xccm(plev)                     ! ccm psh

    real(r8) dummo(plattgcm,plev)           ! dummy array (used for vertical interpolation)

    character*256 locfn                     ! input file name

    real(r8) lontgcmt(plontgcm+1)           ! TIME/GCM longitude array + wrap around point
    real(r8) dummy(plontgcm,plev,plattgcm)  ! dummy array
    real(r8) dummi(plattgcm,plevtgcm)       ! dummy array
    integer  plnout(plat)                   ! contains no. of output longitudes for eacj latitude

    real(r8) :: xccm_mx
    integer :: ktgcm_bt
    real(r8) :: ptgcm(plevtgcm)
    real(r8) :: qsotop(plontgcm,1,plattgcm)


!------------------Allocatable space (because nt is undetermined)--------------------------

    real(r8), allocatable, dimension(:,:,:,:) :: o1tgcmt  ! O1 (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: o2tgcmt  ! O2 (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: o3tgcmt  ! O3 (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: n2tgcmt  ! N2 (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: co2tgcmt ! CO2 (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: notgcmt  ! NO (including wraparound point)
    real(r8), allocatable, dimension(:,:,:,:) :: qtgcmt   ! QSOLAR (including wraparound point)
    real(r8), allocatable, dimension(:,:,:) :: tnt   ! TN (including wraparound point)
    real(r8), allocatable, dimension(:,:,:) :: zaltt   ! ZALTT (including wraparound point)

    real(r8), allocatable, dimension(:,:,:,:) :: o1tgcm   ! O1
    real(r8), allocatable, dimension(:,:,:,:) :: o2tgcm   ! O2
    real(r8), allocatable, dimension(:,:,:,:) :: o3tgcm   ! O3
    real(r8), allocatable, dimension(:,:,:,:) :: n2tgcm   ! N2
    real(r8), allocatable, dimension(:,:,:,:) :: co2tgcm  ! CO2
    real(r8), allocatable, dimension(:,:,:,:) :: notgcm   ! NO
    real(r8), allocatable, dimension(:,:,:,:) :: qtgcm    ! QSOLAR
    real(r8), allocatable, dimension(:,:,:) :: tn    ! TN
    real(r8), allocatable, dimension(:,:,:) :: zalt    ! ZALT

    integer, allocatable, dimension(:,:)     ::  mtime

!----------------------------------------------------------------------------------------

    itgcmcyc    = itgcmcyc_in
    cftgcm      = cftgcm_in

    if (masterproc) then

       write(iulog,*) 'TGCM: Flag for cycling TIME/GCM init fields is',itgcmcyc
       write(iulog,*) 'TGCM: MSS Path name for constituents field', cftgcm

!     Open netCDF file

       write(iulog,*) 'TGCM_INIT: Opening local file '
       call getfil(cftgcm,  locfn, 0)
       write(iulog,*) 'Opening netCDF file ',locfn
       call wrap_open (locfn, 0, ncid_tgcm)

!     Get dimensions id
       write(iulog,*) 'TGCM_INIT: Get dimensions ID'
       call wrap_inq_dimid( ncid_tgcm, 'lon', londimid   )
       call wrap_inq_dimid( ncid_tgcm, 'lev', levdimid   )
       call wrap_inq_dimid( ncid_tgcm, 'lat', latdimid   )
       call wrap_inq_dimid( ncid_tgcm, 'time', tdimid  )
       write(iulog,*) 'TGCM_INIT: lon, lat, lev, time dimIDs:', londimid, latdimid, levdimid, tdimid

!     Get dimensions sizes
       write(iulog,*) 'TGCM_INIT: Get dimensions sizes'
       call wrap_inq_dimlen( ncid_tgcm, londimid, lonsiz   )
       call wrap_inq_dimlen( ncid_tgcm, levdimid, levsiz   )
       call wrap_inq_dimlen( ncid_tgcm, latdimid, latsiz   )
       call wrap_inq_dimlen( ncid_tgcm, tdimid, timesiz   )
       write(iulog,*) 'TGCM_INIT: lon, lat, lev, time dimSizes:', lonsiz, latsiz, levsiz, timesiz

!     Check dimensions sizes are consistent with parameters definition
       write(iulog,*) 'TGCM_INIT: Chck dims are consistent with params def'
       if ((lonsiz-1).ne.plontgcm .or. latsiz.ne.plattgcm .or.  levsiz.ne.plevtgcm) then
          write(iulog,*) '***TGCM_INIT: dimensions sizes do not match***'
          write(iulog,*) 'lonsiz=',lonsiz,'plontgcm=',plontgcm
          write(iulog,*) 'latsiz=',latsiz,'plattgcm=',plattgcm
          write(iulog,*) 'levsiz=',levsiz,'plevtgcm=',plevtgcm
          call endrun
       endif
!
!     Check time dimension is consistent with user provided ITGCMCYC
!     ITGCMCYC == 0 constant input TIME/GCM fields (no lin. interp. required)
!     ITGCMCYC == 1 cycles 12 samples thru a year (lin. interpol. required)
!     ITGCMCYC == 2 only 2 time samples of TIME/GCM fields (lin. interp. between the two is required)
!
       if (itgcmcyc.eq.1) then

          nt=12
          write(iulog,*) 'TGCM_INIT: Found ',timesiz,' time samples'
          if (timesiz.ne.nt) then
             write(iulog,*) 'TGCM_INIT: timesiz=',timesiz,'but expected 12 smpls'
             call endrun
          endif

       elseif (itgcmcyc.eq.0) then

          nt=1
          write(iulog,*) 'TGCM_INIT: Found ',timesiz,' time samples'
          if (timesiz.ne.nt) then
             write(iulog,*) 'TGCM_INIT: timesiz=',timesiz,'but expected 1 smpl'
             call endrun
          endif

       elseif (itgcmcyc.eq.2) then

          nt=2
          write(iulog,*) 'TGCM_INIT: Found ',timesiz,' time samples'
          if (timesiz.ne.nt) then
             write(iulog,*) 'TGCM_INIT: timesiz=',timesiz,'but expected 2 smpls'
             call endrun
          endif

       endif

!     Get variables id
       write(iulog,*) 'TGCM_INIT: Get variables ID'
       call wrap_inq_varid( ncid_tgcm, 'MTIME', mtimeid   )
       if (.not.nlte_use_mo) then
          call wrap_inq_varid( ncid_tgcm, 'O1', o1id   )
          call wrap_inq_varid( ncid_tgcm, 'O2', o2id   )
          call wrap_inq_varid( ncid_tgcm, 'N2', n2id   )
          call wrap_inq_varid( ncid_tgcm, 'CO2', co2id   )
          call wrap_inq_varid( ncid_tgcm, 'NO', noid   )
       endif
       if (use_data_o3) then
          call wrap_inq_varid( ncid_tgcm, 'O3', o3id   )
       endif
       call wrap_inq_varid( ncid_tgcm, 'QSOLAR', qsolarid   )
       call wrap_inq_varid( ncid_tgcm, 'TTOP', tnid   )
       call wrap_inq_varid( ncid_tgcm, 'ZALT', zaltid   )
       call wrap_inq_varid( ncid_tgcm, 'lon', lonid   )
       call wrap_inq_varid( ncid_tgcm, 'lat', latid   )
       call wrap_inq_varid( ncid_tgcm, 'lev', levid   )

!     Check that dimesions are ordered properly
       if (.not.nlte_use_mo) then
          call wrap_inq_vardimid (ncid_tgcm, o1id, dimids)
          write(iulog,*) 'TGCM_INIT: check dimension order'
          if (dimids(1).ne.londimid .or. dimids(2).ne.latdimid .or. dimids(3).ne.levdimid) then
             write(iulog,*)'TGCM_INIT: 1. Data not in expected order'
             write(iulog,*)'dimids:', dimids(1), dimids(2), dimids(3)
             write(iulog,*)'lonid, latid, levid:', londimid, latdimid, levdimid
             call endrun
          end if
       else
          call wrap_inq_vardimid (ncid_tgcm, tnid, dimids(:2))
          write(iulog,*) 'TGCM_INIT: check dimension order'
          if (dimids(1).ne.londimid .or. dimids(2).ne.latdimid) then
             write(iulog,*)'TGCM_INIT: 2. Data not in expected order'
             write(iulog,*)'dimids:', dimids(1), dimids(2)
             write(iulog,*)'lonid, latid:', londimid, latdimid
             call endrun
          end if
       end if
       if (use_data_o3) then
          call wrap_inq_vardimid (ncid_tgcm, o3id, dimids)
          write(iulog,*) 'TGCM_INIT: check dimension order'
          if (dimids(1).ne.londimid .or. dimids(2).ne.latdimid .or. dimids(3).ne.levdimid) then
             write(iulog,*)'TGCM_INIT: 3. Data not in expected order'
             write(iulog,*)'dimids:', dimids(1), dimids(2), dimids(3)
             write(iulog,*)'lonid, latid, levid:', londimid, latdimid, levdimid
             call endrun
          end if
       endif

!
!     Retrieve longitude, latitude and level arrays for subsequent interpolation
!
       write(iulog,*) 'TGCM_INIT: retreive long, lat, level arrays'
       call wrap_get_var_realx (ncid_tgcm, lonid  ,lontgcmt)
       call wrap_get_var_realx (ncid_tgcm, latid  ,lattgcm)
       call wrap_get_var_realx (ncid_tgcm, levid  ,xtgcm)

!     TIME/GCM longitude runs in [-180,+180] (or, [180W,180E]). Remap in to [0,360]
       write(iulog,*) 'TGCM_INIT: re-arrange longitude array'
       lontgcm=cshift(lontgcmt(1:plontgcm),-36,1)
       lontgcm(plontgcm/2+1:plontgcm)=lontgcm(plontgcm/2+1:plontgcm)+360._r8

!     MTIME contains the model (day,hour,minute) triplet: convert to julian date+fraction
       write(iulog,*) 'TGCM_INIT: convert mtime to julian date+fract'
       strt2(1)=1
       strt2(2)=1
       cnt2(1)=3
       cnt2(2)=nt
       allocate (mtime(3,nt))
       call wrap_get_vara_int (ncid_tgcm, mtimeid,strt2,cnt2,mtime)
       caldtgcm(1:ptime)=0.0_r8
       do j=1,nt
          caldtgcm(j)=mtime(1,j)+mtime(2,j)/24._r8+mtime(3,j)/(60._r8*24._r8)
       enddo
       deallocate (mtime)

       strt4(1) = 1
       strt4(2) = 1
       strt4(3) = 1
       strt4(4) = 1
       cnt4(1)  = lonsiz
       cnt4(2)  = latsiz
       cnt4(3)  = levsiz
       cnt4(4)  = nt

       strt3(1) = 1
       strt3(2) = 1
       strt3(3) = 1
       cnt3(1)  = lonsiz
       cnt3(2)  = latsiz
       cnt3(3)  = nt

!     Allocate arrays
       write(iulog,*) 'TGCM_INIT: allocate arrays:'
       write(iulog,*) 'plontgcm=',plontgcm
       write(iulog,*) 'plattgcm=',plattgcm
       write(iulog,*) 'plevtgcm=',plevtgcm
       write(iulog,*) 'nt=',nt

       if (.not.nlte_use_mo) then

          allocate (o1tgcmt    (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'CHETGCMI: memory allocation problem (o1tgcmt)'
             call endrun
          endif

          allocate (o2tgcmt    (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (o2tgcmt)'
             call endrun
          endif

          allocate (n2tgcmt    (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (n2tgcmt)'
             call endrun
          endif

          allocate (co2tgcmt   (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (co2tgcmt)'
             call endrun
          endif

          allocate (notgcmt    (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (notgcmt)'
             call endrun
          endif

       end if
       if (use_data_o3) then
          allocate (o3tgcmt    (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (o3tgcmt)'
             call endrun
          endif
       endif

       allocate (qtgcmt     (plontgcm+1,plattgcm,plevtgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***CHETGCMI: memory allocation problem (qtgcmt)'
          call endrun
       endif

       allocate (tnt     (plontgcm+1,plattgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***CHETGCMI: memory allocation problem (tnt)'
          call endrun
       endif

       allocate (zaltt     (plontgcm+1,plattgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***CHETGCMI: memory allocation problem (zaltt)'
          call endrun
       endif

!     Retrieve the data
       write(iulog,*) 'TGCM_INIT: retreive data'
       if (.not.nlte_use_mo) then
          call wrap_get_vara_realx (ncid_tgcm,o1id    ,strt4,cnt4,o1tgcmt)
          call wrap_get_vara_realx (ncid_tgcm,o2id    ,strt4,cnt4,o2tgcmt)
          call wrap_get_vara_realx (ncid_tgcm,n2id    ,strt4,cnt4,n2tgcmt)
          call wrap_get_vara_realx (ncid_tgcm,co2id   ,strt4,cnt4,co2tgcmt)
          call wrap_get_vara_realx (ncid_tgcm,noid    ,strt4,cnt4,notgcmt)
       endif
       if (use_data_o3) then
          call wrap_get_vara_realx (ncid_tgcm,o3id    ,strt4,cnt4,o3tgcmt)
          call wrap_get_vara_realx (ncid_tgcm,o3id    ,strt4,cnt4,o3tgcmt)
       endif
       call wrap_get_vara_realx (ncid_tgcm,qsolarid,strt4,cnt4,qtgcmt)
       call wrap_get_vara_realx (ncid_tgcm,tnid    ,strt3,cnt3,tnt)
       call wrap_get_vara_realx (ncid_tgcm,zaltid    ,strt3,cnt3,zaltt)


!    Close netCDF file
       call wrap_close(ncid_tgcm)

!     Remap longitude axis for each input field and store
       write(iulog,*) 'TGCM_INIT: remap longitude arrays'

       if (.not.nlte_use_mo) then

          allocate (o1tgcm    (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'CHETGCMI: memory allocation problem (o1tgcm)'
             call endrun
          endif

          allocate (o2tgcm    (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (o2tgcm)'
             call endrun
          endif

          allocate (n2tgcm    (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (n2tgcm)'
             call endrun
          endif

          allocate (co2tgcm   (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (co2tgcm)'
             call endrun
          endif

          allocate (notgcm   (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (notgcm)'
             call endrun
          endif

          o1tgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                     &
               cshift(o1tgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          o2tgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                     &
               cshift(o2tgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          n2tgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                     &
               cshift(n2tgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          co2tgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                    &
               cshift(co2tgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          notgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                    &
               cshift(notgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          deallocate (o1tgcmt,o2tgcmt,n2tgcmt,co2tgcmt,notgcmt, stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory deallocation problem'
             write(iulog,*) '***(o1tgcmt,o2tgcmt,n2tgcmt,co2tgcmt)'
             call endrun
          endif

       end if

       if (use_data_o3) then
          allocate (o3tgcm    (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory allocation problem (o3tgcm)'
             call endrun
          endif

          o3tgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                     &
               cshift(o3tgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

          deallocate (o3tgcmt, stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***CHETGCMI: memory deallocation problem'
             write(iulog,*) '***(o3tgcmt)'
             call endrun
          endif
       endif

       allocate (qtgcm     (plontgcm,plattgcm,plevtgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: memory allocation problem (qtgcm)'
          call endrun
       endif

       allocate (tn     (plontgcm,plattgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: memory allocation problem (tn)'
          call endrun
       endif

       allocate (zalt     (plontgcm,plattgcm,nt),stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: memory allocation problem (zalt)'
          call endrun
       endif

       qtgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                      &
            cshift(qtgcmt(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt),-36,1)

       tn(1:plontgcm,1:plattgcm,1:nt)=                      &
            cshift(tnt(1:plontgcm,1:plattgcm,1:nt),-36,1)

       zalt(1:plontgcm,1:plattgcm,1:nt)=                      &
            cshift(zaltt(1:plontgcm,1:plattgcm,1:nt),-36,1)


       deallocate (qtgcmt,tnt,zaltt,stat=ierr)

       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: memory deallocation problem'
          write(iulog,*) '***(qtgcmt,tnt,zalt)'
          call endrun
       endif

!     transform QSOLAR from K/day to K/sec
       qtgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)=                      &
            qtgcm(1:plontgcm,1:plattgcm,1:plevtgcm,1:nt)/86400._r8

!      
!     Perform linear interpolation to CCM grid
!     Vertical linear interpolation is performed in pressure scale height (psh):
!         
!            X(TGCM,CCM)=LOG(Pref/P)
!
!     where Pref is 5E-7 mbar for TGCM and 1000 mbar for CCM
!
!     It is necessary to transform XTGCM so that it maps as XCCM:
!
!            XTGCM_CCM=XTGCM-LOG(Pref(TGCM)/Pref(CCM))

       do k=1,plevtgcm
          xtgcm_ccm(k)=xtgcm(k)-log(preftgcm/1E5_r8)
          ptgcm(k) = 1e5_r8 * exp (-xtgcm_ccm(k))

          write(iulog,*) xtgcm_ccm(k),ptgcm(k)

       enddo

!     calcualte XCCM 
!     Note that the indexing for xccm is reversed with respect to the CCM grid

       do k=1,plev
          kinv=plev-k+1
          xccm(k)=log(1e5_r8/hypm(kinv))
       enddo

       write(iulog,*) 'TGCM_INIT: space (3D)  interpolation'

       plnout(:)=plon


       if (.not.nlte_use_mo) then

!
!     O1
!
          allocate (o1(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'TGCM_INIT: memory allocation problem (o1)'
             call endrun
          endif

          o1(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm) = o1tgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,            &
                  plev,plev,plattgcm,o1(1,1,1,j),londeg,                      &
                  latdeg,plon,plnout,plev,plat)


          enddo                ! Another time

          write(iulog,*) 'O1: done'
!     deallocate space 
          deallocate (o1tgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of O1TGCM failed***'
             call endrun
          endif

!
!     O2
!
          allocate (o2(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'CHETGCMT: memory allocation problem (o2)'
             call endrun
          endif

          o2(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm) = o2tgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,                    &
                  plev,plev,plattgcm,o2(1,1,1,j),londeg,                              &
                  latdeg,plon,plnout,plev,plat)

          enddo                ! Another time

          write(iulog,*) 'O2: done'
!     deallocate space 
          deallocate (o2tgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of O2TGCM failed***'
             call endrun
          endif

!
!     N2
!
          allocate (n2(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'TGCM_INIT: memory allocation problem (n2)'
             call endrun
          endif

          n2(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm)=n2tgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,       &
                  plev,plev,plattgcm,n2(1,1,1,j),londeg,                 &
                  latdeg,plon,plnout,plev,plat)


          enddo                ! Another time

          write(iulog,*) 'N2: done'
!     deallocate space 
          deallocate (n2tgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of N2TGCM failed***'
             call endrun
          endif

!
!     CO2
!
          allocate (co2(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'TGCM_INIT: memory allocation problem (co2)'
             call endrun
          endif

          co2(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm)=co2tgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,                  &
                  plev,plev,plattgcm,co2(1,1,1,j),londeg,                           &
                  latdeg,plon,plnout,plev,plat)

          enddo                ! Another time

          write(iulog,*) 'CO2: done'
!     deallocate space 
          deallocate (co2tgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of CO2TGCM failed***'
             call endrun
          endif

!
!     NO
!
          allocate (no(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'TGCM_INIT: memory allocation problem (NO)'
             call endrun
          endif

          no(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm)=notgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,                  &
                  plev,plev,plattgcm,no(1,1,1,j),londeg,                           &
                  latdeg,plon,plnout,plev,plat)

          enddo                ! Another time

          write(iulog,*) 'NO: done'
!     deallocate space 
          deallocate (notgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of NOTGCM failed***'
             call endrun
          endif


       endif

       if (use_data_o3) then
!
!     O3
!
          allocate (o3(plon,plev,plat,nt), stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) 'TGCM_INIT: memory allocation problem (o3)'
             call endrun
          endif

          o3(1:plon,1:plev,1:plat,1:nt)=0.0_r8
          do j=1,nt

!     vertical interpolation 
             do i=1,plontgcm

!     collapse indices
                dummi(1:plattgcm,1:plevtgcm)=o3tgcm(i,1:plattgcm,1:plevtgcm,j)
                call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
                dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

             enddo             ! Another longitude

!     longitude/latitude interpolation
             call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,     &
                  plev,plev,plattgcm,o3(1,1,1,j),londeg,               &
                  latdeg,plon,plnout,plev,plat)

          enddo                ! Another time

          write(iulog,*) 'O3: done'
!     deallocate space 
          deallocate (o3tgcm,stat=ierr)
          if (ierr.ne.0) then
             write(iulog,*) '***TGCM_INIT: deallocation of O3TGCM failed***'
             call endrun
          endif

       end if

!
!     QSOLAR
!
       allocate (qsolar(plon,plev,plat,nt), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (qsolar)'
          call endrun
       endif

       qsolar(1:plon,1:plev,1:plat,1:ptime)=0.0_r8
       do j=1,nt

!     vertical interpolation 
          do i=1,plontgcm

!     collapse indices
             dummi(1:plattgcm,1:plevtgcm)=qtgcm(i,1:plattgcm,1:plevtgcm,j)
             call lininterp (dummi,xtgcm_ccm,plattgcm,plevtgcm,dummo,xccm,plev)
             dummy(i,1:plev,1:plattgcm)=transpose(dummo(1:plattgcm,plev:1:-1))

          enddo             ! Another longitude

!     longitude/latitude interpolation
          call bilin (dummy,lontgcm,lattgcm,plontgcm,plontgcm,                &
               plev,plev,plattgcm,qsolar(1,1,1,j),londeg,                      &
               latdeg,plon,plnout,plev,plat)


       enddo        ! Another time

       write(iulog,*) 'QSOLAR: done'
!     deallocate space 
       deallocate (qtgcm,stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: deallocation of QTGCM failed***'
          call endrun
       endif

!
!     TN
!
       t_top_msis = 0.0_r8
       do j=1, nt

          qsotop(:plontgcm,1,:plattgcm)=tn(:plontgcm,:plattgcm,j)

          call bilin (qsotop,lontgcm,lattgcm,plontgcm,plontgcm,    &
               1,1,plattgcm,t_top_msis(1,1,1,j),londeg,        &
               latdeg,plon,plnout,1,plat)


          if (j == 1) then
             write(iulog,*) t_top_msis(:plon,1,5,j)
          endif

       enddo


       write(iulog,*) 'TTOP: done'
!     deallocate space 
       deallocate (tn,stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: deallocation of TTOP failed***'
          call endrun
       endif


!
!     ZALT
!
       z_top_msis = 0.0_r8
       do j=1, nt

          qsotop(:plontgcm,1,:plattgcm)=zalt(:plontgcm,:plattgcm,j)

          call bilin (qsotop,lontgcm,lattgcm,plontgcm,plontgcm,    &
               1,1,plattgcm,z_top_msis(1,1,1,j),londeg,        &
               latdeg,plon,plnout,1,plat)


          if (j == 1) then
             write(iulog,*) z_top_msis(:plon,1,5,j)
          endif

       enddo


       write(iulog,*) 'ZALT: done'
!     deallocate space 
       deallocate (zalt,stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) '***TGCM_INIT: deallocation of ZALT failed***'
          call endrun
       endif

       write(iulog,*) 'TGCM_INIT: Successfully completed'


    endif            ! end masterproc

!  allocate memory 
    if (.not.nlte_use_mo) then

!++bee memory that's referenced on all procs must be allocated on all procs.
!      if (masterproc) then

       allocate (o1c (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (o1c)'
          call endrun
       endif

       allocate (o2c (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (o2c)'
          call endrun
       endif

       allocate (n2c (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (n2c)'
          call endrun
       endif

       allocate (co2c (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (co2c)'
          call endrun
       endif

       allocate (noc (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (noc)'
          call endrun
       endif

!      end if
!--bee

       allocate ( o1_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for O1_TO_CHUNK'
          call endrun
       end if

       allocate ( o2_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for O2_TO_CHUNK'
          call endrun
       end if

       allocate ( n2_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for N2_TO_CHUNK'
          call endrun
       end if

       allocate ( co2_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for CO2_TO_CHUNK'
          call endrun
       end if

       allocate ( no_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for NO_TO_CHUNK'
          call endrun
       end if

    end if

    if (use_data_o3) then

       allocate (o3c (plon,plev,plat), stat=ierr)
       if (ierr.ne.0) then
          write(iulog,*) 'TGCM_INIT: memory allocation problem (o3c)'
          call endrun
       endif

       allocate ( o3_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
       if (ierr /= 0) then
          write(iulog,*) 'TGCM_INIT: Error allocating space for O3_TO_CHUNK'
          call endrun
       end if

    end if


    allocate (qsolarc(plon,plev,plat), stat=ierr)
    if (ierr /= 0) then
       write(iulog,*) 'TGCM_INIT: memory allocation problem (qsolarc)'
       call endrun
    endif

    allocate (qsolar_to_chunk(pcols,plev,begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       write(iulog,*) 'TGCM_INIT: Error allocating space for QSOLAR_TO_CHUNK'
       call endrun
    end if

    allocate (ttop_to_chunk(pcols,1,begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       write(iulog,*) 'TGCM_INIT: Error allocating space for TTOP_TO_CHUNK'
       call endrun
    end if

    allocate (zalt_to_chunk(pcols,1,begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       write(iulog,*) 'TGCM_INIT: Error allocating space for ZALT_TO_CHUNK'
       call endrun
    end if



! end

    return
  end subroutine tgcm_init

!=======================================================================

  subroutine tgcm_timestep_init (nlte_use_mo, use_data_o3)
!
!     Time interpolation of TIME/GCM fields to the current CCM time
!
!------------------------------------------------------------------------

    use time_manager, only: get_curr_calday
    use phys_grid, only: scatter_field_to_chunk
#if ( defined SPMD )
    use mpishorthand
#endif

!------------------------------------------------------------------------
    logical,          intent(in) :: nlte_use_mo
    logical,          intent(in) :: use_data_o3

!---------------------------Local workspace--------------------------------------

    integer k                            ! index
    integer i                            ! index
    integer j                            ! index
    integer kk                           ! index
    integer kkm                          ! index
    integer kkp                          ! index
    integer nt                           ! no. of time samples
    integer ierr                         ! error status
    integer ixshf                        ! number of grid points or x-shift
    integer kinv

    real(r8) cdtgcmp                         ! adjacent (greater) TIME/GCM day to CALDAY
    real(r8) cdtgcmm                         ! adjacent (lower) TIME/GCM day to CALDAY
    real(r8) adtgcm                          ! no. of days between adjacent CALDTGCM()
    real(r8) dm                              ! coefficient for linear interpolation
    real(r8) dp                              ! coefficient for linear interpolation
    real(r8) cday                            ! dummy calendar day
    real(r8) cdayfr                          ! fraction of CALDAY
    real(r8) dummy(plon)                     ! dummy used in cshft (***LEAVE DIM TO PLON***)
    real(r8) :: calday                            ! current calendar day
    real(r8) blah

    if (masterproc) then

       calday = get_curr_calday()
#ifdef TGCM_DIAGS
       write(iulog,*) 'Interpolating waccm bnd data at',calday
#endif

       if (itgcmcyc.eq.0) then

!     cycling the same dataset

          nt=1

          if (.not. nlte_use_mo) then
             o1c (1:plon,1:plev,1:plat) = o1 (1:plon,1:plev,1:plat,nt)
             o2c (1:plon,1:plev,1:plat) = o2 (1:plon,1:plev,1:plat,nt)
             n2c (1:plon,1:plev,1:plat) = n2 (1:plon,1:plev,1:plat,nt)
             co2c(1:plon,1:plev,1:plat) = co2(1:plon,1:plev,1:plat,nt)
             noc (1:plon,1:plev,1:plat) = no (1:plon,1:plev,1:plat,nt)
          endif
          if (use_data_o3) then
             o3c (1:plon,1:plev,1:plat) = o3 (1:plon,1:plev,1:plat,nt)
          endif

          qsolarc(1:plon,1:plev,1:plat)=qsolar(1:plon,1:plev,1:plat,nt)

          t_top_msisc(1:plon,1:plat) = t_top_msis(1:plon,1,1:plat,nt)
          z_top_msisc(1:plon,1:plat) = z_top_msis(1:plon,1,1:plat,nt)

       elseif (itgcmcyc.eq.1) then

!     use 12 datasets
          nt=12

!     Find CALDTGCM() adjacent to CALDAY
          kk=0
          do i=1,nt-1
             kk=merge(i,kk,(calday-caldtgcm(i))*(caldtgcm(i+1)-calday)>=0)
          enddo

          if (kk.ne.0) then

             adtgcm=caldtgcm(kk+1)-caldtgcm(kk)
             dm=(caldtgcm(kk+1)-calday)/adtgcm
             dp=(calday-caldtgcm(kk))/adtgcm

             if (.not.nlte_use_mo) then

                o1c(1:plon,1:plev,1:plat)=                           &
                     o1(1:plon,1:plev,1:plat,kk)*dm                  &
                     +o1(1:plon,1:plev,1:plat,kk+1)*dp

                o2c(1:plon,1:plev,1:plat)=                           &
                     o2(1:plon,1:plev,1:plat,kk)*dm                  &
                     +o2(1:plon,1:plev,1:plat,kk+1)*dp

                n2c(1:plon,1:plev,1:plat)=                           &
                     n2(1:plon,1:plev,1:plat,kk)*dm                  &
                     +n2(1:plon,1:plev,1:plat,kk+1)*dp

                co2c(1:plon,1:plev,1:plat)=                          &
                     co2(1:plon,1:plev,1:plat,kk)*dm                 &
                     +co2(1:plon,1:plev,1:plat,kk+1)*dp

                noc(1:plon,1:plev,1:plat)=                          &
                     no(1:plon,1:plev,1:plat,kk)*dm                 &
                     +no(1:plon,1:plev,1:plat,kk+1)*dp

             endif
             if (use_data_o3) then
                o3c(1:plon,1:plev,1:plat)=                           &
                     o3(1:plon,1:plev,1:plat,kk)*dm                  &
                     +o3(1:plon,1:plev,1:plat,kk+1)*dp
             endif

             qsolarc(1:plon,1:plev,1:plat)=                       &
                  qsolar(1:plon,1:plev,1:plat,kk)*dm              &
                  +qsolar(1:plon,1:plev,1:plat,kk+1)*dp

             t_top_msisc(1:plon,1:plat)=                     &
                  t_top_msis(1:plon,1,1:plat,kk)*dm             &
                  +t_top_msis(1:plon,1,1:plat,kk+1)*dp

             z_top_msisc(1:plon,1:plat)=                     &
                  z_top_msis(1:plon,1,1:plat,kk)*dm             &
                  +z_top_msis(1:plon,1,1:plat,kk+1)*dp

          elseif (kk.eq.0) then

!     KK==0 means interpolation is between i=12(i.e.,december) and i=1(i.e.,january)
             kkm=12
             kkp=1
             cdtgcmp=caldtgcm(kkp)+365._r8
             cdtgcmm=caldtgcm(kkm)
             adtgcm=cdtgcmp-cdtgcmm

!     if CALDAY in january (i.e. .NOT.(335 <= CALDAY <= 366.)) add 365
             if (.not.(calday.ge.335._r8 .and. calday.le.366._r8)) then
!     CALDAY is in january
                cday=calday+365._r8
             else
!     CALDAY is in december
                cday=calday
             endif

             dm=(cdtgcmp-cday)/adtgcm
             dp=(cday-cdtgcmm)/adtgcm

             if (.not.nlte_use_mo) then

                o1c(1:plon,1:plev,1:plat)=                            &
                     o1(1:plon,1:plev,1:plat,kkm)*dm                  &
                     +o1(1:plon,1:plev,1:plat,kkp)*dp

                o2c(1:plon,1:plev,1:plat)=                            &
                     o2(1:plon,1:plev,1:plat,kkm)*dm                  &
                     +o2(1:plon,1:plev,1:plat,kkp)*dp

                n2c(1:plon,1:plev,1:plat)=                            &
                     n2(1:plon,1:plev,1:plat,kkm)*dm                  &
                     +n2(1:plon,1:plev,1:plat,kkp)*dp

                co2c(1:plon,1:plev,1:plat)=                           &
                     co2(1:plon,1:plev,1:plat,kkm)*dm                 &
                     +co2(1:plon,1:plev,1:plat,kkp)*dp

                noc(1:plon,1:plev,1:plat)=                           &
                     no(1:plon,1:plev,1:plat,kkm)*dm                 &
                     +no(1:plon,1:plev,1:plat,kkp)*dp

             endif
             if (use_data_o3) then
                o3c(1:plon,1:plev,1:plat)=                            &
                     o3(1:plon,1:plev,1:plat,kkm)*dm                  &
                     +o3(1:plon,1:plev,1:plat,kkp)*dp
             endif

             qsolarc(1:plon,1:plev,1:plat)=                        &
                  qsolar(1:plon,1:plev,1:plat,kkm)*dm              &
                  +qsolar(1:plon,1:plev,1:plat,kkp)*dp

             t_top_msisc(1:plon,1:plat)=                     &
                  t_top_msis(1:plon,1,1:plat,kkm)*dm             &
                  +t_top_msis(1:plon,1,1:plat,kkp)*dp

             z_top_msisc(1:plon,1:plat)=                     &
                  z_top_msis(1:plon,1,1:plat,kkm)*dm             &
                  +z_top_msis(1:plon,1,1:plat,kkp)*dp

          endif                         ! ends logic KK=0

       endif                            ! ends logic ITGCMCYC

!
!     Adjust for local time. TIME/GCM output is provided 
!     at 0Z, therefore it is necessary to x-shift the constituents + QSOLAR
!     in order to allign with the local time of CCM at 0E
!     The local time difference DHLC at 0E between TIME/GCM and CCM is:
!
!                DHLC = CDAYFR *360
!                            
!     where CDAYFR is the fraction of CALDAY
!     Then the shift of x-grid points  is determined by
!
!                                PLON     
!                IXSHF = DHLC * ------ =  CDAYFR * PLON
!                                360       
!
!     Note: CSHIFT perform a circular shift. The constituents fields
!     are dimesioned to PLOND. DO NOT DO A CIRCULAR SHIFT ON THE ENTIRE 
!     PLOND DIMENSION OTHERWISE ZEROS (OR OTHER MEANINGLESS NUMBERS)
!     WILL APPEAR IN THE ARRAY
!

       cdayfr=calday-int(calday)

       ixshf=int(anint(cdayfr*real(plon,r8)))

! DO NOT SWITCH ORDER OF THE LOOP INDEXES!!!
       do j=1,plat
          do k=1,plev

             if (.not.nlte_use_mo) then

                dummy(1:plon)=o1c(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                o1c(1:plon,k,j)=dummy(1:plon)

                dummy(1:plon)=o2c(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                o2c(1:plon,k,j)=dummy(1:plon)

                dummy(1:plon)=n2c(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                n2c(1:plon,k,j)=dummy(1:plon)

                dummy(1:plon)=co2c(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                co2c(1:plon,k,j)=dummy(1:plon)

                dummy(1:plon)=noc(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                noc(1:plon,k,j)=dummy(1:plon)

             endif
             if (use_data_o3) then
                dummy(1:plon)=o3c(1:plon,k,j)
                dummy=cshift(dummy,ixshf,1)
                o3c(1:plon,k,j)=dummy(1:plon)
             endif

             dummy(1:plon)=qsolarc(1:plon,k,j)
             dummy=cshift(dummy,ixshf,1)
             qsolarc(1:plon,k,j)=dummy(1:plon)

          enddo

          dummy(1:plon)=t_top_msisc(1:plon,j)
          dummy=cshift(dummy,ixshf,1)
          t_top_msisc(1:plon,j)=dummy(1:plon)

          dummy(1:plon)=z_top_msisc(1:plon,j)
          dummy=cshift(dummy,ixshf,1)
          z_top_msisc(1:plon,j)=dummy(1:plon)

       enddo

    endif                  ! end masterproc


! Scatter results
    if (.not.nlte_use_mo) then
       call scatter_field_to_chunk (1,plev,1,plon,o1c ,o1_to_chunk(1,1,begchunk))
       call scatter_field_to_chunk (1,plev,1,plon,o2c ,o2_to_chunk(1,1,begchunk))
       call scatter_field_to_chunk (1,plev,1,plon,n2c ,n2_to_chunk(1,1,begchunk)) 
       call scatter_field_to_chunk (1,plev,1,plon,co2c,co2_to_chunk(1,1,begchunk))  
       call scatter_field_to_chunk (1,plev,1,plon,noc ,no_to_chunk(1,1,begchunk)) 
    endif
    if (use_data_o3) then
       call scatter_field_to_chunk (1,plev,1,plon,o3c ,o3_to_chunk(1,1,begchunk)) 
    endif
    call scatter_field_to_chunk (1,plev,1,plon,qsolarc    ,qsolar_to_chunk(1,1,begchunk)) 
    call scatter_field_to_chunk (1,1,1,plon,t_top_msisc,ttop_to_chunk(1,1,begchunk))
    call scatter_field_to_chunk (1,1,1,plon,z_top_msisc,zalt_to_chunk(1,1,begchunk)) 


    return
  end subroutine tgcm_timestep_init

!================================================================================================

  subroutine tgcm_get_o3 (ncol, lchnk, o3)
!
! Get ozone mass mixing ratios specified from input dataset
!-------------------------------------------------------------------------

! Arguments
    integer,  intent(in)  :: ncol                   ! no. of columns in chunk
    integer,  intent(in)  :: lchnk                  ! chunk identifier
    real(r8), pointer, dimension(:,:) :: o3(:,:)    ! tgcm ozone

! Local workspace
    integer :: k

!------------------------------------------------------------------------

    o3 => o3_to_chunk(:,:,lchnk)

  end subroutine tgcm_get_o3

!================================================================================================

  subroutine tgcm_get_cnst (ncol, lchnk, co2, o1, o2, no, n2)
!
! Get ozone mass mixing ratios specified from input dataset
!-------------------------------------------------------------------------

! Arguments
    integer,  intent(in)  :: ncol                   ! no. of columns in chunk
    integer,  intent(in)  :: lchnk                  ! chunk identifier
    real(r8), pointer, dimension(:,:) :: co2
    real(r8), pointer, dimension(:,:) :: o1
    real(r8), pointer, dimension(:,:) :: o2
    real(r8), pointer, dimension(:,:) :: no
    real(r8), pointer, dimension(:,:) :: n2

! Local workspace
    integer :: k

!------------------------------------------------------------------------

    do k = 1,pver
       co2 => co2_to_chunk(:,:,lchnk)
       o1  => o1_to_chunk (:,:,lchnk)
       o2  => o2_to_chunk (:,:,lchnk)
       no  => no_to_chunk (:,:,lchnk)
       n2  => n2_to_chunk (:,:,lchnk)
    end do

  end subroutine tgcm_get_cnst

!================================================================================================

  subroutine tgcm_get_solar (ncol, lchnk, qrs_mlt)
!
! Get M/LT solar heating rates specified from input dataset
!-------------------------------------------------------------------------

! Arguments
    integer,  intent(in)  :: ncol                   ! no. of columns in chunk
    integer,  intent(in)  :: lchnk                  ! chunk identifier
    real(r8), intent(out) :: qrs_mlt(pcols,pver)    ! M/LT solar heating rates

! Local workspace
    integer :: k

!------------------------------------------------------------------------

    do k = 1,pver
       qrs_mlt (:ncol,k) = qsolar_to_chunk(:ncol,k,lchnk)
    end do

  end subroutine tgcm_get_solar

end module tgcm_forcing
