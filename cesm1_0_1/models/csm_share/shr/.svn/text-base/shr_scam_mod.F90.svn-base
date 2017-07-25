!===============================================================================
! SVN $Id: 
! SVN $URL: 
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_scam_mod.F90 --- Module to handle single column mode share routines.
!
! !DESCRIPTION:
!    Routines needed by drv or several component models for running in single column mode
! 
! !REVISION HISTORY:
!    2007 Sep 14 - B. Kauffman  - svn checkin 
!    2007 Aug 29 - J. Truesdale - first version
!
! !INTERFACE: ------------------------------------------------------------------

module shr_scam_mod

! !USES:

   use shr_kind_mod  ! defines kinds
   use shr_sys_mod   ! system calls
   use shr_file_mod  ! file utilities
   use shr_kind_mod, only : R8=>SHR_KIND_R8,IN=>SHR_KIND_IN,CL=>SHR_KIND_CL
   use shr_log_mod,  only : s_loglev  => shr_log_Level
   use shr_log_mod,  only : s_logunit => shr_log_Unit
   use netcdf

   implicit none

   private           ! By default everything is private to this module

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_scam_getCloseLatLon ! return lat and lon point/index
   public :: shr_scam_checkSurface   ! check grid fraction in focndomain dataset

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP
!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_scam_getCloseLatLon(ncid,targetLat,targetLon, closeLat, closeLon, 
!                                    closeLatIdx, closeLonIdx)
!
! !DESCRIPTION:
!    routine to search in netcdf file and return lat and lon point/index closest to target point
!
! USAGE:
!    call shr_scam_getCloseLatLon(ncid,targetLat,targetLon, closeLat, closeLon, 
!                                 closeLatIdx, closeLonIdx)
!
! !REVISION HISTORY:
!     2007 Aug 29 - J. Truesdale - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_scam_getCloseLatLon(ncid, targetLat,  targetLon, closeLat, closeLon, &
                                    closeLatIdx, closeLonIdx)
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)         :: ncid  ! netcdf id
   real   (R8),intent(in)         :: targetLat  ! find closest latitude to this point
   real   (R8),intent(in)         :: targetLon  ! find closest longitude to this point
   real   (R8),intent(out)        :: closeLat   ! returned close lat
   real   (R8),intent(out)        :: closeLon   ! returned close lon
   integer(IN),intent(out)        :: closeLatIdx  ! index of returned lat point
   integer(IN),intent(out)        :: closeLonIdx  ! index of returned lon point

!EOP

   !----- local variables -----
   real   (R8),allocatable          :: lats(:),lons(:),poslons(:)
   real   (R8)                      :: postargetlon
   integer(IN)                      :: rcode   ! netCDF routine return code
   integer(IN)                      ::  i
   integer(IN)                      ::  len
   integer(IN)                      ::  ndims
   integer(IN)                      ::  nvars
   integer(IN)                      ::  nvarid
   integer(IN)                      ::  ndimid
   integer(IN)                      ::  strt(nf90_max_var_dims),count(nf90_max_var_dims)
   integer(IN)                      ::  nlon,nlat
   integer(IN), dimension(nf90_max_var_dims) :: dimids
   character(len=80)                ::  name,var_name
   character(*),parameter :: subname = "(shr_scam_getCloseLatLon) "

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   !--- Get variable info for search ---

   rcode = nf90_inquire(ncid, nVariables=nvars)
   if (rcode /= nf90_noerr) then
      call shr_sys_abort(subname//"ERROR from nf90_inquire ")
   endif

   !--- Look for/extract lat lon coordinate variables from file ---

   nlat=0
   nlon=0
   nvarid=0

   !--- Loop through all Lat variables until we find lat and lon ---

   do while (nvarid < nvars .and.(nlon.eq.0 .or. nlat.eq.0))
      nvarid=nvarid+1
      rcode = nf90_inquire_variable(ncid, nvarid, var_name, ndims=ndims,dimids = dimids)
      if (rcode /= nf90_noerr) then
         call shr_sys_abort(subname//"ERROR inquiring about variable "//trim(var_name))
      endif

      !--- is this a latitude variable ---

      if ( var_name .eq. 'lat'.or. var_name .eq. 'latixy'.or. var_name .eq. 'yc'.or.var_name.eq.'lsmlat'.or.&
           var_name .eq. 'LAT'.or. var_name .eq. 'LATIXY'.or. var_name .eq. 'YC'.or.var_name.eq.'LSMLAT' )  then

         !--- Loop through all variable dimensions until we find lat and lon ---

         do ndimid =  1,ndims
            rcode = nf90_inquire_dimension(ncid, dimids(ndimid), name, len)
            if (rcode /= nf90_noerr) then
               call shr_sys_abort(subname//"ERROR: Cant read netcdf latitude variable dimension")
            endif
            if ( name .eq. 'lat'.or. name .eq. 'latixy'.or. name .eq. 'nj'.or. name .eq. 'lsmlat' .or. &
                 name .eq. 'LAT'.or. name .eq. 'LATIXY'.or. name .eq. 'NJ'.or. name .eq. 'LSMLAT' )  then
               strt(ndimid) = 1
               count(ndimid) = len
               nlat=len
            else
               strt(ndimid) = 1
               count(ndimid) = 1
            endif
         end do
         if (nlat.eq.0) then
            call shr_sys_abort( subname//"ERROR: Cant find a useable latitude dimension (lat,nj,latixy, or lsmlat")
         endif
         allocate(lats(nlat))
         rcode= nf90_get_var(ncid, nvarid ,lats, start = strt, count = count)
         if (rcode /= nf90_noerr) then
            call shr_sys_abort( subname//"ERROR: Cant read netcdf latitude variable dimension")
         endif
      end if

      !--- is this a longitude variable ---

      if ( var_name .eq. 'lon'.or. var_name .eq. 'longxy'.or. var_name .eq. 'xc'.or.var_name.eq.'lsmlon'.or.&
           var_name .eq. 'LON'.or. var_name .eq. 'LONGXY'.or. var_name .eq. 'XC'.or.var_name.eq.'LSMLON' )  then
         do ndimid =  1,ndims
            rcode = nf90_inquire_dimension(ncid, dimids(ndimid), name, len)
            if (rcode /= nf90_noerr) then
               call shr_sys_abort( subname//"ERROR: Cant read netcdf latitude variable dimension")
            endif
            if ( name .eq. 'lon'.or. name .eq. 'longxy'.or. name .eq. 'ni'.or. name .eq. 'lsmlon' .or. &
                 name .eq. 'LON'.or. name .eq. 'LONGXY'.or. name .eq. 'NI'.or. name .eq. 'LSMLON' )  then
               strt(ndimid) = 1
               count(ndimid) = len
               nlon=len
            else
               strt(ndimid) = 1
               count(ndimid) = 1
            endif
         end do
         if (nlon.eq.0) then
            call shr_sys_abort( subname//"ERROR: Cant find a useable longitude dimension (lon,ni,longxy, or lsmlon")
         endif
         allocate(lons(nlon))
         allocate(poslons(nlon))
         rcode= nf90_get_var(ncid, nvarid ,lons, start = strt, count = count)
         if (rcode /= nf90_noerr) then
            call shr_sys_abort( subname//"ERROR: Cant read netcdf latitude variable dimension")
         endif
      end if
   end do

   !--- Did we get find valid lat and lon coordinate variables ---

   if (nlon.eq.0) then
      call shr_sys_abort( subname//"ERROR: Couldnt find a longitude coordinate variable")
   end if
   if (nlat.eq.0) then
      call shr_sys_abort( subname//"ERROR: Couldnt find a latitude coordinate variable")
   end if

   !--- convert lons array and targetlon to 0,360 ---

   poslons=mod(lons+360._r8,360._r8)
   postargetlon=mod(targetlon+360._r8,360._r8)

   !--- find index of value closest to 0 and set returned values ---

   closelonidx=(MINLOC(abs(poslons-postargetlon),dim=1))
   closelatidx=(MINLOC(abs(lats-targetlat),dim=1))
   closelon=lons(closelonidx)
   closelat=lats(closelatidx)

   !--- if it gets here we need to clean up after ourselves ---

   deallocate(lats)
   deallocate(lons)
   deallocate(poslons)

   return

end subroutine shr_scam_getCloseLatLon

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_scam_checkSurface(scmlon, scmlat, ocn_compid, ocn_mpicom, lnd_present, ocn_present, ice_present)
!
! !DESCRIPTION:
!    routine to check grid fraction from the focndomain dataset
!    and provide information to correctly flag land, ocean or ice for
!    single column mode
!
! !REVISION HISTORY:
!     2007 Aug 29 - J. Truesdale - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_scam_checkSurface(scmlon, scmlat, ocn_compid, ocn_mpicom, lnd_present, ocn_present, ice_present)

! !USES:
   use shr_strdata_mod
   use shr_dmodel_mod    ! shr data model stuff
   use mct_mod
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   real(R8),                     intent(in)  :: scmlon,scmlat ! single column lat lon
   integer(IN),                  intent(in)  :: ocn_compid   ! id for ocean model  
   integer(IN),                  intent(in)  :: ocn_mpicom   ! mpi communicator for ocean  
   logical,            optional, intent(out) :: lnd_present  ! land point
   logical,            optional, intent(out) :: ice_present  ! ice point
   logical,            optional, intent(out) :: ocn_present  ! ocean point

!EOP

   !----- local variables -----
   integer(IN)             :: rcode            ! error code
   integer(IN)             :: ncid_ocn         ! netcdf id for ocn_in
   integer(IN)             :: fracid           !  id for frac variable
   integer(IN)             :: closeLatIdx      ! index of returned lat point
   integer(IN)             :: closeLonIdx      ! index of returned lon point
   integer(IN)             :: unitn            ! io unit
   real   (R8)             :: ocn_frac(1,1)    ! ocean fraction
   real   (R8)             :: closeLat         ! returned close lat
   real   (R8)             :: closeLon         ! returned close lon
   character(len=CL)       :: nrevsn = ' '     ! full path restart file for branch
   character(len=CL)       :: rest_pfile = './rpointer.dom' ! restart pointer file
   character(len=CL)       :: bndtvs           ! sst file
   character(len=CL)       :: focndomain       ! ocn domain file
   logical                 :: sstcyc           ! flag for sst cycling
   logical                 :: docn_exists           ! flag if file exists locally
   logical                 :: ocn_exists            ! flag if file exists locally
   logical                 :: exists            ! flag if file exists locally
   logical                 :: aqua_planet      ! flags
   logical                 :: single_column    ! flags

   !----- formats -----
   character(*),parameter :: subname = "(shr_scam_checkSurface) "
   character(*),parameter :: F00   = "('(shr_scam_checkSurface) ',8a)" 
   type(shr_strdata_type) :: SDAT
   character(len=CL)      :: decomp = '1d' ! restart pointer file
   character(len=CL)      :: restfilm = 'unset' 
   character(len=CL)      :: restfils = 'unset'
   character(len=CL)      :: ocn_in = 'unset'
   integer(IN)   :: nfrac
   namelist /dom_inparm/ sstcyc, nrevsn, rest_pfile, bndtvs, focndomain
   namelist / docn_nml / ocn_in, decomp, restfilm, restfils

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   inquire( file='ocn_in', exist=ocn_exists )
   inquire( file='docn_in', exist=docn_exists )
   if (ocn_exists) then
      !--- read in the ocn_in namelist to get name for focndomain file

      unitn = shr_file_getUnit() ! get an unused unit number
      open( unitn, file='ocn_in', status='old' )
      rcode = 1
      do while ( rcode /= 0 )
         read(unitn, dom_inparm, iostat=rcode)
         if (rcode < 0) then
            call shr_sys_abort( 'shr_scam_checkSurface encountered end-of-file on namelist read' )
         endif
      end do
      close( unitn )
      call shr_file_freeUnit(unitn)

      !--- open the netcdf file ---

      inquire(file=trim(focndomain),exist=exists)
      if (.not.exists) call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(focndomain))
      rcode = nf90_open(focndomain,nf90_nowrite,ncid_ocn)
      if (rCode /= nf90_noerr) call shr_sys_abort(subName//"ERROR opening data file : "//trim(focndomain))
      if (s_loglev > 0) write(s_logunit,F00) 'opened netCDF data file: ',trim(focndomain)

      !--- Extract the fraction for current column ---
   
      call shr_scam_getCloseLatLon(ncid_ocn,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
      rcode = nf90_inq_varid(ncid_ocn, 'frac', fracid)
      if (rcode /= nf90_noerr) then
         call shr_sys_abort(subname//"ERROR getting varid from variable frac in file "//trim(focndomain))
      end if
      rcode = nf90_get_var(ncid_ocn,fracid,ocn_frac,start=(/closelonidx,closelatidx/),count=(/1,1/))
      if (rcode /= nf90_noerr) then
         call shr_sys_abort(subname//"ERROR getting ocean fraction from "//trim(focndomain))
      end if

      !--- Set the appropriate surface flags based on ocean fraction.
   
      if ( present(ocn_present)      )                      ocn_present=(ocn_frac(1,1).gt.0.)
      if ( present(ocn_present).and.present(ice_present))   ice_present=ocn_present
      if ( present(lnd_present))                            lnd_present=(ocn_frac(1,1).lt.1.)
   else if (docn_exists) then
      !--- read in the ocn_in namelist to get name for focndomain file

      unitn = shr_file_getUnit() ! get an unused unit number
      open( unitn, file='docn_in', status='old' )
      rcode = 1
      do while ( rcode /= 0 )
         read (unitn,nml=docn_nml,iostat=rcode)
         if (rcode < 0) then
            call shr_sys_abort( 'shr_scam_checkSurface encountered end-of-file on namelist read' )
         endif
      end do
      close( unitn )
      call shr_file_freeUnit(unitn)
      call shr_strdata_readnml(SDAT,ocn_in)
      call shr_dmodel_readgrid(SDAT%grid,SDAT%gsmap,SDAT%nxg,SDAT%nyg, &
           SDAT%domainfile, ocn_compid, ocn_mpicom, '1d', readfrac=.true., &
           scmmode=.true.,scmlon=scmlon,scmlat=scmlat)
      nfrac = mct_aVect_indexRA(SDAT%grid%data,'frac')
      if ( present(ocn_present)      )                      ocn_present=(SDAT%grid%data%rAttr(nfrac,1).gt.0.)
      if ( present(ocn_present).and.present(ice_present))   ice_present=ocn_present
      if ( present(lnd_present))                            lnd_present=(SDAT%grid%data%rAttr(nfrac,1).lt.1.)
   else
   ! Exit early if no ocn component
      if ( present(ocn_present) ) ocn_present=.false.
      if ( present(ice_present) ) ice_present=.false.
      if ( present(lnd_present) ) lnd_present=.true.
   end if

end subroutine shr_scam_checkSurface

!===============================================================================

end module shr_scam_mod

