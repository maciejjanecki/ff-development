
module cospsimulator_intr
!----------------------------------------------------------------------------------------------------------------------
!Purpose: CAM interface to
!         Name:		CFMIP Observational Simulator Package (COSP)
!         What:		Simulate ISCCP/CloudSat/CALIOP cloud products from GCM inputs
!         Version:	v1.3 released June 2010, updated from v1.1 released May 2009
!         Authors:	Multiple - see http://www.cfmip.net/
!
!Author:  J. Kay (jenkay@ucar.edu) with help from Brian Eaton, John Truesdale, and Y. Zhang/J. Boyle/S. Klein (LLNL)    

! Created: August 2009
! Last modified: August 13, 2010

! NOTES:
! I used ISCCP simulator integration (cloudsimulator_38.F90) as a template.
! This interface is called within radiation. F90.  Andrew says: "clouds and precipitation consistent for the time-step"
! Advice from steve klein on interface design:  "keep the number of changes to COSP core routines to a 
! minimum so that it will be easy to add any updates to the code"
! Keep this interface as self-contained as possible. e.g., get variables from the physics buffer here
! Don't pollute the common space. 

!! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
!! Note: A suggestion was to add all of the CAM variables needed to add to make it possible to run COSP off-line
!! These fields are available and can be called from the namelist though...
!! fincl2 = 'Z3:I','Q:I','T:I','PS:I','CLOUD:I','CONCLD:I','CLDICE:I','CLDLIQ:I','LS_FLXPRC:I','LS_FLXSNW:I', 'ZMFLXPRC:I', 'ZMFLXSNW:I', 'HKFLXPRC:I', 'HKFLXSNW:I',

! FIRST GOAL - GET IT RUNNING WITH CAM4
! I have put "##2check##" next to places in the code where I am unsure what to do or if I have done the right thing.
! random fortran note - spaces don't matter,loops for populating variables from fortran77 not needed here
! Note: These commands start a timer so that you can see how long pieces of the code take to run
!   results in timing file e.g., /ptmp/jenkay/track1_vanilla/ccsm_timing
!   call t_startf ('cospsimulator_intr_run')
!   call t_stopf ('cospsimulator_intr_run')

!----------------------------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use spmd_utils,      only: masterproc
   use ppgrid,          only: pcols, pver, pverp
   use cam_history,     only: addfld, add_default, phys_decomp, outfld, fillvalue
   use perf_mod,        only: t_startf, t_stopf
   use abortutils,      only: endrun
   use cam_pio_utils,   only: max_chars

   implicit none
   private
   save

!Public functions/subroutines

   public :: &
	cospsimulator_intr_readnl,    &
 	cospsimulator_intr_init,    &
	cospsimulator_intr_run

   logical, public :: docosp = .false.  ! whether to do COSP calcs and I/O, default is false
					! if docosp is specified in the atm_in namelist,
					! this value is overwritten and cosp is run

! Private module data


!! Below I set defaults for the CAM and COSP namelists.  Some of the COSP namelist 
!! variables are part of the CAM namelist - they all begin with "cosp_" to keep their 
!! names specific to COSP. I set their CAM namelist defaults here, not in namelist_defaults_cam.xml
!!  Variables identified as namelist variables are defined in
!!  ../models/atm/cam/bld/namelist_files/namelist_definition.xml

!! CAM namelist variable defaults
   logical :: cosp_sample_atrain = .false.    ! CAM namelist variable default, not in COSP namelist
   character(len=256) :: cosp_atrainorbitdata   ! CAM namelist variable, no default, need to specify!
   logical :: cosp_cfmip_3hr = .false.       ! CAM namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_da = .false.        ! CAM namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_off = .false.       ! CAM namelist variable default, not in COSP namelist
   logical :: cosp_cfmip_mon = .false.       ! CAM namelist variable default, not in COSP namelist
   logical :: cosp_lradar_sim = .false.      ! CAM namelist variable default
   logical :: cosp_llidar_sim = .false.      ! CAM namelist variable default
   logical :: cosp_lisccp_sim = .false.      ! CAM namelist variable default
   logical :: cosp_lmisr_sim = .false.       ! CAM namelist variable default
   logical :: cosp_lmodis_sim = .false.      ! CAM namelist variable default (new in v1.3)
   logical :: cosp_lfrac_out = .false.       ! CAM namelist variable default
   integer :: cosp_ncolumns = 50	     ! CAM namelist variable default

!! COSP Namelist variables from cosp_output_nl.txt 
   ! Simulator flags
   logical :: lradar_sim = .false.   		! COSP namelist variable, can be changed from default by CAM namelist
   logical :: llidar_sim = .false.   		! ""
   logical :: lisccp_sim = .false.  		! ""
   logical :: lmisr_sim	 = .false.  		! ""
   logical :: lmodis_sim = .false.  		! ""
   logical :: lrttov_sim = .false.   	     	! not running rttov, always set to .false. (new in v1.3)

   ! Output variables 
   ! All initialized to .false., set to .true. based on CAM namelist in cospsimulator_intr_run
   logical :: lfrac_out = .false.    		! COSP namelist variable, can be changed from default by CAM namelist
   logical :: lalbisccp = .false.
   logical :: latb532 = .false.
   logical :: lboxptopisccp = .false.
   logical :: lboxtauisccp = .false.
   logical :: lcfad_dbze94 = .false.
   logical :: lcfad_lidarsr532 = .false.
   logical :: lclcalipso = .false.
   logical :: lclhcalipso = .false.
   logical :: lclisccp2 = .false.
   logical :: lcllcalipso = .false.
   logical :: lclmcalipso = .false.
   logical :: lcltcalipso = .false.
   logical :: lctpisccp = .false.
   logical :: ldbze94 = .false.
   logical :: ltauisccp = .false.
   logical :: ltclisccp = .false.
   logical :: lparasol_refl = .false.
   logical :: lclmisr = .false.
   logical :: lmeantbisccp = .false.
   logical :: lmeantbclrisccp = .false.
   logical :: lclcalipso2 = .false.
   logical :: lcltlidarradar = .false.
   logical :: lbeta_mol532 = .false.
   logical :: Llongitude = .false.
   logical :: Llatitude =.false.
   logical :: lcltmodis = .false.
   logical :: lclwmodis = .false.
   logical :: lclimodis = .false.
   logical :: lclhmodis = .false.
   logical :: lclmmodis = .false.
   logical :: lcllmodis = .false.
   logical :: ltautmodis = .false.
   logical :: ltauwmodis = .false.
   logical :: ltauimodis = .false.
   logical :: ltautlogmodis = .false.
   logical :: ltauwlogmodis = .false.
   logical :: ltauilogmodis = .false.
   logical :: lreffclwmodis = .false.
   logical :: lreffclimodis = .false.
   logical :: lpctmodis = .false.
   logical :: llwpmodis = .false.
   logical :: liwpmodis = .false.
   logical :: lclmodis = .false.
   logical :: ltbrttov = .false.  !! RTTOV OUTPUT (always set to .false.)

! COSP namelist variables from cosp_input_nl.txt
! Set default values, values discussed with Yuying in ()
! Values from cosp_test.f90 case released with cospv1.1 were used as a template
! Note: Unless otherwise specified, these are parameters that cannot be set by the CAM namelist.

   integer, parameter :: Npoints_it = 10000		! Max # gridpoints to be processed in one iteration (10,000)
   integer :: ncolumns = 50				! Number of subcolumns in SCOPS (50), can be changed from default by CAM namelist
   integer :: nlr = 40          			! Number of levels in statistical outputs 
							! (only used if USE_VGRID=.true.)  (40)
   logical :: use_vgrid	= .true.			! Use fixed vertical grid for outputs? 
							! (if .true. then define # of levels with nlr)	(.true.)
   logical :: csat_vgrid = .true.			! CloudSat vertical grid? 
							! (if .true. then the CloudSat standard grid is used. 
							! If set, overides use_vgrid.) (.true.)
   ! namelist variables for COSP input related to radar simulator
   real(r8) :: radar_freq = 94.0			! CloudSat radar frequency (GHz) (94.0)
   integer :: surface_radar = 0				! surface=1, spaceborne=0 (0)
   integer :: use_mie_tables = 0			! use a precomputed lookup table? yes=1,no=0 (0)
   integer :: use_gas_abs = 1				! include gaseous absorption? yes=1,no=0 (1)
   integer :: do_ray = 0				! calculate/output Rayleigh refl=1, not=0 (0)
   integer :: melt_lay = 0				! melting layer model off=0, on=1 (0)
   real(r8) :: k2 = -1					! |K|^2, -1=use frequency dependent default (-1)
   ! namelist variables for COSP input related to lidar simulator
   integer, parameter :: Nprmts_max_hydro = 12		! Max # params for hydrometeor size distributions (12)
   integer, parameter :: Naero = 1			! Number of aerosol species (Not used) (1)
   integer, parameter :: Nprmts_max_aero = 1		! Max # params for aerosol size distributions (not used) (1)
   integer :: lidar_ice_type = 0			! Ice particle shape in lidar calculations 
							! (0=ice-spheres ; 1=ice-non-spherical) (0)
   integer, parameter :: overlap = 3   		        ! overlap type: 1=max, 2=rand, 3=max/rand (3)

   !! namelist variables for COSP input related to ISCCP simulator
   integer :: isccp_topheight = 1  			! 1 = adjust top height using both a computed infrared 
							! brightness temperature and the visible
                       					! optical depth to adjust cloud top pressure. 
							! Note that this calculation is most appropriate to compare
                       					! to ISCCP data during sunlit hours.
                      					! 2 = do not adjust top height, that is cloud top pressure 
							! is the actual cloud top pressure in the model
                      					! 3 = adjust top height using only the computed infrared 
							! brightness temperature. Note that this calculation is most 
							! appropriate to compare to ISCCP IR only algortihm (i.e. 
							! you can compare to nighttime ISCCP data with this option) (1)
   integer :: isccp_topheight_direction = 2   		! direction for finding atmosphere pressure level with 
							! interpolated temperature equal to the radiance
                                 			! determined cloud-top temperature
                                 			! 1 = find the *lowest* altitude (highest pressure) level 
							! with interpolated temperature 
							! equal to the radiance determined cloud-top temperature
                                 			! 2 = find the *highest* altitude (lowest pressure) level 
							! with interpolated temperature
							! equal to the radiance determined cloud-top temperature
                                 			! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 			! 1 = default setting in COSP v1.1, matches all versions of 
                                 			! ISCCP simulator with versions numbers 3.5.1 and lower
							! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator
   ! RTTOV inputs (not used)
   integer, parameter :: Platform = 1			! satellite platform (1)
   integer, parameter :: Satellite = 15			! satellite (15)
   integer, parameter :: Instrument = 0			! instrument (0)
   integer, parameter :: Nchannels = 8			! Number of channels to be computed (8)
   integer, parameter :: Channels(Nchannels) =  (/1,3,5,6,8,10,11,13/)	
							! Channel numbers (match and supply Nchannels) 
							! (1,3,5,6,8,10,11,13,)
   real(r8), parameter :: Surfem(Nchannels) = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8/)		
							! Surface emissivity (match and supply Nchannels)
							! (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,)
   real(r8), parameter :: ZenAng = 50._r8		! Satellite Zenith Angle (50)
   real(r8), parameter :: co = 	2.098e-07_r8 		! Mixing ratio CO (2.098e-07), not used in cospv1.1
							! set to value from cosp_test.F90's cosp_input_nl.txt
   !per Alejandro: cosp rttov gaseous inputs are also mass mixing ratios
   !values from cosp_test.F90, units (kg/kg)
   !Mixing ratio C02 (5.241e-04), Mixing ratio CH4 (9.139e-07), 
   !Mixing ratio N20 (4.665e-07), Mixing ratio CO (2.098e-07)
   !I get CO2, CH4, N20 from cam radiation interface.

!! Other variables
    integer,parameter :: nhydro = 9                     ! number of COSP hydrometeor classes
    logical :: cosp_runall = .false.     		! flag to run all of the cosp simulator package

!!! Variables read in from atrain orbit file (private module data)
    integer :: norbitdata = 9866324
    real(r8) :: atrainlat(0:9866323) 	        !  A-train orbit latitudes (float in original data file)
    real(r8) :: atrainlon(0:9866323)            !  A-train orbit longitudes
    integer :: atrainday(0:9866323)             !  A-train calendar day (short in data file)
    integer :: atrainhr(0:9866323)              !  A-train hour (byte in data file)
    integer :: atrainmin(0:9866323)             !  A-train minute (byte in data file)
    integer :: atrainsec(0:9866323)             !  A-train minute (byte in data file)
    integer :: idxas = 0			! index to start loop over atrain orbit

CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! subroutine to read namelist variables and run setcospvalues subroutine
! note: cldfrc_readnl is a good template in cloud_fraction.F90
! Make sure that this routine is reading in a namelist.  
! models/atm/cam/bld/build-namelist is the perl script to check

subroutine cospsimulator_intr_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand
   use cosp_share, only: setcospvalues

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input  (nlfile=atm_in)

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cospsimulator_intr_readnl'

    !!! this list should include any variable that you might want to include in the namelist
    !!! philosophy is to not include COSP output flags but just important COSP settings and cfmip controls. 
    namelist /cospsimulator_nl/ docosp, cosp_atrainorbitdata, cosp_cfmip_3hr, cosp_cfmip_da, cosp_cfmip_mon,&
       cosp_cfmip_off, cosp_lfrac_out, cosp_lradar_sim, cosp_llidar_sim, cosp_lisccp_sim, cosp_lmisr_sim,   &
       cosp_lmodis_sim, cosp_ncolumns, cosp_sample_atrain
   !-----------------------------------------------------------------------------

   !! read in the namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )  !! presumably opens the namelist file "nlfile"
      !! position the file to write to the cospsimulator portion of the cam_in namelist
      call find_group_name(unitn, 'cospsimulator_nl', status=ierr)   
      if (ierr == 0) then
	   read(unitn, cospsimulator_nl, iostat=ierr)
	   if (ierr /= 0) then
	      call endrun(subname // ':: ERROR reading namelist')
	   end if
      end if
	   close(unitn)
	   call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(docosp,               1,  mpilog, 0, mpicom)
   call mpibcast(cosp_atrainorbitdata, len(cosp_atrainorbitdata), mpichar, 0, mpicom)
   call mpibcast(cosp_cfmip_3hr,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_da,        1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_mon,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_cfmip_off,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lfrac_out,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lradar_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_llidar_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lisccp_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lmisr_sim,       1,  mpilog, 0, mpicom)
   call mpibcast(cosp_lmodis_sim,      1,  mpilog, 0, mpicom)
   call mpibcast(cosp_ncolumns,        1,  mpiint, 0, mpicom)
   call mpibcast(cosp_sample_atrain,   1,  mpilog, 0, mpicom)
#endif

   !! reset COSP namelist variables based on input from cam namelist variables
   if (cosp_cfmip_3hr) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_da) then
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_off) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
   end if
   if (cosp_cfmip_mon) then
      llidar_sim = .true.
      lisccp_sim = .true.
   end if

   if (cosp_lfrac_out) then
      lfrac_out = .true.
   end if
   if (cosp_lradar_sim) then
      lradar_sim = .true.
   end if
   if (cosp_llidar_sim) then
      llidar_sim = .true.
   end if
   if (cosp_lisccp_sim) then
      lisccp_sim = .true.
   end if
   if (cosp_lmisr_sim) then
      lmisr_sim = .true.
   end if
   if (cosp_lmodis_sim) then
      lmodis_sim = .true.
   end if

   !! reset COSP namelist variables based on input from cam namelist variables
   if (cosp_ncolumns .ne. ncolumns) then
      ncolumns = cosp_ncolumns
   end if

   !! if no simulators are turned on at all and docosp is, set cosp_runall = .true.
   if((docosp) .and. (.not.lradar_sim) .and. (.not.llidar_sim) .and. (.not.lisccp_sim) .and. (.not.lmisr_sim) .and. (.not.lmodis_sim)) then
      cosp_runall = .true.
   end if
   if (cosp_runall) then
      lradar_sim = .true.
      llidar_sim = .true.
      lisccp_sim = .true.
      lmisr_sim = .true.
      lmodis_sim = .true.
      lfrac_out = .true.
   end if

   !! use the namelist variables to overwrite default COSP output variables above (set to .false.)
   if (lradar_sim) then
      !! turn on the outputs for the radar simulator
      lcfad_dbze94 = .true.
      ldbze94 = .true.
      lclcalipso2 = .true.
      lcltlidarradar = .true.
   end if
   if (llidar_sim) then
      !! turn on the outputs for the lidar simulator
      lcllcalipso = .true.
      lclmcalipso = .true.
      lcltcalipso = .true.
      lclcalipso = .true.
      lclhcalipso = .true.
      lcfad_lidarsr532 = .true.
      latb532 = .true.
      lparasol_refl = .true.
      lbeta_mol532 = .true.
   end if
   if (lisccp_sim) then
      !! turn on the outputs for the isccp simulator
      lalbisccp = .true.
      lboxptopisccp = .true.
      lboxtauisccp = .true.
      lclisccp2 = .true.
      lctpisccp = .true.
      ltauisccp = .true.
      ltclisccp = .true.
      lmeantbisccp = .true.
      lmeantbclrisccp = .true.
   end if
   if (lmisr_sim) then
      !! turn on the outputs for the misr simulator
      lclmisr = .true.
   end if
   if (lmodis_sim) then
      !! turn on the outputs for the modis simulator
      lcltmodis = .true.
      lclwmodis = .true.
      lclimodis = .true.
      lclhmodis = .true.
      lclmmodis = .true.
      lcllmodis = .true.
      ltautmodis = .true.
      ltauwmodis = .true.
      ltauimodis = .true.
      ltautlogmodis = .true.
      ltauwlogmodis = .true.
      ltauilogmodis = .true.
      lreffclwmodis = .true.
      lreffclimodis = .true.
      lpctmodis = .true.
      llwpmodis = .true.
      liwpmodis = .true.
      lclmodis = .true.
   end if

   ! Set vertical coordinate and subcolumn cosp options based on namelist inputs
   call setcospvalues(nlr,use_vgrid,csat_vgrid,ncolumns,docosp)

end subroutine cospsimulator_intr_readnl

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine cospsimulator_intr_init

   use cosp_share,          only: nprs_cosp,ntau_cosp,ntau_cosp_modis,&
					nht_cosp,ndbze_cosp,&
					nsr_cosp,nscol_cosp,&
					nhtmisr_cosp,nsza_cosp,&
					nhtml_cosp
   use cam_history,         only: addfld, add_default, phys_decomp

   use mpishorthand
   use netcdf
   use error_messages,      only : handle_ncerr

   integer ncid,latid,lonid,did,hrid,minid,secid
   !------------------------------------------------------------------------------


if (cosp_sample_atrain) then

!!!! READ IN ATRAIN ORBIT DATA FROM INPUT FILE FOR SUB-SAMPLING
!!! What are the line numbers for? ##2check##

   if (masterproc) then
      call handle_ncerr( nf90_open (cosp_atrainorbitdata,NF90_NOWRITE,ncid),&
        'cospsimulator_intr.F90:304')
      call handle_ncerr( nf90_inq_varid (ncid,'lat',latid),&
        'cospsimulator_intr.F90:306')
      call handle_ncerr( nf90_inq_varid (ncid,'lon',lonid),&
        'cospsimulator_intr.F90:308')
      call handle_ncerr( nf90_inq_varid (ncid,'doy365',did),&
        'cospsimulator_intr.F90:310')
      call handle_ncerr( nf90_inq_varid (ncid,'hour',hrid),&
        'cospsimulator_intr.F90:312')
      call handle_ncerr( nf90_inq_varid (ncid,'minute',minid),&
        'cospsimulator_intr.F90:314')
      call handle_ncerr( nf90_inq_varid (ncid,'second',secid),&
        'cospsimulator_intr.F90:316')
      call handle_ncerr( nf90_get_var (ncid,latid,atrainlat),&
        'cospsimulator_intr.F90:318')
      call handle_ncerr( nf90_get_var (ncid,lonid,atrainlon),&
        'cospsimulator_intr.F90:320')
      call handle_ncerr( nf90_get_var (ncid,did,atrainday),&
        'cospsimulator_intr.F90:322')
      call handle_ncerr( nf90_get_var (ncid,hrid,atrainhr),&
        'cospsimulator_intr.F90:324')
      call handle_ncerr( nf90_get_var (ncid,minid,atrainmin),&
        'cospsimulator_intr.F90:326')
      call handle_ncerr( nf90_get_var (ncid,secid,atrainsec),&
        'cospsimulator_intr.F90:328')
      call handle_ncerr( nf90_close (ncid),&
        'cospsimulator_intr.F90:330')
   end if
#if ( defined SPMD )
   call mpibcast (atrainlat,norbitdata , mpir8, 0, mpicom)  !! 2 arg = size, see runtime_opts.F90
   call mpibcast (atrainlon,norbitdata , mpir8, 0, mpicom)
   call mpibcast (atrainday,norbitdata , mpir8, 0, mpicom)
   call mpibcast (atrainhr,norbitdata , mpir8, 0, mpicom)
   call mpibcast (atrainmin,norbitdata , mpir8, 0, mpicom)
   call mpibcast (atrainsec,norbitdata , mpir8, 0, mpicom)
#endif

endif

! ADDFLD ADD_DEFAULT CALLS FOR COSP OUTPUTS
! Notes on addfld/add_default/outfld calls:  
! 1) Dimensions of cosp output should be:  
!	a) ntime=1 (cosp run at every time step)
! 	b) nprofile=ncol
! 2) See cam_history.F90, modified for new cosp output dimensions
! 3) Need conditional logic so that addfld,add_default,and outfld calls are only executed when fields are available.
! 4) nhtml_cosp=height_mlev=pver (height of model levels)
! 5) flag_xyfill=.true. is "non-applicable xy points flagged with fillvalue".  
!  per brian, flag_xyfill should be set to true whenever the fields might contain fillvalues.
!  I think this should be .true. for all COSP outputs.
!  Especially because there are cases where there will definitely be fill values (e.g., tau not calculated when no cloud.)
!  For cosp variables with subcolumns (dimension includes nscol_cosp) I have made the outputs instantaneous by default
!  to get around check_accum failing and aborting run in cam_history.90.  Problem is that the vertical dimension
!  can contain a mix of fillvalue and realvalue (e.g., when there is a cloud in one sub-column but not in another).
!  none of these variables are a part of CFMIP.  Also needed to modify cam_history so that the check_accum is
!  not done when output is instantaneous.
! 6) sampling_seq = radiation timestep.  note: cloudsimulator.F90 does not specify anything special.
! 7) Markings for CFMIP output requirements:
!*cfMon* = CFMIP variable for cfMon - monthly-mean cloud diagnostic fields
!*cfOff* = CFMIP variable for cfOff - monthly-mean offline cloud diagnostic fields
!*cfDa* = CFMIP variable for cfDa - daily-mean cloud diagnostic fields
!*cf3hr* = CFMIP variable for cf3hr - 3-hourly cloud diagnostic fields
! 8) Making it the default to produce a separate CAM history file with COSP outputs.  
! "2" is the history tape index.  I specify 2 here to create a separate history file for cosp output. 
! 9) crash fix: add_default was looking through the master list for input field name and not finding it.
! Solution? I removed all of the spaces in the addfld calls, add_default, and outfld calls.

!!! ISCCP OUTPUTS
   if (lisccp_sim) then
      !! addfld calls for all
      !*cfMon,cfDa* clisccp2 (time,tau,plev,profile), CFMIP wants 7 p bins, 7 tau bins
      call addfld('FISCCP1_COSP','mixed   ',nprs_cosp*ntau_cosp,'A', &
		   'Grid-box fraction covered by each ISCCP D level cloud type',phys_decomp,&
		   flag_xyfill=.true., flag_cospprstaulev=.true.)
      !*cfMon,cfDa* tclisccp (time,profile), CFMIP wants "gridbox mean cloud cover from ISCCP"
      call addfld('CLDTOT_ISCCP','fraction', 1,'A', &
		   'Total Cloud Fraction Calculated by the ISCCP Simulator ',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfDa* albisccp (time,profile)
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)
      call addfld('MEANCLDALB_ISCCP','unitless',1,'A','Mean cloud albedo*CLDTOT_ISCCP',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfDa* ctpisccp (time,profile)
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction (divide by CLDTOT_ISSCP in history file to get weighted average)	
      call addfld('MEANPTOP_ISCCP','mb',1,'A','Mean cloud top pressure*CLDTOT_ISCCP',phys_decomp,flag_xyfill=.true.)
      ! boxtauisccp (time,column,profile)
      call addfld ('TAU_ISCCP','fraction',nscol_cosp,'I','Optical Depth in each Subcolumn',&
		   phys_decomp,flag_xyfill=.true.,flag_cospscol=.true.)
      ! boxptopisccp (time,column,profile)
      call addfld ('CLDPTOP_ISCCP','Pa',nscol_cosp,'I','Cloud Top Pressure in each Subcolumn',&
		   phys_decomp,flag_xyfill=.true.,flag_cospscol=.true.)
      ! tauisccp (time,profile)
      call addfld ('MEANTAU_ISCCP','unitless',1,'A','Mean optical thickness from ISCCP simulator',phys_decomp,flag_xyfill=.true.)
      ! meantbisccp (time,profile), at 10.5 um
      call addfld ('MEANTB_ISCCP','K       ',1,'A','Mean Infrared Tb from ISCCP simulator',phys_decomp,flag_xyfill=.true.)
      ! meantbclrisccp (time,profile)
      call addfld ('MEANTBCLR_ISCCP','K       ',1,'A','Mean Clear-sky Infrared Tb from ISCCP simulator',&
		   phys_decomp,flag_xyfill=.true.)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file
      if (cosp_cfmip_mon.or.cosp_cfmip_da) then
 	   !! add cfmip-requested variables to two separate cam history files
	 if (cosp_cfmip_da) then
	    call add_default ('FISCCP1_COSP',2,' ')
	    call add_default ('CLDTOT_ISCCP',2,' ')
	    call add_default ('MEANCLDALB_ISCCP',2,' ')
	    call add_default ('MEANPTOP_ISCCP',2,' ')
	 end if
	 if (cosp_cfmip_mon) then
	    call add_default ('FISCCP1_COSP',1,' ')
	    call add_default ('CLDTOT_ISCCP',1,' ')
	    call add_default ('MEANCLDALB_ISCCP',1,' ')
	    call add_default ('MEANPTOP_ISCCP',1,' ')
	 end if
      else
 	   !! add all isccp outputs to the history file
	    call add_default ('FISCCP1_COSP',2,' ')
	    call add_default ('CLDTOT_ISCCP',2,' ')
	    call add_default ('MEANCLDALB_ISCCP',2,' ')
	    call add_default ('MEANPTOP_ISCCP',2,' ')
	    call add_default ('TAU_ISCCP',2,' ')
	    call add_default ('CLDPTOP_ISCCP',2,' ')
	    call add_default ('MEANTAU_ISCCP',2,' ')
	    call add_default ('MEANTB_ISCCP',2,' ')
	    call add_default ('MEANTBCLR_ISCCP',2,' ')
      end if
   end if

!!! LIDAR SIMULATOR OUTPUTS
   if (llidar_sim) then
      !! addfld calls for all
      !*cfMon,cfOff,cfDa,cf3hr* cllcalipso (time,profile)
      call addfld('CLDLOW_CAL','fraction',1,'A','Lidar Low-level Cloud Fraction',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfOff,cfDa,cf3hr* clmcalipso (time,profile)
      call addfld('CLDMED_CAL','fraction',1,'A','Lidar Mid-level Cloud Fraction',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfOff,cfDa,cf3hr* clhcalipso (time,profile)
      call addfld('CLDHGH_CAL','fraction',1,'A','Lidar High-level Cloud Fraction',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfOff,cfDa,cf3hr* cltcalipso (time,profile)
      call addfld('CLDTOT_CAL','fraction',1,'A','Lidar Total Cloud Fraction',phys_decomp,flag_xyfill=.true.)
      !*cfMon,cfOff,cfDa,cf3hr* clcalipso (time,height,profile)
      call addfld('CLD_CAL','fraction',nht_cosp,'A','Lidar Cloud Fraction (532 nm)',&
		   phys_decomp, flag_xyfill=.true., flag_cospht=.true.)
      !*cfMon,cfOff,cfDa,cf3hr* parasol_refl (time,sza,profile)
      call addfld ('RFL_PARASOL','fraction',nsza_cosp,'A','PARASOL-like mono-directional reflectance ',&
		   phys_decomp,flag_xyfill=.true.,flag_cospsza=.true.)
      !*cfOff,cf3hr* cfad_lidarsr532 (time,height,scat_ratio,profile), %11%, default is 40 vert levs, 15 SR  bins
      call addfld('CFAD_SR532_CAL','fraction',nht_cosp*nsr_cosp,'A',&
		   'Lidar Scattering Ratio CFAD (532 nm)',phys_decomp,&
		   flag_xyfill=.true., flag_cosphtsrlev=.true.) 
      ! atb532 (time,height_mlev,column,profile)
      call addfld ('ATB532_CAL','no_unit_log10(x)',nhtml_cosp*nscol_cosp,'I', &
		   'Lidar Attenuated Total Backscatter (532 nm) in each Subcolumn',phys_decomp, &
		   flag_xyfill=.true.,flag_cosphtmlscollev=.true.)
      ! beta_mol532 (time,height_mlev,profile)
      call addfld ('MOL532_CAL','m-1sr-1',nhtml_cosp,'A','Lidar Molecular Backscatter (532 nm) ',&
		   phys_decomp,flag_xyfill=.true.,flag_cospht=.true.)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file
      if (cosp_cfmip_mon .or. cosp_cfmip_off .or. cosp_cfmip_da .or. cosp_cfmip_3hr) then
	 if (cosp_cfmip_da) then
	    call add_default ('CLDLOW_CAL',2,' ')
	    call add_default ('CLDMED_CAL',2,' ')
	    call add_default ('CLDHGH_CAL',2,' ')
	    call add_default ('CLDTOT_CAL',2,' ')
	    call add_default ('CLD_CAL',2,' ')
	    call add_default ('RFL_PARASOL',2,' ')
	 end if
	 if (cosp_cfmip_mon.or.cosp_cfmip_off) then
	    call add_default ('CLDLOW_CAL',1,' ')
	    call add_default ('CLDMED_CAL',1,' ')
	    call add_default ('CLDHGH_CAL',1,' ')
	    call add_default ('CLDTOT_CAL',1,' ')
	    call add_default ('CLD_CAL',1,' ')
	    call add_default ('RFL_PARASOL',1,' ')
	 end if
	 if (cosp_cfmip_3hr) then
	    call add_default ('CFAD_SR532_CAL',3,' ')
	    call add_default ('CLDLOW_CAL',3,' ')
	    call add_default ('CLDMED_CAL',3,' ')
	    call add_default ('CLDHGH_CAL',3,' ')
	    call add_default ('CLDTOT_CAL',3,' ')
	    call add_default ('CLD_CAL',3,' ')
	    call add_default ('RFL_PARASOL',3,' ')
	 end if
	 if (cosp_cfmip_off) then
	    call add_default ('CFAD_SR532_CAL',1,' ')
	 end if
      else
 	 !! add all lidar outputs to the history file
	 call add_default ('CLDLOW_CAL',2,' ')
	 call add_default ('CLDMED_CAL',2,' ')
	 call add_default ('CLDHGH_CAL',2,' ')
	 call add_default ('CLDTOT_CAL',2,' ')
	 call add_default ('CLD_CAL',2,' ')
	 call add_default ('RFL_PARASOL',2,' ')
	 call add_default ('CFAD_SR532_CAL',2,' ')
	 call add_default ('ATB532_CAL',2,' ')
	 call add_default ('MOL532_CAL',2,' ')
      end if

   end if

!!! RADAR SIMULATOR OUTPUTS
   if (lradar_sim) then
      !!! addfld calls
      !*cfOff,cf3hr* cfad_dbze94 (time,height,dbze,profile), default is 40 vert levs, 15 dBZ bins 
      call addfld('CFAD_DBZE94_CS','fraction',nht_cosp*ndbze_cosp,'A',&
		   'Radar Reflectivity Factor CFAD (94 GHz)',phys_decomp,&
		   flag_xyfill=.true., flag_cosphtdbzelev=.true.)
      !*cfOff,cf3hr* clcalipso2 (time,height,profile)
      call addfld ('CLD_CAL_NOTCS','fraction',nht_cosp,'A','Cloud occurrence seen by CALIPSO but not CloudSat ',&
		   phys_decomp,flag_xyfill=.true.,flag_cospht=.true.)
      ! dbze94 (time,height_mlev,column,profile),! height_mlevel = height when vgrid_in = .true. (default)
      call addfld ('DBZE_CS','dBZe    ',nhtml_cosp*nscol_cosp,'I',' Radar dBZe (94 GHz) in each Subcolumn',phys_decomp,&
		   flag_xyfill=.true.,flag_cosphtmlscollev=.true.)
      ! cltlidarradar (time,profile)
      call addfld ('CLDTOT_CALCS','fraction',1,'A',' Lidar and Radar Total Cloud Fraction ',phys_decomp,flag_xyfill=.true.)

      !!! add_default calls for CFMIP experiments or else all fields are added to history file
       if (cosp_cfmip_off.or.cosp_cfmip_3hr) then
	  if (cosp_cfmip_3hr) then
	      call add_default ('CFAD_DBZE94_CS',3,' ')
	      call add_default ('CLD_CAL_NOTCS',3,' ')
	  end if
	  if (cosp_cfmip_off) then
	      call add_default ('CFAD_DBZE94_CS',1,' ')
	      call add_default ('CLD_CAL_NOTCS',1,' ')
	  end if
      else
 	 !! add all radar outputs to the history file
	  call add_default ('CFAD_DBZE94_CS',2,' ')
	  call add_default ('CLD_CAL_NOTCS',2,' ')
	  call add_default ('DBZE_CS',2,' ')
	  call add_default ('CLDTOT_CALCS',2,' ')
      end if
   end if

!!! MISR SIMULATOR OUTPUTS
   if (lmisr_sim) then
      ! clMISR (time,tau,CTH_height_bin,profile)
      call addfld ('CLD_MISR','fraction',nhtmisr_cosp*ntau_cosp,'A','Cloud Fraction from MISR Simulator',&
		   phys_decomp,flag_xyfill=.true.,flag_cosphtmisrtaulev=.true.)
      call add_default ('CLD_MISR',2,' ')
   end if

!!! SUB-COLUMN OUTPUT
   if (lfrac_out) then
      ! frac_out (time,height_mlev,column,profile)
      call addfld ('SCOPS_OUT','0=nocld,1=strcld,2=cnvcld',nhtml_cosp*nscol_cosp,'I','SCOPS Subcolumn output',&
		   phys_decomp,flag_xyfill=.true.,flag_cosphtmlscollev=.true.)
      call add_default ('SCOPS_OUT',2,' ')
   end if 

!!! MODIS OUTPUT
   if (lmodis_sim) then
      ! float cltmodis ( time, loc )
      call addfld ('CLTMODIS','%',1,'A','MODIS Total Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLTMODIS',2,' ')
      ! float clwmodis ( time, loc )
      call addfld ('CLWMODIS','%',1,'A','MODIS Liquid Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLWMODIS',2,' ')
      ! float climodis ( time, loc )
      call addfld ('CLIMODIS','%',1,'A','MODIS Ice Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLIMODIS',2,' ')
      ! float clhmodis ( time, loc )
      call addfld ('CLHMODIS','%',1,'A','MODIS High Level Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLHMODIS',2,' ')
      ! float clmmodis ( time, loc )
      call addfld ('CLMMODIS','%',1,'A','MODIS Mid Level Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLMMODIS',2,' ')
      ! float cllmodis ( time, loc )
      call addfld ('CLLMODIS','%',1,'A','MODIS Low Level Cloud Fraction',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('CLLMODIS',2,' ')
      ! float tautmodis ( time, loc )
      call addfld ('TAUTMODIS','1',1,'A','MODIS Total Cloud Optical Thickness',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUTMODIS',2,' ')
      ! float tauwmodis ( time, loc )
      call addfld ('TAUWMODIS','1',1,'A','MODIS Liquid Cloud Optical Thickness',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUWMODIS',2,' ')
      ! float tauimodis ( time, loc )
      call addfld ('TAUIMODIS','1',1,'A','MODIS Ice Cloud Optical Thickness',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUIMODIS',2,' ')
      ! float tautlogmodis ( time, loc )
      call addfld ('TAUTLOGMODIS','1',1,'A','MODIS Total Cloud Optical Thickness (Log10 Mean)',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUTLOGMODIS',2,' ')
      ! float tauwlogmodis ( time, loc )
      call addfld ('TAUWLOGMODIS','1',1,'A','MODIS Liquid Cloud Optical Thickness (Log10 Mean)',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUWLOGMODIS',2,' ')
      ! float tauilogmodis ( time, loc )
      call addfld ('TAUILOGMODIS','1',1,'A','MODIS Ice Cloud Optical Thickness (Log10 Mean)',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('TAUILOGMODIS',2,' ')
      ! float reffclwmodis ( time, loc )
      call addfld ('REFFCLWMODIS','m',1,'A','MODIS Liquid Cloud Particle Size',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('REFFCLWMODIS',2,' ')
      ! float reffclimodis ( time, loc )
      call addfld ('REFFCLIMODIS','m',1,'A','MODIS Ice Cloud Particle Size',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('REFFCLIMODIS',2,' ')
      ! float pctmodis ( time, loc )
      call addfld ('PCTMODIS','Pa',1,'A','MODIS Cloud Top Pressure',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('PCTMODIS',2,' ')
      ! float lwpmodis ( time, loc )
      call addfld ('LWPMODIS','kg m-2',1,'A','MODIS Cloud Liquid Water Path',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('LWPMODIS',2,' ')
      ! float iwpmodis ( time, loc )
      call addfld ('IWPMODIS','kg m-2',1,'A','MODIS Cloud Ice Water Path',&
		   phys_decomp,flag_xyfill=.true.)
      call add_default ('IWPMODIS',2,' ')
      ! float clmodis ( time, plev, tau, loc )
      call addfld ('CLMODIS','%',nprs_cosp*ntau_cosp_modis,'A','MODIS Cloud Area Fraction',&
		   phys_decomp,flag_xyfill=.true., flag_cospprstaumodislev=.true.)
      call add_default ('CLMODIS',2,' ')
   end if

!! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
!! Note: A suggestion was to add all of the CAM variables needed to add to make it possible to run COSP off-line
!! These fields are available and can be called from the namelist thought.

end subroutine cospsimulator_intr_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine cospsimulator_intr_run(state,pbuf,cam_in,emis,cliqwp,cicewp,coszrs)    

   use physics_types,    only: physics_state
   use phys_buffer,	 only: pbuf_fld,pbuf_size_max,pbuf_old_tim_idx,pbuf_get_fld_idx
   use camsrfexch_types, only: cam_in_t
   use constituents, 	 only: cnst_get_ind
   use rad_constituents, only: rad_cnst_get_clim_gas
   use wv_saturation, 	 only: aqsat_water
   use phys_control,     only: phys_getopts
   use interpolate_data, only: lininterp_init,lininterp,lininterp_finish,interp_type    
   use physconst,	 only: pi
   use cam_history,      only: outfld, fillvalue
   use cmparray_mod,     only: CmpDayNite, ExpDayNite
   use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
   use time_manager,     only: get_curr_calday,get_curr_time,get_calendar,get_ref_date

   use cosp_share,       only: nprs_cosp,ntau_cosp,ntau_cosp_modis,&
					nht_cosp,ndbze_cosp,&
					nsr_cosp,nscol_cosp,&
					nhtmisr_cosp,nsza_cosp,&
					nhtml_cosp

#ifdef USE_COSP
   ! cosp simulator package, COSP code was copied into CAM source tree and is compiled as a separate library
   use mod_cosp_constants    !! can i use this to access I_LSCLIQ = 1 etc.
   use mod_cosp_types,	 only: construct_cosp_gridbox,construct_cosp_vgrid,construct_cosp_subgrid, &
				construct_cosp_sgradar,construct_cosp_radarstats,construct_cosp_sglidar, &
				construct_cosp_lidarstats,construct_cosp_isccp,construct_cosp_misr, &
				free_cosp_gridbox,free_cosp_vgrid,free_cosp_subgrid, &
				free_cosp_sgradar,free_cosp_radarstats,free_cosp_sglidar, &
				free_cosp_lidarstats,free_cosp_isccp,free_cosp_misr,&
				cosp_config,cosp_gridbox,cosp_subgrid,cosp_sgradar, &
				cosp_sglidar,cosp_isccp,cosp_misr,cosp_vgrid,cosp_radarstats,&
				cosp_radarstats,cosp_lidarstats
   use mod_cosp_modis_simulator,	 only: construct_cosp_modis,free_cosp_modis,cosp_modis
   use mod_cosp, 	 only: cosp
#endif

! Arguments
   type(physics_state), intent(in) :: state
   type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf
   type(cam_in_t),  intent(in) :: cam_in

   !! vars calculated in subroutine param_cldoptics_calc within param_cldoptics.F90
   real(r8), intent(in) :: cliqwp(pcols,pver)  	   	! in-cloud liquid water path, 
   real(r8), intent(in) :: cicewp(pcols,pver)      	! in-cloud ice water path
   real(r8), intent(in) :: emis(pcols,pver)        	! cloud longwave emissivity
   real(r8), intent(in) :: coszrs(pcols)       		! cosine solar zenith angle (to tell if day or night)

#ifdef USE_COSP
! Local variables
   ! generic
   integer :: lchnk                         		! chunk identifier
   integer :: ncol                          		! number of active atmospheric columns
   integer :: i, k, ip, it, ipt, ih, id, ihd, &
	is, ihs, isc, ihsc, ihm, ihmt, ihml, &
	itim, ifld       				! indices

   ! variables for day/nite and orbital subsetting
    ! Gathered indicies of day and night columns 
    !  chunk_column_index = IdxDay(daylight_column_index)
    integer :: Nday                      ! Number of daylight columns
    integer :: Nno                       ! Number of columns not using for simulator
    integer, dimension(pcols) :: IdxDay  ! Indices of daylight columns
    integer, dimension(pcols) :: IdxNo   ! Indices of columns not using for simulator
    real(r8) :: tmp(pcols)		 ! tempororary variable for array expansion
    real(r8) :: tmp1(pcols,pver)	 ! tempororary variable for array expansion
    real(r8) :: tmp2(pcols,pver)	 ! tempororary variable for array expansion
    real(r8) :: lon_cosp_day(pcols)      ! tempororary variable for sunlit lons
    real(r8) :: lat_cosp_day(pcols)      ! tempororary variable for sunlit lats
    real(r8) :: ptop_day(pcols,pver) 	 ! tempororary variable for sunlit ptop
    real(r8) :: pmid_day(pcols,pver) 	 ! tempororary variable for sunlit pmid
    real(r8) :: ztop_day(pcols,pver) 	 ! tempororary variable for sunlit ztop
    real(r8) :: zmid_day(pcols,pver) 	 ! tempororary variable for sunlit zmid
    real(r8) :: t_day(pcols,pver) 	 ! tempororary variable for sunlit t
    real(r8) :: rh_day(pcols,pver) 	 ! tempororary variable for sunlit rh
    real(r8) :: q_day(pcols,pver) 	 ! tempororary variable for sunlit q
    real(r8) :: concld_day(pcols,pver) 	 ! tempororary variable for sunlit concld
    real(r8) :: cld_day(pcols,pver) 	 ! tempororary variable for sunlit cld
    real(r8) :: ps_day(pcols) 		 ! tempororary variable for sunlit ps
    real(r8) :: ts_day(pcols) 		 ! tempororary variable for sunlit ts
    real(r8) :: landmask_day(pcols) 	 ! tempororary variable for sunlit landmask
    real(r8) :: o3_day(pcols,pver) 	 ! tempororary variable for sunlit o3
    real(r8) :: us_day(pcols) 		 ! tempororary variable for sunlit us
    real(r8) :: vs_day(pcols) 		 ! tempororary variable for sunlit vs
    real(r8) :: cldliq_day(pcols,pver) 	 ! tempororary variable for sunlit cldliq
    real(r8) :: cldice_day(pcols,pver) 	 ! tempororary variable for sunlit cldice
    real(r8) :: mr_ccliq_day(pcols,pver) 	! tempororary variable for sunlit mr_ccliq
    real(r8) :: mr_ccice_day(pcols,pver) 	 ! tempororary variable for sunlit mr_ccice
    real(r8) :: rain_ls_interp_day(pcols,pver) 	 ! tempororary variable for sunlit rain_ls_interp
    real(r8) :: snow_ls_interp_day(pcols,pver) 	 ! tempororary variable for sunlit snow_ls_interp
    real(r8) :: grpl_ls_interp_day(pcols,pver) 	 ! tempororary variable for sunlit grpl_ls_interp
    real(r8) :: rain_cv_interp_day(pcols,pver) 	 ! tempororary variable for sunlit rain_cv_interp
    real(r8) :: snow_cv_interp_day(pcols,pver) 	 ! tempororary variable for sunlit snow_cv_interp
    real(r8) :: reff_cosp_day(pcols,pver,nhydro) ! tempororary variable for sunlit reff_cosp(:,:,:)
    real(r8) :: dtau_s_day(pcols,pver) 	 ! tempororary variable for sunlit dtau_s
    real(r8) :: dtau_c_day(pcols,pver) 	 ! tempororary variable for sunlit dtau_c
    real(r8) :: dem_s_day(pcols,pver) 	 ! tempororary variable for sunlit dem_s
    real(r8) :: dem_c_day(pcols,pver) 	 ! tempororary variable for sunlit dem_c

    !! vars for atrain orbital sub-sampling
    integer :: Natrain                        ! # of atrain columns
    integer, dimension(pcols) :: IdxAtrain                 
    real(r8) :: cam_calday                    ! current CAM calendar day
    real(r8),dimension(norbitdata) :: atrain_calday ! atrain calendar day (real, includes time of day!)
    integer, dimension(pcols) :: atrain_use
    real(r8) :: atrain_latthresh = 1.0    ! absolute lat threshold in degrees
    real(r8) :: atrain_lonthresh = 1.0    ! absolute lon threshold in degrees
    real(r8) :: atrain_daythresh = 0.0208 ! absolute time threshold in fractional day (30 minutes = cam rad timestep)
    integer,parameter :: ntmp = 30000     ! # points to loop over, now ~#points/day in orbit file, could also use atrain_daythresh to set
    real(r8) :: caldaydiff		  ! temp var for comparing calday differences
    real(r8) :: caldaymatch		  ! temp var for comparing caldays 1 = match, 0 = no match
    real(r8) :: caldaymindiff		  ! temp var, min calday difference btwn atrain and cam
    real(r8), dimension(ntmp) :: tmpdiff  ! temp var for comparing lat/lon
    real(r8), dimension(ntmp) :: tmpvar   ! temp var for lat/lon
    integer, dimension(ntmp) :: idxlat
    integer, dimension(ntmp) :: idxlon
    integer, dimension(ntmp) :: idxsum

    real(r8) :: lon_cosp_atrain(pcols)      	! tempororary variable for atrain lons
    real(r8) :: lat_cosp_atrain(pcols)      	! tempororary variable for atrain lats
    real(r8) :: ptop_atrain(pcols,pver) 	! tempororary variable for atrain ptop
    real(r8) :: pmid_atrain(pcols,pver) 	! tempororary variable for atrain pmid
    real(r8) :: ztop_atrain(pcols,pver) 	! tempororary variable for atrain ztop
    real(r8) :: zmid_atrain(pcols,pver) 	! tempororary variable for atrain zmid
    real(r8) :: t_atrain(pcols,pver) 	 	! tempororary variable for atrain t
    real(r8) :: rh_atrain(pcols,pver) 	 	! tempororary variable for atrain rh
    real(r8) :: q_atrain(pcols,pver) 	 	! tempororary variable for atrain q
    real(r8) :: concld_atrain(pcols,pver) 	! tempororary variable for atrain concld
    real(r8) :: cld_atrain(pcols,pver) 	 	! tempororary variable for atrain cld
    real(r8) :: ps_atrain(pcols) 		! tempororary variable for atrain ps
    real(r8) :: ts_atrain(pcols) 		! tempororary variable for atrain ts
    real(r8) :: landmask_atrain(pcols) 	 	! tempororary variable for atrain landmask
    real(r8) :: o3_atrain(pcols,pver) 	 	! tempororary variable for atrain o3
    real(r8) :: us_atrain(pcols) 		! tempororary variable for atrain us
    real(r8) :: vs_atrain(pcols) 		! tempororary variable for atrain vs
    real(r8) :: cldliq_atrain(pcols,pver) 	! tempororary variable for atrain cldliq
    real(r8) :: cldice_atrain(pcols,pver) 	! tempororary variable for atrain cldice
    real(r8) :: mr_ccliq_atrain(pcols,pver) 	! tempororary variable for atrain mr_ccliq
    real(r8) :: mr_ccice_atrain(pcols,pver) 	! tempororary variable for atrain mr_ccice
    real(r8) :: rain_ls_interp_atrain(pcols,pver) ! tempororary variable for atrain rain_ls_interp
    real(r8) :: snow_ls_interp_atrain(pcols,pver) ! tempororary variable for atrain snow_ls_interp
    real(r8) :: grpl_ls_interp_atrain(pcols,pver) ! tempororary variable for atrain grpl_ls_interp
    real(r8) :: rain_cv_interp_atrain(pcols,pver) ! tempororary variable for atrain rain_cv_interp
    real(r8) :: snow_cv_interp_atrain(pcols,pver) ! tempororary variable for atrain snow_cv_interp
    real(r8) :: reff_cosp_atrain(pcols,pver,nhydro) ! tempororary variable for atrain reff_cosp(:,:,:)
    real(r8) :: dtau_s_atrain(pcols,pver) 	! tempororary variable for atrain dtau_s
    real(r8) :: dtau_c_atrain(pcols,pver) 	! tempororary variable for atrain dtau_c
    real(r8) :: dem_s_atrain(pcols,pver) 	! tempororary variable for atrain dem_s
    real(r8) :: dem_c_atrain(pcols,pver) 	! tempororary variable for atrain dem_c

   ! constants for optical depth calculation (from radcswmx.F90)
   real(r8), parameter :: abarl = 2.817e-02_r8    	! A coefficient for extinction optical depth
   real(r8), parameter :: bbarl = 1.305_r8        	! b coefficient for extinction optical depth
   real(r8), parameter :: abari = 3.448e-03_r8    	! A coefficient for extinction optical depth
   real(r8), parameter :: bbari = 2.431_r8        	! b coefficient for extinction optical depth
   real(r8), parameter :: cldmin = 1.0e-80_r8  		! note: cldmin much less than cldmin from cldnrh
   real(r8), parameter :: cldeps = 0.0_r8 

   ! microphysics variables
   integer, parameter :: ncnstmax=4                      ! number of constituents
   character(len=8), dimension(ncnstmax), parameter :: & ! constituent names
     cnst_names = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)
   integer :: ncnst 					! number of constituents (can vary)
   character(len=16) :: microp_scheme  ! microphysics scheme 'RK' or 'MG'
   integer :: &
     ixcldliq,     &					! cloud liquid amount index for state%q
     ixcldice,     &					! cloud ice amount index
     ixnumliq,     &					! cloud liquid number index
     ixnumice						! cloud ice water index

   ! COSP-related local vars
   type(cosp_config) :: cfg   				! Configuration options
   type(cosp_gridbox) :: gbx 				! Gridbox information. Input for COSP
   type(cosp_subgrid) :: sgx     			! Subgrid outputs
   type(cosp_sgradar) :: sgradar 			! Output from radar simulator
   type(cosp_sglidar) :: sglidar 			! Output from lidar simulator
   type(cosp_isccp)   :: isccp   			! Output from ISCCP simulator
   type(cosp_misr)    :: misr    			! Output from MISR simulator
   type(cosp_vgrid)   :: vgrid   			! Information on vertical grid of stats
   type(cosp_radarstats) :: stradar 			! Summary statistics from radar simulator
   type(cosp_lidarstats) :: stlidar 			! Summary statistics from lidar simulator
   type(cosp_modis)   :: modis   			! Output from MODIS simulator (new in cosp v1.3)
   !!!type(cosp_rttov)   :: rttov   			! Output from RTTOV (not using)

   ! COSP input variables that depend on CAM
   ! 1) Npoints = number of gridpoints COSP will process (without subsetting, Npoints=ncol)
   ! 2) Nlevels = number of model levels (Nlevels=pver)
   real(r8), parameter :: time = 1.0_r8 		! time ! Time since start of run [days], set to 1 bc running over single CAM timestep
   real(r8), parameter :: time_bnds(2)=(/0.5_r8,1.5_r8/)	! time_bnds ! Time boundaries - new in cosp v1.3, set following cosp_test.f90 line 121
   integer :: Npoints                          		! Number of gridpoints COSP will process
   integer :: Nlevels	 				! Nlevels
   logical :: use_reff					! True if effective radius to be used by radar simulator 
							! (always used by lidar) (set based on CAM4 or CAM5)
   logical :: use_precipitation_fluxes			! True if precipitation fluxes are input to the algorithm 
							! (set based on CAM4 or CAM5)
   real(r8), parameter :: emsfc_lw = 0.99_r8   		! longwave emissivity of surface at 10.5 microns 
							! set value same as in cloudsimulator.F90

   ! local vars related to calculations to go from CAM input to COSP input
   ! cosp convective value includes both deep and shallow convection
   real(r8) :: ptop(pcols,pver)				! top interface pressure (Pa)
   real(r8) :: ztop(pcols,pver)				! top interface height asl (m)
   real(r8) :: lat_cosp(pcols)    			! lat for cosp (degrees_north)
   real(r8) :: lon_cosp(pcols)	    			! lon for cosp (degrees_east)
   real(r8) :: landmask(pcols)    			! landmask (0 or 1)
   real(r8) :: mr_ccliq(pcols,pver)			! mixing_ratio_convective_cloud_liquid (kg/kg)
   real(r8) :: mr_ccice(pcols,pver)			! mixing_ratio_convective_cloud_ice (kg/kg)
   real(r8) :: rain_cv(pcols,pverp)			! interface flux_convective_cloud_rain (kg m^-2 s^-1)
   real(r8) :: snow_cv(pcols,pverp)			! interface flux_convective_cloud_snow (kg m^-2 s^-1)
   real(r8) :: rain_cv_interp(pcols,pver)		! midpoint flux_convective_cloud_rain (kg m^-2 s^-1)
   real(r8) :: snow_cv_interp(pcols,pver)		! midpoint flux_convective_cloud_snow (kg m^-2 s^-1)
   real(r8) :: grpl_ls_interp(pcols,pver)		! midpoint ls grp flux, should be 0
   real(r8) :: rain_ls_interp(pcols,pver)		! midpoint ls rain flux (kg m^-2 s^-1)
   real(r8) :: snow_ls_interp(pcols,pver)		! midpoint ls snow flux
   real(r8) :: reff_cosp(pcols,pver,nhydro)		! effective radius for cosp input
   real(r8) :: rh(pcols,pver)    			! relative_humidity_liquid_water (%)
   real(r8) :: es(pcols,pver)    			! saturation vapor pressure
   real(r8) :: qs(pcols,pver)    			! saturation mixing ratio (kg/kg), saturation specific humidity
   real(r8) :: tau(pcols,pver)		   	   	! tau - optical depth
   real(r8) :: dtau_s(pcols,pver)  			! dtau_s - Optical depth of stratiform cloud at 0.67 um
   real(r8) :: dtau_c(pcols,pver)   			! dtau_s - Optical depth of convective cloud at 0.67 um
   real(r8) :: dem_s(pcols,pver)   			! dem_s - Longwave emis of stratiform cloud at 10.5 um
   real(r8) :: dem_c(pcols,pver)	   		! dem_c - Longwave emis of convective cloud at 10.5 um

   ! CAM pointers to get variables from radiation interface (get from rad_cnst_get_clim_gas)
   real(r8), pointer, dimension(:,:) :: q 	   	! specific humidity (kg/kg)
   real(r8), pointer, dimension(:,:) :: o3        	! Mass mixing ratio 03
   real(r8), pointer, dimension(:,:) :: co2 	   	! Mass mixing ratio C02
   real(r8), pointer, dimension(:,:) :: ch4        	! Mass mixing ratio CH4
   real(r8), pointer, dimension(:,:) :: n2o 	   	! Mass mixing ratio N20

   ! CAM pointers to get variables from the physics buffer
   real(r8), pointer, dimension(:,:) :: cld        	! cloud fraction, tca - total_cloud_amount (0-1)
   real(r8), pointer, dimension(:,:) :: concld     	! concld fraction, cca - convective_cloud_amount (0-1)
   real(r8), pointer, dimension(:,:) :: rel     	! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:) :: rei     	! ice effective drop size (microns)

   ! More pointers;  pbuff in zm_conv_intr.F90, calc in zm_conv.F90 
   real(r8), pointer, dimension(:,:) :: dp_flxprc   ! deep interface gbm flux_convective_cloud_rain (kg m^-2 s^-1)
   real(r8), pointer, dimension(:,:) :: dp_flxsnw   ! deep interface gbm flux_convective_cloud_snow (kg m^-2 s^-1) 

   ! More pointers;  pbuf in zm_conv_intr.F90, calc in zm_conv.F90, 0 for cam3.5.  value for CAM5?
   real(r8), pointer, dimension(:,:) :: dp_cldliq   ! deep gbm cloud liquid water (kg/kg)
   real(r8), pointer, dimension(:,:) :: dp_cldice   ! deep gmb cloud ice water (kg/kg)

   ! More pointers;  pbuf in convect_shallow.F90, calc in hk_conv.F90/convect_shallow.F90 (CAM4), uwshcu.F90 (CAM5) 
   real(r8), pointer, dimension(:,:) :: sh_flxprc   ! shallow interface gbm flux_convective_cloud_rain (kg m^-2 s^-1) 
   real(r8), pointer, dimension(:,:) :: sh_flxsnw   ! shallow interface gbm flux_convective_cloud_snow (kg m^-2 s^-1)

   ! More pointers;  pbuf in convect_shallow.F90, calc in hk_conv.F90 (should be 0!),uwshcu.F9 
   real(r8), pointer, dimension(:,:) :: sh_cldliq   ! shallow gbm cloud liquid water (kg/kg)
   real(r8), pointer, dimension(:,:) :: sh_cldice   ! shallow gbm cloud ice water (kg/kg)

   ! More pointers;  pbuf in stratiform.F90, getting from pbuf here
   ! a) added as output to pcond subroutine in cldwat.F90 and to nmicro_pcond subroutine in cldwat2m.F90
   real(r8), pointer, dimension(:,:) :: ls_flxprc   ! stratiform interface gbm flux_cloud_rain (kg m^-2 s^-1) 
   real(r8), pointer, dimension(:,:) :: ls_flxsnw   ! stratiform interface gbm flux_cloud_snow (kg m^-2 s^-1)
   ! b) added as output to nmicro_pcond subroutine in cldwat2m.F90
   real(r8), pointer, dimension(:,:) :: ls_mrprc 	! grid-box mean rain mixing ratio (kg/kg)
   real(r8), pointer, dimension(:,:) :: ls_mrsnw 	! grid-box mean rain mixing ratio (kg/kg)

   ! Output CAM variables
   ! Notes:
   ! 1) use pcols (maximum number of columns that code could use, maybe 16)
   ! pcols vs. ncol.  ncol is the number of columns a chunk is actually using, pcols is maximum number
   ! 2) Mixed variables rules/notes, need to collapse because CAM history does not support increased dimensionality
   ! MIXED DIMS: ntau_cosp*nprs_cosp, ndbze_cosp*nht_cosp, nsr_cosp*nht_cosp, nscol_cosp*nhtml_cosp, ntau_cosp*nhtmisr_cosp
   !	a) always making mixed variables VERTICAL*OTHER, e.g., pressure*tau or ht*dbze
   !	b) always collapsing output as V1_1/V2_1...V1_1/V2_N ; V1_2/V2_1 ...V1_2/V2_N etc. to V1_N/V2_1 ... V1_N/V2_N
   ! 	c) here, need vars for both multi-dimensional output from COSP, and two-dimensional output from CAM
   ! 3) ntime=1, nprofile=ncol
   ! 4) dimensions listed in COSP units are from netcdf output from cosp test case, and are not necessarily in the 
   !    correct order.  In fact, most of them are not as I discovered after trying to run COSP in-line.
   !    BE says this could be because FORTRAN and C (netcdf defaults to C) have different conventions.
   ! 5) !! Note: after running COSP, it looks like height_mlev is actually the model levels after all!!

   real(r8) :: clisccp2(pcols,ntau_cosp,nprs_cosp) 	! clisccp2 (time,tau,plev,profile)
   real(r8) :: cfad_dbze94(pcols,ndbze_cosp,nht_cosp)	! cfad_dbze94 (time,height,dbze,profile)
   real(r8) :: cfad_lidarsr532(pcols,nsr_cosp,nht_cosp)	! cfad_lidarsr532 (time,height,scat_ratio,profile)
   real(r8) :: dbze94(pcols,nscol_cosp,nhtml_cosp)	! dbze94 (time,height_mlev,column,profile)
   real(r8) :: atb532(pcols,nscol_cosp,nhtml_cosp)	! atb532 (time,height_mlev,column,profile)
   real(r8) :: clMISR(pcols,ntau_cosp,nhtmisr_cosp)	! clMISR (time,tau,CTH_height_bin,profile)
   real(r8) :: frac_out(pcols,nscol_cosp,nhtml_cosp)  	! frac_out (time,height_mlev,column,profile)

   real(r8) :: fisccp1(pcols,ntau_cosp*nprs_cosp) 	! CAM clisccp2 (time,tau,plev,profile)
   real(r8) :: cldtot_isccp(pcols) 			! CAM tclisccp (time,profile)
   real(r8) :: meancldalb_isccp(pcols)			! CAM albisccp (time,profile)
   real(r8) :: meanptop_isccp(pcols)			! CAM ctpisccp (time,profile)
   real(r8) :: cldlow_cal(pcols)			! CAM cllcalipso (time,profile)
   real(r8) :: cldmed_cal(pcols)			! CAM clmcalipso (time,profile)
   real(r8) :: cldhgh_cal(pcols)			! CAM clhcalipso (time,profile)
   real(r8) :: cldtot_cal(pcols)			! CAM cltcalipso (time,profile)
   real(r8) :: cld_cal(pcols,nht_cosp)			! CAM clcalipso (time,height,profile)
   real(r8) :: cfad_dbze94_cs(pcols,nht_cosp*ndbze_cosp)! CAM cfad_dbze94 (time,height,dbze,profile)
   real(r8) :: cfad_sr532_cal(pcols,nht_cosp*nsr_cosp)	! CAM cfad_lidarsr532 (time,height,scat_ratio,profile)
   real(r8) :: tau_isccp(pcols,nscol_cosp)		! CAM boxtauisccp (time,column,profile)
   real(r8) :: cldptop_isccp(pcols,nscol_cosp)		! CAM boxptopisccp (time,column,profile)
   real(r8) :: meantau_isccp(pcols)			! CAM tauisccp (time,profile)
   real(r8) :: meantb_isccp(pcols)			! CAM meantbisccp (time,profile)
   real(r8) :: meantbclr_isccp(pcols)			! CAM meantbclrisccp (time,profile)	
   real(r8) :: dbze_cs(pcols,nhtml_cosp*nscol_cosp)	! CAM dbze94 (time,height_mlev,column,profile)
   real(r8) :: cldtot_calcs(pcols)			! CAM cltlidarradar (time,profile)
   real(r8) :: cld_cal_notcs(pcols,nht_cosp)     	! CAM clcalipso2 (time,height,profile)
   real(r8) :: atb532_cal(pcols,nhtml_cosp*nscol_cosp)	! CAM atb532 (time,height_mlev,column,profile)
   real(r8) :: mol532_cal(pcols,nhtml_cosp)		! CAM beta_mol532 (time,height_mlev,profile)
   real(r8) :: cld_misr(pcols,nhtmisr_cosp*ntau_cosp)	! CAM clMISR (time,tau,CTH_height_bin,profile)
   real(r8) :: refl_parasol(pcols,nsza_cosp)		! CAM parasol_refl (time,sza,profile)
   real(r8) :: scops_out(pcols,nhtml_cosp*nscol_cosp)  	! CAM frac_out (time,height_mlev,column,profile)

   !! new modis output variables in COSP v1.3
   real(r8) :: cltmodis(pcols)
   real(r8) :: clwmodis(pcols)
   real(r8) :: climodis(pcols)
   real(r8) :: clhmodis(pcols)
   real(r8) :: clmmodis(pcols)
   real(r8) :: cllmodis(pcols)
   real(r8) :: tautmodis(pcols)
   real(r8) :: tauwmodis(pcols)
   real(r8) :: tauimodis(pcols)
   real(r8) :: tautlogmodis(pcols)
   real(r8) :: tauwlogmodis(pcols)
   real(r8) :: tauilogmodis(pcols)
   real(r8) :: reffclwmodis(pcols)
   real(r8) :: reffclimodis(pcols)
   real(r8) :: pctmodis(pcols)
   real(r8) :: lwpmodis(pcols)
   real(r8) :: iwpmodis(pcols)
   real(r8) :: clmodis_cam(pcols,ntau_cosp_modis*nprs_cosp)
   real(r8) :: clmodis(pcols,ntau_cosp_modis,nprs_cosp)

   type(interp_type)  :: interp_wgts
   integer, parameter :: extrap_method = 1  		! sets extrapolation method to boundary value (1)

  !---------------- End of declaration of variables --------------

   !! find the chunk and ncol from the state vector
   lchnk = state%lchnk   !! state variable contains a number of columns, one chunk
   ncol = state%ncol	 !! number of columns in the chunk

   ! initialize temporary variables as fillvalue - need to do this otherwise array expansion puts garbage in history
   ! file for columns over which COSP did make calculations.
   tmp(1:pcols)=fillvalue
   tmp1(1:pcols,1:pver)=fillvalue
   tmp2(1:pcols,1:pver)=fillvalue

   ! Initialize CAM variables as fillvalue, important for history files because it will exclude these from averages
   ! (multi-dimensional output that will be collapsed)
   ! initialize over all pcols, not just ncol.  missing values needed in chunks where ncol<pcols

   clisccp2(1:pcols,1:ntau_cosp,1:nprs_cosp)=fillvalue
   cfad_dbze94(1:pcols,1:ndbze_cosp,1:nht_cosp)=fillvalue
   cfad_lidarsr532(1:pcols,1:nsr_cosp,1:nht_cosp)=fillvalue
   dbze94(1:pcols,1:nscol_cosp,1:nhtml_cosp)=fillvalue
   atb532(1:pcols,1:nscol_cosp,1:nhtml_cosp)=fillvalue
   clMISR(1:pcols,ntau_cosp,1:nhtmisr_cosp)=fillvalue
   frac_out(1:pcols,1:nscol_cosp,1:nhtml_cosp)=fillvalue

   ! (all CAM output variables. including collapsed variables)
   fisccp1(1:pcols,1:ntau_cosp*nprs_cosp)=fillvalue
   cldtot_isccp(1:pcols)=fillvalue
   meancldalb_isccp(1:pcols)=fillvalue
   meanptop_isccp(1:pcols)=fillvalue
   cldlow_cal(1:pcols)=fillvalue
   cldmed_cal(1:pcols)=fillvalue
   cldhgh_cal(1:pcols)=fillvalue
   cldtot_cal(1:pcols)=fillvalue
   cld_cal(1:pcols,1:nht_cosp)=fillvalue
   cfad_dbze94_cs(1:pcols,1:nht_cosp*ndbze_cosp)=fillvalue
   cfad_sr532_cal(1:pcols,1:nht_cosp*nsr_cosp)=fillvalue
   tau_isccp(1:pcols,1:nscol_cosp)=fillvalue
   cldptop_isccp(1:pcols,1:nscol_cosp)=fillvalue
   meantau_isccp(1:pcols)=fillvalue
   meantb_isccp(1:pcols)=fillvalue
   meantbclr_isccp(1:pcols)=fillvalue	
   dbze_cs(1:pcols,1:nhtml_cosp*nscol_cosp)=fillvalue
   cldtot_calcs(1:pcols)=fillvalue
   cld_cal_notcs(1:pcols,1:nht_cosp)=fillvalue
   atb532_cal(1:pcols,1:nhtml_cosp*nscol_cosp)=fillvalue
   mol532_cal(1:pcols,1:nhtml_cosp)=fillvalue
   cld_misr(1:pcols,1:nhtmisr_cosp*ntau_cosp)=fillvalue
   refl_parasol(1:pcols,1:nsza_cosp)=fillvalue
   scops_out(1:pcols,1:nhtml_cosp*nscol_cosp)=fillvalue

   !! new modis output variables in COSP v1.3
   cltmodis(1:pcols)=fillvalue
   clwmodis(1:pcols)=fillvalue
   climodis(1:pcols)=fillvalue
   clhmodis(1:pcols)=fillvalue
   clmmodis(1:pcols)=fillvalue
   cllmodis(1:pcols)=fillvalue
   tautmodis(1:pcols)=fillvalue
   tauwmodis(1:pcols)=fillvalue
   tauimodis(1:pcols)=fillvalue
   tautlogmodis(1:pcols)=fillvalue
   tauwlogmodis(1:pcols)=fillvalue
   tauilogmodis(1:pcols)=fillvalue
   reffclwmodis(1:pcols)=fillvalue
   reffclimodis(1:pcols)=fillvalue
   pctmodis(1:pcols)=fillvalue
   lwpmodis(1:pcols)=fillvalue
   iwpmodis(1:pcols)=fillvalue
   clmodis_cam(1:pcols,1:ntau_cosp_modis*nprs_cosp)=fillvalue
   clmodis(1:pcols,1:ntau_cosp_modis,1:nprs_cosp)=fillvalue

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! POPULATE COSP CONFIGURATION INPUT VARIABLE ("cfg") FROM CAM NAMELIST
! Note: This is a first cut, and these values may get overwritten later.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Define/assign settings that would otherwise be read in using the COSP namelist.
! (e.g., used cosp_test.f90 to identify necessary inputs to COSP and run structure)

! SET COSP RUN OPTIONS
! values come from namelist within atm_in namelist read in by cospsimulator_intr_readnl
! note: some of these will be changed depending on the namelist variables later

   print *, 'Populating cfg for COSP...'
   !!! checked cfg% structure vs. cosp_types.f90 v1.3

   ! Simulator flags
   cfg%Lradar_sim=lradar_sim
   cfg%Llidar_sim=llidar_sim
   cfg%Lisccp_sim=lisccp_sim
   cfg%Lmisr_sim=lmisr_sim
   cfg%Lmodis_sim=lmodis_sim !! new in v1.3
   cfg%Lrttov_sim=lrttov_sim !! new in v1.3

   ! Output variable logicals
   ! Notes on name changes in cfg% from v1.1 to v1.3
   !! cfg%Lctpisccp to cfg%Lpctisccp
   !! cfg%Ltclisccp to cfg%Lcltisccp 
   !! cfg%Lcfad_dbze94 to cfg%Lcfaddbze94
   !! cfg%Lcfad_lidarsr532 to cfg%LcfadLidarsr532
   !! cfg%Lparasol_refl to cfg%LparasolRefl
   !! cfg%Lbeta_mol532 to cfg%LlidarBetaMol532
   !! cfg%lclisccp2 to cfg%Lclisccp

   cfg%Lalbisccp=lalbisccp
   cfg%Latb532=latb532
   cfg%Lboxptopisccp=lboxptopisccp
   cfg%Lboxtauisccp=lboxtauisccp
   cfg%Lcfaddbze94=lcfad_dbze94
   cfg%LcfadLidarsr532=lcfad_lidarsr532
   cfg%Lclcalipso=lclcalipso
   cfg%Lclhcalipso=lclhcalipso
   cfg%Lclisccp=lclisccp2
   cfg%Lcllcalipso=lcllcalipso
   cfg%Lclmcalipso=lclmcalipso
   cfg%Lcltcalipso=lcltcalipso
   cfg%Lpctisccp=lctpisccp
   cfg%Ldbze94=ldbze94
   cfg%Ltauisccp=ltauisccp
   cfg%Lcltisccp=ltclisccp
   cfg%Llongitude=llongitude
   cfg%Llatitude=llatitude
   cfg%LparasolRefl=lparasol_refl
   cfg%LclMISR=lclMISR
   cfg%Lmeantbisccp=lmeantbisccp
   cfg%Lmeantbclrisccp=lmeantbclrisccp
   cfg%Lclcalipso2=lclcalipso2
   cfg%Lcltlidarradar=lcltlidarradar
   cfg%Lfracout=lfrac_out
   cfg%LlidarBetaMol532=lbeta_mol532
   cfg%Lcltmodis=lcltmodis
   cfg%Lclwmodis=lclwmodis
   cfg%Lclimodis=lclimodis
   cfg%Lclhmodis=lclhmodis
   cfg%Lclmmodis=lclmmodis
   cfg%Lcllmodis=lcllmodis
   cfg%Ltautmodis=ltautmodis
   cfg%Ltauwmodis=ltauwmodis
   cfg%Ltauimodis=ltauimodis
   cfg%Ltautlogmodis=ltautlogmodis
   cfg%Ltauwlogmodis=ltauwlogmodis
   cfg%Ltauilogmodis=ltauilogmodis
   cfg%Lreffclwmodis=lreffclwmodis
   cfg%Lreffclimodis=lreffclimodis
   cfg%Lpctmodis=lpctmodis
   cfg%Llwpmodis=llwpmodis
   cfg%Liwpmodis=liwpmodis
   cfg%Lclmodis=lclmodis
   cfg%Ltbrttov=ltbrttov

   ! make COSP radar and lidar stats are calculated
   ! from line 1464 of cosp_io.f90 (cosp_io.f90 code not included here, so need specify this cfg% option)
   if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.

   cfg%Lwrite_output=.false.  !! COSP I/O not used, should always be .false.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! GET CAM GEOPHYSICAL VARIABLES NEEDED FOR COSP INPUT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! 1) state variables (prognostic variables, see physics_types.F90)
! state vars are passed to this subroutine from radiation.F90.   
! I do not need to define these variables.  I can use them as is, e.g., state%t
   !state%lat	! lat (radians) 
   !state%lon	! lon (radians) 
   !state%t	! temperature (K)
   !state%u  	! u_wind zonal wind (m/s)
   !state%v  	! v_wind meridional wind (m/s)
   !state%ps    ! surface pressure (Pa)
   !state%pint  ! p - p_in_full_levels (Pa)
   !state%pmid  ! ph - p_in_half_levels (Pa)
   !state%zm	! geopotential height above surface at midpoints (m), pver
   !state%zi	! geopotential height above surface at interfaces (m), pverp
   !state%q(:,:,ixcldliq) !cldliq or mr_lsliq - mixing_ratio_large_scale_cloud_liquid (kg/kg)
   !state%q(:,:,ixcldice) !cldice or mr_lsice -  mixing_ratio_large_scale_cloud_ice(kg/kg)

   ! need query index for cldliq and cldice
   ! use cnst_get_ind subroutine in constituents.F90.
   ! can also get MG microphysics number from state using similar procedure.
   call cnst_get_ind('CLDLIQ',ixcldliq)  !! replaced cnst_names(1) not setting abort flag which is optional in cnst_get_ind
   call cnst_get_ind(cnst_names(2),ixcldice)

   Npoints = ncol	 ! default is running all columns in the chunk, not pcols = maximum number
   Nlevels = pver

! 2) cam_in variables (see camsrfexch_types.F90)
! I can reference these as is, e.g., cam_in%ts.  
   !cam_in%ts		   	! skt - Skin temperature (K)
   !cam_in%landfrac	   	! land fraction, used to define a landmask (0 or 1) for COSP input

! 3) radiative constituent interface variables:
! specific humidity (q), 03, CH4,C02, N20 mass mixing ratio
! Note: these all have dims (pcol,pver) but the values don't change much for the well-mixed gases.
    call rad_cnst_get_clim_gas('H2O', state, pbuf, q)			
    call rad_cnst_get_clim_gas('O3',  state, pbuf, o3)
    call rad_cnst_get_clim_gas('CH4', state, pbuf, ch4)
    call rad_cnst_get_clim_gas('CO2', state, pbuf, co2)
    call rad_cnst_get_clim_gas('N2O', state, pbuf, n2o)

! 4) variables from variables physics buffer
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('CONCLD')
   concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('REL')
   rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('REI')
   rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)  

   ! Vars added to physics buffer in other interfaces (not radiation.F90)
   ! all "1" at the end ok as is because radiation/intr after when these were added to physics buffer
   ifld = pbuf_get_fld_idx('DP_FLXPRC')
   dp_flxprc  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
   ifld = pbuf_get_fld_idx('DP_FLXSNW')
   dp_flxsnw  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
   ifld = pbuf_get_fld_idx('DP_CLDLIQ')
   dp_cldliq  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('DP_CLDICE')
   dp_cldice  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('SH_FLXPRC')
   sh_flxprc  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
   ifld = pbuf_get_fld_idx('SH_FLXSNW')
   sh_flxsnw  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
   ifld = pbuf_get_fld_idx('SH_CLDLIQ')
   sh_cldliq  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('SH_CLDICE')
   sh_cldice  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('LS_FLXPRC')
   ls_flxprc  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
   ifld = pbuf_get_fld_idx('LS_FLXSNW')
   ls_flxsnw  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

   ! routine used to identify which microphysics scheme is being used.
   ! call phys_getopts(microp_scheme_out=microp_scheme) 	
   ! if using CAM5 - then mg has mixing ratios for prc and snw, preferable
   ! if ( microp_scheme .eq. 'MG' ) then
   !	ifld = pbuf_get_fld_idx('LS_MRPRC')
   !	ls_mrprc  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   !	ifld = pbuf_get_fld_idx('LS_MRSNW')
   !	ls_mrsnw  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ! end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CALCULATE COSP INPUT VARIABLES FROM CAM VARIABLES, done for all columns within chunk
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! 0)create ptop/ztop for gbx%pf and gbx%zlev are for the the interface, 
!  also reverse CAM height/pressure values for input into CSOP
!  CAM state%pint from top to surface, COSP wants surface to top.

   ! initalize
   ptop(1:ncol,1:pver)=0._r8
   ztop(1:ncol,1:pver)=0._r8

   ! assign values from top
   do k=1,pverp-1
	ptop(1:ncol,k)=state%pint(1:ncol,pverp-k)
	ztop(1:ncol,k)=state%zi(1:ncol,pverp-k)
   end do

! 1) lat/lon - convert from radians to cosp input type

   ! initalize lat and lons
   lat_cosp(1:ncol)=0._r8
   lon_cosp(1:ncol)=0._r8
   ! convert from radians to degrees_north and degrees_east 
   lat_cosp=state%lat*180._r8/(pi)  ! needs to go from -90 to +90 degrees north
   lon_cosp=state%lon*180._r8/(pi)  ! needs to go from 0 to 360 degrees east

! 2) rh - relative_humidity_liquid_water (%)

   ! calculate from CAM q and t using CAM built-in functions
   call aqsat_water(state%t, state%pmid, es, qs, pcols, &
	ncol, pver, 1, pver)

   ! initialize rh
   rh(1:ncol,1:pver)=0._r8

   ! calculate rh
   do k=1,pver
   do i=1,ncol
	rh(i,k)=(q(i,k)/qs(i,k))*100
   end do
   end do

! 3) landmask - calculate from cam_in%landfrac

   ! initalize landmask
   landmask(1:ncol)=0._r8
   ! calculate landmask
   do i=1,ncol
	if (cam_in%landfrac(i).gt.0.01_r8) landmask(i)= 1
   end do

! 4) calculate necessary input cloud/precip variables
! CAM4 note: don't take the cloud water from the hack shallow convection scheme or the deep convection.  
! cloud water values for convection are the same as the stratiform value. (Sungsu)
! all precip fluxes are mid points, all values are grid-box mean ("gbm") (Yuying)
! Note: COSP complains when values are below 0, CAM sometimes has rain mixing values below zero.  Is this a problem? ##2check##

   ! initialize local variables
   mr_ccliq(1:ncol,1:pver)=0._r8
   mr_ccice(1:ncol,1:pver)=0._r8
   grpl_ls_interp(1:ncol,1:pver)=0._r8
   rain_ls_interp(1:ncol,1:pver)=0._r8 
   snow_ls_interp(1:ncol,1:pver)=0._r8
   rain_cv(1:ncol,1:pverp)=0._r8
   snow_cv(1:ncol,1:pverp)=0._r8
   rain_cv_interp(1:ncol,1:pver)=0._r8
   snow_cv_interp(1:ncol,1:pver)=0._r8

   reff_cosp(1:ncol,1:pver,1:nhydro)=0._r8	
   ! note: reff_cosp dimensions should be same as cosp (reff_cosp has 9 hydrometeor dimension)
   ! Reff(Npoints,Nlevels,N_HYDRO)

   ! add together deep and shallow convection precipitation
   rain_cv(1:ncol,1:pverp) = sh_flxprc(1:ncol,1:pverp) + dp_flxprc(1:ncol,1:pverp)
   snow_cv(1:ncol,1:pverp) = sh_flxsnw(1:ncol,1:pverp) + dp_flxsnw(1:ncol,1:pverp)

   ! interpolate interface precip fluxes to mid points
   ! rain_cv_interp (pver) from rain_cv (pverp), snow_cv_interp from snow_cv, 
   ! rain_ls_interp from ls_flxprc, snow_ls_interp from ls_flxsnw
   ! interpolate, loop over columns (can also put ":" but this gives less information)
   ! use interpolate_data.F90 in atm/cam/src/control
   do i=1,ncol
	   ! find weights (pressure weighting?)
	   call lininterp_init(state%zi(i,1:pverp),pverp,state%zm(i,1:pver),pver,extrap_method,interp_wgts)
	   ! interpolate  lininterp1d(arrin, nin, arrout, nout, interp_wgts)
	   ! note: lininterp is an interface, contains lininterp1d -- code figures out to use lininterp1d.
	   call lininterp(rain_cv(i,1:pverp),pverp,rain_cv_interp(i,1:pver),pver,interp_wgts)
	   call lininterp(snow_cv(i,1:pverp),pverp,snow_cv_interp(i,1:pver),pver,interp_wgts)
	   call lininterp(ls_flxprc(i,1:pverp),pverp,rain_ls_interp(i,1:pver),pver,interp_wgts)
	   call lininterp(ls_flxsnw(i,1:pverp),pverp,snow_ls_interp(i,1:pver),pver,interp_wgts)
	   call lininterp_finish(interp_wgts)
   end do

   ! if statement here to populate cloud/precip variables depending on model version

   ! currently just use CAM4
   use_precipitation_fluxes = .true.
   use_reff = .false.

   do k=1,pver
   do i=1,ncol
	if (cld(i,k) .gt. 0._r8) then
	   mr_ccliq(i,k) = (state%q(i,k,ixcldliq)/cld(i,k))*concld(i,k)
	   mr_ccice(i,k) = (state%q(i,k,ixcldice)/cld(i,k))*concld(i,k)
	else
	   mr_ccliq(i,k) = 0._r8
	   mr_ccice(i,k) = 0._r8
	end if
   end do
   end do

   !if ( CAM4 ) then !! do this when the trunk contains both track1 and track5
   !end if 
   !if ( CAM5 ) then   !!!
   !	mr_ccliq = sh_cldliq + dp_cldliq
   !	mr_ccice = sh_cldice + dp_cldice
   !	use_precipitation_fluxes = .false.
   ! - effective radius, include both cloud and rain, variable likely does not exist (REI or REL are available)
   !   use_reff = .true.
   !   reff_cosp(1:ncol,1:pver,0)=rel?
   !   reff_cosp(1:ncol,1:pver,1)=rei?
   !  future use size distribution parameters from the MG microphysics as inputs to COSP.
   !end if

! 5) calculate/assign optical depths and emissivities needed for isccp simulator
   ! currently just use CAM4
   ! compute the optical depth (code from radcswmx) -- jek got this code from from cloudsimulator.F90
   do k=1,pver
   do i=1,ncol
	   if(cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
	      tau(i,k) = (abarl + bbarl/rel(i,k)) * cliqwp(i,k) + &
		   (abari + bbari/rei(i,k)) * cicewp(i,k)
	   else
	      tau(i,k) = 0._r8
	   end if
   end do
   end do

   ! CAM4 assumes same radiative properties for stratiform and convective clouds, 
   ! see ISCCP_CLOUD_TYPES subroutine call in cloudsimulator.F90

   ! initialize
   dtau_s(1:ncol,1:pver) = 0._r8
   dtau_c(1:ncol,1:pver) = 0._r8
   dem_s(1:ncol,1:pver) = 0._r8 
   dem_c(1:ncol,1:pver) = 0._r8

   ! assign values
   dtau_s(1:ncol,1:pver) = tau(1:ncol,1:pver)
   dtau_c(1:ncol,1:pver) = tau(1:ncol,1:pver)
   dem_s(1:ncol,1:pver) = emis(1:ncol,1:pver) 
   dem_c(1:ncol,1:pver) = emis(1:ncol,1:pver)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!! END TRANSLATE CAM VARIABLES TO COSP INPUT VARIABLES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! RUN ALL COSP SIMULATORS ON ALL COLUMNS (COSP SETUP/CALL #1)
! This code is invoked when namelist variable "cosp_runall" is set to .true.
! I want to keep this functionality because it is useful for testing
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (cosp_runall) then

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  POPULATE COSP INPUT VARIABLE ("gbx") FROM CAM VARIABLES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Create and Populate the input structure for COSP - "gbx"
! Check cosp_input_nl information for gbx
! contruct_* and free_* routines from /home/jenkay/cosp/cosp.v1.1/cosp_types.f90 
! allocate and deallocate memory within this module.
! gbx %s to set are lon, lat, p, ph, zlev, zlev_half, T, rh, sh, cca, tca, ph, 
   ! skt, landmask, mr_ozone, u_wind, v_wind, sunlit
   ! mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, fl_lsrain, fl_lssnow, fl_lsgrpl, 
   ! fl_ccrain, fl_ccsnow, Reff, dtau_s, dtau_c, dem_s, dem_c 
   ! Make sure that you get cloud information from deep, shallow, and stratiform physics packages.

   print *, 'Allocating memory for gridbox type...'
   call construct_cosp_gridbox(time, &			        ! 1 double precision = real(r8) X
				time_bnds, &			! 1 double precision = real(r8)
				radar_freq, &			! 2 real(r8) X
				surface_radar, &		! 3 integer X
				use_mie_tables, &		! 4 integer X
				use_gas_abs, &			! 5 integer X
				do_ray, &			! 6 integer X 
				melt_lay, &			! 7 integer X 
				k2, &				! 8 real(r8) X
                                Npoints, &			! 9 integer X, same as CAM's ncol
				Nlevels, &			! 10 integer X
				ncolumns,&			! 11 integer X
				nhydro,&			! 12 integer X
				Nprmts_max_hydro,&      	! 13 integer X
				Naero,&				! 14 integer X
				Nprmts_max_aero,&		! 15 integer X
				Npoints_it, &			! 16 integer X
                                lidar_ice_type,&		! 17 integer X
				isccp_topheight,&		! 18 integer X
				isccp_topheight_direction,& 	! 19 integer X
				overlap,&			! 20 integer X
				emsfc_lw, &			! 21 real X
                                use_precipitation_fluxes,& 	! 22 logical X
				use_reff, &			! 23 logical X
                                Platform, &			! 24 integer X
				Satellite, &			! 25 integer X
				Instrument, &			! 26 integer X
				Nchannels, &			! 27 integer X
				ZenAng, &			! 28 real(r8) X
                                Channels(1:Nchannels),&		! 29 integer X
				Surfem(1:Nchannels),&		! 30 real(r8) X
				co2(1,1),&			! 31 real(r8) X
				ch4(1,1),&			! 32 real(r8) X
				n2o(1,1),&			! 33 real(r8) X
				co,&				! 34 real(r8) X
				gbx)				! OUT
    
   print *, 'Populating input structure for COSP...'
   ! Note: GBX expects vertical ordering to be from SURFACE(1) to the TOP(nlev), while
   ! by default CAM uses TOP(1) to SURFACE(nlev).
   ! also gbx is by definition of size ncol, but many variables defined as pcol.  be explicit here, gbx = 1:ncol

   gbx%longitude = lon_cosp(1:ncol)  				! lon (degrees_east)
   gbx%latitude = lat_cosp(1:ncol)   				! lat (degrees_north)
   gbx%p = ptop(1:ncol,1:pver) 					! p_in_full_levels (Pa),upper interface (flipped above)
   gbx%ph = state%pmid(1:ncol,pver:1:-1)			! p_in_half_levels (Pa)
   gbx%zlev = ztop(1:ncol,1:pver)				! height_in_full_levels (m),upper interface,asl (flipped above)
   gbx%zlev_half = state%zm(1:ncol,pver:1:-1) 			! height_in_half_levels (m),asl
   gbx%T = state%t(1:ncol,pver:1:-1)				! air_temperature (K)
   gbx%q = rh(1:ncol,pver:1:-1)					! relative humidity wrt water 
								! note: it is confusing that it is called "q" within cosp,
								! but it is rh according to Yuying
   gbx%sh = q(1:ncol,pver:1:-1)					! specific humidity (kg/kg), range [0,0.019ish]
   gbx%cca = concld(1:ncol,pver:1:-1)				! convective_cloud_amount (0-1)
   gbx%tca = cld(1:ncol,pver:1:-1)				! total_cloud_amount (0-1)
   gbx%psfc = state%ps(1:ncol)					! surface pressure (Pa)
   gbx%skt  = cam_in%ts(1:ncol)					! skin temperature (K)
   gbx%land = landmask(1:ncol)					! landmask (0 or 1)
   gbx%mr_ozone  = o3(1:ncol,pver:1:-1)				! ozone mass mixing ratio (kg/kg)
   gbx%u_wind  = state%u(1:ncol,pver)				! surface u_wind (m/s)
   gbx%v_wind  = state%v(1:ncol,pver)				! surface v_wind (m/s)
   gbx%sunlit  = 1						!  ??? what does this mean??
   gbx%mr_hydro(:,:,I_LSCLIQ) = state%q(1:ncol,pver:1:-1,ixcldliq)	! mr_lsliq, mixing_ratio_large_scale_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_LSCICE) = state%q(1:ncol,pver:1:-1,ixcldice)	! mr_lsice - mixing_ratio_large_scale_cloud_ice (kg/kg)
   gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq(1:ncol,pver:1:-1)		! mr_ccliq - mixing_ratio_convective_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice(1:ncol,pver:1:-1)		! mr_ccice - mixing_ratio_convective_cloud_ice (kg/kg)
   gbx%rain_ls = rain_ls_interp(1:ncol,pver:1:-1)			! midpoint gbm
   gbx%snow_ls = snow_ls_interp(1:ncol,pver:1:-1)			! midpoint gbm			
   gbx%grpl_ls = grpl_ls_interp(1:ncol,pver:1:-1)			! midpoint gbm
   gbx%rain_cv = rain_cv_interp(1:ncol,pver:1:-1)			! midpoint gbm
   gbx%snow_cv = snow_cv_interp(1:ncol,pver:1:-1)			! midpoint gbm
   gbx%Reff = reff_cosp(1:ncol,pver:1:-1,:)
   gbx%dtau_s   = dtau_s(1:ncol,pver:1:-1)
   gbx%dtau_c   = dtau_c(1:ncol,pver:1:-1)
   gbx%dem_s    = dem_s(1:ncol,pver:1:-1)
   gbx%dem_c    = dem_c(1:ncol,pver:1:-1)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! STARTUP RELATED TO COSP OUTPUT (see cosp.F90)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Note: COSP output variables are sgx (sub-grid outputs), sgradar (radar outputs), sglidar (lidar outputs), 
! isccp (isccp outputs), misr (misr simulator outputs), vgrid (vertical grid info), stradar 
! (summary statistics radar simulator), stlidar (summary statistics lidar simulator)

! Define new vertical grid
   print *, 'Defining new vertical grid...'
   call construct_cosp_vgrid(gbx,nlr,use_vgrid,csat_vgrid,vgrid)
        
! Allocate memory for output
   print *, 'Allocating memory for other types...'
   call construct_cosp_subgrid(Npoints, ncolumns, Nlevels, sgx)
   call construct_cosp_sgradar(cfg,Npoints,ncolumns,Nlevels,nhydro,sgradar)
   call construct_cosp_radarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,stradar)
   call construct_cosp_sglidar(cfg,Npoints,ncolumns,Nlevels,nhydro,PARASOL_NREFL,sglidar)
   call construct_cosp_lidarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,PARASOL_NREFL,stlidar)
   call construct_cosp_isccp(cfg,Npoints,ncolumns,Nlevels,isccp)
   call construct_cosp_misr(cfg,Npoints,misr)
   !! new for cosp v1.3
   call construct_cosp_modis(cfg,Npoints,modis)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  RUN COSP on all columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   print *, 'Calling the COSP simulator and running on all columns...'
   call cosp(overlap,ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  TRANSLATE COSP OUTPUT INTO INDIVIDUAL VARIABLES FOR OUTPUT (see nc_write_cosp_1d in cosp_io.f90)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Populate CAM vars with COSP output (the reasoning for this is based on cosp_io.f90)
   ! (1d variables)
   cldlow_cal(1:ncol) = stlidar%cldlayer(:,1)		! CAM version of cllcalipso (time,profile)	
   cldmed_cal(1:ncol) = stlidar%cldlayer(:,2)		! CAM version of clmcalipso (time,profile)
   cldhgh_cal(1:ncol) = stlidar%cldlayer(:,3)		! CAM version of clhcalipso (time,profile)
   cldtot_cal(1:ncol) = stlidar%cldlayer(:,4)		! CAM version of cltcalipso (time,profile)
   cldtot_calcs(1:ncol) = stradar%radar_lidar_tcc	! CAM version of cltlidarradar (time,profile)
   cldtot_isccp(1:ncol) = isccp%totalcldarea		! CAM version of tclisccp (time,profile)
   meanptop_isccp(1:ncol) = isccp%meanptop		! CAM version of ctpisccp (time,profile)
   meantau_isccp(1:ncol) = isccp%meantaucld		! CAM version of meantbisccp (time,profile)
   meancldalb_isccp(1:ncol) = isccp%meanalbedocld	! CAM version of albisccp (time,profile)
   meantb_isccp(1:ncol) = isccp%meantb			! CAM version of meantbisccp (time,profile)
   meantbclr_isccp(1:ncol) = isccp%meantbclr		! CAM version of meantbclrisccp (time,profile)
   ! (2d variables)
   cld_cal(1:ncol,1:nht_cosp) = stlidar%lidarcld			! CAM version of clcalipso (time,height,profile)
   cld_cal_notcs(1:ncol,1:nht_cosp) = stradar%lidar_only_freq_cloud	! CAM version of clcalipso2 (time,height,profile)
   mol532_cal(1:ncol,1:nhtml_cosp) = sglidar%beta_mol			! CAM version of beta_mol532 (time,height_mlev,profile)
   tau_isccp(1:ncol,1:nscol_cosp) = isccp%boxtau			! CAM version of boxtauisccp (time,column,profile)
   cldptop_isccp(1:ncol,1:nscol_cosp) = isccp%boxptop			! CAM version of boxptopisccp (time,column,profile)
   refl_parasol(1:ncol,1:nsza_cosp) = stlidar%parasolrefl		! CAM version of parasolrefl (time,sza,profile)
   ! (3d variables)
   dbze94(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sgradar%Ze_tot		! dbze94 (time,height_mlev,column,profile)
   atb532(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sglidar%beta_tot		! atb532 (time,height_mlev,column,profile)
   frac_out(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sgx%frac_out		! frac_out (time,height_mlev,column,profile)
   cfad_dbze94(1:ncol,1:ndbze_cosp,1:nht_cosp) = stradar%cfad_ze	! cfad_dbze94 (time,height,dbze,profile)
   cfad_lidarsr532(1:ncol,1:nsr_cosp,1:nht_cosp)= stlidar%cfad_sr	! cfad_lidarsr532 (time,height,scat_ratio,profile)
   clisccp2(1:ncol,1:ntau_cosp,1:nprs_cosp) = isccp%fq_isccp		! clisccp2 (time,tau,plev,profile)
   clMISR(1:ncol,1:ntau_cosp,1:nhtmisr_cosp) = misr%fq_MISR		! clMISR (time,tau,CTH_height_bin,profile)
   !! modis variables
   cltmodis(1:ncol)=modis%Cloud_Fraction_Total_Mean
   clwmodis(1:ncol)=modis%Cloud_Fraction_Water_Mean
   climodis(1:ncol)=modis%Cloud_Fraction_Ice_Mean
   clhmodis(1:ncol)=modis%Cloud_Fraction_High_Mean
   clmmodis(1:ncol)=modis%Cloud_Fraction_Mid_Mean
   cllmodis(1:ncol)=modis%Cloud_Fraction_Low_Mean
   tautmodis(1:ncol)=modis%Optical_Thickness_Total_Mean
   tauwmodis(1:ncol)=modis%Optical_Thickness_Water_Mean
   tauimodis(1:ncol)=modis%Optical_Thickness_Ice_Mean
   tautlogmodis(1:ncol)=modis%Optical_Thickness_Total_LogMean
   tauwlogmodis(1:ncol)=modis%Optical_Thickness_Water_LogMean
   tauilogmodis(1:ncol)=modis%Optical_Thickness_Ice_LogMean
   reffclwmodis(1:ncol)=modis%Cloud_Particle_Size_Water_Mean
   reffclimodis(1:ncol)=modis%Cloud_Particle_Size_Ice_Mean
   pctmodis(1:ncol)=modis%Cloud_Top_Pressure_Total_Mean
   lwpmodis(1:ncol)=modis%Liquid_Water_Path_Mean
   iwpmodis(1:ncol)=modis%Ice_Water_Path_Mean
   clmodis(1:ncol,1:ntau_cosp_modis,1:nprs_cosp)=modis%Optical_Thickness_vs_Cloud_Top_Pressure

   print *, 'Done mapping COSP output to local variables ...'

   ! Use high-dimensional output to populate CAM collapsed output variables
   ! see cosp_share.F90 for mixed dimension definitions
   ! i am using the convention of starting vertical coordinates at the surface, up to down, COSP convention, not CAM.

   do i=1,ncol
	! CAM clisccp2 (time,tau,plev,profile)
	do ip=1,nprs_cosp
	do it=1,ntau_cosp
	   ipt=(ip-1)*ntau_cosp+it
	   fisccp1(i,ipt) = clisccp2(i,it,ip)   		! fisccp1(pcols,ntau_cosp*nprs_cosp) 
	end do
	end do
	! CAM cfad_dbze94 (time,height,dbze,profile) 
	do ih=1,nht_cosp
	do id=1,ndbze_cosp
	   ihd=(ih-1)*ndbze_cosp+id			
	   cfad_dbze94_cs(i,ihd) = cfad_dbze94(i,id,ih)  	! cfad_dbze94_cs(pcols,nht_cosp*ndbze_cosp)
	end do
	end do
	! CAM cfad_lidarsr532 (time,height,scat_ratio,profile)
	do ih=1,nht_cosp
	do is=1,nsr_cosp
	   ihs=(ih-1)*nsr_cosp+is			
	   cfad_sr532_cal(i,ihs) = cfad_lidarsr532(i,is,ih)  ! cfad_sr532_cal(pcols,nht_cosp*nsr_cosp)
	end do
	end do
	! CAM dbze94 (time,height_mlev,column,profile)
	do ihml=1,nhtml_cosp
	do isc=1,nscol_cosp
	   ihsc=(ihml-1)*nscol_cosp+isc			
	   dbze_cs(i,ihsc) = dbze94(i,isc,ihml)   		! dbze_cs(pcols,pver*nscol_cosp) 
	end do
	end do
	! CAM clMISR (time,tau,CTH_height_bin,profile)
	do ihm=1,nhtmisr_cosp
	do it=1,ntau_cosp
	   ihmt=(ihm-1)*ntau_cosp+it			
	   cld_misr(i,ihmt) = clMISR(i,it,ihm)   		! cld_misr(pcols,nhtmisr_cosp*ntau_cosp)
	end do
	end do
	! CAM atb532 (time,height_mlev,column,profile)  FIX
	do ihml=1,nhtml_cosp
	do isc=1,nscol_cosp
	   ihsc=(ihml-1)*nscol_cosp+isc			
	   atb532_cal(i,ihsc) = atb532(i,isc,ihml)  	 	! atb532_cal(pcols,nht_cosp*nscol_cosp)
	end do
	end do
	! CAM frac_out (time,height_mlev,column,profile)  FIX
	do ihml=1,nhtml_cosp
	do isc=1,nscol_cosp
	   ihsc=(ihml-1)*nscol_cosp+isc			
	   scops_out(i,ihsc) = frac_out(i,isc,ihml)    		! scops_out(pcols,nht_cosp*nscol_cosp)
	end do
	end do
	! CAM clmodis
	do ip=1,nprs_cosp
	do it=1,ntau_cosp_modis
	   ipt=(ip-1)*ntau_cosp_modis+it
	   clmodis_cam(i,ipt) = clmodis(i,it,ip)
	end do
	end do
   end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DEALLOCATE MEMORY for running with all columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call free_cosp_gridbox(gbx)
   call free_cosp_vgrid(vgrid)
   call free_cosp_subgrid(sgx)
   call free_cosp_sgradar(sgradar)
   call free_cosp_radarstats(stradar)
   call free_cosp_sglidar(sglidar)
   call free_cosp_lidarstats(stlidar)
   call free_cosp_isccp(isccp)
   call free_cosp_misr(misr) 
   call free_cosp_modis(modis)

end if !! if cosp_runall

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  RUN COSP ONLY ON SUNLIT COLUMNS FOR ISCCP and/or MISR and/or MODIS SIMULATORS (COSP SETUP/CALL #2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ((.not.cosp_runall) .and. (lisccp_sim.or.lmisr_sim.or.lmodis_sim)) then  !! RUNNING ON SUNLIT

   cfg%Lalbisccp=lalbisccp
   cfg%Latb532=latb532
   cfg%Lboxptopisccp=lboxptopisccp
   cfg%Lboxtauisccp=lboxtauisccp
   cfg%Lcfaddbze94=lcfad_dbze94
   cfg%LcfadLidarsr532=lcfad_lidarsr532
   cfg%Lclcalipso=lclcalipso
   cfg%Lclhcalipso=lclhcalipso
   cfg%Lclisccp=lclisccp2
   cfg%Lcllcalipso=lcllcalipso
   cfg%Lclmcalipso=lclmcalipso
   cfg%Lcltcalipso=lcltcalipso
   cfg%Lpctisccp=lctpisccp
   cfg%Ldbze94=ldbze94
   cfg%Ltauisccp=ltauisccp
   cfg%Lcltisccp=ltclisccp
   cfg%Llongitude=llongitude
   cfg%Llatitude=llatitude
   cfg%LparasolRefl=lparasol_refl
   cfg%LclMISR=lclMISR
   cfg%Lmeantbisccp=lmeantbisccp
   cfg%Lmeantbclrisccp=lmeantbclrisccp
   cfg%Lclcalipso2=lclcalipso2
   cfg%Lcltlidarradar=lcltlidarradar
   cfg%Lfracout=lfrac_out
   cfg%LlidarBetaMol532=lbeta_mol532

! Only run the isccp and misr simulators
if (lisccp_sim.or.lmodis_sim) then
   cfg%Lisccp_sim= .true.
   cfg%Lalbisccp= .true.
   cfg%Latb532= .true.
   cfg%Lboxptopisccp= .true.
   cfg%Lboxtauisccp= .true.
   cfg%Lpctisccp= .true.
   cfg%Ltauisccp= .true.
   cfg%Lcltisccp= .true.
   cfg%Lmeantbisccp= .true.
   cfg%Lmeantbclrisccp= .true.
end if
if (lmisr_sim) then
   cfg%Lmisr_sim= .true.
   cfg%LclMISR= .true.
end if
if (lmodis_sim) then
   cfg%Lmodis_sim= .true.
   cfg%Lcltmodis= .true.
   cfg%Lclwmodis= .true.
   cfg%Lclimodis= .true.
   cfg%Lclhmodis= .true.
   cfg%Lclmmodis= .true.
   cfg%Lcllmodis= .true.
   cfg%Ltautmodis= .true.
   cfg%Ltauwmodis= .true.
   cfg%Ltauimodis= .true.
   cfg%Ltautlogmodis= .true.
   cfg%Ltauwlogmodis= .true.
   cfg%Ltauilogmodis= .true.
   cfg%Lreffclwmodis= .true.
   cfg%Lreffclimodis= .true.
   cfg%Lpctmodis= .true.
   cfg%Llwpmodis= .true.
   cfg%Liwpmodis= .true.
   cfg%Lclmodis= .true.
end if

! Don't run the radar/lidar simulators by setting COSP input structure cfg
cfg%Lradar_sim= .false.
cfg%Llidar_sim= .false.

! Also set radar/lidar/modis simulator outputs to false in COSP input structure cfg
cfg%Latb532= .false.
cfg%Lcfaddbze94= .false.
cfg%LcfadLidarsr532= .false.
cfg%Lclcalipso= .false.
cfg%Lclhcalipso= .false.
cfg%Lcllcalipso= .false.
cfg%Lclmcalipso= .false.
cfg%Lcltcalipso= .false.
cfg%Ldbze94= .false.
cfg%LparasolRefl= .false.
cfg%Lclcalipso2= .false.
cfg%Lcltlidarradar= .false.
cfg%LlidarBetaMol532= .false.

! Don't run cosp_stats.f90, only for radar and lidar
cfg%Lstats = .false.

! Subsetting using the TOD restriction
    ! Gather night/day column indices. (code from radsw.F90)
    Nday = 0
    Nno = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
          Nday = Nday + 1
          IdxDay(Nday) = i
       else
          Nno = Nno + 1
          IdxNo(Nno) = i
       end if
    end do

! Reset the number of points for input into COSP
Npoints = Nday

if (Nday .gt. 0) then  !! only run COSP if there are daylight points

!!! rearrange input arrays - compress arrays so that COSP only runs on sunlit columns
!e.g., a la
!   call CmpDayNite(E_pmid, pmid,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
! (function from ..atm/cam/src/physics/cam/cmparray_mod.F90, see also radsw.F90)
!! CmpDayNite_2d_R_Copy(InArray, OutArray...

call CmpDayNite(lon_cosp,lon_cosp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(lat_cosp,lat_cosp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(ptop, ptop_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%pmid, pmid_day, 	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(ztop,ztop_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%zm, zmid_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%t,t_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rh, rh_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(q, q_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(concld, concld_day, 	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(cld, cld_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%ps ,ps_day, 	Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(cam_in%ts ,ts_day, 	Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(landmask,landmask_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(o3, o3_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%u(1:ncol,pver),us_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(state%v(1:ncol,pver),vs_day,		Nday, IdxDay, Nno, IdxNo, 1, pcols)
call CmpDayNite(state%q(:,:,ixcldliq),cldliq_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%q(:,:,ixcldice),cldice_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(mr_ccliq, mr_ccliq_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(mr_ccice, mr_ccice_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rain_ls_interp, rain_ls_interp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(snow_ls_interp, snow_ls_interp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(grpl_ls_interp, grpl_ls_interp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rain_cv_interp,rain_cv_interp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)			
call CmpDayNite(snow_cv_interp,snow_cv_interp_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dtau_s,dtau_s_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dtau_c,dtau_c_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dem_s,dem_s_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dem_c,dem_c_day,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)

!! need to loop over nhydro for reff_cosp_day
do i=1,nhydro
   call CmpDayNite(reff_cosp(:,:,i),tmp1,		Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
   reff_cosp_day(:,:,i)=tmp1
   tmp1(1:pcols,1:pver)=fillvalue
end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  POPULATE COSP INPUT VARIABLE ("gbx") FROM CAM VARIABLES for sunlit columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   print *, 'Allocating memory for gridbox type for isccp sub-columns only...'
   call construct_cosp_gridbox(time, &			        ! 1 double precision = real(r8) X
				time_bnds, &			! 1 double precision = real(r8)
				radar_freq, &			! 2 real(r8) X
				surface_radar, &		! 3 integer X
				use_mie_tables, &		! 4 integer X
				use_gas_abs, &			! 5 integer X
				do_ray, &			! 6 integer X 
				melt_lay, &			! 7 integer X 
				k2, &				! 8 real(r8) X
                                Npoints, &			! 9 integer X
				Nlevels, &			! 10 integer X
				ncolumns,&			! 11 integer X
				nhydro,&			! 12 integer X
				Nprmts_max_hydro,&      	! 13 integer X
				Naero,&				! 14 integer X
				Nprmts_max_aero,&		! 15 integer X
				Npoints_it, &			! 16 integer X
                                lidar_ice_type,&		! 17 integer X
				isccp_topheight,&		! 18 integer X
				isccp_topheight_direction,& 	! 19 integer X
				overlap,&			! 20 integer X
				emsfc_lw, &			! 21 real X
                                use_precipitation_fluxes,& 	! 22 logical X
				use_reff, &			! 23 logical X
                                Platform, &			! 24 integer X
				Satellite, &			! 25 integer X
				Instrument, &			! 26 integer X
				Nchannels, &			! 27 integer X
				ZenAng, &			! 28 real(r8) X
                                Channels(1:Nchannels),&		! 29 integer X
				Surfem(1:Nchannels),&		! 30 real(r8) X
				co2(1,1),&			! 31 real(r8) X
				ch4(1,1),&			! 32 real(r8) X
				n2o(1,1),&			! 33 real(r8) X
				co,&				! 34 real(r8) X
				gbx)				! OUT
    
   print *, 'Populating input structure for COSP running on day columns only...'
   !! be very explicit about gbx structure sizes.
   !! flip vertical dimension from COSP convention to CAM convention
   !! variables initialzed as having 1:pcols dimension, but gbx only wants 1:Nday.

   gbx%longitude = lon_cosp_day(1:Nday)  			! lon (degrees_east)
   gbx%latitude = lat_cosp_day(1:Nday)    			! lat (degrees_north)
   gbx%p = ptop_day(1:Nday,1:pver) 				! p_in_full_levels (Pa),upper interface (flipping done above)
   gbx%ph = pmid_day(1:Nday,pver:1:-1)				! p_in_half_levels (Pa)
   gbx%zlev = ztop_day(1:Nday,1:pver)				! height_in_full_levels (m),upper interface,asl (flipping done above)
   gbx%zlev_half = zmid_day(1:Nday,pver:1:-1)			! height_in_half_levels (m),asl
   gbx%T = t_day(1:Nday,pver:1:-1)				! air_temperature (K)
   gbx%q = rh_day(1:Nday,pver:1:-1)				! relative humidity wrt water 
								! note: it is confusing that it is called "q" within cosp,
								! but it is rh according to Yuying
   gbx%sh = q_day(1:Nday,pver:1:-1)				! specific humidity (kg/kg), range [0,0.019ish]
   gbx%cca = concld_day(1:Nday,pver:1:-1)			! convective_cloud_amount (0-1)
   gbx%tca = cld_day(1:Nday,pver:1:-1)				! total_cloud_amount (0-1)
   gbx%psfc = ps_day(1:Nday)					! surface pressure (Pa)
   gbx%skt  = ts_day(1:Nday)					! skin temperature (K)
   gbx%land = landmask_day(1:Nday)				! landmask (0 or 1)
   gbx%mr_ozone  = o3_day(1:Nday,pver:1:-1)			! ozone mass mixing ratio (kg/kg)
   gbx%u_wind  = us_day(1:Nday)					! surface u_wind (m/s)
   gbx%v_wind  = vs_day(1:Nday)					! surface v_wind (m/s)
   gbx%sunlit  = 1
   gbx%mr_hydro(:,:,I_LSCLIQ) = cldliq_day(1:Nday,pver:1:-1)		! mr_lsliq, mixing_ratio_large_scale_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_LSCICE) = cldice_day(1:Nday,pver:1:-1)		! mr_lsice - mixing_ratio_large_scale_cloud_ice (kg/kg)
   gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq_day(1:Nday,pver:1:-1)	! mr_ccliq - mixing_ratio_convective_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice_day(1:Nday,pver:1:-1)	! mr_ccice - mixing_ratio_convective_cloud_ice (kg/kg)
   gbx%rain_ls = rain_ls_interp_day(1:Nday,pver:1:-1)		! midpoint gbm
   gbx%snow_ls = snow_ls_interp_day(1:Nday,pver:1:-1)		! midpoint gbm			
   gbx%grpl_ls = grpl_ls_interp_day(1:Nday,pver:1:-1)		! midpoint gbm
   gbx%rain_cv = rain_cv_interp_day(1:Nday,pver:1:-1)		! midpoint gbm
   gbx%snow_cv = snow_cv_interp_day(1:Nday,pver:1:-1)		! midpoint gbm
   gbx%Reff = reff_cosp_day(1:Nday,pver:1:-1,:)
   gbx%dtau_s   = dtau_s_day(1:Nday,pver:1:-1)
   gbx%dtau_c   = dtau_c_day(1:Nday,pver:1:-1)
   gbx%dem_s    = dem_s_day(1:Nday,pver:1:-1)
   gbx%dem_c    = dem_c_day(1:Nday,pver:1:-1)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! STARTUP RELATED TO COSP OUTPUT (see cosp.F90) for sunlit columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Note: COSP output variables are sgx (sub-grid outputs), sgradar (radar outputs), sglidar (lidar outputs), 
! isccp (isccp outputs), misr (misr simulator outputs), vgrid (vertical grid info), stradar 
! (summary statistics radar simulator), stlidar (summary statistics lidar simulator)

! Define new vertical grid
   print *, 'Defining new vertical grid...'
   call construct_cosp_vgrid(gbx,nlr,use_vgrid,csat_vgrid,vgrid)
        
! Allocate memory for output
   print *, 'Allocating memory for other types...'
   call construct_cosp_subgrid(Npoints, ncolumns, Nlevels, sgx)
   call construct_cosp_sgradar(cfg,Npoints,ncolumns,Nlevels,nhydro,sgradar)
   call construct_cosp_radarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,stradar)
   call construct_cosp_sglidar(cfg,Npoints,ncolumns,Nlevels,nhydro,PARASOL_NREFL,sglidar)
   call construct_cosp_lidarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,PARASOL_NREFL,stlidar)
   call construct_cosp_isccp(cfg,Npoints,ncolumns,Nlevels,isccp)
   call construct_cosp_misr(cfg,Npoints,misr)
   !! new for cosp v1.3
   call construct_cosp_modis(cfg,Npoints,modis)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  RUN COSP only on sunlit columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   print *, 'Calling the COSP simulator and running on sunlit columns...'
   call cosp(overlap,ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  TRANSLATE COSP OUTPUT INTO INDIVIDUAL VARIABLES FOR OUTPUT (see nc_write_cosp_1d in cosp_io.f90)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Populate CAM vars with COSP output (the reasoning for this is based on cosp_io.f90)
   !! steps 1) set as fill value, 2) fill in the columns calculated by COSP, 3) expand array
   !! only do this for the isccp-related and misr output
   !!! populate variables with fill values where COSP didn't run.
   ! Note: To expand to output CAM arrays, use ExpDayNite (function from ..atm/cam/src/physics/cam/cmparray_mod.F90, see also radsw.F90)
   ! e.g.,    call ExpDayNite(solin,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
   !   call ExpDayNite(pmxrgn,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pverp)

if (lisccp_sim) then
   ! (1d variables)
   cldtot_isccp(1:Nday) = isccp%totalcldarea
   call ExpDayNite(cldtot_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   meanptop_isccp(1:Nday) = isccp%meanptop
   call ExpDayNite(meanptop_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   meantau_isccp(1:Nday) = isccp%meantaucld
   call ExpDayNite(meantau_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   meancldalb_isccp(1:Nday) = isccp%meanalbedocld
   call ExpDayNite(meancldalb_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   meantb_isccp(1:Nday) = isccp%meantb
   call ExpDayNite(meantb_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   meantbclr_isccp(1:Nday) = isccp%meantbclr
   call ExpDayNite(meantbclr_isccp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

   ! (2d variables), loop over subcolumns to do expansion.
   do i=1,nscol_cosp
      tmp(1:Nday) = isccp%boxtau(1:Nday,i)	! CAM version of boxtauisccp (time,column,profile)
      call ExpDayNite(tmp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
      tau_isccp(1:ncol,i) = tmp(1:ncol)
      tmp(1:pcols)=fillvalue
   end do

   do i=1,nscol_cosp
      tmp(1:Nday) = isccp%boxptop(1:Nday,i)	! CAM version of boxtauisccp (time,column,profile)
      call ExpDayNite(tmp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
      cldptop_isccp(1:ncol,i) = tmp(1:ncol)
      tmp(1:pcols)=fillvalue
   end do

   ! (3d variables), loop over tau and nprs_cosp/htmisr_cosp to do expansion.
   !!  I think there might be a problem here...
   do it=1,ntau_cosp
   do ip=1,nprs_cosp
      tmp(1:Nday) = isccp%fq_isccp(1:Nday,it,ip)
      call ExpDayNite(tmp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)	
      clisccp2(1:ncol,it,ip) = tmp(1:ncol)	! clisccp2 (time,tau,plev,profile)
      tmp(1:pcols)=fillvalue
   end do
   end do

   ! Use high-dimensional output to populate CAM collapsed output variables
   do i=1,ncol
	! CAM clisccp2 (time,tau,plev,profile)
	do ip=1,nprs_cosp
	do it=1,ntau_cosp
	   ipt=(ip-1)*ntau_cosp+it
	   fisccp1(i,ipt) = clisccp2(i,it,ip)   		! fisccp1(pcols,ntau_cosp*nprs_cosp) 
	end do
	end do
   end do
end if

if (lmisr_sim) then
   do it=1,ntau_cosp
   do ihm=1,nhtmisr_cosp
      tmp(1:Nday) = misr%fq_MISR(1:Nday,it,ihm)
      call ExpDayNite(tmp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)	
      clMISR(1:ncol,it,ihm) = tmp(1:ncol)			! clMISR (time,tau,CTH_height_bin,profile)
      tmp(1:pcols)=fillvalue
   end do
   end do

   ! Use high-dimensional output to populate CAM collapsed output variables
   do i=1,ncol
	! CAM clMISR (time,tau,CTH_height_bin,profile)
	do ihm=1,nhtmisr_cosp
	do it=1,ntau_cosp
	   ihmt=(ihm-1)*ntau_cosp+it			
	   cld_misr(i,ihmt) = clMISR(i,it,ihm)   		! cld_misr(pcols,nhtmisr_cosp*ntau_cosp)
	end do
	end do
   end do
end if

if (lfrac_out) then
   frac_out = sgx%frac_out				! frac_out (time,height_mlev,column,profile)

   ! Use high-dimensional output to populate CAM collapsed output variables
   ! CAM frac_out (time,height_mlev,column,profile)  FIX
   do ihml=1,nhtml_cosp
   do isc=1,nscol_cosp
	ihsc=(ihml-1)*nscol_cosp+isc			
	scops_out(i,ihsc) = frac_out(i,isc,ihml)    	! scops_out(pcols,nht_cosp*nscol_cosp)
   end do
   end do
end if
 
if (lmodis_sim) then
      ! (1d variables from modis simulator)
      cltmodis(1:Nday) = modis%Cloud_Fraction_Total_Mean
      call ExpDayNite(cltmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      clwmodis(1:Nday) = modis%Cloud_Fraction_Water_Mean
      call ExpDayNite(clwmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      climodis(1:Nday) = modis%Cloud_Fraction_Ice_Mean
      call ExpDayNite(climodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tautmodis(1:Nday) = modis%Cloud_Fraction_High_Mean
      call ExpDayNite(tautmodis,Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauwmodis(1:Nday) = modis%Cloud_Fraction_Mid_Mean
      call ExpDayNite(tauwmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauimodis(1:Nday) = modis%Cloud_Fraction_Low_Mean
      call ExpDayNite(tauimodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tautmodis(1:Nday) = modis%Optical_Thickness_Total_Mean
      call ExpDayNite(tautmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauwmodis(1:Nday) = modis%Optical_Thickness_Water_Mean
      call ExpDayNite(tauwmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauimodis(1:Nday) = modis%Optical_Thickness_Ice_Mean
      call ExpDayNite(tauimodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tautlogmodis(1:Nday) = modis%Optical_Thickness_Total_LogMean
      call ExpDayNite(tautlogmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauwlogmodis(1:Nday) = modis%Optical_Thickness_Water_LogMean
      call ExpDayNite(tauwlogmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      tauilogmodis(1:Nday) = modis%Optical_Thickness_Ice_LogMean
      call ExpDayNite(tauilogmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      reffclwmodis(1:Nday) = modis%Cloud_Particle_Size_Water_Mean
      call ExpDayNite(reffclwmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      reffclimodis(1:Nday) = modis%Cloud_Particle_Size_Ice_Mean
      call ExpDayNite(reffclimodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      pctmodis(1:Nday) = modis%Cloud_Top_Pressure_Total_Mean
      call ExpDayNite(pctmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      lwpmodis(1:Nday) = modis%Liquid_Water_Path_Mean
      call ExpDayNite(lwpmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      iwpmodis(1:Nday) = modis%Ice_Water_Path_Mean
      call ExpDayNite(iwpmodis,	Nday, IdxDay, Nno, IdxNo, 1, pcols)

      ! (3d variable), need to loop 
      do it=1,ntau_cosp_modis
      do ip=1,nprs_cosp
      	 tmp(1:Nday) = modis%Optical_Thickness_vs_Cloud_Top_Pressure(1:Nday,it,ip)
      	 call ExpDayNite(tmp,	Nday, IdxDay, Nno, IdxNo, 1, pcols)	
      	 clmodis(1:ncol,it,ip) = tmp(1:ncol)
         tmp(1:pcols)=fillvalue
      end do
      end do

   ! Use high-dimensional output to populate CAM collapsed output variables
      do i=1,ncol
	   ! CAM clmodis
	   do ip=1,nprs_cosp
	   do it=1,ntau_cosp_modis
	      ipt=(ip-1)*ntau_cosp_modis+it
	      clmodis_cam(i,ipt) = clmodis(i,it,ip)
	   end do
	   end do
      end do
end if

   print *, 'Done mapping COSP output to local variables for running only ISCCP/MISR/MODIS on sunlit columns ...'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DEALLOCATE MEMORY for running only daylit columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call free_cosp_gridbox(gbx)
   call free_cosp_vgrid(vgrid)
   call free_cosp_subgrid(sgx)
   call free_cosp_sgradar(sgradar)
   call free_cosp_radarstats(stradar)
   call free_cosp_sglidar(sglidar)
   call free_cosp_lidarstats(stlidar)
   call free_cosp_isccp(isccp)
   call free_cosp_misr(misr) 
   call free_cosp_modis(modis)

end if !! if there are columns to run

end if  !!! END RUNNING COSP ONLY ON SUNLIT COLUMNS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  RUN COSP FOR ONLY RADAR AND/OR LIDAR SIMULATORS (COSP SETUP/CALL #3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ((.not.cosp_runall) .and. (lradar_sim.or.llidar_sim)) then !! RUNNING RADAR/LIDAR  

! Only run the lidar and radar simulators and their outputs
if (lradar_sim) then
   cfg%Lradar_sim= .true.
   cfg%Lcfaddbze94= .true.
   cfg%Ldbze94= .true.
   cfg%Lclcalipso2= .true.
   cfg%Lcltlidarradar= .true.
end if
if (llidar_sim) then
   cfg%Llidar_sim= .true.
   cfg%Latb532= .true.
   cfg%LcfadLidarsr532= .true.
   cfg%Lclcalipso= .true.
   cfg%Lclhcalipso= .true.
   cfg%Lclisccp= .true.
   cfg%Lcllcalipso= .true.
   cfg%Lclmcalipso= .true.
   cfg%Lcltcalipso= .true.
   cfg%LparasolRefl= .true.
   cfg%LlidarBetaMol532= .true.
end if

! Run cosp_stats.f90 for radar and lidar output
cfg%Lstats = .true.

! Don't run the iscc/misr/modis or misr simulators by changing COSP input structure cfg
cfg%Lisccp_sim= .false.
cfg%Lmisr_sim= .false.
cfg%Lmodis_sim= .false.

! Also set isccp,misr,modis output flags to false in COSP input structure cfg
cfg%Lalbisccp= .false.
cfg%Latb532= .false.
cfg%Lboxptopisccp= .false.
cfg%Lboxtauisccp= .false.
cfg%Lpctisccp= .false.
cfg%Ltauisccp= .false.
cfg%Lcltisccp= .false.
cfg%Lmeantbisccp= .false.
cfg%Lmeantbclrisccp= .false.
cfg%LclMISR= .false.
cfg%Lcltmodis= .false.
cfg%Lclwmodis= .false.
cfg%Lclimodis= .false.
cfg%Lclhmodis= .false.
cfg%Lclmmodis= .false.
cfg%Lcllmodis= .false.
cfg%Ltautmodis= .false.
cfg%Ltauwmodis= .false.
cfg%Ltauimodis= .false.
cfg%Ltautlogmodis= .false.
cfg%Ltauwlogmodis= .false.
cfg%Ltauilogmodis= .false.
cfg%Lreffclwmodis= .false.
cfg%Lreffclimodis= .false.
cfg%Lpctmodis= .false.
cfg%Llwpmodis= .false.
cfg%Liwpmodis= .false.
cfg%Lclmodis= .false.

! Subsetting using an orbital curtain for radar and lidar simulators
! Note: This is a very basic implementation and should be checked with other groups.

if (.not.cosp_sample_atrain) then
   Natrain = 0
   Npoints = ncol !!! running on all columns
end if

if (cosp_sample_atrain) then

! Subsetting using the orbit data file.
! Subsetting CAM columns that will have simulator run on them using time and geographic constraint
! 1) time - use fractional calendar day (1-365), eliminates the need for time of day.
! for day - CAM default is a NO_LEAP calendar
!    calendar = get_calendar()  !! 'NO_LEAP'
! always use first 365 days of atrain orbit for 2008, even though 2008 has 366 days.
! for time of day, use fractional value of calendar day.
! note that you can use this to find exactly where the Atrain orbit has been during the CAM timestep

   !! OK, apply time constraint by using these values to find a beginning index to loop over
   !!! FIND WHERE THE ATRAIN HAS BEEN DURING THIS CAM TIMESTEP
   !!  This will enable you to loop over idxas,idxae instead of looping over 1, norbitdata
   !!! first time look through whole array using a plane old do loop, then
   !!! keep track of last time index to look for new time in the next time step
   !!! use information on time sequence in cloudsat orbit - map calendar day to index in big list
   !!! idxas will be private data stored at the module data level,  index of the last time match, it will be overwritten every time you call.
   !!idxas = beginning of cam timestep

! 2) location - using lat/lon degree distance.
! is it better to use a distance approach?  because of convergence at the poles, this will lead to un-uniform sampling
! lat/lon only works for a rectangular grid.  what is the lat/lon threshold?  it will be resolution dependent!
!!! programming help notes
!! look in models/atm/cam/src/utils/time_manager.F90 for CAM time and date functions
!! get_curr_time used in models/atm/cam/src/control/cam_history.F90
!! for location - use lat_cosp, lon_cosp in the degrees_north and degrees_east already

!! initialize variable to decide if CAM column on A-train orbit
   atrain_use(1:ncol) = 0     !! array to indicate if column is in atrain orbit (1=yes,0=no)

!! TIME CONSTRAINT

   !! find idxas using calendar day from cam (calday) and from atrain (atrain_calday)
   cam_calday = get_curr_calday()    !returns fractional calendar day at end of current timestep
   atrain_calday = atrainday + (atrainhr*60*60.0 + atrainmin*60.0 + atrainsec)/86400.0

   !!print *, cam_calday
   !!print *, atrain_calday(1:200)
 
   !! if first time through or after failed time search, start at the beginning of the file
   if (idxas .eq. 0) then
     idxas = 1
   end if

   !! loop to find idxas
   caldaydiff = abs(cam_calday - atrain_calday(idxas))
   caldaymindiff = caldaydiff
   caldaymatch = 0

   do while (caldaydiff .gt. atrain_daythresh)
      idxas = idxas + 1
      caldaydiff = abs(cam_calday - atrain_calday(idxas))
      if (caldaydiff .lt. atrain_daythresh) then
	 !print *, 'TIME MATCH FOUND ...'
	 !print *, cam_calday
	 !print *, atrain_calday(idxas)
	 !print *, caldaydiff
         caldaymatch = 1
      end if
      if (caldaydiff .lt. caldaymindiff) then
          caldaymindiff = caldaydiff
      end if

      if (idxas .ge. norbitdata) then
	print *, 'No time match found for A-train orbit data, going back to beginning'
	idxas = 1
	caldaydiff = abs(cam_calday - atrain_calday(idxas))
	do while (caldaydiff .gt. atrain_daythresh)
	   idxas = idxas + 1
	   caldaydiff = abs(cam_calday - atrain_calday(idxas))
	   if (caldaydiff .lt. atrain_daythresh) then
	      print *, 'TIME MATCH FOUND AFTER GOING BACK TO BEGINNING...'
	      print *, cam_calday
	      print *, atrain_calday(idxas)
	      print *, caldaydiff
	      caldaymatch = 1
	   end if
	   if (caldaydiff .lt. caldaymindiff) then
	       caldaymindiff = caldaydiff
	   end if
	   if (idxas .ge. norbitdata) then
	     print *, 'No time match found for A-train orbit data after search from beginning ERROR!'
	     print *, 'Need to increase atrain_daythresh for overpass. Min difference in days is:'
             print *, caldaymindiff
	     !! go back to the beginning for next time search
	     idxas = 0
	     !! exit the inner while
	     exit
 	   end if
        end do
      end if
      !! no match found, exit the while
      if (idxas .eq. 0) then
	exit
      end if
   end do

!! GEOGRAPHIC CONSTRAINT
   !! starting with a lat/lon degree criteria
   !! already have lat_cosp and lon_cosp in the right units.

if (caldaymatch .eq. 1) then
   
   !! initialize 
   tmpvar(1:ntmp) = 0.0
   tmpdiff(1:ntmp) = 0.0
   idxlat(1:ntmp) = 0.0
   idxlon(1:ntmp) = 0.0
   idxsum(1:ntmp) = 0.0

   !! test by forcing a location match
   !!print *, atrainlat(idxas:idxas+10)
   !!print *, atrainlon(idxas:idxas+10)
   !!lat_cosp(1)=atrainlat(idxas)
   !!lon_cosp(1)=atrainlon(idxas)

   do i = 1, ncol
       !! right latitude? atrainlat vs. lat_cosp(i)
       !! 1) create a vector of CAM latitude
       tmpvar(1:ntmp) = lat_cosp(i)
       !! 2) difference the CAM latitude with the Atrain orbit latitudes only over the right time!
       !! lehey compains "Two entities must be in shape conformance"
       tmpdiff = atrainlat(idxas:idxas+ntmp-1) - tmpvar
       !! 3) populate array to indicate if a matching latitude was found
       idxlat(1:ntmp) = 0
       do k = 1,ntmp
	  if (abs(tmpdiff(k)) .lt. atrain_latthresh) then
	     !!print *, 'LAT MATCH FOUND ...'
	     idxlat(k) = 1
	  end if
       end do

       !! right longitude? atrainlon vs. lon_cosp(i) -- same thing as for latitude
       tmpvar(1:ntmp) = lon_cosp(i)
       tmpdiff = atrainlon(idxas:idxas+ntmp-1) - tmpvar
       idxlon(1:ntmp) = 0		         !! array to indicate if a matching longitude was found
       do k = 1,ntmp
	  if (abs(tmpdiff(k)) .lt. atrain_lonthresh) then
	     !!print *, 'LON MATCH FOUND ...'
	     idxlon(k) = 1
	  end if
       end do

       !! sum the two query arrays and populate array to indicate if lat/lon criteria was met
       idxsum = idxlat + idxlon
       !!print *, maxval(idxsum)
       if (maxval(idxsum) .eq. 2) then
	  !!print *, 'LOCATION AND CALDAY MATCH FOUND ...'
	  !!atrain_use(i) = 1
       end if
   end do

end if !! if calday match

!!! assign indices for compression (note: code for day/nite from radsw.F90)
    Natrain = 0
    Nno = 0
    do i = 1, ncol
       !! figure out if part of atrain orbit
	if (atrain_use(i) .eq. 1) then
          Natrain = Natrain + 1
          IdxAtrain(Natrain) = i
       else
          Nno = Nno + 1
          IdxNo(Nno) = i
       end if
    end do

    Npoints = Natrain

if (Natrain .gt. 0) then  !! only run if you need to, otherwise default values used.

!!! rearrange input arrays - compress arrays so that COSP only runs on atrain columns
!e.g., a la
!   call CmpDayNite(E_pmid, pmid,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pver)
! (function from ..atm/cam/src/physics/cam/cmparray_mod.F90, see also radsw.F90)
!! CmpDayNite_2d_R_Copy(InArray, OutArray

call CmpDayNite(lon_cosp, lon_cosp_atrain, 	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(lat_cosp, lat_cosp_atrain, 	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(ptop, ptop_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%pmid, pmid_atrain,         Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(ztop, ztop_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%zm, zmid_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%t,t_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rh, rh_atrain,			Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(q, q_atrain,			Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(concld, concld_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(cld, cld_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%ps, ps_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(cam_in%ts ,ts_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(landmask,landmask_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(o3, o3_atrain,			Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%u(1:ncol,pver),us_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(state%v(1:ncol,pver),vs_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
call CmpDayNite(state%q(:,:,ixcldliq),cldliq_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(state%q(:,:,ixcldice),cldice_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(mr_ccliq, mr_ccliq_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(mr_ccice, mr_ccice_atrain,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rain_ls_interp, rain_ls_interp_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(snow_ls_interp, snow_ls_interp_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(grpl_ls_interp, grpl_ls_interp_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(rain_cv_interp,rain_cv_interp_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)			
call CmpDayNite(snow_cv_interp,snow_cv_interp_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dtau_s,dtau_s_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dtau_c,dtau_c_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dem_s,dem_s_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
call CmpDayNite(dem_c,dem_c_atrain,		Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)

!! need to loop over nhydro for reff_cosp_atrain
do i=1,nhydro
   call CmpDayNite(reff_cosp(:,:,i),tmp1,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols, 1, pver)
   reff_cosp_atrain(:,:,i)=tmp1
   tmp1(1:pcols,1:pver)=fillvalue
end do

else

   print *, 'No columns met Atrain sampling requirement...'

end if ! if there are atrain columns to run

else  !! not atrain sampling, 

!! running lidar or radar simulator on all columns
!! populate gbx inputs with all CAM columns
  lon_cosp_atrain=lon_cosp
  lat_cosp_atrain=lat_cosp
  ptop_atrain=ptop
  pmid_atrain=state%pmid
  ztop_atrain=ztop
  zmid_atrain=state%zm
  t_atrain=state%t
  rh_atrain=rh
  q_atrain=q
  concld_atrain=concld
  cld_atrain=cld
  ps_atrain=state%ps
  ts_atrain=cam_in%ts
  landmask_atrain=landmask
  o3_atrain=o3
  us_atrain(1:ncol)=state%u(1:ncol,pver)
  vs_atrain(1:ncol)=state%v(1:ncol,pver)
  cldliq_atrain=state%q(:,:,ixcldliq)
  cldice_atrain=state%q(:,:,ixcldice)
  mr_ccliq_atrain=mr_ccliq
  mr_ccice_atrain=mr_ccice
  rain_ls_interp_atrain=rain_ls_interp
  snow_ls_interp_atrain=snow_ls_interp
  grpl_ls_interp_atrain=grpl_ls_interp
  rain_cv_interp_atrain=rain_cv_interp			
  snow_cv_interp_atrain=snow_cv_interp
  dtau_s_atrain=dtau_s
  dtau_c_atrain=dtau_c
  dem_s_atrain=dem_s
  dem_c_atrain=dem_c
  reff_cosp_atrain=reff_cosp

end if ! if atrain sampling

!! if not atrain sampling or there are points to sample within atrain
if ((.not.cosp_sample_atrain) .or. ((cosp_sample_atrain) .and. Natrain .gt. 0)) then

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  POPULATE COSP INPUT VARIABLE ("gbx") FROM CAM VARIABLES for atrain columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   print *, 'Allocating memory for gridbox type for isccp sub-columns only...'
   call construct_cosp_gridbox(time, &			        ! 1 double precision = real(r8) X
				time_bnds, &			! 1 double precision = real(r8)
				radar_freq, &			! 2 real(r8) X
				surface_radar, &		! 3 integer X
				use_mie_tables, &		! 4 integer X
				use_gas_abs, &			! 5 integer X
				do_ray, &			! 6 integer X 
				melt_lay, &			! 7 integer X 
				k2, &				! 8 real(r8) X
                                Npoints, &			! 9 integer X
				Nlevels, &			! 10 integer X
				ncolumns,&			! 11 integer X
				nhydro,&			! 12 integer X
				Nprmts_max_hydro,&      	! 13 integer X
				Naero,&				! 14 integer X
				Nprmts_max_aero,&		! 15 integer X
				Npoints_it, &			! 16 integer X
                                lidar_ice_type,&		! 17 integer X
				isccp_topheight,&		! 18 integer X
				isccp_topheight_direction,& 	! 19 integer X
				overlap,&			! 20 integer X
				emsfc_lw, &			! 21 real X
                                use_precipitation_fluxes,& 	! 22 logical X
				use_reff, &			! 23 logical X
                                Platform, &			! 24 integer X
				Satellite, &			! 25 integer X
				Instrument, &			! 26 integer X
				Nchannels, &			! 27 integer X
				ZenAng, &			! 28 real(r8) X
                                Channels(1:Nchannels),&		! 29 integer X
				Surfem(1:Nchannels),&		! 30 real(r8) X
				co2(1,1),&			! 31 real(r8) X
				ch4(1,1),&			! 32 real(r8) X
				n2o(1,1),&			! 33 real(r8) X
				co,&				! 34 real(r8) X
				gbx)				! OUT
    
   print *, 'Populating input structure for COSP for running lidar/radar simulators only...'
   !! be very explicit about gbx structure sizes.
   !! flip vertical dimension from COSP convention to CAM convention
   !! Npoints can be ncol or Natrain - so here I just specify dimensions as Npoints
   gbx%longitude = lon_cosp_atrain(1:Npoints)  			! lon (degrees_east)
   gbx%latitude = lat_cosp_atrain(1:Npoints)   			! lat (degrees_north)
   gbx%p = ptop_atrain(1:Npoints,1:pver) 			! p_in_full_levels (Pa),upper interface (flipping done above)
   gbx%ph = pmid_atrain(1:Npoints,pver:1:-1)			! p_in_half_levels (Pa)
   gbx%zlev = ztop_atrain(1:Npoints,1:pver)			! height_in_full_levels (m),upper interface,asl (flipping done above)
   gbx%zlev_half = zmid_atrain(1:Npoints,pver:1:-1)		! height_in_half_levels (m),asl
   gbx%T = t_atrain(1:Npoints,pver:1:-1)			! air_temperature (K)
   gbx%q = rh_atrain(1:Npoints,pver:1:-1)			! relative humidity wrt water 
								! note: it is confusing that it is called "q" within cosp,
								! but it is rh according to Yuying
   gbx%sh = q_atrain(1:Npoints,pver:1:-1)			! specific humidity (kg/kg), range [0,0.019ish]
   gbx%cca = concld_atrain(1:Npoints,pver:1:-1)			! convective_cloud_amount (0-1)
   gbx%tca = cld_atrain(1:Npoints,pver:1:-1)			! total_cloud_amount (0-1)
   gbx%psfc = ps_atrain(1:Npoints)				! surface pressure (Pa)
   gbx%skt  = ts_atrain(1:Npoints)				! skin temperature (K)
   gbx%land = landmask_atrain(1:Npoints)			! landmask (0 or 1)
   gbx%mr_ozone  = o3_atrain(1:Npoints,pver:1:-1)		! ozone mass mixing ratio (kg/kg)
   gbx%u_wind  = us_atrain(1:Npoints)				! surface u_wind (m/s)
   gbx%v_wind  = vs_atrain(1:Npoints)				! surface v_wind (m/s)
   gbx%sunlit  = 1
   gbx%mr_hydro(:,:,I_LSCLIQ) = cldliq_atrain(1:Npoints,pver:1:-1)	! mr_lsliq, mixing_ratio_large_scale_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_LSCICE) = cldice_atrain(1:Npoints,pver:1:-1)	! mr_lsice - mixing_ratio_large_scale_cloud_ice (kg/kg)
   gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq_atrain(1:Npoints,pver:1:-1)	! mr_ccliq - mixing_ratio_convective_cloud_liquid (kg/kg)
   gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice_atrain(1:Npoints,pver:1:-1)	! mr_ccice - mixing_ratio_convective_cloud_ice (kg/kg)
   gbx%rain_ls = rain_ls_interp_atrain(1:Npoints,pver:1:-1)		! midpoint gbm
   gbx%snow_ls = snow_ls_interp_atrain(1:Npoints,pver:1:-1)		! midpoint gbm			
   gbx%grpl_ls = grpl_ls_interp_atrain(1:Npoints,pver:1:-1)		! midpoint gbm
   gbx%rain_cv = rain_cv_interp_atrain(1:Npoints,pver:1:-1)		! midpoint gbm
   gbx%snow_cv = snow_cv_interp_atrain(1:Npoints,pver:1:-1)		! midpoint gbm
   gbx%Reff = reff_cosp_atrain(1:Npoints,pver:1:-1,:)
   gbx%dtau_s   = dtau_s_atrain(1:Npoints,pver:1:-1)
   gbx%dtau_c   = dtau_c_atrain(1:Npoints,pver:1:-1)
   gbx%dem_s    = dem_s_atrain(1:Npoints,pver:1:-1)
   gbx%dem_c    = dem_c_atrain(1:Npoints,pver:1:-1)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! STARTUP RELATED TO COSP OUTPUT (see cosp.F90) for atrain columns
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Note: COSP output variables are sgx (sub-grid outputs), sgradar (radar outputs), sglidar (lidar outputs), 
! isccp (isccp outputs), misr (misr simulator outputs), vgrid (vertical grid info), stradar 
! (summary statistics radar simulator), stlidar (summary statistics lidar simulator)

! Define new vertical grid
   print *, 'Defining new vertical grid...'
   call construct_cosp_vgrid(gbx,nlr,use_vgrid,csat_vgrid,vgrid)
        
! Allocate memory for output
   print *, 'Allocating memory for other types...'
   call construct_cosp_subgrid(Npoints, ncolumns, Nlevels, sgx)
   call construct_cosp_sgradar(cfg,Npoints,ncolumns,Nlevels,nhydro,sgradar)
   call construct_cosp_radarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,stradar)
   call construct_cosp_sglidar(cfg,Npoints,ncolumns,Nlevels,nhydro,PARASOL_NREFL,sglidar)
   call construct_cosp_lidarstats(cfg,Npoints,ncolumns,vgrid%Nlvgrid,nhydro,PARASOL_NREFL,stlidar)
   call construct_cosp_isccp(cfg,Npoints,ncolumns,Nlevels,isccp)
   call construct_cosp_misr(cfg,Npoints,misr)
   !! new for cosp v1.3
   call construct_cosp_modis(cfg,Npoints,modis)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  RUN COSP - RADAR/LIDAR only
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   print *, 'Calling the COSP simulator and running on all/atrain columns...'
   call cosp(overlap,ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  TRANSLATE COSP OUTPUT INTO INDIVIDUAL VARIABLES FOR OUTPUT (see nc_write_cosp_1d in cosp_io.f90)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Map from COSP output to CAM output, only for simulators that were run!
   ! Populate CAM vars with COSP output (the reasoning for this is based on cosp_io.f90)
   !! steps 1) set as fill value, 2) fill in the columns calculated by COSP, 3) expand array
   ! Note: To expand to output CAM arrays, use ExpDayNite (function from ..atm/cam/src/physics/cam/cmparray_mod.F90, see also radsw.F90)
   ! e.g.,    call ExpDayNite(solin,	Nday, IdxDay, Nno, IdxNo, 1, pcols)
   !   call ExpDayNite(pmxrgn,	Nday, IdxDay, Nno, IdxNo, 1, pcols, 1, pverp)

if (lradar_sim) then

   if (cosp_sample_atrain) then
      ! (1d variables from radar simulator)
      cldtot_calcs(1:Natrain)=stradar%radar_lidar_tcc	! CAM version of cltlidarradar (time,profile)
      call ExpDayNite(cldtot_calcs,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)

      ! (2d variables from radar simulator), loop over nht_cosp to do expansion.
      do i=1,nht_cosp
	 tmp(1:Natrain) = stradar%lidar_only_freq_cloud(1:Natrain,i)	! CAM version of clcalipso2 (time,height,profile)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
	 cld_cal_notcs(1:ncol,i)=tmp(1:ncol)
         tmp(1:pcols)=fillvalue
      end do

      ! (3d variables from radar simulator), loop over tau and plev to do expansion.
      do isc=1,nscol_cosp
      do ihml=1,nhtml_cosp
	 tmp(1:Natrain) = sgradar%Ze_tot(1:Natrain,isc,ihml)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)	
	 dbze94(1:ncol,isc,ihml)=tmp(1:ncol)		! dbze94 (time,height_mlev,column,profile)	
         tmp(1:pcols)=fillvalue	
      end do
      end do
      do id=1,ndbze_cosp
      do ih=1,nht_cosp
	 tmp(1:Natrain) = stradar%cfad_ze(1:Natrain,id,ih)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)	
	 cfad_dbze94(1:ncol,id,ih)=tmp(1:ncol)		! cfad_dbze94 (time,height,dbze,profile)
         tmp(1:pcols)=fillvalue			
      end do
      end do

   else
      ! (1d variables)
      cldtot_calcs(1:ncol) = stradar%radar_lidar_tcc		! CAM version of cltlidarradar (time,profile)
      ! (2d variables)
      cld_cal_notcs(1:ncol,1:nht_cosp) = stradar%lidar_only_freq_cloud	! CAM version of clcalipso2 (time,height,profile)
      ! (3d variables)
      dbze94(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sgradar%Ze_tot	 ! dbze94 (time,height_mlev,column,profile)
      cfad_dbze94(1:ncol,1:ndbze_cosp,1:nht_cosp) = stradar%cfad_ze	! cfad_dbze94 (time,height,dbze,profile)
   end if

   ! Use high-dimensional output to populate CAM collapsed output variables
   do i=1,ncol
	! CAM dbze94 (time,height_mlev,column,profile)
	do ihml=1,nhtml_cosp
	do isc=1,nscol_cosp
	   ihsc=(ihml-1)*nscol_cosp+isc			
	   dbze_cs(i,ihsc) = dbze94(i,isc,ihml)   		! dbze_cs(pcols,pver*nscol_cosp) 
	end do
	end do
	! CAM cfad_dbze94 (time,height,dbze,profile) 
	do ih=1,nht_cosp
	do id=1,ndbze_cosp
	   ihd=(ih-1)*ndbze_cosp+id			
	   cfad_dbze94_cs(i,ihd) = cfad_dbze94(i,id,ih)  	! cfad_dbze94_cs(pcols,nht_cosp*ndbze_cosp)
	end do
	end do
   end do

end if

if (llidar_sim) then

   if (cosp_sample_atrain) then

      ! (1d variables from lidar simulator)
      cldlow_cal(1:Natrain) = stlidar%cldlayer(:,1)	! CAM version of cllcalipso (time,profile)
      call ExpDayNite(cldlow_cal,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
      cldmed_cal(1:Natrain) = stlidar%cldlayer(:,2)	! CAM version of clmcalipso (time,profile)
      call ExpDayNite(cldmed_cal,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
      cldhgh_cal(1:Natrain) = stlidar%cldlayer(:,3)	! CAM version of clhcalipso (time,profile)
      call ExpDayNite(cldhgh_cal,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
      cldtot_cal(1:Natrain) = stlidar%cldlayer(:,4)	! CAM version of cltcalipso (time,profile)
      call ExpDayNite(cldtot_cal,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)

      ! (2d variables from lidar simulator), loop over 2nd dimension to do expansion.
      do i=1,nht_cosp
	 tmp(1:Natrain) = stlidar%lidarcld(1:Natrain,i)	! CAM version of clcalipso (time,height,profile)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
	 cld_cal(1:ncol,i) = tmp(1:ncol)
         tmp(1:pcols)=fillvalue
      end do

      do i=1,nhtml_cosp
	 tmp(1:Natrain) = sglidar%beta_mol(1:Natrain,i)	! CAM version of beta_mol532 (time,height_mlev,profile) 
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
	 mol532_cal(1:ncol,i) = tmp(1:ncol)
         tmp(1:pcols)=fillvalue
      end do

      do i=1,nsza_cosp
	 tmp(1:Natrain) = stlidar%parasolrefl(1:Natrain,i)	! CAM version of parasolrefl (time,sza,profile)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)
	 refl_parasol(1:ncol,i) = tmp(1:ncol)
         tmp(1:pcols)=fillvalue
      end do

      ! (3d variables from lidar simulator), loop over 2nd and 3rd dimensions to do expansion.
      do isc=1,nscol_cosp
      do ihml=1,nhtml_cosp
	 tmp(1:Natrain) = sglidar%beta_tot(1:Natrain,isc,ihml)
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)	
	 atb532(1:ncol,isc,ihml) = tmp(1:ncol)		! atb532 (time,height_mlev,column,profile)
         tmp(1:pcols)=fillvalue	
      end do
      end do

      do is=1,nsr_cosp
      do ih=1,nht_cosp
	 tmp(1:Natrain) = stlidar%cfad_sr(1:Natrain,is,ih)	
	 call ExpDayNite(tmp,	Natrain, IdxAtrain, Nno, IdxNo, 1, pcols)	
	 cfad_lidarsr532(1:ncol,is,ih) = tmp(1:ncol)	! cfad_lidarsr532 (time,height,scat_ratio,profile)
         tmp(1:pcols)=fillvalue		
      end do
      end do

   else
      ! (1d variables)
      cldlow_cal(1:ncol) = stlidar%cldlayer(:,1)		! CAM version of cllcalipso (time,profile)	
      cldmed_cal(1:ncol) = stlidar%cldlayer(:,2)		! CAM version of clmcalipso (time,profile)
      cldhgh_cal(1:ncol) = stlidar%cldlayer(:,3)		! CAM version of clhcalipso (time,profile)
      cldtot_cal(1:ncol) = stlidar%cldlayer(:,4)		! CAM version of cltcalipso (time,profile)
      ! (2d variables)
      cld_cal(1:ncol,1:nht_cosp) = stlidar%lidarcld		! CAM version of clcalipso (time,height,profile)
      mol532_cal(1:ncol,1:nhtml_cosp) = sglidar%beta_mol	! CAM version of beta_mol532 (time,height_mlev,profile)
      refl_parasol(1:ncol,1:nsza_cosp) = stlidar%parasolrefl	! CAM version of parasolrefl (time,sza,profile)
      ! (3d variables)
      atb532(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sglidar%beta_tot	! atb532 (time,height_mlev,column,profile)
      cfad_lidarsr532(1:ncol,1:nsr_cosp,1:nht_cosp)= stlidar%cfad_sr	! cfad_lidarsr532 (time,height,scat_ratio,profile)
   end if

   ! Use high-dimensional output to populate CAM collapsed output variables
   do i=1,ncol
	! CAM atb532 (time,height_mlev,column,profile)  FIX
	do ihml=1,nhtml_cosp
	do isc=1,nscol_cosp
	   ihsc=(ihml-1)*nscol_cosp+isc			
	   atb532_cal(i,ihsc) = atb532(i,isc,ihml)  	 	! atb532_cal(pcols,nht_cosp*nscol_cosp)
	end do
	end do
	! CAM cfad_lidarsr532 (time,height,scat_ratio,profile)
	do ih=1,nht_cosp
	do is=1,nsr_cosp
	   ihs=(ih-1)*nsr_cosp+is			
	   cfad_sr532_cal(i,ihs) = cfad_lidarsr532(i,is,ih)  ! cfad_sr532_cal(pcols,nht_cosp*nsr_cosp)
	end do
	end do
   end do

end if

if ((lfrac_out) .and. (.not.lradar_sim) .and. (.not.llidar_sim)) then
   frac_out(1:ncol,1:nscol_cosp,1:nhtml_cosp) = sgx%frac_out		! frac_out (time,height_mlev,column,profile)

   ! Use high-dimensional output to populate CAM collapsed output variables
   ! CAM frac_out (time,height_mlev,column,profile)  FIX
   do ihml=1,nhtml_cosp
   do isc=1,nscol_cosp
	ihsc=(ihml-1)*nscol_cosp+isc			
	scops_out(i,ihsc) = frac_out(i,isc,ihml)    	! scops_out(pcols,nht_cosp*nscol_cosp)
   end do
   end do
end if

   print *, 'Done mapping COSP output to local variables for running only radar/lidar ...'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DEALLOCATE MEMORY for running only radar/lidar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call free_cosp_gridbox(gbx)
   call free_cosp_vgrid(vgrid)
   call free_cosp_subgrid(sgx)
   call free_cosp_sgradar(sgradar)
   call free_cosp_radarstats(stradar)
   call free_cosp_sglidar(sglidar)
   call free_cosp_lidarstats(stlidar)
   call free_cosp_isccp(isccp)
   call free_cosp_misr(misr) 
   call free_cosp_modis(modis)

end if !! if cosp_sample_atrain = false or there are atrain columns to run

end if  !!! END RUNNING COSP ONLY RADAR/LIDAR

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OUTFLD CALLS
! Note: I replace COSP missing values with CAM missing values (fillvalue=1e+36)
! otherwise the averaging done by cam_history.F90 will be screwed up.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!! ISCCP OUTPUTS
   if (lisccp_sim) then
      where (fisccp1(:ncol,:ntau_cosp*nprs_cosp) .eq. -1.00000E+30_r8)
          fisccp1(:ncol,:ntau_cosp*nprs_cosp) = fillvalue
      end where   
      call outfld('FISCCP1_COSP',fisccp1, pcols,lchnk)

      where (cldtot_isccp(:ncol) .eq. -1.00000E+30_r8)
           cldtot_isccp(:ncol) = fillvalue
      end where  
      call outfld('CLDTOT_ISCCP',cldtot_isccp    ,pcols,lchnk)

      where (meancldalb_isccp(:ncol) .eq. -1.00000E+30_r8)
           meancldalb_isccp(:ncol) = fillvalue
      end where
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction
      where ((meancldalb_isccp(:ncol) .ne. fillvalue)  .and. (cldtot_isccp(:ncol) .ne. fillvalue))
      	   meancldalb_isccp(:ncol) = meancldalb_isccp(:ncol)*cldtot_isccp(:ncol)
      end where
      where (cldtot_isccp(:ncol) .eq. fillvalue)
      	   meancldalb_isccp(:ncol) = fillvalue
      end where
      call outfld('MEANCLDALB_ISCCP',meancldalb_isccp,pcols,lchnk)

      where (meanptop_isccp(:ncol) .eq. -1.00000E+30_r8)
           meanptop_isccp(:ncol) = fillvalue
      end where  
      !!! Per CFMIP request - weight by ISCCP Total Cloud Fraction
      where ((meanptop_isccp(:ncol) .ne. fillvalue)  .and. (cldtot_isccp(:ncol) .ne. fillvalue))
      	   meanptop_isccp(:ncol) = meanptop_isccp(:ncol)*cldtot_isccp(:ncol)
      end where
      where (cldtot_isccp(:ncol) .eq. fillvalue)
      	   meanptop_isccp(:ncol) = fillvalue
      end where
      call outfld('MEANPTOP_ISCCP',meanptop_isccp    ,pcols,lchnk)

      where (tau_isccp(:ncol,:nscol_cosp) .eq. -1.00000E+30_r8)
            tau_isccp(:ncol,:nscol_cosp) = fillvalue
      end where  
      call outfld('TAU_ISCCP',tau_isccp    ,pcols,lchnk) !! fails check_accum if 'A'

      where (cldptop_isccp(:ncol,:nscol_cosp) .eq. -1.00000E+30_r8)
            cldptop_isccp(:ncol,:nscol_cosp) = fillvalue
      end where  
      call outfld('CLDPTOP_ISCCP',cldptop_isccp    ,pcols,lchnk) !! fails check_accum if 'A'

      where (meantau_isccp(:ncol) .eq. -1.00000E+30_r8)
            meantau_isccp(:ncol) = fillvalue
      end where  
      call outfld('MEANTAU_ISCCP',meantau_isccp  ,pcols,lchnk)

      where (meantb_isccp(:ncol) .eq. -1.00000E+30_r8)
            meantb_isccp(:ncol) = fillvalue
      end where  
      call outfld('MEANTB_ISCCP',     meantb_isccp   ,pcols,lchnk)

      where (meantbclr_isccp(:ncol) .eq. -1.00000E+30_r8)
            meantbclr_isccp(:ncol) = fillvalue
      end where  
      call outfld('MEANTBCLR_ISCCP',  meantbclr_isccp ,pcols,lchnk)

   end if

!!! LIDAR SIMULATOR OUTPUTS
   if (llidar_sim) then
      where (cldlow_cal(:ncol) .eq. -1.00000E+30_r8)
            cldlow_cal(:ncol) = fillvalue
      end where  
      call outfld('CLDLOW_CAL',cldlow_cal    ,pcols,lchnk)

      where (cldmed_cal(:ncol) .eq. -1.00000E+30_r8)
            cldmed_cal(:ncol) = fillvalue
      end where  
      call outfld('CLDMED_CAL',cldmed_cal    ,pcols,lchnk)

      where (cldhgh_cal(:ncol) .eq. -1.00000E+30_r8)
            cldhgh_cal(:ncol) = fillvalue
      end where  
      call outfld('CLDHGH_CAL',cldhgh_cal    ,pcols,lchnk)

      where (cldtot_cal(:ncol) .eq. -1.00000E+30_r8)
            cldtot_cal(:ncol) = fillvalue
      end where  
      call outfld('CLDTOT_CAL',cldtot_cal    ,pcols,lchnk)

      where (cld_cal(:ncol,:nht_cosp) .eq. -1.00000E+30_r8)
	    !! setting missing values to 0 (clear air).  
	    !! I'm not sure why COSP produces a mix of fillvalue and realvalue in the nht_cosp dimension.
            cld_cal(:ncol,:nht_cosp) = 0.0
      end where  
      call outfld('CLD_CAL',cld_cal    ,pcols,lchnk)  !! fails check_accum if 'A'

      where (atb532_cal(:ncol,:nhtml_cosp*nscol_cosp) .eq. -1.00000E+30_r8)
            atb532_cal(:ncol,:nhtml_cosp*nscol_cosp) = fillvalue
      end where  
      call outfld('ATB532_CAL',atb532_cal    ,pcols,lchnk) !! fails check_accum if 'A'

      where (mol532_cal(:ncol,:nhtml_cosp) .eq. -1.00000E+30_r8)
            mol532_cal(:ncol,:nhtml_cosp) = fillvalue
      end where  
      call outfld('MOL532_CAL',mol532_cal    ,pcols,lchnk)

      where (cfad_sr532_cal(:ncol,:nht_cosp*nsr_cosp) .eq. -1.00000E+30_r8)
            cfad_sr532_cal(:ncol,:nht_cosp*nsr_cosp) = fillvalue
      end where  
      call outfld('CFAD_SR532_CAL',cfad_sr532_cal    ,pcols,lchnk)

      where (refl_parasol(:ncol,:nsza_cosp) .eq. -1.00000E+30_r8)
	    !! setting missing values to 0 (clear air).  
	    !! I'm not sure why COSP produces a mix of fillvalue and realvalue in the nht_cosp dimension.
            refl_parasol(:ncol,:nsza_cosp) = 0
      end where  
      call outfld('RFL_PARASOL',refl_parasol   ,pcols,lchnk) !! fails check_accum if 'A'
   end if

!!! RADAR SIMULATOR OUTPUTS
   if (lradar_sim) then
      where (cfad_dbze94_cs(:ncol,:nht_cosp*ndbze_cosp) .eq. -1.00000E+30_r8)
            cfad_dbze94_cs(:ncol,:nht_cosp*ndbze_cosp) = fillvalue
      end where  
      call outfld('CFAD_DBZE94_CS',cfad_dbze94_cs    ,pcols,lchnk)

      where (dbze_cs(:ncol,:nhtml_cosp*nscol_cosp) .eq. -1.00000E+30_r8)
            dbze_cs(:ncol,:nhtml_cosp*nscol_cosp) = fillvalue
      end where  
      call outfld('DBZE_CS',dbze_cs   ,pcols,lchnk) !! fails check_accum if 'A'

      where (cldtot_calcs(1:ncol) .eq. -1.00000E+30_r8)
            cldtot_calcs(1:ncol) = fillvalue
      end where  
      call outfld('CLDTOT_CALCS',cldtot_calcs    ,pcols,lchnk)

      where (cld_cal_notcs(:ncol,:nht_cosp) .eq. -1.00000E+30_r8)
            cld_cal_notcs(:ncol,:nht_cosp)= fillvalue
      end where  
      call outfld('CLD_CAL_NOTCS',cld_cal_notcs    ,pcols,lchnk)
   end if

!!! MISR SIMULATOR OUTPUTS
   if (lmisr_sim) then
      where (cld_misr(:ncol,:nhtmisr_cosp*ntau_cosp)  .eq. -1.00000E+30_r8)
            cld_misr(:ncol,:nhtmisr_cosp*ntau_cosp) = fillvalue
      end where  
      call outfld('CLD_MISR',cld_misr    ,pcols,lchnk)
   end if

!!! SUB-COLUMN OUTPUT
   if (lfrac_out) then
      where (scops_out(:ncol,:nhtml_cosp*nscol_cosp) .eq. -1.00000E+30_r8)
            scops_out(:ncol,:nhtml_cosp*nscol_cosp) = fillvalue
      end where  
      call outfld('SCOPS_OUT',scops_out   ,pcols,lchnk)!!!-1.00000E+30 !! fails check_accum if 'A'
   end if

!!! MODIS SIMULATOR OUTPUTS
   if (lmodis_sim) then

      where (cltmodis(:ncol) .eq. -1.00000E+30_r8)
            cltmodis(:ncol) = fillvalue
      end where  
      call outfld('CLTMODIS',cltmodis    ,pcols,lchnk)

      where (clwmodis(:ncol) .eq. -1.00000E+30_r8)
            clwmodis(:ncol) = fillvalue
      end where  
      call outfld('CLWMODIS',clwmodis    ,pcols,lchnk)

      where (climodis(:ncol) .eq. -1.00000E+30_r8)
            climodis(:ncol) = fillvalue
      end where  
      call outfld('CLIMODIS',climodis    ,pcols,lchnk)

      where (clhmodis(:ncol) .eq. -1.00000E+30_r8)
            clhmodis(:ncol) = fillvalue
      end where  
      call outfld('CLHMODIS',clhmodis    ,pcols,lchnk)

      where (clmmodis(:ncol) .eq. -1.00000E+30_r8)
            clmmodis(:ncol) = fillvalue
      end where  
      call outfld('CLMMODIS',clmmodis    ,pcols,lchnk)

      where (cllmodis(:ncol) .eq. -1.00000E+30_r8)
            cllmodis(:ncol) = fillvalue
      end where  
      call outfld('CLLMODIS',cllmodis    ,pcols,lchnk)

      where (tautmodis(:ncol) .eq. -1.00000E+30_r8)
            tautmodis(:ncol) = fillvalue
      end where  
      call outfld('TAUTMODIS',tautmodis    ,pcols,lchnk)

      where (tauwmodis(:ncol) .eq. -1.00000E+30_r8)
            tauwmodis(:ncol) = fillvalue
      end where  
      call outfld('TAUWMODIS',tauwmodis    ,pcols,lchnk)

      where (tauimodis(:ncol) .eq. -1.00000E+30_r8)
            tauimodis(:ncol) = fillvalue
      end where  
      call outfld('TAUIMODIS',tauimodis    ,pcols,lchnk)

      where (tautlogmodis(:ncol) .eq. -1.00000E+30_r8)
            tautlogmodis(:ncol) = fillvalue
      end where  
      call outfld('TAUTLOGMODIS',tautlogmodis    ,pcols,lchnk)

      where (tauwlogmodis(:ncol) .eq. -1.00000E+30_r8)
            tauwlogmodis(:ncol) = fillvalue
      end where  
      call outfld('TAUWLOGMODIS',tauwlogmodis    ,pcols,lchnk)

      where (tauilogmodis(:ncol) .eq. -1.00000E+30_r8)
            tauilogmodis(:ncol) = fillvalue
      end where  
      call outfld('TAUILOGMODIS',tauilogmodis    ,pcols,lchnk)

      where (reffclwmodis(:ncol) .eq. -1.00000E+30_r8)
            reffclwmodis(:ncol) = fillvalue
      end where  
      call outfld('REFFCLWMODIS',reffclwmodis    ,pcols,lchnk)

      where (reffclimodis(:ncol) .eq. -1.00000E+30_r8)
            reffclimodis(:ncol) = fillvalue
      end where  
      call outfld('REFFCLIMODIS',reffclimodis    ,pcols,lchnk)

      where (pctmodis(:ncol) .eq. -1.00000E+30_r8)
            pctmodis(:ncol) = fillvalue
      end where  
      call outfld('PCTMODIS',pctmodis    ,pcols,lchnk)

      where (lwpmodis(:ncol) .eq. -1.00000E+30_r8)
            lwpmodis(:ncol) = fillvalue
      end where  
      call outfld('LWPMODIS',lwpmodis    ,pcols,lchnk)

      where (iwpmodis(:ncol) .eq. -1.00000E+30_r8)
            iwpmodis(:ncol) = fillvalue
      end where  
      call outfld('IWPMODIS',iwpmodis    ,pcols,lchnk)

      where (clmodis_cam(:ncol,:nprs_cosp*ntau_cosp_modis) .eq. -1.00000E+30_r8)
            clmodis_cam(:ncol,:nprs_cosp*ntau_cosp_modis) = fillvalue
      end where  
      call outfld('CLMODIS',clmodis_cam    ,pcols,lchnk)

   end if

#endif
! USE_COSP

end subroutine cospsimulator_intr_run

!#######################################################################

end module cospsimulator_intr
