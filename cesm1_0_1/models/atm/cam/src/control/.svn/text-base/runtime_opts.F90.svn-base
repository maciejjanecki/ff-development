
module runtime_opts

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for reading CAM namelist cam_inparm 
!          and broadcasting namelist values if needed.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, September 2003
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8, SHR_KIND_CL
use spmd_utils,      only: masterproc
use namelist_utils,  only: find_group_name
use pmgrid,          only: plat, plev, plon
use cam_history
use cam_control_mod
use cam_diagnostics, only: inithist_all
use cam_logfile,     only: iulog
use pspect
!use shr_orb_mod
use units
use constituents,    only: pcnst, readtrace
use tracers,         only: tracers_flag
use time_manager,    only: dtime
use filenames,       only: ncdata, bnd_topo, &
                           absems_data, modal_optics, &
                           caseid, isccpdata, &
                           brnch_retain_casename
use cloudsimulator,  only: doisccp
use cloudsimulator_38, only: doisccp_38
use dycore,          only: dycore_is
use abortutils,      only: endrun
use fv_control_mod,  only: nsplit, nspltrac, nspltvrm, iord, jord, kord, dyn_conservative,    &
                           filtcw, ct_overlap, trac_decomp, fft_flt, div24del2flag, del2coef
use rayleigh_friction, only: rayk0, raykrange, raytau0

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
implicit none
private
save


!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
public read_namelist        ! Set and/or get all runtime options

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

character(len=SHR_KIND_CL), private :: nlfilename = 'atm_in' ! Namelist filename

!-----------------------------------------------------------------------
!
! SOMEWHAT ALPHABETICAL listing of variables in the cam_inparm namelist:
!
! variable                description
! --------             -----------------
!
! bnd_topo             Path and filename of topography dataset
! 
! absems_data          Dataset with absorption and emissivity factors.
!
! modal_optics         Dataset with coefficients for modal aerosol optics.
!
!
! dif2 = nnn.n,        del2 horizontal diffusion coeff. Default value 
!                      defined in module comhd.  
real(r8) :: dif2
! dif4 = nnn.n,        del4 horizontal diffusion coeff. Default value 
!                      defined in module comhd.  
real(r8) :: dif4
! 
! kmxhdc = nn          number of levels (starting from model top) to
!                      apply Courant limiter.  Default value defined 
!                      in module comhd.  
integer :: kmxhdc
! 
! divdampn = 0.        Number of days (from nstep 0) to run divergence
!                      damper
!
! dtime = nnnn,        Model time step in seconds. Default is dycore dependent.
! 
! eps = nnn.n,         time filter coefficient. Defaults to 0.06.
! 
! fincl1 = 'field1', 'field2',...
!                      List of fields to add to the primary history file.
! fincl1lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl1 fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      single character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl1 fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fincl[2..6] = 'field1', 'field2',...
!                      List of fields to add to the auxiliary history file.
!
! fincl2..6]lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl[2..6] fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      singel character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl[2..6] fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fexcl1 = 'field1','field2',... 
!                      List of field names to exclude from default
!                      primary history file (default fields on the 
!                      Master Field List).
! 
! fexcl[2..6] = 'field1','field2',... 
!                      List of field names to exclude from
!                      auxiliary history files.
! 
! fhstpr1 = 'field1', 'field2',...
!                      List of fields to change buffer size in
!                      primary history file
!
! fhstpr[2..6] = 'field1', 'field2',...
!                      List of fields to change buffer size in auxiliary files
!
! fwrtpr1 = 'field1', 'field2',...
!                      List of fields to change output data type in
!                      primary history file
!
! fwrtpr[2..6] = 'field1', 'field2',...
!                      List of fields to change output data type in
!                      auxiliary files
!
! mfilt = nn,nn,nn     Array containing the maximum number of time 
!                      samples per disk history file. Defaults to 5.
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! ncdata               Path and filename of initial condition dataset.
! 
! nhtfrq = nn,nn,nn,.. Output history frequency for each tape
!
!                      If = 0 : monthly average
!                      If > 0 : output every nhtfrq time steps.
!                      If < 0 : output every abs(nhtfrq) hours.
! 
! nlvdry = nn,         Number of layers over which to do dry
!                      adjustment. Defaults to 3.
! 
! cam_branch_file      Filepath of restart file to branch from (nsrest=3)
!                      Full pathname required.
character(len=256) :: cam_branch_file = ' '
!
! use_64bit_nc         True if new 64-bit netCDF formit, false otherwise (default false)
! 
!------------------------------------------------------------------
! The following are specific to f-v dynamics (see dynpkg for info)
!------------------------------------------------------------------
! nsplit               Intended number of dynamics timesteps per physics timestep
! nspltrac             Intended number of tracer advection timesteps per physics timestep
!                      Nsplit is partitioned into nspltrac and nsplit/nspltrac,
!                      with the latter being the number of dynamics timesteps per
!                      tracer timestep, possibly rounded upward; after initialization,
!                      the code quantity nsplit is redefined to be the number of
!                      dynamics timesteps per tracer timestep.
! nspltvrm             Number of vertical re-mapping timesteps per physics timestep
!                      Nspltrac is partitioned into nspltvrm and nspltrac/nspltvrm,
!                      with the latter being the number of tracer timesteps per
!                      re-mapping timestep, possibly rounded upward; after initialization,
!                      the code quantity nspltrac is redefined to be the number of
!                      tracer timesteps per re-mapping timestep.
! iord                 scheme to be used for E-W transport (default: 4)
! jord                 scheme to be used for N-S transport (default: 4)
! kord                 scheme to be used for vertical mapping (default: 4)
! filtcw               filter c-grid winds if nonzero
! ct_overlap           Nonzero for overlapping cd_core and trac2d, 0 otherwise
! trac_decomp          Size of trac2d domain decomposition (default 1)
! fft_flt              1 for FFT filter always, 0 for algebraic/FFT filter
!                      (default: 0 for cam3, cam3_5; 1 otherwise)
! ldiv2                .true. : CAM3, CAM3.5 default divergence damping
!                      .false.: 4th-order divergence damping except model top
!                               layers where 2nd-order Laplacian damping is used
!------------------------------------------------------------------
!
!------------------------------------------------------------------
! The following 29 are specific to f-v decomposition and transposes 
! (see spmd_dyn for info)
!------------------------------------------------------------------
! npr_yz(4)            yz and xy decompositions
integer :: npr_yz(4)
! geopktrans           geopotential method (routine geopk, geopk16, 
!                       or geopk_d)
integer :: geopktrans
! geopkblocks          number of stages to use in Z-serial non-transpose
!                       geopotential method (routine geopk_d)
integer :: geopkblocks
! force_2d             option to force transpose computation for 1D decomp.
integer :: force_2d
! modcomm_transpose    mod_comm transpose method
integer :: modcomm_transpose
! modcomm_geopk        mod_comm geopk method
integer :: modcomm_geopk
! modcomm_gatscat      mod_comm gather/scatter method
integer :: modcomm_gatscat
! modc_sw_dynrun       mod_comm irregular underlying communication method for dyn_run/misc
integer :: modc_sw_dynrun
! modc_hs_dynrun       mod_comm irregular communication handshaking for dyn_run/misc
logical :: modc_hs_dynrun
! modc_send_dynrun     mod_comm irregular communication send protocol for dyn_run/misc
logical :: modc_send_dynrun
! modc_mxreq_dynrun   mod_comm irregular communication nonblocking request throttle for dyn_run/misc
integer :: modc_mxreq_dynrun
! modc_sw_cdcore       mod_comm irregular underlying communication method for cd_core/geopk
integer :: modc_sw_cdcore
! modc_hs_cdcore       mod_comm irregular communication handshaking for cd_core/geopk
logical :: modc_hs_cdcore
! modc_send_cdcore     geopk_d and mod_comm irregular communication send protocol for cd_core/geopk
logical :: modc_send_cdcore
! modc_mxreq_cdcore   mod_comm irregular communication nonblocking request throttle for cd_core/geopk
integer :: modc_mxreq_cdcore
! modc_sw_gather      mod_comm irregular underlying communication method for gather
integer :: modc_sw_gather
! modc_hs_gather      mod_comm irregular communication handshaking for gather
logical :: modc_hs_gather
! modc_send_gather    mod_comm irregular communication send protocol for gather
logical :: modc_send_gather
! modc_mxreq_gather   mod_comm irregular communication nonblocking request throttle for gather
integer :: modc_mxreq_gather
! modc_sw_scatter      mod_comm irregular underlying communication method for scatter
integer :: modc_sw_scatter
! modc_hs_scatter      mod_comm irregular communication handshaking for scatter
logical :: modc_hs_scatter
! modc_send_scatter    mod_comm irregular communication send protocol for scatter
logical :: modc_send_scatter
! modc_mxreq_scatter   mod_comm irregular communication nonblocking request throttle for scatter
integer :: modc_mxreq_scatter
! modc_sw_tracer       mod_comm irregular underlying communication method for multiple tracers
integer :: modc_sw_tracer
! modc_hs_tracer       mod_comm irregular communication handshaking for multiple tracers
logical :: modc_hs_tracer
! modc_send_tracer     mod_comm irregular communication send protocol for multiple tracers
logical :: modc_send_tracer
! modc_mxreq_tracer   mod_comm irregular communication nonblocking request throttle for multiple tracers
integer :: modc_mxreq_tracer
! modc_onetwo          one or two simultaneous mod_comm irregular communications (excl. tracers)
integer :: modc_onetwo
! modc_tracers         max number of tracers for simultaneous mod_comm irregular communications
!                      if 0, then use mp_sendirr/mp_recvirr to do one at a time
integer :: modc_tracers
!------------------------------------------------------------------

!------------------------------------------------------------------
! The following 5 are specific to eul/sld communication algorithms
! (see { eul | sld }/spmd_dyn for info)
!------------------------------------------------------------------
! dyn_alltoall         dynamics transpose option.
integer :: dyn_alltoall
! dyn_allgather        dynamics gather option.
integer :: dyn_allgather
! dyn_equi_by_col      dynamics load balancing option.
logical :: dyn_equi_by_col
! dyn_npes             number of processes assigned to dynamics (EUL and SLD dycores)
integer :: dyn_npes
! dyn_npes_stride      stride for dynamics processes (EUL and SLD dycores)
!                      (e.g., if stride=2, assign every second process to the dynamics)
integer :: dyn_npes_stride
!------------------------------------------------------------------
! The following 3 are specific to the spmd_utils module, used
! in the point-to-point implementations of eul/sld and physics
! communication algorithms and in the flow-controlled gather collective
! (see spmd_utils for info)
!------------------------------------------------------------------
! swap_comm_protocol   Performance tuning option for swap communication.
integer :: swap_comm_protocol
! swap_comm_maxreq     Performance tuning option for swap communication.
integer :: swap_comm_maxreq
! fc_gather_flow_cntl  Tuning option for flow-controlled gather, used
!                      to improve robustness at high process count and
!                      to improve performance
integer :: fc_gather_flow_cntl
!------------------------------------------------------------------
! The following 3 are specific to Rayleigh friction
! integer rayk0         vertical level at which rayleigh friction term is centered
! real(r8) raykrange    range of rayleigh friction profile; if 0, range is set automatically
! real(r8) raytau0      approximate value of decay time at model top (days);
!                       if 0., no rayleigh friction is applied
!------------------------------------------------------------------
!
!
! hfilename_spec       Flexible filename specifier for history files
!
! 
! pertlim = n.n        Max size of perturbation to apply to initial
!                      temperature field.
!
! phys_alltoall        Dynamics/physics transpose option. See phys_grid module.
!
integer :: phys_alltoall
! 
! phys_loadbalance     Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_loadbalance
! 
! phys_twin_algorithm  Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_twin_algorithm
! 
! phys_chnk_per_thd    Performance tuning option for physics chunks.  See 
!                      phys_grid module.  
integer :: phys_chnk_per_thd
!
! repro_sum_use_ddpdd  Flag indicating that the double-double summation 
!                      algorithm should be used instead of the repro 
!                      fixed precision algorithm
logical :: repro_sum_use_ddpdd
!
! repro_sum_rel_diff_max Maximum permissible relative difference between
!                      repro and nonrepro scalable algorithms. 
!                      See repro_sum_mod  module.
real(r8) :: repro_sum_rel_diff_max
!
! repro_sum_recompute  Flag indicating that an alternative, serial, algorithm
!                      should be used when the relative difference between
!                      the repro and nonrepro scalable algorithms exceeds
!                      the tolerance specified by repro_sum_rel_diff_max.
!                      See repro_sum_mod  module.
logical :: repro_sum_recompute
! 
! tracers_flag = .F.    If true, implement tracer test code. Number of tracers determined
!                      in tracers_suite.F90 must agree with PCNST
!
! readtrace = .T.      If true, tracer initial conditions obtained from 
!                      initial file. 
!
! inithist             Generate initial dataset as auxillary history file
!                      can be set to '6-HOURLY', 'DAILY', 'MONTHLY', 'YEARLY' or 'NONE'. 
!                      default: 'YEARLY'
!
! empty_htapes         true => no fields by default on history tapes
!
! print_step_cost      true => print per timestep cost info
!
! avgflag_pertape      A, I, X, or M means avg, instantaneous, max or min for all fields on
!                      that tape
!
! doisccp              whether to do ISCCP calcs (version 3.4) and history output (default false)
! doisccp_38           whether to do ISCCP calcs (version 3.8) and history output (default false)
!
!   logical indirect     
!                    ! true => include indirect radiative effects of
!                    ! sulfate aerosols.  Default is false.
!
! inithist_all         .false.:  include only REQUIRED fields on IC file
!                      .true. :  include required AND optional fields on IC file
!                      default:  .false.
!
! met_data_file        name of file that contains the offline meteorology data
!
! met_filenames_list   name of file that contains names of the offline 
!                      meteorology data files
!
! met_remove_file      true => the offline meteorology file will be removed
!
! met_cell_wall_winds  true => the offline meteorology winds are defined on the model
!                      grid cell walls
! Physics buffer
logical :: pbuf_global_allocate       ! allocate all buffers as global (default: .true.)


! Diagnostics options

character(len=8) :: diag_cnst_conv_tend ! output constituent tendencies due to convection
                                        ! 'none', 'q_only' or 'all'

! Conservation checks

logical            :: print_energy_errors ! switch for diagnostic output from check_energy module

! Radiative heating rate calculation options

integer :: iradsw        ! freq. of shortwave radiation calc in time steps (positive)
                         ! or hours (negative).  Default: -1
integer :: iradlw        ! frequency of longwave rad. calc. in time steps (positive)
                         ! or hours (negative).  Default: -1
integer :: iradae        ! frequency of absorp/emis calc in time steps (positive)
                         ! or hours (negative).  Default: -12
integer :: irad_always   ! Specifies length of time in timesteps (positive)
                         ! or hours (negative) SW/LW radiation will be run continuously
                         ! from the start of an initial run.  Default: 0

#if (defined WACCM_PHYS)
! iondrag / efield
character(len=256) :: efield_lflux_file
character(len=256) :: efield_hflux_file
character(len=256) :: efield_wei96_file
! waccm qbo data variables
character(len=256) :: qbo_forcing_file
logical            :: qbo_use_forcing
logical            :: qbo_cyclic
#endif

! Upper atmosphere radiative processes (waccm phys)
logical :: nlte_use_mo              ! Determines which constituents are used from NLTE calculations
                                    !  = .true. uses MOZART constituents
                                    !  = .false. uses constituents from bnd dataset cftgcm
integer :: itgcmcyc                 ! flag for cycling TIME/GCM input dataset:
                                    !  = 0 : one time sample only
                                    !  = 1 : full annual cycle (default)
                                    !  = 2 : two time samples
character(len=256) :: cftgcm        ! Pathname of time-variant TIME/GCM output

! SCM Options
logical  :: single_column
real(r8) :: scmlat,scmlon
integer, parameter :: max_chars = 128
character(len=max_chars) iopfile
logical  :: scm_iop_srf_prop
logical  :: scm_relaxation
logical  :: scm_diurnal_avg
logical  :: scm_crm_mode

#if ( defined OFFLINE_DYN )
logical :: met_remove_file
logical :: met_cell_wall_winds
character(len=256) :: met_data_file
character(len=256) :: met_filenames_list
real(r8) :: met_rlx_top ! (km) top of relaxation region of winds for offline waccm
real(r8) :: met_rlx_bot ! (km) bottom of relaxation region of winds for offline waccm
real(r8) :: met_max_rlx ! maximum of vertical relaxation function in bottom portion (default is 1.0)
logical  :: met_fix_mass ! switch to turn on/off mass fixer for offline driver (default is TRUE)
character(len=16)  :: met_shflx_name ! srf heat flux field name in met data file
character(len=16)  :: met_qflx_name  ! water vapor flux field name in met data file
real(r8)           :: met_shflx_factor ! multiplication factor for srf heat flux 
real(r8)           :: met_qflx_factor  ! multiplication factor for water vapor flux 
#endif

contains

!=======================================================================

  subroutine read_namelist(single_column_in, scmlon_in, scmlat_in, nlfilename_in )

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Read data from namelist cam_inparm to define the run. Process some of the
   ! namelist variables to determine history and restart/branch file path 
   ! names.  Check input namelist variables for validity and print them
   ! to standard output. 
   ! 
   ! Method: 
   ! Important Note for running on SUN systems: "implicit automatic (a-z)"
   ! will not work because namelist data must be static.
   !
   ! Author: 
   ! Original version:  CCM1
   ! Standardized:      L. Bath, June 1992
   !                    T. Acker, March 1996
   !     
   !-----------------------------------------------------------------------

   ! Note that the following interfaces are prototypes proposed by Henderson 
   ! and Eaton.  They minimize coupling with other modules.  Design of these 
   ! interfaces should be refined via review by other CAM developers.  
   ! Interface *_defaultopts() gets default values from the responsible 
   ! module (Expert) prior to namelist read.  
   ! Interface *_setopts() sends values to the responsible module (Expert) 
   ! after namelist read.  Erroneous values are handled by Experts.  
   ! TBH  9/8/03 
   !
   use phys_grid,        only: phys_grid_defaultopts, phys_grid_setopts
   use repro_sum_mod,    only: repro_sum_defaultopts, repro_sum_setopts
   use phys_buffer,      only: pbuf_defaultopts, pbuf_setopts
#if ( defined SPMD )
   use spmd_utils,       only: swap_comm_defaultopts, swap_comm_setopts, &
                               fc_gather_defaultopts, fc_gather_setopts
   use spmd_dyn,         only: spmd_dyn_defaultopts, spmd_dyn_setopts
#endif
#if (defined WACCM_PHYS)
   use iondrag,          only: iondrag_defaultopts, iondrag_setopts
   use qbo,              only: qbo_defaultopts, qbo_setopts
#endif

   use chem_surfvals,    only: chem_surfvals_readnl
   use cam_diagnostics,  only: diag_defaultopts, diag_setopts
   use check_energy,     only: check_energy_defaultopts, check_energy_setopts
   use comhd,            only: comhd_defaultopts, comhd_setopts
   use radiation,        only: radiation_defaultopts, radiation_setopts, radiation_printopts
   use cam_restart,      only: restart_defaultopts, restart_setopts, restart_printopts
   use radheat,          only: radheat_defaultopts, radheat_setopts
#if ( defined OFFLINE_DYN )
   use metdata,          only: offline_met_defaultopts, offline_met_setopts
#endif
   use co2_cycle,        only: co2_cycle_readnl
   use shr_string_mod,   only: shr_string_toUpper
   use scamMod,          only: scam_setopts,scam_default_opts

   ! Some modules read their own namelist input.
   use physconst,           only: physconst_readnl
   use phys_control,        only: phys_ctl_readnl
   use cam3_aero_data,      only: cam3_aero_data_readnl
   use cam3_ozone_data,     only: cam3_ozone_data_readnl
   use cloud_fraction,      only: cldfrc_readnl
   use cldwat,              only: cldwat_readnl
   use zm_conv,             only: zmconv_readnl
   use hk_conv,             only: hkconv_readnl
   use pkg_cld_sediment,    only: cld_sediment_readnl
   use gw_drag,             only: gw_drag_readnl
   use phys_debug_util,     only: phys_debug_readnl
   use rad_constituents,    only: rad_cnst_readnl
   use radiation_data,      only: rad_data_readnl
   use chemistry,           only: chem_readnl
   use prescribed_volcaero, only: prescribed_volcaero_readnl
   use aerodep_flx,         only: aerodep_flx_readnl
   use solar_data,          only: solar_data_readnl
   use tropopause,          only: tropopause_readnl
   use aoa_tracers,         only: aoa_tracers_readnl
   use prescribed_ozone,    only: prescribed_ozone_readnl
   use prescribed_aero,     only: prescribed_aero_readnl
   use prescribed_ghg,      only: prescribed_ghg_readnl
   use aircraft_emit,       only: aircraft_emit_readnl
   use cospsimulator_intr,  only: cospsimulator_intr_readnl

!---------------------------Arguments-----------------------------------

   logical , intent(in), optional :: single_column_in 
   real(r8), intent(in), optional :: scmlon_in
   real(r8), intent(in), optional :: scmlat_in
   character(len=*)    , optional :: nlfilename_in
!-----------------------------------------------------------------------

   include 'netcdf.inc'

!---------------------------Local variables-----------------------------
   character(len=*), parameter ::  subname = "read_namelist"
! 
   character ctemp*8      ! Temporary character strings
   integer ntspdy         ! number of timesteps per day
   integer t              ! history tape index
   integer lastchar       ! index to last char of a char variable
   integer ierr           ! error code
   integer unitn          ! namelist unit number

   integer f, i
   integer, parameter :: fieldname_len = 16
   integer, parameter :: fieldname_lenp2 = fieldname_len + 2
   integer, parameter :: max_chars = 128

   character(len=fieldname_lenp2) fincl1(pflds)
   character(len=fieldname_lenp2) fincl2(pflds)
   character(len=fieldname_lenp2) fincl3(pflds)
   character(len=fieldname_lenp2) fincl4(pflds)
   character(len=fieldname_lenp2) fincl5(pflds)
   character(len=fieldname_lenp2) fincl6(pflds)

   character(len=max_chars) fincl1lonlat(pflds)
   character(len=max_chars) fincl2lonlat(pflds)
   character(len=max_chars) fincl3lonlat(pflds)
   character(len=max_chars) fincl4lonlat(pflds)
   character(len=max_chars) fincl5lonlat(pflds)
   character(len=max_chars) fincl6lonlat(pflds)

   character(len=fieldname_len) fexcl1(pflds)
   character(len=fieldname_len) fexcl2(pflds)
   character(len=fieldname_len) fexcl3(pflds)
   character(len=fieldname_len) fexcl4(pflds)
   character(len=fieldname_len) fexcl5(pflds)
   character(len=fieldname_len) fexcl6(pflds)

   character(len=fieldname_lenp2) fhstpr1(pflds)
   character(len=fieldname_lenp2) fhstpr2(pflds)
   character(len=fieldname_lenp2) fhstpr3(pflds)
   character(len=fieldname_lenp2) fhstpr4(pflds)
   character(len=fieldname_lenp2) fhstpr5(pflds)
   character(len=fieldname_lenp2) fhstpr6(pflds)

   character(len=fieldname_lenp2) fwrtpr1(pflds)
   character(len=fieldname_lenp2) fwrtpr2(pflds)
   character(len=fieldname_lenp2) fwrtpr3(pflds)
   character(len=fieldname_lenp2) fwrtpr4(pflds)
   character(len=fieldname_lenp2) fwrtpr5(pflds)
   character(len=fieldname_lenp2) fwrtpr6(pflds)

!
! Define the cam_inparm namelist
!
! TBH:  NOTE:  Move the definition of cam_inparm outside of this routine 
! TBH:  as soon as common blocks in comctl.h, comtfc.h, 
! TBH:  comsol.h, comadj.h, and perturb.h have been converted to 
! TBH:  modules.  
!        
! ***NOTE*** If a namelist option is not described in the CAM Users Guide,
!            it is not supported.  

  namelist /cam_inparm/ ncdata, bnd_topo, &
                    cam_branch_file  ,nhstpr  ,ndens   ,nhtfrq  , &
                    mfilt   ,absems_data ,modal_optics, &
                    fincl1  ,fincl2  ,fincl3  ,fincl4  ,fincl5  , &
                    fincl1lonlat,fincl2lonlat,fincl3lonlat, &
                    fincl4lonlat  ,fincl5lonlat  , &
                    fincl6  ,fexcl1  ,fexcl2  ,fexcl3  ,fexcl4  , &
                    fexcl5  ,fexcl6  ,hfilename_spec, &
                    fhstpr1 ,fhstpr2 ,fhstpr3 ,fhstpr4 ,fhstpr5 ,fhstpr6 , &
                    fwrtpr1 ,fwrtpr2 ,fwrtpr3, fwrtpr4 ,fwrtpr5 ,fwrtpr6 , &
                    dtime,  eps     ,dif2    ,dif4    ,kmxhdc  , &
                    nlvdry,  &
                    pertlim ,divdampn, &
                    readtrace, rayk0, raykrange, raytau0, &
                    tracers_flag, &
                    inithist, indirect, nsplit, nspltrac, nspltvrm, filtcw, fft_flt, &
                    iord, jord, kord, dyn_conservative, div24del2flag, del2coef, &
                    npr_yz, geopktrans, geopkblocks, &
                    force_2d, ct_overlap, trac_decomp, &
                    modcomm_transpose, modcomm_geopk, modcomm_gatscat, &
                    modc_sw_dynrun, modc_hs_dynrun, modc_send_dynrun, modc_mxreq_dynrun, &
                    modc_sw_cdcore, modc_hs_cdcore, modc_send_cdcore, modc_mxreq_cdcore, &
                    modc_sw_gather, modc_hs_gather, modc_send_gather, modc_mxreq_gather, &
                    modc_sw_scatter, modc_hs_scatter, modc_send_scatter, modc_mxreq_scatter, &
                    modc_sw_tracer, modc_hs_tracer, modc_send_tracer, modc_mxreq_tracer, &
                    modc_onetwo, modc_tracers, &
                    dyn_alltoall, dyn_allgather, dyn_equi_by_col, &
                    dyn_npes, dyn_npes_stride, &
                    swap_comm_protocol, swap_comm_maxreq, &
                    fc_gather_flow_cntl, &
                    empty_htapes, use_64bit_nc, &
                    print_step_cost, avgflag_pertape, &
                    doisccp, doisccp_38, isccpdata, &
                    phys_alltoall, phys_loadbalance, phys_twin_algorithm, &
                    phys_chnk_per_thd, &
                    repro_sum_use_ddpdd, repro_sum_rel_diff_max, repro_sum_recompute, &
                    inithist_all

  ! physics buffer
  namelist /cam_inparm/ pbuf_global_allocate

  ! diagnostic options
  namelist /cam_inparm/ diag_cnst_conv_tend

  ! conservation checks
  namelist /cam_inparm/ print_energy_errors

  ! radiative heating calculation options
  namelist /cam_inparm/ iradsw, iradlw, iradae, irad_always

#if (defined WACCM_PHYS)
  ! iondrag / efield options
  namelist /cam_inparm/ efield_lflux_file, efield_hflux_file, efield_wei96_file
  ! waccm qbo namelist variables
  namelist /cam_inparm/ qbo_use_forcing, qbo_forcing_file, qbo_cyclic
#endif

  ! upper atmosphere radiative processes
  namelist /cam_inparm/ nlte_use_mo, itgcmcyc, cftgcm

  ! scam
  namelist /cam_inparm/ iopfile,scm_iop_srf_prop,scm_relaxation, &
                        scm_diurnal_avg,scm_crm_mode

#if ( defined OFFLINE_DYN )
  ! offline meteorology parameters
  namelist /cam_inparm/ met_data_file, met_remove_file, met_cell_wall_winds, &
                        met_filenames_list, met_rlx_top, met_rlx_bot, met_max_rlx, &
                        met_fix_mass, met_shflx_name, met_qflx_name, &
                        met_shflx_factor, met_qflx_factor
#endif
! 
!-----------------------------------------------------------------------
  if (present(nlfilename_in)) then
     nlfilename = nlfilename_in
  end if
!
! Determine preset values (this is currently being phased out)
!
   call preset ()
!
! Finite volume code only: Set Lagrangian time splits.  A default of zero indicates the number
! should be automatically computed unless the user enters something.
!
   ! This is a hack to set nsplit for the waccm model at resolutions of 2x2.5 and finer.
   if (plev > 65 .and. plat > 90) then
      nsplit = 8
   end if
!
! Preset sulfate aerosol related variables

   indirect  = .false.
! 
! Get default values of runtime options for spmd_dyn
!
#if ( defined SPMD )
   if ( dycore_is ('UNSTRUCTURED') ) then
      call spmd_dyn_defaultopts()
   else if ( dycore_is ('LR') ) then
      call spmd_dyn_defaultopts(                       &
             npr_yz_out            =npr_yz,            &
             geopktrans_out        =geopktrans,        &
             geopkblocks_out       =geopkblocks,       &
             force_2d_out          =force_2d,          &
             modcomm_transpose_out =modcomm_transpose, &
             modcomm_geopk_out     =modcomm_geopk,     &
             modcomm_gatscat_out   =modcomm_gatscat,   &
             modc_sw_dynrun_out    =modc_sw_dynrun,    &
             modc_hs_dynrun_out    =modc_hs_dynrun,    &
             modc_send_dynrun_out  =modc_send_dynrun,  &
             modc_mxreq_dynrun_out =modc_mxreq_dynrun, &
             modc_sw_cdcore_out    =modc_sw_cdcore,    &
             modc_hs_cdcore_out    =modc_hs_cdcore,    &
             modc_send_cdcore_out  =modc_send_cdcore,  &
             modc_mxreq_cdcore_out =modc_mxreq_cdcore, &
             modc_sw_gather_out    =modc_sw_gather,    &
             modc_hs_gather_out    =modc_hs_gather,    &
             modc_send_gather_out  =modc_send_gather,  &
             modc_mxreq_gather_out =modc_mxreq_gather, &
             modc_sw_scatter_out   =modc_sw_scatter,   &
             modc_hs_scatter_out   =modc_hs_scatter,   &
             modc_send_scatter_out =modc_send_scatter, &
             modc_mxreq_scatter_out=modc_mxreq_scatter,&
             modc_sw_tracer_out    =modc_sw_tracer,    &
             modc_hs_tracer_out    =modc_hs_tracer,    &
             modc_send_tracer_out  =modc_send_tracer,  &
             modc_mxreq_tracer_out =modc_mxreq_tracer, &
             modc_onetwo_out       =modc_onetwo,       &
             modc_tracers_out      =modc_tracers       )
   else if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call spmd_dyn_defaultopts(                  &
             dyn_alltoall_out    =dyn_alltoall,   &
             dyn_allgather_out   =dyn_allgather,  &
             dyn_equi_by_col_out =dyn_equi_by_col,&
             dyn_npes_out        =dyn_npes,       &
             dyn_npes_stride_out =dyn_npes_stride )
   endif

   ! Get default values of runtime options for swap communication
   call swap_comm_defaultopts(                   &
      swap_comm_protocol_out=swap_comm_protocol, &
      swap_comm_maxreq_out=swap_comm_maxreq      )

   ! Get default values of runtime options for flow-controlled
   ! gather collective
   call fc_gather_defaultopts(                    &
      fc_gather_flow_cntl_out=fc_gather_flow_cntl  )

#endif

   ! restart write interval
   call restart_defaultopts( &
      cam_branch_file_out          =cam_branch_file            )

   ! Get default values of runtime options for physics chunking.
   call phys_grid_defaultopts(                      &
      phys_loadbalance_out    =phys_loadbalance,    &
      phys_twin_algorithm_out =phys_twin_algorithm, &
      phys_alltoall_out       =phys_alltoall,       &
      phys_chnk_per_thd_out   =phys_chnk_per_thd    )

   ! Get default values of runtime options for reproducible sum
   ! calculations.
   call repro_sum_defaultopts(                           &
      repro_sum_use_ddpdd_out=repro_sum_use_ddpdd,       &
      repro_sum_rel_diff_max_out=repro_sum_rel_diff_max, &
      repro_sum_recompute_out=repro_sum_recompute        )

   ! physics buffer
   call pbuf_defaultopts( &
      pbuf_global_allocate_out = pbuf_global_allocate )

   ! diagnostics
   call diag_defaultopts( &
      diag_cnst_conv_tend_out = diag_cnst_conv_tend )

   ! conservation
   call check_energy_defaultopts( &
      print_energy_errors_out = print_energy_errors )

   ! radiative heating calcs
   call radiation_defaultopts( &
      iradsw_out      = iradsw,     &
      iradlw_out      = iradlw,     &
      iradae_out      = iradae,     &
      irad_always_out = irad_always )

#if (defined WACCM_PHYS)
   ! iondrag / efield
   call iondrag_defaultopts( &
      efield_lflux_file_out =efield_lflux_file, &
      efield_hflux_file_out =efield_hflux_file, &
      efield_wei96_file_out =efield_wei96_file )
   ! qbo forcing
   call qbo_defaultopts( &
      qbo_use_forcing_out  = qbo_use_forcing, &
      qbo_forcing_file_out = qbo_forcing_file,&
      qbo_cyclic_out       = qbo_cyclic       )
#endif

   ! Upper atmosphere radiative processes
   call radheat_defaultopts( &
      nlte_use_mo_out =nlte_use_mo, &
      itgcmcyc_out    =itgcmcyc,    &
      cftgcm_out      =cftgcm       )

! 
! Get default values of runtime options for comhd
!
   call comhd_defaultopts(dif2_out  =dif2, &
                          dif4_out  =dif4, &
                          kmxhdc_out=kmxhdc)

#if ( defined OFFLINE_DYN )
!
! Get runtime defualts for the metdata module
!
   call offline_met_defaultopts( met_data_file_out = met_data_file, &
                                 met_remove_file_out = met_remove_file, &
                                 met_cell_wall_winds_out = met_cell_wall_winds, &
                                 met_filenames_list_out = met_filenames_list, &
                                 met_rlx_top_km_out = met_rlx_top, &
                                 met_rlx_bot_km_out = met_rlx_bot, &
                                 met_max_rlx_out = met_max_rlx, &
                                 met_fix_mass_out = met_fix_mass, &
                                 met_shflx_name_out = met_shflx_name, &
                                 met_qflx_name_out = met_qflx_name, &
                                 met_shflx_factor_out = met_shflx_factor, &
                                 met_qflx_factor_out = met_qflx_factor )
#endif

   if (present(single_column_in)) then
      call scam_default_opts(scmlat_out=scmlat,scmlon_out=scmlon, &
        single_column_out=single_column, &
        scm_iop_srf_prop_out=scm_iop_srf_prop,&
        scm_relaxation_out=scm_relaxation, &
        scm_diurnal_avg_out=scm_diurnal_avg, &
        scm_crm_mode_out=scm_crm_mode)
   end if

   do f = 1, pflds
      fincl1(f) = ' '         
      fincl2(f) = ' '         
      fincl3(f) = ' '         
      fincl4(f) = ' '         
      fincl5(f) = ' '         
      fincl6(f) = ' '         
      fincl1lonlat(f) = ' '
      fincl2lonlat(f) = ' '
      fincl3lonlat(f) = ' '
      fincl4lonlat(f) = ' '
      fincl5lonlat(f) = ' '
      fincl6lonlat(f) = ' '
      fexcl1(f) = ' '
      fexcl2(f) = ' '
      fexcl3(f) = ' '
      fexcl4(f) = ' '
      fexcl5(f) = ' '
      fexcl6(f) = ' '
      fhstpr1(f) = ' '
      fhstpr2(f) = ' '
      fhstpr3(f) = ' '
      fhstpr4(f) = ' '
      fhstpr5(f) = ' '
      fhstpr6(f) = ' '
      fwrtpr1(f) = ' '
      fwrtpr2(f) = ' '
      fwrtpr3(f) = ' '
      fwrtpr4(f) = ' '
      fwrtpr5(f) = ' '
      fwrtpr6(f) = ' '
   enddo

   ! Read in the cam_inparm namelist from input filename

   if (masterproc) then
      write(iulog,*) 'Read in cam_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for cam_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'cam_inparm', status=ierr)
      if (ierr == 0) then  ! found cam_inparm
         read(unitn, cam_inparm, iostat=ierr)  ! read the cam_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                          ' error condition for cam_inparm' )
         end if
      else
         call endrun(subname // ':: can''t find cam_inparm in file ' // trim(nlfilename))
      end if
      close( unitn )
      call freeunit( unitn )
      !
      ! Check CASE namelist variable
      !
      if (caseid==' ') then
         call endrun ('READ_NAMELIST: Namelist variable CASEID must be set')
      end if

      lastchar = len(caseid)
      if (caseid(lastchar:lastchar) /= ' ') then
         write(iulog,*)'READ_NAMELIST: CASEID must not exceed ', len(caseid)-1, ' characters'
         call endrun
      end if

      do f=1, pflds
         fincl(f, 1) = fincl1(f)
         fincl(f, 2) = fincl2(f)
         fincl(f, 3) = fincl3(f)
         fincl(f, 4) = fincl4(f)
         fincl(f, 5) = fincl5(f)
         fincl(f, 6) = fincl6(f)
         
         fincllonlat(f, 1) = fincl1lonlat(f)
         fincllonlat(f, 2) = fincl2lonlat(f)
         fincllonlat(f, 3) = fincl3lonlat(f)
         fincllonlat(f, 4) = fincl4lonlat(f)
         fincllonlat(f, 5) = fincl5lonlat(f)
         fincllonlat(f, 6) = fincl6lonlat(f)
         if(dycore_is('UNSTRUCTURED') ) then
            do i=1,6
               if (fincllonlat(f,i) .ne. ' ') then
                  call endrun('READ_NAMELIST: Column output is not supported in Unstructered Grids')
               end if
            end do
         end if


         fexcl(f, 1) = fexcl1(f)
         fexcl(f, 2) = fexcl2(f)
         fexcl(f, 3) = fexcl3(f)
         fexcl(f, 4) = fexcl4(f)
         fexcl(f, 5) = fexcl5(f)
         fexcl(f, 6) = fexcl6(f)

         fhstpr(f, 1) = fhstpr1(f)
         fhstpr(f, 2) = fhstpr2(f)
         fhstpr(f, 3) = fhstpr3(f)
         fhstpr(f, 4) = fhstpr4(f)
         fhstpr(f, 5) = fhstpr5(f)
         fhstpr(f, 6) = fhstpr6(f)

         fwrtpr(f, 1) = fwrtpr1(f)
         fwrtpr(f, 2) = fwrtpr2(f)
         fwrtpr(f, 3) = fwrtpr3(f)
         fwrtpr(f, 4) = fwrtpr4(f)
         fwrtpr(f, 5) = fwrtpr5(f)
         fwrtpr(f, 6) = fwrtpr6(f)
      enddo
   end if
!
! Scatter namelist data to all processes
#if ( defined SPMD )
   call distnl ( )
#endif
!
! Auxiliary history files:
! Store input auxf values in array aux (from common block /comhst/).
!
! If generate an initial conditions history file as an auxillary tape:
!
   ctemp = shr_string_toUpper(inithist) 
   inithist = trim(ctemp)
   if (inithist /= '6-HOURLY' .and. inithist /= 'DAILY' .and. &
       inithist /= 'MONTHLY'  .and. inithist /= 'YEARLY' .and. &
       inithist /= 'CAMIOP'   .and. inithist /= 'ENDOFRUN') then
      inithist = 'NONE'
   endif
!
! Ensure that monthly averages have not been specified for aux. tapes
!
   do t=2,ptapes
      if (nhtfrq(t) == 0) then
         call endrun ('READ_NAMELIST: Only the primary history file may be monthly averaged')
      end if
   end do
! 
! History file write up times
! Convert write freq. of hist files from hours to timesteps if necessary.
! 
   do t=1,ptapes
      if (nhtfrq(t) < 0) then
         nhtfrq(t) = nint((-nhtfrq(t)*3600._r8)/dtime)
      end if
   end do
!
! Initialize the filename specifier if not already set
! This is the format for the history filenames:
! %c= caseid, %t=tape no., %y=year, %m=month, %d=day, %s=second, %%=%
! See the filenames module for more information
!
   do t = 1, ptapes
      if ( len_trim(hfilename_spec(t)) == 0 )then
         if ( nhtfrq(t) == 0 )then
            hfilename_spec(t) = '%c.cam2.h%t.%y-%m.nc'        ! Monthly files
         else
            hfilename_spec(t) = '%c.cam2.h%t.%y-%m-%d-%s.nc'
         end if
      end if
   end do
!
! Only one time sample allowed per monthly average file
! 
   if (nhtfrq(1) == 0) mfilt(1) = 1

   ! Print per-tape averaging flags
   if (masterproc) then
      do t=1,ptapes
         if (avgflag_pertape(t) /= ' ') then
            write(iulog,*)'Unless overridden by namelist input on a per-field basis (FINCL),'
            write(iulog,*)'All fields on history file ',t,' will have averaging flag ',avgflag_pertape(t)
         end if
      end do
   end if

#if ( defined SPMD )
! 
! Set runtime options for spmd_dyn
!
   if ( dycore_is ('HOMME') ) then
	call spmd_dyn_setopts()
   else if ( dycore_is ('LR') ) then
      call spmd_dyn_setopts(                          &
             npr_yz_in            =npr_yz,            &
             geopktrans_in        =geopktrans,        &
             geopkblocks_in       =geopkblocks,       &
             force_2d_in          =force_2d,          &
             modcomm_transpose_in =modcomm_transpose, &
             modcomm_geopk_in     =modcomm_geopk,     &
             modcomm_gatscat_in   =modcomm_gatscat,   &
             modc_sw_dynrun_in    =modc_sw_dynrun,    &
             modc_hs_dynrun_in    =modc_hs_dynrun,    &
             modc_send_dynrun_in  =modc_send_dynrun,  &
             modc_mxreq_dynrun_in =modc_mxreq_dynrun, &
             modc_sw_cdcore_in    =modc_sw_cdcore,    &
             modc_hs_cdcore_in    =modc_hs_cdcore,    &
             modc_send_cdcore_in  =modc_send_cdcore,  &
             modc_mxreq_cdcore_in =modc_mxreq_cdcore, &
             modc_sw_gather_in    =modc_sw_gather,    &
             modc_hs_gather_in    =modc_hs_gather,    &
             modc_send_gather_in  =modc_send_gather,  &
             modc_mxreq_gather_in =modc_mxreq_gather, &
             modc_sw_scatter_in   =modc_sw_scatter,   &
             modc_hs_scatter_in   =modc_hs_scatter,   &
             modc_send_scatter_in =modc_send_scatter, &
             modc_mxreq_scatter_in=modc_mxreq_scatter,&
             modc_sw_tracer_in    =modc_sw_tracer,    &
             modc_hs_tracer_in    =modc_hs_tracer,    &
             modc_send_tracer_in  =modc_send_tracer,  &
             modc_mxreq_tracer_in =modc_mxreq_tracer, &
             modc_onetwo_in       =modc_onetwo,       &
             modc_tracers_in      =modc_tracers       )
   else if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call spmd_dyn_setopts(                     &
             dyn_alltoall_in    =dyn_alltoall,   &
             dyn_allgather_in   =dyn_allgather,  &
             dyn_equi_by_col_in =dyn_equi_by_col,&
             dyn_npes_in        =dyn_npes,       &
             dyn_npes_stride_in =dyn_npes_stride )
   endif
! 
! Set runtime options for swap communications.
!
   call swap_comm_setopts(                          &
          swap_comm_protocol_in=swap_comm_protocol, &
          swap_comm_maxreq_in=swap_comm_maxreq)
! 
! Set runtime options for flow-controlled gather collective.
!
   call fc_gather_setopts(                           &
          fc_gather_flow_cntl_in=fc_gather_flow_cntl  )
#endif
   ! restart write interval
   call restart_setopts( nsrest,            &
      cam_branch_file_in          =cam_branch_file            )


   ! Set runtime options for physics chunking.
   call phys_grid_setopts(                          &
       phys_loadbalance_in    =phys_loadbalance,    &
       phys_twin_algorithm_in =phys_twin_algorithm, &
       phys_alltoall_in       =phys_alltoall,       &
       phys_chnk_per_thd_in   =phys_chnk_per_thd    )

   ! Set runtime options for reproducible sum calculations.
   call repro_sum_setopts(                              &
      repro_sum_use_ddpdd_in=repro_sum_use_ddpdd,       &
      repro_sum_rel_diff_max_in=repro_sum_rel_diff_max, &
      repro_sum_recompute_in=repro_sum_recompute,       &
      repro_sum_master=masterproc,                      &
      repro_sum_logunit=iulog                           )

   ! physics buffer
   call pbuf_setopts( &
      pbuf_global_allocate_in = pbuf_global_allocate )

   ! diagnostics
   call diag_setopts( &
      diag_cnst_conv_tend_in = diag_cnst_conv_tend )

   ! conservation
   call check_energy_setopts( &
      print_energy_errors_in = print_energy_errors )

   call radiation_setopts( dtime, nhtfrq(1), &
      iradsw_in      = iradsw,     &
      iradlw_in      = iradlw,     &
      iradae_in      = iradae,     &
      irad_always_in = irad_always )

#if (defined WACCM_PHYS)
   ! iondrag / efield
   call iondrag_setopts( &
        efield_lflux_file_in =efield_lflux_file, &
        efield_hflux_file_in =efield_hflux_file, &
        efield_wei96_file_in =efield_wei96_file)
   ! qbo forcing
   call qbo_setopts( &
        qbo_use_forcing_in  = qbo_use_forcing, &
        qbo_forcing_file_in = qbo_forcing_file,&
        qbo_cyclic_in       = qbo_cyclic       )
#endif

   ! Upper atmosphere radiative processes
   call radheat_setopts( &
      nlte_use_mo_in =nlte_use_mo, &
      itgcmcyc_in    =itgcmcyc,    &
      cftgcm_in      =cftgcm       )

! 
! Set runtime options for comhd
!
   call comhd_setopts( dif2_in  =dif2, &
                       dif4_in  =dif4, &
                       kmxhdc_in=kmxhdc)
! 
! Set runtime options for single column mode
!
   if (present(single_column_in) .and. present(scmlon_in) .and. present(scmlat_in)) then 
      if (single_column_in) then
         single_column = single_column_in
         scmlon = scmlon_in
         scmlat = scmlat_in
         call scam_setopts( scmlat_in=scmlat,scmlon_in=scmlon, &
                            iopfile_in=iopfile,single_column_in=single_column,&
                            scm_iop_srf_prop_in=scm_iop_srf_prop,&
                            scm_relaxation_in=scm_relaxation, &
                            scm_diurnal_avg_in=scm_diurnal_avg, &
                            scm_crm_mode_in=scm_crm_mode)
      end if
   endif
#if ( defined OFFLINE_DYN )
! 
! Set runtime options for comhd
!
   call offline_met_setopts( met_data_file_in = met_data_file, &
                             met_remove_file_in = met_remove_file, &
                             met_cell_wall_winds_in = met_cell_wall_winds, &
                             met_filenames_list_in = met_filenames_list, &
                             met_rlx_top_km_in = met_rlx_top, &
                             met_rlx_bot_km_in = met_rlx_bot, &
                             met_max_rlx_in = met_max_rlx, &
                             met_fix_mass_in = met_fix_mass, &
                             met_shflx_name_in = met_shflx_name, &
                             met_qflx_name_in = met_qflx_name, &
                             met_shflx_factor_in = met_shflx_factor, &
                             met_qflx_factor_in = met_qflx_factor )
#endif

   ! Call subroutines for modules to read their own namelist.
   ! In some cases namelist default values may depend on settings from
   ! other modules, so there may be an order dependence in the following
   ! calls.
   ! ***N.B.*** In particular, physconst_readnl should be called before
   !            the other readnl methods in case that method is used to set
   !            physical constants, some of which are set at runtime
   !            by the physconst_readnl method.
   ! Modules that read their own namelist are responsible for making sure
   ! all processes receive the values.

   call physconst_readnl(nlfilename)
   call chem_surfvals_readnl(nlfilename)
   call phys_ctl_readnl(nlfilename)
   call cam3_aero_data_readnl(nlfilename)
   call cam3_ozone_data_readnl(nlfilename)
   call cldfrc_readnl(nlfilename)
   call zmconv_readnl(nlfilename)
   call cldwat_readnl(nlfilename)
   call hkconv_readnl(nlfilename)
   call cld_sediment_readnl(nlfilename)
   call gw_drag_readnl(nlfilename)
   call phys_debug_readnl(nlfilename)
   call rad_cnst_readnl(nlfilename)
   call rad_data_readnl(nlfilename)
   call chem_readnl(nlfilename)
   call prescribed_volcaero_readnl(nlfilename)
   call solar_data_readnl(nlfilename)
   call tropopause_readnl(nlfilename)
   call aoa_tracers_readnl(nlfilename)
   call aerodep_flx_readnl(nlfilename)
   call prescribed_ozone_readnl(nlfilename)
   call prescribed_aero_readnl(nlfilename)
   call prescribed_ghg_readnl(nlfilename)
   call co2_cycle_readnl(nlfilename)
   call aircraft_emit_readnl(nlfilename)
   call cospsimulator_intr_readnl(nlfilename)

! 
! Print cam_inparm input variables to standard output
! 
   if (masterproc) then
      write(iulog,*)' ------------------------------------------'
      write(iulog,*)'     *** INPUT VARIABLES (CAM_INPARM) ***'
      write(iulog,*)' ------------------------------------------'
      if (nsrest/=0) then
         write(iulog,*) '  Continuation of an earlier run'
      else
         write(iulog,*) '         Initial run'
      end if
      write(iulog,*) ' ********** CASE = ',trim(caseid),' **********'
      write(iulog,'(1x,a)') ctitle
      if (len_trim(ncdata) > 0) then
         write(iulog,*) 'Initial dataset is: ',trim(ncdata)
      end if
      write(iulog,*)'Topography dataset is: ', trim(bnd_topo)
      write(iulog,*)'Time-invariant (absorption/emissivity) factor dataset is: ', trim(absems_data)
#ifdef MODAL_AERO
      write(iulog,*)'Time-invariant modal aerosol optics datset is: ', trim(modal_optics)
#endif

      ! Type of run
      write(iulog,*)'Run type flag (NSREST) 0=initial, 1=restart, 3=branch ',nsrest

      call restart_printopts()

   end if
!
! History file info 
!
   if (masterproc) then
      if (inithist == '6-HOURLY' ) then
         write(iulog,*)'Initial conditions history files will be written 6-hourly.'
      else if (inithist == 'DAILY' ) then
         write(iulog,*)'Initial conditions history files will be written daily.'
      else if (inithist == 'MONTHLY' ) then
         write(iulog,*)'Initial conditions history files will be written monthly.'
      else if (inithist == 'YEARLY' ) then
         write(iulog,*)'Initial conditions history files will be written yearly.'
      else if (inithist == 'CAMIOP' ) then
         write(iulog,*)'Initial conditions history files will be written for IOP.'
      else if (inithist == 'ENDOFRUN' ) then
         write(iulog,*)'Initial conditions history files will be written at end of run.'
      else
         write(iulog,*)'Initial conditions history files will not be created'
      end if
!
! Write physics variables from namelist cam_inparm to std. output
!
      write(iulog,9108) eps,dif2,dif4,kmxhdc,nlvdry
9108 format(' Time filter coefficient (EPS)                 ',f10.3,/,&
            ' DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3/, &
            ' DEL4 Horizontal diffusion coefficient (DIF4)  ',e10.3/, &
            ' Number of levels Courant limiter applied      ',i10/,   &
            ' Lowest level for dry adiabatic adjust (NLVDRY)',i10)

      call radiation_printopts()

      if (dycore_is ('EUL') .or. dycore_is ('SLD')) then
         if (divdampn > 0._r8) then
            write(iulog,*) 'Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
         elseif (divdampn < 0._r8) then
            call endrun ('READ_NAMELIST: divdampn must be a positive number')
         else
            write(iulog,*) 'Divergence damper for spectral dycore NOT invoked'
         endif
      end if

      if ( (adiabatic .and. ideal_phys) .or. (adiabatic .and. aqua_planet) .or. &
           (ideal_phys .and. aqua_planet) ) then
         call endrun ('READ_NAMELIST: Only one of ADIABATIC, IDEAL_PHYS, or AQUA_PLANET can be .true.')
      end if

#ifdef COUP_SOM
      if (adiabatic .or. ideal_phys .or. aqua_planet )then
         call endrun ('READ_NAMELIST: adiabatic, ideal_phys or aqua_planet can not be used with SOM')
      end if
#else
      if (adiabatic)   write(iulog,*) 'Model will run ADIABATICALLY (i.e. no physics)'
      if (ideal_phys)  write(iulog,*) 'Run ONLY the "idealized" dynamical core of the ', &
                                  'model  (dynamics + Held&Suarez-specified physics)'
      if (aqua_planet) write(iulog,*) 'Run model in "AQUA_PLANET" mode'
#endif
   end if

   ! set public data in cam_control_mod
   moist_physics = (.not. adiabatic) .and. (.not. ideal_phys)

#ifdef PERGRO
   if (masterproc) then
      write(iulog,*)'pergro for cloud water is true'
   end if
#endif

   ntspdy = nint(86400._r8/dtime) ! no. timesteps per day

   if (masterproc) then
      if (doisccp) then
         write(iulog,*)'ISCCP calcs (version 3.4) and history IO will be done'
      else if (doisccp_38) then
         write(iulog,*)'ISCCP calcs (version 3.8) and history IO will be done'
      else
         write(iulog,*)'ISCCP calcs and history IO will NOT be done'
      end if
   end if

end subroutine read_namelist


!=======================================================================

#ifdef SPMD
subroutine distnl
!-----------------------------------------------------------------------
!     
! Purpose:     
! Distribute namelist data to all processors.
!
! The cpp SPMD definition provides for the funnelling of all program i/o
! through the master processor. Processor 0 either reads restart/history
! data from the disk and distributes it to all processors, or collects
! data from all processors and writes it to disk.
!     
!---------------------------Code history-------------------------------
!
! Original version:  CCM2
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
!
!-----------------------------------------------------------------------
   use mpishorthand
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
! 
   call mpibcast (dtime,       1,mpiint,0,mpicom)

   call mpibcast (nhstpr  ,ptapes,mpiint,0,mpicom)
   call mpibcast (ndens   ,ptapes,mpiint,0,mpicom)
   call mpibcast (nhtfrq  ,ptapes,mpiint,0,mpicom)
   call mpibcast (mfilt   ,ptapes,mpiint,0,mpicom)
   call mpibcast (nsrest  ,1,mpiint,0,mpicom)
   call mpibcast (kmxhdc  ,1,mpiint,0,mpicom)
   call mpibcast (nlvdry  ,1,mpiint,0,mpicom)

   ! f-v dynamics specific
   call mpibcast (nsplit  ,1,mpiint,0,mpicom)
   call mpibcast (nspltrac,1,mpiint,0,mpicom)
   call mpibcast (nspltvrm,1,mpiint,0,mpicom)
   call mpibcast (iord    ,1,mpiint,0,mpicom)
   call mpibcast (jord    ,1,mpiint,0,mpicom)
   call mpibcast (kord    ,1,mpiint,0,mpicom)
   call mpibcast (dyn_conservative,1,mpilog,0,mpicom)
   call mpibcast (filtcw  ,1,mpiint,0,mpicom)
   call mpibcast (ct_overlap  ,1,mpiint,0,mpicom)
   call mpibcast (trac_decomp ,1,mpiint,0,mpicom)
   call mpibcast (fft_flt ,1,mpiint,0,mpicom)
   call mpibcast (div24del2flag ,1,mpiint,0,mpicom)
   call mpibcast (del2coef ,1,mpir8,0,mpicom)

   call mpibcast (rayk0    ,1,mpiint,0,mpicom)
   call mpibcast (raykrange,1,mpir8,0,mpicom)
   call mpibcast (raytau0  ,1,mpir8,0,mpicom)

   call mpibcast (divdampn,1,mpir8,0,mpicom)
   call mpibcast (eps     ,1,mpir8,0,mpicom)
   call mpibcast (dif2    ,1,mpir8,0,mpicom)
   call mpibcast (dif4    ,1,mpir8,0,mpicom)

   call mpibcast (tracers_flag,1,mpilog,0,mpicom)
   call mpibcast (readtrace   ,1,mpilog,0,mpicom)
   call mpibcast (adiabatic   ,1,mpilog,0,mpicom)
   call mpibcast (ideal_phys  ,1,mpilog,0,mpicom)
   call mpibcast (aqua_planet ,1,mpilog,0,mpicom)

   call mpibcast (empty_htapes,1,mpilog,0,mpicom)
   call mpibcast (use_64bit_nc,1,mpilog,0,mpicom)
   call mpibcast (print_step_cost,1,mpilog,0,mpicom)
   call mpibcast (inithist_all   ,1,mpilog,0,mpicom)
   call mpibcast (doisccp     ,1,mpilog,0,mpicom)
   call mpibcast (doisccp_38  ,1,mpilog,0,mpicom)
   call mpibcast (pertlim     ,1, mpir8,  0, mpicom )

   call mpibcast (caseid  ,len(caseid) ,mpichar,0,mpicom)
   call mpibcast (avgflag_pertape, ptapes, mpichar,0,mpicom)
   call mpibcast (ctitle  ,len(ctitle),mpichar,0,mpicom)
   call mpibcast (ncdata  ,len(ncdata) ,mpichar,0,mpicom)
   call mpibcast (bnd_topo  ,len(bnd_topo) ,mpichar,0,mpicom)
   call mpibcast (absems_data,len(absems_data),mpichar,0,mpicom)
#ifdef MODAL_AERO
   call mpibcast (modal_optics,len(modal_optics),mpichar,0,mpicom)
#endif
   call mpibcast (cam_branch_file  ,len(cam_branch_file) ,mpichar,0,mpicom)
   call mpibcast (inithist,len(inithist)  ,mpichar,0,mpicom)
   call mpibcast (hfilename_spec, len(hfilename_spec(1))*ptapes, mpichar, 0, mpicom)
   call mpibcast (fincl   ,len(fincl (1,1))*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fexcl   ,len(fexcl (1,1))*pflds*ptapes,mpichar,0,mpicom)

   call mpibcast (fincllonlat   ,len(fincllonlat (1,1))*pflds*ptapes,mpichar,0,mpicom)

   call mpibcast (fhstpr  ,len(fhstpr(1,1))*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fwrtpr  ,len(fwrtpr(1,1))*pflds*ptapes,mpichar,0,mpicom)
#if (defined OFFLINE_DYN)
!
! Offline dynamics parameters
!
   call mpibcast (met_data_file  ,len(met_data_file) ,mpichar,0,mpicom)
   call mpibcast (met_remove_file    ,1 ,mpilog, 0, mpicom )
   call mpibcast (met_cell_wall_winds,1 ,mpilog, 0, mpicom )
   call mpibcast (met_filenames_list ,len(met_filenames_list),mpichar,0,mpicom)
   call mpibcast (met_rlx_top,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_rlx_bot,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_max_rlx,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_fix_mass,       1 ,mpilog, 0, mpicom )
   call mpibcast (met_qflx_name      ,len(met_qflx_name),     mpichar,0,mpicom)
   call mpibcast (met_shflx_name     ,len(met_shflx_name),    mpichar,0,mpicom)
   call mpibcast (met_qflx_factor    ,1, mpir8,  0, mpicom )
   call mpibcast (met_shflx_factor   ,1, mpir8,  0, mpicom )
#endif
!
! Orbital stuff
!
   call mpibcast (eccen   ,1  ,mpir8 ,0,mpicom)
   call mpibcast (obliqr  ,1  ,mpir8 ,0,mpicom)
   call mpibcast (mvelpp  ,1  ,mpir8 ,0,mpicom)
   call mpibcast (lambm0  ,1  ,mpir8 ,0,mpicom)
   call mpibcast (indirect     , 1 ,mpilog, 0,mpicom)

   ! spmd_dyn
   if ( dycore_is ('LR') ) then
      call mpibcast (npr_yz            ,4,mpiint,0,mpicom)
      call mpibcast (geopktrans        ,1,mpiint,0,mpicom)
      call mpibcast (geopkblocks       ,1,mpiint,0,mpicom)
      call mpibcast (force_2d          ,1,mpiint,0,mpicom)
      call mpibcast (modcomm_transpose ,1,mpiint,0,mpicom)
      call mpibcast (modcomm_geopk     ,1,mpiint,0,mpicom)
      call mpibcast (modcomm_gatscat   ,1,mpiint,0,mpicom)
      call mpibcast (modc_sw_dynrun    ,1,mpiint,0,mpicom)
      call mpibcast (modc_hs_dynrun    ,1,mpilog,0,mpicom)
      call mpibcast (modc_send_dynrun  ,1,mpilog,0,mpicom)
      call mpibcast (modc_mxreq_dynrun ,1,mpiint,0,mpicom)
      call mpibcast (modc_sw_cdcore    ,1,mpiint,0,mpicom)
      call mpibcast (modc_hs_cdcore    ,1,mpilog,0,mpicom)
      call mpibcast (modc_send_cdcore  ,1,mpilog,0,mpicom)
      call mpibcast (modc_mxreq_cdcore ,1,mpiint,0,mpicom)
      call mpibcast (modc_sw_gather    ,1,mpiint,0,mpicom)
      call mpibcast (modc_hs_gather    ,1,mpilog,0,mpicom)
      call mpibcast (modc_send_gather  ,1,mpilog,0,mpicom)
      call mpibcast (modc_mxreq_gather ,1,mpiint,0,mpicom)
      call mpibcast (modc_sw_scatter   ,1,mpiint,0,mpicom)
      call mpibcast (modc_hs_scatter   ,1,mpilog,0,mpicom)
      call mpibcast (modc_send_scatter ,1,mpilog,0,mpicom)
      call mpibcast (modc_mxreq_scatter,1,mpiint,0,mpicom)
      call mpibcast (modc_sw_tracer    ,1,mpiint,0,mpicom)
      call mpibcast (modc_hs_tracer    ,1,mpilog,0,mpicom)
      call mpibcast (modc_send_tracer  ,1,mpilog,0,mpicom)
      call mpibcast (modc_mxreq_tracer ,1,mpiint,0,mpicom)
      call mpibcast (modc_onetwo       ,1,mpiint,0,mpicom)
      call mpibcast (modc_tracers      ,1,mpiint,0,mpicom)
   endif
   if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call mpibcast (dyn_alltoall   ,1,mpiint,0,mpicom)
      call mpibcast (dyn_allgather  ,1,mpiint,0,mpicom)
      call mpibcast (dyn_equi_by_col,1,mpilog,0,mpicom)
      call mpibcast (dyn_npes       ,1,mpiint,0,mpicom)
      call mpibcast (dyn_npes_stride,1,mpiint,0,mpicom)
   endif

   ! Physics chunk tuning
   call mpibcast (phys_loadbalance   ,1,mpiint,0,mpicom)
   call mpibcast (phys_twin_algorithm,1,mpiint,0,mpicom)
   call mpibcast (phys_alltoall      ,1,mpiint,0,mpicom)
   call mpibcast (phys_chnk_per_thd  ,1,mpiint,0,mpicom)

   !  Reproducible sum options
   call mpibcast (repro_sum_use_ddpdd,1,mpilog,0,mpicom)
   call mpibcast (repro_sum_rel_diff_max,1,mpir8,0,mpicom)
   call mpibcast (repro_sum_recompute,1,mpilog,0,mpicom)

   ! Interprocessor communication tuning
   call mpibcast (swap_comm_protocol,1,mpiint,0,mpicom)
   call mpibcast (swap_comm_maxreq,1,mpiint,0,mpicom)
   call mpibcast (fc_gather_flow_cntl,1,mpiint,0,mpicom)

   ! Physics buffer
   call mpibcast (pbuf_global_allocate, 1, mpilog, 0, mpicom)

   ! Diagnostic options
   call mpibcast (diag_cnst_conv_tend, len(diag_cnst_conv_tend), mpichar, 0, mpicom)

   ! Conservation
   call mpibcast (print_energy_errors, 1, mpilog, 0, mpicom)

   ! Radiative heating calculation
   call mpibcast (iradsw,     1, mpiint, 0, mpicom)
   call mpibcast (iradlw,     1, mpiint, 0, mpicom)
   call mpibcast (iradae,     1, mpiint, 0, mpicom)
   call mpibcast (irad_always,1, mpiint, 0, mpicom)

#if (defined WACCM_PHYS)
   ! iondrag / efield options
   call mpibcast (efield_lflux_file, len(efield_lflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_hflux_file, len(efield_hflux_file), mpichar, 0, mpicom)
   call mpibcast (efield_wei96_file, len(efield_wei96_file), mpichar, 0, mpicom)
   ! qbo variables
   call mpibcast (qbo_forcing_file,  len(qbo_forcing_file ), mpichar, 0, mpicom)
   call mpibcast (qbo_use_forcing,   1,                      mpilog,  0, mpicom)
   call mpibcast (qbo_cyclic,        1,                      mpilog,  0, mpicom)
#endif

end subroutine distnl
#endif



subroutine preset
!----------------------------------------------------------------------- 
! 
! Purpose: Preset namelist CAM_INPARM input variables and initialize some other variables
! 
! Method: Hardwire the values
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use cam_history,  only: fincl, fexcl, fhstpr, fwrtpr,fincllonlat
   use rgrid
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Preset character history variables here because module initialization of character arrays
! does not work on all machines
! $$$ TBH:  is this still true?  12/14/03
!
   fincl(:,:)  = ' '
   fincllonlat(:,:)  = ' '
   fexcl(:,:)  = ' '
   fhstpr(:,:) = ' '
   fwrtpr(:,:) = ' '
!
! Flags
!
   print_step_cost = .false.   ! print per timestep cost info
!
! Numerical scheme default values
!
   eps    = 0.06_r8
   nlvdry = 3
!
! No divergence damping
!
   divdampn = 0._r8
!
! rgrid: set default to full grid
!
   nlon(:) = plon
!!
!! Unit numbers: set to invalid
!!
!   ncid_ini = -1
!   ncid_sst = -1
!   ncid_trc = -1
!
! /perturb/
!
  pertlim = 0.0_r8

   return
end subroutine preset

end module runtime_opts
