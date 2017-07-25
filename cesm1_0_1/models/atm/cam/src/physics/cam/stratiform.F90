
  #undef DEBUG

  module stratiform

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macro and microphysics
  !
  ! Author: Byron Boville  Sept 04, 2002
  !         modified by D.B. Coleman May 2004
  !         modified by Sungsu Park. Dec.2009
  !         modified by J. Kay Jan. 2010 to add COSP simulator info from RK microphysics to physics buffer
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice
  use abortutils,    only: endrun
  use chemistry,     only: chem_is
  use phys_control,  only: phys_getopts

  use perf_mod
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: stratiform_register, stratiform_init_cnst, stratiform_implements_cnst
  public :: stratiform_init, conv_water_4rad
  public :: stratiform_tend

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! Choose either 'intermediate' ('inter') or complete ('compl') cloud microphysics 
  ! inter : Microphysics assumes 'liquid stratus frac = ice stratus frac = max( liquid stratus frac, ice stratus frac )'.
  ! compl : Microphysics explicitly treats 'liquid stratus frac .ne. ice stratus frac'  
  ! for CAM5, only 'inter' is functional

    character(len=5), private, parameter :: micro_treatment = 'inter' 

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus 
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus, 
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.

    logical,          private, parameter :: cu_det_st  = .false.  

  ! -------------------------------- !
  ! End of Private Module Parameters !
  ! -------------------------------- !

  integer, parameter :: ncnstmax = 4                    ! Number of constituents
  integer            :: ncnst 	  		        ! Number of constituents (can vary)
  character(len=8), dimension(ncnstmax), parameter &    ! Constituent names
                     :: cnst_names = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)
  logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc
  character(len=16)  :: microp_scheme                   ! Microphysics scheme

  integer :: &
     cldo_idx,     &! old cld index in physics buffer
     kvh_idx,      &! kvh index in physics buffer
#ifdef MODAL_AERO
     qqcw_idx,     &! qqcw index in physics buffer
#endif
     ixcldliq,     &! cloud liquid amount index
     ixcldice,     &! cloud ice amount index
     ixnumliq,     &! cloud liquid number index
     ixnumice,     &! cloud ice water index
     qcwat_idx,    &! qcwat index in physics buffer
     lcwat_idx,    &! lcwat index in physics buffer
     iccwat_idx,   &! iccwat index in physics buffer
     nlwat_idx,    &! nlwat index in physics buffer
     niwat_idx,    &! niwat index in physics buffer
     tcwat_idx,    &! tcwat index in physics buffer
     CC_T_idx,     &!
     CC_qv_idx,    &!
     CC_ql_idx,    &!
     CC_qi_idx,    &!
     CC_nl_idx,    &!
     CC_ni_idx,    &!
     CC_qlst_idx,  &!
     cld_idx,      &! cld index in physics buffer
     ast_idx,      &! stratiform cloud fraction index in physics buffer
     aist_idx,     &! ice stratiform cloud fraction index in physics buffer
     alst_idx,     &! liquid stratiform cloud fraction index in physics buffer
     qist_idx,     &! ice stratiform in-cloud IWC 
     qlst_idx,     &! liquid stratiform in-cloud LWC  
     concld_idx,   &! concld index in physics buffer
     rhdfda_idx,   &! rhdfda index in physics buffer
     rhu00_idx,    &! rhu00 index in physics buffer
     rel2_idx,     &! rel2 index in physics buffer
     rei2_idx,     &! rei2 index in physics buffer
     ls_flxprc_idx,&
     ls_flxsnw_idx


  contains

  ! ===============================================================================

  subroutine stratiform_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add

    integer idx

  !-----------------------------------------------------------------------

    call phys_getopts(microp_scheme_out=microp_scheme)

    if ( microp_scheme .eq. 'MG' ) then
	ncnst = 4
    else if ( microp_scheme .eq. 'RK' ) then
	ncnst = 2
    end if

  ! Register cloud water and determine index.

    call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount')
    call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount')

    if ( microp_scheme .eq. 'MG' ) then
       call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, ixnumliq, &
          longname='Grid box averaged cloud liquid number')
       call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, ixnumice, &
          longname='Grid box averaged cloud ice number')
    end if

  ! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add('QCWAT',   'global',  1, pver, pbuf_times,   qcwat_idx)
    call pbuf_add('LCWAT',   'global',  1, pver, pbuf_times,   lcwat_idx)
    call pbuf_add('ICCWAT',  'global',  1, pver, pbuf_times,  iccwat_idx)
    call pbuf_add('NLWAT',   'global',  1, pver, pbuf_times,   nlwat_idx)
    call pbuf_add('NIWAT',   'global',  1, pver, pbuf_times,   niwat_idx)
    call pbuf_add('CC_T',    'global',  1, pver, pbuf_times,    CC_T_idx)
    call pbuf_add('CC_qv',   'global',  1, pver, pbuf_times,   CC_qv_idx)
    call pbuf_add('CC_ql',   'global',  1, pver, pbuf_times,   CC_ql_idx)
    call pbuf_add('CC_qi',   'global',  1, pver, pbuf_times,   CC_qi_idx)
    call pbuf_add('CC_nl',   'global',  1, pver, pbuf_times,   CC_nl_idx)
    call pbuf_add('CC_ni',   'global',  1, pver, pbuf_times,   CC_ni_idx)
    call pbuf_add('CC_qlst', 'global',  1, pver, pbuf_times, CC_qlst_idx)
    call pbuf_add('TCWAT',   'global',  1, pver, pbuf_times,   tcwat_idx)
    call pbuf_add('CLD',     'global',  1, pver, pbuf_times,     cld_idx)
    call pbuf_add('CLDO',    'global',  1, pver, pbuf_times,    cldo_idx)
    call pbuf_add('AST',     'global',  1, pver, pbuf_times,     ast_idx)
    call pbuf_add('AIST',    'global',  1, pver, pbuf_times,    aist_idx)
    call pbuf_add('ALST',    'global',  1, pver, pbuf_times,    alst_idx)
    call pbuf_add('QIST',    'global',  1, pver, pbuf_times,    qist_idx)
    call pbuf_add('QLST',    'global',  1, pver, pbuf_times,    qlst_idx)
    call pbuf_add('CONCLD',  'global',  1, pver, pbuf_times,  concld_idx)
    call pbuf_add('RHDFDA',  'global',  1, pver, pbuf_times,  rhdfda_idx)
    call pbuf_add('RHU00',   'global',  1, pver, pbuf_times,   rhu00_idx)
    call pbuf_add('REL2',    'global',  1, pver, pbuf_times,    rel2_idx)
    call pbuf_add('REI2',    'global',  1, pver, pbuf_times,    rei2_idx)

  ! Physics buffer variables for convective cloud properties.

    call pbuf_add('CONCLDQL',   'physpkg', 1, pver, 1, idx) 
    call pbuf_add('FICE',       'physpkg', 1, pver, 1, idx) 
    call pbuf_add('SH_FRAC',    'physpkg', 1, pver, 1, idx) 
    call pbuf_add('DP_FRAC',    'physpkg', 1, pver, 1, idx) 

    call pbuf_add('QINI'      , 'physpkg', 1,pver, 1, idx)
    call pbuf_add('CLDLIQINI' , 'physpkg', 1,pver, 1, idx)
    call pbuf_add('CLDICEINI' , 'physpkg', 1,pver, 1, idx)
    call pbuf_add('TINI'      , 'physpkg', 1,pver, 1, idx)

    call pbuf_add('QME',        'physpkg', 1, pver, 1, idx)
    call pbuf_add('PRAIN' ,     'physpkg', 1, pver, 1, idx)
    call pbuf_add('NEVAPR' ,    'physpkg', 1, pver, 1, idx)

    call pbuf_add('WSEDL',      'physpkg', 1, pver, 1, idx)

    call pbuf_add('REI',        'physpkg', 1, pver, 1, idx)
    call pbuf_add('REL',        'physpkg', 1, pver, 1, idx)
    call pbuf_add('REL_FN',     'physpkg', 1, pver, 1, idx)          ! REL at fixed number for indirect rad forcing

    call pbuf_add('DEI',        'physpkg', 1, pver, 1, idx)          ! Mitchell ice effective diameter for radiation
    call pbuf_add('MU',         'physpkg', 1, pver, 1, idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('LAMBDAC',    'physpkg', 1, pver, 1, idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('ICIWP',      'physpkg', 1, pver, 1, idx)          ! In cloud ice water path for radiation
    call pbuf_add('ICLWP',      'physpkg', 1, pver, 1, idx)          ! In cloud liquid water path for radiation

    call pbuf_add('DEICONV',    'physpkg', 1, pver, 1, idx)          ! Convective ice effective diameter for radiation
    call pbuf_add('MUCONV',     'physpkg', 1, pver, 1, idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('LAMBDACONV', 'physpkg', 1, pver, 1, idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('ICIWPST',    'physpkg', 1, pver, 1, idx)          ! Stratiform only in cloud ice water path for radiation
    call pbuf_add('ICLWPST',    'physpkg', 1, pver, 1, idx)          ! Stratiform in cloud liquid water path for radiation
    call pbuf_add('ICIWPCONV',  'physpkg', 1, pver, 1, idx)          ! Convective only in cloud ice water path for radiation
    call pbuf_add('ICLWPCONV',  'physpkg', 1, pver, 1, idx)          ! Convective in cloud liquid water path for radiation

    call pbuf_add('DES',        'physpkg', 1, pver, 1, idx)          ! Snow effective diameter for radiation
    call pbuf_add('ICSWP',      'physpkg', 1, pver, 1, idx)          ! In cloud snow water path for radiation
    call pbuf_add('CLDFSNOW',   'physpkg', 1, pver ,pbuf_times, idx) ! Cloud fraction for liquid drops + snow

#ifdef MODAL_AERO
    call pbuf_add('QQCW',       'global',  1, pver, pcnst, qqcw_idx)
    call pbuf_add('RATE1_CW2PR_ST', 'physpkg', 1, pver, 1, idx)   ! rce 2010/05/01
#endif

    call pbuf_add('LS_FLXPRC',  'physpkg', 1, pverp, 1, ls_flxprc_idx)
    call pbuf_add('LS_FLXSNW',  'physpkg', 1, pverp, 1, ls_flxsnw_idx)

  end subroutine stratiform_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  function stratiform_implements_cnst(name)

  !----------------------------------------------------------------------------- ! 
  !                                                                              !    
  ! Purpose: return true if specified constituent is implemented by this package !
  !                                                                              !
  ! Author: B. Eaton                                                             !
  !                                                                              ! 
  !----------------------------------------------------------------------------- !
     implicit none
  !-----------------------------Arguments---------------------------------
     character(len=*), intent(in) :: name      ! constituent name
     logical :: stratiform_implements_cnst     ! return value
  !---------------------------Local workspace-----------------------------
     integer :: m
  !-----------------------------------------------------------------------

     stratiform_implements_cnst = .false.

     do m = 1, ncnst
        if (name == cnst_names(m)) then
           stratiform_implements_cnst = .true.
           return
        end if
     end do
  end function stratiform_implements_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_init_cnst(name, q, gcid)

  !----------------------------------------------------------------------- !
  !                                                                        !
  ! Initialize the cloud water mixing ratios (liquid and ice), if they are !
  ! not read from the initial file                                         ! 
  !                                                                        !
  !----------------------------------------------------------------------- !
    implicit none
  !---------------------------- Arguments ---------------------------------
    character(len=*), intent(in)  :: name     ! constituent name
    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id
  !-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMICE' ) then
       q = 0.0_r8
       return
    end if

  end subroutine stratiform_init_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_init

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !

    use cldwat,          only: inimc
    use microp_aero,     only: ini_microp_aero
    use cldwat2m_micro,  only: ini_micro
    use cldwat2m_macro,  only: ini_macro
    use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst
    use cam_history,     only: addfld, add_default, phys_decomp
    use physconst,       only: tmelt, rh2o, rhodair
    use convect_shallow, only: convect_shallow_use_shfrc
    use dycore,          only: dycore_is
    use phys_control,    only: cam_physpkg_is
#ifdef MODAL_AERO
    use ndrop,           only: activate_init
    use cam_history,     only: fieldname_len
    use spmd_utils,      only: masterproc
    use modal_aero_data, only: cnst_name_cw, &
                               lmassptr_amode, lmassptrcw_amode, &
                               nspec_amode, ntot_amode, numptr_amode, numptrcw_amode, maxd_amode
#endif

    integer              :: m, mm
    logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_microphysics ! Output variables for microphysics diagnostics package
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields
#ifdef MODAL_AERO
    integer                        :: l, lphase, lspec
    character(len=fieldname_len)   :: tmpname
    character(len=fieldname_len+3) :: fieldname
    character(128)                 :: long_name
    character(8)                   :: unit
#endif

  !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol      , &
                       history_microphysics_out   = history_microphysics , & 
                       history_budget_out         = history_budget       , &
                       history_budget_histfile_num_out = history_budget_histfile_num)

  ! Initialization routine for cloud macrophysics and microphysics

    if( microp_scheme .eq. 'MG' ) then
	call ini_micro
	call ini_macro
        call ini_microp_aero
    elseif( microp_scheme .eq. 'RK' ) then
        call inimc(tmelt, rhodair/1000.0_r8, gravit, rh2o)
    endif
#ifdef MODAL_AERO
    call activate_init
#endif

  ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if( convect_shallow_use_shfrc() ) then
        use_shfrc = .true.
    else 
        use_shfrc = .false.
    endif

  ! Register history variables

    do m = 1, ncnst
       call cnst_get_ind( cnst_names(m), mm )
       call addfld( cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp )
       call addfld( sflxnam  (mm), 'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp )
       call add_default( cnst_name(mm), 1, ' ' )
       call add_default( sflxnam  (mm), 1, ' ' )
    enddo

    call addfld (apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
    call addfld (apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
    call addfld (bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
    call addfld (bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

    if( history_budget) then
       call add_default (cnst_name(ixcldliq), history_budget_histfile_num, 'I')
       call add_default (cnst_name(ixcldice), history_budget_histfile_num, 'I')
       call add_default (apcnst   (ixcldliq), history_budget_histfile_num, 'I')
       call add_default (apcnst   (ixcldice), history_budget_histfile_num, 'I')
       call add_default (bpcnst   (ixcldliq), history_budget_histfile_num, 'I')
       call add_default (bpcnst   (ixcldice), history_budget_histfile_num, 'I')
    end if

    call addfld ('FWAUT    ', 'fraction', pver, 'A', 'Relative importance of liquid autoconversion'            ,phys_decomp)
    call addfld ('FSAUT    ', 'fraction', pver, 'A', 'Relative importance of ice autoconversion'               ,phys_decomp)
    call addfld ('FRACW    ', 'fraction', pver, 'A', 'Relative importance of rain accreting liquid'            ,phys_decomp)
    call addfld ('FSACW    ', 'fraction', pver, 'A', 'Relative importance of snow accreting liquid'            ,phys_decomp)
    call addfld ('FSACI    ', 'fraction', pver, 'A', 'Relative importance of snow accreting ice'               ,phys_decomp)
    call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
    call addfld ('CMEICE   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of ice within the cloud'               ,phys_decomp)
    call addfld ('CMELIQ   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of liq within the cloud'               ,phys_decomp)
    call addfld ('ICE2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of ice to precip'                     ,phys_decomp)
    call addfld ('LIQ2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of liq to precip'                     ,phys_decomp)
    call addfld ('ZMDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from ZM convection'               ,phys_decomp)
    call addfld ('DPDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from deep convection'             ,phys_decomp)
    call addfld ('DPDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from deep convection'                      ,phys_decomp)
    call addfld ('SHDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from shallow convection'          ,phys_decomp)
    call addfld ('SHDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from shallow convection'                   ,phys_decomp)
    call addfld ('DPDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to deep convective detrainment'           ,phys_decomp)
    call addfld ('SHDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to shallow convective detrainment'        ,phys_decomp)
    call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
    call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
    call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
    call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
    call addfld ('HCME     ', 'W/kg'    , pver, 'A', 'Heating from cond-evap within the cloud'                 ,phys_decomp)
    call addfld ('HEVAP    ', 'W/kg'    , pver, 'A', 'Heating from evaporation of falling precip'              ,phys_decomp)
    call addfld ('HFREEZ   ', 'W/kg'    , pver, 'A', 'Heating rate due to freezing of precip'                  ,phys_decomp)
    call addfld ('HMELT    ', 'W/kg'    , pver, 'A', 'Heating from snow melt'                                  ,phys_decomp)
    call addfld ('HREPART  ', 'W/kg'    , pver, 'A', 'Heating from cloud ice/liquid repartitioning'            ,phys_decomp)
    call addfld ('REPARTICE', 'kg/kg/s' , pver, 'A', 'Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('REPARTLIQ', 'kg/kg/s' , pver, 'A', 'Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
    call addfld ('ICWMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud water mixing ratio'                  ,phys_decomp)
    call addfld ('ICIMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud ice mixing ratio'                    ,phys_decomp)
    call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
    call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)
    call addfld ('PCSNOW   ', 'm/s     ', 1   , 'A', 'Snow fall from prognostic clouds'                        ,phys_decomp)
  ! MG microphysics diagnostics
    call addfld ('QCSEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling cloud water'              ,phys_decomp)
    call addfld ('QISEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of sublimation of falling cloud ice'                ,phys_decomp)
    call addfld ('QVRES    ', 'kg/kg/s ', pver, 'A', 'Rate of residual condensation term'                      ,phys_decomp)
    call addfld ('CMEIOUT  ', 'kg/kg/s ', pver, 'A', 'Rate of deposition/sublimation of cloud ice'             ,phys_decomp)
    call addfld ('VTRMC    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud water fallspeed'                     ,phys_decomp)
    call addfld ('VTRMI    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud ice fallspeed'                       ,phys_decomp)
    call addfld ('QCSEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud water mixing ratio tendency from sedimentation'    ,phys_decomp)
    call addfld ('QISEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud ice mixing ratio tendency from sedimentation'      ,phys_decomp)
    call addfld ('PRAO     ', '1/s     ', pver, 'A', 'Accretion of cloud water by rain'                        ,phys_decomp)
    call addfld ('PRCO     ', '1/s     ', pver, 'A', 'Autoconversion of cloud water'                           ,phys_decomp)
    call addfld ('MNUCCCO  ', '1/s     ', pver, 'A', 'Immersion freezing of cloud water'                       ,phys_decomp)
    call addfld ('MNUCCTO  ', '1/s     ', pver, 'A', 'Contact freezing of cloud water'                         ,phys_decomp)
    call addfld ('MSACWIO  ', '1/s     ', pver, 'A', 'Conversion of cloud water from rime-splintering'         ,phys_decomp)
    call addfld ('PSACWSO  ', '1/s     ', pver, 'A', 'Accretion of cloud water by snow'                        ,phys_decomp)
    call addfld ('BERGSO   ', '1/s     ', pver, 'A', 'Conversion of cloud water to snow from bergeron'         ,phys_decomp)
    call addfld ('BERGO    ', '1/s     ', pver, 'A', 'Conversion of cloud water to cloud ice from bergeron'    ,phys_decomp)
    call addfld ('MELTO    ', '1/s     ', pver, 'A', 'Melting of cloud ice'                                    ,phys_decomp)
    call addfld ('HOMOO    ', '1/s     ', pver, 'A', 'Homogeneous freezing of cloud water'                     ,phys_decomp)
    call addfld ('QCRESO   ', '1/s     ', pver, 'A', 'Residual condensation term for cloud water'              ,phys_decomp)
    call addfld ('PRCIO    ', '1/s     ', pver, 'A', 'Autoconversion of cloud ice'                             ,phys_decomp)
    call addfld ('PRAIO    ', '1/s     ', pver, 'A', 'Accretion of cloud ice by rain'                          ,phys_decomp)
    call addfld ('QIRESO   ', '1/s     ', pver, 'A', 'Residual deposition term for cloud ice'                  ,phys_decomp)
    call addfld ('MNUCCRO  ', '1/s     ', pver, 'A', 'Heterogeneous freezing of rain to snow'                  ,phys_decomp)
    call addfld ('PRACSO   ', '1/s     ', pver, 'A', 'Accretion of rain by snow'                               ,phys_decomp)
    call addfld ('MELTSDT  ', 'W/kg    ', pver, 'A', 'Latent heating rate due to melting of snow'              ,phys_decomp)
    call addfld ('FRZRDT   ', 'W/kg    ', pver, 'A', 'Latent heating rate due to homogeneous freezing of rain' ,phys_decomp)
  ! Convective cloud water variables.
    call addfld ('ICIMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud ice mixing ratio '                   ,phys_decomp)
    call addfld ('ICLMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud liquid mixing ratio '                ,phys_decomp)	
    call addfld ('ICIMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud ice mixing ratio '                        ,phys_decomp)
    call addfld ('ICLMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud liquid mixing ratio '                     ,phys_decomp)
    call addfld ('ICWMRSH  ', 'kg/kg   ', pver, 'A', 'Shallow Convection in-cloud water mixing ratio '         ,phys_decomp)
    call addfld ('ICWMRDP  ', 'kg/kg   ', pver, 'A', 'Deep Convection in-cloud water mixing ratio '            ,phys_decomp)

    call addfld ('DQSED    ', 'kg/kg/s ', pver, 'A', 'Water vapor tendency from cloud sedimentation'           ,phys_decomp)
    call addfld ('DLSED    ', 'kg/kg/s ', pver, 'A', 'Cloud liquid tendency from sedimentation'                ,phys_decomp)
    call addfld ('DISED    ', 'kg/kg/s ', pver, 'A', 'Cloud ice tendency from sedimentation'                   ,phys_decomp)
    call addfld ('HSED     ', 'W/kg    ', pver, 'A', 'Heating from cloud sediment evaporation'                 ,phys_decomp)
    call addfld ('SNOWSED  ', 'm/s     ', 1   , 'A', 'Snow from cloud ice sedimentation'                       ,phys_decomp)
    call addfld ('RAINSED  ', 'm/s     ', 1   , 'A', 'Rain from cloud liquid sedimentation'                    ,phys_decomp)
    call addfld ('PRECSED  ', 'm/s     ', 1   , 'A', 'Precipitation from cloud sedimentation'                  ,phys_decomp)

    call add_default ('FICE    ', 1, ' ')

    call addfld ('CNVCLD   ', 'fraction', 1,    'A', 'Vertically integrated convective cloud amount'           ,phys_decomp)
    call addfld ('CLDST    ', 'fraction', pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
    call addfld ('CONCLD   ', 'fraction', pver, 'A', 'Convective cloud cover'                                  ,phys_decomp)
    call addfld ('SH_CLD   ', 'fraction', pver, 'A', 'Shallow convective cloud cover'                          ,phys_decomp)
    call addfld ('DP_CLD   ', 'fraction', pver, 'A', 'Deep convective cloud cover'                             ,phys_decomp)
	
    call add_default ('CONCLD  ', 1, ' ')

    call addfld ('AST','fraction',pver, 'A','Stratus cloud fraction',phys_decomp)

  ! History variables for CAM5 macro-microphysics
    call addfld ('MPDT     ', 'W/kg    ', pver, 'A', 'Heating tendency - Morrison microphysics'                ,phys_decomp)
    call addfld ('MACPDT   ', 'W/kg    ', pver, 'A', 'Heating tendency - Revised  macrophysics'                ,phys_decomp)
    call addfld ('MPDQ     ', 'kg/kg/s ', pver, 'A', 'Q tendency - Morrison microphysics'                      ,phys_decomp)
    call addfld ('MACPDQ   ', 'kg/kg/s ', pver, 'A', 'Q tendency - Revised macrophysics'                       ,phys_decomp)
    call addfld ('MPDLIQ   ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Morrison microphysics'                 ,phys_decomp)
    call addfld ('MACPDLIQ ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Revised macrophysics'                  ,phys_decomp)
    call addfld ('MPDICE   ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Morrison microphysics'                 ,phys_decomp)
    call addfld ('MACPDICE ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Revised macrophysics'                  ,phys_decomp)
    call addfld ('MPDW2V   ', 'kg/kg/s ', pver, 'A', 'Water <--> Vapor tendency - Morrison microphysics'       ,phys_decomp)
    call addfld ('MPDW2I   ', 'kg/kg/s ', pver, 'A', 'Water <--> Ice tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDW2P   ', 'kg/kg/s ', pver, 'A', 'Water <--> Precip tendency - Morrison microphysics'      ,phys_decomp)
    call addfld ('MPDI2V   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Vapor tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDI2W   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Water tendency - Morrison microphysics'         ,phys_decomp)
    call addfld ('MPDI2P   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Precip tendency - Morrison microphysics'        ,phys_decomp)
    call addfld ('CLDVAPADJ', 'kg/kg/s ', pver, 'A', 'Q tendency associated with liq/ice adjustment - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQADJ', 'kg/kg/s ', pver, 'A', 'CLDLIQ adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDICEADJ', 'kg/kg/s ', pver, 'A', 'CLDICE adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDLIQDET', 'kg/kg/s ', pver, 'A', 'Detrainment of conv cld liq into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDICEDET', 'kg/kg/s ', pver, 'A', 'Detrainment of conv cld ice into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQLIM', 'kg/kg/s ', pver, 'A', 'CLDLIQ limiting tendency - Revised macrophysics'         ,phys_decomp)
    call addfld ('CLDICELIM', 'kg/kg/s ', pver, 'A', 'CLDICE limiting tendency - Revised macrophysics'         ,phys_decomp)
    call addfld ('LIQCLDF  ', 'fraction', pver, 'A', 'Stratus Liquid cloud fraction'                           ,phys_decomp)
    call addfld ('ICECLDF  ', 'fraction', pver, 'A', 'Stratus ICE cloud fraction'                              ,phys_decomp)
    call addfld ('IWC      ', 'kg/m3   ', pver, 'A', 'Grid box average ice water content'                      ,phys_decomp)
    call addfld ('LWC      ', 'kg/m3   ', pver, 'A', 'Grid box average liquid water content'                   ,phys_decomp)
    call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
    call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
    call addfld ('EFFLIQ   ', 'Micron  ', pver, 'A', 'Prognostic droplet effective radius'                     ,phys_decomp)
    call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
    call addfld ('EFFICE   ', 'Micron  ', pver, 'A', 'Prognostic ice effective radius'                         ,phys_decomp)
    call addfld ('WSUB     ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity'                   ,phys_decomp)
    call addfld ('WSUBI    ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity for ice'           ,phys_decomp)
    call addfld ('CDNUMC   ', '#/m2    ', 1,    'A', 'Vertically-integrated droplet concentration'             ,phys_decomp)

    if ( history_budget ) then

       call add_default ('EVAPSNOW ', history_budget_histfile_num, ' ')
       call add_default ('EVAPPREC ', history_budget_histfile_num, ' ')
       call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then

          call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')
          call add_default ('CME      ', history_budget_histfile_num, ' ')
          call add_default ('DQSED    ', history_budget_histfile_num, ' ')
          call add_default ('DISED    ', history_budget_histfile_num, ' ')
          call add_default ('DLSED    ', history_budget_histfile_num, ' ')
          call add_default ('HSED     ', history_budget_histfile_num, ' ')
          call add_default ('CMEICE   ', history_budget_histfile_num, ' ')
          call add_default ('LIQ2PR   ', history_budget_histfile_num, ' ')
          call add_default ('ICE2PR   ', history_budget_histfile_num, ' ')
          call add_default ('HCME     ', history_budget_histfile_num, ' ')
          call add_default ('HEVAP    ', history_budget_histfile_num, ' ')
          call add_default ('HFREEZ   ', history_budget_histfile_num, ' ')
          call add_default ('HMELT    ', history_budget_histfile_num, ' ')
          call add_default ('HREPART  ', history_budget_histfile_num, ' ')
          call add_default ('HPROGCLD ', history_budget_histfile_num, ' ')
          call add_default ('REPARTLIQ', history_budget_histfile_num, ' ')
          call add_default ('REPARTICE', history_budget_histfile_num, ' ')

       elseif( cam_physpkg_is('cam5') ) then

          call add_default ('QVRES    ', history_budget_histfile_num, ' ')
          call add_default ('QISEVAP  ', history_budget_histfile_num, ' ')
          call add_default ('QCSEVAP  ', history_budget_histfile_num, ' ')
          call add_default ('QISEDTEN ', history_budget_histfile_num, ' ')
          call add_default ('QCSEDTEN ', history_budget_histfile_num, ' ')
          call add_default ('QIRESO   ', history_budget_histfile_num, ' ')
          call add_default ('QCRESO   ', history_budget_histfile_num, ' ')
          call add_default ('PSACWSO  ', history_budget_histfile_num, ' ')
          call add_default ('PRCO     ', history_budget_histfile_num, ' ')
          call add_default ('PRCIO    ', history_budget_histfile_num, ' ')
          call add_default ('PRAO     ', history_budget_histfile_num, ' ')
          call add_default ('PRAIO    ', history_budget_histfile_num, ' ')
          call add_default ('PRACSO   ', history_budget_histfile_num, ' ')
          call add_default ('MSACWIO  ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2V   ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2P   ', history_budget_histfile_num, ' ')
          call add_default ('MPDW2I   ', history_budget_histfile_num, ' ')
          call add_default ('MPDT     ', history_budget_histfile_num, ' ')
          call add_default ('MPDQ     ', history_budget_histfile_num, ' ')
          call add_default ('MPDLIQ   ', history_budget_histfile_num, ' ')
          call add_default ('MPDICE   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2W   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2V   ', history_budget_histfile_num, ' ')
          call add_default ('MPDI2P   ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCTO  ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCRO  ', history_budget_histfile_num, ' ')
          call add_default ('MNUCCCO  ', history_budget_histfile_num, ' ')
          call add_default ('MELTSDT  ', history_budget_histfile_num, ' ')
          call add_default ('MELTO    ', history_budget_histfile_num, ' ')
          call add_default ('MACPDT   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDQ   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDLIQ ', history_budget_histfile_num, ' ')
          call add_default ('MACPDICE ', history_budget_histfile_num, ' ')
          call add_default ('HOMOO    ', history_budget_histfile_num, ' ')
          call add_default ('FRZRDT   ', history_budget_histfile_num, ' ')
          call add_default ('CMEIOUT  ', history_budget_histfile_num, ' ')
          call add_default ('CLDVAPADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQLIM', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQDET', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDICELIM', history_budget_histfile_num, ' ')
          call add_default ('CLDICEDET', history_budget_histfile_num, ' ')
          call add_default ('CLDICEADJ', history_budget_histfile_num, ' ')
          call add_default ('BERGSO   ', history_budget_histfile_num, ' ')
          call add_default ('BERGO    ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFT   ', history_budget_histfile_num, ' ')

       end if

    end if

  ! Averaging for cloud particle number and size
    call addfld ('AWNC     ', 'm-3     ', pver, 'A', 'Average cloud water number conc'                         ,phys_decomp)
    call addfld ('AWNI     ', 'm-3     ', pver, 'A', 'Average cloud ice number conc'                           ,phys_decomp)
    call addfld ('AREL     ', 'Micron  ', pver, 'A', 'Average droplet effective radius'                        ,phys_decomp)
    call addfld ('AREI     ', 'Micron  ', pver, 'A', 'Average ice effective radius'                            ,phys_decomp)
  ! Frequency arrays for above
    call addfld ('FREQL    ', 'fraction', pver, 'A', 'Fractional occurance of liquid'                          ,phys_decomp)
    call addfld ('FREQI    ', 'fraction', pver, 'A', 'Fractional occurance of ice'                             ,phys_decomp)

    if( history_microphysics) then
        call add_default ('CDNUMC   ', 1, ' ')
        call add_default ('IWC      ', 1, ' ')
        call add_default ('WSUB     ', 1, ' ')
        call add_default ('FREQL    ', 1, ' ')
        call add_default ('FREQI    ', 1, ' ')
        call add_default ('AREI     ', 1, ' ')
        call add_default ('AREL     ', 1, ' ')
        call add_default ('AWNC     ', 1, ' ')
        call add_default ('AWNI     ', 1, ' ')
    endif

  ! Average cloud top particle size and number (liq, ice) and frequency
    call addfld ('ACTREL   ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet effective radius'              ,phys_decomp)
    call addfld ('ACTREI   ', 'Micron  ', 1,    'A', 'Average Cloud Top ice effective radius'                  ,phys_decomp)
    call addfld ('ACTNL    ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet number'                        ,phys_decomp)
    call addfld ('ACTNI    ', 'Micron  ', 1,    'A', 'Average Cloud Top ice number'                            ,phys_decomp)

    call addfld ('FCTL     ', 'fraction', 1,    'A', 'Fractional occurance of cloud top liquid'                ,phys_decomp)
    call addfld ('FCTI     ', 'fraction', 1,    'A', 'Fractional occurance of cloud top ice'                   ,phys_decomp)

    call add_default ('ICWMR', 1, ' ')
    call add_default ('ICIMR', 1, ' ')

    call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain flux', phys_decomp)
    call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)

#ifdef MODAL_AERO
! Add dropmixnuc tendencies for all modal aerosol species
    do m = 1, ntot_amode
    do lphase = 1, 2
    do lspec = 0, nspec_amode(m)+1   ! loop over number + chem constituents + water
       unit = 'kg/m2/s'
       if (lspec == 0) then   ! number
          unit = '#/m2/s'
          if (lphase == 1) then
             l = numptr_amode(m)
          else
             l = numptrcw_amode(m)
          endif
       else if (lspec <= nspec_amode(m)) then   ! non-water mass
          if (lphase == 1) then
             l = lmassptr_amode(lspec,m)
          else
             l = lmassptrcw_amode(lspec,m)
          endif
       else   ! water mass
          cycle
       end if
       if (lphase == 1) then
          tmpname = cnst_name(l)
       else
          tmpname = cnst_name_cw(l)
       end if

       fieldname = trim(tmpname) // '_mixnuc1'
       long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
       call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
       if ( history_aerosol ) then 
          call add_default( fieldname, 1, ' ' )
          if ( masterproc ) write(*,'(2a)') 'stratiform_init addfld - ', fieldname
       endif
       
    end do   ! lspec
    end do   ! lphase
    end do   ! m
#endif
    return
  end subroutine stratiform_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_tend(                             &
             state, ptend_all, dtime, icefrac, landfrac,  &
             ocnfrac, landm, snowh,                       &
#ifdef MODAL_AERO
             cflx,                                        &
#endif
             dlf, dlf2, rliq, cmfmc, cmfmc2, ts,          &
             sst, zdu, prec_str, snow_str, prec_sed,      &
             snow_sed, prec_pcw, snow_pcw, pbuf, state_eq )

  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to sedimentation, detrain, cloud fraction and !
  !        cloud macro - microphysics subroutines           !
  !                                                         ! 
  ! Author: D.B. Coleman                                    !
  ! Date: Apr 2004                                          !
  !                                                         !
  !-------------------------------------------------------- !

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid
  use cloud_fraction,   only: cldfrc
  use physics_types,    only: physics_state, physics_ptend, physics_tend
  use physics_types,    only: physics_ptend_init, physics_update, physics_tend_init
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use cam_history,      only: outfld
  use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
  use constituents,     only: cnst_get_ind, pcnst
  use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
  use cldwat,           only: pcond, cldwat_fice
  use cldwat2m_micro,   only: mmicro_pcond
  use microp_aero,      only: microp_aero_ts 
  use cldwat2m_macro,   only: mmacro_pcond
  use physconst,        only: cpair
  use rad_constituents, only: rad_cnst_get_clim_info, rad_cnst_get_clim_aer
  use time_manager,     only: is_first_step, get_nstep
  use pkg_cldoptics,    only: cldefr
#ifdef MODAL_AERO
  use modal_aero_data
#endif
  ! Debug
    use phys_debug_util,  only: phys_debug_col
  ! Debug

  implicit none

  ! Debug
    integer icol
  ! Debug

!
! Parameters
!
  real(r8) pnot                  ! Reference pressure
  parameter (pnot = 1.e5_r8)

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend_all   ! Package tendencies
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

  real(r8), intent(in)  :: dtime                    ! Timestep
  real(r8), intent(in)  :: icefrac (pcols)          ! Sea ice fraction (fraction)
  real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
  real(r8), intent(in)  :: landm(pcols)             ! Land fraction ramped over water
  real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)
#ifdef MODAL_AERO
  real(r8), intent(in)  :: cflx(pcols,pcnst)        ! Constituent flux from surface
#endif

  real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
  real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
  real(r8), intent(in)  :: rliq(pcols)              ! Vertical integral of liquid not yet in q(ixcldliq)
  real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
  real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

  real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
  real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
  real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection

  real(r8), intent(out) :: prec_str(pcols)          ! [Total] Sfc flux of precip from stratiform [ m/s ] 
  real(r8), intent(out) :: snow_str(pcols)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
  real(r8), intent(out) :: prec_sed(pcols)          ! Surface flux of total cloud water from sedimentation
  real(r8), intent(out) :: snow_sed(pcols)          ! Surface flux of cloud ice from sedimentation
  real(r8), intent(out) :: prec_pcw(pcols)          ! Sfc flux of precip from microphysics [ m/s ]
  real(r8), intent(out) :: snow_pcw(pcols)          ! Sfc flux of snow from microphysics [ m/s ]

  ! Equilibrium state variables at the end of macrophysics
  ! Below 'state_eq' is for future use as the input of radiation'PBL scheme

  type(physics_state), intent(out) :: state_eq   ! Equilibrium state variables at the end of macrophysics

  !
  ! Local variables
  !

  type(physics_state)   :: state1                   ! Local copy of the state variable
  type(physics_tend )   :: tend                     ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend)   :: ptend_loc                ! Package tendencies

  integer i,k,m
  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns
  integer :: conv_water_in_rad

  ! Physics buffer fields

  integer itim, ifld
  real(r8), pointer, dimension(:,:) :: rhdfda       !
  real(r8), pointer, dimension(:,:) :: rhu00        ! 
  real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
  real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
  real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: CC_T         ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qv        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ql        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qi        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_nl        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ni        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qlst      ! In-liquid stratus microphysical tendency
  real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
  real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
  real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
  real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
  real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
  real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
  real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
  real(r8), pointer, dimension(:,:) :: qme
  real(r8), pointer, dimension(:,:) :: prain        ! Total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: nevapr       ! Evaporation of total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: rel          ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei          ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: rel2         ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei2         ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: cldo         ! Old cloud fraction
  real(r8), pointer, dimension(:,:) :: kkvh         ! Vertical eddy diffusivity
  real(r8), pointer, dimension(:,:) :: wsedl        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
#ifdef MODAL_AERO
  real(r8), pointer, dimension(:,:,:) :: qqcw       ! Cloud-borne aerosol
  real(r8), pointer, dimension(:,:,:) :: dgnumwet   ! Number mode diameter
  real(r8), pointer, dimension(:,:,:) :: dgnum      ! Number mode diameter
  real(r8), pointer, dimension(:,:) :: rate1ord_cw2pr_st   ! 1st order rate for direct conversion of
                                                    ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
#endif
  real(r8) :: shfrc(pcols,pver)                     ! Cloud fraction from shallow convection scheme
  real(r8), pointer, dimension(:,:) :: rel_fn       ! Ice effective drop size at fixed number (indirect effect) (microns)

  real(r8)  rate1cld(pcols,pver)                    ! array to hold rate1ord_cw2pr_st from microphysics

  ! physics buffer fields for COSP simulator (RK only)

  real(r8), pointer, dimension(:,:) :: rkflxprc     ! RK grid-box mean flux_large_scale_cloud_rain at interfaces (kg/m2/s)
  real(r8), pointer, dimension(:,:) :: rkflxsnw     ! RK grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)

  ! physics buffer fields for radiation

  real(r8), pointer, dimension(:,:) :: dei          ! Ice effective diameter (meters) (AG: microns?)
  real(r8), pointer, dimension(:,:) :: mu           ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdac      ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwp        ! In-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwp        ! In-cloud liquid water path for radiation
  
  ! For rrtm optics. specificed distribution.

  real(r8) :: mucon                                 ! Convective size distribution shape parameter
  real(r8) :: dcon                                  ! Convective size distribution effective radius (meters)
  real(r8) :: lamcon                                ! Convective size distribution slope parameter (meters-1)
  real(r8) :: deicon                                ! Convective ice effective diameter (meters)

  ! Physics buffer fields

  real(r8), pointer, dimension(:,:) :: deiconv      ! Ice effective diameter (microns)
  real(r8), pointer, dimension(:,:) :: muconv       ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdaconv   ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwpst      ! Stratiform in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpst      ! Stratiform in-cloud liquid water path for radiation
  real(r8), pointer, dimension(:,:) :: iciwpconv    ! Convective in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpconv    ! Convective in-cloud liquid water path for radiation

  real(r8), pointer, dimension(:,:) :: tke          ! TKE from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: turbtype     ! Turbulence type from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: smaw         ! Instability function of momentum from the moist PBL scheme

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: concld_ql    ! Convective cloud
  real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

  ! Local variables for in-cloud water quantities adjusted for convective water

  real(r8)  allcld_ice (pcols,pver)                 ! All-cloud cloud ice
  real(r8)  allcld_liq (pcols,pver)                 ! All-cloud liquid

  ! Snow

  real(r8), pointer, dimension(:,:) :: cldfsnow     ! Cloud fraction for liquid+snow
  real(r8), pointer, dimension(:,:) :: icswp        ! In-cloud snow water path
  real(r8), pointer, dimension(:,:) :: des          ! Snow effective diameter (m)
  real(r8)  qsout(pcols,pver)                       ! Snow mixing ratio
 
  ! Local variables for stratiform_sediment

  real(r8)  rain(pcols)                             ! Surface flux of cloud liquid
  real(r8)  pvliq(pcols,pver+1)                     ! Vertical velocity of cloud liquid drops (Pa/s)
  real(r8)  pvice(pcols,pver+1)                     ! Vertical velocity of cloud ice particles (Pa/s)

  ! Local variables for cldfrc

  real(r8)  cldst(pcols,pver)                       ! Stratus cloud fraction
  real(r8)  rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
  real(r8)  rhcloud2(pcols,pver)                    ! Relative humidity cloud (perturbation)
  real(r8)  clc(pcols)                              ! Column convective cloud amount
  real(r8)  relhum(pcols,pver)                      ! RH, output to determine drh/da
  real(r8)  rhu002(pcols,pver)                      ! Same as rhu00 but for perturbed rh 
  real(r8)  cld2(pcols,pver)                        ! Same as cld but for perturbed rh
  real(r8)  concld2(pcols,pver)                     ! Same as concld but for perturbed rh 
  real(r8)  cldst2(pcols,pver)                      ! Same as cldst but for perturbed rh 
  real(r8)  relhum2(pcols,pver)                     ! RH after  perturbation            
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf_out(pcols,pver)                 ! Ice cloud fraction
  real(r8)  liqcldf_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf2(pcols,pver)                    ! Ice cloud fraction
  real(r8)  liqcldf2(pcols,pver)                    ! Liquid cloud fraction (combined into cloud)

  ! Local variables for microphysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
  real(r8)  evapheat(pcols,pver)                    ! Heating rate due to evaporation of precip
  real(r8)  evapsnow(pcols,pver)                    ! Local evaporation of snow
  real(r8)  prfzheat(pcols,pver)                    ! Heating rate due to freezing of precip (W/kg)
  real(r8)  meltheat(pcols,pver)                    ! Heating rate due to phase change of precip
  real(r8)  cmeheat (pcols,pver)                    ! Heating rate due to phase change of precip
  real(r8)  prodsnow(pcols,pver)                    ! Local production of snow
  real(r8)  totcw(pcols,pver)                       ! Total cloud water mixing ratio
  real(r8)  fice(pcols,pver)                        ! Fractional ice content within cloud
  real(r8)  fsnow(pcols,pver)                       ! Fractional snow production
  real(r8)  repartht(pcols,pver)                    ! Heating rate due to phase repartition of input precip
  real(r8)  icimr(pcols,pver)                       ! In cloud ice mixing ratio
  real(r8)  icwmr(pcols,pver)                       ! In cloud water mixing ratio
  real(r8)  icimrst(pcols,pver)                     ! In stratus ice mixing ratio
  real(r8)  icwmrst(pcols,pver)                     ! In stratus water mixing ratio
  real(r8)  icimrst_out(pcols,pver)                 ! In stratus ice mixing ratio
  real(r8)  icwmrst_out(pcols,pver)                 ! In stratus water mixing ratio
  real(r8)  fwaut(pcols,pver)              
  real(r8)  fsaut(pcols,pver)              
  real(r8)  fracw(pcols,pver)              
  real(r8)  fsacw(pcols,pver)              
  real(r8)  fsaci(pcols,pver)              
  real(r8)  cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
  real(r8)  cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud
  real(r8)  ice2pr(pcols,pver)                      ! Rate of conversion of ice to precip
  real(r8)  liq2pr(pcols,pver)                      ! Rate of conversion of liquid to precip
  real(r8)  liq2snow(pcols,pver)                    ! Rate of conversion of liquid to snow
  real(r8)  temp(pcols)
  real(r8)  res(pcols,pver)
  real(r8)  droprad                                 ! Radius of droplets detrained from cumulus (m)
  real(r8)  invdropmass                             ! Inverse of mean droplet mass (#/kg)

! MG micro diagnostics

  real(r8)  qcsevap(pcols,pver)                     ! Evaporation of falling cloud water
  real(r8)  qisevap(pcols,pver)                     ! Sublimation of falling cloud ice
  real(r8)  qvres(pcols,pver)                       ! Residual condensation term to remove excess saturation
  real(r8)  cmeiout(pcols,pver)                     ! Deposition/sublimation rate of cloud ice
  real(r8)  vtrmc(pcols,pver)                       ! Mass-weighted cloud water fallspeed
  real(r8)  vtrmi(pcols,pver)                       ! Mass-weighted cloud ice fallspeed
  real(r8)  qcsedten(pcols,pver)                    ! Cloud water mixing ratio tendency from sedimentation
  real(r8)  qisedten(pcols,pver)                    ! Cloud ice mixing ratio tendency from sedimentation

  real(r8)  prao(pcols,pver)  
  real(r8)  prco(pcols,pver)  
  real(r8)  mnuccco(pcols,pver)  
  real(r8)  mnuccto(pcols,pver)  
  real(r8)  msacwio(pcols,pver)  
  real(r8)  psacwso(pcols,pver)  
  real(r8)  bergso(pcols,pver)  
  real(r8)  bergo(pcols,pver)  
  real(r8)  melto(pcols,pver)  
  real(r8)  homoo(pcols,pver)  
  real(r8)  qcreso(pcols,pver)  
  real(r8)  prcio(pcols,pver)  
  real(r8)  praio(pcols,pver)  
  real(r8)  qireso(pcols,pver)
  real(r8)  ftem(pcols,pver)
  real(r8)  mnuccro(pcols,pver) 
  real(r8)  pracso (pcols,pver) 
  real(r8)  meltsdt(pcols,pver) 
  real(r8)  frzrdt (pcols,pver) 
  real(r8)  dpdlfliq(pcols,pver)
  real(r8)  dpdlfice(pcols,pver)
  real(r8)  shdlfliq(pcols,pver)
  real(r8)  shdlfice(pcols,pver)
  real(r8)  dpdlft  (pcols,pver)
  real(r8)  shdlft  (pcols,pver)

#ifdef MODAL_AERO
  integer l, lnum, lnumcw, lmass, lmasscw
#endif

  ! Variables for MG microphysics

  real(r8)  dum1,dum2
  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)
  real(r8)  icinc(pcols,pver)                       ! In cloud ice number conc
  real(r8)  cdnumc(pcols)                           ! Vertically-integrated droplet concentration
  real(r8)  icwnc(pcols,pver)                       ! In cloud water number conc
  real(r8)  iwc(pcols,pver)                         ! Grid box average ice water content
  real(r8)  lwc(pcols,pver)                         ! Grid box average liquid water content  
  real(r8)  effliq(pcols,pver)                      ! In cloud liq eff rad
  real(r8)  effice(pcols,pver)                      ! In cloud ice eff rad
  real(r8)  effliq_fn(pcols,pver)                   ! In cloud liq eff rad at fixed number concentration	
  real(r8)  wsub(pcols,pver)                        ! Sub-grid vertical velocity (m/s)
  real(r8)  wsubi(pcols,pver)                       ! Sub-grid vertical velocity for ice (m/s)

  ! Output from mmicro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)
  real(r8)  effc(pcols,pver)
  real(r8)  effc_fn(pcols,pver)                     ! Liquid effective radius at fixed number (for indirect calc)
  real(r8)  effi(pcols,pver)
  real(r8)  prect(pcols)
  real(r8)  preci(pcols)

  ! Output from microp_aero_ts for aerosol actication
  real(r8)  naai(pcols,pver)      !ice nucleation number
  real(r8)  npccn(pcols,pver)     !liquid activation number
  real(r8)  rndst(pcols,pver,4)
  real(r8)  nacon(pcols,pver,4)

  ! Output from mmacro_pcond

  real(r8)  qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
  real(r8)  qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
  real(r8)  qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
  real(r8)  qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
  real(r8)  qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

  ! Averaging arrays for effective radius and number....

  real(r8)  efiout(pcols,pver)
  real(r8)  efcout(pcols,pver)
  real(r8)  ncout(pcols,pver)
  real(r8)  niout(pcols,pver)

  real(r8)  freqi(pcols,pver)
  real(r8)  freql(pcols,pver)

  ! Average cloud top radius & number

  real(r8)  ctrel(pcols)
  real(r8)  ctrei(pcols)
  real(r8)  ctnl(pcols)
  real(r8)  ctni(pcols)
  real(r8)  fcti(pcols)
  real(r8)  fctl(pcols)

  ! Gather mass mixing ratio for all aerosols affecting the climate

  integer               :: naer_all
  real(r8), pointer     :: aermmr1(:,:)
  real(r8), allocatable :: aer_mmr(:,:,:)           ! Aerosol mass mixing ratio

  ! For revised macophysics, mmacro_pcond

  real(r8)  itend(pcols,pver)
  real(r8)  lmitend(pcols,pver)
  real(r8)  zeros(pcols,pver)
  real(r8)  t_inout(pcols,pver)
  real(r8)  qv_inout(pcols,pver)
  real(r8)  ql_inout(pcols,pver)
  real(r8)  qi_inout(pcols,pver)
  real(r8)  prsed(pcols,pver)
  real(r8)  pssed(pcols,pver)
  real(r8)  ersed(pcols,pver)
  real(r8)  essed(pcols,pver)
  real(r8)  alst_mic(pcols,pver)
  real(r8)  aist_mic(pcols,pver)
  real(r8)  concld_old(pcols,pver)

  real(r8)  nl_inout(pcols,pver)
  real(r8)  ni_inout(pcols,pver)
  real(r8)  dum1D(pcols)
  real(r8)  nltend(pcols,pver)
  real(r8)  nitend(pcols,pver)

  real(r8)  zero1D(pcols)
  real(r8)  t_out(pcols,pver)
  real(r8)  qv_out(pcols,pver)
  real(r8)  ql_out(pcols,pver)
  real(r8)  qi_out(pcols,pver)
  real(r8)  nl_out(pcols,pver)
  real(r8)  ni_out(pcols,pver)
  real(r8)  QQw(pcols,pver)
  real(r8)  QQi(pcols,pver)
  real(r8)  QQnl(pcols,pver)
  real(r8)  QQni(pcols,pver)

  ! For detraining cumulus condensate into the 'stratus' without evaporation
  ! This is for use in mmacro_pcond

  real(r8)  dlf_T(pcols,pver)
  real(r8)  dlf_qv(pcols,pver)
  real(r8)  dlf_ql(pcols,pver)
  real(r8)  dlf_qi(pcols,pver)
  real(r8)  dlf_nl(pcols,pver)
  real(r8)  dlf_ni(pcols,pver)

  real(r8)  rel_detcu
  real(r8)  rei_detcu

  ! ======================================================================

  lchnk = state%lchnk
  ncol  = state%ncol

  call phys_getopts( conv_water_in_rad_out = conv_water_in_rad )

  call physics_state_copy(state,state1)             ! Copy state to local state1.
  call physics_ptend_init(ptend_loc)                ! Initialize local ptend type
  call physics_ptend_init(ptend_all)                ! Initialize output ptend type
  call physics_tend_init(tend)                      ! tend here is just a null place holder

  ! Associate pointers with physics buffer fields

  itim = pbuf_old_tim_idx()

  ifld = pbuf_get_fld_idx('QCWAT')
  qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('RHDFDA')
  rhdfda => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('RHU00')
  rhu00 => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('TCWAT')
  tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('LCWAT')
  lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('ICCWAT')
  iccwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('NLWAT')
  nlwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('NIWAT')
  niwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_T')
  CC_T => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_qv')
  CC_qv => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_ql')
  CC_ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_qi')
  CC_qi => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_nl')
  CC_nl => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_ni')
  CC_ni => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CC_qlst')
  CC_qlst => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CLD')
  cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('AST')
  ast => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('AIST')
  aist => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('ALST')
  alst => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('QIST')
  qist => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('QLST')
  qlst => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CONCLD')  
  concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('CLDO')
  cldo => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('REL2')
  rel2 => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ifld = pbuf_get_fld_idx('REI2')
  rei2 => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

#ifdef MODAL_AERO

  ifld = pbuf_get_fld_idx('QQCW')
  qqcw => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:pcnst)

  ifld = pbuf_get_fld_idx('DGNUMWET')
  dgnumwet => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:maxd_amode)

  ifld = pbuf_get_fld_idx('DGNUM' )
  dgnum => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:maxd_amode)

  ifld = pbuf_get_fld_idx('RATE1_CW2PR_ST')                             
  rate1ord_cw2pr_st => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)    

  if( is_first_step() ) then
      qqcw(1:pcols,1:pver,1:pcnst)          = 1.e-38_r8
      dgnumwet(1:pcols,1:pver,1:maxd_amode) = 0.0_r8
      dgnum(1:pcols,1:pver,1:maxd_amode)    = 0.0_r8
      rate1ord_cw2pr_st(1:pcols,1:pver)     = 0.0_r8   
  endif

#endif

! For purposes of convective ql.

  ifld = pbuf_get_fld_idx('CONCLDQL')                        
  concld_ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)                                      

  ifld = pbuf_get_fld_idx('FICE')                        
  fice_ql => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  ifld = pbuf_get_fld_idx('QME')
  qme  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('PRAIN')
  prain  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('NEVAPR')
  nevapr  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('REL')
  rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('REI')
  rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('REL_FN')
  rel_fn  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('kvh')
  kkvh => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  ifld = pbuf_get_fld_idx('DEI')
  dei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('MU')
  mu  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('LAMBDAC')
  lambdac  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICIWP')
  iciwp  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICLWP')
  iclwp  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('DEICONV')
  deiconv  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('MUCONV')
  muconv  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('LAMBDACONV')
  lambdaconv  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICIWPST')
  iciwpst  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICLWPST')
  iclwpst  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICIWPCONV')
  iciwpconv  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICLWPCONV')
  iclwpconv  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('DES')
  des  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('ICSWP')
  icswp  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ifld = pbuf_get_fld_idx('CLDFSNOW')
  cldfsnow  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, itim)

  ifld = pbuf_get_fld_idx('tke')
  tke => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  ifld = pbuf_get_fld_idx('turbtype')
  turbtype => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  ifld = pbuf_get_fld_idx('smaw')
  smaw => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  ifld = pbuf_get_fld_idx('WSEDL')
  wsedl  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ! If first timestep, initialize heatflux....in pbuf at all time levels.

  if( is_first_step() ) then
      kkvh(:,:)     = 0._r8
      tke(:,:)      = 0._r8
      turbtype(:,:) = 0._r8
      smaw(:,:)     = 0._r8
  endif

  ! Assign default size distribution parameters for no-stratiform clouds (convection only)
  ! Also put into physics buffer for possible separate use by radiation

  dcon   = 25.e-6_r8
  mucon  = 5.3_r8
  deicon = 50._r8

  muconv(:,:)     = mucon
  lambdaconv(:,:) = (mucon + 1._r8)/dcon
  deiconv(:,:)    = deicon

  ! Initialize convective detrainment tendency

  dlf_T(:,:)  = 0._r8
  dlf_qv(:,:) = 0._r8
  dlf_ql(:,:) = 0._r8
  dlf_qi(:,:) = 0._r8
  dlf_nl(:,:) = 0._r8
  dlf_ni(:,:) = 0._r8

   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !

   ! ------------- !
   ! Sedimentation !
   ! ------------- !

   if( microp_scheme .eq. 'RK' ) then

     ! Allow the cloud liquid drops and ice particles to sediment.
     ! This is done before adding convectively detrained cloud water, 
     ! because the phase of the detrained water is unknown.

       call t_startf('stratiform_sediment')

       ptend_loc%name         = 'pcwsediment'
       ptend_loc%ls           = .TRUE.
       ptend_loc%lq(1)        = .TRUE.
       ptend_loc%lq(ixcldice) = .TRUE.
       ptend_loc%lq(ixcldliq) = .TRUE.

       call cld_sediment_vel( ncol,                                                           &
                              icefrac, landfrac, ocnfrac, state1%pmid, state1%pdel, state1%t, &
                              cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice),            & 
                              pvliq, pvice, landm, snowh )

       wsedl(:ncol,:pver) = pvliq(:ncol,:pver)/gravit/(state1%pmid(:ncol,:pver)/(287.15_r8*state1%t(:ncol,:pver)))

       call cld_sediment_tend( ncol, dtime ,                                                             &
                               state1%pint, state1%pmid, state1%pdel, state1%t,                          &
                               cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice), pvliq, pvice,        &
                               ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,1), &
                               ptend_loc%s, rain, snow_sed )

     ! Convert rain and snow fluxes at the surface from [kg/m2/s] to [m/s]
     ! Compute total precipitation flux at the surface in [m/s]

       snow_sed(:ncol) = snow_sed(:ncol)/1000._r8
       rain(:ncol)     = rain(:ncol)/1000._r8
       prec_sed(:ncol) = rain(:ncol) + snow_sed(:ncol)

     ! Record history variables
       lchnk = state1%lchnk
       call outfld( 'DQSED'   ,ptend_loc%q(:,:,1)       , pcols,lchnk )
       call outfld( 'DISED'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk )
       call outfld( 'DLSED'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk )
       call outfld( 'HSED'    ,ptend_loc%s              , pcols,lchnk )
       call outfld( 'PRECSED' ,prec_sed                 , pcols,lchnk )
       call outfld( 'SNOWSED' ,snow_sed                 , pcols,lchnk )
       call outfld( 'RAINSED' ,rain                     , pcols,lchnk )

     ! Add tendency from this process to tend from other processes here
       call physics_ptend_sum( ptend_loc, ptend_all, state )

     ! Update physics state type state1 with ptend_loc 
       call physics_update( state1, tend, ptend_loc, dtime )
       call physics_ptend_init( ptend_loc )

       call t_stopf('stratiform_sediment')

     ! Accumulate prec and snow flux at the surface [ m/s ]
       prec_str(:ncol) = prec_sed(:ncol)
       snow_str(:ncol) = snow_sed(:ncol)

   endif  ! End of 'Sediment'

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   ptend_loc%name = 'pcwdetrain'

   if( microp_scheme .eq. 'MG' ) then

       ptend_loc%lq(ixcldliq) = .TRUE.
       ptend_loc%lq(ixcldice) = .TRUE.
       ptend_loc%lq(ixnumliq) = .TRUE.
       ptend_loc%lq(ixnumice) = .TRUE.
       ptend_loc%ls           = .TRUE.

     ! Procedures :
     ! (1) Partition detrained convective cloud water into liquid and ice based on T.
     !     This also involves heating.
     !     If convection scheme can handle this internally, this step is not necssary.
     ! (2) Assuming a certain effective droplet radius, computes number concentration
     !     of detrained convective cloud liquid and ice.
     ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into 
     !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
     !     not involve any macrophysical evaporation while the latter does. This is
     !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded 
     !     by qcst_min and qcst_max in mmacro_pcond.
     ! (4) In contrast to liquid, convective ice is detrained into the environment 
     !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
     ! This is the key procesure generating upper-level cirrus clouds.
     ! The unit of dlf : [ kg/kg/s ]

       do k = 1, pver
       do i = 1, state1%ncol
          if( state1%t(i,k) > 268.15_r8 ) then
              dum1 = 0.0_r8
          elseif( state1%t(i,k) < 238.15_r8 ) then
              dum1 = 1.0_r8
          else
              dum1 = ( 268.15_r8 - state1%t(i,k) ) / 30._r8
          endif
          ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
          ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
        ! dum2                      = dlf(i,k) * ( 1._r8 - dum1 )
          ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) / (4._r8*3.14_r8* 8.e-6_r8**3*997._r8) + & ! Deep    Convection
                                      3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) / (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection 
        ! dum2                      = dlf(i,k) * dum1
          ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) / (4._r8*3.14_r8*25.e-6_r8**3*500._r8) + & ! Deep    Convection
                                      3._r8 * (                         dlf2(i,k)    *  dum1 ) / (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
          ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice
        ! Targetted detrainment of convective liquid water either directly into the
        ! existing liquid stratus or into the environment. 
          if( cu_det_st ) then
              dlf_T(i,k)  = ptend_loc%s(i,k)/cpair
              dlf_qv(i,k) = 0._r8
              dlf_ql(i,k) = ptend_loc%q(i,k,ixcldliq)
              dlf_qi(i,k) = ptend_loc%q(i,k,ixcldice)
              dlf_nl(i,k) = ptend_loc%q(i,k,ixnumliq)
              dlf_ni(i,k) = ptend_loc%q(i,k,ixnumice)         
              ptend_loc%q(i,k,ixcldliq) = 0._r8
              ptend_loc%q(i,k,ixcldice) = 0._r8
              ptend_loc%q(i,k,ixnumliq) = 0._r8
              ptend_loc%q(i,k,ixnumice) = 0._r8
              ptend_loc%s(i,k)          = 0._r8
              dpdlfliq(i,k)             = 0._r8
              dpdlfice(i,k)             = 0._r8
              shdlfliq(i,k)             = 0._r8
              shdlfice(i,k)             = 0._r8
              dpdlft  (i,k)             = 0._r8
              shdlft  (i,k)             = 0._r8
           else
              dpdlfliq(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( 1._r8 - dum1 )
              dpdlfice(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( dum1 )
              shdlfliq(i,k) = dlf2(i,k) * ( 1._r8 - dum1 )
              shdlfice(i,k) = dlf2(i,k) * ( dum1 )
              dpdlft  (i,k) = ( dlf(i,k) - dlf2(i,k) ) * dum1 * latice/cpair
              shdlft  (i,k) = dlf2(i,k) * dum1 * latice/cpair
          endif
       end do
       end do

       call outfld( 'DPDLFLIQ ', dpdlfliq, pcols, lchnk )
       call outfld( 'DPDLFICE ', dpdlfice, pcols, lchnk )
       call outfld( 'SHDLFLIQ ', shdlfliq, pcols, lchnk )
       call outfld( 'SHDLFICE ', shdlfice, pcols, lchnk )
       call outfld( 'DPDLFT   ', dpdlft  , pcols, lchnk )
       call outfld( 'SHDLFT   ', shdlft  , pcols, lchnk )

   else if ( microp_scheme .eq. 'RK' ) then

     ! Put all of the detraining cloud water from convection into the large scale cloud.
     ! It all goes in liquid for the moment.
     ! Strictly speaking, this approach is detraining all the cconvective water into 
     ! the environment, not the large-scale cloud.

       ptend_loc%lq(ixcldliq) = .TRUE.
       do k = 1, pver
       do i = 1, state1%ncol
          ptend_loc%q(i,k,ixcldliq) = dlf(i,k)
       end do
       end do

   end if

   call outfld( 'ZMDLF', dlf, pcols, state1%lchnk )

 ! Add hie detrainment tendency to tend from the other prior processes

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc )

 ! Accumulate prec and snow, reserved liquid has now been used.
 ! For MG, this is performed later.

   if( microp_scheme .eq. 'RK' ) then
       prec_str(:ncol) = prec_str(:ncol) - rliq(:ncol)  ! ( snow contribution is zero )
   endif

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !  
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- ! 

   concld_old(:ncol,:pver) = concld(:ncol,:pver)

   if( use_shfrc ) then
       ifld = pbuf_get_fld_idx('shfrc')
       shfrc(:pcols,:pver) = pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   else 
       shfrc(:,:) = 0._r8
   endif

   ! CAM5 only uses 'concld' output from the below subroutine. 
   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 
   ! Note 'shfrc' and' deep convective cloud fraction' will be saved into the 
   ! physical buffer (SH_FRAC,DP_FRAC) within cldfrc.

   call t_startf("cldfrc")
   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld, rhcloud, clc, state1%pdel,                                    &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                relhum, 0 )    

   ! Re-calculate cloud with perturbed rh add call cldfrc  

   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld2, rhcloud2, clc, state1%pdel,                                  &
                cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                relhum2, 1 )              

   call t_stopf("cldfrc")

   ! Add following to estimate rhdfda. Below block is only for CAM4                       

   rhu00(:ncol,1) = 2.0_r8

   do k = 1, pver
   do i = 1, ncol
      if( relhum(i,k) < rhu00(i,k) ) then
          rhdfda(i,k) = 0.0_r8
      elseif( relhum(i,k) >= 1.0_r8 ) then
          rhdfda(i,k) = 0.0_r8
      else
         ! Under certain circumstances, rh+ cause cld not to changed
         ! when at an upper limit, or w/ strong subsidence
         if( ( cld2(i,k) - cld(i,k) ) < 1.e-4_r8 ) then
               rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8  
         else
               rhdfda(i,k) = 0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
         endif
      endif
   enddo
   enddo

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   call t_startf('stratiform_microphys')

   lchnk  = state1%lchnk
   ncol   = state1%ncol
   rdtime = 1._r8/dtime

 ! Define fractional amount of stratus condensate and precipitation in ice phase.
 ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
 ! The ramp within convective cloud may be different

   call cldwat_fice( ncol, state1%t, fice, fsnow )

   if( microp_scheme .eq. 'RK' ) then

     ! Perform repartitioning of stratiform condensate.    
     ! Corresponding heating tendency will be added later. 

       totcw(:ncol,:pver)     = state1%q(:ncol,:pver,ixcldice) + state1%q(:ncol,:pver,ixcldliq)
       repartht(:ncol,:pver)  = state1%q(:ncol,:pver,ixcldice)
       ptend_loc%q(:ncol,:pver,ixcldice) = rdtime * ( totcw(:ncol,:pver)*fice(:ncol,:pver)          - state1%q(:ncol,:pver,ixcldice) )
       ptend_loc%q(:ncol,:pver,ixcldliq) = rdtime * ( totcw(:ncol,:pver)*(1.0_r8-fice(:ncol,:pver)) - state1%q(:ncol,:pver,ixcldliq) )

       call outfld( 'REPARTICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk )
       call outfld( 'REPARTLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk )

       ptend_loc%name         = 'cldwat-repartition'
       ptend_loc%lq(ixcldice) = .true.
       ptend_loc%lq(ixcldliq) = .true.

       call physics_ptend_sum( ptend_loc, ptend_all, state )
       call physics_update( state1, tend, ptend_loc, dtime )
       call physics_ptend_init( ptend_loc )

   endif

   ptend_loc%name         = 'cldwat'
   ptend_loc%ls           = .true.
   ptend_loc%lq(1)        = .true.
   ptend_loc%lq(ixcldice) = .true.
   ptend_loc%lq(ixcldliq) = .true.

   if( microp_scheme .eq. 'MG' ) then
       ptend_loc%lq(ixnumliq) = .true.
       ptend_loc%lq(ixnumice) = .true.
#ifdef MODAL_AERO
       do m = 1, ntot_amode
          lnum = numptr_amode(m)
          if( lnum > 0 ) then
              ptend_loc%lq(lnum)= .true.
          endif
          do l = 1, nspec_amode(m)
             lmass = lmassptr_amode(l,m)
             ptend_loc%lq(lmass)= .true.
         enddo
      enddo
#endif
   end if

   if( microp_scheme .eq. 'RK' ) then

     ! Determine repartition heating from change in cloud ice.

       repartht(:ncol,:pver) = (latice/dtime) * ( state1%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver) )

     ! Non-micro and non-macrophysical external advective forcings to compute net condensation rate. 
     ! Note that advective forcing of condensate is aggregated into liquid phase.

       qtend(:ncol,:pver) = ( state1%q(:ncol,:pver,1) - qcwat(:ncol,:pver) ) * rdtime
       ttend(:ncol,:pver) = ( state1%t(:ncol,:pver)   - tcwat(:ncol,:pver) ) * rdtime
       ltend(:ncol,:pver) = ( totcw   (:ncol,:pver)   - lcwat(:ncol,:pver) ) * rdtime

     ! Compute Stratiform Macro-Microphysical Tendencies

       ! Add rain and snow fluxes as output variables from pcond, and into physics buffer
       rkflxprc  => pbuf(ls_flxprc_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
       rkflxsnw  => pbuf(ls_flxsnw_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)

       call t_startf('pcond')
       call pcond( lchnk, ncol,                                                &
                   state1%t, ttend, state1%q(1,1,1), qtend, state1%omega,      &
                   totcw, state1%pmid , state1%pdel, cld, fice, fsnow,         &
                   qme, prain, prodsnow, nevapr, evapsnow, evapheat, prfzheat, &
                   meltheat, prec_pcw, snow_pcw, dtime, fwaut,                 &
                   fsaut, fracw, fsacw, fsaci, ltend,                          &
                   rhdfda, rhu00, icefrac, state1%zi, ice2pr, liq2pr,          &
                   liq2snow, snowh, rkflxprc, rkflxsnw)
       call t_stopf('pcond')

   elseif( microp_scheme .eq. 'MG' ) then

     ! ------------------------------ !
     ! Liquid Stratiform Macrophysics !
     ! ------------------------------ !

       call t_startf('mmacro_pcond')

       zeros(:ncol,:pver)  = 0._r8
       qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
       qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
       nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
       ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)

     ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... ) 
     ! are separately provided into the prognostic stratiform macrophysics scheme. This is an
     ! attempt to resolve in-cloud and out-cloud forcings. 

       if( get_nstep() .le. 1 ) then
           tcwat(:ncol,:pver)   = state1%t(:ncol,:pver)
           qcwat(:ncol,:pver)   = state1%q(:ncol,:pver,1)
           lcwat(:ncol,:pver)   = qc(:ncol,:pver) + qi(:ncol,:pver)
           iccwat(:ncol,:pver)  = qi(:ncol,:pver)
           nlwat(:ncol,:pver)   = nc(:ncol,:pver)
           niwat(:ncol,:pver)   = ni(:ncol,:pver)
           ttend(:ncol,:pver)   = 0._r8
           qtend(:ncol,:pver)   = 0._r8
           ltend(:ncol,:pver)   = 0._r8
           itend(:ncol,:pver)   = 0._r8
           nltend(:ncol,:pver)  = 0._r8
           nitend(:ncol,:pver)  = 0._r8
           CC_T(:ncol,:pver)    = 0._r8
           CC_qv(:ncol,:pver)   = 0._r8
           CC_ql(:ncol,:pver)   = 0._r8
           CC_qi(:ncol,:pver)   = 0._r8
           CC_nl(:ncol,:pver)   = 0._r8
           CC_ni(:ncol,:pver)   = 0._r8
           CC_qlst(:ncol,:pver) = 0._r8
       else
           ttend(:ncol,:pver)   = ( state1%t(:ncol,:pver)   -  tcwat(:ncol,:pver)) * rdtime -   CC_T(:ncol,:pver) 
           qtend(:ncol,:pver)   = ( state1%q(:ncol,:pver,1) -  qcwat(:ncol,:pver)) * rdtime -  CC_qv(:ncol,:pver)
           ltend(:ncol,:pver)   = ( qc(:ncol,:pver) + qi(:ncol,:pver)  -  & 
                                         lcwat(:ncol,:pver) ) * rdtime - (CC_ql(:ncol,:pver) + CC_qi(:ncol,:pver))
           itend(:ncol,:pver)   = ( qi(:ncol,:pver)         - iccwat(:ncol,:pver)) * rdtime -  CC_qi(:ncol,:pver)
           nltend(:ncol,:pver)  = ( nc(:ncol,:pver)         -  nlwat(:ncol,:pver)) * rdtime -  CC_nl(:ncol,:pver)
           nitend(:ncol,:pver)  = ( ni(:ncol,:pver)         -  niwat(:ncol,:pver)) * rdtime -  CC_ni(:ncol,:pver)
       endif
       lmitend(:ncol,:pver) = ltend(:ncol,:pver) - itend(:ncol,:pver)

       t_inout(:ncol,:pver)  =  tcwat(:ncol,:pver) 
       qv_inout(:ncol,:pver) =  qcwat(:ncol,:pver)
       ql_inout(:ncol,:pver) =  lcwat(:ncol,:pver) - iccwat(:ncol,:pver)
       qi_inout(:ncol,:pver) = iccwat(:ncol,:pver)
       nl_inout(:ncol,:pver) =  nlwat(:ncol,:pver)
       ni_inout(:ncol,:pver) =  niwat(:ncol,:pver)

     ! Liquid Stratiform Macrophysics.
     ! The main roles of this subroutines are
     ! (1) compute net condensation rate of strayiform liquid ( cmeliq )
     ! (2) compute liquid stratus and ice stratus fractions. 
     ! Note 'ttend...' are advective tendencies except microphysical process while
     !      'CC...'    are microphysical tendencies. 

       call mmacro_pcond( lchnk, ncol, dtime, state1%pmid, state1%pdel,              &
                          t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout, &                  
                          ttend, qtend, lmitend, itend, nltend, nitend,              &
                          CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst,          & 
                          dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni,             &
                          concld_old, concld, landfrac, snowh,                       &
                          tlat, qvlat, qcten, qiten, ncten, niten,                   &
                          cmeliq, qvadj, qladj, qiadj, qllim, qilim,                 &
                          cld, alst, aist, qlst, qist ) 

     ! Compute net stratus fraction using maximum over-lapping assumption

       do k = 1, pver
       do i = 1, ncol
          ast(i,k) = max( alst(i,k), aist(i,k) )      
       enddo
       enddo

       call t_stopf('mmacro_pcond')

       do k = 1, pver
       do i = 1, ncol
          ptend_loc%s(i,k)          =  tlat(i,k)
          ptend_loc%q(i,k,1)        = qvlat(i,k)
          ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
          ptend_loc%q(i,k,ixcldice) = qiten(i,k)
          ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
          ptend_loc%q(i,k,ixnumice) = niten(i,k)
       end do
       end do   

       call outfld( 'MACPDT   ', tlat ,  pcols, lchnk )
       call outfld( 'MACPDQ   ', qvlat,  pcols, lchnk )
       call outfld( 'MACPDLIQ ', qcten,  pcols, lchnk )
       call outfld( 'MACPDICE ', qiten,  pcols, lchnk )
       call outfld( 'CLDVAPADJ', qvadj,  pcols, lchnk )
       call outfld( 'CLDLIQADJ', qladj,  pcols, lchnk )
       call outfld( 'CLDICEADJ', qiadj,  pcols, lchnk )
       call outfld( 'CLDLIQDET', dlf_ql, pcols, lchnk )
       call outfld( 'CLDICEDET', dlf_qi, pcols, lchnk )
       call outfld( 'CLDLIQLIM', qllim,  pcols, lchnk )
       call outfld( 'CLDICELIM', qilim,  pcols, lchnk )

     ! Here 'state_eq' is the equlibrium state after macrophysics for potential
     ! use in the radiation scheme later.

       call physics_ptend_sum( ptend_loc, ptend_all, state )
       call physics_update( state1, tend, ptend_loc, dtime )
       call physics_ptend_init( ptend_loc ) 
       call physics_state_copy( state1, state_eq )

     ! ----------------------- !
     ! Stratiform Microphysics !
     ! ----------------------- !

       ptend_loc%name         = 'cldmicro'
       ptend_loc%ls           = .true.
       ptend_loc%lq(1)        = .true.
       ptend_loc%lq(ixcldliq) = .true.
       ptend_loc%lq(ixcldice) = .true.
       ptend_loc%lq(ixnumliq) = .true.
       ptend_loc%lq(ixnumice) = .true.

#ifndef MODAL_AERO
       call rad_cnst_get_clim_info( naero = naer_all )
       allocate( aer_mmr( pcols, pver, naer_all ) )
       do m = 1, naer_all
          call rad_cnst_get_clim_aer( m, state1, pbuf, aermmr1 )
          aer_mmr(:ncol,:,m) = aermmr1(:ncol,:)
       enddo
#endif

#ifdef MODAL_AERO
       do m = 1, ntot_amode
          lnum = numptr_amode(m)
          if( lnum > 0 ) then
              ptend_loc%lq(lnum)= .true.
          endif
          do l = 1, nspec_amode(m)
             lmass = lmassptr_amode(l,m)
             ptend_loc%lq(lmass)= .true.
          enddo
       enddo
#endif

      call t_startf('mmicro_pcond')

      qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
      qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
      nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
      ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)

      if( micro_treatment .eq. 'inter' ) then
          alst_mic(:ncol,:pver) = ast(:ncol,:pver)
          aist_mic(:ncol,:pver) = ast(:ncol,:pver)
      elseif( micro_treatment .eq. 'compl' ) then
          alst_mic(:ncol,:pver) = alst(:ncol,:pver)
          aist_mic(:ncol,:pver) = aist(:ncol,:pver)
      endif

      ! calculate aerosol activation (naai for ice, npccn for liquid) and 
      ! dust size (rndst) and number (nacon) for contact nucleation

      call microp_aero_ts ( lchnk, ncol, dtime, state1%t, zeros,                       &
                         state1%q(1,1,1), qc, qi,                     &
                         nc, ni, state1%pmid, state1%pdel, ast,                     &
                         alst_mic, aist_mic,                                        &
	                 cldo, state1%pint, state1%rpdel, state1%zm, state1%omega,  &
#ifdef MODAL_AERO
                         state1%q, cflx, ptend_loc%q, qqcw, dgnumwet, dgnum,        &
#else
                         aer_mmr,                                                   &
#endif
                         kkvh, tke, turbtype, smaw, wsub, wsubi,                    &
                         naai, npccn, rndst,nacon)

      ! Call MG Microphysics

      call mmicro_pcond( lchnk, ncol, dtime, state1%t, zeros,                       &
                         state1%q(1,1,1), zeros, zeros, qc, qi,                     &
                         nc, ni, state1%pmid, state1%pdel, ast,                     &
                         alst_mic, aist_mic,                                        &
	                 cldo, state1%pint, state1%rpdel, state1%zm, state1%omega,  &
                          rate1cld,                                                 & 
                         naai, npccn, rndst,nacon,                                  &
                         rhdfda, rhu00, fice, tlat, qvlat,                          &
                         qcten, qiten, ncten, niten, effc,                          &
                         effc_fn, effi, prect, preci,                               & 
                         nevapr, evapsnow,                                          &
                         prain, prodsnow, cmeice, dei, mu,                          &
                         lambdac, qsout, des,                                       &   
                         qcsevap, qisevap, qvres, cmeiout,                          &
                         vtrmc, vtrmi, qcsedten, qisedten,                          &
                         prao, prco, mnuccco, mnuccto, msacwio, psacwso,            &
                         bergso, bergo, melto, homoo, qcreso, prcio, praio, qireso, &
                         mnuccro, pracso, meltsdt, frzrdt )

    ! Reassign rate1 if modal aerosols
#ifdef MODAL_AERO
      rate1ord_cw2pr_st(1:ncol,1:pver)=rate1cld(1:ncol,1:pver)
#endif
    ! Sedimentation velocity for liquid stratus cloud droplet

      wsedl(:ncol,:pver) = vtrmc(:ncol,:pver)

    ! Nominal values for no stratiform (convective only) cloud.
    ! Convert snow mixing ratio to microns

      do k = 1, pver
      do i = 1, ncol 
         des(i,k) = des(i,k) * 1.e6_r8
         if( ast(i,k) .lt. 1.e-4_r8 ) then
             mu(i,k) = mucon
             lambdac(i,k) = (mucon + 1._r8)/dcon
             dei(i,k) = deicon
         endif
      end do
      end do 

    ! Microphysical tendencies for use in the macrophysics at the next time step

      CC_T(:ncol,:pver)    =  tlat(:ncol,:pver)/cpair
      CC_qv(:ncol,:pver)   = qvlat(:ncol,:pver)
      CC_ql(:ncol,:pver)   = qcten(:ncol,:pver)
      CC_qi(:ncol,:pver)   = qiten(:ncol,:pver)
      CC_nl(:ncol,:pver)   = ncten(:ncol,:pver)
      CC_ni(:ncol,:pver)   = niten(:ncol,:pver)
      CC_qlst(:ncol,:pver) = qcten(:ncol,:pver)/max(0.01_r8,alst_mic(:ncol,:pver))

    ! Net stratiform condensation rate

      qme(:ncol,:pver) = cmeliq(:ncol,:pver) + cmeiout(:ncol,:pver) 

      call t_stopf('mmicro_pcond')

#ifndef MODAL_AERO
      deallocate(aer_mmr) 
#endif

   end if

   if( microp_scheme .eq. 'RK' ) then

       do k = 1, pver
       do i = 1, ncol
          ptend_loc%s(i,k)          =   qme(i,k)*( latvap + latice*fice(i,k) ) + &
                                        evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)
          ptend_loc%q(i,k,1)        = - qme(i,k) + nevapr(i,k)
          ptend_loc%q(i,k,ixcldice) =   qme(i,k)*fice(i,k)         - ice2pr(i,k)
          ptend_loc%q(i,k,ixcldliq) =   qme(i,k)*(1._r8-fice(i,k)) - liq2pr(i,k)
       end do
       end do
 
       do k = 1, pver
       do i = 1, ncol
          aist(i,k)  = cld(i,k)
          alst(i,k)  = cld(i,k)
          ast(i,k)   = cld(i,k)
          icimr(i,k) = (state1%q(i,k,ixcldice) + dtime*ptend_loc%q(i,k,ixcldice)) / max(0.01_r8,aist(i,k))
          icwmr(i,k) = (state1%q(i,k,ixcldliq) + dtime*ptend_loc%q(i,k,ixcldliq)) / max(0.01_r8,alst(i,k))
       end do
       end do

    ! Convert precipitation from [ kg/m2 ] to [ m/s ]
      snow_pcw(:ncol) = snow_pcw(:ncol)/1000._r8
      prec_pcw(:ncol) = prec_pcw(:ncol)/1000._r8

      do k = 1, pver
      do i = 1, ncol
         cmeheat(i,k) = qme(i,k) * ( latvap + latice*fice(i,k) )
         cmeice (i,k) = qme(i,k) *   fice(i,k)
         cmeliq (i,k) = qme(i,k) * ( 1._r8 - fice(i,k) )
      end do
      end do

    ! Record history variables

      call outfld( 'FWAUT'   , fwaut,       pcols, lchnk )
      call outfld( 'FSAUT'   , fsaut,       pcols, lchnk )
      call outfld( 'FRACW'   , fracw,       pcols, lchnk )
      call outfld( 'FSACW'   , fsacw,       pcols, lchnk )
      call outfld( 'FSACI'   , fsaci,       pcols, lchnk )

      call outfld( 'PCSNOW'  , snow_pcw,    pcols, lchnk )
      call outfld( 'FICE'    , fice,        pcols, lchnk )
      call outfld( 'CMEICE'  , cmeice,      pcols, lchnk )
      call outfld( 'CMELIQ'  , cmeliq,      pcols, lchnk )
      call outfld( 'ICE2PR'  , ice2pr,      pcols, lchnk )
      call outfld( 'LIQ2PR'  , liq2pr,      pcols, lchnk )
      call outfld( 'HPROGCLD', ptend_loc%s, pcols, lchnk )
      call outfld( 'HEVAP   ', evapheat,    pcols, lchnk )
      call outfld( 'HMELT'   , meltheat,    pcols, lchnk )
      call outfld( 'HCME'    , cmeheat ,    pcols, lchnk )
      call outfld( 'HFREEZ'  , prfzheat,    pcols, lchnk )
      call outfld( 'HREPART' , repartht,    pcols, lchnk )
      call outfld('LS_FLXPRC', rkflxprc,    pcols, lchnk )
      call outfld('LS_FLXSNW', rkflxsnw,    pcols, lchnk )

   elseif( microp_scheme .eq. 'MG' ) then

      do k = 1, pver
      do i = 1, ncol
         ptend_loc%s(i,k)          =  tlat(i,k)
         ptend_loc%q(i,k,1)        = qvlat(i,k)
         ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
         ptend_loc%q(i,k,ixcldice) = qiten(i,k)
         ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
         ptend_loc%q(i,k,ixnumice) = niten(i,k)
      enddo
      enddo
    
    ! For precip, accumulate only total precip in prec_pwc and snow_pwc variables.
    ! Other precip output varirables are set to 0
      prec_pcw(:ncol) = prect(:ncol)
      snow_pcw(:ncol) = preci(:ncol)
      prec_sed(:ncol) = 0._r8
      snow_sed(:ncol) = 0._r8
      prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol) - rliq(:ncol)
      snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol) 

   endif

   ! ------------------------------- !
   ! Update microphysical tendencies !
   ! ------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   ptend_all%name = 'stratiform'
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc ) 

   ! Below 'cldfrc' block should not be performed for MG since 
   ! 'cld' has already been computed from the mmacro_pcond.

   if ( microp_scheme .ne. 'MG' ) then

   call t_startf("cldfrc")
   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld, rhcloud, clc, state1%pdel,                                    &
                cmfmc, cmfmc2, landfrac, snowh, concld, cldst,                     &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                relhum, 0 )    
   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld2, rhcloud2, clc, state1%pdel,                                  &
                cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                relhum2, 1  )              

   call t_stopf("cldfrc")

   do k = 1, pver
   do i = 1, ncol
      if( relhum(i,k) < rhu00(i,k) ) then
          rhdfda(i,k)=0.0_r8
      elseif( relhum(i,k) >= 1.0_r8 ) then
          rhdfda(i,k)=0.0_r8
      else
        ! Under certain circumstances, rh+ causes cld not to be changed
        ! when at an upper limit, or w/ strong subsidence
          if( ( rhcloud2(i,k) - rhcloud(i,k) ) < 1.e-4_r8 ) then
                rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8 
          else
                rhdfda(i,k) = 0.01_r8*relhum(i,k)/(rhcloud2(i,k)-rhcloud(i,k))
          endif
      endif
   enddo
   enddo

   endif

 ! Copy of concld/fice to put in physics buffer
 ! Below are used only for convective cloud.

   concld_ql(:ncol,:pver) = concld(:ncol,:pver)
   fice_ql(:ncol,:pver)   = fice(:ncol,:pver)

   if( micro_treatment .eq. 'inter' ) then
       icecldf(:ncol,:pver) = ast(:ncol,:pver)
       liqcldf(:ncol,:pver) = ast(:ncol,:pver)
   elseif( micro_treatment .eq. 'compl' ) then
       icecldf(:ncol,:pver) = aist(:ncol,:pver)
       liqcldf(:ncol,:pver) = alst(:ncol,:pver)
   endif

   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )
   call outfld( 'CNVCLD  ', clc,    pcols, lchnk )

   call outfld( 'ICECLDF ', aist,   pcols, lchnk )
   call outfld( 'LIQCLDF ', alst,   pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )   

   if( microp_scheme .eq. 'MG' ) then

     ! ------------------------------------------------- !
     ! Save equilibrium state variables for macrophysics !        
     ! at the next time step                             !
     ! ------------------------------------------------- !
       do k = 1, pver
          tcwat(:ncol,k)  = state_eq%t(:ncol,k) 
          qcwat(:ncol,k)  = state_eq%q(:ncol,k,1)
          lcwat(:ncol,k)  = state_eq%q(:ncol,k,ixcldliq) + state_eq%q(:ncol,k,ixcldice)
          iccwat(:ncol,k) = state_eq%q(:ncol,k,ixcldice)
          nlwat(:ncol,k)  = state_eq%q(:ncol,k,ixnumliq)
          niwat(:ncol,k)  = state_eq%q(:ncol,k,ixnumice)         
       end do
     ! Effective droplet radius
       rel(:ncol,:pver)        = effc(:ncol,:pver)
       rel_fn(:ncol,:pver)     = effc_fn(:ncol,:pver)	
       rei(:ncol,:pver)        = effi(:ncol,:pver)
       rel2(:ncol,:pver)       = rel(:ncol,:pver) * 0.9071_r8 ! Convect to effective volume radius assuming pgam = 8
       rei2(:ncol,:pver)       = rei(:ncol,:pver) * 0.6057_r8 ! Convect to effective volume radius at pgam = 0 for ice 
     ! ----------------------------------------------------------- ! 
     ! Adjust in-cloud water values to take account of convective  !
     ! in-cloud water. It is used to calculate the values of       !
     ! icwlp and iciwp to pass to the radiation.                   ! 
     ! ----------------------------------------------------------- !
       allcld_ice(:ncol,:pver) = 0._r8 ! Grid-avg all cloud liquid
       allcld_liq(:ncol,:pver) = 0._r8 ! Grid-avg all cloud ice
       if( conv_water_in_rad /= 0 ) then
	   call conv_water_4rad( lchnk, ncol, pbuf, conv_water_in_rad, rei, state1%pdel, &
	                         state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice),         &
                                 allcld_liq, allcld_ice )
	else
	   allcld_liq(:ncol,:) = state1%q(:ncol,:,ixcldliq)  ! Grid-ave all cloud liquid
	   allcld_ice(:ncol,:) = state1%q(:ncol,:,ixcldice)  !           "        ice 
	end if
      ! ------------------------------------------------------------ !
      ! Compute in cloud ice and liquid mixing ratios                !
      ! Note that 'iclwp, iciwp' are used for radiation computation. !
      ! ------------------------------------------------------------ !
        do k = 1, pver
        do i = 1, ncol
            ! Limits for in-cloud mixing ratios consistent with MG microphysics
            ! in-cloud mixing ratio 0.0001 to 0.005 kg/kg
              icimr(i,k)     = min( allcld_ice(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
              icwmr(i,k)     = min( allcld_liq(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
              icimrst(i,k)   = min( state1%q(i,k,ixcldice) / max(0.0001_r8,icecldf(i,k)),0.005_r8 )
              icwmrst(i,k)   = min( state1%q(i,k,ixcldliq) / max(0.0001_r8,liqcldf(i,k)),0.005_r8 )
              icinc(i,k)     = state1%q(i,k,ixnumice) / max(0.0001_r8,icecldf(i,k)) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
              icwnc(i,k)     = state1%q(i,k,ixnumliq) / max(0.0001_r8,liqcldf(i,k)) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
 	      iwc(i,k)       = allcld_ice(i,k) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
	      lwc(i,k)       = allcld_liq(i,k) * state1%pmid(i,k) / (287.15_r8*state1%t(i,k))
              effliq(i,k)    = effc(i,k)
              effliq_fn(i,k) = effc_fn(i,k)
              effice(i,k)    = effi(i,k)
            ! Calculate total cloud water paths in each layer
              iciwp(i,k)     = icimr(i,k) * state1%pdel(i,k) / gravit
              iclwp(i,k)     = icwmr(i,k) * state1%pdel(i,k) / gravit
            ! Calculate stratiform cloud water paths in each layer
            ! Note: uses stratiform cloud fraction!
              iciwpst(i,k)   = min(state1%q(i,k,ixcldice)/max(0.0001_r8,ast(i,k)),0.005_r8) * state1%pdel(i,k) / gravit
              iclwpst(i,k)   = min(state1%q(i,k,ixcldliq)/max(0.0001_r8,ast(i,k)),0.005_r8) * state1%pdel(i,k) / gravit
            ! Calculate convective in-cloud LWP.
              iclwpconv(i,k) = max(allcld_liq(i,k) - state%q(i,k,ixcldliq),0._r8)/max(0.0001_r8,concld(i,k)) 
              iciwpconv(i,k) = max(allcld_ice(i,k) - state%q(i,k,ixcldice),0._r8)/max(0.0001_r8,concld(i,k)) 
            ! ------------------------------ !
            ! Adjust cloud fraction for snow !
            ! ------------------------------ !
              cldfsnow(i,k) = cld(i,k)
            ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
              if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. & 
                  ( concld(i,k)   .lt. 1.e-4_r8 ) .and. & 
                  ( state1%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
                    cldfsnow(i,k) = 0._r8
              endif
            ! If no cloud and snow, then set to 0.25
              if( ( cldfsnow(i,k) .lt. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then 
                    cldfsnow(i,k) = 0.25_r8
              endif
            ! Calculate in-cloud snow water path
              icswp(i,k) = qsout(i,k) / max( 0.0001_r8, cldfsnow(i,k) ) * state1%pdel(i,k) / gravit
        enddo 
        enddo

      ! --------------------- !
      ! History Output Fields !
      ! --------------------- !

      ! Column droplet concentration

        do i = 1, ncol
           cdnumc(i) = 0._r8
           do k = 1, pver
              cdnumc(i) = cdnumc(i) + state1%q(i,k,ixnumliq)*state1%pdel(i,k)/gravit
           end do
        end do

      ! Averaging for new output fields

        efcout(:,:)      = 0._r8
        efiout(:,:)      = 0._r8
        ncout(:,:)       = 0._r8
        niout(:,:)       = 0._r8	
        freql(:,:)       = 0._r8
        freqi(:,:)       = 0._r8
        liqcldf_out(:,:) = 0._r8
        icecldf_out(:,:) = 0._r8
        icwmrst_out(:,:) = 0._r8
        icimrst_out(:,:) = 0._r8
        do k = 1, pver
        do i = 1, ncol
           if( liqcldf(i,k) .gt. 0.01_r8 .and. icwmrst(i,k) .gt. 5.e-5_r8 ) then
               efcout(i,k) = effc(i,k)
               ncout(i,k)  = icwnc(i,k)
               freql(i,k)  = 1._r8
               liqcldf_out(i,k) = liqcldf(i,k)
               icwmrst_out(i,k) = icwmrst(i,k)
           endif
           if( icecldf(i,k) .gt. 0.01_r8 .and. icimrst(i,k) .gt. 1.e-6_r8 ) then
               efiout(i,k) = effi(i,k)
               niout(i,k)  = icinc(i,k)
               freqi(i,k)  = 1._r8
               icecldf_out(i,k) = icecldf(i,k)
               icimrst_out(i,k) = icimrst(i,k)
            endif
        end do
        end do

        call outfld( 'AREL' , efcout,  pcols, lchnk )
        call outfld( 'AREI' , efiout,  pcols, lchnk )
        call outfld( 'AWNC' , ncout,   pcols, lchnk )
        call outfld( 'AWNI' , niout,   pcols, lchnk )
        call outfld( 'FREQL', freql,   pcols, lchnk )
        call outfld( 'FREQI', freqi,   pcols, lchnk )

      ! Cloud top effective radius and number.

        fcti(:)  = 0._r8
        fctl(:)  = 0._r8
        ctrel(:) = 0._r8
        ctrei(:) = 0._r8
        ctnl(:)  = 0._r8
        ctni(:)  = 0._r8
        do i = 1, ncol
        do k = 1, pver
           if( liqcldf(i,k) .gt. 0.01_r8 .and. icwmrst(i,k) .gt. 1.e-7_r8 ) then
               ctrel(i) = effc(i,k)
               ctnl(i)  = icwnc(i,k)
               fctl(i)  = 1._r8
               exit
           endif
           if( icecldf(i,k) .gt. 0.01_r8 .and. icimrst(i,k) .gt. 1.e-7_r8 ) then
               ctrei(i) = effi(i,k)
               ctni(i)  = icinc(i,k)
               fcti(i)  = 1._r8
               exit
           endif
        enddo
        enddo

        call outfld( 'ACTREL'     , ctrel,     pcols, lchnk )
        call outfld( 'ACTREI'     , ctrei,     pcols, lchnk )
        call outfld( 'ACTNL'      , ctnl,      pcols, lchnk )
        call outfld( 'ACTNI'      , ctni,      pcols, lchnk )
        call outfld( 'FCTL'       , fctl,      pcols, lchnk )
        call outfld( 'FCTI'       , fcti,      pcols, lchnk )

        call outfld( 'MPDT'       , tlat,      pcols, lchnk )
        call outfld( 'MPDQ'       , qvlat,     pcols, lchnk )
        call outfld( 'MPDLIQ'     , qcten,     pcols, lchnk )
        call outfld( 'MPDICE'     , qiten,     pcols, lchnk )
        call outfld( 'ICINC'      , icinc,     pcols, lchnk )
        call outfld( 'ICWNC'      , icwnc,     pcols, lchnk )
        call outfld( 'EFFLIQ'     , effliq,    pcols, lchnk )
        call outfld( 'EFFLIQ_IND' , effliq_fn, pcols, lchnk )
        call outfld( 'EFFICE'     , effice,    pcols, lchnk )
        call outfld( 'WSUB'       , wsub,      pcols, lchnk )
        call outfld( 'WSUBI'      , wsubi,     pcols, lchnk )
        call outfld( 'CDNUMC'     , cdnumc,    pcols, lchnk )

   elseif( microp_scheme .eq. 'RK' ) then

     do k = 1, pver
     do i = 1, ncol
	iwc(i,k)   = state1%q(i,k,ixcldice)*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
	lwc(i,k)   = state1%q(i,k,ixcldliq)*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
        icimr(i,k) = state1%q(i,k,ixcldice) / max(0.01_r8,rhcloud(i,k))
        icwmr(i,k) = state1%q(i,k,ixcldliq) / max(0.01_r8,rhcloud(i,k))
     end do
     end do

   endif ! RK,MG microphysics

   ! --------------------------------------------- !
   ! Common outfield calls for either microphysics !
   ! --------------------------------------------- !

   call outfld( 'IWC'      , iwc,         pcols, lchnk )
   call outfld( 'LWC'      , lwc,         pcols, lchnk )
   call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
   call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )
   call outfld( 'ICIMRST'  , icimrst_out, pcols, lchnk )
   call outfld( 'ICWMRST'  , icwmrst_out, pcols, lchnk )
   call outfld( 'CME'      , qme,         pcols, lchnk )
   call outfld( 'PRODPREC' , prain,       pcols, lchnk )
   call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
   call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )
   call outfld( 'QCSEVAP'  , qcsevap,     pcols, lchnk )
   call outfld( 'QISEVAP'  , qisevap,     pcols, lchnk )
   call outfld( 'QVRES'    , qvres,       pcols, lchnk )
   call outfld( 'CMEIOUT'  , cmeiout,     pcols, lchnk )
   call outfld( 'CMELIQ'   , cmeliq,      pcols, lchnk )
   call outfld( 'VTRMC'    , vtrmc,       pcols, lchnk )
   call outfld( 'VTRMI'    , vtrmi,       pcols, lchnk )
   call outfld( 'QCSEDTEN' , qcsedten,    pcols, lchnk )
   call outfld( 'QISEDTEN' , qisedten,    pcols, lchnk )
   call outfld( 'PRAO'     , prao,        pcols, lchnk )
   call outfld( 'PRCO'     , prco,        pcols, lchnk )
   call outfld( 'MNUCCCO'  , mnuccco,     pcols, lchnk )
   call outfld( 'MNUCCTO'  , mnuccto,     pcols, lchnk )
   call outfld( 'MSACWIO'  , msacwio,     pcols, lchnk )
   call outfld( 'PSACWSO'  , psacwso,     pcols, lchnk )
   call outfld( 'BERGSO'   , bergso,      pcols, lchnk )
   call outfld( 'BERGO'    , bergo,       pcols, lchnk )
   call outfld( 'MELTO'    , melto,       pcols, lchnk )
   call outfld( 'HOMOO'    , homoo,       pcols, lchnk )
   call outfld( 'QCRESO'   , qcreso,      pcols, lchnk )
   call outfld( 'PRCIO'    , prcio,       pcols, lchnk )
   call outfld( 'PRAIO'    , praio,       pcols, lchnk )
   call outfld( 'QIRESO'   , qireso,      pcols, lchnk )
   call outfld( 'MNUCCRO'  , mnuccro,     pcols, lchnk )
   call outfld( 'PRACSO'   , pracso ,     pcols, lchnk )
   call outfld( 'MELTSDT'  , meltsdt,     pcols, lchnk )
   call outfld( 'FRZRDT'   , frzrdt ,     pcols, lchnk )

   if( microp_scheme .eq. 'MG' ) then
       ftem(:ncol,:pver) =  qcreso(:ncol,:pver)
       call outfld( 'MPDW2V', ftem, pcols, lchnk )
       ftem(:ncol,:pver) =  melto(:ncol,:pver) - mnuccco(:ncol,:pver) - mnuccto(:ncol,:pver) - &
                            bergo(:ncol,:pver) - homoo  (:ncol,:pver) - msacwio(:ncol,:pver)
       call outfld( 'MPDW2I', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -prao(:ncol,:pver) - prco(:ncol,:pver) - psacwso(:ncol,:pver) - &
                            bergso(:ncol,:pver)
       call outfld( 'MPDW2P', ftem, pcols, lchnk )
       ftem(:ncol,:pver) =  cmeiout(:ncol,:pver) + qireso (:ncol,:pver)
       call outfld( 'MPDI2V', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -melto(:ncol,:pver) + mnuccco(:ncol,:pver) + mnuccto(:ncol,:pver) + &
                            bergo(:ncol,:pver) + homoo  (:ncol,:pver) + msacwio(:ncol,:pver)
       call outfld( 'MPDI2W', ftem, pcols, lchnk )
       ftem(:ncol,:pver) = -prcio(:ncol,:pver) - praio  (:ncol,:pver)
       call outfld( 'MPDI2P', ftem, pcols, lchnk )
   endif

   call t_stopf('stratiform_microphys')

   if( microp_scheme .eq. 'RK' ) then

      prec_str(:ncol) = prec_str(:ncol) + prec_pcw(:ncol)
      snow_str(:ncol) = snow_str(:ncol) + snow_pcw(:ncol)

      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)

    ! Save variables for use in the macrophysics at the next time step

      do k = 1, pver
         qcwat(:ncol,k) = state1%q(:ncol,k,1)
         tcwat(:ncol,k) = state1%t(:ncol,k)
         lcwat(:ncol,k) = state1%q(:ncol,k,ixcldice) + state1%q(:ncol,k,ixcldliq)
      end do
  
    ! Cloud water and ice particle sizes, saved in physics buffer for radiation

      call cldefr( lchnk, ncol, landfrac, state1%t, rel, rei, state1%ps, state1%pmid, landm, icefrac, snowh )
      rel2(:ncol,:pver) = rel(:ncol,:pver)      
      rei2(:ncol,:pver) = rei(:ncol,:pver)      

   end if

   end subroutine stratiform_tend

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_1(state1,ptend,i,k, &
        dtime,qme,fice,snow_pcw,prec_pcw, &
        prain,nevapr,prodsnow, evapsnow, &
        ice2pr,liq2pr,liq2snow)

     use physics_types, only: physics_state, physics_ptend
     use physconst,     only: tmelt

     implicit none
     
     integer, intent(in) :: i,k
     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     type(physics_ptend), intent(in) :: ptend  ! local copy of the ptend variable
     real(r8), intent(in)  :: dtime                ! timestep
     real(r8), intent(in) :: qme(pcols,pver)          ! local condensation - evaporation of cloud water

     real(r8), intent(in) :: prain(pcols,pver)          ! local production of precipitation
     real(r8), intent(in) :: nevapr(pcols,pver)          ! local evaporation of precipitation
     real(r8), intent(in) :: prodsnow(pcols,pver)          ! local production of snow
     real(r8), intent(in) :: evapsnow(pcols,pver)          ! local evaporation of snow
     real(r8), intent(in) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
     real(r8), intent(in) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
     real(r8), intent(in) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
     real(r8), intent(in) :: fice    (pcols,pver)          ! Fractional ice content within cloud
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: prec_pcw(pcols)

     real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
     real(r8) w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

     pr1 = 0
     hs1 = 0
     qv1 = 0
     ql1 = 0
     qi1 = 0
     qs1 = 0
     qr1 = 0
     w1 = 0
     wl = 0
     wv = 0
     wi = 0
     wlf = 0
     wvf = 0 
     wif = 0


     write(iulog,*) 
     write(iulog,*) ' input state, t, q, l, i ', k, state1%t(i,k), state1%q(i,k,1), state1%q(i,k,ixcldliq),  state1%q(i,k,ixcldice)
     write(iulog,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip before accumulation                      ', k, pr1

     wv = wv + state1%q(i,k,1       )*state1%pdel(i,k)/gravit
     wl = wl + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     wi = wi + state1%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

     qvf = state1%q(i,k,1) + ptend%q(i,k,1)*dtime
     qlf = state1%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dtime
     qif = state1%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dtime

     if (qvf.lt.0._r8) then
        write(iulog,*) ' qvf is negative *******', qvf
     endif
     if (qlf.lt.0._r8) then
        write(iulog,*) ' qlf is negative *******', qlf
     endif
     if (qif.lt.0._r8) then
        write(iulog,*) ' qif is negative *******', qif
     endif
     write(iulog,*) ' qvf, qlf, qif ', qvf, qlf, qif

     wvf = wvf + qvf*state1%pdel(i,k)/gravit
     wlf = wlf + qlf*state1%pdel(i,k)/gravit
     wif = wif + qif*state1%pdel(i,k)/gravit

     hs1 = hs1 + ptend%s(i,k)*state1%pdel(i,k)/gravit
     pr1 = pr1 + state1%pdel(i,k)/gravit*(prain(i,k)-nevapr(i,k))
     qv1 = qv1 - (qme(i,k)-nevapr(i,k))*state1%pdel(i,k)/gravit    ! vdot
     w1  = w1  + (qme(i,k)-prain(i,k))*state1%pdel(i,k)/gravit    ! cdot
     qi1 = qi1 + ((qme(i,k))*fice(i,k)        -ice2pr(i,k) )*state1%pdel(i,k)/gravit   ! idot
     ql1 = ql1 + ((qme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state1%pdel(i,k)/gravit   ! ldot

     qr1 = qr1 &
          + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
          -(nevapr(i,k)-evapsnow(i,k)) &     ! rain evaporation
          )*state1%pdel(i,k)/gravit
     qs1 = qs1 &
          + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
          -evapsnow(i,k)               &     ! snow evaporation
          )*state1%pdel(i,k)/gravit

     if (state1%t(i,k).gt.tmelt) then
        qr1 = qr1 + qs1
        qs1 = 0._r8
     endif
     write(iulog,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip after accumulation      ', k, pr1
     write(iulog,*)
     write(iulog,*) ' layer prain, nevapr, pdel ', prain(i,k), nevapr(i,k), state1%pdel(i,k)
     write(iulog,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
     write(iulog,*) ' layer prain-prodsnow, liq2pr-liq2snow ', prain(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
     write(iulog,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), nevapr(i,k)-evapsnow(i,k)
     write(iulog,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
     write(iulog,*) ' layer ice2pr+liq2pr, prain ', ice2pr(i,k)+liq2pr(i,k), prain(i,k)
     write(iulog,*)
     write(iulog,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
     write(iulog,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
     write(iulog,*) ' condensate produce after accum                   ', k, w1
     write(iulog,*) ' liq+ice tends accum                              ', k, ql1+qi1
     write(iulog,*) ' change in total water after accum                ', k, qv1+ql1+qi1
     write(iulog,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
     write(iulog,*) ' fice at this lev ', fice(i,k)
     write(iulog,*)

     res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36_r8))
     write(iulog,*) ' relative residual in column method 1             ', k, res

     write(iulog,*) ' relative residual in column method 2             ',&
	 k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36_r8))
     !            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
     if (res.gt.1.e-14_r8) then
        call endrun ('STRATIFORM_TEND')
     endif

     !             w3  = qme(i,k) * (latvap + latice*fice(i,k)) &
     !               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

     res = qs1+qr1-pr1
     w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
     if (w4.gt.0._r8)  then
        if (res/w4.gt.1.e-14_r8) then
           write(iulog,*) ' imbalance in precips calculated two ways '
           write(iulog,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                res/w4, pr1, qr1, qs1, qr1+qs1
           !                   call endrun()
        endif
     endif
     if (k.eq.pver) then
        write(iulog,*) ' pcond returned precip, rain and snow rates ', prec_pcw(i), prec_pcw(i)-snow_pcw(i), snow_pcw(i)
        write(iulog,*) ' I calculate ', pr1, qr1, qs1
        !               call endrun
        write(iulog,*) ' byrons water check ', wv+wl+wi-pr1*dtime, wvf+wlf+wif
     endif
     write(iulog,*)


   end subroutine debug_microphys_1

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_2(state1,&
        snow_pcw,fsaut,fsacw ,fsaci, meltheat)

     use ppgrid,        only: pver
     use physconst,     only: tmelt
     use physics_types, only: physics_state
     
     implicit none

     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: fsaut(pcols,pver)              
     real(r8), intent(in) :: fsacw(pcols,pver)              
     real(r8), intent(in) :: fsaci(pcols,pver)              
     real(r8), intent(in) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip


     integer  i,ncol,lchnk


     ncol = state1%ncol
     lchnk = state1%lchnk
     
     do i = 1,ncol
        if (snow_pcw(i) .gt. 0.01_r8/8.64e4_r8  .and.  state1%t(i,pver) .gt. tmelt) then
           write(iulog,*) ' stratiform: snow, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
        
        if (snow_pcw(i)*8.64e4_r8 .lt. -1.e-5_r8) then
           write(iulog,*) ' neg snow ', snow_pcw(i)*8.64e4_r8
           write(iulog,*) ' stratiform: snow_pcw, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
     end do
     
   end subroutine debug_microphys_2

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine conv_water_4rad( lchnk, ncol, pbuf, conv_water_mode, &
                               rei, pdel, ls_liq, ls_ice, totg_liq, totg_ice )

   ! --------------------------------------------------------------------- ! 
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       ! 
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !                                                                       !
   !---------------------------------------------------------------------- !

   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx 
   use cam_history,     only: outfld
   use phys_control,    only: phys_getopts
   use phys_debug_util, only: phys_debug_col
   
   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !

   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf
   integer,  intent(in) :: lchnk
   integer,  intent(in) :: ncol
   integer,  intent(in) :: conv_water_mode
   real(r8), intent(in) :: rei(pcols,pver)        ! Ice effective drop size (microns)
   real(r8), intent(in) :: pdel(pcols,pver)       ! Moist pressure difference across layer
   real(r8), intent(in) :: ls_liq(pcols,pver)     ! Large-scale contributions to GBA cloud liq      
   real(r8), intent(in) :: ls_ice(pcols,pver)     ! Large-scale contributions to GBA cloud ice 

   real(r8), intent(out):: totg_ice(pcols,pver)   ! Total GBA in-cloud ice
   real(r8), intent(out):: totg_liq(pcols,pver)   ! Total GBA in-cloud liquid

   ! --------------- !
   ! Local Workspace !
   ! --------------- !

   ! Physics buffer fields

   real(r8), pointer, dimension(:,:) ::  ast      ! Physical liquid+ice stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  cu_frac  ! Final convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction

   real(r8), pointer, dimension(:,:) ::  alst     ! Physical liquid stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  aist     ! Physical ice    stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  qlst     ! Physical in-stratus LWC [kg/kg]
   real(r8), pointer, dimension(:,:) ::  qist     ! Physical in-stratus IWC [kg/kg]

   real(r8), pointer, dimension(:,:) ::  dp_icwmr ! Deep conv. cloud water
   real(r8), pointer, dimension(:,:) ::  sh_icwmr ! Shallow conv. cloud water
   real(r8), pointer, dimension(:,:) ::  fice     ! Ice partitioning ratio

   ! Local Variables

   real(r8) :: conv_ice(pcols,pver)               ! Convective contributions to IC cloud ice
   real(r8) :: conv_liq(pcols,pver)               ! Convective contributions to IC cloud liquid
   real(r8) :: tot_ice(pcols,pver)                ! Total IC ice
   real(r8) :: tot_liq(pcols,pver)                ! Total IC liquid

   integer  :: i,k,itim,ifld                      ! Lon, lev indices buff stuff.
   real(r8) :: cu_icwmr                           ! Convective  water for this grid-box.   
   real(r8) :: ls_icwmr                           ! Large-scale water for this grid-box. 
   real(r8) :: tot_icwmr                          ! Large-scale water for this grid-box.  
   real(r8) :: ls_frac                            ! Large-scale cloud frac for this grid-box. 
   real(r8) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac 
   real(r8) :: kabs, kabsi, kabsl, alpha, dp0, sh0, ic_limit, frac_limit  
   real(r8) :: wrk1         

   ! --------- !
   ! Parameter !
   ! --------- !

   parameter( kabsl = 0.090361_r8, frac_limit = 0.01_r8, ic_limit = 1.e-12_r8 )

 ! Get microphysics option

   character(len=16) :: microp_scheme 
   call phys_getopts( microp_scheme_out = microp_scheme )

 ! Get convective in-cloud water and ice/water temperature partitioning.

   ifld = pbuf_get_fld_idx('ICWMRSH')
   sh_icwmr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('ICWMRDP')
   dp_icwmr => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('FICE')
   fice => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

 ! Get convective in-cloud fraction    

   ifld = pbuf_get_fld_idx('SH_FRAC')
   sh_frac => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('DP_FRAC')
   dp_frac => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('CONCLDQL')
   cu_frac => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)  

   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('AST')
   ast => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 

   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('ALST')
   alst => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('AIST')
   aist => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QLST')
   qlst => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QIST')
   qist => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim) 

   ! --------------------------------------------------------------- !
   ! Loop through grid-boxes and determine:                          !
   ! 1. Effective mean in-cloud convective ice/liquid (deep+shallow) !
   ! 2. Effective mean in-cloud total ice/liquid (ls+convective)     !
   ! --------------------------------------------------------------- !

   do k = 1, pver
   do i = 1, ncol

      if( sh_frac(i,k) <= frac_limit .or. sh_icwmr(i,k) <= ic_limit ) then
          sh0_frac = 0._r8
      else
          sh0_frac = sh_frac(i,k)
      endif
      if( dp_frac(i,k) <= frac_limit .or. dp_icwmr(i,k) <= ic_limit ) then
          dp0_frac = 0._r8
      else
          dp0_frac = dp_frac(i,k)
      endif
      cu0_frac = sh0_frac + dp0_frac

    ! For the moment calculate the emissivity based upon the ls clouds ice fraction

      wrk1 = min(1._r8,max(0._r8, ls_ice(i,k)/(ls_ice(i,k)+ls_liq(i,k)+1.e-36_r8)))

      if( ( cu0_frac < frac_limit ) .or. ( ( sh_icwmr(i,k) + dp_icwmr(i,k) ) < ic_limit ) ) then

            cu0_frac = 0._r8
            cu_icwmr = 0._r8
         
            ls_frac = ast(i,k)
            if( ls_frac < frac_limit ) then
                ls_frac  = 0._r8
                ls_icwmr = 0._r8
            else
                ls_icwmr = ( ls_liq(i,k) + ls_ice(i,k) )/max(frac_limit,ls_frac) ! Convert to IC value.
            end if

            tot0_frac = ls_frac
            tot_icwmr = ls_icwmr
           
      else

          ! Select radiation constants (effective radii) for emissivity averaging.
            
            if( microp_scheme .eq. 'MG' ) then
                kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(i,k)),130._r8)
            elseif( microp_scheme .eq. 'RK' ) then
                kabsi = 0.005_r8 + 1._r8/rei(i,k)
            endif
            kabs  = kabsl * ( 1._r8 - wrk1 ) + kabsi * wrk1
            alpha = -1.66_r8*kabs*pdel(i,k)/gravit*1000.0_r8

          ! Selecting cumulus in-cloud water.            

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted arithmetic average
               cu_icwmr = ( sh0_frac * sh_icwmr(i,k) + dp0_frac*dp_icwmr(i,k))/max(frac_limit,cu0_frac)
            case (2)
               sh0 = exp(alpha*sh_icwmr(i,k))
               dp0 = exp(alpha*dp_icwmr(i,k))               
               cu_icwmr = log((sh0_frac*sh0+dp0_frac*dp0)/max(frac_limit,cu0_frac))
               cu_icwmr = cu_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

          ! Selecting total in-cloud water. 
          ! Attribute large-scale/convective area fraction differently from default.

            ls_frac   = ast(i,k) 
            ls_icwmr  = (ls_liq(i,k) + ls_ice(i,k))/max(frac_limit,ls_frac) ! Convert to IC value.
            tot0_frac = (ls_frac + cu0_frac) 

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted 'arithmetic in emissivity' average
               tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit,tot0_frac)
            case (2)
               tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr)+cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit,tot0_frac))
               tot_icwmr = tot_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

      end if

    ! Repartition convective cloud water into liquid and ice phase.
    ! Currently, this partition is made using the ice fraction of stratus condensate.
    ! In future, we should use ice fraction explicitly computed from the convection scheme.

      conv_ice(i,k) = cu_icwmr * wrk1
      conv_liq(i,k) = cu_icwmr * (1._r8-wrk1)

      tot_ice(i,k)  = tot_icwmr * wrk1
      tot_liq(i,k)  = tot_icwmr * (1._r8-wrk1)

      totg_ice(i,k) = tot0_frac * tot_icwmr * wrk1
      totg_liq(i,k) = tot0_frac * tot_icwmr * (1._r8-wrk1)

   end do
   end do

  ! Output convective IC WMRs
   
   call outfld( 'ICLMRCU ', conv_liq  , pcols, lchnk )
   call outfld( 'ICIMRCU ', conv_ice  , pcols, lchnk )
   call outfld( 'ICWMRSH ', sh_icwmr  , pcols, lchnk )
   call outfld( 'ICWMRDP ', dp_icwmr  , pcols, lchnk ) 
   call outfld( 'ICLMRTOT', tot_liq   , pcols, lchnk )
   call outfld( 'ICIMRTOT', tot_ice   , pcols, lchnk )
   call outfld( 'SH_CLD  ', sh_frac   , pcols, lchnk )
   call outfld( 'DP_CLD  ', dp_frac   , pcols, lchnk )

  end subroutine conv_water_4rad

  end module stratiform
