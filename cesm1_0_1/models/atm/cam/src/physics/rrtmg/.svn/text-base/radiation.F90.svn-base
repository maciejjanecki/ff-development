module radiation

!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to RRTMG
!
! Revision history:
! May  2004, D. B. Coleman,  Initial version of interface module.
! July 2004, B. Eaton,       Use interfaces from new shortwave, longwave, and ozone modules.
! Feb  2005, B. Eaton,       Add namelist variables and control of when calcs are done.
! May  2008, Mike Iacono     Initial version for RRTMG
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8
use spmd_utils,      only: masterproc
use ppgrid,          only: pcols, pver, pverp
use physics_types,   only: physics_state, physics_ptend
use physconst,       only: cpair, cappa
use time_manager,    only: get_nstep
use abortutils,      only: endrun
use error_messages,  only: handle_err
use cam_control_mod, only: lambm0, obliqr, mvelpp, eccen
use scamMod,         only: scm_crm_mode, single_column,have_cld,cldobs,&
                           have_clwp,clwpobs,have_tg,tground
use perf_mod,        only: t_startf, t_stopf
use cam_logfile,     only: iulog

use rad_constituents, only: N_DIAG
use radconstants,     only: rrtmg_sw_cloudsim_band, rrtmg_lw_cloudsim_band, nswbands

implicit none
private
save

public :: &
   radiation_register,    &! registers radiation physics buffer fields
   radiation_defaultopts, &! set default values of namelist variables in runtime_opts
   radiation_setopts,     &! set namelist values from runtime_opts
   radiation_printopts,   &! print namelist values to log
   radiation_get,         &! provide read access to private module data
   radiation_nextsw_cday, &! calendar day of next radiation calculation
   radiation_do,          &! query which radiation calcs are done this timestep
   radiation_init,        &! calls radini
   radiation_tend          ! moved from radctl.F90

! Private module data

! Default values for namelist variables

integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).

integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run

character(len=16) :: microp_scheme  ! microphysics scheme

character(len=4) :: diag(10) =(/'_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

logical :: dohirs = .false. ! diagnostic  brightness temperatures at the top of the
                            ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
                            ! channels (1,2,3,4).
integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

!===============================================================================
contains
!===============================================================================

  subroutine radiation_register
!-----------------------------------------------------------------------
! 
! Register radiation fields in the physics buffer
!
!-----------------------------------------------------------------------

    use phys_buffer,  only: pbuf_times, pbuf_add

    integer idx
    call pbuf_add('QRS' , 'global', 1,pver,1, idx) ! shortwave radiative heating rate 
    call pbuf_add('QRL' , 'global', 1,pver,1, idx) ! longwave  radiative heating rate 

  end subroutine radiation_register

!================================================================================================

subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, irad_always_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = -999
   if ( present(irad_always_out) ) irad_always_out = irad_always

end subroutine radiation_defaultopts

!================================================================================================

subroutine radiation_setopts(dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, &
   irad_always_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! *** NOTE *** This routine needs information about dtime (init by dycore) 
!              and nhtfrq (init by history) to do its checking.  Being called
!              from runtime_opts provides these values possibly before they
!              have been set in the modules responsible for them.
!-----------------------------------------------------------------------

   integer, intent(in)           :: dtime           ! timestep size (s)
   integer, intent(in)           :: nhtfrq          ! output frequency of primary history file
   integer, intent(in), optional :: iradsw_in
   integer, intent(in), optional :: iradlw_in
   integer, intent(in), optional :: iradae_in
   integer, intent(in), optional :: irad_always_in

   ! Local
   integer :: ntspdy   ! no. timesteps per day
   integer :: nhtfrq1  ! local copy of input arg nhtfrq
   integer :: iradae   ! not used by RRTMG
!-----------------------------------------------------------------------

   if ( present(iradsw_in) )      iradsw = iradsw_in
   if ( present(iradlw_in) )      iradlw = iradlw_in
   if ( present(iradae_in) )      iradae = iradae_in
   if ( present(irad_always_in) ) irad_always = irad_always_in

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Has user specified iradae?
   if (iradae /= -999) then
      call endrun('radiation_setopts: iradae not used by RRTMG.')
   end if

end subroutine radiation_setopts

!===============================================================================

subroutine radiation_get(iradsw_out, iradlw_out, iradae_out, irad_always_out)
!----------------------------------------------------------------------- 
! Purpose: Provide access to private module data.  (This should be eliminated.)
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = -999
   if ( present(irad_always_out) ) irad_always_out = irad_always

end subroutine radiation_get

!================================================================================================

subroutine radiation_printopts
!----------------------------------------------------------------------- 
! Purpose: Print runtime options to log.
!-----------------------------------------------------------------------


   if(irad_always /= 0) write(iulog,10) irad_always
   write(iulog,20) iradsw,iradlw
10 format(' Execute SW/LW radiation continuously for the first ',i5, ' timestep(s) of this run')
20 format(' Frequency of Shortwave Radiation calc. (IRADSW)     ',i5/, &
          ' Frequency of Longwave Radiation calc. (IRADLW)      ',i5)

end subroutine radiation_printopts

!================================================================================================

function radiation_do(op, timestep)
!----------------------------------------------------------------------- 
! Purpose: Returns true if the specified operation is done this timestep.
!-----------------------------------------------------------------------

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   !-----------------------------------------------------------------------

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradsw == 1                     &
                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('lw') ! do a longwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradlw == 1                     &
                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
      ! for RRTMG there is no abs/ems restart file
      radiation_do = .false.
         
   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
!----------------------------------------------------------------------- 
! Purpose: Returns calendar day of next sw radiation calculation
!-----------------------------------------------------------------------

   use time_manager, only: get_curr_calday, get_nstep, get_step_size

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

  subroutine radiation_init()
!-----------------------------------------------------------------------
!
! Initialize the radiation parameterization, add fields to the history buffer
! 
!-----------------------------------------------------------------------

    use cam_history,   only: addfld, add_default, phys_decomp
    use physconst,     only: gravit, cpair, epsilo, stebol, &
                             pstd, mwdry, mwco2, mwo3
    use phys_control,  only: phys_getopts
    use radsw,         only: radsw_init
    use radlw,         only: radlw_init
    use radae,         only: radae_init
    use hirsbt,        only: hirsbt_init
    use hirsbtpar,     only: hirsname, msuname
    use cloudsimulator_38, only: doisccp_38, cloudsimulator_init_38=>cloudsimulator_init
    use radiation_data, only: init_rad_data
    use phys_control,   only: phys_getopts

#ifdef MODAL_AERO
    use modal_aer_opt, only: modal_aer_opt_init
    use radsw,         only: wavmin, wavmax
    use radlw,         only: wavenumber1_longwave,wavenumber2_longwave
#endif
    integer :: idiag
    integer :: nstep                       ! current timestep number
    logical :: history_budget              ! output tendencies and state variables for CAM4
                                           ! temperature, water vapor, cloud ice and cloud
                                           ! liquid budgets.
    integer :: history_budget_histfile_num ! output history file number for budget fields
    !-----------------------------------------------------------------------

    call init_rad_data() ! initialize output fields for offline driver

    call phys_getopts(microp_scheme_out=microp_scheme)

    call radsw_init()
    call radlw_init()
#ifdef MODAL_AERO
    call modal_aer_opt_init(wavmin, wavmax, wavenumber1_longwave, wavenumber2_longwave)
#endif

    call hirsbt_init()

    ! "irad_always" is number of time steps to execute radiation continuously from start of
    ! initial OR restart run

    nstep = get_nstep()
    if ( irad_always > 0) then
       nstep       = get_nstep()
       irad_always = irad_always + nstep
    end if

    if (doisccp_38) call cloudsimulator_init_38

    ! Shortwave radiation
    call addfld ('SOLIN   ','W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLL    ','W/m2    ',1,    'A','Solar downward near infrared direct  to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLS    ','W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('SOLLD   ','W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('SOLSD   ','W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRS     ','K/s     ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRSC    ','K/s     ',pver, 'A','Clearsky solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNS    ','W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNT    ','W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOA  ','W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSUTOA  ','W/m2    ',1,    'A','Upwelling solar flux at top of atmosphere',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTOAC ','W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
                                                                                                          sampling_seq='rad_lwsw')
    call addfld ('FSN200  ','W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSN200C ','W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNTC   ','W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNSC   ','W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSDSC   ','W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FSDS    ','W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUS     ','W/m2    ',pverp,'I','Shortwave upward flux',phys_decomp)
    call addfld ('FDS     ','W/m2    ',pverp,'I','Shortwave downward flux',phys_decomp)
    call addfld ('FUSC    ','W/m2    ',pverp,'I','Shortwave clear-sky upward flux',phys_decomp)
    call addfld ('FDSC    ','W/m2    ',pverp,'I','Shortwave clear-sky downward flux',phys_decomp)
    call addfld ('FSNIRTOA','W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAC','W/m2    ',1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
            phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FSNRTOAS','W/m2    ',1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('SWCF    ','W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')

    call add_default ('SOLIN   ', 1, ' ')
    call add_default ('QRS     ', 1, ' ')
    call add_default ('FSNS    ', 1, ' ')
    call add_default ('FSNT    ', 1, ' ')
    call add_default ('FSNTOA  ', 1, ' ')
    call add_default ('FSUTOA  ', 1, ' ')
    call add_default ('FSNTOAC ', 1, ' ')
    call add_default ('FSNTC   ', 1, ' ')
    call add_default ('FSNSC   ', 1, ' ')
    call add_default ('FSDSC   ', 1, ' ')
    call add_default ('FSDS    ', 1, ' ')
    call add_default ('SWCF    ', 1, ' ')
    if (single_column.and.scm_crm_mode) then
       call add_default ('FUS     ', 1, ' ')
       call add_default ('FUSC    ', 1, ' ')
       call add_default ('FDS     ', 1, ' ')
       call add_default ('FDSC    ', 1, ' ')
    endif

    do idiag = 1, N_DIAG
       ! Shortwave radiation
       call addfld ('SOLIN'//diag(idiag),'W/m2    ',1,    'A','Solar insolation',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('SOLL'//diag(idiag),'W/m2    ',1,    'A','Solar downward near infrared direct  to surface', &
	       phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('SOLS'//diag(idiag),'W/m2    ',1,    'A','Solar downward visible direct  to surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('SOLLD'//diag(idiag),'W/m2    ',1,    'A','Solar downward near infrared diffuse to surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('SOLSD'//diag(idiag),'W/m2    ',1,    'A','Solar downward visible diffuse to surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('QRS'//diag(idiag),'K/s     ',pver, 'A','Solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('QRSC'//diag(idiag),'K/s     ',pver, 'A','Clearsky solar heating rate',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSNS'//diag(idiag),'W/m2    ',1,    'A','Net solar flux at surface',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSNT'//diag(idiag),'W/m2    ',1,    'A','Net solar flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSNTOA'//diag(idiag),'W/m2    ',1,    'A','Net solar flux at top of atmosphere',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSNTOAC'//diag(idiag),'W/m2    ',1,    'A','Clearsky net solar flux at top of atmosphere',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSN200'//diag(idiag),'W/m2    ',1,    'A','Net shortwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSN200C'//diag(idiag),'W/m2    ',1,    'A','Clearsky net shortwave flux at 200 mb',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSNTC'//diag(idiag),'W/m2    ',1,    'A','Clearsky net solar flux at top of model',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSNSC'//diag(idiag),'W/m2    ',1,    'A','Clearsky net solar flux at surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSDSC'//diag(idiag),'W/m2    ',1,    'A','Clearsky downwelling solar flux at surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FSDS'//diag(idiag),'W/m2    ',1,    'A','Downwelling solar flux at surface',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FUS'//diag(idiag),'W/m2    ',pverp,'I','Shortwave upward flux',phys_decomp)
       call addfld ('FDS'//diag(idiag),'W/m2    ',pverp,'I','Shortwave downward flux',phys_decomp)
       call addfld ('FUSC'//diag(idiag),'W/m2    ',pverp,'I','Shortwave clear-sky upward flux',phys_decomp)
       call addfld ('FDSC'//diag(idiag),'W/m2    ',pverp,'I','Shortwave clear-sky downward flux',phys_decomp)
       call addfld ('FSNIRTOA'//diag(idiag),'W/m2    ',1,    'A','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
               phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSNRTOAC'//diag(idiag),'W/m2    ', &
               1,    'A','Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
               phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FSNRTOAS'//diag(idiag),'W/m2    ', &
               1,    'A','Net near-infrared flux (>= 0.7 microns) at top of atmosphere',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('SWCF'//diag(idiag),'W/m2    ',1,    'A','Shortwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    enddo

    ! Longwave radiation
    call addfld ('QRL     ','K/s     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('QRLC    ','K/s     ',pver, 'A','Clearsky longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNS    ','W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLDS    ','W/m2    ',1,    'A','Downwelling longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNT    ','W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUT    ','W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLUTC   ','W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
                                                                                                        sampling_seq='rad_lwsw')
    call addfld ('FLNTC   ','W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200  ','W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLN200C ','W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FLNSC   ','W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('LWCF    ','W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
    call addfld ('FUL     ','W/m2    ',pverp,'I','Longwave upward flux',phys_decomp)
    call addfld ('FDL     ','W/m2    ',pverp,'I','Longwave downward flux',phys_decomp)
    call addfld ('FULC    ','W/m2    ',pverp,'I','Longwave clear-sky upward flux',phys_decomp)
    call addfld ('FDLC    ','W/m2    ',pverp,'I','Longwave clear-sky downward flux',phys_decomp)
    call add_default ('QRL     ', 1, ' ')
    call add_default ('FLNS    ', 1, ' ')
    call add_default ('FLDS    ', 1, ' ')
    call add_default ('FLNT    ', 1, ' ')
    call add_default ('FLUT    ', 1, ' ')
    call add_default ('FLUTC   ', 1, ' ')
    call add_default ('FLNTC   ', 1, ' ')
    call add_default ('FLNSC   ', 1, ' ')
    call add_default ('LWCF    ', 1, ' ')
    if (single_column.and.scm_crm_mode) then
       call add_default ('FUL     ', 1, ' ')
       call add_default ('FULC    ', 1, ' ')
       call add_default ('FDL     ', 1, ' ')
       call add_default ('FDLC    ', 1, ' ')
    endif

    do idiag = 1, N_DIAG
       call addfld ('QRL'//diag(idiag),'K/s     ',pver, 'A','Longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('QRLC'//diag(idiag),'K/s     ',pver, 'A','Clearsky longwave heating rate',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FLNS'//diag(idiag),'W/m2    ',1,    'A','Net longwave flux at surface',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FLNT'//diag(idiag),'W/m2    ',1,    'A','Net longwave flux at top of model',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FLUT'//diag(idiag),'W/m2    ',1,    'A','Upwelling longwave flux at top of model',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FLUTC'//diag(idiag),'W/m2    ',1,    'A','Clearsky upwelling longwave flux at top of model',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FLNTC'//diag(idiag),'W/m2    ',1,    'A','Clearsky net longwave flux at top of model',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('LWCF'//diag(idiag),'W/m2    ',1,    'A','Longwave cloud forcing',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FLN200'//diag(idiag),'W/m2    ',1,    'A','Net longwave flux at 200 mb',phys_decomp, sampling_seq='rad_lwsw')
       call addfld ('FLN200C'//diag(idiag),'W/m2    ',1,    'A','Clearsky net longwave flux at 200 mb',phys_decomp, &
               sampling_seq='rad_lwsw')
       call addfld ('FLNSC'//diag(idiag),'W/m2    ',1,    'A','Clearsky net longwave flux at surface',phys_decomp, &
               sampling_seq='rad_lwsw')
   
    enddo

    ! HIRS/MSU diagnostic brightness temperatures
    if (dohirs) then
       call addfld (hirsname(1),'K       ',1,'A','HIRS CH2 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(2),'K       ',1,'A','HIRS CH4 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(3),'K       ',1,'A','HIRS CH6 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(4),'K       ',1,'A','HIRS CH8 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(5),'K       ',1,'A','HIRS CH10 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(6),'K       ',1,'A','HIRS CH11 infra-red brightness temperature',phys_decomp)
       call addfld (hirsname(7),'K       ',1,'A','HIRS CH12 infra-red brightness temperature',phys_decomp)
       call addfld (msuname(1),'K       ',1,'A','MSU CH1 microwave brightness temperature',phys_decomp)
       call addfld (msuname(2),'K       ',1,'A','MSU CH2 microwave brightness temperature',phys_decomp)
       call addfld (msuname(3),'K       ',1,'A','MSU CH3 microwave brightness temperature',phys_decomp)
       call addfld (msuname(4),'K       ',1,'A','MSU CH4 microwave brightness temperature',phys_decomp)
       call add_default (hirsname(1), 1, ' ')
       call add_default (hirsname(2), 1, ' ')
       call add_default (hirsname(3), 1, ' ')
       call add_default (hirsname(4), 1, ' ')
       call add_default (hirsname(5), 1, ' ')
       call add_default (hirsname(6), 1, ' ')
       call add_default (hirsname(7), 1, ' ')
       call add_default (msuname(1), 1, ' ')
       call add_default (msuname(2), 1, ' ')
       call add_default (msuname(3), 1, ' ')
       call add_default (msuname(4), 1, ' ')
    end if

    ! Heating rate needed for d(theta)/dt computation
    call addfld ('HR      ','K/s     ',pver, 'A','Heating rate needed for d(theta)/dt computation',phys_decomp)

    ! (Almost) net radiative flux at surface, does not have lwup.
    call addfld ('SRFRAD  ','W/m2    ',1,    'A','Net radiative flux at surface',phys_decomp)
    call add_default ('SRFRAD  ', 1, ' ')

    call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)

    if ( history_budget .and. history_budget_histfile_num > 1 ) then
       call add_default ('QRL     ', history_budget_histfile_num, ' ')
       call add_default ('QRS     ', history_budget_histfile_num, ' ')
    end if

  end subroutine radiation_init

!===============================================================================
  
  subroutine radiation_tend(state,ptend,pbuf, &
       cam_out, cam_in, &
       landfrac,landm,icefrac,snowh, &
       fsns,    fsnt, flns,    flnt,  &
       fsds, net_flx)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Driver for radiation computation.
    ! 
    ! Method: 
    ! Radiation uses cgs units, so conversions must be done from
    ! model fields to radiation fields.
    !
    ! Revision history:
    ! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
    ! 2004-08-09  B. Eaton         Add pointer variables for constituents.
    ! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
    ! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
    ! 2007-11-05  M. Iacono        Install rrtmg_lw and sw as radiation model.
    ! 2007-12-27  M. Iacono        Modify to use CAM cloud optical properties with rrtmg.
    !-----------------------------------------------------------------------


    use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
    use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
    use physics_types,   only: physics_state, physics_ptend
    use time_manager,    only: get_curr_calday
    use camsrfexch_types,only: cam_out_t, cam_in_t
    use cam_history,     only: outfld
    use cloudsimulator_38, only: doisccp_38, cloudsimulator_run_38=>cloudsimulator_run
    use parrrtm,         only: nbndlw
    use parrrsw,         only: nbndsw
    use hirsbt,          only: hirsrtm
    use hirsbtpar,       only: pnb_hirs, pnf_msu, hirsname, msuname
    use radheat,         only: radheat_tend
    use ppgrid
    use pspect
    use physconst,        only: cpair, stebol
    use radconstants,     only: nlwbands
    use radsw,            only: rad_rrtmg_sw
    use radlw,            only: rad_rrtmg_lw
    use rad_constituents, only: rad_cnst_get_clim_gas, rad_cnst_out, N_DIAG, rad_cnst_get_diag_list, oldcldoptics, &
               liqcldoptics, icecldoptics
    use aer_rad_props,    only: aer_rad_props_get_clim_sw, aer_rad_props_get_clim_lw
    use interpolate_data, only: vertinterp
    use cloud_rad_props,  only: get_ice_optics_sw, get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
               ice_cloud_get_rad_props_lw, cloud_rad_props_get_lw, snow_cloud_get_rad_props_lw, get_snow_optics_sw
    use slingo,           only: slingo_liq_get_rad_props_lw, slingo_liq_optics_sw
    use ebert_curry,      only: ec_ice_optics_sw, ec_ice_get_rad_props_lw
    use rad_solar_var,    only: get_variability
    use radiation_data,   only: output_rad_data

    ! Arguments
    real(r8), intent(in)    :: landfrac(pcols)  ! land fraction
    real(r8), intent(in)    :: landm(pcols)     ! land fraction ramp
    real(r8), intent(in)    :: icefrac(pcols)   ! land fraction
    real(r8), intent(in)    :: snowh(pcols)     ! Snow depth (liquid water equivalent)
    real(r8), intent(inout) :: fsns(pcols)      ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)      ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)      ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)      ! Surface solar down flux
    real(r8), intent(inout) :: net_flx(pcols)

    type(physics_state), intent(in), target :: state
    type(physics_ptend), intent(out)        :: ptend
    type(pbuf_fld),      intent(inout)      :: pbuf(pbuf_size_max)
    type(cam_out_t),     intent(inout)      :: cam_out
    type(cam_in_t),      intent(in)         :: cam_in

    ! Local variables

    logical :: dosw, dolw
    integer nstep                       ! current timestep number
    real(r8) britemp(pcols,pnf_msu)     ! Microwave brightness temperature
    real(r8) tb_ir(pcols,pnb_hirs)      ! Infrared brightness temperature
    real(r8) ts(pcols)                  ! surface temperature
    real(r8) pintmb(pcols,pverp)        ! Model interface pressures (hPa)
    real(r8) oro(pcols)                 ! Land surface flag, sea=0, land=1

    integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
    real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
                                               !    maximally overlapped region.
                                               !    0->pmxrgn(i,1) is range of pressure for
                                               !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
                                               !    2nd region, etc
    real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
    real(r8) :: cldtau(pcols,pver)             ! Cloud longwave optical depth
    real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
    real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
    real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
    real(r8) cllow(pcols)                      !       "     low  cloud cover
    real(r8) clmed(pcols)                      !       "     mid  cloud cover
    real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
    real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables

    ! combined cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: c_cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    real(r8) :: c_cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    real(r8) :: c_cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    real(r8) :: c_cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    real(r8) :: c_cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    real(r8) :: cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    real(r8) :: cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    real(r8) :: cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    real(r8) :: cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: ice_tau    (nbndsw,pcols,pver) ! ice extinction optical depth
    real(r8) :: ice_tau_w  (nbndsw,pcols,pver) ! ice single scattering albedo * tau
    real(r8) :: ice_tau_w_g(nbndsw,pcols,pver) ! ice assymetry parameter * tau * w
    real(r8) :: ice_tau_w_f(nbndsw,pcols,pver) ! ice forward scattered fraction * tau * w
    real(r8) :: ice_lw_abs (nbndlw,pcols,pver)   ! ice absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: snow_tau    (nbndsw,pcols,pver) ! snow extinction optical depth
    real(r8) :: snow_tau_w  (nbndsw,pcols,pver) ! snow single scattering albedo * tau
    real(r8) :: snow_tau_w_g(nbndsw,pcols,pver) ! snow assymetry parameter * tau * w
    real(r8) :: snow_tau_w_f(nbndsw,pcols,pver) ! snow forward scattered fraction * tau * w
    real(r8) :: snow_lw_abs (nbndlw,pcols,pver)   ! snow absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: liq_tau    (nbndsw,pcols,pver) ! liquid extinction optical depth
    real(r8) :: liq_tau_w  (nbndsw,pcols,pver) ! liquid single scattering albedo * tau
    real(r8) :: liq_tau_w_g(nbndsw,pcols,pver) ! liquid assymetry parameter * tau * w
    real(r8) :: liq_tau_w_f(nbndsw,pcols,pver) ! liquid forward scattered fraction * tau * w
    real(r8) :: liq_lw_abs (nbndlw,pcols,pver) ! liquid absorption optics depth (LW)

    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: rel      ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rel_fn   ! liquid effective drop radius at fixed number
                                                  ! for indirect effect (microns)	
    real(r8), pointer, dimension(:,:) :: rei      ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cldfsnow ! cloud fraction of just "snow clouds- whatever they are"
    real(r8) :: cldfprime(pcols,pver)             ! combined cloud fraction (snow plus regular)
    real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction
    real(r8), pointer, dimension(:,:) :: qrs      ! shortwave radiative heating rate 
    real(r8), pointer, dimension(:,:) :: qrl      ! longwave  radiative heating rate 
    real(r8) :: qrsc(pcols,pver)                  ! clearsky shortwave radiative heating rate 
    real(r8) :: qrlc(pcols,pver)                  ! clearsky longwave  radiative heating rate 

    integer lchnk, ncol, lw
    real(r8) :: calday                        ! current calendar day
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
    logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps

    ! Local variables from radctl
    integer i, k                  ! index
    integer :: istat
    real(r8) solin(pcols)         ! Solar incident flux
    real(r8) fsntoa(pcols)        ! Net solar flux at TOA
    real(r8) fsutoa(pcols)        ! Upwelling solar flux at TOA
    real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
    real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
    real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
    real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
    real(r8) flut(pcols)          ! Upward flux at top of model
    real(r8) lwcf(pcols)          ! longwave cloud forcing
    real(r8) swcf(pcols)          ! shortwave cloud forcing
    real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
    real(r8) flntc(pcols)         ! Clear sky lw flux at model top
    real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
    real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) fns(pcols,pverp)     ! net shortwave flux
    real(r8) fcns(pcols,pverp)    ! net clear-sky shortwave flux
    real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
    real(r8) fsn200c(pcols)       ! fcns interpolated to 200 mb
    real(r8) fnl(pcols,pverp)     ! net longwave flux
    real(r8) fcnl(pcols,pverp)    ! net clear-sky longwave flux

    real(r8) pbr(pcols,pver)      ! Model mid-level pressures (dynes/cm2)
    real(r8) pnm(pcols,pverp)     ! Model interface pressures (dynes/cm2)
    real(r8) eccf                 ! Earth/sun distance factor
    real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

    real(r8) dy                   ! Temporary layer pressure thickness
    real(r8) tint(pcols,pverp)    ! Model interface temperature
    real(r8) :: sfac(1:nswbands)  ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band

    real(r8), pointer, dimension(:,:) :: n2o    ! nitrous oxide mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4    ! methane mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o2     ! O2 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc11  ! cfc11 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc12  ! cfc12 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
    real(r8), dimension(pcols) :: co2_col_mean  ! co2 column mean mmr
    real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity

    ! Aerosol radiative properties
    real(r8) :: aer_tau    (pcols,0:pver,nbndsw) ! aerosol extinction optical depth
    real(r8) :: aer_tau_w  (pcols,0:pver,nbndsw) ! aerosol single scattering albedo * tau
    real(r8) :: aer_tau_w_g(pcols,0:pver,nbndsw) ! aerosol assymetry parameter * w * tau
    real(r8) :: aer_tau_w_f(pcols,0:pver,nbndsw) ! aerosol forward scattered fraction * w * tau
    real(r8) :: aer_lw_abs (pcols,pver,nbndlw)   ! aerosol absorption optics depth (LW)

    ! Gathered indicies of day and night columns 
    !  chunk_column_index = IdxDay(daylight_column_index)
    integer :: Nday                      ! Number of daylight columns
    integer :: Nnite                     ! Number of night columns
    integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
    integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

    integer :: idiag  ! index through diagnostic radiation calls
    logical :: diagnosticcall(N_DIAG)

    character(*), parameter :: name = 'radiation_tend'
!----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol = state%ncol

    calday = get_curr_calday()

    itim = pbuf_old_tim_idx()
    ifld = pbuf_get_fld_idx('CLDFSNOW')
    cldfsnow => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    ifld = pbuf_get_fld_idx('CLD')
    cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    ifld = pbuf_get_fld_idx('CONCLD')
    concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    ifld = pbuf_get_fld_idx('QRS')
    qrs => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    ifld = pbuf_get_fld_idx('QRL')
    qrl => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
    ifld = pbuf_get_fld_idx('REL')
    rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    ifld = pbuf_get_fld_idx('REL_FN')
    rel_fn  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    ifld = pbuf_get_fld_idx('REI')
    rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    
!  For CRM, make cloud equal to input observations:
    if (single_column.and.scm_crm_mode.and.have_cld) then
       do k = 1,pver
          cld(:ncol,k)= cldobs(k)
       enddo
    endif

    !
    ! Cosine solar zenith angle for current time step
    !
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)

    call output_rad_data( pbuf, state, cam_in, landm, coszrs )

    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
          Nday = Nday + 1
          IdxDay(Nday) = i
       else
          Nnite = Nnite + 1
          IdxNite(Nnite) = i
       end if
    end do

    dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
    dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?

    if (dosw .or. dolw) then

       ! For CRM, make cloud liquid water path equal to input observations
       if(single_column.and.scm_crm_mode.and.have_clwp)then
          call endrun('cloud water path must be passed through radiation interface')
          !do k=1,pver
          !   cliqwp(:ncol,k) = clwpobs(k)
          !end do
       endif

       call t_startf('cldoptics')

       if (dosw) then
          if(oldcldoptics) then
             call ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.false.)
             call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.false.)
          else
             select case (icecldoptics)
             case ('ebertcurry')
                call  ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.false.)
             case ('mitchell')
                call get_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
             case default
                call endrun('iccldoptics must be one either ebertcurry or mitchell')
             end select
             select case (liqcldoptics)
             case ('slingo')
                call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.false.)
             case ('gammadist')
                call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
             case default
                call endrun('liqcldoptics must be either slingo or gammadist')
             end select
          endif
          cld_tau    (:,1:ncol,:) =  liq_tau    (:,1:ncol,:) + ice_tau    (:,1:ncol,:)
          cld_tau_w  (:,1:ncol,:) =  liq_tau_w  (:,1:ncol,:) + ice_tau_w  (:,1:ncol,:)
          cld_tau_w_g(:,1:ncol,:) =  liq_tau_w_g(:,1:ncol,:) + ice_tau_w_g(:,1:ncol,:)
          cld_tau_w_f(:,1:ncol,:) =  liq_tau_w_f(:,1:ncol,:) + ice_tau_w_f(:,1:ncol,:)
 
          ! add in snow
          call get_snow_optics_sw(state, pbuf, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)
          do i=1,ncol
            do k=1,pver
              cldfprime(i,k)=max(cld(i,k),cldfsnow(i,k))
              if(cldfprime(i,k) > 0.)then
                c_cld_tau    (1:nbndsw,i,k)= &
			(cldfsnow(i,k)*snow_tau    (1:nbndsw,i,k) + cld(i,k)*cld_tau    (1:nbndsw,i,k))/cldfprime(i,k)
                c_cld_tau_w  (1:nbndsw,i,k)= &
			(cldfsnow(i,k)*snow_tau_w  (1:nbndsw,i,k) + cld(i,k)*cld_tau_w  (1:nbndsw,i,k))/cldfprime(i,k)
                c_cld_tau_w_g(1:nbndsw,i,k)= &
			(cldfsnow(i,k)*snow_tau_w_g(1:nbndsw,i,k) + cld(i,k)*cld_tau_w_g(1:nbndsw,i,k))/cldfprime(i,k)
                c_cld_tau_w_f(1:nbndsw,i,k)= &
			(cldfsnow(i,k)*snow_tau_w_f(1:nbndsw,i,k) + cld(i,k)*cld_tau_w_f(1:nbndsw,i,k))/cldfprime(i,k)
              else
                c_cld_tau    (1:nbndsw,i,k)= 0._r8
                c_cld_tau_w  (1:nbndsw,i,k)= 0._r8
                c_cld_tau_w_g(1:nbndsw,i,k)= 0._r8
                c_cld_tau_w_f(1:nbndsw,i,k)= 0._r8
              endif
            enddo
          enddo
       endif

       if (dolw) then
          if(oldcldoptics) then
             call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldcloud=.true.)

          else
             select case (icecldoptics)
             case ('ebertcurry')
                call    ec_ice_get_rad_props_lw(state, pbuf, ice_lw_abs, oldicewp=.false.)
             case ('mitchell')
                call ice_cloud_get_rad_props_lw(state, pbuf, ice_lw_abs)
             case default
                call endrun('iccldoptics must be one either ebertcurry or mitchell')
             end select
             select case (liqcldoptics)
             case ('slingo')
                call   slingo_liq_get_rad_props_lw(state, pbuf, liq_lw_abs, oldliqwp=.false.)
             case ('gammadist')
                call liquid_cloud_get_rad_props_lw(state, pbuf, liq_lw_abs)
             case default
                call endrun('liqcldoptics must be either slingo or gammadist')
             end select
             cld_lw_abs(:,1:ncol,:) = liq_lw_abs(:,1:ncol,:) + ice_lw_abs(:,1:ncol,:)
          endif
          !call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldliq=.true., oldice=.true.)
          !call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldcloud=.true.)
          !call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldliq=.true., oldice=.true.)

          ! add in snow
          call snow_cloud_get_rad_props_lw(state, pbuf, snow_lw_abs)
          do i=1,ncol
            do k=1,pver
              cldfprime(i,k)=max(cld(i,k),cldfsnow(i,k))
              if(cldfprime(i,k) > 0.)then
                c_cld_lw_abs(1:nbndlw,i,k)= &
			(cldfsnow(i,k)*snow_lw_abs(1:nbndlw,i,k) + cld(i,k)*cld_lw_abs(1:nbndlw,i,k))/cldfprime(i,k)
              else
                c_cld_lw_abs(1:nbndlw,i,k)= 0._r8
              endif
            enddo
          enddo
       endif

       call t_stopf('cldoptics')

       ! construct cgs unit reps of pmid and pint and get "eccf" - earthsundistancefactor
       call radinp(ncol, state%pmid, state%pint, pbr, pnm, eccf)

       ! Calculate interface temperatures (following method
       ! used in radtpl for the longwave), using surface upward flux and
       ! stebol constant in mks units
       do i = 1,ncol
          tint(i,1) = state%t(i,1)
          tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i)/stebol))
          do k = 2,pver
             dy = (state%lnpint(i,k) - state%lnpmid(i,k)) / (state%lnpmid(i,k-1) - state%lnpmid(i,k))
             tint(i,k) = state%t(i,k) - dy * (state%t(i,k) - state%t(i,k-1))
          end do
       end do

       ! Solar radiation computation

       if (dosw) then

          call get_variability(sfac)

          ! These are diagnostic shortwave calculations
          call rad_cnst_get_diag_list(diagnosticcall)
          do idiag = 1, N_DIAG
              if (diagnosticcall(idiag)) then
                  ! Get specific humidity
                  call rad_cnst_get_clim_gas('H2O', state, pbuf, sp_hum, diagnosticindex=idiag)
                  ! Get ozone mass mixing ratio.
                  call rad_cnst_get_clim_gas('O3',  state, pbuf, o3, diagnosticindex=idiag)
                  ! Get CO2 mass mixing ratio and compute column mean values
                  call rad_cnst_get_clim_gas('CO2', state, pbuf, co2, diagnosticindex=idiag)
                  ! Get N2O mass mixing ratio
                  call rad_cnst_get_clim_gas('N2O', state, pbuf, n2o, diagnosticindex=idiag)
                  ! Get CH4 mass mixing ratio
                  call rad_cnst_get_clim_gas('CH4', state, pbuf, ch4, diagnosticindex=idiag)
                  ! Get O2  mass mixing ratio
                  call rad_cnst_get_clim_gas('O2', state, pbuf, o2, diagnosticindex=idiag)

                  call aer_rad_props_get_clim_sw(state, pbuf, nnite, idxnite, &
                     aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f, diagnosticindex=idiag)

                  call rad_rrtmg_sw( &
                       lchnk,        ncol,                                                   &
                       state%pint,   state%pmid,   tint,         state%t,                    &
                       sp_hum,       co2,          o3,           ch4,          o2, n2o,          &
                       cldfprime,                                                                  &
                       aer_tau,      aer_tau_w,    aer_tau_w_g,  aer_tau_w_f,                &
                       eccf,         coszrs,       solin,        sfac,         &
                       cam_in%asdir, cam_in%asdif, cam_in%aldir, cam_in%aldif,               &
                       qrs,          qrsc,         fsnt,         fsntc,        fsntoa, fsutoa, &
                       fsntoac,      fsnirt,       fsnrtc,       fsnirtsq,     fsns,           &
                       fsnsc,        fsdsc,        fsds,         cam_out%sols, cam_out%soll,   &
                       cam_out%solsd, cam_out%solld, fns,        fcns,                         &
                       Nday,         Nnite,        IdxDay,       IdxNite,                      &
                       E_cld_tau=c_cld_tau, E_cld_tau_w=c_cld_tau_w, E_cld_tau_w_g=c_cld_tau_w_g, E_cld_tau_w_f=c_cld_tau_w_f, &
                       old_convert = .false.)

                  !  Output net fluxes at 200 mb
                  call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, fsn200c)
                  call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns, fsn200)

                  do i=1,ncol
                     swcf(i)=fsntoa(i) - fsntoac(i)
                  end do
        
                  ! Dump shortwave radiation information to history tape buffer (diagnostics)
                  ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
                  call outfld('QRS'//diag(idiag),ftem  ,pcols,lchnk)
                  ftem(:ncol,:pver) = qrsc(:ncol,:pver)/cpair
                  call outfld('QRSC'//diag(idiag),ftem  ,pcols,lchnk)
                  call outfld('SOLIN'//diag(idiag),solin ,pcols,lchnk)
                  call outfld('FSDS'//diag(idiag),fsds  ,pcols,lchnk)
                  call outfld('FSNIRTOA'//diag(idiag),fsnirt,pcols,lchnk)
                  call outfld('FSNRTOAC'//diag(idiag),fsnrtc,pcols,lchnk)
                  call outfld('FSNRTOAS'//diag(idiag),fsnirtsq,pcols,lchnk)
                  call outfld('FSNT'//diag(idiag),fsnt  ,pcols,lchnk)
                  call outfld('FSNS'//diag(idiag),fsns  ,pcols,lchnk)
                  call outfld('FSNTC'//diag(idiag),fsntc ,pcols,lchnk)
                  call outfld('FSNSC'//diag(idiag),fsnsc ,pcols,lchnk)
                  call outfld('FSDSC'//diag(idiag),fsdsc ,pcols,lchnk)
                  call outfld('FSNTOA'//diag(idiag),fsntoa,pcols,lchnk)
                  call outfld('FSNTOAC'//diag(idiag),fsntoac,pcols,lchnk)
                  call outfld('SOLS'//diag(idiag),cam_out%sols  ,pcols,lchnk)
                  call outfld('SOLL'//diag(idiag),cam_out%soll  ,pcols,lchnk)
                  call outfld('SOLSD'//diag(idiag),cam_out%solsd ,pcols,lchnk)
                  call outfld('SOLLD'//diag(idiag),cam_out%solld ,pcols,lchnk)
                  call outfld('FSN200'//diag(idiag),fsn200,pcols,lchnk)
                  call outfld('FSN200C'//diag(idiag),fsn200c,pcols,lchnk)
                  call outfld('SWCF'//diag(idiag),swcf  ,pcols,lchnk)

              endif !(diagnosticcall(idiag))
          enddo !idiag = 1, N_DIAG

          ! This shortwave calculation forces the climate.
          call t_startf('rad_sw')

          ! Get specific humidity
          call rad_cnst_get_clim_gas('H2O', state, pbuf, sp_hum)
          ! Get ozone mass mixing ratio.
          call rad_cnst_get_clim_gas('O3',  state, pbuf, o3)
          ! Get CO2 mass mixing ratio and compute column mean values
          call rad_cnst_get_clim_gas('CO2', state, pbuf, co2)
          ! Get N2O mass mixing ratio
          call rad_cnst_get_clim_gas('N2O', state, pbuf, n2o)
          ! Get CH4 mass mixing ratio
          call rad_cnst_get_clim_gas('CH4', state, pbuf, ch4)
          ! Get O2  mass mixing ratio
          call rad_cnst_get_clim_gas('O2', state, pbuf, o2)

          call t_startf('aero_optics_sw')
          call aer_rad_props_get_clim_sw(state, pbuf, nnite, idxnite, &
               aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)
          call t_stopf('aero_optics_sw')

          call rad_rrtmg_sw( &
               lchnk,        ncol,                                                   &
               state%pint,   state%pmid,   tint,         state%t,                    &
               sp_hum,       co2,          o3,           ch4,          o2, n2o,          &
               cldfprime,                                                                  &
               aer_tau,      aer_tau_w,    aer_tau_w_g,  aer_tau_w_f,                &
               eccf,         coszrs,       solin,        sfac,         &
               cam_in%asdir, cam_in%asdif, cam_in%aldir, cam_in%aldif,               &
               qrs,          qrsc,         fsnt,         fsntc,        fsntoa, fsutoa, &
               fsntoac,      fsnirt,       fsnrtc,       fsnirtsq,     fsns,           &
               fsnsc,        fsdsc,        fsds,         cam_out%sols, cam_out%soll,   &
               cam_out%solsd, cam_out%solld, fns,        fcns,                         &
               Nday,         Nnite,        IdxDay,       IdxNite,                      &
               E_cld_tau=c_cld_tau, E_cld_tau_w=c_cld_tau_w, E_cld_tau_w_g=c_cld_tau_w_g, E_cld_tau_w_f=c_cld_tau_w_f,   &
               old_convert =.false.)!, ancientmethod=.true.)

          call t_stopf('rad_sw')

          !  Output net fluxes at 200 mb
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, fsn200c)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns, fsn200)

          do i=1,ncol
             swcf(i)=fsntoa(i) - fsntoac(i)
          end do

          ! Dump shortwave radiation information to history tape buffer (diagnostics)
          ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
          call outfld('QRS     ',ftem  ,pcols,lchnk)
          ftem(:ncol,:pver) = qrsc(:ncol,:pver)/cpair
          call outfld('QRSC    ',ftem  ,pcols,lchnk)
          call outfld('SOLIN   ',solin ,pcols,lchnk)
          call outfld('FSDS    ',fsds  ,pcols,lchnk)
          call outfld('FSNIRTOA',fsnirt,pcols,lchnk)
          call outfld('FSNRTOAC',fsnrtc,pcols,lchnk)
          call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
          call outfld('FSNT    ',fsnt  ,pcols,lchnk)
          call outfld('FSNS    ',fsns  ,pcols,lchnk)
          call outfld('FSNTC   ',fsntc ,pcols,lchnk)
          call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
          call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
          call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
          call outfld('FSUTOA  ',fsutoa,pcols,lchnk)
          call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
          call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
          call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
          call outfld('SOLSD   ',cam_out%solsd ,pcols,lchnk)
          call outfld('SOLLD   ',cam_out%solld ,pcols,lchnk)
          call outfld('FSN200  ',fsn200,pcols,lchnk)
          call outfld('FSN200C ',fsn200c,pcols,lchnk)
          call outfld('SWCF    ',swcf  ,pcols,lchnk)

       end if   ! dosw

       ! Output aerosol mmr
       call rad_cnst_out(state, pbuf)

       ! Longwave radiation computation

       if (dolw) then
          !
          ! Convert upward longwave flux units to CGS
          !
          do i=1,ncol
             lwupcgs(i) = cam_in%lwup(i)*1000._r8
             if(single_column.and.scm_crm_mode.and.have_tg) &
                  lwupcgs(i) = 1000*stebol*tground(1)**4
          end do

          call rad_cnst_get_diag_list(diagnosticcall)
          do idiag = 1, N_DIAG
              if (diagnosticcall(idiag)) then

                  ! Get specific humidity
                  call rad_cnst_get_clim_gas('H2O', state, pbuf, sp_hum, diagnosticindex=idiag)
                  ! Get oxygen mass mixing ratio.
                  call rad_cnst_get_clim_gas('O2',  state, pbuf, o2, diagnosticindex=idiag)
                  ! Get ozone mass mixing ratio.
                  call rad_cnst_get_clim_gas('O3',  state, pbuf, o3, diagnosticindex=idiag)
                  ! Get CO2 mass mixing ratio
                  call rad_cnst_get_clim_gas('CO2', state, pbuf, co2, diagnosticindex=idiag)
                  ! Get N2O mass mixing ratio
                  call rad_cnst_get_clim_gas('N2O', state, pbuf, n2o, diagnosticindex=idiag)
                  ! Get CH4 mass mixing ratio
                  call rad_cnst_get_clim_gas('CH4', state, pbuf, ch4, diagnosticindex=idiag)
                  ! Get CFC mass mixing ratios
                  call rad_cnst_get_clim_gas('CFC11', state, pbuf, cfc11, diagnosticindex=idiag)
                  call rad_cnst_get_clim_gas('CFC12', state, pbuf, cfc12, diagnosticindex=idiag)

                  call aer_rad_props_get_clim_lw(state, pbuf, aer_lw_abs, diagnosticindex=idiag)

                  call rad_rrtmg_lw( &
                         lchnk,        ncol,                                                    &
                         tint,         state%t,      sp_hum,       co2,    o3,                  &
                         state%pmid,   state%pint,   aer_lw_abs,                                &
                         n2o,          ch4,          o2,           cfc11,        cfc12,         &
                         cldfprime,    c_cld_lw_abs,                                            &
                         qrl,          qrlc,                                                    &
                         flns,         flnt,         flnsc,        flntc,        cam_out%flwds, &
                         flut,         flutc,        fnl,          fcnl)

                  do i=1,ncol
                     lwcf(i)=flutc(i) - flut(i)
                  end do

                  !  Output fluxes at 200 mb
                  call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl, fln200)
                  call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, fln200c)
                  ! Dump longwave radiation information to history tape buffer (diagnostics)
                  call outfld('QRL'//diag(idiag),qrl (:ncol,:)/cpair,ncol,lchnk)
                  call outfld('QRLC'//diag(idiag),qrlc(:ncol,:)/cpair,ncol,lchnk)
                  call outfld('FLNT'//diag(idiag),flnt  ,pcols,lchnk)
                  call outfld('FLUT'//diag(idiag),flut  ,pcols,lchnk)
                  call outfld('FLUTC'//diag(idiag),flutc ,pcols,lchnk)
                  call outfld('FLNTC'//diag(idiag),flntc ,pcols,lchnk)
                  call outfld('FLNS'//diag(idiag),flns  ,pcols,lchnk)
                  call outfld('FLNSC'//diag(idiag),flnsc ,pcols,lchnk)
                  call outfld('LWCF'//diag(idiag),lwcf  ,pcols,lchnk)
                  call outfld('FLN200'//diag(idiag),fln200,pcols,lchnk)
                  call outfld('FLN200C'//diag(idiag),fln200c,pcols,lchnk)

              endif ! diagnosticcall
          enddo ! idiag = 1, N_DIAG

          call t_startf('rad_lw')

          ! Get specific humidity
          call rad_cnst_get_clim_gas('H2O', state, pbuf, sp_hum)
          ! Get ozone mass mixing ratio.
          call rad_cnst_get_clim_gas('O2',  state, pbuf, o2)
          ! Get ozone mass mixing ratio.
          call rad_cnst_get_clim_gas('O3',  state, pbuf, o3)
          ! Get CO2 mass mixing ratio
          call rad_cnst_get_clim_gas('CO2', state, pbuf, co2)
          ! Get N2O mass mixing ratio
          call rad_cnst_get_clim_gas('N2O', state, pbuf, n2o)
          ! Get CH4 mass mixing ratio
          call rad_cnst_get_clim_gas('CH4', state, pbuf, ch4)
          ! Get CFC mass mixing ratios
          call rad_cnst_get_clim_gas('CFC11', state, pbuf, cfc11)
          call rad_cnst_get_clim_gas('CFC12', state, pbuf, cfc12)

          call t_startf('aero_optics_lw')
          call aer_rad_props_get_clim_lw(state, pbuf, aer_lw_abs)
          call t_stopf('aero_optics_lw')

          call rad_rrtmg_lw( &
               lchnk,        ncol,                                                    &
               tint,         state%t,      sp_hum,       co2, o3,                     &
               state%pmid,   state%pint,   aer_lw_abs,                                &
               n2o,          ch4,          o2,           cfc11,        cfc12,         &
               cldfprime,    c_cld_lw_abs,                                            &
               qrl,          qrlc,                                                    &
               flns,         flnt,         flnsc,        flntc,        cam_out%flwds, &
               flut,         flutc,        fnl,          fcnl)

          call t_stopf('rad_lw')

          !  Output fluxes at 200 mb
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl, fln200)
          call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, fln200c)

          do i=1,ncol
             lwcf(i)=flutc(i) - flut(i)
          end do

          ! Dump longwave radiation information to history tape buffer (diagnostics)
          call outfld('QRL     ',qrl (:ncol,:)/cpair,ncol,lchnk)
          call outfld('QRLC    ',qrlc(:ncol,:)/cpair,ncol,lchnk)
          call outfld('FLNT    ',flnt  ,pcols,lchnk)
          call outfld('FLUT    ',flut  ,pcols,lchnk)
          call outfld('FLUTC   ',flutc ,pcols,lchnk)
          call outfld('FLNTC   ',flntc ,pcols,lchnk)
          call outfld('FLNS    ',flns  ,pcols,lchnk)
          call outfld('FLDS    ',cam_out%flwds ,pcols,lchnk)
          call outfld('FLNSC   ',flnsc ,pcols,lchnk)
          call outfld('LWCF    ',lwcf  ,pcols,lchnk)
          call outfld('FLN200  ',fln200,pcols,lchnk)
          call outfld('FLN200C ',fln200c,pcols,lchnk)

       end if  !dolw

       ! mji/hirsrtm - Add call to HIRSRTM package
       ! HIRS brightness temperature calculation in 7 infra-red channels and 4 microwave
       ! channels as a diagnostic to compare to TOV/MSU satellite data.
       ! Done if dohirs set to .true. at time step frequency ihirsfq

       nstep = get_nstep()

       if ( dohirs .and. (mod(nstep-1,ihirsfq) .eq. 0) ) then

          do i= 1, ncol
             ts(i) = sqrt(sqrt(cam_in%lwup(i)/stebol))
             ! Set oro (land/sea flag) for compatibility with landfrac/icefrac/ocnfrac
             ! oro=0 (sea or ice); oro=1 (land)
             if (landfrac(i).ge.0.001) then
                oro(i)=1.
             else
                oro(i)=0.
             endif
             ! Convert pressure from Pa to hPa
             do k = 1, pver
                pintmb(i,k) = state%pint(i,k)*1.e-2_r8        
             end do
             pintmb(i,pverp) = state%pint(i,pverp)*1.e-2_r8 
          end do

          call calc_col_mean(state, co2, co2_col_mean)
          call hirsrtm(lchnk  ,ncol , &
             pintmb ,state%t  ,sp_hum ,co2_col_mean, &
             o3     ,ts       ,oro    ,tb_ir  ,britemp )

          do i = 1, pnb_hirs
             call outfld(hirsname(i),tb_ir(1,i),pcols,lchnk)
          end do
          do i = 1, pnf_msu
             call outfld(msuname(i),britemp(1,i),pcols,lchnk)
          end do

       end if

       ! Run ISCCP cloud simulator
       if (doisccp_38) then
          emis(:ncol,:) = 1._r8 - exp(-cld_lw_abs(rrtmg_lw_cloudsim_band,:ncol,:))
          call cloudsimulator_run_38(state, cam_in%ts, concld, cld, &
                  cld_tau(rrtmg_sw_cloudsim_band,:,:), emis, coszrs  )
       end if

    else   !  if (dosw .or. dolw) then

       ! convert radiative heating rates from Q*dp to Q for energy conservation
       if (conserve_energy) then
!DIR$ CONCURRENT
          do k =1 , pver
!DIR$ CONCURRENT
             do i = 1, ncol
                qrs(i,k) = qrs(i,k)/state%pdel(i,k)
                qrl(i,k) = qrl(i,k)/state%pdel(i,k)
             end do
          end do
       end if

    end if   !  if (dosw .or. dolw) then

    ! Compute net radiative heating tendency
    call radheat_tend(state, pbuf, ptend, qrl, qrs, fsns, &
                      fsnt, flns, flnt, cam_in%asdir, net_flx)

    ! Compute heating rate for dtheta/dt 
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
       end do
    end do
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
!DIR$ CONCURRENT
       do k =1 , pver
!DIR$ CONCURRENT
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)*state%pdel(i,k)
             qrl(i,k) = qrl(i,k)*state%pdel(i,k)
          end do
       end do
    end if
 
    ! Compute net surface radiative flux for use by surface temperature code.
    ! Note that units have already been converted to mks in RADCTL.  Since
    ! fsns and flwds are in the buffer, array values will be carried across
    ! timesteps when the radiation code is not invoked.
    cam_out%srfrad(:ncol) = fsns(:ncol) + cam_out%flwds(:ncol)
    call outfld('SRFRAD  ',cam_out%srfrad,pcols,lchnk)

 end subroutine radiation_tend

!===============================================================================

subroutine radinp(ncol, pmid, pint, pmidrd, pintrd, eccf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
! Convert model pressures to cgs.
! 
! Author: CCM1, CMS Contact J. Kiehl
!-----------------------------------------------------------------------
   use shr_orb_mod
   use time_manager, only: get_curr_calday

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)
!
! Output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor

!
!---------------------------Local variables-----------------------------
!
   integer i                ! Longitude loop index
   integer k                ! Vertical loop index

   real(r8) :: calday       ! current calendar day
   real(r8) :: delta        ! Solar declination angle
!-----------------------------------------------------------------------
!
   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

!
! Convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0_r8
         pintrd(i,k) = pint(i,k)*10.0_r8
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0_r8
   end do

end subroutine radinp

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)
!----------------------------------------------------------------------- 
! 
! Compute the column mean mass mixing ratio.  
!
!-----------------------------------------------------------------------

   use cam_logfile,  only: iulog

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

!===============================================================================

end module radiation

