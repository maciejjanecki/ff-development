module inital
!-----------------------------------------------------------------------
!BOP
! !MODULE:  inital --- Define initial conditions for first run of case
!
!----------------------------------------------------------------------- 
! !USES:

   implicit none

   private   ! By default everything private to this module

! !PUBLIC MEMBER FUNCTIONS:

   public cam_initial   ! Cam initialization (formally inital)
!
! !DESCRIPTION: Module for CAM initialization
!
! !REVISION HISTORY:
!   05.08.11   Kluzek     Creation
!   05.11.10   Sawyer     Now using dyn_import/export_t containers
!   06.04.13   Sawyer     Removed dependency on prognostics
!
!
!EOP
!-----------------------------------------------------------------------


contains

!
!----------------------------------------------------------------------- 
!

!-----------------------------------------------------------------------
!BOP
!ROUTINE:  inital --- Define initial conditions for first run of case
!
! !INTERFACE:

subroutine cam_initial( dyn_in, dyn_out, NLFileName )

! !USES:
   use shr_kind_mod,         only : r8 => shr_kind_r8
   use physconst,             only : omega, rearth,         &
                                    rair, cpair,         &
                                    zvir,  pi
   use pmgrid,               only : plon, plat, plev, plevp,                   &
                                    beglonxy, endlonxy, beglatxy, endlatxy,    &
                                    beglat,   endlat,   beglev,   endlev,      &
                                    npr_y, npr_z, nprxy_x, nprxy_y,            &
                                    myid_y, myid_z, myidxy_x, myidxy_y,        &
                                    twod_decomp, mod_geopk, mod_transpose,     &
                                    mod_gatscat

   use dyn_comp,             only : dyn_import_t, dyn_export_t, dyn_init
   use constituents,         only : pcnst
   use phys_grid,            only : phys_grid_init
   use chem_surfvals,        only : chem_surfvals_init
   use time_manager,         only : get_step_size
   use startup_initialconds, only : setup_initial, initial_conds
   use hycoef,               only : hyai, hybi
   use dynamics_vars,        only : T_FVDYCORE_STATE   
   use dyn_internal_state,   only : get_dyn_state
   use fv_control_mod,       only : kord, jord, iord, nsplit, nspltrac, nspltvrm, dyn_conservative, filtcw
!-----------------------------------------------------------------------
!
! Arguments
!
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out
   character(len=*), intent(in) :: NLFileName
!------------------------------Parameters-------------------------------
! !DESCRIPTION:
!
!   Define initial conditions for first run of case
!
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.03.01      Acker         Modifications
!   96.04.01      Boville       Reviewed
!   01.06.17      Sawyer        Added call to dynamics_init
!   01.07.12      Sawyer        Added arguments to dynamics_init
!   05.11.10      Sawyer        Added dyn_in, dyn_out, dyn_state to init
!   06.04.13      Sawyer        Removed call to initial_prognostics
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   real(r8), parameter ::  D0_0                  =   0.0_r8
   real(r8), parameter ::  D1E5                  =   1.0e5_r8
   type (T_FVDYCORE_STATE), pointer :: dyn_state
   real(r8) :: dtime          ! timestep size
   real(r8), allocatable :: ak(:), bk(:)
   integer :: k               ! Index
   integer :: ks              ! transition index to pressure coordinates
!
!-----------------------------------------------------------------------
   call setup_initial()
   !
   !
   dtime = get_step_size()

   allocate( ak(plev+1) )
   allocate( bk(plev+1) )
   do k = 1, plev+1
      ak(k) = hyai(k) * D1E5
      bk(k) = hybi(k)
      if( bk(k) == D0_0 ) ks = k-1
   end do

   !
   ! Initialize dynamics
   !
   dyn_state => get_dyn_state()

   call dyn_init(  "NOT USED", "NOT_USED", plon, plat, plev, dtime, NSPLIT,&
                   NSPLTRAC, NSPLTVRM, IORD, JORD, KORD, 0, AK, BK, pcnst, pcnst, KS,& ! "0" temporary
                   filtcw, beglonxy, endlonxy, beglatxy, endlatxy,         &
                   beglat, endlat, beglev, endlev, pi,           &
                   omega, cpair, rair,     &
                   rearth, rair/cpair,     &
                   zvir,  (/ npr_y, npr_z, nprxy_x, nprxy_y /),  &
                   twod_decomp, mod_geopk, mod_transpose, mod_gatscat,     &
                   dyn_conservative, dyn_state, dyn_in, dyn_out, NLFileName )
   deallocate( ak, bk )
   !
   ! Initialize dynamics grid
   !
   call initcom
   !
   ! Define physics data structures
   !
   call phys_grid_init
   !
   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   !
   call chem_surfvals_init()
   call initial_conds( dyn_in )
!EOC
end subroutine cam_initial

!
!----------------------------------------------------------------------- 
!

end module inital
