module atm_comp_esmf

! !USES:

  use ESMF_MOD
  use esmfshr_mod
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_esmf
  public :: atm_run_esmf
  public :: atm_final_esmf
  public :: atm_register_esmf

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine atm_register_esmf(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In atm register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETINIT, &
      atm_init_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETRUN, &
      atm_run_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETFINAL, &
      atm_final_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)


end subroutine

!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_init_esmf
!
! !DESCRIPTION:
!     initialize dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_init_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

    rc = ESMF_SUCCESS

    ! Set flag to specify dead components
!    call ESMF_AttributeSet(export_state, name="atm_present", value=.false., rc=rc)
!    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(export_state, name="atm_prognostic", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine atm_init_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_run_esmf
!
! !DESCRIPTION:
!     run method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_run_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

   rc = ESMF_SUCCESS

end subroutine atm_run_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_final_esmf
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_final_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
 

   rc = ESMF_SUCCESS

end subroutine atm_final_esmf
!===============================================================================

end module atm_comp_esmf
