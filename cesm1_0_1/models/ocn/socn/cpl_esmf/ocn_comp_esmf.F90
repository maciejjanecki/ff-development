module ocn_comp_esmf

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

  public :: ocn_init_esmf
  public :: ocn_run_esmf
  public :: ocn_final_esmf
  public :: ocn_register_esmf

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ocn_register_esmf(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In ocn register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETINIT, &
      ocn_init_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETRUN, &
      ocn_run_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_SETFINAL, &
      ocn_final_esmf, ESMF_SINGLEPHASE, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)


end subroutine

!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_init_esmf
!
! !DESCRIPTION:
!     initialize dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ocn_init_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

    rc = ESMF_SUCCESS

    ! Set flag to specify dead components
    call ESMF_AttributeSet(export_state, name="ocn_present", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(export_state, name="ocn_prognostic", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

    call ESMF_AttributeSet(export_state, name="ocnrof_prognostic", value=.false., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine ocn_init_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_run_esmf
!
! !DESCRIPTION:
!     run method for dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ocn_run_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc


!EOP

   rc = ESMF_SUCCESS

end subroutine ocn_run_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_final_esmf
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ocn_final_esmf(comp, import_state, export_state, EClock, rc)

! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
 

   rc = ESMF_SUCCESS

end subroutine ocn_final_esmf
!===============================================================================

end module ocn_comp_esmf
