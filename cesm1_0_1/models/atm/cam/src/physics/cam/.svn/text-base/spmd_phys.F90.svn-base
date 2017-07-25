
module spmd_phys

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM for physics.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid, only: begchunk, endchunk
   implicit none

   private
   public spmdinit_phys

CONTAINS

!========================================================================

   subroutine spmdinit_phys ()
      use dyn_grid, only: get_dyn_grid_parm
      begchunk = get_dyn_grid_parm('beglat')
      endchunk = get_dyn_grid_parm('endlat')

      return
   end subroutine spmdinit_phys

end module spmd_phys
