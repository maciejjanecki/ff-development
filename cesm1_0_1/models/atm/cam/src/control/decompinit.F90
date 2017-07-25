
subroutine decompinit
!----------------------------------------------------------------------- 
! 
! Purpose: Set up dynamics and physics decompositions
! 
! Method: Ask dynamics and physics packages to do it.  
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

#if ( defined SPMD )
   use spmd_dyn,  only: spmdinit_dyn
   use spmd_phys, only: spmdinit_phys
#endif

   implicit none

#if ( defined SPMD )
!
! Currently spmdinit_dyn must be called before spmdinit_phys because the latter just copies
! in data computed in the former
!
   call spmdinit_dyn ()
   call spmdinit_phys ()
#endif   
   return
end subroutine decompinit

