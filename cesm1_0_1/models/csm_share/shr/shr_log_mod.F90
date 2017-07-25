!===============================================================================
! SVN $Id: shr_const_mod.F90 6354 2007-09-11 22:49:33Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk/shr/shr_const_mod.F90 $
!===============================================================================

MODULE shr_log_mod

   use shr_kind_mod

   !----------------------------------------------------------------------------
   ! low-level shared variables for logging, these may not be parameters
   !----------------------------------------------------------------------------
   public

   integer(SHR_KIND_IN) :: shr_log_Level = 1
   integer(SHR_KIND_IN) :: shr_log_Unit  = 6

END MODULE shr_log_mod
