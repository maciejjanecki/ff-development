
      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine     nlnmat( mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)





















      call     nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      subroutine     nlnmat_finit( mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(   1) = lmat(   1)
         mat(   3) = lmat(   3)
         mat(   4) = lmat(   4)
         mat(   5) = lmat(   5)
         mat(   6) = lmat(   6)
         mat(   2) = 0._r8
         mat(   7) = 0._r8
         mat(   8) = 0._r8
         mat(   9) = 0._r8
         mat(  10) = 0._r8
         mat(  11) = 0._r8
         mat(  12) = 0._r8
         mat(  13) = 0._r8
         mat(  14) = 0._r8
         mat(  15) = 0._r8
         mat(  16) = 0._r8
         mat(  17) = 0._r8
         mat(  18) = 0._r8
         mat(  19) = 0._r8
         mat(  20) = 0._r8
         mat(  21) = 0._r8
         mat(  22) = 0._r8
         mat(   1) = mat(   1) - dti
         mat(   2) = -dti
         mat(   4) = mat(   4) - dti
         mat(   6) = mat(   6) - dti
         mat(   7) = -dti
         mat(   8) = -dti
         mat(   9) = -dti
         mat(  10) = -dti
         mat(  11) = -dti
         mat(  12) = -dti
         mat(  13) = -dti
         mat(  14) = -dti
         mat(  15) = -dti
         mat(  16) = -dti
         mat(  17) = -dti
         mat(  18) = -dti
         mat(  19) = -dti
         mat(  20) = -dti
         mat(  21) = -dti
         mat(  22) = -dti

      end subroutine nlnmat_finit

      end module mo_nln_matrix
