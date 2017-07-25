






      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)


      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)



!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(1) = ( + rxt(1) + rxt(3) + het_rates(1))* y(1)
         prod(1) = 0._r8
         loss(2) = 0._r8
         prod(2) =rxt(4)*y(3)
         loss(3) = ( + rxt(4) + het_rates(2))* y(3)
         prod(3) = (rxt(5) +.500_r8*rxt(6) +rxt(7))*y(4)
         loss(4) = ( + rxt(5) + rxt(6) + rxt(7))* y(4)
         prod(4) = 0._r8
         loss(5) = ( + rxt(8))* y(5)
         prod(5) = 0._r8
         loss(6) = 0._r8
         prod(6) = 0._r8
         loss(7) = 0._r8
         prod(7) = 0._r8
         loss(8) = 0._r8
         prod(8) = 0._r8
         loss(9) = 0._r8
         prod(9) = 0._r8
         loss(10) = 0._r8
         prod(10) = 0._r8
         loss(11) = 0._r8
         prod(11) = 0._r8
         loss(12) = 0._r8
         prod(12) = 0._r8
         loss(13) = 0._r8
         prod(13) = 0._r8
         loss(14) = 0._r8
         prod(14) = 0._r8
         loss(15) = 0._r8
         prod(15) = 0._r8
         loss(16) = 0._r8
         prod(16) = 0._r8
         loss(17) = 0._r8
         prod(17) = 0._r8
         loss(18) = 0._r8
         prod(18) = 0._r8
         loss(19) = 0._r8
         prod(19) = 0._r8
         loss(20) = 0._r8
         prod(20) = 0._r8
         loss(21) = 0._r8
         prod(21) = 0._r8
         loss(22) = 0._r8
         prod(22) = 0._r8
         loss(23) = 0._r8
         prod(23) = 0._r8
         loss(24) = 0._r8
         prod(24) = 0._r8
         loss(25) = 0._r8
         prod(25) = 0._r8
         loss(26) = 0._r8
         prod(26) = 0._r8
         loss(27) = 0._r8
         prod(27) = 0._r8
         loss(28) = 0._r8
         prod(28) = 0._r8
         loss(29) = 0._r8
         prod(29) = 0._r8
         loss(30) = 0._r8
         prod(30) = 0._r8
         loss(31) = 0._r8
         prod(31) = 0._r8
         loss(32) = 0._r8
         prod(32) = 0._r8
         loss(33) = 0._r8
         prod(33) = 0._r8
         loss(34) = 0._r8
         prod(34) = 0._r8
         loss(35) = 0._r8
         prod(35) = 0._r8
         loss(36) = 0._r8
         prod(36) = 0._r8
         loss(37) = 0._r8
         prod(37) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss
