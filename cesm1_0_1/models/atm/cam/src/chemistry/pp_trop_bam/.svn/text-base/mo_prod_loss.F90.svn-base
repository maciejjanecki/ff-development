      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid,       only : pver

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)    ::  y(:,:,:)
      real(r8), intent(in)    ::  rxt(:,:,:)
      real(r8), intent(in)    ::  het_rates(:,:,:)


      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid,       only : pver

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)    ::  y(:)
      real(r8), intent(in)    ::  rxt(:)
      real(r8), intent(in)    ::  het_rates(:)



!--------------------------------------------------------------------
!       ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(1) = ( + het_rates(1))* y(1)
         prod(1) = 0._r8
         loss(2) = ( + het_rates(2))* y(2)
         prod(2) = 0._r8
         loss(3) = 0._r8
         prod(3) = 0._r8
         loss(4) = 0._r8
         prod(4) = 0._r8
         loss(5) = ( + rxt(1))* y(7)
         prod(5) = 0._r8
         loss(6) = 0._r8
         prod(6) =rxt(1)*y(7)
         loss(7) = ( + rxt(2))* y(5)
         prod(7) = 0._r8
         loss(8) = 0._r8
         prod(8) =rxt(2)*y(5)
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

      end subroutine imp_prod_loss

      end module mo_prod_loss
