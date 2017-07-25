module sld_control_mod
  use shr_kind_mod, only : r8=>shr_kind_r8
  use pmgrid, only: plat, plon, plev

! from comqfl.h
!-----------------------------------------------------------------------
!
! Purpose: Global integrals for moisture and mass conservation and geopotential height
!
!-----------------------------------------------------------------------
!
   real(r8) ,public ::  tmass(plat)  ! Mass integral for each latitude pair
   real(r8) ,public ::  tmass0       ! Specified dry mass of atmosphere
   real(r8) ,public ::  tmassf       ! Global mass integral
   real(r8) ,public ::  qmassf       ! Global moisture integral
   real(r8) ,public ::  fixmas       ! Proportionality factor for ps in dry mass fixer
   real(r8) ,public ::  qmass1       ! Contribution to global moisture integral (mass
                         !  weighting is based upon the "A" part of the hybrid grid)
   real(r8) ,public ::  qmass2       ! Contribution to global moisture integral (mass
                         !  weighting is based upon the "B" part of the hybrid grid)
   real(r8) ,public ::  pdela(plon,plev)   ! pressure difference between interfaces (pressure
                         !  defined using the "A" part of hybrid grid only)
   real(r8) ,public ::  zgsint       ! global integral of geopotential height

! from comfft.h

      integer pcray          ! length of vector register (words)
      parameter (pcray=64)
!
      real(r8) trig (3*plon/2+1,plat)         ! trigonometric funct values used by fft
      integer ifax(19,plat)           ! fft factorization of plon/2





end module sld_control_mod
