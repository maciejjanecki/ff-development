
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,45) = 2.2e-10_r8
      rate(:,:,46) = 1.1e-10_r8
      rate(:,:,57) = 6.7e-11_r8
      rate(:,:,58) = 4.9e-11_r8
      rate(:,:,76) = 1.5e-10_r8
      rate(:,:,83) = 9.e-12_r8
      rate(:,:,88) = 1.e-14_r8
      rate(:,:,93) = 2.e-13_r8
      rate(:,:,94) = 6.8e-14_r8
      rate(:,:,108) = 1e-12_r8
      rate(:,:,121) = 5.4e-11_r8
      rate(:,:,123) = 3.5e-12_r8
      rate(:,:,127) = 6.8e-13_r8
      rate(:,:,133) = 3.e-12_r8
      rate(:,:,134) = 1.e-11_r8
      rate(:,:,138) = 1.1e-11_r8
      rate(:,:,142) = 2.4e-12_r8
      rate(:,:,146) = 1.4e-11_r8
      rate(:,:,153) = 2.4e-12_r8
      rate(:,:,156) = 1.4e-11_r8
      rate(:,:,159) = 5.e-12_r8
      rate(:,:,172) = 7.e-13_r8
      rate(:,:,175) = 2.4e-12_r8
      rate(:,:,179) = 4.5e-11_r8
      rate(:,:,183) = 2.4e-12_r8
      rate(:,:,192) = 4.e-14_r8
      rate(:,:,193) = 3.e-12_r8
      rate(:,:,194) = 1.e-11_r8
      rate(:,:,195) = 2.1e-6_r8
      rate(:,:,196) = 7.1e-6_r8
      rate(:,:,202) = 7.1e-6_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,42) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,43) = 2.1e-11_r8 * exp( 115._r8 * itemp(:,:) )
      rate(:,:,44) = 3.2e-11_r8 * exp( 70._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,47) = 5.5e-12_r8 * exp_fac(:,:)
      rate(:,:,140) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,48) = 2.2e-11_r8 * exp( 120._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,49) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,81) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,95) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,101) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,115) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,120) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,126) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,137) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,144) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,157) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,171) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,50) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,51) = 1.e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,53) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,54) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,59) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,55) = 4.2e-12_r8 * exp( -240._r8 * itemp(:,:) )
      rate(:,:,60) = 3e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,61) = 5.6e-12_r8 * exp_fac(:,:)
      rate(:,:,87) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,99) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,112) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,122) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,124) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,129) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,135) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,141) = 2.2e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,62) = 1.2e-13_r8 * exp( -2450._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 170._r8 * itemp(:,:) )
      rate(:,:,63) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,68) = 1.5e-11_r8 * exp_fac(:,:)
      rate(:,:,70) = 1.3e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,75) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,77) = 2.8e-12_r8 * exp_fac(:,:)
      rate(:,:,117) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,78) = 5.e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,79) = 1.9e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,80) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,82) = 6.0e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,86) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,89) = 1.6e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,90) = 8.7e-12_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,91) = 2.6e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,92) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,100) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,113) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,118) = 8.6e-13_r8 * exp_fac(:,:)
      rate(:,:,125) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,130) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,136) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,143) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,154) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,170) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,176) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,184) = 8.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,97) = 6.5e-15_r8 * exp_fac(:,:)
      rate(:,:,103) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,98) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,102) = 5.6e-12_r8 * exp_fac(:,:)
      rate(:,:,104) = 8.1e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,106) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,160) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,107) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,110) = 2.5e-12_r8 * exp_fac(:,:)
      rate(:,:,119) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,161) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,111) = 1.0e-11_r8 * exp( -660._r8 * itemp(:,:) )
      rate(:,:,114) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,128) = 2.3e-12_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,132) = 1.7e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,139) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,145) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,155) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,185) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,147) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,148) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,149) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,181) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,150) = 4.4e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,151) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,152) = 1.3e-13_r8 * exp_fac(:,:)
      rate(:,:,158) = 5.3e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,182) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,162) = 4.6e-12_r8 * exp_fac(:,:)
      rate(:,:,163) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,166) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,167) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,168) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,173) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,177) = 8.4e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,178) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,186) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,187) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,189) = 7.3e-12_r8 * exp( -620._r8 * itemp(:,:) )
      rate(:,:,190) = 6.9e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,198) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,200) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,201) = 1.7e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 6.9e-31_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,56), m, .6_r8, ko, kinf, n )

      ko(:,:) = 2.e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**.7_r8
      call jpl( rate(1,1,64), m, .6_r8, ko, kinf, n )

      ko(:,:) = 2.0e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.5e-11_r8
      call jpl( rate(1,1,66), m, .6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 4.7e-12_r8 * itemp(:,:)**1.4_r8
      call jpl( rate(1,1,69), m, .6_r8, ko, kinf, n )

      ko(:,:) = 1.e-28_r8 * itemp(:,:)**.8_r8
      kinf(:,:) = 8.8e-12_r8
      call jpl( rate(1,1,85), m, .6_r8, ko, kinf, n )

      ko(:,:) = 8.e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.e-11_r8
      call jpl( rate(1,1,96), m, .5_r8, ko, kinf, n )

      ko(:,:) = 8.5e-29_r8 * itemp(:,:)**6.5_r8
      kinf(:,:) = 1.1e-11_r8 * itemp(:,:)
      call jpl( rate(1,1,105), m, .6_r8, ko, kinf, n )

      ko(:,:) = 8.e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.e-11_r8
      call jpl( rate(1,1,191), m, .5_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
