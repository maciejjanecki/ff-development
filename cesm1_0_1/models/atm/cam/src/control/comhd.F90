
module comhd
!----------------------------------------------------------------------- 
! 
! Purpose: horizontal diffusion module
!
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use infnan,       only: inf
   use pspect, only: pnmax
   use pmgrid, only: plev
   use abortutils, only: endrun

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public comhd_defaultopts     ! get default runtime options
   public comhd_setopts         ! set runtime options

!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
!
   real(r8), public :: cnfac  ! Courant num factor(multiply by max |V|)
   real(r8), public :: cnlim  ! Maximum allowable courant number
   real(r8), public :: dif2 = inf      ! Del^2 diffusion coeff. (also namelist 
                                       ! variable)
   real(r8), public :: dif4 = inf      ! Del^4 diffusion coeff. (also namelist 
                                       ! variable)
   real(r8), public :: hdfsd2(pnmax)   ! Del^2 mult. for each wave 
                                       ! (vorticity-divergence)
   real(r8), public :: hdfst2(pnmax)   ! Del^2 multiplier for each wave (t-q)
   real(r8), public :: hdfsd4(pnmax)   ! Del^4 mult. for each wave 
                                       ! (vorticity-divergence)
   real(r8), public :: hdfst4(pnmax)   ! Del^4 multiplier for each wave (t-q)
   real(r8), public :: hdiftq(pnmax,plev)   ! Temperature-tracer diffusion 
                                            ! factors
   real(r8), public :: hdifzd(pnmax,plev)   ! Vorticity-divergence diffusion 
                                            ! factors
!
   integer, public :: kmnhd4    ! Top level for del^4 diffusion
   integer, public :: kmxhd2    ! Bottom level for increased del^2 diffusion
   ! The Courant limiter is only used in eul dycore.
   ! Number of levels to apply Courant limiter (also in namelist)
   integer, private, parameter :: def_kmxhdc = 5    ! default
   integer, public :: kmxhdc = def_kmxhdc

   integer, public :: nindex(plev)    ! Starting index for spectral truncation
   integer, public :: nmaxhd    ! Maximum two dimensional wave number
!

contains


   subroutine comhd_defaultopts(dif2_out, &
                                dif4_out, &
                                kmxhdc_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
     use dycore,       only: dycore_is, get_resolution
!------------------------------Arguments--------------------------------
     ! Del^2 diffusion coeff.
     real(r8), intent(out), optional :: dif2_out
     ! Del^4 diffusion coeff.
     real(r8), intent(out), optional :: dif4_out
     ! Number of levels to apply Courant limiter
     integer, intent(out), optional :: kmxhdc_out
!-----------------------------------------------------------------------
!
! Horizontal diffusion is used in both eul and sld dycores, but for sld
! the parameters are zero in the production version.
! Horizontal diffusion is not used in the fv dycore.
!
! Set dif2 and dif4 values for known resolutions, otherwise require in 
! namelist
!
     if ( present(dif2_out) ) then
       if (dycore_is ('EUL')) then
          select case ( get_resolution() )
          case ( 'T5' )
             dif2_out = 2.5e7_r8
          case ( 'T21' )
             dif2_out = 2.5e5_r8
          case ( 'T31' )
             dif2_out = 2.5e5_r8
          case ( 'T42' )
             dif2_out = 2.5e5_r8
          case ( 'T85' )
             dif2_out = 2.5e5_r8
          case ( 'T170' )
             dif2_out = 2.5e5_r8
          case default
             dif2_out = 2.5e5_r8
          end select
       else
          dif2_out = 0._r8
       endif
     endif
     if ( present(dif4_out) ) then
       if (dycore_is ('EUL')) then
          select case ( get_resolution() )
          case ( 'T5' )
             dif4_out = 1.0e18_r8
          case ( 'T21' )
             dif4_out = 2.0e16_r8
          case ( 'T31' )
             dif4_out = 2.0e16_r8
          case ( 'T42' )
             dif4_out = 1.0e16_r8
          case ( 'T85' )
             dif4_out = 1.0e15_r8
          case ( 'T170' )
             dif4_out = 1.5e14_r8
          case default
             dif4_out = inf
          end select
       else
          dif4_out = 0._r8
       endif
     endif
     if ( present(kmxhdc_out) ) then
       kmxhdc_out = def_kmxhdc
     endif
   end subroutine comhd_defaultopts


   subroutine comhd_setopts(dif2_in, &
                            dif4_in, &
                            kmxhdc_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     ! Del^2 diffusion coeff.
     real(r8), intent(in), optional :: dif2_in
     ! Del^4 diffusion coeff.
     real(r8), intent(in), optional :: dif4_in
     ! Number of levels to apply Courant limiter
     integer, intent(in), optional :: kmxhdc_in
!-----------------------------------------------------------------------
!
! Horizontal diffusion is used in both eul and sld dycores, but for sld
! the parameters are zero in the production version.
! Horizontal diffusion is not used in the fv dycore.
!
! Set dif2 and dif4 values for known resolutions, otherwise require in 
! namelist
!
     if ( present(dif2_in) ) then
       dif2 = dif2_in
     endif
     if ( present(dif4_in) ) then
       dif4 = dif4_in
       if ( dif4 == inf ) then
          call endrun ('COMHD_SETOPTS: dif4 must be set in the namelist for this resolution.')
       end if
     endif
     if ( present(kmxhdc_in) ) then
       kmxhdc = kmxhdc_in
       if (kmxhdc >= plev .or. kmxhdc < 0) then
          call endrun ('COMHD_SETOPTS:  ERROR:  KMXHDC must be between 0 and plev-1')
       end if
     endif
   end subroutine comhd_setopts


end module comhd

