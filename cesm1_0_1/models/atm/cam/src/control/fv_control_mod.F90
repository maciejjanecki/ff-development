module fv_control_mod
  use shr_kind_mod, only : r8=> shr_kind_r8

  public
  real(r8) :: tmass0
  real(r8) :: zgsint
  ! f-v dynamics specific
  ! _ord = 1: first order upwind
  ! _ord = 2: 2nd order van Leer (Lin et al 1994)
  ! _ord = 3: standard PPM 
  ! _ord = 4: enhanced PPM (default)

  integer :: nsplit = 0                  ! Lagrangian time splits (Lin-Rood only)
  integer :: nspltrac = 0                ! Tracer time splits (Lin-Rood only)
  integer :: nspltvrm = 0                ! Vertical re-mapping time splits (Lin-Rood only)
  integer :: iord = 4                    ! scheme to be used in E-W direction
  integer :: jord = 4                    ! scheme to be used in N-S direction
  integer :: kord = 4                    ! scheme to be used for vertical mapping
  logical :: dyn_conservative = .false.  ! Flag indicating whether the dynamics is conservative
  integer :: filtcw = 0                  ! flag for filtering c-grid winds
  integer :: ct_overlap = 0              ! nonzero for overlap of cd_core and trac2d, 0 otherwise
  integer :: trac_decomp = 1             ! size of tracer domain decomposition for trac2d
  integer :: fft_flt = 1                 ! 0 => FFT/algebraic filter; 1 => FFT filter
  integer :: div24del2flag = 2           ! 2 for 2nd order div damping, 4 for 4th order div damping,
                                         ! 42 for 4th order div damping plus 2nd order velocity damping
  real(r8):: del2coef = 3.e5_r8          ! strength of 2nd order velocity damping

end module fv_control_mod
