
module mo_setinv

  use cam_logfile, only: iulog

  implicit none

  save

  integer :: id_o, id_o2, id_h
  integer :: m_ndx, o2_ndx, n2_ndx, h2o_ndx, o3_ndx
  logical :: has_o2, has_n2, has_h2o, has_o3, has_var_o2

  private
  public :: setinv_inti, setinv, has_h2o, o2_ndx, h2o_ndx, n2_ndx

contains

  subroutine setinv_inti
    !-----------------------------------------------------------------
    !        ... initialize the module
    !-----------------------------------------------------------------

    use mo_chem_utls, only : get_inv_ndx, get_spc_ndx
    use spmd_utils,   only : masterproc

    implicit none

    m_ndx   = get_inv_ndx( 'M' )
    n2_ndx  = get_inv_ndx( 'N2' )
    o2_ndx  = get_inv_ndx( 'O2' )
    h2o_ndx = get_inv_ndx( 'H2O' )
    o3_ndx  = get_inv_ndx( 'O3' )

    id_o  = get_spc_ndx('O')
    id_o2 = get_spc_ndx('O2')
    id_h  = get_spc_ndx('H')

    has_var_o2 = id_o2>0 .and. id_o>0 .and. id_h>0

    has_n2  = n2_ndx > 0
    has_o2  = o2_ndx > 0
    has_h2o = h2o_ndx > 0
    has_o3  = o3_ndx > 0

    if (masterproc) write(iulog,*) 'setinv_inti: m,n2,o2,h2o ndx = ',m_ndx,n2_ndx,o2_ndx,h2o_ndx

  end subroutine setinv_inti

  subroutine setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk )
    !-----------------------------------------------------------------
    !        ... set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use shr_kind_mod,  only : r8 => shr_kind_r8
    use ppgrid,        only : pcols, pver
    use chem_mods,     only : nfs, gas_pcnst
    use mo_constants,  only : boltz=>boltzmann
    use tracer_cnst,   only : num_tracer_cnst, tracer_cnst_flds, get_cnst_data
    use mo_chem_utls,  only : get_inv_ndx

    implicit none

!!$      real(r8), parameter ::  boltz   = 1.38044e-16_r8         ! erg/K

    !-----------------------------------------------------------------
    !        ... dummy arguments
    !-----------------------------------------------------------------
    real(r8), intent(in)  ::      tfld(pcols,pver)          ! temperature
    real(r8), intent(in)  ::      h2ovmr(ncol,pver)         ! water vapor vmr
    real(r8), intent(in)  ::      pmid(pcols,pver)          ! pressure
    integer,  intent(in)  ::      ncol                      ! chunk column count
    integer,  intent(in)  ::      lchnk                     ! chunk number
    real(r8), intent(in)  ::      vmr(ncol,pver,gas_pcnst)  ! vmr
    real(r8), intent(out) ::      invariants(ncol,pver,nfs) ! invariant array

    real(r8) :: cnst_offline( ncol, pver )

    !-----------------------------------------------------------------
    !        .. local variables
    !-----------------------------------------------------------------
    integer :: k, i, ndx
    real(r8), parameter ::  Pa_xfac = 10._r8                 ! Pascals to dyne/cm^2
    real(r8) :: sum1(ncol)

    !-----------------------------------------------------------------
    !        note: invariants are in cgs density units.
    !              the pmid array is in pascals and must be
    !	       mutiplied by 10. to yield dynes/cm**2.
    !-----------------------------------------------------------------
    invariants(:,:,:) = 0._r8
    !-----------------------------------------------------------------
    !	... set m, n2, o2, and h2o densities
    !-----------------------------------------------------------------
    do k = 1,pver
       invariants(:ncol,k,m_ndx) = Pa_xfac * pmid(:ncol,k) / (boltz*tfld(:ncol,k))
    end do
    if( has_n2 ) then
       if ( has_var_o2 ) then
          do k = 1,pver
             sum1(:ncol) = (vmr(:ncol,k,id_o) + vmr(:ncol,k,id_o2) + vmr(:ncol,k,id_h))
             invariants(:ncol,k,n2_ndx) = (1._r8 - sum1(:)) * invariants(:ncol,k,m_ndx)
          end do
       else
          do k = 1,pver
             invariants(:ncol,k,n2_ndx) = .79_r8 * invariants(:ncol,k,m_ndx)
          end do
       endif
    end if
    if( has_o2 ) then
       do k = 1,pver
          invariants(:ncol,k,o2_ndx) = .21_r8 * invariants(:ncol,k,m_ndx)
       end do
    end if
    if( has_h2o ) then
       do k = 1,pver
          invariants(:ncol,k,h2o_ndx) = h2ovmr(:ncol,k) * invariants(:ncol,k,m_ndx)
       end do
    end if

    do i = 1,num_tracer_cnst

       call get_cnst_data( tracer_cnst_flds(i), cnst_offline,  ncol, lchnk )
       ndx =  get_inv_ndx( tracer_cnst_flds(i) )

       do k = 1,pver
          invariants(:ncol,k,ndx) = cnst_offline(:ncol,k)*invariants(:ncol,k,m_ndx)
       enddo

    enddo

  end subroutine setinv

end module mo_setinv
