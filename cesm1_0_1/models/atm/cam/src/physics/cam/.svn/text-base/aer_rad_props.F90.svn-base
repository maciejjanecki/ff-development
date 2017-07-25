module aer_rad_props

!------------------------------------------------------------------------------------------------
! Converts aerosol masses to bulk optical properties for sw and lw radiation
! computations.  
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: rga
use physics_types,    only: physics_state
use phys_buffer,      only: pbuf_size_max, pbuf_fld
use radconstants,     only: nrh, nswbands, nlwbands, idx_sw_diag, ot_length
use rad_constituents, only: rad_cnst_get_clim_info, rad_cnst_get_clim_aer, &
                            rad_cnst_get_clim_aer_props
use wv_saturation,    only: aqsat
use cam_history,      only: fieldname_len, addfld, phys_decomp, outfld, fillvalue, add_default

use abortutils,       only: endrun
use constituents,     only: pcnst

implicit none
private
save

public :: &
   aer_rad_props_init,        &
   aer_rad_props_get_clim_sw, & ! return SW optical props of total bulk aerosols
   aer_rad_props_get_clim_lw    ! return LW optical props of total bulk aerosols

! Private data
character(len=fieldname_len), pointer :: odv_names(:)  ! outfld names for visible OD


!==============================================================================
contains
!==============================================================================

subroutine aer_rad_props_init()


   integer                    :: i
   integer                    :: numaerosols  ! number of aerosols
   character(len=64), pointer :: aernames(:)  ! aerosol names
   !----------------------------------------------------------------------------


   call addfld ('AEROD_v ', '1', 1, 'A', &
      'Total Aerosol Optical Depth in visible band', phys_decomp, flag_xyfill=.true.)
   call add_default ('AEROD_v', 1, ' ')

   ! Contributions to AEROD_v from individual aerosols (climate species).

   ! number of aerosols in climate list
   call rad_cnst_get_clim_info(naero=numaerosols)

   ! get names of aerosols
   allocate(aernames(numaerosols))
   call rad_cnst_get_clim_info(aernames=aernames)

   ! create outfld names for visible OD and call addfld
   allocate(odv_names(numaerosols))
   do i = 1, numaerosols
      odv_names(i) = 'ODV_'//trim(aernames(i))
      call addfld (odv_names(i), '1', 1, 'A', &
         trim(aernames(i))//' optical depth in visible band', phys_decomp, flag_xyfill=.true.)
      call add_default (odv_names(i), 1, ' ')
   end do
   deallocate(aernames)

end subroutine aer_rad_props_init

!==============================================================================

subroutine aer_rad_props_get_clim_sw(state, pbuf, nnite, idxnite, &
                                     tau, tau_w, tau_w_g, tau_w_f,&
                                     diagnosticindex)

! return totaled (across all species) layer tau, omega, g, f 
! for all spectral interval for aerosols affecting the climate

   ! this use should be removed when the relative humidity calc is fixed
   use physconst, only: epsilo
#ifdef MODAL_AERO
   use modal_aero_data,   only: maxd_amode
   use modal_aer_opt,     only: modal_size_parameters, modal_aero_sw
   use modal_aer_opt,     only: ncoef
   use phys_buffer,       only: pbuf_get_fld_idx
#endif

   ! Arguments
   type(physics_state), intent(in) :: state
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   integer,             intent(in) :: nnite                ! number of night columns
   integer,             intent(in) :: idxnite(:)           ! local column indices of night columns
   integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information

   real(r8), intent(out) :: tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8), intent(out) :: tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8), intent(out) :: tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * tau * w
   real(r8), intent(out) :: tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * tau * w

   ! Local variables

   integer :: ncol
   integer :: lchnk
   integer :: k, i    ! lev and daycolumn indices
   integer :: iswband ! sw band indices

   ! optical props for each aerosol
   ! hygroscopic
   real(r8), pointer :: h_ext(:,:)
   real(r8), pointer :: h_ssa(:,:)
   real(r8), pointer :: h_asm(:,:)
   ! non-hygroscopic
   real(r8), pointer :: n_ext(:)
   real(r8), pointer :: n_ssa(:)
   real(r8), pointer :: n_asm(:)
   real(r8), pointer :: n_scat(:)
   real(r8), pointer :: n_ascat(:)
   ! radius-dependent
   real(r8), pointer :: r_ext(:,:)    ! radius-dependent mass-specific extinction
   real(r8), pointer :: r_scat(:,:)
   real(r8), pointer :: r_ascat(:,:)
   real(r8), pointer :: r_mu(:)       ! log(radius) domain variable for r_ext, r_scat, r_ascat

   ! radiative properties for each aerosol
   real(r8) :: ta (pcols,pver,nswbands)
   real(r8) :: tw (pcols,pver,nswbands)
   real(r8) :: twf(pcols,pver,nswbands)
   real(r8) :: twg(pcols,pver,nswbands)

   ! aerosol masses
   real(r8), pointer :: aermmr(:,:)    ! mass mixing ratio of aerosols
   real(r8) :: mmr_to_mass(pcols,pver) ! conversion factor for mmr to mass
   real(r8) :: aermass(pcols,pver)     ! mass of aerosols
#ifdef MODAL_AERO
   real(r8), pointer, dimension(:,:,:) :: dgnumwet ! number mode diameter
   real(r8), pointer, dimension(:,:,:) :: qaerwat ! aerosol water
!   real(r8) :: raer(pcols, pver, pcnst)
   real(r8) :: radsurf(pcols,pver,maxd_amode) ! aerosol surface mode radius
   real(r8) :: logradsurf(pcols,pver,maxd_amode) ! log(aerosol surface mode radius)
   real(r8) :: cheb(ncoef,maxd_amode,pcols,pver)
   integer :: ifld
#endif

   ! for table lookup into rh grid
   real(r8) :: esat(pcols,pver)     ! saturation vapor pressure
   real(r8) :: qsat(pcols,pver)     ! saturation specific humidity
   real(r8) :: rh(pcols,pver)
   real(r8) :: rhtrunc(pcols,pver)
   real(r8) :: wrh(pcols,pver)
   integer  :: krh(pcols,pver)
 
   integer  :: iaerosol        ! aerosol index
   integer  :: numaerosols     ! number of aerosols in climate list
   integer  :: diag_idx = 0    ! idx for climate call
   character(len=ot_length) :: opticstype       ! hygro or nonhygro
   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   ! compute mixing ratio to mass conversion
   do k = 1, pver
      mmr_to_mass(:ncol,k) = rga * state%pdeldry(:ncol,k)
   enddo

   ! initialize to conditions that would cause failure
   tau     (:,:,:) = -100._r8
   tau_w   (:,:,:) = -100._r8
   tau_w_g (:,:,:) = -100._r8
   tau_w_f (:,:,:) = -100._r8

   ! top layer (ilev = 0) has no aerosol (ie tau = 0)
   ! also initialize rest of layers to accumulate od's
   tau    (1:ncol,:,:) = 0._r8
   tau_w  (1:ncol,:,:) = 0._r8
   tau_w_g(1:ncol,:,:) = 0._r8
   tau_w_f(1:ncol,:,:) = 0._r8

   ! calculate relative humidity for table lookup into rh grid
   call aqsat(state%t, state%pmid, esat, qsat, pcols, &
               ncol, pver, 1, pver)
   rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver)

   rhtrunc(1:ncol,1:pver) = min(rh(1:ncol,1:pver),1._r8)
   krh(1:ncol,1:pver) = min(floor( rhtrunc(1:ncol,1:pver) * nrh ) + 1, nrh - 1) ! index into rh mesh
   wrh(1:ncol,1:pver) = rhtrunc(1:ncol,1:pver) * nrh - krh(1:ncol,1:pver)       ! (-) weighting on left side values

   ! number of aerosols in climate list
   if (present(diagnosticindex)) then
       call rad_cnst_get_clim_info(naero=numaerosols, diagnosticindex=diagnosticindex)
       diag_idx = diagnosticindex
   else
       call rad_cnst_get_clim_info(naero=numaerosols)
       diag_idx = 0
   endif

#ifdef MODAL_AERO
   if (present(diagnosticindex)) then
      tau    (1:ncol,:,:) = 0._r8
      tau_w  (1:ncol,:,:) = 0._r8
      tau_w_g(1:ncol,:,:) = 0._r8
      tau_w_f(1:ncol,:,:) = 0._r8
   ! Loop over aerosols in climate list
      do iaerosol = 1, numaerosols
      ! get aerosol type
         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype, diagnosticindex=diagnosticindex)
         if (trim(opticstype) == 'volcanic_radius') then

      ! get aerosol mass mixing ratio

            call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr, diagnosticindex=diagnosticindex)
            aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

     ! get optical properties for volcanic aerosols

            call rad_cnst_get_clim_aer_props(iaerosol, r_sw_ext=r_ext, r_sw_scat=r_scat, r_sw_ascat=r_ascat, mu=r_mu, diagnosticindex=diagnosticindex)
            call get_volcanic_radius_rad_props(lchnk, ncol, aermass, pbuf, r_ext, r_scat, r_ascat, r_mu, ta, tw, twg, twf)
            tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
            tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
            tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
            tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)
            call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaerosol, ta(:,:,idx_sw_diag), diag_idx)
         end if
      end do
   else
      ifld = pbuf_get_fld_idx('DGNUMWET')
      dgnumwet => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:maxd_amode)
      ifld  = pbuf_get_fld_idx( 'QAERWAT' )
      if ( associated(pbuf(ifld)%fld_ptr) ) then
         qaerwat => pbuf(ifld)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:maxd_amode )
      else
         call endrun( 'pbuf for QAERWAT not allocated in modal_size_parameters' )
      end if

      call modal_size_parameters(ncol, dgnumwet, radsurf, logradsurf, cheb)
!      raer(:,:,:)=state%q(:,:,:) ! for starters. later use rad_cnst_get_clim_aer
      do iswband=1,nswbands
         call modal_aero_sw(ncol, lchnk, iswband, nnite, idxnite, tau(1,0,iswband), tau_w(1,0,iswband), &
                  tau_w_g(1,0,iswband), tau_w_f(1,0,iswband), state%q, state, qaerwat, radsurf, logradsurf, cheb )
      end do

   ! Loop over aerosols in climate list
      do iaerosol = 1, numaerosols

      ! get aerosol type

         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype)

         if (trim(opticstype) == 'volcanic_radius') then

      ! get aerosol mass mixing ratio

            call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr)
            aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

     ! get optical properties for volcanic aerosols
            call rad_cnst_get_clim_aer_props(iaerosol, r_sw_ext=r_ext, r_sw_scat=r_scat, r_sw_ascat=r_ascat, mu=r_mu )
            call get_volcanic_radius_rad_props(lchnk, ncol, aermass, pbuf, r_ext, r_scat, r_ascat, r_mu, ta, tw, twg, twf)
            tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
            tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
            tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
            tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)
         ! diagnostic output of individual aerosol optical properties
         ! currently implemented for climate list only
            call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaerosol, ta(:,:,idx_sw_diag), diag_idx)
         end if
      end do
   end if
#else

   ! Loop over aerosols in climate list
   do iaerosol = 1, numaerosols

      ! get aerosol mass mixing ratio
      if (present(diagnosticindex)) then
          call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr, diagnosticindex=diagnosticindex)
      else
          call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr)
      endif
      aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

      ! get aerosol type
      if (present(diagnosticindex)) then
          call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype, diagnosticindex=diagnosticindex)
      else
          call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype)
      endif

      select case (trim(opticstype))
      case('hygro')
         ! get optical properties for hygroscopic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, sw_hygro_ext=h_ext, sw_hygro_ssa=h_ssa, sw_hygro_asm=h_asm, &
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, sw_hygro_ext=h_ext, sw_hygro_ssa=h_ssa, sw_hygro_asm=h_asm)
         endif
         call get_hygro_rad_props(ncol, krh, wrh, aermass, h_ext, h_ssa, h_asm, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('hygroscopic','hygroscopi')
         ! get optical properties for hygroscopic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, sw_hygro_ext=h_ext, sw_hygro_ssa=h_ssa, sw_hygro_asm=h_asm, &
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, sw_hygro_ext=h_ext, sw_hygro_ssa=h_ssa, sw_hygro_asm=h_asm)
         endif
         call get_hygro_rad_props(ncol, krh, wrh, aermass, h_ext, h_ssa, h_asm, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('nonhygro')
         ! get optical properties for non-hygroscopic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_ssa=n_ssa, sw_nonhygro_asm=n_asm, &
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_ssa=n_ssa, sw_nonhygro_asm=n_asm)
         endif

         call get_nonhygro_rad_props(ncol, aermass, n_ext, n_ssa, n_asm, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('insoluble ')
         ! get optical properties for non-hygroscopic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_ssa=n_ssa, sw_nonhygro_asm=n_asm, &
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_ssa=n_ssa, sw_nonhygro_asm=n_asm)
         endif

         call get_nonhygro_rad_props(ncol, aermass, n_ext, n_ssa, n_asm, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('volcanic')
         ! get optical properties for volcanic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_scat=n_scat, sw_nonhygro_ascat=n_ascat,&
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_scat=n_scat, sw_nonhygro_ascat=n_ascat)
         endif

         call get_volcanic_rad_props(ncol, aermass, n_ext, n_scat, n_ascat, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('volcanic_radius')
         ! get optical properties for volcanic aerosols
         if (present(diagnosticindex)) then
             call rad_cnst_get_clim_aer_props(iaerosol, r_sw_ext=r_ext, r_sw_scat=r_scat, r_sw_ascat=r_ascat, mu=r_mu,&
                                              diagnosticindex=diagnosticindex)
         else
             call rad_cnst_get_clim_aer_props(iaerosol, r_sw_ext=r_ext, r_sw_scat=r_scat, r_sw_ascat=r_ascat, mu=r_mu )
         endif

         call get_volcanic_radius_rad_props(lchnk, ncol, aermass, pbuf, r_ext, r_scat, r_ascat, r_mu, ta, tw, twg, twf)
         tau    (1:ncol,1:pver,:) = tau    (1:ncol,1:pver,:) + ta (1:ncol,:,:)
         tau_w  (1:ncol,1:pver,:) = tau_w  (1:ncol,1:pver,:) + tw (1:ncol,:,:)
         tau_w_g(1:ncol,1:pver,:) = tau_w_g(1:ncol,1:pver,:) + twg(1:ncol,:,:)
         tau_w_f(1:ncol,1:pver,:) = tau_w_f(1:ncol,1:pver,:) + twf(1:ncol,:,:)

      case('zero')
         ! no effect of "zero" aerosols, so update nothing
      case default
         call endrun('aer_rad_props_get_clim_sw: unsupported opticstype :'//trim(opticstype)//':')
      end select

      ! diagnostic output of individual aerosol optical properties
      ! currently implemented for climate list only
      call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaerosol, ta(:,:,idx_sw_diag), diag_idx)

   enddo
#endif

   ! diagnostic output of total aerosol optical properties
   ! currently implemented for climate list only
   call aer_vis_diag_out(lchnk, ncol, nnite, idxnite, 0, tau(:,:,idx_sw_diag), diag_idx)

end subroutine aer_rad_props_get_clim_sw

!==============================================================================

subroutine aer_rad_props_get_clim_lw(state, pbuf, odap_aer, diagnosticindex)

   use physconst,     only: epsilo
   use radconstants,  only: ot_length
   use phys_buffer,   only: pbuf_get_fld_idx
#ifdef MODAL_AERO
   use modal_aer_opt, only: modal_aero_lw
#endif

! Purpose: Compute aerosol transmissions needed in absorptivity/
!    emissivity calculations

! lw extinction is the same representation for all 
! species.  If this changes, this routine will need to do something
! similar to the sw with routines like get_hygro_lw_abs

   ! Arguments
   type(physics_state), intent(in)  :: state
   type(pbuf_fld),      intent(in)  :: pbuf(pbuf_size_max)
   real(r8),            intent(out) :: odap_aer(pcols,pver,nlwbands) ! [fraction] absorption optical depth, per layer
   integer, optional,   intent(in)  :: diagnosticindex

   ! Local variables

   integer :: bnd_idx     ! LW band index
   integer :: i           ! column index
   integer :: k           ! lev index
   integer :: ncol        ! number of columns
   integer :: lchnk       ! chunk index
   integer :: iaerosol    ! index into aerosol list
   integer :: numaerosols ! index into aerosol list
   character(len=ot_length) :: opticstype       ! hygro or nonhygro


   ! optical props for each aerosol
   real(r8), pointer :: lw_abs(:)
   real(r8), pointer :: lw_hygro_abs(:,:)
   real(r8), pointer :: geometric_radius(:,:)
 
   ! volcanic lookup table
   real(r8), pointer :: r_lw_abs(:,:)  ! radius dependent mass-specific absorption coefficient
   real(r8), pointer :: r_mu(:)        ! log(geometric_mean_radius) domain samples of r_lw_abs(:,:)
   integer  :: idx                     ! index to pbuf for geometric radius
   real(r8) :: mu(pcols,pver)          ! log(geometric_radius)
   real(r8) :: r_mu_min, r_mu_max, wmu, mutrunc
   integer  :: nmu, kmu

   ! for table lookup into rh grid
   real(r8) :: esat(pcols,pver)     ! saturation vapor pressure
   real(r8) :: qsat(pcols,pver)     ! saturation specific humidity
   real(r8) :: rh(pcols,pver)
   real(r8) :: rhtrunc(pcols,pver)
   real(r8) :: wrh(pcols,pver)
   integer  :: krh(pcols,pver)

   ! aerosol (vertical) mass path and extinction
   ! aerosol masses
   real(r8), pointer :: aermmr(:,:)    ! mass mixing ratio of aerosols
   real(r8) :: mmr_to_mass(pcols,pver) ! conversion factor for mmr to mass
   real(r8) :: aermass(pcols,pver)     ! mass of aerosols
#ifdef MODAL_AERO
!   real(r8) :: raer(pcols, pver, pcnst)
#endif
   !-----------------------------------------------------------------------------

   ncol = state%ncol
   lchnk = state%lchnk

#ifdef MODAL_AERO
   if (present(diagnosticindex)) then
      odap_aer = 0._r8
! hack to get volcanic aerosols into mam
! number of aerosols in climate list
      call rad_cnst_get_clim_info(naero=numaerosols, diagnosticindex=diagnosticindex)

! Loop over aerosols in climate list

      do iaerosol = 1, numaerosols

         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype, diagnosticindex=diagnosticindex)

         if (trim(opticstype) == 'volcanic_radius') then

! compute mixing ratio to mass conversion
            do k = 1, pver
               mmr_to_mass(:ncol,k) = rga * state%pdeldry(:ncol,k)
            end do
            
! compute optical depths odap_aer (summed over all aerosols)
! total optical depths from top (odap_aer_ttl)

! calculate relative humidity for table lookup into rh grid
            call aqsat(state%t, state%pmid, esat, qsat, pcols, &
                 ncol, pver, 1, pver)
            rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver)
            
            rhtrunc(1:ncol,1:pver) = min(rh(1:ncol,1:pver),1._r8)
            krh(1:ncol,1:pver) = min(floor( rhtrunc(1:ncol,1:pver) * nrh ) + 1, nrh - 1) ! index into rh mesh
            wrh(1:ncol,1:pver) = rhtrunc(1:ncol,1:pver) * nrh - krh(1:ncol,1:pver)       ! (-) weighting on left side values

      ! get aerosol mass mixing ratio

            call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr, diagnosticindex=diagnosticindex)
            aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

          ! get optical properties for hygroscopic aerosols
            call rad_cnst_get_clim_aer_props(iaerosol, r_lw_abs=r_lw_abs, mu=r_mu, diagnosticindex=diagnosticindex)

          ! get microphysical properties for volcanic aerosols
            idx = pbuf_get_fld_idx('VOLC_RAD_GEOM')
            geometric_radius => pbuf(idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1)

          ! interpolate in radius
          ! caution: clip the table with no warning when outside bounds
            nmu = size(r_mu)
            r_mu_max = r_mu(nmu)
            r_mu_min = r_mu(1)
            do i = 1, ncol
               do k = 1, pver
                  if(geometric_radius(i,k) > 0.) then
                     mu(i,k) = log(geometric_radius(i,k))
                  else
                     mu(i,k) = 0.
                  endif
                  mutrunc = max(min(mu(i,k),r_mu_max),r_mu_min)
                  kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
                  wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
                  do bnd_idx = 1, nlwbands
                     odap_aer(i,k,bnd_idx) = odap_aer(i,k,bnd_idx) + &
                          aermass(i,k) * &
                          ((1._r8 - wmu) * r_lw_abs(bnd_idx, kmu  ) + &
                          (wmu) * r_lw_abs(bnd_idx, kmu+1))
                  end do
               enddo
            enddo
         end if
      enddo
   else
!      raer(:,:,:)=state%q(:,:,:) ! for starters. later use rad_cnst_get_clim_aer
      call modal_aero_lw(ncol,lchnk, pbuf, state%q, odap_aer, state%pdeldry)

! hack to get volcanic aerosols into mam
! number of aerosols in climate list
      call rad_cnst_get_clim_info(naero=numaerosols)

! Loop over aerosols in climate list

      do iaerosol = 1, numaerosols

         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype)

         if (trim(opticstype) == 'volcanic_radius') then

! compute mixing ratio to mass conversion
            do k = 1, pver
               mmr_to_mass(:ncol,k) = rga * state%pdeldry(:ncol,k)
            end do
            
! compute optical depths odap_aer (summed over all aerosols)
! total optical depths from top (odap_aer_ttl)

! calculate relative humidity for table lookup into rh grid
            call aqsat(state%t, state%pmid, esat, qsat, pcols, &
                 ncol, pver, 1, pver)
            rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver)
            
            rhtrunc(1:ncol,1:pver) = min(rh(1:ncol,1:pver),1._r8)
            krh(1:ncol,1:pver) = min(floor( rhtrunc(1:ncol,1:pver) * nrh ) + 1, nrh - 1) ! index into rh mesh
            wrh(1:ncol,1:pver) = rhtrunc(1:ncol,1:pver) * nrh - krh(1:ncol,1:pver)       ! (-) weighting on left side values

      ! get aerosol mass mixing ratio

            call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr)
            aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

          ! get optical properties for hygroscopic aerosols
            call rad_cnst_get_clim_aer_props(iaerosol, r_lw_abs=r_lw_abs, mu=r_mu)

          ! get microphysical properties for volcanic aerosols
            idx = pbuf_get_fld_idx('VOLC_RAD_GEOM')
            geometric_radius => pbuf(idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1)

          ! interpolate in radius
          ! caution: clip the table with no warning when outside bounds
            nmu = size(r_mu)
            r_mu_max = r_mu(nmu)
            r_mu_min = r_mu(1)
            do i = 1, ncol
               do k = 1, pver
                  if(geometric_radius(i,k) > 0.) then
                     mu(i,k) = log(geometric_radius(i,k))
                  else
                     mu(i,k) = 0.
                  endif
                  mutrunc = max(min(mu(i,k),r_mu_max),r_mu_min)
                  kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
                  wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
                  do bnd_idx = 1, nlwbands
                     odap_aer(i,k,bnd_idx) = odap_aer(i,k,bnd_idx) + &
                          aermass(i,k) * &
                          ((1._r8 - wmu) * r_lw_abs(bnd_idx, kmu  ) + &
                          (wmu) * r_lw_abs(bnd_idx, kmu+1))
                  end do
               enddo
            enddo
         end if
      enddo
   end if 
#else
   ! compute mixing ratio to mass conversion
   do k = 1, pver
      mmr_to_mass(:ncol,k) = rga * state%pdeldry(:ncol,k)
   end do

   ! compute optical depths odap_aer (summed over all aerosols)
   ! total optical depths from top (odap_aer_ttl)
   odap_aer = 0._r8

   ! calculate relative humidity for table lookup into rh grid
   call aqsat(state%t, state%pmid, esat, qsat, pcols, &
               ncol, pver, 1, pver)
   rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver)

   rhtrunc(1:ncol,1:pver) = min(rh(1:ncol,1:pver),1._r8)
   krh(1:ncol,1:pver) = min(floor( rhtrunc(1:ncol,1:pver) * nrh ) + 1, nrh - 1) ! index into rh mesh
   wrh(1:ncol,1:pver) = rhtrunc(1:ncol,1:pver) * nrh - krh(1:ncol,1:pver)       ! (-) weighting on left side values


   ! number of aerosols in climate list
   if (present(diagnosticindex)) then
      call rad_cnst_get_clim_info(naero=numaerosols, diagnosticindex=diagnosticindex)
   else
      call rad_cnst_get_clim_info(naero=numaerosols)
   endif

   ! Loop over aerosols in climate list
   do iaerosol = 1, numaerosols

      ! get aerosol mass mixing ratio
      if (present(diagnosticindex)) then
         call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr, diagnosticindex=diagnosticindex)
      else
         call rad_cnst_get_clim_aer(iaerosol, state, pbuf, aermmr)
      endif
      aermass(1:ncol,1:pver) = aermmr(1:ncol,1:pver) * mmr_to_mass(1:ncol,1:pver)

      if (present(diagnosticindex)) then
         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype, diagnosticindex=diagnosticindex)
      else
         call rad_cnst_get_clim_aer_props(iaerosol, opticstype=opticstype)
      endif
      select case (trim(opticstype))
      case('hygroscopic')
          ! get optical properties for hygroscopic aerosols
          if (present(diagnosticindex)) then
              call rad_cnst_get_clim_aer_props(iaerosol, lw_hygro_ext=lw_hygro_abs, diagnosticindex=diagnosticindex)
          else
              call rad_cnst_get_clim_aer_props(iaerosol, lw_hygro_ext=lw_hygro_abs)
          endif
          do bnd_idx = 1, nlwbands
             do k = 1, pver
                do i = 1, ncol
                   odap_aer(i, k, bnd_idx) = odap_aer(i, k, bnd_idx) + &
                         aermass(i, k) * &
                            ((1 + wrh(i,k)) * lw_hygro_abs(krh(i,k)+1,bnd_idx) &
                                - wrh(i,k)  * lw_hygro_abs(krh(i,k),  bnd_idx))
                enddo
             enddo
          enddo
      case('insoluble','nonhygro','hygro','volcanic')
          ! get optical properties for hygroscopic aerosols
          if (present(diagnosticindex)) then
              call rad_cnst_get_clim_aer_props(iaerosol, lw_ext=lw_abs, diagnosticindex=diagnosticindex)
          else
              call rad_cnst_get_clim_aer_props(iaerosol, lw_ext=lw_abs)
          endif
          do bnd_idx = 1, nlwbands
             do k = 1, pver          
                do i = 1, ncol
                   odap_aer(i,k,bnd_idx) = odap_aer(i,k,bnd_idx) + lw_abs(bnd_idx)*aermass(i,k)
                end do
             end do
          end do

      case('volcanic_radius')
          ! get optical properties for hygroscopic aerosols
          if (present(diagnosticindex)) then
              call rad_cnst_get_clim_aer_props(iaerosol, r_lw_abs=r_lw_abs, mu=r_mu, diagnosticindex=diagnosticindex)
          else
              call rad_cnst_get_clim_aer_props(iaerosol, r_lw_abs=r_lw_abs, mu=r_mu)
          endif
          ! get microphysical properties for volcanic aerosols
          idx = pbuf_get_fld_idx('VOLC_RAD_GEOM')
          geometric_radius => pbuf(idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1)

          ! interpolate in radius
          ! caution: clip the table with no warning when outside bounds
          nmu = size(r_mu)
          r_mu_max = r_mu(nmu)
          r_mu_min = r_mu(1)
          do i = 1, ncol
             do k = 1, pver
                if(geometric_radius(i,k) > 0.) then
                   mu(i,k) = log(geometric_radius(i,k))
                else
                   mu(i,k) = 0.
                endif
                mutrunc = max(min(mu(i,k),r_mu_max),r_mu_min)
                kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
                wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
                do bnd_idx = 1, nlwbands
                   odap_aer(i,k,bnd_idx) = odap_aer(i,k,bnd_idx) + &
                      aermass(i,k) * &
                      ((1._r8 - wmu) * r_lw_abs(bnd_idx, kmu  ) + &
                      (wmu) * r_lw_abs(bnd_idx, kmu+1))
                end do
             enddo
          enddo

      case('zero')
          ! zero aerosols types have no optical effect, so do nothing.
      case default
         call endrun('aer_rad_props_get_clim_lw: unsupported opticstype: '//trim(opticstype))
      end select
   end do
#endif

end subroutine aer_rad_props_get_clim_lw

!==============================================================================
! Private methods
!==============================================================================

subroutine get_hygro_rad_props(ncol, krh, wrh, mass, ext, ssa, asm, &
                               tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   integer,  intent(in) :: krh(pcols,pver)  ! index for linear interpolation of optics on rh
   real(r8), intent(in) :: wrh(pcols,pver)  ! weight for linear interpolation of optics on rh
   real(r8), intent(in) :: mass(pcols,pver)
   real(r8), intent(in) :: ext(:,:)
   real(r8), intent(in) :: ssa(:,:)
   real(r8), intent(in) :: asm(:,:)

   real(r8), intent(out) :: tau    (pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w  (pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w_g(pcols,pver,nswbands)
   real(r8), intent(out) :: tau_w_f(pcols,pver,nswbands)

   ! Local variables
   real(r8) :: ext1, ssa1, asm1
   integer :: icol, ilev, iswband
   !-----------------------------------------------------------------------------

   do iswband = 1, nswbands
      do icol = 1, ncol
         do ilev = 1, pver
            ext1 = (1 + wrh(icol,ilev)) * ext(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * ext(krh(icol,ilev),  iswband)
            ssa1 = (1 + wrh(icol,ilev)) * ssa(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * ssa(krh(icol,ilev),  iswband)
            asm1 = (1 + wrh(icol,ilev)) * asm(krh(icol,ilev)+1,iswband) &
                      - wrh(icol,ilev)  * asm(krh(icol,ilev),  iswband)
  
            tau    (icol, ilev, iswband) = mass(icol, ilev) * ext1
            tau_w  (icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1
            tau_w_g(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1
            tau_w_f(icol, ilev, iswband) = mass(icol, ilev) * ext1 * ssa1 * asm1 * asm1
         enddo
      enddo
   enddo

end subroutine get_hygro_rad_props 

!==============================================================================
    
subroutine get_nonhygro_rad_props(ncol, mass, ext, ssa, asm, &
                                  tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   real(r8), intent(in) :: ext(:)
   real(r8), intent(in) :: ssa(:)
   real(r8), intent(in) :: asm(:)

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands) 

   ! Local variables
   integer  :: iswband
   real(r8) :: ext1, ssa1, asm1
   !-----------------------------------------------------------------------------
   
   do iswband = 1, nswbands
      ext1 = ext(iswband)
      ssa1 = ssa(iswband)
      asm1 = asm(iswband)
      tau    (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1
      tau_w  (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1
      tau_w_g(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1 * asm1
      tau_w_f(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext1 * ssa1 * asm1 * asm1
   enddo

end subroutine get_nonhygro_rad_props

!==============================================================================
    
subroutine get_volcanic_radius_rad_props(lchnk, ncol, mass, pbuf, r_ext, r_scat, r_ascat, r_mu, &
                                  tau, tau_w, tau_w_g, tau_w_f)

   use phys_buffer,       only: pbuf_get_fld_idx, pbuf_size_max

   ! Arguments
   integer,  intent(in) :: lchnk
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   real(r8), intent(in) :: r_ext(:,:)
   real(r8), intent(in) :: r_scat(:,:)
   real(r8), intent(in) :: r_ascat(:,:)
   real(r8), intent(in) :: r_mu(:) ! log(radius) domain of mass-specific optics

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands) 

   ! Local variables
   integer  :: iswband
   real(r8) :: g

   integer  :: idx                             ! index to radius in physics buffer
   real(r8), pointer :: geometric_radius(:,:)  ! geometric mean radius of volcanic aerosol
   real(r8) :: mu(pcols,pver)                  ! log(geometric mean radius of volcanic aerosol)
   integer  :: imu                             ! index into table values of log radius
   integer  :: kmu, nmu
   real(r8) :: wmu, mutrunc, r_mu_max, r_mu_min
 
   ! interpolated values from table
   real(r8) :: ext(nswbands)
   real(r8) :: scat(nswbands)
   real(r8) :: ascat(nswbands)

   integer :: i, k ! column level iterator
   !-----------------------------------------------------------------------------

   tau    =0._r8                 
   tau_w  =0._r8                 
   tau_w_g=0._r8                 
   tau_w_f=0._r8                  

   ! get microphysical properties for volcanic aerosols
   idx = pbuf_get_fld_idx('VOLC_RAD_GEOM')
   geometric_radius => pbuf(idx)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1 )

   ! interpolate in radius
   ! caution: clip the table with no warning when outside bounds
   nmu = size(r_mu)
   r_mu_max = r_mu(nmu)
   r_mu_min = r_mu(1)
   do i = 1, ncol
      do k = 1, pver
         if(geometric_radius(i,k) > 0.) then
            mu(i,k) = log(geometric_radius(i,k))
         else
            mu(i,k) = 0.
         endif
         mutrunc = max(min(mu(i,k),r_mu_max),r_mu_min)
         kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
         wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
         do iswband = 1, nswbands
            ext(iswband) =  &
               ((1._r8 - wmu) * r_ext(iswband, kmu  ) + &
               (wmu) * r_ext(iswband, kmu+1))
            scat(iswband) =  &
               ((1._r8 - wmu) * r_scat(iswband, kmu  ) + &
               (wmu) * r_scat(iswband, kmu+1))
            ascat(iswband) =  &
               ((1._r8 - wmu) * r_ascat(iswband, kmu  ) + &
               (wmu) * r_ascat(iswband, kmu+1))
            if (scat(iswband).gt.0.) then
               g = ascat(iswband)/scat(iswband)
            else
               g=0._r8
            endif
            tau    (i,k,iswband) = mass(i,k) * ext(iswband)  
            tau_w  (i,k,iswband) = mass(i,k) * scat(iswband)  
            tau_w_g(i,k,iswband) = mass(i,k) * ascat(iswband)  
            tau_w_f(i,k,iswband) = mass(i,k) * g * ascat(iswband)  
         end do
      enddo
   enddo

end subroutine get_volcanic_radius_rad_props

!==============================================================================
    
subroutine get_volcanic_rad_props(ncol, mass, ext, scat, ascat, &
                                  tau, tau_w, tau_w_g, tau_w_f)

   ! Arguments
   integer,  intent(in) :: ncol
   real(r8), intent(in) :: mass(pcols, pver)
   real(r8), intent(in) :: ext(:)
   real(r8), intent(in) :: scat(:)
   real(r8), intent(in) :: ascat(:)

   real(r8), intent(out) :: tau    (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w  (pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_g(pcols, pver, nswbands)
   real(r8), intent(out) :: tau_w_f(pcols, pver, nswbands) 

   ! Local variables
   integer  :: iswband
   real(r8) :: g
   !-----------------------------------------------------------------------------
   
   do iswband = 1, nswbands
      if (scat(iswband).gt.0.) then
         g = ascat(iswband)/scat(iswband)
      else
         g=0._r8
      endif
      tau    (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ext(iswband) 
      tau_w  (1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * scat(iswband) 
      tau_w_g(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * ascat(iswband) 
      tau_w_f(1:ncol,1:pver,iswband) = mass(1:ncol,1:pver) * g * ascat(iswband) 
   enddo

end subroutine get_volcanic_rad_props

!==============================================================================

subroutine aer_vis_diag_out(lchnk, ncol, nnite, idxnite, iaer, tau, diag_idx)

   ! output aerosol optical depth for the visible band

   integer,          intent(in) :: lchnk
   integer,          intent(in) :: ncol           ! number of columns
   integer,          intent(in) :: nnite          ! number of night columns
   integer,          intent(in) :: idxnite(:)     ! local column indices of night columns
   integer,          intent(in) :: iaer           ! aerosol index -- if 0 then tau is a total for all aerosols
   real(r8),         intent(in) :: tau(:,:)       ! aerosol optical depth for the visible band
   integer,          intent(in) :: diag_idx       ! identifies whether the aerosol optics
                                                  ! is for the climate calc or a diagnostic calc
 
   ! Local variables
   integer  :: i
   real(r8) :: tmp(pcols)
   !-----------------------------------------------------------------------------

   ! currently only implemented for climate calc
   if (diag_idx > 0) return

   ! compute total column aerosol optical depth
   tmp(1:ncol) = sum(tau(1:ncol,:), 2)
   ! use fillvalue to indicate night columns
   do i = 1, nnite
      tmp(idxnite(i)) = fillvalue
   end do

   if (iaer > 0) then
      call outfld(odv_names(iaer), tmp, pcols, lchnk)
   else
      call outfld('AEROD_v', tmp, pcols, lchnk)
   end if

end subroutine aer_vis_diag_out

!==============================================================================

end module aer_rad_props