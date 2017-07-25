
module rad_constituents

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Provide constituent distributions to the radiation routines.
! 
! By default the subroutine that returns constituent mixing ratios returns
! the distribution to be used for the interactive calculation.  That method
! also provides an optional argument for indicating that the requested
! distribution is passive and only used in the diagnostic radiative forcing
! calculation.
! 
! The logic to control which constituent distribution is interactive and
! which is passive is contained in this module.  By default, if a prognostic
! version of a constituent is found (by looking in the constituent array),
! then it is used for the interactive calculation, and a prescribed
! distribution will be used for the diagnostic calculation.
!
! Revision history:
! 2004-08-28  B. Eaton      Original version
! 2005-02-01  B. Eaton      Add functionality for non-constant CO2
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: masterproc
use error_messages, only: alloc_err   
use abortutils,     only: endrun
use cam_logfile,    only: iulog
use ppgrid,         only: pcols, pver
use physconst,      only: rga
use physics_types,  only: physics_state
use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_get_fld_name
use constituents,   only: cnst_name, cnst_get_ind
use radconstants,   only: gasnamelength, nradgas, gaslist, rad_gas_index, ot_length
use cam_history,    only: addfld, fieldname_len, phys_decomp, add_default, outfld

implicit none
private
save

! Public interfaces

public :: &
   rad_cnst_readnl,             &! read namelist values and parse
   rad_cnst_init,               &! find optics files and all constituents
   rad_cnst_get_clim_info,      &! return info about climate lists
   rad_cnst_get_clim_gas,       &! return pointer to mmr for gasses in climate list
   rad_cnst_get_clim_aer,       &! return pointer to mmr for aerosols in climate list
   rad_cnst_get_clim_aer_props, &! return physical properties for aerosols in the climate list
   rad_cnst_out,                &! output constituent diagnostics (mass per layer and column burden)
   rad_cnst_get_diag_list        ! return list of diagnostic calls to radiation

! Private module data

integer, parameter :: N_RAD_CNST = 100
integer, public, parameter :: N_DIAG = 10
integer, parameter :: cs1 = 256
character(len=cs1),public :: iceopticsfile, liqopticsfile
character(len=32),public :: icecldoptics,liqcldoptics
logical,public :: oldcldoptics = .false.

! Namelist variables
character(len=cs1), dimension(N_RAD_CNST) :: rad_climate = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_1  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_2  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_3  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_4  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_5  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_6  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_7  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_8  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_9  = ' '
character(len=cs1), dimension(N_RAD_CNST) :: rad_diag_10 = ' '

type :: rad_cnst_namelist_t
   integer :: ncnst
   character(len=  1), pointer :: source(:)  ! 'D' if in pbuf or 'P' if in state
   character(len= 64), pointer :: camname(:) ! name registered in pbuf or constituents
   character(len=cs1), pointer :: radname(:) ! radname is the name as identfied in radiation,
                                             ! must be one of (rgaslist if a gas) or
                                             ! (/fullpath/filename.nc if an aerosol)
   character(len=  1), pointer :: type(:)    ! 'A' if aerosol or 'G' if gas
end type rad_cnst_namelist_t

! Storage for parsed namelist variables
type(rad_cnst_namelist_t) :: clim_namelist         ! list of constituents interacting with climate
type(rad_cnst_namelist_t) :: diag_namelist(N_DIAG) ! constituents for diagnostic calculations
logical :: calldiagnostic(N_DIAG)     ! was nth diagnostic call to radiation specified?

type :: gas_t
   character(len=1)  :: source ='Z'  ! P is state, D is pbuf, Z is as near to zero as rad allows
   character(len=64) :: camname      ! name of constituent in physics state or buffer
   character(len=32) :: mass_name    ! name for mass per layer field in history output
   integer           :: idx          ! index from constituents or from pbuf
end type gas_t

type :: gaslist_t
   integer                :: ngas
   character(len=2)       :: list_id  ! set to "  " for climate list, or two character integer
                                      ! (include leading zero) to identify diagnostic list
   type(gas_t), pointer   :: gas(:)   ! dimension(ngas) where ngas = nradgas is from radconstants
end type gaslist_t

! Storage for gas identifiers
type(gaslist_t), target :: clim_gaslist          ! gasses interacting with climate
type(gaslist_t), target :: diag_gaslist(N_DIAG)  ! gasses used in diagnostic calculations


type :: aerosol_t
   character(len=1)   :: source         ! (numaersols) P is state, D is pbuf.
   character(len=64)  :: camname        ! name of constituent in physics state or buffer
   character(len=32)  :: mass_name      ! name for mass per layer field in history output
   integer            :: idx            ! index of constituent in physics state or buffer
   integer            :: physprop_id    ! ID used to access physical properties from phys_prop module
end type aerosol_t

type :: aerlist_t
  integer                  :: numaerosols  ! number of aerosols
  character(len=2)         :: list_id      ! set to "  " for climate list, or two character integer
                                           ! (include leading zero) to identify diagnostic list
  type(aerosol_t), pointer :: aer(:)       ! dimension(numaerosols)
end type aerlist_t

! Storage for aerosol identifiers and associated physical properties
type(aerlist_t), target :: clim_aerosollist         ! list of aerosols interacting with climate
type(aerlist_t), target :: diag_aerosollist(N_DIAG) ! list of aerosols used in diagnostic calcs

! mmr values for gasses required by radiation but for which no source data is specified
real(r8), allocatable, target :: zerommr(:,:) 

! produce a bunch of output about what is being associated by this module
logical :: debugthiscode = .true.

!==============================================================================
contains
!==============================================================================

subroutine rad_cnst_readnl(nlfile)

   ! Read rad_cnst_nl namelist group.  Parse input.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=*), parameter :: subname = 'rad_cnst_readnl'

   namelist /rad_cnst_nl/ rad_climate,  &
                          rad_diag_1, &
                          rad_diag_2, &
                          rad_diag_3, &
                          rad_diag_4, &
                          rad_diag_5, &
                          rad_diag_6, &
                          rad_diag_7, &
                          rad_diag_8, &
                          rad_diag_9, &
                          rad_diag_10,&
                          iceopticsfile, &
                          liqopticsfile, &
                          icecldoptics, &
                          liqcldoptics, &
                          oldcldoptics

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rad_cnst_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rad_cnst_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (rad_climate, len(rad_climate(1))*N_RAD_CNST,   mpichar, 0, mpicom)
   call mpibcast (rad_diag_1,  len(rad_diag_1(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_2,  len(rad_diag_2(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_3,  len(rad_diag_3(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_4,  len(rad_diag_4(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_5,  len(rad_diag_5(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_6,  len(rad_diag_6(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_7,  len(rad_diag_7(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_8,  len(rad_diag_8(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_9,  len(rad_diag_9(1))*N_RAD_CNST,    mpichar, 0, mpicom)
   call mpibcast (rad_diag_10, len(rad_diag_10(1))*N_RAD_CNST,   mpichar, 0, mpicom)
   call mpibcast (iceopticsfile, len(iceopticsfile),   mpichar, 0, mpicom)
   call mpibcast (liqopticsfile, len(liqopticsfile),   mpichar, 0, mpicom)
   call mpibcast (liqcldoptics, len(liqcldoptics),   mpichar, 0, mpicom)
   call mpibcast (icecldoptics, len(icecldoptics),   mpichar, 0, mpicom)
   call mpibcast (oldcldoptics, 1,   mpilog , 0, mpicom)
#endif

   ! Parse the namelist input.
   call parse_rad_specifier(rad_climate, clim_namelist)
   do i = 1,N_DIAG
      select case (i)
      case (1)
         call parse_rad_specifier(rad_diag_1, diag_namelist(i))
      case (2)
         call parse_rad_specifier(rad_diag_2, diag_namelist(i))
      case (3)
         call parse_rad_specifier(rad_diag_3, diag_namelist(i))
      case (4)
         call parse_rad_specifier(rad_diag_4, diag_namelist(i))
      case (5)
         call parse_rad_specifier(rad_diag_5, diag_namelist(i))
      case (6)
         call parse_rad_specifier(rad_diag_6, diag_namelist(i))
      case (7)
         call parse_rad_specifier(rad_diag_7, diag_namelist(i))
      case (8)
         call parse_rad_specifier(rad_diag_8, diag_namelist(i))
      case (9)
         call parse_rad_specifier(rad_diag_9, diag_namelist(i))
      case (10)
         call parse_rad_specifier(rad_diag_10, diag_namelist(i))
      end select
   enddo

   ! were there any constituents specified for the nth diagnostic call?
   ! if so, radiation will make a call with those consituents
   calldiagnostic(:) = (diag_namelist(:)%ncnst > 0)

end subroutine rad_cnst_readnl

!================================================================================================

subroutine rad_cnst_init(pbuf, phys_state)

   use phys_prop, only: physprop_accum_unique_files, physprop_init

   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   type(physics_state), pointer    :: phys_state(:) 

   integer :: num_aerosols
   integer :: i
   character(len=2) :: suffix
   logical, parameter :: stricttest = .true.
   !-----------------------------------------------------------------------------

   ! The list_id component can be used to identify whether a given list is the
   ! active (climate) constituents or one of the lists for diagnostic calculations
   clim_gaslist%list_id = "  "
   clim_aerosollist%list_id = "  "

   ! Create a list of the unique set of filenames containing property data
   call physprop_accum_unique_files(clim_namelist%radname, clim_namelist%type, num_aerosols)
   clim_aerosollist%numaerosols = num_aerosols

   do i = 1,N_DIAG
      if ( diag_namelist(i)%ncnst > 0 ) then
         call physprop_accum_unique_files(diag_namelist(i)%radname, diag_namelist(i)%type, num_aerosols)
         write ( suffix, fmt = '(i2.2)' ) i
         diag_aerosollist(i)%numaerosols = num_aerosols
         diag_aerosollist(i)%list_id = suffix
         diag_gaslist(i)%list_id = suffix
      else
         diag_aerosollist(i)%numaerosols = 0
      endif
   enddo

   ! Allocate storage for the physical properties of each aerosol; read properties from
   ! the data files.
   call physprop_init()

   ! Start checking that specified radiative constituents are present
   if (masterproc) write(iulog,*) 'rad_cnst_init: checking for radiative constituents'

   ! Initialize the gas and aerosol lists for the constituents affecting the climate.
   call init_lists(clim_namelist, clim_gaslist, clim_aerosollist, pbuf, phys_state)

   ! Check that all gases supported by the radiative transfer code have been specified.
   if(stricttest) then
      do i = 1, nradgas
         if (clim_gaslist%gas(i)%source .eq. 'Z' ) then
            call endrun("list of radiative gasses must include all radiation gasses for the climate specication")
         endif
      enddo
   endif

   if (masterproc .and. debugthiscode) call debug_nl(clim_gaslist, clim_aerosollist)

   ! Initialize history output of climate diagnostic quantities
   call rad_gas_diag_init(clim_gaslist)
   call rad_aer_diag_init(clim_aerosollist)

   ! Initialize the gas and aerosol lists for the constituents used in diagnostic calculations
   do i = 1,N_DIAG
      if (diag_namelist(i)%ncnst > 0) then

         call init_lists(diag_namelist(i), diag_gaslist(i), diag_aerosollist(i), pbuf, phys_state)

         if (masterproc .and. debugthiscode) call debug_nl(diag_gaslist(i), diag_aerosollist(i))

      endif
   enddo

   ! memory to point to if no mass is specified
   allocate(zerommr(pcols,pver))
   zerommr = 0._r8

end subroutine rad_cnst_init

!================================================================================================

subroutine rad_cnst_get_clim_gas(gasname, state, pbuf, mmr, diagnosticindex)

! Return pointer to mass mixing ratio for the the instance of the gas
! specified to affect the climate simulation.

   ! Arguments
   character(len=*),            intent(in) :: gasname
   type(physics_state), target, intent(in) :: state
   type(pbuf_fld),              intent(in) :: pbuf(pbuf_size_max)
   real(r8),                    pointer    :: mmr(:,:)
   integer, optional,           intent(in) :: diagnosticindex

   ! Local variables
   integer :: lchnk
   integer :: igas
   integer :: idx
   character(len=1) :: source
   type(gaslist_t),pointer :: gaslist
   !-----------------------------------------------------------------------------

   gaslist => clim_gaslist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           gaslist => diag_gaslist(diagnosticindex)
       endif
   endif
      
   lchnk = state%lchnk

   ! Get index of gas in internal arrays.  rad_gas_index will abort if the 
   ! specified gasname is not recognized by the radiative transfer code.
   igas = rad_gas_index(trim(gasname))
 
   ! Get data source from the climate gaslist
   source = gaslist%gas(igas)%source
   idx    = gaslist%gas(igas)%idx
   select case( source )
   case ('P')
      mmr => state%q(:,:,idx)
   case ('D')
      mmr => pbuf(idx)%fld_ptr(1,:,:,lchnk,1)
   case ('Z')
      mmr => zerommr(:,:)
   end select

end subroutine rad_cnst_get_clim_gas

!================================================================================================

subroutine rad_cnst_get_clim_aer(list_idx, state, pbuf, mmr, diagnosticindex)

! Return pointer to mass mixing ratio for the the instance of the aerosol
! specified to affect the climate simulation.

! This method is used to access the aerosols in the climate list sequentially
! without specifying aerosol names.

   ! Arguments
   integer,                     intent(in) :: list_idx
   type(physics_state), target, intent(in) :: state
   type(pbuf_fld),              intent(in) :: pbuf(pbuf_size_max)
   real(r8),                    pointer    :: mmr(:,:)
   integer, optional,           intent(in) :: diagnosticindex

   ! Local variables
   integer :: lchnk
   integer :: idx
   character(len=1) :: source
   type(aerlist_t),pointer :: aerlist
   !-----------------------------------------------------------------------------

   aerlist => clim_aerosollist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           aerlist => diag_aerosollist(diagnosticindex)
       endif
   endif

   lchnk = state%lchnk

   ! Check for valid input index
   if (list_idx < 1  .or.  list_idx > aerlist%numaerosols) then
      write(iulog,*) 'list_idx= ', list_idx, '  numaerosols= ', aerlist%numaerosols
      call endrun('rad_cnst_get_clim_aer: aerosol list index out of range')
   end if

   ! Get data source from the climate aerosollist
   source = aerlist%aer(list_idx)%source
   idx    = aerlist%aer(list_idx)%idx
   select case( source )
   case ('P')
      mmr => state%q(:,:,idx)
   case ('D')
      mmr => pbuf(idx)%fld_ptr(1,:,:,lchnk,1)
   end select

end subroutine rad_cnst_get_clim_aer

!================================================================================================

subroutine rad_cnst_get_clim_info( naero, aernames, aersources, aerindices, &
                                   ngas,  gasnames, gassources, gasindices, &
                                   use_data_o3, diagnosticindex )

! Return info about aerosol climate list

   ! Arguments
   integer,           optional, intent(out) :: naero
   integer,           optional, intent(out) :: ngas
   character(len=64), optional, intent(out) :: aernames(:)
   character(len=64), optional, intent(out) :: gasnames(:)
   character(len=1),  optional, intent(out) :: aersources(:)
   character(len=1),  optional, intent(out) :: gassources(:)
   integer,           optional, intent(out) :: aerindices(:)
   integer,           optional, intent(out) :: gasindices(:)
   logical,           optional, intent(out) :: use_data_o3
   integer,           optional, intent(in)  :: diagnosticindex

   ! Local variables
   type(gaslist_t), pointer :: gaslist ! local pointer to gas list of interest
   type(aerlist_t), pointer :: aerlist ! local pointer to aerosol list of interest

   integer          :: i
   integer          :: arrlen  ! length of assumed shape array
   integer          :: gaslen  ! length of assumed shape array
   integer          :: igas    ! index of a gas in the gas list
   character(len=1) :: source  ! P is state, D is pbuf, Z is as near to zero as rad allows

   character(len=*), parameter :: subname = 'rad_cnst_get_clim_info'
   !-----------------------------------------------------------------------------

   gaslist => clim_gaslist
   aerlist => clim_aerosollist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           gaslist => diag_gaslist(diagnosticindex)
           aerlist => diag_aerosollist(diagnosticindex)
       endif
   endif

   ! number of aerosols in list
   if (present(naero)) then
      naero = aerlist%numaerosols
   endif
   ! number of gases in list
   if (present(ngas)) then
      ngas = gaslist%ngas
   endif

   ! names of aerosols in list
   if (present(aernames)) then

      ! check that output array is long enough
      arrlen = size(aernames)
      if (arrlen < aerlist%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', aerlist%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aernames too short')
      end if

      do i = 1, aerlist%numaerosols
         aernames(i) = aerlist%aer(i)%camname
      end do

   end if

   ! sources of aerosols in list
   if (present(aersources)) then

      ! check that output array is long enough
      arrlen = size(aersources)
      if (arrlen < aerlist%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', aerlist%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aersources too short')
      end if

      do i = 1, aerlist%numaerosols
         aersources(i) = aerlist%aer(i)%source
      end do

   end if

   ! indices of aerosols in list
   if (present(aerindices)) then

      ! check that output array is long enough
      arrlen = size(aerindices)
      if (arrlen < aerlist%numaerosols) then
         write(iulog,*) subname//': ERROR: naero=', aerlist%numaerosols, '  arrlen=', arrlen
         call endrun(subname//': ERROR: aerindices too short')
      end if

      do i = 1, aerlist%numaerosols
         aerindices(i) = aerlist%aer(i)%idx
      end do

   end if

   ! names of gas in list
   if (present(gasnames)) then

      ! check that output array is long enough
      gaslen = size(gasnames)
      if (gaslen < gaslist%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', gaslist%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasnames too short')
      end if

      do i = 1, gaslist%ngas
         gasnames(i) = gaslist%gas(i)%camname
      end do

   end if

   ! sources of gas in list
   if (present(gassources)) then

      ! check that output array is long enough
      gaslen = size(gassources)
      if (gaslen < gaslist%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', gaslist%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gassources too short')
      end if

      do i = 1, gaslist%ngas
         gassources(i) = gaslist%gas(i)%source
      end do

   end if

   ! indices of gas in list
   if (present(gasindices)) then

      ! check that output array is long enough
      gaslen = size(gasindices)
      if (gaslen < gaslist%ngas) then
         write(iulog,*) subname//': ERROR: ngas=', gaslist%ngas, '  gaslen=', gaslen
         call endrun(subname//': ERROR: gasindices too short')
      end if

      do i = 1, gaslist%ngas
         gasindices(i) = gaslist%gas(i)%idx
      end do

   end if

   ! Does the climate calculation use data ozone?
   if (present(use_data_o3)) then

      ! get index of O3 in gas list
      igas = rad_gas_index('O3')
 
      ! Get data source from the climate gaslist
      source = gaslist%gas(igas)%source

      use_data_o3 = .false.
      if (source == 'D') use_data_o3 = .true.
   endif

end subroutine rad_cnst_get_clim_info

!================================================================================================

subroutine rad_cnst_get_diag_list(diagnosticcalllist)

! Return info about aerosol climate list

   ! Arguments
   logical, intent(out) :: diagnosticcalllist(N_DIAG)
   !-----------------------------------------------------------------------------

   diagnosticcalllist(:) = calldiagnostic(:)

end subroutine rad_cnst_get_diag_list

!================================================================================================

subroutine rad_cnst_get_clim_aer_props( &
   list_idx, diagnosticindex, opticstype, &
   sw_hygro_ext, sw_hygro_ssa, sw_hygro_asm, lw_hygro_ext, &
   sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm, &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_ext, &
   refindex_real_aer_sw, refindex_im_aer_sw, refindex_real_aer_lw, refindex_im_aer_lw, &
   refindex_real_water_sw, refindex_im_water_sw, refindex_real_water_lw, refindex_im_water_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, num_to_mass_aer)

   ! Return requested properties for specified climate aerosol.

   use phys_prop, only: physprop_get

   ! Arguments
   integer,                     intent(in)  :: list_idx
   integer,           optional, intent(in)  :: diagnosticindex ! index to rad diagnostic call
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:) 
   real(r8),          optional, pointer     :: lw_hygro_ext(:,:)         
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_ext(:)         
   real(r8),          optional, pointer     :: refindex_real_aer_sw(:)
   real(r8),          optional, pointer     :: refindex_im_aer_sw(:)
   real(r8),          optional, pointer     :: refindex_real_aer_lw(:)
   real(r8),          optional, pointer     :: refindex_im_aer_lw(:)
   real(r8),          optional, pointer     :: refindex_real_water_sw(:)
   real(r8),          optional, pointer     :: refindex_im_water_sw(:)
   real(r8),          optional, pointer     :: refindex_real_water_lw(:)
   real(r8),          optional, pointer     :: refindex_im_water_lw(:)
   character(len=20), optional, intent(out) :: aername           
   real(r8),          optional, intent(out) :: density_aer
   real(r8),          optional, intent(out) :: hygro_aer
   real(r8),          optional, intent(out) :: dryrad_aer        
   real(r8),          optional, intent(out) :: dispersion_aer    
   real(r8),          optional, intent(out) :: num_to_mass_aer   

   real(r8),          optional, pointer     :: r_sw_ext(:,:)         
   real(r8),          optional, pointer     :: r_sw_scat(:,:)         
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)         
   real(r8),          optional, pointer     :: r_lw_abs(:,:)         
   real(r8),          optional, pointer     :: mu(:)         

   ! Local variables
   integer :: id
   character(len=*), parameter :: subname = 'rad_cnst_get_clim_aer_props'
   type(aerlist_t), pointer :: aerlist
   !------------------------------------------------------------------------------------

   aerlist => clim_aerosollist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           aerlist => diag_aerosollist(diagnosticindex)
       endif
   endif

   if (list_idx <= 0 .or. list_idx > aerlist%numaerosols) then
      if (.not. present(diagnosticindex)) then
         write(iulog,*) subname//': climate aerosol list index out of range: ', list_idx
      else
         write(iulog,*) subname//': aerosol list index out of range: ', list_idx ,' diag call: ',diagnosticindex
      endif
      call endrun(subname//': list_idx out of range')
   end if

   id = aerlist%aer(list_idx)%physprop_id

   if (present(opticstype))        call physprop_get(id, opticstype=opticstype)
   if (present(sw_hygro_ext))      call physprop_get(id, sw_hygro_ext=sw_hygro_ext)
   if (present(sw_hygro_ssa))      call physprop_get(id, sw_hygro_ssa=sw_hygro_ssa)
   if (present(sw_hygro_asm))      call physprop_get(id, sw_hygro_asm=sw_hygro_asm)
   if (present(lw_hygro_ext))      call physprop_get(id, lw_hygro_abs=lw_hygro_ext)
   if (present(sw_nonhygro_ext))   call physprop_get(id, sw_nonhygro_ext=sw_nonhygro_ext)
   if (present(sw_nonhygro_ssa))   call physprop_get(id, sw_nonhygro_ssa=sw_nonhygro_ssa)
   if (present(sw_nonhygro_asm))   call physprop_get(id, sw_nonhygro_asm=sw_nonhygro_asm)
   if (present(sw_nonhygro_scat))  call physprop_get(id, sw_nonhygro_scat=sw_nonhygro_scat)
   if (present(sw_nonhygro_ascat)) call physprop_get(id, sw_nonhygro_ascat=sw_nonhygro_ascat)
   if (present(lw_ext))            call physprop_get(id, lw_abs=lw_ext)
   if (present(refindex_real_aer_sw))   call physprop_get(id, refindex_real_aer_sw=refindex_real_aer_sw)
   if (present(refindex_im_aer_sw))     call physprop_get(id, refindex_im_aer_sw=refindex_im_aer_sw)
   if (present(refindex_real_aer_lw))   call physprop_get(id, refindex_real_aer_lw=refindex_real_aer_lw)
   if (present(refindex_im_aer_lw))     call physprop_get(id, refindex_im_aer_lw=refindex_im_aer_lw)
   if (present(refindex_real_water_sw)) call physprop_get(id, refindex_real_water_sw=refindex_real_water_sw)
   if (present(refindex_im_water_sw))   call physprop_get(id, refindex_im_water_sw=refindex_im_water_sw)
   if (present(refindex_real_water_lw)) call physprop_get(id, refindex_real_water_lw=refindex_real_water_lw)
   if (present(refindex_im_water_lw))   call physprop_get(id, refindex_im_water_lw=refindex_im_water_lw)
   if (present(aername))           call physprop_get(id, aername=aername)
   if (present(density_aer))       call physprop_get(id, density_aer=density_aer)
   if (present(hygro_aer))         call physprop_get(id, hygro_aer=hygro_aer)
   if (present(dryrad_aer))        call physprop_get(id, dryrad_aer=dryrad_aer)
   if (present(dispersion_aer))    call physprop_get(id, dispersion_aer=dispersion_aer)
   if (present(num_to_mass_aer))   call physprop_get(id, num_to_mass_aer=num_to_mass_aer)

   if (present(r_lw_abs))          call physprop_get(id, r_lw_abs=r_lw_abs)
   if (present(r_sw_ext))          call physprop_get(id, r_sw_ext=r_sw_ext)
   if (present(r_sw_scat))         call physprop_get(id, r_sw_scat=r_sw_scat)
   if (present(r_sw_ascat))        call physprop_get(id, r_sw_ascat=r_sw_ascat)
   if (present(mu))                call physprop_get(id, mu=mu)

end subroutine rad_cnst_get_clim_aer_props
!================================================================================================

subroutine rad_cnst_out(state, pbuf, diagnosticindex)

! Output the mass per layer, and total column burdens for gas and aerosol
! constituents in either the climate or diagnostic lists

   ! Arguments
   type(physics_state), target, intent(in) :: state
   type(pbuf_fld),              intent(in) :: pbuf(pbuf_size_max)
   integer, optional,           intent(in) :: diagnosticindex

   ! Local variables
   integer :: i, naer, ngas, lchnk, ncol
   integer :: idx
   character(len=1)  :: source
   character(len=32) :: name, cbname
   real(r8)          :: mass(pcols,pver)
   real(r8)          :: cb(pcols)
   real(r8), pointer :: mmr(:,:)
   type(aerlist_t), pointer :: aerlist
   type(gaslist_t),pointer :: gaslist
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointer with requested aerosol list
   aerlist => clim_aerosollist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           aerlist => diag_aerosollist(diagnosticindex)
       endif
   endif

   naer = aerlist%numaerosols
   do i = 1, naer

      source = aerlist%aer(i)%source
      idx    = aerlist%aer(i)%idx
      name   = aerlist%aer(i)%mass_name
      ! construct name for column burden field by replacing the 'm_' prefix by 'cb_'
      cbname = 'cb_' // name(3:len_trim(name))

      select case( source )
      case ('P')
         mmr => state%q(:,:,idx)
      case ('D')
         mmr => pbuf(idx)%fld_ptr(1,:,:,lchnk,1)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

   ! Associate pointer with requested gas list
   gaslist => clim_gaslist
   if (present(diagnosticindex)) then
       if (diagnosticindex<N_DIAG) then
           gaslist => diag_gaslist(diagnosticindex)
       endif
   endif
   ngas = gaslist%ngas
   do i = 1, ngas

      source = gaslist%gas(i)%source
      idx    = gaslist%gas(i)%idx
      name   = gaslist%gas(i)%mass_name
      cbname = 'cb_' // name(3:len_trim(name))
      select case( source )
      case ('P')
         mmr => state%q(:,:,idx)
      case ('D')
         mmr => pbuf(idx)%fld_ptr(1,:,:,lchnk,1)
      end select

      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga
      call outfld(trim(name), mass, pcols, lchnk)

      cb(:ncol) = sum(mass(:ncol,:),2)
      call outfld(trim(cbname), cb, pcols, lchnk)

   end do

end subroutine rad_cnst_out

!================================================================================================
! Private methods
!================================================================================================

subroutine init_lists(namelist, gaslist, aerlist, pbuf, phys_state)

   ! Initialize the gas and aerosol lists with the constituents specified in the
   ! namelist.

   use ppgrid,        only: pcols, pver, begchunk, endchunk
   use phys_buffer,   only: pbuf_get_fld_idx, pbuf_size_max, pbuf_fld
   use physics_types, only: physics_state
   use phys_prop,     only: physprop_get_id

   type(rad_cnst_namelist_t), intent(in) :: namelist            ! namelist input for climate or diagnostics
   type(pbuf_fld),            intent(in) :: pbuf(pbuf_size_max)
   type(physics_state),       pointer    :: phys_state(:) 

   type(gaslist_t),        intent(inout) :: gaslist
   type(aerlist_t),        intent(inout) :: aerlist


   ! Local variables
   integer :: constindx, idx
   integer :: igas, ifileindex, aeridx
   integer :: istat
   integer, parameter :: pbuffailure = -1 ! if pbuf fails, return this value
   !-----------------------------------------------------------------------------

   ! nradgas is set by the radiative transfer code
   gaslist%ngas = nradgas

   ! index into list of aerosols
   aeridx = 1

   ! allocate storage for the aerosol and gas lists
   allocate(aerlist%aer(aerlist%numaerosols), stat=istat)
   call alloc_err(istat, 'init_lists', 'aerlist%aer', aerlist%numaerosols)
   allocate(gaslist%gas(gaslist%ngas), stat=istat)
   call alloc_err(istat, 'rad_constituents', 'gaslist%gas', gaslist%ngas)

   ! Loop over the constituents specified in the namelist
   do constindx = 1, namelist%ncnst

      if (masterproc .and. debugthiscode) &
         write(iulog,*) "rad namelist spec: "// trim(namelist%source(constindx)) &
         //"_"//trim(namelist%camname(constindx))//":"//trim(namelist%radname(constindx))

      ! Check that the source specifier is legal.
      if (namelist%source(constindx) /= 'D' .and. namelist%source(constindx) /= 'P') then
         call endrun("init_lists: constituent source must either be D or P:"//&
                     " illegal specifier in namelist input: "//namelist%source(constindx))
      end if

      ! Add constituent to either the aerosol or the gas list.
      if (namelist%type(constindx) == 'A') then 

         ! Add to aerosol list

         aerlist%aer(aeridx)%source  = namelist%source(constindx)
         aerlist%aer(aeridx)%camname = namelist%camname(constindx)

         ! get the physprop_id from the phys_prop module
         aerlist%aer(aeridx)%physprop_id = physprop_get_id(namelist%radname(constindx))

         ! locate in the pbuf or state
         if (namelist%source(constindx) == 'D') then

            idx = pbuf_get_fld_idx(namelist%camname(constindx), pbuffailure)
            if (idx == pbuffailure) then
               call endrun('init_lists: data constituent not found in pbuf: '//namelist%camname(constindx))
            end if
            aerlist%aer(aeridx)%idx = idx

         else if (namelist%source(constindx) == 'P') then

            call cnst_get_ind(namelist%camname(constindx), idx, abort=.false.)
            if (idx < 0) then
               call endrun('init_lists: prognostic constituent not found in state: '//namelist%camname(constindx))
            end if
            aerlist%aer(aeridx)%idx = idx

         endif

         aeridx = aeridx + 1

      else 

         ! Add to gas list

         ! The radiative transfer code requires the input of a specific set of gases
         ! which is hardwired into the code.  The CAM interface to the RT code uses
         ! the names in the radconstants module to refer to these gases.  The user
         ! interface (namelist) also uses these names to identify the gases treated
         ! by the RT code.  We use the index order set in radconstants for convenience
         ! only.

         ! First check that the gas name specified by the user is allowed.
         ! rad_gas_index will abort on illegal names.
         igas = rad_gas_index(namelist%radname(constindx))

         ! Set values in the igas index
         gaslist%gas(igas)%source  = namelist%source(constindx)
         gaslist%gas(igas)%camname = namelist%camname(constindx)

         ! locate in pbuf or state
         if (namelist%source(constindx) == 'D') then
            idx = pbuf_get_fld_idx(namelist%camname(constindx), pbuffailure)
            if (idx == pbuffailure) then 
               call endrun('init_lists: data constituent not found in pbuf: '//namelist%camname(constindx))
            end if
            gaslist%gas(igas)%idx = idx

         else if (namelist%source(constindx) == 'P') then
            call cnst_get_ind(namelist%camname(constindx), idx, abort=.false.)
            if (idx < 0) then
               call endrun('init_lists: prognostic constituent not found in state: '//namelist%camname(constindx))
            end if
            gaslist%gas(igas)%idx = idx

         endif

      endif ! (gas vs aerosol)

   enddo    ! loop over constituents

end subroutine init_lists

!================================================================================================

subroutine rad_gas_diag_init(glist)

! Add diagnostic fields to the master fieldlist.

   type(gaslist_t), intent(inout) :: glist

   integer :: i, ngas
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   ngas = glist%ngas
   if (ngas == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = glist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, ngas

      ! construct names for mass per layer diagnostics
      name = 'm_' // trim(glist%gas(i)%camname) // trim(suffix)
      glist%gas(i)%mass_name = name
      long_name = trim(glist%gas(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), 'kg/m^2', pver, 'A', trim(long_name), phys_decomp)

      ! construct names for column burden diagnostics
      name = 'cb_' // trim(glist%gas(i)%camname) // trim(suffix)
      long_name = trim(glist%gas(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), 'kg/m^2', 1, 'A', trim(long_name), phys_decomp)

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_gas_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_gas_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_gas_diag_init

!================================================================================================

subroutine rad_aer_diag_init(alist)

! Add diagnostic fields to the master fieldlist.

   type(aerlist_t), intent(inout) :: alist

   integer :: i, naer
   character(len=64) :: name
   character(len=2)  :: list_id
   character(len=4)  :: suffix
   character(len=128):: long_name
   character(len=32) :: long_name_description
   !-----------------------------------------------------------------------------

   naer = alist%numaerosols
   if (naer == 0) return

   ! Determine whether this is a climate or diagnostic list.
   list_id = alist%list_id
   if (len_trim(list_id) == 0) then
      suffix = '_c'
      long_name_description = ' used in climate calculation'
   else
      suffix = '_d' // list_id
      long_name_description = ' used in diagnostic calculation'
   end if

   do i = 1, naer

      ! construct names for mass per layer diagnostic fields
      name = 'm_' // trim(alist%aer(i)%camname) // trim(suffix)
      alist%aer(i)%mass_name = name
      long_name = trim(alist%aer(i)%camname)//' mass per layer'//long_name_description
      call addfld(trim(name), 'kg/m^2', pver, 'A', trim(long_name), phys_decomp)

      ! construct names for column burden diagnostic fields
      name = 'cb_' // trim(alist%aer(i)%camname) // trim(suffix)
      long_name = trim(alist%aer(i)%camname)//' column burden'//long_name_description
      call addfld(trim(name), 'kg/m^2', 1, 'A', trim(long_name), phys_decomp)

      ! error check for name length
      if (len_trim(name) > fieldname_len) then
         write(iulog,*) 'rad_aer_diag_init: '//trim(name)//' longer than ', fieldname_len, ' characters'
         call endrun('rad_aer_diag_init: name too long: '//trim(name))
      end if

   end do

end subroutine rad_aer_diag_init


!================================================================================================

subroutine parse_rad_specifier(specifier, namelist_data)

!-----------------------------------------------------------------------------
! Private method for parsing the radiation namelist specifiers.  The specifiers
! are of the form 'source_camname:radname' where:
! source  -- either 'D' for data or 'P' for prognostic
! camname -- the name of a constituent that must be found in the constituent
!            component of the state when source=P or in the physics buffer
!            when source=D
! radname -- For gases this is a name that identifies the constituent to the
!            radiative transfer codes.  These names are contained in the
!            radconstants module.  For aerosols this is a filename, which is
!            identified by a ".nc" suffix.  The file contains optical and 
!            other physical properties of the aerosol.
!
! This code also identifies whether the constituent is a gas or an aerosol
! and adds that info to a structure that stores the parsed data.
!-----------------------------------------------------------------------------

    character(len=*), dimension(:), intent(in) :: specifier
    type(rad_cnst_namelist_t),   intent(inout) :: namelist_data

    ! Local variables
    character(len=1),   dimension(N_RAD_CNST) :: source
    character(len=64),  dimension(N_RAD_CNST) :: camname
    character(len=cs1), dimension(N_RAD_CNST) :: radname
    character(len=1),   dimension(N_RAD_CNST) :: type
    integer :: number, i,j,k, astat
    !-------------------------------------------------------------------------
  
    number = 0

    parse_loop: do i = 1, N_RAD_CNST
      if ( len_trim(specifier(i)) == 0 ) then 
         exit parse_loop
      endif

      ! Locate the '-' separating source from camname.  This is the first
      ! occurance of '_' and allows for the possibility that camname or radname
      ! may also contain underscores.
      j = scan( specifier(i),'_')
      source(i) = trim(adjustl( specifier(i)(:j-1) ))

      ! locate the ':' separating camname from radname
      k = scan( specifier(i),':')
 
      camname(i) = trim(adjustl( specifier(i)(j+1:k-1) ))
      radname(i) = trim(adjustl( specifier(i)(k+1:) ))

      ! determine the type of constituent
      if (index(radname(i),".nc") .gt. 0) then
         type(i) = 'A'
      else
         type(i) = 'G'
      end if

      number = number+1    
    end do parse_loop

    namelist_data%ncnst = number

    if (number == 0) return

    allocate(namelist_data%source (number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%source')
    allocate(namelist_data%camname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%camname')
    allocate(namelist_data%radname(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%radname')
    allocate(namelist_data%type(number), stat=astat)
    if( astat/= 0 ) call endrun('parse_rad_specifier: not able to allocate namelist_data%type')

    namelist_data%source(:namelist_data%ncnst)  = source (:namelist_data%ncnst)
    namelist_data%camname(:namelist_data%ncnst) = camname(:namelist_data%ncnst)
    namelist_data%radname(:namelist_data%ncnst) = radname(:namelist_data%ncnst)
    namelist_data%type(:namelist_data%ncnst)    = type(:namelist_data%ncnst)

end subroutine parse_rad_specifier

!================================================================================================

subroutine debug_nl(dgaslist,daerosollist)

   use phys_prop,    only: physprop_get

   type(aerlist_t), intent(in) :: daerosollist
   type(gaslist_t), intent(in) :: dgaslist

   integer :: iconst, id
   character(len=256) :: filename

   write(iulog,*) '   gaslist'
   do iconst = 1,nradgas
      if (dgaslist%gas(iconst)%source .eq. 'D') then
         write(iulog,*) dgaslist%gas(iconst)%source//'_'//gaslist(iconst)//' has pbuf name:'//&
                        pbuf_get_fld_name(dgaslist%gas(iconst)%idx)
      else if (dgaslist%gas(iconst)%source .eq. 'P') then
         write(iulog,*) dgaslist%gas(iconst)%source//'_'//gaslist(iconst)//' has constituents name:'//&
                        cnst_name(dgaslist%gas(iconst)%idx)
      endif
   enddo

   write(iulog,*) '   aerosollist'
   do iconst = 1,daerosollist%numaerosols
      id = daerosollist%aer(iconst)%physprop_id
      call physprop_get(id, sourcefile=filename)
      write(iulog,*) trim(daerosollist%aer(iconst)%source)//'_'//trim(daerosollist%aer(iconst)%camname)//&
                     ' optics and phys props in :'//trim(filename)
   enddo

end subroutine debug_nl

!================================================================================================

end module rad_constituents
