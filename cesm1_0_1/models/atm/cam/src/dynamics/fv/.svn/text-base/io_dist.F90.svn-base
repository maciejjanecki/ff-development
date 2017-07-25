module io_dist

   use shr_kind_mod,       only : r8 => shr_kind_r8, r4 => shr_kind_r4
   use decompmodule,       only : decomptype
   use dynamics_vars,      only : T_FVDYCORE_GRID
   use dyn_internal_state, only : get_dyn_state_grid
   use abortutils,         only : endrun
   use cam_logfile,        only : iulog

#if ( defined SPMD )
   use mod_comm,           only: mp_sendirr, mp_recvirr, mp_sendirr_r4, mp_recvirr_r4, &
                                 mp_sendirr_i4, mp_recvirr_i4
   use parutilitiesmodule, only: parpatterntype
#endif

   implicit none
   private
   save

   public :: &
      fv_scatter_r4,&
      fv_gather_r4, &
      fv_read_r4,   &
      fv_write_r4,  &
      fv_scatter_r8,&
      fv_gather_r8, &
      fv_read_r8,   &
      fv_write_r8,  &
      fv_scatter_i4,&
      fv_gather_i4, &
      fv_read_i4

!====================================================================================
CONTAINS
!====================================================================================

   subroutine fv_scatter_r8(pattern_type, bufres, lenarr, arr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      real(r8) bufres(*)                                 ! Array to scatter
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      real(r8), intent(out)              :: arr(*)       ! Array to be gathered
#else
      real(r8), intent(out)              :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r8' )
!
!  Should check if this is a "scatter" pattern
!
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                      modc=grid%modc_scatter )
      CALL mp_recvirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                      modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy
#else 
      arr(1:lenarr) = bufres(1:lenarr)
#endif
      return
   end subroutine fv_scatter_r8

!====================================================================================

   subroutine fv_gather_r8( pattern_type, arr, lenarr, bufres)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r8) arr(*)            ! Array to be gathered
#else
      real(r8) arr(lenarr)       ! Array (SMP-only)
#endif
      real(r8), intent(out)              :: bufres(*)    ! Gathered array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
!
!  Should check if this is a "gather" pattern
!
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r8' )
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                      modc=grid%modc_gather )
      CALL mp_recvirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                      modc=grid%modc_gather )
     endif  !  grid%iam .lt. grid%npes_xy
#else
      bufres(1:lenarr) = arr(1:lenarr)
#endif
      return
   end subroutine fv_gather_r8

!====================================================================================

   subroutine fv_read_r8(iu, pattern_type, arr, lenarr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      integer, intent(in)                :: iu           ! Logical unit
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      real(r8), intent(out)              :: arr(*)       ! Array to be gathered
#else
      real(r8), intent(out)              :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
      real(r8), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r8' )
!
!  Should check if this is a "scatter" pattern
!

      if (grid%iam == 0) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write(iulog,*) 'FV_READ_R8 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                      modc=grid%modc_scatter )
      CALL mp_recvirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                      modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy

      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write(iulog,*) 'FV_READ_R8 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine fv_read_r8


!====================================================================================

   subroutine fv_write_r8(iu, pattern_type, arr, lenarr)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      integer, intent(in)                :: iu           ! Logical unit
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer, intent(in)                :: lenarr       ! Global length of array
#if defined( SPMD )
      real(r8) arr(*)                                    ! Array to be gathered
#else
      real(r8) arr(lenarr)                               ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
      real(r8), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
!
!  Should check if this is a "gather" pattern
!
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r8' )

      if ( grid%iam == 0 ) then
         allocate( bufres(lenarr) ) 
      else
         allocate( bufres(1) )
      endif
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                      modc=grid%modc_gather )
      CALL mp_recvirr(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                      modc=grid%modc_gather )
     endif  !  grid%iam .lt. grid%npes_xy

      if (grid%iam == 0) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write(iulog,*) 'FV_WRITE_R8 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun ('FV_WRITE_R8')
         end if
      endif
      deallocate( bufres )
#else
      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write(iulog,*) 'fv_write_r8 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun ('FV_WRITE_R8')
      end if
#endif
      return
   end subroutine fv_write_r8

!====================================================================================

   subroutine fv_scatter_r4(pattern_type, bufres, lenarr, arr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      real(r4) bufres(*)                                 ! Array to scatter
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      real(r4), intent(out)              :: arr(*)       ! Array to be gathered
#else
      real(r4), intent(out)              :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r4' )
!
!  Should check if this is a "scatter" pattern
!
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
      CALL mp_recvirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy
#else 
      arr(1:lenarr) = bufres(1:lenarr)
#endif
      return
   end subroutine fv_scatter_r4

!====================================================================================

   subroutine fv_gather_r4(pattern_type, arr, lenarr, bufres)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r4) arr(*)            ! Array to be gathered
#else
      real(r4) arr(lenarr)       ! Array (SMP-only)
#endif
      real(r4), intent(out)              :: bufres(*)    ! Gathered array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
!
!  Should check if this is a "gather" pattern
!
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r4' )
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
      CALL mp_recvirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
     endif  !  grid%iam .lt. grid%npes_xy
#else
      bufres(1:lenarr) = arr(1:lenarr)
#endif
      return
   end subroutine fv_gather_r4

!====================================================================================

   subroutine fv_read_r4(iu, pattern_type, arr, lenarr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      integer, intent(in)                :: iu           ! Logical unit
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      real(r4), intent(out)              :: arr(*)       ! Array to be gathered
#else
      real(r4), intent(out)              :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
      real(r4), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r4' )
!
!  Should check if this is a "scatter" pattern
!

      if (grid%iam == 0) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write(iulog,*) 'FV_READ_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
      CALL mp_recvirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy

      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write(iulog,*) 'FV_READ_R4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine fv_read_r4


!====================================================================================

   subroutine fv_write_r4(iu, pattern_type, arr, lenarr)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      integer, intent(in)                :: iu           ! Logical unit
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer, intent(in)                :: lenarr       ! Global length of array
#if defined( SPMD )
      real(r4) arr(*)                                    ! Array to be gathered
#else
      real(r4) arr(lenarr)                               ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
      real(r4), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
!
!  Should check if this is a "gather" pattern
!
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'r4' )

      if ( grid%iam == 0 ) then
         allocate( bufres(lenarr) ) 
      else
         allocate( bufres(1) )
      endif
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
      CALL mp_recvirr_r4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
     endif  !  grid%iam .lt. grid%npes_xy

      if (grid%iam == 0) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write(iulog,*) 'FV_WRITE_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun ('FV_WRITE_R4')
         end if
      endif
      deallocate( bufres )
#else
      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write(iulog,*) 'fv_write_r4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun ('FV_WRITE_R4')
      end if
#endif
      return
   end subroutine fv_write_r4

!====================================================================================

   subroutine fv_scatter_i4(pattern_type, bufres, lenarr, arr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer bufres(*)                                 ! Array to scatter
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      integer, intent(out)              :: arr(*)       ! Array to be gathered
#else
      integer, intent(out)              :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'i4' )
!
!  Should check if this is a "scatter" pattern
!
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
      CALL mp_recvirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy
#else 
      arr(1:lenarr) = bufres(1:lenarr)
#endif
      return
   end subroutine fv_scatter_i4

!====================================================================================

   subroutine fv_gather_i4(pattern_type, arr, lenarr, bufres)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer lenarr             ! Global size of array
#if defined( SPMD )
      integer arr(*)             ! Array to be gathered
#else
      integer arr(lenarr)        ! Array (SMP-only)
#endif
      integer, intent(out)               :: bufres(*)    ! Gathered array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
!
!  Should check if this is a "gather" pattern
!
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'i4' )
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
      CALL mp_recvirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, arr, bufres,   &
                         modc=grid%modc_gather )
     endif  !  grid%iam .lt. grid%npes_xy
#else
      bufres(1:lenarr) = arr(1:lenarr)
#endif
      return
   end subroutine fv_gather_i4

!====================================================================================

   subroutine fv_read_i4(iu, pattern_type, arr, lenarr)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      integer, intent(in)                :: iu           ! Logical unit
      character(len=*), intent(in)       :: pattern_type ! Type of comm pattern
      integer, intent(in)                :: lenarr       ! Global size of array
#if defined( SPMD )
      integer, intent(out)               :: arr(*)       ! Array to be gathered
#else
      integer, intent(out)               :: arr(lenarr)  ! Array (SMP-only)
#endif
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      type (parpatterntype):: pattern ! Pattern descriptor
      integer, allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      grid => get_dyn_state_grid()
      pattern =  get_pattern( grid, pattern_type, 'i4' )
!
!  Should check if this is a "scatter" pattern
!

      if (grid%iam == 0) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write(iulog,*) 'FV_READ_I4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
     if (grid%iam .lt. grid%npes_xy) then
      CALL mp_sendirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
      CALL mp_recvirr_i4(pattern%Comm, pattern%SendDesc, pattern%RecvDesc, bufres, arr,   &
                         modc=grid%modc_scatter )
     endif  !  grid%iam .lt. grid%npes_xy

      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write(iulog,*) 'FV_READ_I4 ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine fv_read_i4

!====================================================================================

   function get_decomp( grid, type )
      ! this utility function probably belongs where the decomptype objects are stored
      type (T_FVDYCORE_GRID), pointer :: grid  ! grid information
      character(len=*), intent(in) :: type
      type (decomptype)            :: get_decomp ! Decomposition descriptor
      select case (type)
      case ('2d')
         get_decomp = grid%strip2d
      case ('3dxzy')
         get_decomp = grid%strip3dxzy
      case ('3dxzyp')
         get_decomp = grid%strip3dxzyp
      case ('3dxyz')
         get_decomp = grid%strip3dxyz
      case default
         write(iulog,*)'get_decomp: invalid number decomposition type=', type
         call endrun ()
      end select
   end function get_decomp

!====================================================================================
#if defined( SPMD )
   function get_pattern( grid, type, data_type )

      ! this utility function probably belongs where the decomptype objects are stored
      type (T_FVDYCORE_GRID) :: grid  ! grid information
      character(len=*), intent(in) :: type
      character(len=*), intent(in) :: data_type
      type (ParPatternType) :: get_pattern ! grid information
      select case (type)
      case ('s_2dxy')
         select case (data_type)
         case ('r8')
            get_pattern = grid%s_2dxy_r8
         case ('r4')
            get_pattern = grid%s_2dxy_r4
         case ('i4')
            get_pattern = grid%s_2dxy_i4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select

      case ('s_3dxyz')
         select case (data_type)
         case ('r8')
            get_pattern = grid%s_3dxyz_r8
         case ('r4')
            get_pattern = grid%s_3dxyz_r4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select

      case ('s_3dxyzp')
         select case (data_type)
         case ('r8')
            get_pattern = grid%s_3dxyzp_r8
         case ('r4')
            get_pattern = grid%s_3dxyzp_r4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select

      case ('g_2dxy')
         select case (data_type)
         case ('r8')
            get_pattern = grid%g_2dxy_r8
         case ('r4')
            get_pattern = grid%g_2dxy_r4
         case ('i4')
            get_pattern = grid%g_2dxy_i4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select
      case ('g_3dxyz')
         select case (data_type)
         case ('r8')
            get_pattern = grid%g_3dxyz_r8
         case ('r4')
            get_pattern = grid%g_3dxyz_r4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select
      case ('g_3dxyzp')
         select case (data_type)
         case ('r8')
            get_pattern = grid%g_3dxyzp_r8
         case ('r4')
            get_pattern = grid%g_3dxyzp_r4
         case default
            write(iulog,*)'get_pattern: ', data_type, ' not supported in decomposition ', type
            call endrun ()
         end select
      case default
         write(iulog,*)'get_pattern: invalid number decomposition type=', type
         call endrun ()
      end select
   end function get_pattern
#endif
end module io_dist
