
module spmd_dyn
!BOP
!
! !MODULE: Subroutines to initialize SPMD implementation of CAM
!

#if (defined SPMD)

!
! !USES:
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: iam, masterproc, npes
   use pmgrid,             only: plat, plon, numbnd, &
                                 numlats, beglat, endlat, &
                                 plev, beglev, endlev, endlevp1, &
                                 endlevp, myid_y, myid_z, npr_y, npr_z, plevp, &
                                 myidxy_x, myidxy_y, nprxy_x, nprxy_y, &
                                 beglonxy, endlonxy, beglatxy, endlatxy, &
                                 twod_decomp, spmd_on, mod_transpose, mod_geopk, &
                                 mod_gatscat
   use mpishorthand,       only: mpir8, mpicom, mpiint, mpi_success
   use decompmodule,       only: decomptype, decompcreate
   use ghostmodule,        only: ghosttype
   use parutilitiesmodule, only: parpatterntype
   use fv_control_mod,     only: ct_overlap, trac_decomp
   use infnan,             only: inf
   use abortutils,         only: endrun
   use cam_logfile,        only: iulog

   implicit none

! !PUBLIC MEMBER FUNCTIONS:

   public spmdinit_dyn, decomp_wavenumbers
   public compute_gsfactors, spmdbuf

! !PUBLIC DATA MEMBERS:

   logical :: local_dp_map=.false.    ! flag indicates that mapping between dynamics 
                                      !  and physics decompositions does not require 
                                      !  interprocess communication
   integer :: block_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in dynamics decomposition (including level 0)
   integer :: chunk_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in physics decomposition (including level 0)

   integer :: proc(plat)              ! processor id associated with a given lat.
   integer, allocatable :: cut(:,:)   ! partition for MPI tasks
   integer, allocatable :: nlat_p(:)  ! number of latitudes per subdomain
   integer, allocatable :: kextent(:) ! number of levels per subdomain

   integer comm_y            ! communicator in latitude
   integer comm_z            ! communicator in vertical
   integer commxy_x          ! communicator in longitude (xy second. decomp.)
   integer commxy_y          ! communicator in latitude (xy second. decomp.)
   integer mpicom_yz         ! communicator for yz decomposition
   integer mpicom_nyz        ! communicator for multiple yz decomposition
   integer mpicom_xy         ! communicator for xy decomposition
   integer npes_yz           ! number of processes for yz decomposition
   integer npes_xy           ! number of processes for xy decomposition
   integer, allocatable :: lonrangexy(:,:)   ! global xy-longitude subdomain index
   integer, allocatable :: latrangexy(:,:)   ! global xy-latitude subdomain index
   logical geopkdist         ! use a distributed method for geopotential calculation 
                             !  with 2D decomp.
   logical geopk16byte       ! use Z-parallel distributed method for geopotential 
                             !  calculation with 2D decomp.; otherwise, use Z-serial 
                             !  pipeline algorithm
   integer geopkblocks       ! number of stages to use in Z-serial non-transpose
                             !  geopotential method (routine geopk_d)
                             !  with 2D decomp.
   integer modc_sw_dynrun    ! mod_comm irregular underlying communication method for dyn_run/misc
                             !  0 for original mp_sendirr/mp_recvirr
                             !  1 for mp_swapirr and point-to-point communications
                             !  2 for mp_swapirr and all-to-all communications
   logical modc_hs_dynrun    ! true for mod_comm irregular communication handshaking for dyn_run/misc
   logical modc_send_dynrun  ! true for mod_comm irregular communication blocking send for
                             !  dyn_run/misc, false for nonblocking send
   integer modc_mxreq_dynrun ! maximum number of nonblocking communication requests to allow
                             !  when using mp_swapirr and point-to-point communications for
                             !  dyn_run/misc
                             !  < 0 implies no limits
   integer modc_sw_cdcore    ! mod_comm irregular underlying communication method for cd_core/geopk
                             !  0 for original mp_sendirr/mp_recvirr
                             !  1 for mp_swapirr and point-to-point communications
                             !  2 for mp_swapirr and all-to-all communications
   logical modc_hs_cdcore    ! true for mod_comm irregular communication handshaking for cd_core/geopk
   logical modc_send_cdcore  ! true for geopk_d or mod_comm irregular communication blocking send for
                             !  cd_core/geopk, false for nonblocking send
   integer modc_mxreq_cdcore ! maximum number of nonblocking communication requests to allow
                             !  when using mp_swapirr and point-to-point communications for
                             !  cd_core/geopk
                             !  < 0 implies no limits
   integer modc_sw_gather    ! mod_comm irregular underlying communication method for gather
                             !  0 for original mp_sendirr/mp_recvirr
                             !  1 for mp_swapirr and point-to-point communications
                             !  2 for mp_swapirr and all-to-all communications
   logical modc_hs_gather    ! true for mod_comm irregular communication handshaking for gather
   logical modc_send_gather  ! true for mod_comm irregular communication blocking send for
                             !  gather, false for nonblocking send
   integer modc_mxreq_gather ! maximum number of nonblocking communication requests to allow
                             !  when using mp_swapirr and point-to-point communications for
                             !  gather
                             !  < 0 implies no limits
   integer modc_sw_scatter   ! mod_comm irregular underlying communication method for scatter
                             !  0 for original mp_sendirr/mp_recvirr
                             !  1 for mp_swapirr and point-to-point communications
                             !  2 for mp_swapirr and all-to-all communications
   logical modc_hs_scatter   ! true for mod_comm irregular communication handshaking for scatter
   logical modc_send_scatter ! true for mod_comm irregular communication blocking send for
                             !  scatter, false for nonblocking send
   integer modc_mxreq_scatter! maximum number of nonblocking communication requests to allow
                             !  when using mp_swapirr and point-to-point communications for
                             !  scatter
                             !  < 0 implies no limits
   integer modc_sw_tracer    ! mod_comm irregular underlying communication method for multiple tracers
                             !  0 for original mp_sendirr/mp_recvirr
                             !  1 for mp_swapirr and point-to-point communications
                             !  2 for mp_swapirr and all-to-all communications
   logical modc_hs_tracer    ! true for mod_comm irregular communication handshaking for multiple tracers
   logical modc_send_tracer  ! true for mod_comm irregular communication blocking send for
                             !  multiple tracers, false for nonblocking send
   integer modc_mxreq_tracer ! maximum number of nonblocking communication requests to allow
                             !  when using mp_swapirr and point-to-point communications for
                             !  multiple tracers
                             !  < 0 implies no limits
   integer modc_onetwo       ! one or two simultaneous mod_comm irregular communications (excl. tracers)
   integer modc_tracers      ! max number of tracers for simultaneous mod_comm irregular communications 
                             !  0 for original mp_sendirr/mp_recvirr communications
                             !  positive for special tracer routines

   type (ghosttype), save  :: ghostpe_yz, ghostpe1_yz
   type (parpatterntype)   :: ikj_xy_to_yz, ijk_yz_to_xy, ijk_xy_to_yz, &
                              pexy_to_pe, pkxy_to_pkc
!
! !DESCRIPTION: 
!   {\bf Purpose:} Subroutines to initialize SPMD implementation of CAM
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Alterations for LR SPMD mode
!   01.05.09  Mirin              2-D yz decomposition
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.12.20  Sawyer             Changed index order of Q3 decomposition
!   02.12.11  Sawyer             Use parbegin/endtransfer for transposes
!   03.05.07  Sawyer             Removed unneeded decompositions
!   06.03.01  Sawyer             Removed tracertrans-related variables
!
!EOP
!-----------------------------------------------------------------------

contains

  subroutine spmd_dyn_defaultopts(npr_yz_out, geopktrans_out,    &
               geopkblocks_out,                                  &
               force_2d_out, modcomm_transpose_out,              &
               modcomm_geopk_out, modcomm_gatscat_out,           &
               dyn_alltoall_out, dyn_allgather_out,              &
               dyn_equi_by_col_out,                              &
               dyn_npes_out, dyn_npes_stride_out,                &
               modc_sw_dynrun_out, modc_hs_dynrun_out,           &
               modc_send_dynrun_out, modc_mxreq_dynrun_out,      &
               modc_sw_cdcore_out, modc_hs_cdcore_out,           &
               modc_send_cdcore_out, modc_mxreq_cdcore_out,      &
               modc_sw_gather_out, modc_hs_gather_out,           &
               modc_send_gather_out, modc_mxreq_gather_out,      &
               modc_sw_scatter_out, modc_hs_scatter_out,         &
               modc_send_scatter_out, modc_mxreq_scatter_out,    &
               modc_sw_tracer_out, modc_hs_tracer_out,           &
               modc_send_tracer_out, modc_mxreq_tracer_out,      &
               modc_onetwo_out, modc_tracers_out                 )
!----------------------------------------------------------------------
! Purpose: Return default runtime options
! Author: Art Mirin
!----------------------------------------------------------------------
!------------------------------Arguments-------------------------------
     ! yz and xy decompositions
     integer, intent(out), optional :: npr_yz_out(4)
     ! geopotential method (routine geopk, geopk16, or geopk_d)
     integer, intent(out), optional :: geopktrans_out
     ! number of stages to use in geopotential method geopk_d
     integer, intent(out), optional :: geopkblocks_out
     ! option to force transpose computation for 1D decomp.
     integer, intent(out), optional :: force_2d_out
! Original mod_comm irregular communication options
     ! mod_comm transpose method
     integer, intent(out), optional :: modcomm_transpose_out
     ! mod_comm geopk method
     integer, intent(out), optional :: modcomm_geopk_out
     ! mod_comm gather/scatter method
     integer, intent(out), optional :: modcomm_gatscat_out
! EUL/SLD-only arguments
     integer, intent(out), optional :: dyn_alltoall_out
     integer, intent(out), optional :: dyn_allgather_out
     logical, intent(out), optional :: dyn_equi_by_col_out
     integer, intent(out), optional :: dyn_npes_out
     integer, intent(out), optional :: dyn_npes_stride_out
! Additional mod_comm irregular communication options
     integer, intent(out), optional :: modc_sw_dynrun_out
     logical, intent(out), optional :: modc_hs_dynrun_out
     logical, intent(out), optional :: modc_send_dynrun_out
     integer, intent(out), optional :: modc_mxreq_dynrun_out
     integer, intent(out), optional :: modc_sw_cdcore_out
     logical, intent(out), optional :: modc_hs_cdcore_out
     logical, intent(out), optional :: modc_send_cdcore_out
     integer, intent(out), optional :: modc_mxreq_cdcore_out
     integer, intent(out), optional :: modc_sw_gather_out
     logical, intent(out), optional :: modc_hs_gather_out
     logical, intent(out), optional :: modc_send_gather_out
     integer, intent(out), optional :: modc_mxreq_gather_out
     integer, intent(out), optional :: modc_sw_scatter_out
     logical, intent(out), optional :: modc_hs_scatter_out
     logical, intent(out), optional :: modc_send_scatter_out
     integer, intent(out), optional :: modc_mxreq_scatter_out
     integer, intent(out), optional :: modc_sw_tracer_out
     logical, intent(out), optional :: modc_hs_tracer_out
     logical, intent(out), optional :: modc_send_tracer_out
     integer, intent(out), optional :: modc_mxreq_tracer_out
     integer, intent(out), optional :: modc_onetwo_out
     integer, intent(out), optional :: modc_tracers_out

!----------------------------------------------------------------------
     if (present(npr_yz_out) ) then
        npr_yz_out(1) = npes
        npr_yz_out(2) = 1
        npr_yz_out(3) = 1
        npr_yz_out(4) = npes
     endif
     if (present(geopktrans_out) ) then
        geopktrans_out = 0
     endif
     if (present(geopkblocks_out) ) then
        geopkblocks_out = 1
     endif
     if (present(force_2d_out) ) then
        force_2d_out = 0
     endif
     if (present(modcomm_transpose_out) ) then
        modcomm_transpose_out = 0
     endif
     if (present(modcomm_geopk_out) ) then
        modcomm_geopk_out = 0
     endif
     if (present(modcomm_gatscat_out) ) then
        modcomm_gatscat_out = 0
     endif

     ! dynrun: handshaking and send
     if (present(modc_sw_dynrun_out) ) then
        modc_sw_dynrun_out = 0
     endif
     if (present(modc_hs_dynrun_out) ) then
        modc_hs_dynrun_out = .true.
     endif
     if (present(modc_send_dynrun_out) ) then
        modc_send_dynrun_out = .true.
     endif
     if (present(modc_mxreq_dynrun_out) ) then
        modc_mxreq_dynrun_out = -1
     endif

     ! cd_core: handshaking and send
     if (present(modc_sw_cdcore_out) ) then
        modc_sw_cdcore_out = 0
     endif
     if (present(modc_hs_cdcore_out) ) then
        modc_hs_cdcore_out = .true.
     endif
     if (present(modc_send_cdcore_out) ) then
        modc_send_cdcore_out = .true.
     endif
     if (present(modc_mxreq_cdcore_out) ) then
        modc_mxreq_cdcore_out = -1
     endif

     ! gather: handshaking and mxreq (and swap)
     if (present(modc_sw_gather_out) ) then
        modc_sw_gather_out = 1
     endif
     if (present(modc_hs_gather_out) ) then
        modc_hs_gather_out = .true.
     endif
     if (present(modc_send_gather_out) ) then
        modc_send_gather_out = .true.
     endif
     if (present(modc_mxreq_gather_out) ) then
        modc_mxreq_gather_out = 64
     endif

     ! scatter: no restrictions
     if (present(modc_sw_scatter_out) ) then
        modc_sw_scatter_out = 0
     endif
     if (present(modc_hs_scatter_out) ) then
        modc_hs_scatter_out = .false.
     endif
     if (present(modc_send_scatter_out) ) then
        modc_send_scatter_out = .true.
     endif
     if (present(modc_mxreq_scatter_out) ) then
        modc_mxreq_scatter_out = -1
     endif

     ! tracer: handshaking and send
     if (present(modc_sw_tracer_out) ) then
        modc_sw_tracer_out = 0
     endif
     if (present(modc_hs_tracer_out) ) then
        modc_hs_tracer_out = .true.
     endif
     if (present(modc_send_tracer_out) ) then
        modc_send_tracer_out = .true.
     endif
     if (present(modc_mxreq_tracer_out) ) then
        modc_mxreq_tracer_out = -1
     endif
     if (present(modc_onetwo_out) ) then
        modc_onetwo_out = 2
     endif
     if (present(modc_tracers_out) ) then
        modc_tracers_out = 3
     endif

     return

  end subroutine spmd_dyn_defaultopts

  subroutine spmd_dyn_setopts(npr_yz_in, geopktrans_in,       &
               geopkblocks_in,                                &
               force_2d_in, modcomm_transpose_in,             &
               modcomm_geopk_in, modcomm_gatscat_in,          &
               dyn_alltoall_in, dyn_allgather_in,             &
               dyn_equi_by_col_in,                            &
               dyn_npes_in, dyn_npes_stride_in,               &
               modc_sw_dynrun_in, modc_hs_dynrun_in,          &
               modc_send_dynrun_in, modc_mxreq_dynrun_in,     &
               modc_sw_cdcore_in, modc_hs_cdcore_in,          &
               modc_send_cdcore_in, modc_mxreq_cdcore_in,     &
               modc_sw_gather_in, modc_hs_gather_in,          &
               modc_send_gather_in, modc_mxreq_gather_in,     &
               modc_sw_scatter_in, modc_hs_scatter_in,        &
               modc_send_scatter_in, modc_mxreq_scatter_in,   &
               modc_sw_tracer_in, modc_hs_tracer_in,          &
               modc_send_tracer_in, modc_mxreq_tracer_in,     &
               modc_onetwo_in, modc_tracers_in                )
!----------------------------------------------------------------------
! Purpose: Set runtime options
! Author: Art Mirin
!----------------------------------------------------------------------
!------------------------------Arguments-------------------------------
! yz and xy decompositions (npr_y, npr_z, nprxy_x, nprxy_y)
     integer, intent(in), optional :: npr_yz_in(4)
! geopotential method (routines geopk, geopk16, and geopk_d)
! 0 for transpose method, 1 for method using semi-global z communication 
!   with optional 16-byte arithmetic, 2 for method using local
!   z communication; method 0, method 1 with 16-byte arithmetic and 
!   method 2 are all bit-for-bit across decompositions; method 0
!   scales better than method 1 with npr_z, and method 1 is superior 
!   to method 0 for small npr_z. The optimum speed is attained either 
!   using method 1 with 8-byte arithmetic (standard for geopk16) or 
!   method 2 when utilizing the optimal value for the associated 
!   parameter geopkblocks; see geopk.F90.
     integer, intent(in), optional :: geopktrans_in
! number of stages to use in geopotential method geopk_d
     integer, intent(in), optional :: geopkblocks_in
! option to force transpose computation for 1D decomp.
! the only purpose for invoking this option is debugging
     integer, intent(in), optional :: force_2d_in
! mod_comm transpose/geopk/gatscat method
!   0 for temporary contiguous buffers
!   1 for mpi derived types
     integer, intent(in), optional :: modcomm_transpose_in, modcomm_geopk_in, &
                                      modcomm_gatscat_in
! Additional mod_comm irregular communication options
     integer, intent(in), optional :: modc_sw_dynrun_in
     logical, intent(in), optional :: modc_hs_dynrun_in
     logical, intent(in), optional :: modc_send_dynrun_in
     integer, intent(in), optional :: modc_mxreq_dynrun_in
     integer, intent(in), optional :: modc_sw_cdcore_in
     logical, intent(in), optional :: modc_hs_cdcore_in
     logical, intent(in), optional :: modc_send_cdcore_in
     integer, intent(in), optional :: modc_mxreq_cdcore_in
     integer, intent(in), optional :: modc_sw_gather_in
     logical, intent(in), optional :: modc_hs_gather_in
     logical, intent(in), optional :: modc_send_gather_in
     integer, intent(in), optional :: modc_mxreq_gather_in
     integer, intent(in), optional :: modc_sw_scatter_in
     logical, intent(in), optional :: modc_hs_scatter_in
     logical, intent(in), optional :: modc_send_scatter_in
     integer, intent(in), optional :: modc_mxreq_scatter_in
     integer, intent(in), optional :: modc_sw_tracer_in
     logical, intent(in), optional :: modc_hs_tracer_in
     logical, intent(in), optional :: modc_send_tracer_in
     integer, intent(in), optional :: modc_mxreq_tracer_in
     integer, intent(in), optional :: modc_onetwo_in
     integer, intent(in), optional :: modc_tracers_in
! EUL/SLD-only arguments
     integer, intent(in), optional :: dyn_alltoall_in
     integer, intent(in), optional :: dyn_allgather_in
     logical, intent(in), optional :: dyn_equi_by_col_in
     integer, intent(in), optional :: dyn_npes_in
     integer, intent(in), optional :: dyn_npes_stride_in
!----------------------------------------------------------------------
     integer omp_get_num_threads
     integer color, ierror, ntemp
!----------------------------------------------------------------------
     npes_yz = npes
     npes_xy = npes
     if (present(npr_yz_in) ) then
        npr_y   = npr_yz_in(1)
        npr_z   = npr_yz_in(2)
        nprxy_x = npr_yz_in(3)
        nprxy_y = npr_yz_in(4)
        npes_yz = npr_y*npr_z
        npes_xy = nprxy_x*nprxy_y
        if (masterproc) then
           write(iulog,*) 'npr_y = ', npr_y, '  npr_z = ', npr_z
           write(iulog,*) 'nprxy_x = ', nprxy_x, '  nprxy_y = ', nprxy_y
           write(iulog,*) 'npes = ', npes, '  npes_yz= ', npes_yz, '  npes_xy = ', npes_xy
        endif
        if (npes_yz > npes) then
           call endrun ('SPMD_DYN_SET : incorrect yz domain decomposition - aborting')
        endif
        if (npes_xy > npes) then
           call endrun ('SPMD_DYN_SET : incorrect xy domain decomposition - aborting')
        endif
        if (npes_xy < npes) then
           if (masterproc) then
              write(iulog,*) 'WARNING - proceeding with auxiliary dynamics processes'
           endif
        endif
        if (npes_yz < npes_xy) then
           if (masterproc) then
              write(iulog,*) 'WARNING - proceeding with smaller yz decomposition'
           endif
        endif
     else
        npr_y   = npes
        npr_z   = 1
        nprxy_x = 1
        nprxy_y = npes
        if (masterproc) then
           write(iulog,*) 'WARNING : npr_yz not present - using 1-D domain decomposition'
        endif
        npes_yz = npes
        npes_xy = npes
     endif
     if (ct_overlap .ne. 0) then
        if (npes .lt. 2*npes_yz) then
           call endrun ('SPMDINIT_SETOPTS: Not enough processes to overlap cd_core and trac2d')
        else
           if (masterproc) then
              write(iulog,*) 'Overlapping tracer and dynamics subcycles'
           endif
        endif
     endif
     if (trac_decomp .le. 0) then
        call endrun ('SPMDINIT_SETOPTS: trac_decomp improperly initialized')
     endif
     if (npes .lt. trac_decomp*npes_yz) then
        call endrun ('SPMDINIT_SETOPTS: Not enough processes to decompose tracers ')
     else
        if (masterproc) then
           write(iulog,*) 'Decomposing tracers into ', trac_decomp, ' groups'
        endif
     endif
     if (ct_overlap .gt. 0 .and. trac_decomp .gt. 1) then
        call endrun ('SPMDINIT_SETOPTS: Cannot simultaneously overlap cd_core/trac2d and decompose tracers')
     endif
     myid_z   = iam/npr_y
     myid_y   = iam - myid_z*npr_y
     color = iam/npes_yz
     call mpi_comm_split(mpicom, color, iam, mpicom_yz, ierror)
     if (ierror /= mpi_success) then
        write(iulog,*) 'SPMD_DYN_SETOPTS:  ERROR:  mpi_comm_split_yz failed with IER=', ierror
        call endrun
     endif
     call mpi_comm_size(mpicom_yz, ntemp, ierror)
     if (masterproc .and. ntemp .ne. npes_yz) then
        write(iulog,*) 'SPMD_DYN_SETOPTS:  ERROR:  mpicom_yz has incorrect size of ', ntemp
     endif
     if (ct_overlap .gt. 0 .or. trac_decomp .gt. 1) then
! These are mutually exclusive options
        if ((ct_overlap .gt. 0 .and. iam .lt. 2*npes_yz) .or.         &
           (trac_decomp .gt. 1 .and. iam .lt. trac_decomp*npes_yz)) then
              color = 1
        else
              color = 0
        endif
        call mpi_comm_split(mpicom, color, iam, mpicom_nyz, ierror)
        if (ierror /= mpi_success) then
           write (iulog,*) 'SPMD_DYN_SETOPTS:  ERROR:  mpi_comm_split_nyz failed with IER=', ierror
           call endrun
        endif
     else
        mpicom_nyz = mpicom_yz
     endif
     myidxy_y = iam/nprxy_x
     myidxy_x = iam - myidxy_y*nprxy_x
     color = iam/npes_xy
     call mpi_comm_split(mpicom, color, iam, mpicom_xy, ierror)
     if (ierror /= mpi_success) then
        write(iulog,*) 'SPMD_DYN_SETOPTS:  ERROR:  mpi_comm_split_xy failed with IER=', ierror
        call endrun
     endif
     call mpi_comm_size(mpicom_xy, ntemp, ierror)
     if (ntemp .ne. npes_xy) then
        write(iulog,*) 'SPMD_DYN_SETOPTS:  ERROR:  mpicom_xy has incorrect size of ', ntemp
     endif

     geopkdist   = .false.
     geopk16byte = .false.
     if (present(geopktrans_in) ) then
        if (geopktrans_in .ne. 0) geopkdist   = .true.
        if (geopktrans_in .eq. 1) geopk16byte = .true.
#ifdef NO_CRAY_POINTERS
        if (geopk16byte) then
           call endrun ('SPMD_DYN_SET : cannot use geopk16 unless compiler supports cray pointers')
        end if
#endif
        if (masterproc) then
           write(iulog,*) 'non-transpose geopk communication method = ', geopkdist
           write(iulog,*) 'Z-parallel non-transpose geopk communication method = ', geopk16byte
        endif
     else
        if (masterproc) then
           write(iulog,*) 'WARNING : geopktrans not present - using transpose method'
        endif
     endif

     if (present(geopkblocks_in) ) then
        geopkblocks = max(1,geopkblocks_in)
     else
        geopkblocks = 1
     endif
     if ((masterproc) .and. (geopkdist) .and. (.not. geopk16byte)) then
        write(iulog,*) 'number of stages in Z-serial non-transpose geopk method = ', geopkblocks
     endif

     twod_decomp = 1

     if (present(force_2d_in) ) then
        if (npr_z .eq. 1 .and. nprxy_x .eq. 1 .and. force_2d_in .eq. 0) then
           twod_decomp = 0
           if (masterproc) then
              write(iulog,*) 'decomposition is effectively 1D - skipping transposes'
           endif
        else
           if (masterproc) then
              write(iulog,*) 'using multi-2d decomposition methodology'
           endif
        endif
     else
        if (npr_z .eq. 1 .and. nprxy_x .eq. 1 ) twod_decomp = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : force_2d not present - defaulting'
        endif
     endif

     if (present(modcomm_transpose_in) ) then
        mod_transpose = modcomm_transpose_in
        if (masterproc) then
           write(iulog,*) 'modcomm transpose method = ', mod_transpose
        endif
     else
        mod_transpose = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modcomm_transpose not present - defaulting'
        endif
     endif

     if (present(modcomm_geopk_in) ) then
        mod_geopk = modcomm_geopk_in
        if (masterproc) then
           write(iulog,*) 'modcomm geopk method = ', mod_geopk
        endif
     else
        mod_geopk = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modcomm_geopk not present - defaulting'
        endif
     endif

     if (present(modcomm_gatscat_in) ) then
        mod_gatscat = modcomm_gatscat_in
        if (masterproc) then
           write(iulog,*) 'modcomm gatscat method = ', mod_gatscat
        endif
     else
        mod_gatscat= 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modcomm_gatscat not present - defaulting'
        endif
     endif

     if (present(modc_sw_dynrun_in) ) then
        modc_sw_dynrun = modc_sw_dynrun_in
        if (masterproc) then
           write(iulog,*) 'modc_sw_dynrun = ', modc_sw_dynrun
        endif
        if (modc_sw_dynrun .lt. 0 .or. modc_sw_dynrun .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_dynrun')
        endif
        if (modc_sw_dynrun .gt. 0 .and. mod_transpose .gt. 0) then
           modc_sw_dynrun = 0
           if (masterproc) then
              write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_dynrun reset to 0 for consistency'
           endif
        endif
     else
        modc_sw_dynrun = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_sw_dynrun not present - defaulting'
        endif
     endif

     if (present(modc_hs_dynrun_in) ) then
        modc_hs_dynrun = modc_hs_dynrun_in
        if (masterproc) then
           write(iulog,*) 'modc_hs_dynrun = ', modc_hs_dynrun
        endif
     else
        modc_hs_dynrun = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_hs_dynrun not present - defaulting'
        endif
     endif

     if (present(modc_send_dynrun_in) ) then
        modc_send_dynrun = modc_send_dynrun_in
        if (masterproc) then
           write(iulog,*) 'modc_send_dynrun = ', modc_send_dynrun
        endif
     else
        modc_send_dynrun = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_send_dynrun not present - defaulting'
        endif
     endif

     if (present(modc_mxreq_dynrun_in) ) then
        modc_mxreq_dynrun = modc_mxreq_dynrun_in
        if (masterproc) then
           write(iulog,*) 'modc_mxreq_dynrun = ', modc_mxreq_dynrun
        endif
     else
        modc_mxreq_dynrun = -1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_mxreq_dynrun not present - defaulting'
        endif
     endif

     if (present(modc_sw_cdcore_in) ) then
        modc_sw_cdcore = modc_sw_cdcore_in
        if (masterproc) then
           write(iulog,*) 'modc_sw_cdcore = ', modc_sw_cdcore
        endif
        if (modc_sw_cdcore .lt. 0 .or. modc_sw_cdcore .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_cdcore')
        endif
        if (modc_sw_cdcore .gt. 0 .and. (mod_transpose .gt. 0 .or. (mod_geopk .gt. 0 .and. geopk16byte))) then
           modc_sw_cdcore = 0
           if (masterproc) then
              write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_cdcore reset to 0 for consistency'
           endif
        endif
     else
        modc_sw_cdcore = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_sw_cdcore not present - defaulting'
        endif
     endif

     if (present(modc_hs_cdcore_in) ) then
        modc_hs_cdcore = modc_hs_cdcore_in
        if (masterproc) then
           write(iulog,*) 'modc_hs_cdcore = ', modc_hs_cdcore
        endif
     else
        modc_hs_cdcore = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_hs_cdcore not present - defaulting'
        endif
     endif

     if (present(modc_send_cdcore_in) ) then
        modc_send_cdcore = modc_send_cdcore_in
        if (masterproc) then
           write(iulog,*) 'modc_send_cdcore = ', modc_send_cdcore
        endif
     else
        modc_send_cdcore = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_send_cdcore not present - defaulting'
        endif
     endif

     if (present(modc_mxreq_cdcore_in) ) then
        modc_mxreq_cdcore = modc_mxreq_cdcore_in
        if (masterproc) then
           write(iulog,*) 'modc_mxreq_cdcore = ', modc_mxreq_cdcore
        endif
     else
        modc_mxreq_cdcore = -1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_mxreq_cdcore not present - defaulting'
        endif
     endif

     if (present(modc_sw_gather_in) ) then
        modc_sw_gather = modc_sw_gather_in
        if (masterproc) then
           write(iulog,*) 'modc_sw_gather = ', modc_sw_gather
        endif
        if (modc_sw_gather .lt. 0 .or. modc_sw_gather .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_gather')
        endif
        if (modc_sw_gather .gt. 0 .and. mod_gatscat .gt. 0) then
           modc_sw_gather = 0
           if (masterproc) then
              write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_gather reset to 0 for consistency'
           endif
        endif
     else
        modc_sw_gather = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_sw_gather not present - defaulting'
        endif
     endif

     if (present(modc_hs_gather_in) ) then
        modc_hs_gather = modc_hs_gather_in
        if (masterproc) then
           write(iulog,*) 'modc_hs_gather = ', modc_hs_gather
        endif
     else
        modc_hs_gather = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_hs_gather not present - defaulting'
        endif
     endif

     if (present(modc_send_gather_in) ) then
        modc_send_gather = modc_send_gather_in
        if (masterproc) then
           write(iulog,*) 'modc_send_gather = ', modc_send_gather
        endif
     else
        modc_send_gather = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_send_gather not present - defaulting'
        endif
     endif

     if (present(modc_mxreq_gather_in) ) then
        modc_mxreq_gather = modc_mxreq_gather_in
        if (masterproc) then
           write(iulog,*) 'modc_mxreq_gather = ', modc_mxreq_gather
        endif
     else
        modc_mxreq_gather = -1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_mxreq_gather not present - defaulting'
        endif
     endif

     if (present(modc_sw_scatter_in) ) then
        modc_sw_scatter = modc_sw_scatter_in
        if (masterproc) then
           write(iulog,*) 'modc_sw_scatter = ', modc_sw_scatter
        endif
        if (modc_sw_scatter .lt. 0 .or. modc_sw_scatter .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_scatter')
        endif
        if (modc_sw_scatter .gt. 0 .and. mod_gatscat .gt. 0) then
           modc_sw_scatter = 0
           if (masterproc) then
              write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_scatter reset to 0 for consistency'
           endif
        endif
     else
        modc_sw_scatter = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_sw_scatter not present - defaulting'
        endif
     endif

     if (present(modc_hs_scatter_in) ) then
        modc_hs_scatter = modc_hs_scatter_in
        if (masterproc) then
           write(iulog,*) 'modc_hs_scatter = ', modc_hs_scatter
        endif
     else
        modc_hs_scatter = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_hs_scatter not present - defaulting'
        endif
     endif

     if (present(modc_send_scatter_in) ) then
        modc_send_scatter = modc_send_scatter_in
        if (masterproc) then
           write(iulog,*) 'modc_send_scatter = ', modc_send_scatter
        endif
     else
        modc_send_scatter = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_send_scatter not present - defaulting'
        endif
     endif

     if (present(modc_mxreq_scatter_in) ) then
        modc_mxreq_scatter = modc_mxreq_scatter_in
        if (masterproc) then
           write(iulog,*) 'modc_mxreq_scatter = ', modc_mxreq_scatter
        endif
     else
        modc_mxreq_scatter = -1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_mxreq_scatter not present - defaulting'
        endif
     endif

     if (present(modc_sw_tracer_in) ) then
        modc_sw_tracer = modc_sw_tracer_in
        if (masterproc) then
           write(iulog,*) 'modc_sw_tracer = ', modc_sw_tracer
        endif
        if (modc_sw_tracer .lt. 0 .or. modc_sw_tracer .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_tracer')
        endif
        if (modc_sw_tracer .gt. 0 .and. mod_transpose .gt. 0) then
           modc_sw_tracer = 0
           if (masterproc) then
              write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_tracer reset to 0 for consistency'
           endif
        endif
     else
        modc_sw_tracer = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_sw_tracer not present - defaulting'
        endif
     endif

     if (present(modc_hs_tracer_in) ) then
        modc_hs_tracer = modc_hs_tracer_in
        if (masterproc) then
           write(iulog,*) 'modc_hs_tracer = ', modc_hs_tracer
        endif
     else
        modc_hs_tracer = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_hs_tracer not present - defaulting'
        endif
     endif

     if (present(modc_send_tracer_in) ) then
        modc_send_tracer = modc_send_tracer_in
        if (masterproc) then
           write(iulog,*) 'modc_send_tracer = ', modc_send_tracer
        endif
     else
        modc_send_tracer = .true.
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_send_tracer not present - defaulting'
        endif
     endif

     if (present(modc_mxreq_tracer_in) ) then
        modc_mxreq_tracer = modc_mxreq_tracer_in
        if (masterproc) then
           write(iulog,*) 'modc_mxreq_tracer = ', modc_mxreq_tracer
        endif
     else
        modc_mxreq_tracer = -1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_mxreq_tracer not present - defaulting'
        endif
     endif

     if (present(modc_onetwo_in) ) then
        modc_onetwo = modc_onetwo_in
        if (masterproc) then
           write(iulog,*) 'modc_onetwo = ', modc_onetwo
        endif
        if (modc_onetwo .lt. 1 .or. modc_onetwo .gt. 2) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_onetwo')
        endif
     else
        modc_onetwo = 1
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_onetwo not present - defaulting'
        endif
     endif

     if (present(modc_tracers_in) ) then
        modc_tracers = modc_tracers_in
        if (masterproc) then
           write(iulog,*) 'modc_tracers = ', modc_tracers
        endif
        if (modc_tracers .lt. 0) then
           call endrun ('SPMD_DYN_SET : inadmissable value of modc_tracers')
        endif
     else
        modc_tracers = 0
        if (masterproc) then
           write(iulog,*) 'WARNING : modc_tracers not present - defaulting'
        endif
     endif

     return

  end subroutine spmd_dyn_setopts

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: spmdinit_dyn --- SPMD initialization for dynamics
!
! !INTERFACE:

   subroutine spmdinit_dyn ()

! !USES:
      use parutilitiesmodule, only : parinit, parsplit
      use decompmodule, only : decompcreate

! !DESCRIPTION:
!
!   SPMD initialization routine: get number of cpus, processes, tids, etc
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Added LR-specific initialization
!   01.03.26  Sawyer             Added ProTeX documentation
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.10.16  Sawyer             Added Y at each Z decompositions
!   03.07.22  Sawyer             Removed decomps used by highp2
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer procid    ! processor id
      integer procids   ! processor id SH
      integer procidn   ! processor id NH
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer workleft  ! amount of work still to be parcelled out
      integer actual    ! actual amount of work parcelled out
      integer ideal     ! ideal amt of work to parcel out
      integer pesleft   ! number of procs still to be given work
      integer isum      ! running total of work parcelled out
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer xdist(1)  ! number of lons per subdomain
      integer, allocatable :: ydist(:) ! number of lats per subdomain
      integer, allocatable :: zdist(:) ! number of levels per subdomain
      integer, allocatable :: zdistq(:) ! number of levels per subdomain for Q3
      integer ier       ! error flag
      integer rank_y, size_y   !  rank and size wrt y-communicator
      integer rank_z, size_z   !  rank and size wrt z-communicator
      integer rankxy_x, sizexy_x   !  rank and size wrt xy x-communicator
      integer rankxy_y, sizexy_y   !  rank and size wrt xy y-communicator
      integer zdist1(1) ! used for misc. decomposition definitions
      integer, allocatable :: xdistxy(:) ! number of xy-longs per subdomain
      integer, allocatable :: ydistxy(:) ! number of xy-lats per subdomain
      integer, allocatable :: ydistqxy(:) ! number of xy tracer/lats per subdomain
      integer zdistxy(1)  ! number of xy-verts per subdomain
      integer j, k, vert, lonn
      integer ydistk(1)
      integer mod_maxirr

      spmd_on = 1

! Default 2D decomposition
      beglev = 1
      endlev = plev
      endlevp1 = plev + 1
      endlevp = plev + 1
      mod_maxirr = max(modc_onetwo, modc_tracers)
!
! Addition for LR dynamical core to initialize PILGRIM library
!
      call parinit(comm=mpicom, &
                   npryzxy = (/ npr_y, npr_z, nprxy_x, nprxy_y /), &
                   mod_method  = mod_transpose, &
                   mod_geopk   = mod_geopk,     &
                   mod_maxirr  = mod_maxirr,    &
                   mod_gatscat = mod_gatscat )
!
! Form separate communicators
!
      call parsplit(mpicom, myid_z, iam, comm_y, rank_y, size_y)
      call parsplit(mpicom, myid_y, iam, comm_z, rank_z, size_z)
      call parsplit(mpicom, myidxy_y, iam, commxy_x, rankxy_x, sizexy_x)
      call parsplit(mpicom, myidxy_x, iam, commxy_y, rankxy_y, sizexy_y)

!
!-----------------------------------------------------------------------
!
! Compute y decomposition
!
      allocate (ydist  (npr_y))
      allocate (nlat_p (0:npes-1))
      allocate (cut    (2,0:npes-1))

      ydist(:) = 0
      nlat_p(:) = 0
      cut(1,:) = -1
      cut(2,:) = -2

      lat = plat / npr_y
      workleft = plat - lat * npr_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 latitudes per subdomain')
      endif
!
! Be careful:  ydist is 1-based.  NCARs arrays, e.g., cut,  are 0-based
!
      do procid=1,npr_y
         ydist(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_y
            ydist(procids) = ydist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydist(procidn) = ydist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydist) /= plat ) then
         write(iulog,*)'SPMDINIT_DYN:', ydist,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(y) not zero.  Value is ',workleft
         call endrun
      end if

! Set the NCAR data structures

      lat  = 0
      do procid=0,npr_y-1
         cut(1,procid) = lat+1
         lat = lat + ydist(procid+1)
         cut(2,procid) = lat
         nlat_p(procid) = ydist(procid+1)

         if (masterproc) then
            write(iulog,*) 'nlat_p(',procid,') = ', nlat_p(procid)
         end if

         if (myid_y == procid) then
            beglat  = cut(1,myid_y)
            endlat  = cut(2,myid_y)
            numlats = ydist(procid+1)
         end if
      enddo

      do k = 1, npr_z-1
         do j = 0, npr_y-1
            procid = j + k*npr_y
            cut(1,procid) = cut(1,j)
            cut(2,procid) = cut(2,j)
            nlat_p(procid) = nlat_p(j)
         enddo
      enddo
!
! Compute z decomposition
!
      allocate (zdist ((npes-1)/npr_y+1))
      allocate (zdistq(npr_z))
      allocate (kextent(npr_z))

      zdist(:) = 0

      vert = plev / npr_z
      workleft = plev - vert * npr_z
      if ( vert < 1 ) then
         call endrun ('SPMDINIT_DYN: less than 1 verticals per subdomain')
      endif

      do procid=1,npr_z
         zdist(procid) = vert
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_z+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_z
            zdist(procids) = zdist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               zdist(procidn) = zdist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(zdist) /= plev ) then
         write(iulog,*)'SPMDINIT_DYN:', zdist,' does not add up to ', plev
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(z) not zero.  Value is ',workleft
         call endrun
      end if

! kextent is global, zdist is local to this module
      kextent(:) = zdist(:)

! Compute local limits

      beglev = 1
      endlev = zdist(1)
      do procid = 1, myid_z
         beglev = endlev + 1
         endlev = beglev + zdist(procid+1) - 1
      enddo
      endlevp1 = endlev + 1
      endlevp = endlev
      if (myid_z == npr_z-1) endlevp = endlev + 1

      if (iam .ge. npes_yz) then
! Auxiliary processes only
        beglat = 1
        endlat = 0
        numlats = 0
        beglev = 1
        endlev = 0
        endlevp = endlev + 1
        endlevp1 = endlev + 1
      endif

!
! Compute x secondary decomposition
!
      allocate (xdistxy (nprxy_x))

      xdistxy(:) = 0

      lonn = plon / nprxy_x
      workleft = plon - lonn * nprxy_x
      if ( lonn < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-longitudes per subdomain')
      endif

      do procid=1,nprxy_x
         xdistxy(procid) = lonn
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_x+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_x
            xdistxy(procids) = xdistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               xdistxy(procidn) = xdistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(xdistxy) /= plon ) then
         write(iulog,*)'SPMDINIT_DYN:', xdistxy,' does not add up to ', plon
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(xy-x) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglonxy = 1
      endlonxy = xdistxy(1)
      do procid = 1, myidxy_x
         beglonxy = endlonxy + 1
         endlonxy = beglonxy + xdistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (lonrangexy(2,nprxy_x))
      lonrangexy(1,1) = 1
      lonrangexy(2,1) = xdistxy(1)
      do procid = 2, nprxy_x
         lonrangexy(1,procid) = lonrangexy(2,procid-1) + 1
         lonrangexy(2,procid) = lonrangexy(1,procid) + xdistxy(procid) - 1
      enddo
!
! Compute y secondary decomposition
!
      allocate (ydistxy ((npes-1)/nprxy_x+1))

      ydistxy(:) = 0

      lat = plat / nprxy_y
      workleft = plat - lat * nprxy_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-latitudes per subdomain')
      endif

      do procid=1,nprxy_y
         ydistxy(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_y
            ydistxy(procids) = ydistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydistxy(procidn) = ydistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydistxy) /= plat ) then
         write(iulog,*)'SPMDINIT_DYN:', ydistxy,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(xy-y) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglatxy = 1
      endlatxy = ydistxy(1)
      do procid = 1, myidxy_y
         beglatxy = endlatxy + 1
         endlatxy = beglatxy + ydistxy(procid+1) - 1
      enddo

      if (iam .ge. npes_xy) then
! Auxiliary processes only
        beglonxy = 1
        endlonxy = 0
        beglatxy = 1
        endlatxy = 0
      endif

! Compute global table

      allocate (latrangexy(2,nprxy_y))
      latrangexy(1,1) = 1
      latrangexy(2,1) = ydistxy(1)
      do procid = 2, nprxy_y
         latrangexy(1,procid) = latrangexy(2,procid-1) + 1
         latrangexy(2,procid) = latrangexy(1,procid) + ydistxy(procid) - 1
      enddo

!
! Do generic NCAR decomposition
!
      proc(:) = 0
      do procid=0,npr_y*npr_z-1
         if (iam == 0) then
            write(iulog,*)'procid ',procid,' assigned ', &
                 cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                 cut(1,procid),' through ',cut(2,procid)
         endif
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
      end do

      nmostlat = plat
      smostlat = 1
      if (iam .lt. npes_yz) then

! Primary processes only
!
! Number of neighbor processors needed for boundary communication.  North
! first.
!
        nmostlat = 0
        isum = 0
        do procid=myid_y+1,npr_y-1
           nmostlat = cut(2,procid)
           isum = isum + cut(2,procid) - cut(1,procid) + 1
           if (isum >= numbnd) goto 20
        end do
20      if (myid_y /= npr_y-1 .and. isum < numbnd .and. nmostlat /= plat)then
           call endrun ('SPMDINIT_DYN: Something wrong in computation of northern neighbors')
        end if

        smostlat = 0
        isum = 0
        do procid=myid_y-1,0,-1
           smostlat = cut(1,procid)
           isum = isum + cut(2,procid) - cut(1,procid) + 1
           if (isum >= numbnd) goto 30
        end do
30      if (myid_y /= 0 .and. isum < numbnd .and. smostlat /= 1) then
           call endrun ('SPMDINIT_DYN: Something wrong in computation of southern neighbors')
        end if

!        write(iulog,*)'-----------------------------------------'
!        write(iulog,*)'Number of lats passed north & south = ',numbnd
!        write(iulog,*)'Node  Partition'
!        write(iulog,*)'-----------------------------------------'
!        do procid=0,npes-1
!           write(iulog,200) procid,cut(1,procid),cut(2,procid)
!        end do
!        write(iulog,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!        write(iulog,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      endif

      deallocate (ydist)
      deallocate (zdist)

      return
!
! Formats
!
200   format(i3,4x,i3,'-',i3,7x,i3,'-',i3)

!EOC
   end subroutine spmdinit_dyn

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      call endrun ('decomp_wavenumbers() should never be called in LR dynamics')

   end subroutine decomp_wavenumbers

   subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: placeholder for buffer allocation routine 
! 
! Method: 
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      return

   end subroutine spmdbuf

   subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer, intent(in) :: numperlat            ! number of elements per latitude
!
! Output arguments
!
      integer, intent(out) :: numtot               ! total number of elements (to send or recv)
      integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
      integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
      integer :: p                    ! index
   
      numtot = numperlat*numlats
   
      do p=0,npes-1
         numperproc(p) = numperlat*nlat_p(p)
      end do
     
      displs(:) = 0
      do p=1,npr_y-1
         displs(p) = displs(p-1) + numperproc(p-1)
      end do

      if (npr_z > 1) then
         do p=1,npr_z-1
            displs(p*npr_y:(p+1)*npr_y-1) = displs(0:npr_y-1)
         enddo
      endif

    end subroutine compute_gsfactors

#endif

end module spmd_dyn

