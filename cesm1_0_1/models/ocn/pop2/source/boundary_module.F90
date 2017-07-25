!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module boundary_module !lateral boundary conditions

!BOP
! !MODULE: prognostic
! !DESCRIPTION:
!  This module contains all lateral boundary conditions for the
!  Baltic Sea
!
! !REVISION HISTORY:
!  CVS:$Id: lbc_module.F90,v 1.0 2012/09/17 jjakacki
!  jjakacki@iopan.gda.pl
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use blocks
   use communicate
   use distribution
   use constants
   use domain_size
   use exit_mod
   use prognostic
   use POP_HaloMod
   use POP_FieldMod
   use POP_KindsMod
   use POP_GridHorzMod
   use POP_ErrorMod
   use gather_scatter
   use broadcast
   use io
   use grid
   use time_management

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_lbc, &
             lbc_orlanski_north, &
             lbc_nmn_west

!! PUBLIC DATA
!   logical(log_kind), public :: lbc
!! LOCAL VARIABLES
!! this is beta version, works only with orlanski_north
!   logical(log_kind), public :: orlanski_north,orlanski_south, &
!                        orlanski_east,orlanski_west
!   integer(int_kind),public :: lbc_option !current only orlanski and neumann are implemented
!   integer(int_kind),dimension(nx_block,ny_block,max_blocks_clinic),target :: &
!!========= T & U points are on Arakawa B-grid
!                     orlanski_north_U_MASK, orlanski_north_T_MASK, &
!                     orlanski_south_U_MASK, orlanski_south_T_MASK, &
!                     orlanski_east_U_MASK, orlanski_east_T_MASK, &
!                     orlanski_west_U_MASK, orlanski_west_T_MASK
 contains

!***********************************************************************
!BOP
! !IROUTINE: init_grid1
! !INTERFACE:

 subroutine init_lbc

! !DESCRIPTION:
!  Initializes only grid quantities necessary for running
!  lateral boundary routines
!  
!
! !REVISION HISTORY:
!  same as module
!==============================================
! lbc_option = 1 - orlanski_north
!              2 - orlanski_south
!              3 - orlanski_east
!              4 - orlanski_west
!==============================================
 integer (i4), dimension(:,:), allocatable :: MASK_G
 real(r8), dimension(:,:), allocatable :: MASK_GR 
 integer (int_kind) :: nml_error           ! namelist i/o error flag
 character (char_len) :: north_bfile, south_bfile, &
                         east_bfile, west_bfile
 character (char_len) :: bfile_fmt, lbc_type
 integer (int_kind) :: ii,jj,kk

 namelist /lbc_nml/ orlanski_north, orlanski_south, &
                      orlanski_east, orlanski_west, nmn_north, &
                      nmn_south,nmn_east, nmn_west, &
                      north_bfile, south_bfile, &
                      east_bfile, west_bfile, &
                      bfile_fmt
    !!!!!!
    !!!!!! this should read from pop2_in, but it will be done later

   orlanski_north = .false.
   orlanski_south = .false.
   orlanski_east  = .false.
   orlanski_west  = .false.
   nmn_north = .false.
   nmn_south = .false.
   nmn_east  = .false.
   nmn_west  = .false.
   north_bfile    = 'unknown_north_bfile'
   south_bfile    = 'unknown_south_bfile'
   east_bfile     = 'unknown_east_bfile'
   west_bfile     = 'unknown_west_bfile'
   bfile_fmt      = 'bin'
   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=lbc_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR: reading lbc_nml')
   endif
!   call broadcast_scalar(lbc_type,  master_task)
   call broadcast_scalar(orlanski_north,  master_task)
   call broadcast_scalar(orlanski_south,  master_task)
   call broadcast_scalar(orlanski_east,   master_task)
   call broadcast_scalar(orlanski_west,   master_task)
   call broadcast_scalar(nmn_north,  master_task)
   call broadcast_scalar(nmn_south,  master_task)
   call broadcast_scalar(nmn_east,   master_task)
   call broadcast_scalar(nmn_west,   master_task)
   call broadcast_scalar(north_bfile,     master_task)
   call broadcast_scalar(south_bfile,     master_task)
   call broadcast_scalar(east_bfile,      master_task)
   call broadcast_scalar(west_bfile,      master_task)
   call broadcast_scalar(bfile_fmt,       master_task)
   if (my_task==master_task) then
      write(stdout,*) 'orlanski N/S/E/W = ',orlanski_north, orlanski_south, &
                orlanski_east,orlanski_west
      write(stdout,*) 'nmn N/S/E/W = ',nmn_north, nmn_south, &
                nmn_east,nmn_west
   endif

!       orlanski_north = .true.

   if (orlanski_north) lbc_option = 1
   if (orlanski_south) lbc_option = 2
   if (orlanski_east)  lbc_option = 3
   if (orlanski_west)  lbc_option = 4

! lbc_option = 1
 ORL_RESTORING_MASK=1._r8
 VVEL_RESTORING_MASK=1._r8
 NMN_RESTORING_MASK=1._r8
 NMN_VEL_MASK = 1._r8
 allocate(MASK_G(nx_global,ny_global))
 allocate(MASK_GR(nx_global,ny_global))
! orlanski north boundary
!-- create U-mask for the boundary
 if (orlanski_north) then
     orl_north_first_call = .true.

     ORL_AVE_SPEED_TEMP = 0._r8
     ORL_AVE_SPEED_SALT = 0._r8
     ORL_AVE_SPEED_VVEL = 0._r8
     ORL_AVE_SPEED_VBTP = 0._r8
     ORL_AVE_SPEED_SHGT = 0._r8

     ORL_RESTORING_MASK = 0._r8

     MASK_G(:,:) = 0_i4
     MASK_G(110:160,219) = 1_i4
     call scatter_global(orlanski_north_U_MASK, MASK_G, master_task, distrb_clinic, &
                             field_loc_center, field_type_scalar)
 ! orlanski_north_U_MASK
     MASK_G(:,:) = 0_i4
     MASK_G(110:160,220) = 1_i4
     call scatter_global(orlanski_north_T_MASK, MASK_G, master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)
 ! create  mask 
 ! or read it from the file
    MASK_GR(:,:) = 0._r8
    MASK_GR(80:180,200:222) = 1._r8
    kk = 1
    do jj = 199,180,-1
       MASK_GR(80:180,jj)=1._r8-real(kk,r8)/20._r8
       if (my_task == master_task) then
          write(stdout,*) 'MASK_GR(145,',jj,') = ',MASK_GR(145,jj)
       endif
       kk = kk + 1
    enddo
    where (MASK_GR <= 0._r8)
       MASK_GR = 0._r8
    endwhere
    call scatter_global(VVEL_RESTORING_MASK, MASK_GR, master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)

  endif
    
  if (nmn_west) then


    !---- interior restoring mask
    ! identical like in orlanski
    MASK_GR(:,:) = 0._r8
    MASK_GR(1:75,50:320)=1._r8
    MASK_GR(74:180,190:320)=1._r8
    kk = 1
    do jj = 189,170,-1
       MASK_GR(80:180,jj)=1._r8-real(kk,r8)/20._r8
       if (my_task == master_task) then
          write(stdout,*) 'MASK_GR(145,',jj,') = ',MASK_GR(145,jj)
       endif
       kk = kk + 1
    enddo
    where (MASK_GR <= 0._r8)
       MASK_GR = 0._r8
    endwhere
! ====================================================
!   old version
!     MASK_GR(:,:) = 0._r8
!     MASK_GR(1:60,1:300) = 1._r8
!     MASK_GR(60:80,200:300) = 1._r8
!     kk = 1
!     do ii = 110,81,-1
!        MASK_GR(ii,210:310)=1._r8-real(kk,r8)/30._r8
!        if (my_task == master_task) then
!           write(stdout,*) 'MASK_GR(',ii,'275) = ',MASK_GR(ii,275)
!        endif
!        kk = kk + 1
!     enddo
!    
!     where (MASK_GR <= 0._r8)
!        MASK_GR = 0._r8
!     endwhere
! ====================================================
     call scatter_global(NMN_RESTORING_MASK, MASK_GR, master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)
     if (my_task == master_task) then
        open(88,file='nmn_restoring_mask.bin', access='direct',form = 'unformatted', &
                     recl=600*640*8, status='unknown');
        write(88,rec=1) MASK_GR
        close(88)
     endif
     MASK_GR(:,:) = 1._r8
     MASK_GR(1:30,1:300) = 1._r8
     kk = 1
     do ii = 60,1,-1
!jj        MASK_GR(ii,1:300)=1._r8-real(kk,r8)/65._r8
!jj        MASK_GR(ii,1:300)=1._r8-real(kk,r8)/85._r8
!jj for tides 03 and 04        MASK_GR(ii,1:300)=1._r8-real(kk,r8)/75._r8
        MASK_GR(ii,1:300)=1._r8-real(kk,r8)/95._r8
        if (my_task == master_task) then
           write(stdout,*) 'MASK_GR(',ii,'150) = ',MASK_GR(ii,150)
        endif
        kk = kk + 1
     enddo
     call scatter_global(NMN_VEL_MASK, MASK_GR, master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)
     NMN_VEL_MASK = 1._r8

     MASK_G(:,:) = 0_i4
     MASK_G(3,:) = 1_i4
     call scatter_global(nmn_MASK, MASK_G, master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)
  endif
 deallocate(MASK_G)
 deallocate(MASK_GR)
 !only temporary
! ORL_RESTORING_MASK=1._r8
! VVEL_RESTORING_MASK=1._r8
! NMN_RESTORING_MASK=1._r8
! ORL_RESTORING_MASK = NMN_RESTORING_MASK
 end subroutine init_lbc

 subroutine lbc_orlanski_north(k,field_type,prognostic_field,iblock)


 integer(i4), intent(in) :: k,iblock  !level number and local block id 
 character(len=6), intent(in) :: field_type
 character(4), intent(in) :: prognostic_field
! field_type = 'vector' or 'scalar' 
! vector for baroclinic and barotropic velocities
! scalar for tracers and sea level
! 
! equation is the same for all prognostic variables
! routine will use only different mask
!
! routine should be called for each level separately
! local variables

 real(r8) :: R(nx_block,ny_block,max_blocks_clinic) ! = c*delta_t/delta_x

! arrays balow are temporary and will work for velocities,
! tracers and surface height
 
 real(r8),dimension(:,:,:),allocatable :: &
         NP1_I,N_I,N_IM1,NP1_IM1,N_IM2
 integer(i4),dimension(:,:,:),allocatable :: &
         current_mask
 integer(i4) :: prognostic_field_option, errorCode, ii,jj,kk
 real(r8), parameter :: & 
!      cfl_limit = 0.5    !tides02_02
   cfl_limit = 1.         !tides02_03

 !!=====================================================
 !! prognostic_field =
 !!        'temp'
 !!        'salt'
 !!        'vvel'
 !!        'vbtp'
 !!        'shgt'
 !!=====================================================
!jj bid = this_block%local_id
 allocate(current_mask(nx_block,ny_block,max_blocks_clinic))
 current_mask = 0_i4
 prognostic_field_option = 0_i4
!xx write(*,*) 'lbc north 001'
 select case(trim(prognostic_field))
   case('temp')
       prognostic_field_option = 1_i4
       current_mask = orlanski_north_T_MASK
   case('salt')
       prognostic_field_option = 2_i4
       current_mask = orlanski_north_T_MASK
   case('uvel') 
       prognostic_field_option = 3_i4
       current_mask = orlanski_north_U_MASK
   case('vvel') 
       prognostic_field_option = 4_i4
       current_mask = orlanski_north_U_MASK
   case('vbtp') 
       prognostic_field_option = 5_i4
       current_mask = orlanski_north_U_MASK
   case('shgt')
       prognostic_field_option = 6_i4
       current_mask = orlanski_north_T_MASK
   case default
     call exit_POP(sigAbort, &
       'lbc_north: ERROR: boundary/prognostic_field_option')
   end select
!xx   write(*,*) 'lbc north 002'
! if (trim(prognostic_field) == 'temp') prognostic_field_option = 1_i4
! if (trim(prognostic_field) == 'salt') prognostic_field_option = 2_i4
! if (trim(prognostic_field) == 'uvel') prognostic_field_option = 3_i4
! if (trim(prognostic_field) == 'vvel') prognostic_field_option = 4_i4
! if (trim(prognostic_field) == 'vbtp') prognostic_field_option = 5_i4
! if (trim(prognostic_field) == 'shgt') prognostic_field_option = 6_i4
! if (prognostic_field_option == 0) &
!     call exit_POP(sigAbort,'ERROR: boundary/prognostic_field_option') 
! if (k == 0 & prognostic_field_option /= 5)  &
!     call exit_POP(sigAbort,'ERROR: boundary/prognostic_field_option/shgt')
! first test is for baroclinic velocity, barotropic should be calculated based on 
! baroclinic, tracers and sea level should be changed based on velocity in the next step
! this is the first interpolation, will see what will happen 
! k = 1
! iblock = 1
 !!!  ===  initialize
 allocate(NP1_I(nx_block,ny_block,max_blocks_clinic), &
            N_I(nx_block,ny_block,max_blocks_clinic), &
          N_IM1(nx_block,ny_block,max_blocks_clinic), &
        NP1_IM1(nx_block,ny_block,max_blocks_clinic), &
          N_IM2(nx_block,ny_block,max_blocks_clinic))
! allocate(current_mask(nx_block,ny_block,max_blocks_clinic))
! current_mask = 0_i4
!xx write(*,*) 'lbc north 003'
 !if     (field_type == 'scalar' .and. prognostic_field_option == 1) then
 !   current_mask = orlanski_north_T_MASK
 !elseif (field_type == 'scalar' .and. prognostic_field_option == 2) then
 !   current_mask = orlanski_north_T_MASK
 !elseif (field_type == 'vector' .and. prognostic_field_option == 4) then
 !   current_mask = orlanski_north_U_MASK
 !elseif (field_type == 'vector' .and. prognostic_field_option == 5) then
 !   current_mask = orlanski_north_U_MASK
 !elseif (field_type == 'scalar' .and. prognostic_field_option == 6) then
 !   current_mask = orlanski_north_T_MASK
 !else
 !   current_mask = 0_i4 ! lbc will work but with no effect
 !endif
!  NP1_I  ---->>> X(N+1,I) - target 
if (orlanski_north) then 
!
!   real (r8), dimension(nx_block,ny_block,km,nt,3,max_blocks_clinic), &
!      target :: &
!      TRACER  
!xx   write(*,*) 'lbc north 004'
   select case (prognostic_field_option)
   case (1)     !***** TEMP *****
!xx      write(*,*) 'lbc north 005'
      N_I(:,:,iblock)    = TRACER(:,:,k,1,curtime,iblock)                               ! X(N,I)
      N_IM1(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,1,curtime,iblock),SHIFT = -1, DIM = 2)   ! X(N,I-1)
      NP1_IM1(:,:,iblock)= EOSHIFT(TRACER(:,:,k,1,newtime,iblock),SHIFT = -1, DIM = 2)   ! X(N+1,I-1)
      N_IM2(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,1,curtime,iblock),SHIFT = -2, DIM = 2)   ! X(N,I-2)
!      call POP_HaloUpdate(N_IM1,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(NP1_IM1,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(N_IM2,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
   case (2)     !***** SALT *****
!xx      write(*,*) 'lbc north 006'
      N_I(:,:,iblock)    = TRACER(:,:,k,2,curtime,iblock)                               ! X(N,I)
      N_IM1(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,2,curtime,iblock),SHIFT = -1, DIM = 2)   ! X(N,I-1)
      NP1_IM1(:,:,iblock)= EOSHIFT(TRACER(:,:,k,2,newtime,iblock),SHIFT = -1, DIM = 2)   ! X(N+1,I-1)
      N_IM2(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,2,curtime,iblock),SHIFT = -2, DIM = 2)   ! X(N,I-2)
!      call POP_HaloUpdate(N_IM1,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(NP1_IM1,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(N_IM2,POP_haloClinic, POP_gridHorzLocCenter, &
!                            POP_fieldKindScalar, errorCode,fillValue = 0.0_POP_r8)
   case (4)     !***** VVEL *****
!xx      write(*,*) 'lbc north 007'
      N_I(:,:,iblock)    = VVEL(:,:,k,curtime,iblock)                               ! X(N,I)
      N_IM1(:,:,iblock)  = EOSHIFT(VVEL(:,:,k,curtime,iblock),SHIFT = -1, DIM = 2)   ! X(N,I-1)
      NP1_IM1(:,:,iblock)= EOSHIFT(VVEL(:,:,k,newtime,iblock),SHIFT = -1, DIM = 2)   ! X(N+1,I-1)
      N_IM2(:,:,iblock)  = EOSHIFT(VVEL(:,:,k,curtime,iblock),SHIFT = -2, DIM = 2)   ! X(N,I-2)
!      call POP_HaloUpdate(N_IM1,POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(NP1_IM1,POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(N_IM2,POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8)
   case (5)
!xx      write(*,*) 'lbc north 008'
      N_I(:,:,iblock)    = VBTROP(:,:,curtime,iblock)                               ! X(N,I)
      N_IM1(:,:,iblock)  = EOSHIFT(VBTROP(:,:,curtime,iblock),SHIFT = -1, DIM = 2)   ! X(N,I-1)
      NP1_IM1(:,:,iblock)= EOSHIFT(VBTROP(:,:,newtime,iblock),SHIFT = -1, DIM = 2)   ! X(N+1,I-1)
      N_IM2(:,:,iblock)  = EOSHIFT(VBTROP(:,:,curtime,iblock),SHIFT = -2, DIM = 2)   ! X(N,I-2)
!      call POP_HaloUpdate(N_IM1, POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(NP1_IM1, POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(N_IM2, POP_haloClinic,POP_gridHorzLocNECorner, &
!                          POP_fieldKindVector, errorCode, fillValue = 0.0_POP_r8) 
   case (6)     !***** PSURF *****
!xx      write(*,*) 'lbc north 009'
      N_I(:,:,iblock)    = PSURF(:,:,curtime,iblock)                               ! X(N,I)
      N_IM1(:,:,iblock)  = EOSHIFT(PSURF(:,:,curtime,iblock),SHIFT = -1, DIM = 2)   ! X(N,I-1)
      NP1_IM1(:,:,iblock)= EOSHIFT(PSURF(:,:,newtime,iblock),SHIFT = -1, DIM = 2)   ! X(N+1,I-1)
      N_IM2(:,:,iblock)  = EOSHIFT(PSURF(:,:,curtime,iblock),SHIFT = -2, DIM = 2)   ! X(N,I-2)
!      call POP_HaloUpdate(N_IM1,POP_haloClinic, POP_gridHorzLocCenter,  &
!                          POP_fieldKindScalar, errorCode, fillValue = 0.0_POP_r8) 
!      call POP_HaloUpdate(NP1_IM1,POP_haloClinic, POP_gridHorzLocCenter,  &
!                          POP_fieldKindScalar, errorCode, fillValue = 0.0_POP_r8)
!      call POP_HaloUpdate(N_IM2,POP_haloClinic, POP_gridHorzLocCenter,  &
!                          POP_fieldKindScalar, errorCode, fillValue = 0.0_POP_r8)
   case default
     call exit_POP(sigAbort, &
       'lbc_north: invalid prognostic_field_option')
   end select
   R(:,:,iblock) = 0._r8
   R(:,:,iblock) = N_IM1(:,:,iblock) - N_IM2(:,:,iblock)
   where (R(:,:,iblock) /= 0._r8)
      R(:,:,iblock) = (NP1_IM1(:,:,iblock)-N_IM1(:,:,iblock))/ &
                      (N_IM1(:,:,iblock) - N_IM2(:,:,iblock))
   end where
   select case (prognostic_field_option)
   case (1)       !*** TEMP ***
      where(R(:,:,iblock) >= 0._r8)
         R(:,:,iblock) = 0._r8
      endwhere
      where( R(:,:,iblock) <= -cfl_limit)
         R(:,:,iblock) = -cfl_limit
      endwhere
      ORL_AVE_SPEED_TEMP(:,:,k,iblock,2:4) = &
               ORL_AVE_SPEED_TEMP(:,:,k,iblock,1:3)
      ORL_AVE_SPEED_TEMP(:,:,k,iblock,1) = R(:,:,iblock)
      ORL_AVE_SPEED_TEMP(:,:,k,iblock,5) = &
             ( ORL_AVE_SPEED_TEMP(:,:,k,iblock,1) + &
               ORL_AVE_SPEED_TEMP(:,:,k,iblock,2) + &
               ORL_AVE_SPEED_TEMP(:,:,k,iblock,3) + &
               ORL_AVE_SPEED_TEMP(:,:,k,iblock,4) )/4._r8
      R(:,:,iblock) = ORL_AVE_SPEED_TEMP(:,:,k,iblock,5)
      NP1_I = N_I + R*(N_I-N_IM1)
      where(current_mask(:,:,iblock) == 1_i4)
            TRACER(:,:,k,1,newtime,iblock) = (NP1_I(:,:,iblock)+ &
                   NP1_IM1(:,:,iblock))/2
      endwhere
   case (2)       !*** SALT ***
     where(R(:,:,iblock) >= 0._r8)
         R(:,:,iblock) = 0._r8
      endwhere
      where( R(:,:,iblock) <= -cfl_limit)
         R(:,:,iblock) = -cfl_limit
      endwhere
      ORL_AVE_SPEED_SALT(:,:,k,iblock,2:4) = &
               ORL_AVE_SPEED_SALT(:,:,k,iblock,1:3)
      ORL_AVE_SPEED_SALT(:,:,k,iblock,1) = R(:,:,iblock)
      ORL_AVE_SPEED_SALT(:,:,k,iblock,5) = &
             ( ORL_AVE_SPEED_SALT(:,:,k,iblock,1) + &
               ORL_AVE_SPEED_SALT(:,:,k,iblock,2) + &
               ORL_AVE_SPEED_SALT(:,:,k,iblock,3) + &
               ORL_AVE_SPEED_SALT(:,:,k,iblock,4) )/4._r8
      R(:,:,iblock) = ORL_AVE_SPEED_SALT(:,:,k,iblock,5)

      NP1_I = N_I + R*(N_I-N_IM1)
      where(current_mask(:,:,iblock) == 1_i4)
            TRACER(:,:,k,2,newtime,iblock) = (NP1_I(:,:,iblock)+ &
                   NP1_IM1(:,:,iblock))/2
      endwhere
   case (4)       !*** VVEL ***
!      ORL_RESTORING_MASK(:,:,k,iblock) = 0._r8
      where(R(:,:,iblock) >= 0._r8)
         R(:,:,iblock) = 0._r8
      endwhere
      where( R(:,:,iblock) <= -cfl_limit)
         R(:,:,iblock) = -cfl_limit
      endwhere
         ORL_AVE_SPEED_VVEL(:,:,k,iblock,2:4) = &
               ORL_AVE_SPEED_VVEL(:,:,k,iblock,1:3)
         ORL_AVE_SPEED_VVEL(:,:,k,iblock,1) = R(:,:,iblock)
      ORL_AVE_SPEED_VVEL(:,:,k,iblock,5) = &
             ( ORL_AVE_SPEED_VVEL(:,:,k,iblock,1) + &
               ORL_AVE_SPEED_VVEL(:,:,k,iblock,2) + &
               ORL_AVE_SPEED_VVEL(:,:,k,iblock,3) + &
               ORL_AVE_SPEED_VVEL(:,:,k,iblock,4) )/4._r8
      R(:,:,iblock) = ORL_AVE_SPEED_VVEL(:,:,k,iblock,5)
      NP1_I = N_I + R*(N_I-N_IM1)
      where(current_mask(:,:,iblock) == 1_i4)
!         VVEL(:,:,k,newtime,iblock) = (NP1_I(:,:,iblock)+ &
!                   NP1_IM1(:,:,iblock))/2
         VVEL(:,:,k,newtime,iblock) = NP1_I(:,:,iblock)
      endwhere
!      where (VVEL(:,:,k,newtime,iblock) < 0._r8)
!         ORL_RESTORING_MASK(:,:,k,iblock) = 1._r8
!      endwhere
!      ORL_RESTORING_MASK(:,:,k,iblock) = &
!              VVEL_RESTORING_MASK(:,:,iblock)
   case (5)       !*** VBTP ***
      where(R(:,:,iblock) >= 0._r8)
         R(:,:,iblock) = 0._r8
      endwhere
      where( R(:,:,iblock) <= -cfl_limit)
         R(:,:,iblock) = -cfl_limit
      endwhere
         ORL_AVE_SPEED_VBTP(:,:,iblock,2:4) = &
               ORL_AVE_SPEED_VBTP(:,:,iblock,1:3)
         ORL_AVE_SPEED_VBTP(:,:,iblock,1) = R(:,:,iblock)
      ORL_AVE_SPEED_VBTP(:,:,iblock,5) = &
             ( ORL_AVE_SPEED_VBTP(:,:,iblock,1) + &
               ORL_AVE_SPEED_VBTP(:,:,iblock,2) + &
               ORL_AVE_SPEED_VBTP(:,:,iblock,3) + &
               ORL_AVE_SPEED_VBTP(:,:,iblock,4) )/4._r8
      R(:,:,iblock) = ORL_AVE_SPEED_VBTP(:,:,iblock,5)
      NP1_I = N_I + R*(N_I-N_IM1)
      where(current_mask(:,:,iblock) == 1_i4)
!         VBTROP(:,:,newtime,iblock) = (NP1_I(:,:,iblock) + &
!                   NP1_IM1(:,:,iblock))/2
         VBTROP(:,:,newtime,iblock) = NP1_I(:,:,iblock)
      endwhere
!      ORL_RESTORING_MASK(:,:,1,iblock) = 0._r8
!      where(VBTROP(:,:,newtime,iblock) <=0._r8 .and. &
!                 VVEL_RESTORING_MASK(:,:,iblock) > 0._r8)
!         ORL_RESTORING_MASK(:,:,1,iblock) = 1._r8
!      endwhere
     do ii = 1,km
        ORL_RESTORING_MASK(:,:,ii,iblock) = VVEL_RESTORING_MASK(:,:,iblock) !0._r8  !no restoring at the boundary
     enddo
   case (6)
      where(R(:,:,iblock) >= 0._r8)
         R(:,:,iblock) = 0._r8
      endwhere
      where( R(:,:,iblock) <= -cfl_limit)
         R(:,:,iblock) = -cfl_limit
      endwhere
         ORL_AVE_SPEED_SHGT(:,:,iblock,2:4) = &
               ORL_AVE_SPEED_SHGT(:,:,iblock,1:3)
         ORL_AVE_SPEED_SHGT(:,:,iblock,1) = R(:,:,iblock)
      ORL_AVE_SPEED_SHGT(:,:,iblock,5) = &
             ( ORL_AVE_SPEED_SHGT(:,:,iblock,1) + &
               ORL_AVE_SPEED_SHGT(:,:,iblock,2) + &
               ORL_AVE_SPEED_SHGT(:,:,iblock,3) + &
               ORL_AVE_SPEED_SHGT(:,:,iblock,4) )/4._r8
      R(:,:,iblock) = ORL_AVE_SPEED_SHGT(:,:,iblock,5)

      NP1_I = N_I + R*(N_I-N_IM1)
      where(current_mask(:,:,iblock) == 1_i4)
!         PSURF(:,:,newtime,iblock) = (NP1_I(:,:,iblock) + &
!                   NP1_IM1(:,:,iblock))/2
          PSURF(:,:,newtime,iblock) = NP1_IM1(:,:,iblock) ! zero gradient
      endwhere
   case default
     call exit_POP(sigAbort, &
       'lbc_north: invalid prognostic_field_option')
   end select
endif
 deallocate(NP1_I,N_I,N_IM1,NP1_IM1,N_IM2) 
 end subroutine lbc_orlanski_north

 subroutine lbc_nmn_west(k,prognostic_field,iblock)
  
! in the west direction all can be calculated in identical way
! so field type is not important

 integer(i4), intent(in) :: k,iblock  !level number and local block id 
 character(4), intent(in) :: prognostic_field
 real(r8),dimension(:,:,:),allocatable :: &
         N_IM1

 allocate(N_IM1(nx_block,ny_block,max_blocks_clinic))

 select case(trim(prognostic_field))
   case('temp')

!      where(current_mask(:,:,iblock) == 1_i4)
!            TRACER(:,:,k,1,newtime,iblock) = (NP1_I(:,:,iblock)+ &
!                   NP1_IM1(:,:,iblock))/2
!      endwhere
       ! e = eoshift(b, SHIFT=1, DIM=1) 
       N_IM1(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,1,newtime,iblock),SHIFT = 1, DIM = 1)
       where (nmn_MASK(:,:,iblock) == 1_int_kind)
          TRACER(:,:,k,1,newtime,iblock) = N_IM1(:,:,iblock)
       endwhere
   case('salt')
       N_IM1(:,:,iblock)  = EOSHIFT(TRACER(:,:,k,2,newtime,iblock),SHIFT = 1, DIM = 1)
       where (nmn_MASK(:,:,iblock) == 1) 
          TRACER(:,:,k,2,newtime,iblock) = N_IM1(:,:,iblock)
       endwhere
   case('uvel') 
       N_IM1(:,:,iblock)  = EOSHIFT(UVEL(:,:,k,newtime,iblock),SHIFT = 1, DIM = 1)
       where (nmn_MASK(:,:,iblock) == 1)
           UVEL(:,:,k,newtime,iblock) = N_IM1(:,:,iblock)
       endwhere
       UVEL(:,:,k,newtime,iblock) = UVEL(:,:,k,newtime,iblock)*NMN_VEL_MASK(:,:,iblock)
   case('vvel')
!jj       N_IM1(:,:,iblock)  = EOSHIFT(VVEL(:,:,k,newtime,iblock),SHIFT = 1, DIM = 1)
!jj       where (nmn_MASK(:,:,iblock) == 1)
!jj           VVEL(:,:,k,newtime,iblock) = N_IM1(:,:,iblock)
!jj       endwhere
       VVEL(:,:,k,newtime,iblock) = VVEL(:,:,k,newtime,iblock)*NMN_VEL_MASK(:,:,iblock)
   case('vbtp') 
!jj       N_IM1(:,:,iblock)  = EOSHIFT(VBTROP(:,:,newtime,iblock),SHIFT = 1, DIM = 1)
!jj       where (nmn_MASK(:,:,iblock) == 1)
!jj            VBTROP(:,:,newtime,iblock) = N_IM1(:,:,iblock)
!jj       endwhere
       VBTROP(:,:,newtime,iblock)=VBTROP(:,:,newtime,iblock)*NMN_VEL_MASK(:,:,iblock)
   case('ubtp')
       N_IM1(:,:,iblock)  = EOSHIFT(UBTROP(:,:,newtime,iblock),SHIFT = 1, DIM = 1)
       where (nmn_MASK(:,:,iblock) == 1)
            UBTROP(:,:,newtime,iblock) = N_IM1(:,:,iblock)
       endwhere
       UBTROP(:,:,newtime,iblock) = UBTROP(:,:,newtime,iblock)*NMN_VEL_MASK(:,:,iblock)
   case('shgt')
       N_IM1(:,:,iblock)  = EOSHIFT(PSURF(:,:,newtime,iblock),SHIFT = 1, DIM = 1)
       where (nmn_MASK(:,:,iblock) == 1)
            PSURF(:,:,newtime,iblock) = N_IM1(:,:,iblock)
       endwhere
!       PSURF(:,:,newtime,iblock) = PSURF(:,:,newtime,iblock)*NMN_VEL_MASK(:,:,iblock)
   case default
     call exit_POP(sigAbort, &
       'lbc_north: ERROR: boundary/prognostic_field_option')
   end select
 deallocate(N_IM1) 
 end subroutine lbc_nmn_west
   
 end module boundary_module




