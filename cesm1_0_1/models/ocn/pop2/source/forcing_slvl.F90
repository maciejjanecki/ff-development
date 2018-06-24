!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_slvl

!BOP
! !MODULE: forcing_slvl
!
! !DESCRIPTION:
!  Contains routines and variables used for determining the
!  surface atmopheric pressure for using atmospheric pressure
!  forcing at the surface.
!
! !REVISION HISTORY:
!  SVN:$Id: forcing_slvl.F90 12674 2008-10-31 22:21:32Z njn01 $
!  Jaromir Jakacki @ 25.05.2018
!
! !USES:

   use kinds_mod
   use domain
   use constants
   use broadcast
   use io
   use forcing_tools
   use time_management
   use grid
   use exit_mod
   use prognostic, only: slvl_assimilation
   use gather_scatter
   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public:: init_slvl,     &
            set_slvl

! !PUBLIC DATA MEMBERS:

   real (r8), public :: & ! public for use in restart
      slvl_interp_last      ! time when last interpolation was done

   character (char_len), public :: &! needed in barotropic
      slvl_data_type     !  keyword for the period of the forcing data

   real (r8), dimension(:,:,:), allocatable, public :: &
          SLVL_MASK,SLVL_OUT_MASK

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:,:), allocatable :: &
      SEA_LEVEL_DATA !  external sea level data at multiple times that
                     !  may be interpolated in time to get SEA_LEVEL

   real (r8), dimension(12) :: &
      slvl_data_time  !  time (in hours) corresponding to sea level
                    !  in SEA_LEVEL_DATA

   real (r8), dimension(20) :: &
      slvl_data_renorm   ! scale factors for changing input units
                       ! to model units
   logical (log_kind) :: use_slvl_mask

   real (r8) ::       &
      slvl_data_inc,    &! time increment between values of forcing data
      slvl_data_next,   &! time to be used for next forcing data
      slvl_data_update, &! time when new forcing value added to interp set
      slvl_interp_inc,  &! time increment between interpolation
      slvl_interp_next   ! time when next interpolation will be done

   integer (int_kind) ::   &
      slvl_interp_order,     &! order of temporal interpolation
      slvl_data_time_min_loc  ! index of 3rd dimension of SEA_LEVEL_DATA,
                            !   containing the minimum forcing time

   character (char_len) :: &
      slvl_filename,    &!  name of file conainting forcing data
      slvl_file_fmt,    &!  data format of forcing file (bin or nc)
      slvl_interp_freq, &!  keyword for period of temporal interpolation
      slvl_interp_type, &
      slvl_data_label,  &
      slvl_formulation, &
      slvl_mask_file, &
      slvl_mask_file_fmt

   character (char_len), dimension(:), allocatable :: &
      slvl_data_names    ! field names for getting ap data

   integer (int_kind), dimension(:), allocatable :: &
      slvl_bndy_loc,    &! location and field type for ghost cell
      slvl_bndy_type     !   updates

   type (datafile) :: &
      slvl_data_file    ! io file type for ap data file

   type (io_field_desc) :: &
      slvl_data_in      ! io field type for reading ap data

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_slvl
! !INTERFACE:

 subroutine init_slvl(SEA_LEVEL)

! !DESCRIPTION:
!  Initializes sea level assimilation by either calculating or reading
!  in the sea level.  Also does initial book-keeping concerning
!  when new data is needed for the temporal interpolation and
!  when the forcing will need to be updated.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      SEA_LEVEL        !  external sea level data at current timestep

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::     &
      n, iblock, nu,         &! dummy loop index and unit number
      nml_error               ! namelist i/o error flag

   character (char_len) ::   &
      forcing_filename        ! name of file containing forcing data

   type (datafile) :: &
      slvl_data_file  ! data file descriptor for interior pot temp data

   type (io_field_desc) :: &
      slvl_data_in          ! io field descriptor for input pot temp data

   type (block) :: &
      this_block    ! block information for current block

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptor for horiz dims
      month_dim      ! dimension descriptor for monthly fields

   real (r8), dimension(:,:),allocatable :: &
      WORK

   real (r8), dimension(:,:,:,:), allocatable :: &
      TEMP_DATA   ! temp array for reading monthly data

   namelist /forcing_slvl_nml/ slvl_data_type,   slvl_data_inc,    &
                             slvl_interp_type, slvl_interp_freq, &
                             slvl_interp_inc,  slvl_filename,    &
                             slvl_data_renorm, slvl_file_fmt,    &
                             use_slvl_mask, slvl_mask_file,slvl_mask_file_fmt

!-----------------------------------------------------------------------
!
!  read atmospheric pressure namelist input after setting
!  default values.
!
!-----------------------------------------------------------------------

   slvl_data_type   = 'none'
   slvl_data_inc    = 1.e20_r8
   slvl_interp_type = 'nearest'
   slvl_interp_freq = 'never'
   slvl_interp_inc  = 1.e20_r8
   slvl_filename    = 'unknown-ap'
   slvl_file_fmt    = 'bin'
   slvl_data_renorm = c1
!  slvl_data_renorm = c10  ! convert from Pa to dynes/cm**2

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
         read(nml_in, nml=forcing_slvl_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading forcing_slvl namelist')
   endif

   call broadcast_scalar(slvl_data_type,     master_task)
   call broadcast_scalar(slvl_data_inc,      master_task)
   call broadcast_scalar(slvl_interp_type,   master_task)
   call broadcast_scalar(slvl_interp_freq,   master_task)
   call broadcast_scalar(slvl_interp_inc,    master_task)
   call broadcast_scalar(slvl_filename,      master_task)
   call broadcast_scalar(slvl_file_fmt,      master_task)
   call broadcast_array (slvl_data_renorm,   master_task)
   call broadcast_scalar(use_slvl_mask,      master_task)
   call broadcast_scalar(slvl_mask_file_fmt, master_task)
   call broadcast_scalar(slvl_mask_file,     master_task)




   if (use_slvl_mask .and. slvl_data_type .ne. 'none') then
      slvl_assimilation = 1
   endif

   slvl_formulation = char_blank

!-----------------------------------------------------------------------
!
!  convert data_type to 'monthly-calendar' if input is 'monthly'
!
!-----------------------------------------------------------------------

   if (slvl_data_type == 'monthly') slvl_data_type = 'monthly-calendar'

!-----------------------------------------------------------------------
!
!  convert interp_type to corresponding integer value.
!
!-----------------------------------------------------------------------

   select case (slvl_interp_type)

   case ('nearest')
      slvl_interp_order = 1

   case ('linear')
      slvl_interp_order = 2

   case ('4point')
      slvl_interp_order = 4

   case default
      call exit_POP(sigAbort, &
                    'init_slvl: Unknown value for slvl_interp_type')

   end select

!-----------------------------------------------------------------------
!
!  set values of atm pressure arrays (SEA_LEVEL or SEA_LEVEL_DATA)
!  depending on the type of the atm pressure data.
!
!-----------------------------------------------------------------------

   select case (slvl_data_type)

   case ('none')

      !***  no atm pressure, therefore no interpolation in time needed
      !***  (slvl_interp_freq = 'none'), nor are there any new values
      !***  to be used (slvl_data_next = slvl_data_update = never).

      SEA_LEVEL = c0
      slvl_data_next = never
      slvl_data_update = never
      slvl_interp_freq = 'never'

   case ('analytic')

      !*** simple analytic atm pressure that is constant in time,
      !*** therefore no interpolation in time is needed
      !*** (slvl_interp_freq = 'none'), nor are there any new values
      !*** to be used (slvl_data_next = slvl_data_update = never).

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,WORK)

      do iblock=1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         WORK = (ULAT(:,:,iblock)*radian -  25.0_r8)**2 + &
                (ULON(:,:,iblock)*radian - 180.0_r8)**2

         where (CALCT(:,:,iblock) .and. WORK < 100.0_r8**2)
            SEA_LEVEL(:,:,iblock) = grav*exp(-WORK*p5/15.0_r8**2)
            SEA_LEVEL(:,:,iblock) = WORK - &
                   sum(WORK*TAREA(:,:,iblock),mask=CALCT(:,:,iblock))/ &
                   sum(     TAREA(:,:,iblock),mask=CALCT(:,:,iblock))
         elsewhere
            SEA_LEVEL(:,:,iblock) = c0
         end where

         if (slvl_data_renorm(1) /= c1) SEA_LEVEL(:,:,iblock) = &
                                      SEA_LEVEL(:,:,iblock)*  &
                                      slvl_data_renorm(1)

      end do
      !$OMP END PARALLEL DO

      slvl_data_next = never
      slvl_data_update = never
      slvl_interp_freq = 'never'

   case ('annual')

      !*** annual mean climatological atm pressure (read in from a file)
      !*** constant in time, therefore no interpolation in time needed
      !*** (slvl_interp_freq = 'none'), nor are there any new values
      !***  to be used (slvl_data_next = slvl_data_update = never).

      allocate(slvl_data_names(1), &
               slvl_bndy_loc  (1), &
               slvl_bndy_type (1))
      slvl_data_names(1) = 'EXTERNAL SEA LEVEL'
      slvl_bndy_loc  (1) = field_loc_center
      slvl_bndy_type (1) = field_type_scalar

      slvl_data_file = construct_file(slvl_file_fmt,                 &
                                    full_name=trim(slvl_filename), &
                                    record_length=rec_type_dbl,  &
                                    recl_words=nx_global*ny_global)

      call data_set(slvl_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      slvl_data_in = construct_io_field(trim(slvl_data_names(1)),  &
                     dim1=i_dim, dim2=j_dim,                   &
                     field_loc = slvl_bndy_loc(1),               &
                     field_type = slvl_bndy_type(1),             &
                     d2d_array = SEA_LEVEL)

      call data_set (slvl_data_file, 'define', slvl_data_in)
      call data_set (slvl_data_file, 'read',   slvl_data_in)
      call data_set (slvl_data_file, 'close')
      call destroy_io_field(slvl_data_in)
      call destroy_file(slvl_data_file)

      if (slvl_data_renorm(1) /= c1) SEA_LEVEL = &
                                   SEA_LEVEL*slvl_data_renorm(1)

      slvl_data_next = never
      slvl_data_update = never
      slvl_interp_freq = 'never'

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a22,a)') ' Sea Level Annual file read: ', &
                                 trim(forcing_filename)
      endif

   case ('monthly-equal','monthly-calendar')

      !*** monthly mean climatological atm pressure. all 12 months
      !*** are read in from a file. interpolation order
      !*** (slvl_interp_order) may be specified with namelist input.

      allocate(SEA_LEVEL_DATA(nx_block,ny_block,max_blocks_clinic, &
                                                          1,0:12), &
               TEMP_DATA(nx_block,ny_block,12,max_blocks_clinic))

      allocate(slvl_data_names(1), &
               slvl_bndy_loc  (1), &
               slvl_bndy_type (1))

      SEA_LEVEL_DATA = c0
      slvl_data_names(1) = 'EXTERNAL SEA LEVEL'
      slvl_bndy_loc  (1) = field_loc_center
      slvl_bndy_type (1) = field_type_scalar

      call find_forcing_times(slvl_data_time, slvl_data_inc,            &
                              slvl_interp_type, slvl_data_next,         &
                              slvl_data_time_min_loc, slvl_data_update, &
                              slvl_data_type)

      slvl_data_file = construct_file(slvl_file_fmt,                 &
                                    full_name=trim(slvl_filename), &
                                    record_length=rec_type_dbl,  &
                                    recl_words=nx_global*ny_global)

      call data_set(slvl_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)
      month_dim = construct_io_dim('month',12)

      slvl_data_in = construct_io_field(trim(slvl_data_names(1)),  &
                     dim1=i_dim, dim2=j_dim, dim3=month_dim,   &
                     field_loc = slvl_bndy_loc(1),               &
                     field_type = slvl_bndy_type(1),             &
                     d3d_array = TEMP_DATA)

      call data_set (slvl_data_file, 'define', slvl_data_in)
      call data_set (slvl_data_file, 'read',   slvl_data_in)
      call destroy_io_field(slvl_data_in)

      call data_set (slvl_data_file, 'close')
      call destroy_file(slvl_data_file)

      !$OMP PARALLEL DO PRIVATE(iblock,n)
      do iblock=1,nblocks_clinic
      do n=1,12
         SEA_LEVEL_DATA(:,:,iblock,1,n) = TEMP_DATA(:,:,n,iblock)

         if (slvl_data_renorm(1) /= c1) &
            SEA_LEVEL_DATA(:,:,iblock,1,n) = &
            SEA_LEVEL_DATA(:,:,iblock,1,n)*slvl_data_renorm(1)
      end do
      end do
      !$OMP END PARALLEL DO

      deallocate(TEMP_DATA)

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,'(a23,a)') ' Monthly Sea Level file read: ', &
                                 trim(forcing_filename)
      endif

   case ('n-hour')

      !*** atm pressure specified every n-hours, where the n-hour
      !*** increment should be specified with namelist input
      !*** (slvl_data_inc).  Only as many times as are necessary based
      !*** on the order of the temporal interpolation scheme
      !*** (slvl_interp_order) reside in memory at any given time.

      allocate (SEA_LEVEL_DATA(nx_block,ny_block,max_blocks_clinic, &
                                              1,0:slvl_interp_order))
      allocate(slvl_data_names(1), &
               slvl_bndy_loc  (1), &
               slvl_bndy_type (1))

      SEA_LEVEL_DATA = c0
      slvl_data_names(1) = 'EXTERNAL SEA LEVEL'
      slvl_bndy_loc  (1) = field_loc_center
      slvl_bndy_type (1) = field_type_scalar
!      write(*,*) 'SEA_LEVEL_DATA/slvl_data_names(1)/slvl_bndy_loc  (1)/slvl_bndy_type (1) = ',&
!              slvl_data_names(1),slvl_bndy_loc  (1),slvl_bndy_type (1)

      call find_forcing_times(slvl_data_time,         slvl_data_inc,    &
                              slvl_interp_type,       slvl_data_next,   &
                              slvl_data_time_min_loc, slvl_data_update, &
                              slvl_data_type)

      do n = 1, slvl_interp_order
         call get_forcing_filename(forcing_filename, slvl_filename, &
                                   slvl_data_time(n), slvl_data_inc)

         slvl_data_file = construct_file(slvl_file_fmt,                 &
                                  full_name=trim(forcing_filename), &
                                  record_length=rec_type_dbl,       &
                                  recl_words=nx_global*ny_global)

         call data_set(slvl_data_file, 'open_read')

         i_dim = construct_io_dim('i',nx_global)
         j_dim = construct_io_dim('j',ny_global)

         slvl_data_in = construct_io_field(trim(slvl_data_names(1)),  &
                        dim1=i_dim, dim2=j_dim,                   &
                        field_loc = slvl_bndy_loc(1),               &
                        field_type = slvl_bndy_type(1),             &
                        d2d_array = SEA_LEVEL_DATA(:,:,:,1,n))

         call data_set (slvl_data_file, 'define', slvl_data_in)
         call data_set (slvl_data_file, 'read',   slvl_data_in)
         call data_set (slvl_data_file, 'close')
         call destroy_io_field(slvl_data_in)
         call destroy_file(slvl_data_file)

         if (my_task == master_task) then
            write(stdout,blank_fmt)
            write(stdout,'(A,A)') ' Sea Level n-hour file read: ', forcing_filename
         endif
      enddo

      if (slvl_data_renorm(1) /= c1) SEA_LEVEL_DATA = &
                                   SEA_LEVEL_DATA*slvl_data_renorm(1)

   case default

      call exit_POP(sigAbort,'init_slvl: Unknown value for slvl_data_type')

   end select

!-----------------------------------------------------------------------
!
!  now check interpolation period (slvl_interp_freq) to set the
!    time for the next temporal interpolation (slvl_interp_next).
!
!  if no interpolation is to be done, set next interpolation time
!    to a large number so the atm pressure update test in routine
!    set_surface_forcing will always be false.
!
!  if interpolation is to be done every n-hours, find the first
!    interpolation time greater than the current time.
!
!  if interpolation is to be done every timestep, set next interpolation
!    time to a large negative number so the atm pressure update test in
!    routine set_surface_forcing will always be true.
!
!-----------------------------------------------------------------------

   select case (slvl_interp_freq)

   case ('never')

      slvl_interp_next = never
      slvl_interp_last = never
      slvl_interp_inc  = c0

   case ('n-hour')

      call find_interp_time(slvl_interp_inc, slvl_interp_next)

   case ('every-timestep')

      slvl_interp_next = always
      slvl_interp_inc  = c0

   case default

      call exit_POP(sigAbort, &
                    'init_slvl: Unknown value for slvl_interp_freq')

   end select

   if (nsteps_total == 0) slvl_interp_last = thour00

!-----------------------------------------------------------------------
!
!  echo forcing options to stdout.
!
!-----------------------------------------------------------------------

   slvl_data_label = 'Sea level assimilation'
   call echo_forcing_options(slvl_data_type,   slvl_formulation,  &
                             slvl_data_inc,    slvl_interp_freq,  &
                             slvl_interp_type, slvl_interp_inc,   &
                             slvl_data_label)

   if (my_task == master_task) then
     write(stdout,'(A)') 'sea level mask options: '
     write(stdout,*) 'sea level mask: ',use_slvl_mask
     write(stdout,'(A,A)') 'sea level mask file: ',trim(slvl_mask_file)
   endif

!-----------------------------------------------------------------------
!EOC


   if (use_slvl_mask) then
      allocate(SLVL_MASK(nx_block,ny_block,max_blocks_clinic))
      allocate(SLVL_OUT_MASK(nx_block,ny_block,max_blocks_clinic))
      forcing_filename = slvl_mask_file

      slvl_data_file = construct_file(slvl_mask_file_fmt,  &
                                    full_name=trim(forcing_filename),  &
                                    record_length=rec_type_dbl,        &
                                    recl_words=nx_global*ny_global)

      call data_set(slvl_data_file, 'open_read')

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      slvl_data_in = construct_io_field('SEA_LEVEL_MASK',          &
                             dim1=i_dim, dim2=j_dim,                   &
                             field_loc = field_loc_center,             &
                             field_type = field_type_scalar,           &
                             d2d_array = SLVL_MASK)

      call data_set (slvl_data_file, 'define', slvl_data_in)
      call data_set (slvl_data_file, 'read',   slvl_data_in)
      call data_set (slvl_data_file, 'close')
      call destroy_file(slvl_data_file)
      if (my_task == master_task) then
        allocate(WORK(nx_global,ny_global))
        WORK = 1._r8
        WORK(:,510:ny_global) = 0._r8
        
      endif
       
      call scatter_global(SLVL_OUT_MASK,WORK, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      if (my_task == master_task) deallocate(WORK)
      
  else
      allocate(SLVL_MASK(nx_block,ny_block,max_blocks_clinic))
      SLVL_MASK = 1._r8
      SLVL_OUT_MASK = 1._r8
  endif


 end subroutine init_slvl

!***********************************************************************
!BOP
! !IROUTINE: set_slvl
! !INTERFACE:

 subroutine set_slvl(SEA_LEVEL)

! !DESCRIPTION:
!  Updates the current value of the sea level array (ATM\_PRESS) by
!  interpolating to the current time.  It may also be necessary to
!  use (read) new data values in the interpolation.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      SEA_LEVEL        !  external sea level at current timestep

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  check if new data is necessary for interpolation.  if yes, then
!  shuffle indices in SEA_LEVEL_DATA and slvl_data_time arrays.
!  and read in new data if necessary ('n-hour' case).  also
!  increment values of slvl_data_time_min_loc, slvl_data_next and
!  slvl_data_update.  then do the temporal interpolation.
!
!-----------------------------------------------------------------------

   if (thour00 >= slvl_interp_next .or. nsteps_run == 0) then

   select case(slvl_data_type)
   case ('monthly-equal','monthly-calendar')

      slvl_data_label = 'Sea Level Monthly'
      if (thour00 >= slvl_data_update)                                  &
         call update_forcing_data(slvl_data_time, slvl_data_time_min_loc, &
                                  slvl_interp_type, slvl_data_next,       &
                                  slvl_data_update, slvl_data_type,       &
                                  slvl_data_inc,                        &
                                  SEA_LEVEL_DATA(:,:,:,:,1:12),       &
                                  slvl_data_renorm,                     &
                                  slvl_data_label, slvl_data_names,       &
                                  slvl_bndy_loc, slvl_bndy_type,          &
                                  slvl_filename, slvl_file_fmt)

      call interpolate_forcing(SEA_LEVEL_DATA(:,:,:,:,0),          &
                               SEA_LEVEL_DATA(:,:,:,:,1:12),       &
                               slvl_data_time, slvl_interp_type,       &
                               slvl_data_time_min_loc,               &
                               slvl_interp_freq, slvl_interp_inc,      &
                               slvl_interp_next, slvl_interp_last,     &
                               nsteps_run)

      SEA_LEVEL = SEA_LEVEL_DATA(:,:,:,1,0)

   case ('n-hour')


    

      slvl_data_label = 'Sea Level n-hour'
      if (thour00 >= slvl_data_update)                                  &
         call update_forcing_data(slvl_data_time, slvl_data_time_min_loc, &
                                  slvl_interp_type, slvl_data_next,       &
                                  slvl_data_update, slvl_data_type,       &
                                  slvl_data_inc,                        &
                           SEA_LEVEL_DATA(:,:,:,:,1:slvl_interp_order), &
                                  slvl_data_renorm,                     &
                                  slvl_data_label, slvl_data_names,       &
                                  slvl_bndy_loc, slvl_bndy_type,          &
                                  slvl_filename, slvl_file_fmt)

!    write(*,*) 'SEA_LEVEL_DATA/slvl_data_names(1)/slvl_bndy_loc  (1)/slvl_bndy_type (1) = ',&
!              slvl_data_names,slvl_bndy_loc ,slvl_bndy_type 


      call interpolate_forcing(SEA_LEVEL_DATA(:,:,:,:,0),             &
                           SEA_LEVEL_DATA(:,:,:,:,1:slvl_interp_order), &
                               slvl_data_time, slvl_interp_type,          &
                               slvl_data_time_min_loc,                  &
                               slvl_interp_freq, slvl_interp_inc,         &
                               slvl_interp_next, slvl_interp_last,        &
                               nsteps_run)

      SEA_LEVEL = SEA_LEVEL_DATA(:,:,:,1,0)

   end select

   if (nsteps_run /= 0) slvl_interp_next = slvl_interp_next + slvl_interp_inc

   endif ! thour00 > slvl_interp_next

!-----------------------------------------------------------------------
!EOC

 end subroutine set_slvl

!***********************************************************************

 end module forcing_slvl

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
