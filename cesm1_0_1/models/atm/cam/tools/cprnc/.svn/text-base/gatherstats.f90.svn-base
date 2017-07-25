module gatherstats

   use prec,          only: r4, r8
   use specialvalues, only: uninit_int, uninit_r8
   use cntlvars,      only: iprs, ipre, jprs, jpre, kprs, kpre
   use utils,         only: endrun

   implicit none

PRIVATE

   include 'netcdf.inc'

   public :: varinfo                   ! var stats datatype
   public :: diffinfo                  ! diff datatype
   public :: zero_stats                ! gather non-comparison stats (max, min, avg, nfill)
   public :: zero_comparison_stats     ! gather non-comparison stats (max, min, avg, nfill)
   public :: gather_stats              ! gather non-comparison stats (max, min, avg, nfill)
   public :: gather_comparison_stats   ! gather comparison stats (max diff, max relative diff, ...)
   public :: print_nullstats           ! print var stats
   public :: print_stats               ! print var stats
   public :: print_comparison_stats    ! print compares

   type varinfo
      character(len=NF_MAX_NAME) :: name          ! variable name
      character(len=NF_MAX_NAME) :: dimnames(4)   ! spatial dimension names
 
      integer :: dimsizes(4)   ! total dimension sizes
      integer :: npossible     ! number of field values
      integer :: nfill         ! number of fill values
      integer :: imax          ! 1st dimension index where field max occurs
      integer :: jmax          ! 2nd dimension index where field max occurs
      integer :: kmax          ! 3rd dimension index where field max occurs
      integer :: imin          ! 1st dimension index where field min occurs
      integer :: jmin          ! 2nd dimension index where field min occurs
      integer :: kmin          ! 3rd dimension index where field min occurs
      integer :: lidx          ! time index for this case

      real(r8) :: accum        ! accumulated sum
      integer  :: nvalid       ! number to average over
      real(r8) :: avg          ! field average value (accum/navg)
      real(r8) :: arrmax       ! field max
      real(r8) :: arrmin       ! field min
   end type varinfo

   type diffinfo
      character(len=NF_MAX_NAME) :: name          ! variable name
      character(len=NF_MAX_NAME) :: dimnames(4)   ! spatial dimension names

      real(r8) :: dmax         ! max diff
      real(r8) :: rdmax        ! max relative diff
      real(r8) :: dvals(2)     ! field values where max diff occurs
      real(r8) :: rdvals(2)    ! field values where max relative diff occurs
      real(r8) :: rms          ! rms difference
      real(r8) :: mwrms        ! mass-weighted rms diff (unimplemented)
      real(r8) :: rdbar        ! mean relative difference
      real(r8) :: rdlnbar      ! used in computation of number of matching digits

      integer :: npossible     ! number of non-fillvalue locations
      integer :: ndiff         ! number of diffs
      integer :: code          ! error code
      integer :: imax          ! 1st dimension index where biggest diff occurs
      integer :: jmax          ! 2nd dimension index where biggest diff occurs
      integer :: kmax          ! 3rd dimension index where biggest diff occurs
      integer :: irmax         ! 1st dimension index where biggest relative diff occurs
      integer :: jrmax         ! 2nd dimension index where biggest relative diff occurs
      integer :: krmax         ! 3rd dimension index where biggest relative diff occurs
   end type diffinfo

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine zero_stats (varstats)
!-------------------------------------------------------------------------------------------
! Purpose: Zero out varstats type for future gathering and accumulating
!
! Method:  This method zeros out the varstats type
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(varinfo), intent(out) :: varstats                       ! stats will be put into this struct
!
! Local workspace
!
      real(r8) :: dummy         ! dummy

!
! Zero varstats
!
      dummy = 1.0_r8   ! value not important, may not even need a value

      varstats%name        = ' '
      varstats%dimnames(:) = ' '
      varstats%dimsizes(:) = 0
      varstats%npossible   = 0
      varstats%nfill       = 0

      varstats%imax      = 0
      varstats%jmax      = 0
      varstats%kmax      = 0
      varstats%imin      = 0
      varstats%jmin      = 0
      varstats%kmin      = 0
      varstats%lidx      = 0
      varstats%accum     = 0.0_r8
      varstats%nvalid    = 0
      varstats%avg       = 0.0_r8
      varstats%arrmax    = -huge(dummy)
      varstats%arrmin    =  huge(dummy)

      return
   end subroutine zero_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine zero_comparison_stats (diffstats,code)
!-------------------------------------------------------------------------------------------
! Purpose: Zero out diffstats type for future gathering and accumulating
!
! Method:  This method zeros out the diffstats type
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(diffinfo), intent(out) :: diffstats                       ! stats will be put into this struct
      integer,optional,intent(in) :: code
!
! Local workspace
!
      real(r8) :: dummy         ! dummy

!
! Zero diffstats
!
      dummy = 1.0_r8   ! value not important, may not even need a value

      diffstats%name        = ' '
      diffstats%dimnames(:) = ' '
      diffstats%dmax        = -huge(dummy)
      diffstats%rdmax       = -huge(dummy)
      diffstats%npossible   = 0
      diffstats%ndiff       = 0
      diffstats%code        = 0
      diffstats%imax        = uninit_int
      diffstats%jmax        = uninit_int
      diffstats%kmax        = uninit_int
      diffstats%irmax       = uninit_int
      diffstats%jrmax       = uninit_int
      diffstats%krmax       = uninit_int
      diffstats%dvals(:)    = uninit_r8
      diffstats%rdvals(:)   = uninit_r8
      diffstats%rms         = 0._r8
      diffstats%rdbar       = 0._r8
      diffstats%rdlnbar     = 0._r8

      if (present(code)) then
         diffstats%code = code
      endif

      return

   end subroutine zero_comparison_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine gather_stats (name, ndims, arrsiz1, arrsiz2, arrsiz3, indxt, &
                            arr, fillvalue, doprint, dimnames, varstats)
!-------------------------------------------------------------------------------------------
! Purpose: Gather max, min, avg, and associated position info for this variable
!
! Method:  Check only those points where _FillValue does not apply. Copy output into
!          varstats.  This subroutine is now accumulating, varstats is in/out.
!          Accumulates on top of what's already in varstats.  To zero this out, 
!          call zero_stats.
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name                         ! variable name
      integer, intent(in) :: ndims                                 ! number of dimensions
      integer, intent(in) :: arrsiz1                               ! 1st non-time dimension size
      integer, intent(in) :: arrsiz2                               ! 2nd non-time dimension size (1 if absent)
      integer, intent(in) :: arrsiz3                               ! 3rd non-time dimension size (1 if absent)
      integer, intent(in) :: indxt                                 ! time index
      double precision, intent(in) :: arr(arrsiz1,arrsiz2,arrsiz3) ! array to be analyzed
      real(r8), intent(in) :: fillvalue                            ! fillvalue (if applicable)
      logical, intent(in) :: doprint                               ! whether or not to print requested data
      character(len=*), intent(in) :: dimnames(:)                  ! dimension names

      type(varinfo), intent(inout) :: varstats                       ! stats will be put into this struct
!
! Local workspace
!
      integer :: i, j, k        ! spatial indices

      logical :: chkfillvalue   ! whether to check for fillvalue
      logical :: validv         ! whether valid (i.e. non-fill) values exist


      call zero_stats(varstats)
!
! Fillvalue will contain uninit_r8 if it does not apply to this variable
!
      chkfillvalue = fillvalue /= uninit_r8

      do k=1,arrsiz3
         do j=1,arrsiz2
            do i=1,arrsiz1
!
! First print data if so requested
!
               if (doprint) then
                  if (i >= iprs .and. i <= ipre .and. &
                      (arrsiz2 == 1 .or. (j >= jprs .and. j <= jpre)) .and. &
                      (arrsiz3 == 1 .or. (k >= kprs .and. k <= kpre))) then
                     write(6,100)i, j, k, arr(i,j,k)
100                  format(' i,j,k=',3i5,' arr=',1p,e23.15)
                  end if
               end if

               validv = .not. chkfillvalue .or. (chkfillvalue .and. arr(i,j,k) /= fillvalue)
               if (validv) then
                  varstats%accum = varstats%accum + abs (arr(i,j,k))
                  varstats%nvalid = varstats%nvalid + 1
                  if (arr(i,j,k) > varstats%arrmax) then
                     varstats%arrmax = arr(i,j,k)
                     varstats%imax   = i
                     varstats%jmax   = j
                     varstats%kmax   = k
                  end if
               
                  if (arr(i,j,k) < varstats%arrmin) then
                     varstats%arrmin = arr(i,j,k)
                     varstats%imin   = i
                     varstats%jmin   = j
                     varstats%kmin   = k
                  end if
               else
                  varstats%nfill = varstats%nfill + 1
               end if
            end do
         end do
      end do
!
! Copy output to varstats
!
      varstats%name         = name
      varstats%dimnames(:)  = dimnames(:)
      varstats%npossible    = varstats%npossible + arrsiz1*arrsiz2*arrsiz3
      varstats%lidx         = indxt
      if (varstats%nvalid > 0) then
         varstats%avg       = varstats%accum / varstats%nvalid
      else
         varstats%avg       = fillvalue
      end if

      return
   end subroutine gather_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine gather_comparison_stats (name, ndims, dimsiz1, dimsiz2, dimsiz3, &
                                       arr, fillvalue, dimnames, diffstats, dodiffer)
!-------------------------------------------------------------------------------------------
! Purpose: Gather difference stats for the 2 instances of this variable. Stats include
!          number of diffs, largest absolute diff, largest relative diff, position 
!          info, and rms diff.
!
! Method:  Check only those points where _FillValue does not apply. Copy output into
!          diffstats. Also, set dodiffer to indicate whether any diffs were found
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*), intent(in) :: name                           ! variable name
      integer, intent(in) :: ndims                                   ! number of dimensions                     
      integer, intent(in) :: dimsiz1                                 ! 1st non-time dimension size              
      integer, intent(in) :: dimsiz2                                 ! 2nd non-time dimension size (1 if absent)
      integer, intent(in) :: dimsiz3                                 ! 3rd non-time dimension size (1 if absent)
!
! Declare arrays as double precision not r8 for netcdf convenience
!
      double precision, intent(in) :: arr(dimsiz1,dimsiz2,dimsiz3,2) ! array to be analyzed
      real(r8), intent(in) :: fillvalue(2)                           ! fillvalue (if applicable)
      character(len=*), intent(in) :: dimnames(4)                    ! dimension names
      real(r8) :: diff(dimsiz1)                ! difference
      real(r8) :: rdiff(dimsiz1)               ! relative difference

      type(diffinfo), intent(out) :: diffstats                       ! difference stats will be put into this struct
      logical, intent(out) :: dodiffer                               ! flag indicates whether any differences exist
!
! Local workspace
!
      integer :: i, ii, j, k                   ! spatial indices
      integer :: isub(1)                       ! i-index
      integer :: indx(dimsiz1)                 ! indices where there are diffs
      integer :: nval                          ! intermediate number of diffs

      real(r8) :: meanvertwt                   ! mean vertwt between 2 cases
      real(r8) :: denom                        ! denominator of expression

      logical :: chkfillvalue1, chkfillvalue2  ! whether to check for fillvalue
      logical :: validv1, validv2              ! whether an element equals fillvalue or not

      call zero_comparison_stats(diffstats)

      chkfillvalue1 = fillvalue(1) /= uninit_r8
      chkfillvalue2 = fillvalue(2) /= uninit_r8

      do k=1,dimsiz3
         do j=1,dimsiz2
            nval = 0
            do i=1,dimsiz1
!
! First print data if so requested
!
               if (i >= iprs .and. i <= ipre .and. j >= jprs .and. j <= jpre .and. &
                   (dimsiz3 == 1 .or. k >= kprs .and. k <= kpre)) then
                  write(6,100)i, j, k, arr(i,j,k,1), arr(i,j,k,2), arr(i,j,k,1) - arr(i,j,k,2)
100               format(' i,j,k=',3i5,' arr1,arr2,diff=',1p,3e23.15)
               end if

               diffstats%npossible = diffstats%npossible + 1
               validv1 = .not. chkfillvalue1 .or. (chkfillvalue1 .and. arr(i,j,k,1) /= fillvalue(1))
               validv2 = .not. chkfillvalue2 .or. (chkfillvalue2 .and. arr(i,j,k,2) /= fillvalue(2))

               if (validv1 .and. validv2) then
                  diff(i) = abs (arr(i,j,k,1) - arr(i,j,k,2))
                  diffstats%rms = diffstats%rms + diff(i)**2
                  if (diff(i) /= 0.) then
                     nval = nval + 1
                     indx(nval) = i
                  end if
               else
!
! Setting diff to zero here means locations with fillvalue in one file and valid data
! in the other file will be treated as not different.  This is probably not the best
! way to handle things.
!
                  diff(i) = 0.
               end if
            end do
!
! Save max diff info
! 
            if (nval > 0) then                     ! at least 1 difference
               diffstats%ndiff = diffstats%ndiff + nval
               isub = maxloc (diff(:))
               if (diff(isub(1)) > diffstats%dmax) then      ! Save max diff info
                  diffstats%dmax = diff(isub(1))
                  diffstats%dvals(1) = arr(isub(1),j,k,1)
                  diffstats%dvals(2) = arr(isub(1),j,k,2)
                  diffstats%imax = isub(1)
                  diffstats%jmax = j
                  diffstats%kmax = k
               end if
!
! Save max relative diff info
!
               rdiff(:) = 0.
               do ii=1,nval                              ! Compute relative diffs
                  i = indx(ii)
                  denom = max (abs (arr(i,j,k,1)), abs (arr(i,j,k,2)))
                  rdiff(i) = diff(i)/(2.*denom)
                  diffstats%rdbar = diffstats%rdbar + rdiff(i)
                  diffstats%rdlnbar = diffstats%rdlnbar - log10 (rdiff(i))
               end do
               isub = maxloc (rdiff(:))
            
               if (rdiff(isub(1)) > diffstats%rdmax) then          ! Save max relative diff info
                  diffstats%rdmax = rdiff(isub(1))
                  diffstats%rdvals(1) = arr(isub(1),j,k,1)         ! Save values & indices
                  diffstats%rdvals(2) = arr(isub(1),j,k,2)
                  diffstats%irmax = isub(1)
                  diffstats%jrmax = j
                  diffstats%krmax = k
               end if
            end if          ! nval > 0
         end do             ! j=1,dimsiz2
      end do                ! k=1,dimsiz3
!
! Final normalization
!
      diffstats%rms   = sqrt (diffstats%rms/diffstats%npossible)
      diffstats%rdbar = diffstats%rdbar/diffstats%npossible
!
! Copy values into struct
!
      diffstats%name        = name
      diffstats%dimnames(:) = dimnames(:)

      dodiffer = diffstats%ndiff > 0

      return
   end subroutine gather_comparison_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_nullstats (name,code)
!-------------------------------------------------------------------------------------------
! Purpose: Print various statistics (max, min, etc.) contained in input
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      character(len=*),intent(in) :: name
      integer :: code

      if (code == 1) then
         write(6,803) trim(name),' NOT ANALYZED'
         write(6,*) ' '
      elseif (code == 2) then
!         write(6,803) trim(name),' ARE ON BOTH FILES BUT DIFFERENT SIZES'
!         write(6,*) ' '
      endif

      return

803   format(1x,a10,1x,a)
   end subroutine print_nullstats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_stats (varstats)
!-------------------------------------------------------------------------------------------
! Purpose: Print various statistics (max, min, etc.) contained in input
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(varinfo), intent(in) :: varstats   ! derived type containing all the info
      character(len=NF_MAX_NAME) :: dnames(4)
      integer :: n,m,len

      do n = 1,4
         dnames(n) = trim(varstats%dimnames(n))
         len = len_trim(dnames(n))
      enddo

      write(6,809)  trim(varstats%name), &
                    (trim(dnames(n)),n=1,4), varstats%lidx
      write(6,810) varstats%npossible, &
                   varstats%imax, varstats%jmax, varstats%kmax, &
                   varstats%imin, varstats%jmin, varstats%kmin
      write(6,803) varstats%nvalid, varstats%arrmax, varstats%arrmin
      write(6,805) varstats%avg
      
      return

803   format(12x,i8,1x,1p2e23.15)
805   format(12x,'avg abs field values: ',1pe23.15,/)
809   format(1x,a,3x,'(',a,',',a,',',a,',',a,')',4x,'t_index = ',i6)
810   format(12x,i8,2x,'(',i6,',',i6,',',i6,') (',i6,',',i6,',',i6,')')
   end subroutine print_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_comparison_stats (varstats, diffstats)
!-------------------------------------------------------------------------------------------
! Purpose: Print various statistics (number of diffs, biggest diff, etc.) contained in input
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      type(varinfo), intent(in) :: varstats(2)  ! derived type for single-variable stats
      type(diffinfo), intent(in) :: diffstats   ! derived type for difference stats
!
! Local workspace
!
      real(r8) :: digbar      ! mean number of digits of agreement
      real(r8) :: digworst    ! number of digits of agreement for worst poin in the domain
      character(len=NF_MAX_NAME) :: dnames(4)
      integer :: n,m,len
      character(len=64) :: message
      real(r8) :: rmsout

      message = ''
      rmsout = diffstats%rms
      if (diffstats%code == -1) then
         message = 'comparison not possible, field shape/size difference'
         rmsout = huge(1.0_r8)
      endif

      do n = 1,4
         dnames(n) = trim(varstats(1)%dimnames(n))
         len = len_trim(dnames(n))
      enddo

      write(6,809) trim(varstats(1)%name), &
                   (trim(dnames(n)),n=1,4), varstats(1)%lidx, varstats(2)%lidx

      if (diffstats%ndiff > 0) then
         write(6,810)diffstats%ndiff, varstats(1)%npossible, &
                     varstats(1)%imax, varstats(1)%jmax, varstats(1)%kmax, &
                     varstats(1)%imin, varstats(1)%jmin, varstats(1)%kmin, &
                     diffstats%imax,  diffstats%jmax,  diffstats%kmax, &
                     diffstats%irmax, diffstats%jrmax, diffstats%krmax
         write(6,803) varstats(1)%nvalid, varstats(1)%arrmax, varstats(1)%arrmin, &
                     diffstats%dmax, diffstats%dvals(1), diffstats%rdmax, diffstats%rdvals(1), &
                                     varstats(2)%nvalid, varstats(2)%arrmax, varstats(2)%arrmin, &
                                     diffstats%dvals(2),                  diffstats%rdvals(2)
         write(6,811)varstats(2)%npossible, &
                     varstats(2)%imax, varstats(2)%jmax, varstats(2)%kmax, &
                     varstats(2)%imin, varstats(2)%jmin, varstats(2)%kmin
!
! Compute # digits accuracy for worst case & avg differences
!
         digbar = diffstats%rdlnbar / diffstats%ndiff
         digworst = log10 (1. / diffstats%rdmax)

         write(6,805) varstats(1)%avg, diffstats%rms, diffstats%rdbar, &
                      varstats(2)%avg, digbar, digworst
      else
         write(6,810)diffstats%ndiff, varstats(1)%npossible, &
                     varstats(1)%imax, varstats(1)%jmax, varstats(1)%kmax, &
                     varstats(1)%imin, varstats(1)%jmin, varstats(1)%kmin
         write(6,814) varstats(1)%nvalid, varstats(1)%arrmax, varstats(1)%arrmin, &
                     trim(message), &
                                       varstats(2)%nvalid, varstats(2)%arrmax, varstats(2)%arrmin
         write(6,811)varstats(2)%npossible, &
                     varstats(2)%imax, varstats(2)%jmax, varstats(2)%kmax, &
                     varstats(2)%imin, varstats(2)%jmin, varstats(2)%kmin
         write(6,815)varstats(1)%avg, varstats(2)%avg
      end if
!
! The following is strictly for test-model: Actual mass-weighting not yet implemented
!
      write(6,'(a,a8,1pe11.4,/)') ' RMS ', varstats(1)%name, rmsout
      return
!
! formats
!
803   format(12x,      i8,1x,1p2e23.15,2x,e8.1,e23.15,e8.1,e23.15,/, &
             12x,      i8,1x,  2e23.15,2x,8x,  e23.15,8x,  e23.15)
814   format(12x,      i8,1x,1p2e23.15,10x,a,/, &
             12x,      i8,1x,  2e23.15)
809   format(1x,a,3x,'(',a,',',a,',',a,',',a,')',4x,'t_index = ',2i6)
810   format(3x,i8,1x,i8,2x,'(',i6,',',i6,',',i6,') (',i6,',',i6,',',i6,')', &
             11x,'(',i6,',',i6,',',i6,')',9x,'(',i6,',',i6,',',i6,')')
811   format(12x,     i8,2x,'(',i6,',',i6,',',i6,') (',i6,',',i6,',',i6,')', &
             11x,'(',i6,',',i6,',',i6,')',9x,'(',i6,',',i6,',',i6,')')

805   format(10x,'avg abs field values:  ',1pe23.15,4x,'rms diff:',e8.1, &
              3x,'avg rel diff(npos): ',e8.1,/, &
             10x,'                       ',  e23.15,24x, &
             'avg decimal digits(ndif): ',0p,f4.1,' worst: ',f4.1)
815   format(10x,'avg abs field values:  ',1pe23.15,/, &
             10x,'                       ',  e23.15)
   end subroutine print_comparison_stats

end module gatherstats
