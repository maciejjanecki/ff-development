module driver

   use prec,          only: r4, r8
   use cntlvars,      only: matchts, lnorm, constps
   use netcdfids,     only: timevid, timedid
   use specialvalues, only: uninit_r8
   use gatherstats,   only: varinfo, diffinfo, gather_stats, gather_comparison_stats, &
                            print_stats, print_comparison_stats, print_nullstats, &
                            zero_comparison_stats
   use nfwrappers,    only: wrap_inq_var, wrap_get_vara_double, wrap_inq_dimlen, &
                            wrap_inq_dimname, wrap_get_var_double
   use printers,      only: prhddiff, prheader
   use utils,         only: endrun

   implicit none

PRIVATE

   include 'netcdf.inc'

   public :: timeloop       ! driver routine which loops through time and calls stats-gathering routines
   private :: kosher        ! logical function determines whether field should be analyzed
   private :: get_fillvalue ! retrieve _FillValue attribute if present
   private :: setsizes      ! sets start, kount arrays for netcdf calls, and total array size
   integer,private :: debug = 1

CONTAINS

   subroutine timeloop (nfid, nvars, numcases, ntime)
!-------------------------------------------------------------------------------------------
! Purpose: Loop through time samples, computing and printing stats for each variable
!
! Method:  For each time index, loop through the variables and print stats (max, min, etc.).
!          Print comparison stats for fields which are on both files
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)             ! file id(s)
      integer, intent(in) :: nvars(2)            ! number of variables on each file
      integer, intent(in) :: numcases            ! number of files (1 or 2)
      integer, intent(in) :: ntime(2)            ! length of time dimension
!
! Local workspace
!
      real(r8), parameter :: timeepsilon = 1.e-9 ! time diff less than this considered insignificant
 
      character(len=NF_MAX_NAME) :: name         ! variable name
      character(len=NF_MAX_NAME) :: dimnames(4)  ! dimension names

      integer :: t(2)                   ! time indices
      integer :: t1, t2
      equivalence (t(1),t1)             ! equivalence for convenience
      equivalence (t(2),t2)

      integer :: n,nn                   ! file indices
      integer :: v                      ! variable index
      integer :: xtype                  ! variable type from netcdf
      integer :: numfld1                ! number of fields analyzed on file 1
      integer :: numfsk1                ! number of fields skipped on file 1
      integer :: numfld2                ! number of fields analyzed on file 2
      integer :: numfcpr                ! number of fields compared (total)
      integer :: numfdif                ! number of fields which differed (total)
      integer :: numfsha                ! number of fields which differed in shape
      integer :: ndims(2)               ! number of dimensions in variable
      integer :: dimsizes(4)            ! sizes of non-time (probably spatial) dimensions
      integer :: length(2)              ! size needed for buf allocation

      integer :: varid(2)               ! variable ids
      integer :: idum                   ! required arg to subroutine that is unused
      integer :: start(4,2)             ! shape input to get_vara routines
      integer :: kount(4,2)             ! shape input to get_vara routines
      integer :: dimids(4,2)            ! dimension ids
      integer :: natts                  ! number of attributes (required for netcdf calls)
      
      logical :: twocases               ! comparison or single file analysis?
      logical :: dodiffer               ! flag indicates differences were found
      logical :: dodshape               ! flag indicates different shape/size were found

      real(r8) :: fillvalue(2)          ! _FillValue (if applicable).  Points are ignored.
      real(r8) :: time1(ntime(1))       ! time variable file 1
      real(r8) :: time2(ntime(2))       ! time variable file 2
      real(r8),allocatable :: time1x(:),time2x(:)         ! current time1, time2
      integer  :: ntm,ktm               ! number of matching time indices, counter
      integer, allocatable :: tm1(:)    ! matching time indices from file 1
      integer, allocatable :: tm2(:)    ! matching time indices from file 2
      logical  :: done                  ! flag used in time index search
      logical  :: analyze               ! flag used to determine analyze mode
      real(r8) :: timediff              ! time difference between samples when there are 2 files
!
! Declare data array as double precision instead of r8 for netcdf consistency
!
      double precision, allocatable :: buf(:,:)    ! field values

      type (varinfo) :: varstats(numcases)         ! variable statistics
      type (diffinfo) :: diffstats                 ! difference statistics
!
! Retrieve time variable(s) and print header diffs if 2 files are being analyzed
!
      twocases = numcases > 1
      numfld1 = 0
      numfsk1 = 0
      numfld2 = 0
      numfcpr = 0
      numfdif = 0
      numfsha = 0

      ntm = 0
      if (twocases .and. timevid(1) > 0 .and. timevid(2) > 0) then
         call prhddiff (nfid, numcases)
         call wrap_get_var_double (nfid(1), timevid(1), time1)
         call wrap_get_var_double (nfid(2), timevid(2), time2)
         ntm = min(ntime(1),ntime(2))
         allocate(tm1(ntm),tm2(ntm))
         !--- Default is just compare by index
         do ktm = 1,ntm
            tm1(ktm) = ktm
            tm2(ktm) = ktm
         enddo
         if (matchts) then
            ntm = 0
            t1 = 1
            t2 = 1
            done = .false.
            do while (.not. done)
               timediff = abs (time1(t1) - time2(t2))
               if (timediff < timeepsilon) then
                  ntm = ntm + 1
                  tm1(ntm) = t1
                  tm2(ntm) = t2
                  t1 = t1 + 1
                  t2 = t2 + 1
               elseif (time1(t1) < time2(t2)) then
                  write(6,*) 'Skipping a time sample on file 1'
                  t1 = t1 + 1
                else
                  write(6,*) 'Skipping a time sample on file 2'
                  t2 = t2 + 1
               end if
               if (t1 > ntime(1)) done = .true.
               if (t2 > ntime(2)) done = .true.
            enddo
         endif
         write(6,*) 'found ',ntm,' common timesteps',tm1(1:ntm),' and ',tm2(1:ntm)
         allocate(time1x(ntime(1)),time2x(ntime(2)))
         time1x = time1
         time2x = time2
      elseif (.not.twocases .and. timevid(1) > 0) then
         call wrap_get_var_double (nfid(1), timevid(1), time1)
         ntm = ntime(1)
         allocate(tm1(ntm))
         do ktm = 1,ntm
            tm1(ktm) = ktm
         enddo
         write(6,*) 'file has ',ntm,' timesteps',tm1(1:ntm)
         allocate(time1x(ntm),time2x(ntm),tm2(ntm))
         time1x = time1
         time2x = 0.
         tm2 = 1
      else  ! pretend time doesn't exist, still have to loop through time index
         ntm = 1
         allocate(tm1(ntm),tm2(ntm))
         tm1 = 1
         tm2 = 1
         allocate(time1x(ntm),time2x(ntm))
         time1x = 0.
         time2x = 0.
         write(6,*) 'file does not have a time dimension'
      end if

!
! Top of loop over time samples: first match up times
!
      do ktm = 1,ntm
      t1 = tm1(ktm)
      if (twocases) then
         t2 = tm2(ktm)

         call prheader (nfid, numcases, t, ntime)
         write(6,800) t1, t2, time1x(t1), time2x(t2)
         write(6,801)
         write(6,802)
         write(6,803)
         write(6,804)
         write(6,805)
         write(6,*)
!
! Loop over file 1 variables, print comparison stats for those fields which are also on file 2
!
         do v=1,nvars(1)
            varid(1) = v
            dimids(:,1) = -1
            call wrap_inq_var (nfid(1), varid(1), name, xtype, ndims(1), dimids(:,1), natts)
!
! Skip to the next variable if type and shape criteria are not met (floating point, something by time)
!
            if (.not. kosher (ktm, nfid(1), name, ndims(1), dimids(:,1), timedid(1), xtype)) then
               numfsk1 = numfsk1 + 1
               call print_nullstats(name,1)
               cycle
            end if
!
! Compare only if field also exists on file 2, and has the same characteristics as the file 1 field
!
            if (nf_inq_varid (nfid(2), name, varid(2)) == NF_NOERR) then
               dimids(:,2) = -1
               call wrap_inq_var (nfid(2), varid(2), name, xtype, ndims(2), dimids(:,2), natts)
               if (kosher (ktm, nfid(2), name, ndims(2), dimids(:,2), timedid(2), xtype)) then
!
! Determine sizes, perform sanity check, then allocate buffer
!
                  do n=1,numcases
                     call setsizes (nfid(n), t(n), ndims(n), dimids(:,n), timedid(n), &
                                    start(:,n), kount(:,n), length(n), dimsizes, dimnames)
                  end do

                  if (shapesmatch (nfid, ndims, dimids) .and. length(1) == length(2)) then

                     dodshape = .false.
                     allocate (buf(length(1),numcases))
!
! Compute non-comparison stats for this variable on both files
!
                     do n=1,numcases
                        call wrap_get_vara_double (nfid(n), varid(n), start(:,n), kount(:,n), buf(:,n))
                        fillvalue(n) = get_fillvalue (nfid(n), varid(n))
                        call gather_stats (name, ndims(n), kount(1,n), kount(2,n), kount(3,n), t(n), &
                                           buf(:,n), fillvalue(n), .false., dimnames, varstats(n))
                     end do
!
! Gather comparison stats, then release buffer space
!
                     call gather_comparison_stats (name, ndims(1), kount(1,1), kount(2,1), kount(3,1), &
                                                   buf, fillvalue, dimnames, diffstats, dodiffer)
                     deallocate (buf)
!
! Accumulate stats on number of fields compared and number that differ
! Then print the comparison stats
!
                  else
!
! Compute non-comparison stats for this variable on both files
!
                     dodshape = .true.
                     do n=1,numcases
                        allocate (buf(length(n),numcases))
                        call wrap_get_vara_double (nfid(n), varid(n), start(:,n), kount(:,n), buf(:,n))
                        fillvalue(n) = get_fillvalue (nfid(n), varid(n))
                        call gather_stats (name, ndims(n), kount(1,n), kount(2,n), kount(3,n), t(n), &
                                           buf(:,n), fillvalue(n), .false., dimnames, varstats(n))
                        deallocate (buf)
                     end do
                     call zero_comparison_stats(diffstats,-1)
                     dodiffer = .true.
                  end if    ! shapes match
!
! Compute overall summary stats
!
                  numfcpr = numfcpr + 1
                  if (dodiffer) then
                     numfdif = numfdif + 1
                  end if
                  if (dodshape) then
                     numfsha = numfsha + 1
                  end if
                  call print_comparison_stats (varstats, diffstats)

               endif     ! kosher 2
            end if       ! variable exists on both files
         end do          ! v=1,nvars(1)
      end if             ! twocases
!
! Gather and print non-comparison stats (e.g. when analyzing only 1 file, or variable does
! not exist on both files)
!
      do n=1,numcases
         nn = mod(n,2) + 1
         if (n == 1) then
            write(6,900) n, t(n), time1x(t(n))
         else if (n == 2) then
            write(6,900) n, t(n), time2x(t(n))
         end if
         write(6,811)
         write(6,802)
         write(6,904)
         write(6,*)

         do v=1,nvars(n)
            varid(n) = v
            dimids(:,n) = -1
            call wrap_inq_var (nfid(n), varid(n), name, xtype, ndims(n), dimids(:,n), natts)
!
! If data are floating point and shape is something by time, then analyze if field is not 
! also on other the file. If it is on the other file, stats have already been printed in 
! the "do v=1,nvars(1)" loop above
!
            analyze = .true.
            if (.not.kosher (ktm, nfid(n), name, ndims(n), dimids(:,n), timedid(n), xtype)) analyze=.false.
            if (twocases .and. nf_inq_varid (nfid(nn), name, varid(nn)) == NF_NOERR) then
               dimids(:,nn) = -1
               call wrap_inq_var (nfid(nn), varid(nn), name, xtype, ndims(nn), dimids(:,nn), natts)
               if (kosher (ktm, nfid(nn), name, ndims(nn), dimids(:,nn), timedid(nn), xtype)) analyze=.false.
            endif

            if (analyze) then
               if (n == 1) numfld1 = numfld1 + 1
               if (n == 2) numfld2 = numfld2 + 1
!
! Determine sizes, allocate buffer, read in the variable, then gather and print stats
!                               
               call setsizes (nfid(n), t(n), ndims(n), dimids(:,n), timedid(n), &
                              start(:,n), kount(:,n), length(n), dimsizes, dimnames)

               allocate (buf(length(n),1))
               call wrap_get_vara_double (nfid(n), varid(n), start(:,n), kount(:,n), buf)
               fillvalue(n) = get_fillvalue (nfid(n), varid(n))

               call gather_stats (name, ndims(n), kount(1,n), kount(2,n), kount(3,n), t(n), &
                                  buf, fillvalue(n), .true., dimnames, varstats(n))
               deallocate (buf)

               call print_stats (varstats(n))
            end if
         end do     ! v=1,nvars(n)
      end do        ! n=1,numcases

      write(6,806)
!
! Move to next time index in file(s)
!
      enddo   ! ktm loop

!
! Summarize results
!
      write(6,*) ' '
      write(6,700) 'SUMMARY of cprnc:'
      if (twocases .and. (ntm < ntime(1) .or. ntm < ntime(2))) then
         do n=1,numcases
            if (ntm < ntime(n)) then
               write(6,701) '  A total number of ',ntime(n)-ntm, ' time samples on file ', n
            end if
         end do
      end if
      if (twocases) then
      write(6,700) '  A total number of ',numfcpr,' fields were compared'
      write(6,700) '          of which  ',numfdif,' has non-zero differences'
      write(6,700) '          including ',numfsha,' that had inconsistent size/shape'
      endif
      write(6,700) '  A total number of ',numfld1,' fields in file 1 were analyzed (non-compare mode)'
      write(6,700) '  A total number of ',numfsk1,' fields in file 1 could not be analyzed'
      if (twocases) then
      write(6,700) '  A total number of ',numfld2,' fields in file 2 were analyzed (non-compare mode)'
      if (numfcpr > 0 .and. numfdif == 0) then
         write(6,700) '  diff_test: the two files seem to be IDENTICAL '
      else
         write(6,700) '  diff_test: the two files seem to be DIFFERENT '
      endif
      endif
      write(6,*) ' '

700   format(a,i6,a)
701   format(a,i6,a,i2) 
800   format( 'SUMMARY of matching fields time samples ', 2i6, ' times ', 1p, 2e15.8, ':')
801   format( '  FIELD       (  dim1,  dim2,  dim3,  dim4)    t_index =   tid1  tid2 ')
811   format( '                      (  dim1,  dim2,  dim3,  dim4)')
802   format( '     NDIFFS   ARRSIZ1 ( indx1, indx2, indx3) file 1')
803   format( '              NVALID1',10x,'MAX1',18x,'MIN1',12x,'DIFFMAX',2x,'VALUES',16x,'RDIFMAX',2x,'VALUES')
804   format( '              NVALID2',10x,'MAX2',18x,'MIN2')
805   format( '              ARRSIZ2 ( indx1, indx2, indx3) file 2')
900   format( 'SUMMARY of fields only on file ', i1, ' time sample', i6, ' time ', 1p, e15.8, ':')
904   format( ' FIELD        NVALID ',10x,'MAX ',18x,'MIN')
806   format(132('*'))

   end subroutine timeloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setsizes (nfid, t, ndims, dimids, timedid, start, &
                        kount, length, dimsizes, dimnames)
!-------------------------------------------------------------------------------------------
! Purpose: Determine dimension sizes and total length, and set start and kount arrays
!          for input to nf_get_vara_double
!
! Method:  Call appropriate netcdf inquiry functions
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid                  ! file id
      integer, intent(in) :: t                     ! time index
      integer, intent(in) :: ndims                 ! number of dimensions in variable
      integer, intent(in) :: dimids(:)             ! dimension ids of variable
      integer, intent(in) :: timedid               ! time dimension id

      integer, intent(out) :: start(:)             ! for input to nf_get_vara routines
      integer, intent(out) :: kount(:)             ! for input to nf_get_vara routines
      integer, intent(out) :: length               ! total array size (excluding time index)
      integer, intent(out) :: dimsizes(:)          ! individual non-time dimension sizes

      character(len=*), intent(out) :: dimnames(:) ! dimension names
!
! Local workspace
!
      integer :: n                                 ! dimension index

!
! ndims > 4 should have been skipped and a msg printed by function kosher
!
      if (ndims > 4) then
         call endrun ('setsizes: ndims too big')
      end if

      start(:) = 1
      kount(:) = 1
      length = 1
!
! Loop through dimensions, setting start, kount, dimnames, dimsizes, and length for output
!
      do n=1,ndims
         call wrap_inq_dimlen (nfid, dimids(n), dimsizes(n))
         call wrap_inq_dimname (nfid, dimids(n), dimnames(n))
         if (dimids(n) == timedid) then
            start(n) = t
            kount(n) = 1
         else
            start(n) = 1
            kount(n) = dimsizes(n)
            length   = length * kount(n)
         endif
      end do
!
! Unlimited dimension: want 1 slice at time index t
!
!
! Set inapplicable dimension sizes to 1 (convenience for gather_stats routines)
! Set inapplicable dimension names to dashes for later printing (print_stats routines)
!
      dimsizes(ndims+1:4) = 1
      dimnames(ndims+1:4) = '-----'

      return
   end subroutine setsizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function kosher (kn, nfid, name, ndims, dimids, timedid_loc, xtype)
!-------------------------------------------------------------------------------------------
! Purpose: Determine whether input field has appropriate characterisics for analysis
!
! Method:  Examine input vars and make appropriate netcdf calls to see if the field
!          is floating point, dimensioned something by time, where "something" may
!          be 1, 2, or 3 dimensions
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: kn                  ! times through the loop
      integer, intent(in) :: nfid                ! file id
      character(len=*), intent(in) :: name       ! variable name
      integer, intent(in) :: ndims               ! number of dimensions associated with variable
      integer, intent(in) :: dimids(:)           ! dimension ids
      integer, intent(in) :: timedid_loc         ! _loc distinguishes from global variable
      integer, intent(in) :: xtype               ! variable type (netcdf)

      integer :: n                               ! dimension index
      integer :: dimlen                          ! dimension length

      if (xtype < 0) then
         call endrun ('kosher: bad input')
      end if
!
! No more than 3 dimensions plus time allowed
!
      if (ndims < 0 .or. ndims > 4) then
         kosher = .false.
         if (debug > 1) write(6,*)'kosher: variable ', name, ' must be 1-4 dimensions'
         return
      end if
!
! Ensure that all dimensions have length > 0
!
      do n=1,ndims
         call wrap_inq_dimlen (nfid, dimids(n), dimlen)
         if (dimlen == 0) then
            if (debug > 1) write(6,*) 'kosher: dimlen == 0',n
            kosher = .false.
            return
         end if
      end do

!
! Ensure that time dim is last if it exists
!
      do n=1,ndims-1
         if (dimids(n) == timedid_loc) then
            if (debug > 1) write(6,*) 'kosher: time dim not last',n
            kosher = .false.
            return
         end if
      end do
!
! Ensure that if kn > 1, then field has a time dim
!
      if (kn > 1) then
         kosher = .false.
         do n=1,ndims
            if (dimids(n) == timedid_loc) then
               kosher = .true.
            end if
         end do
         if (.not.kosher) then
           if (debug > 2) write(6,*) 'kosher: no time dim'
           return
         endif
      endif
!
! Variable must be of type float, and dimensioned something by time
!
      if (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) then
         kosher = .true.
      else
         if (debug > 1) write(6,*) 'kosher: type'
         kosher = .false.
         return
      end if

      return
   end function kosher

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function shapesmatch (nfid, ndims, dimids)
!-------------------------------------------------------------------------------------------
! Purpose: Determine whether shapes of 2 netcdf variables match
!
! Method:  Make netcdf calls to get dimension names and lengths. Both names AND sizes
!          must all match in order for shapes to match
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)               ! file ids
      integer, intent(in) :: ndims(2)              ! number of dimensions
      integer, intent(in) :: dimids(4,2)           ! dimension ids
!
! Local workspace
!
      integer :: n                                 ! dimension index
      integer :: dimlen1, dimlen2                  ! dimension lengths
      
      character(len=NF_MAX_NAME) :: dimname1, dimname2 ! dimension names
!
! Check on sanity of input
!
      if (ndims(1) /= ndims(2)) then
         shapesmatch = .false.
         return
      end if

      do n=1,ndims(1)
         call wrap_inq_dimname (nfid(1), dimids(n,1), dimname1)
         call wrap_inq_dimname (nfid(2), dimids(n,2), dimname2)
!
! If dimension names do not match return false
! One might allow the unlimited dimension name to differ, but we do not
!
         if (dimname1 /= dimname2) then
            shapesmatch = .false.
            return
         end if
!
! If dimension lengths other than the unlimited dimension do not match return false
!
         call wrap_inq_dimlen (nfid(1), dimids(n,1), dimlen1)
         call wrap_inq_dimlen (nfid(2), dimids(n,2), dimlen2)
         if (dimlen1 /= dimlen2) then
            if (dimids(n,1) /= timedid(1) .or. dimids(n,2) /= timedid(2)) then
               shapesmatch = .false.
               return
            end if
         end if
      end do

      shapesmatch = .true.
      return
   end function shapesmatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real(r8) function get_fillvalue (nfid, varid)
!-------------------------------------------------------------------------------------------
! Purpose: Determine _FillValue attribute for a given variable (if present)
!
! Method:  Make appropriate netcdf calls.  If _FillValue is not present, return uninit_r8
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid   ! file id
      integer, intent(in) :: varid  ! variable id
!
! Local workspace
!
      integer :: xtype               ! variable type (netcdf)
      integer :: numelem             ! number of elements in attribute
      integer :: ret                 ! return code

      double precision :: fillvalue  ! value of _FillValue attribute
!
! If fillvalue is there, determine its value and return it. Otherwise return uninit_r8
!
      get_fillvalue = uninit_r8
      if (nf_inq_att (nfid, varid, '_FillValue', xtype, numelem) == NF_NOERR) then
         if ((xtype == NF_REAL .or. xtype == NF_DOUBLE) .and. numelem == 1) then
            ret = nf_get_att_double (nfid, varid, '_FillValue', fillvalue)
            get_fillvalue = fillvalue
         end if
      end if

      return
   end function get_fillvalue
end module driver
