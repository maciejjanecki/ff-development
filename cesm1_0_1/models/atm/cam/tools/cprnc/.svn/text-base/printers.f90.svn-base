module printers

   use prec,          only: r8
   use specialvalues, only: uninit_int
   use nfwrappers,    only: wrap_inq_varid, wrap_get_var_int, wrap_get_vara_text, wrap_get_vara_int
   use netcdfids,     only: dateid, datesecid, nstephid, date_writtenid, time_writtenid

   implicit none

   include 'netcdf.inc'

PRIVATE

   public :: prheader                 ! print header info
   public :: prhddiff                 ! print header diffs

CONTAINS

   subroutine prheader (nfid, numcases, t, ntime)
!-------------------------------------------------------------------------------------------
! Purpose: Print header info for 1 or 2 cases
!
! Method:  Make netcdf calls to retrieve the values or print the value of the input argument
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)       ! file id(s)
      integer, intent(in) :: numcases      ! number of files
      integer, intent(in) :: t(2)          ! time index
      integer, intent(in) :: ntime(2)      ! number of time samples
!
! Local workspace
!
      character(len=8) :: date_written(2)  ! date stamp on the time sample
      character(len=8) :: time_written(2)  ! time of day stamp

      integer :: date(2)                   ! simulation date
      integer :: datesec(2)                ! simulation time of day
      integer :: nsteph(2)                 ! time step
      integer :: startc(2)                 ! start array for nf_get_vara_text calls
      integer :: kountc(2)                 ! kount array for nf_get_vara_text calls
      integer :: start(1)                  ! start array for nf_get_vara_int calls
      integer :: kount(1)                  ! kount array for nf_get_vara_int calls
      integer :: n                         ! file index
!
! Initialize variables in case they do not exist on the file
!
      date(:)    = uninit_int
      datesec(:) = uninit_int
      nsteph(:)  = uninit_int
      date_written(:) = 'N/A'
      time_written(:) = 'N/A'

      write(6,*)'SUMMARY OF HEADER INFO TIME SAMPLES', t(1), t(2),':'
      do n=1,numcases
         startc(:) = (/1,t(n)/)
         kountc(:) = (/8,1/)
         if (date_writtenid(n) > 0) call wrap_get_vara_text (nfid(n), date_writtenid(n), startc, kountc, date_written(n))
         if (time_writtenid(n) > 0) call wrap_get_vara_text (nfid(n), time_writtenid(n), startc, kountc, time_written(n))

         start(1) = t(n)
         kount(1) = 1
         if (nstephid(n) > 0)  call wrap_get_vara_int (nfid(n), nstephid(n), start, kount, nsteph(n))
         if (dateid(n) > 0)    call wrap_get_vara_int (nfid(n), dateid(n), start, kount, date(n))
         if (datesecid(n) > 0) call wrap_get_vara_int (nfid(n), datesecid(n), start, kount, datesec(n))
      end do

      if (numcases > 1) then
         write(6,800)date_written(1), date_written(2), &
                     time_written(1), time_written(2), &
                     t(1), t(2), &
                     ntime(1), ntime(2), &
                     nsteph(1), nsteph(2), &
                     date(1), date(2), &
                     datesec(1),  datesec(2)
      else
         write(6,801)date_written(1), &
                     time_written(1), &
                     t(1), &
                     ntime(1), &
                     nsteph(1), &
                     date(1), &
                     datesec(1)
      end if
      write(6,*)

      return

800   format(' date_written:                ',2(A8,1X),/, &
             ' time_written:                ',2(A8,1X),/, &
             ' time sample:                 ',2I9,/,      &
             ' total possible time samples: ',2I9,/,      &
             ' NSTEPH:                      ',2I9,/,      &
             ' DATE:                        ',2I9,/,      &
             ' DATESEC:                     ',2I9)

801   format(' date_written:                ',1(A8,1X),/, &
             ' time_written:                ',1(A8,1X),/, &
             ' time sample:                 ',1I9,/,      &
             ' total possible time samples: ',1I9,/,      &
             ' NSTEPH:                      ',1I9,/,      &
             ' DATE:                        ',1I9,/,      &
             ' DATESEC:                     ',1I9)
   end subroutine prheader

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine prhddiff (nfid, numcases)
!-------------------------------------------------------------------------------------------
! Purpose: Print header difference info between 2 files
!
! Method:  Make netcdf calls to retrieve the values. Print difference info if they do not
!          match
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)        ! file id(s)
      integer, intent(in) :: numcases       ! number of files
!
! Local workspace
!
      integer :: n                          ! file index
      integer :: ntrmid(2)                  ! netcdf id for ntrm
      integer :: ntrnid(2)                  ! netcdf id for ntrn
      integer :: ntrkid(2)                  ! netcdf id for ntrk
      integer :: ndbaseid(2)                ! netcdf id for ndbase
      integer :: nsbaseid(2)                ! netcdf id for nsbase
      integer :: nbdateid(2)                ! netcdf id for nbdate
      integer :: nbsecid(2)                 ! netcdf id for nbsec
      integer :: mdtid(2)                   ! netcdf id for mdt

      integer :: ntrm(2)       = uninit_int ! spectral truncation parameter M
      integer :: ntrn(2)       = uninit_int ! spectral truncation parameter N
      integer :: ntrk(2)       = uninit_int ! spectral truncation parameter K
      integer :: ndbase(2)     = uninit_int ! base day
      integer :: nsbase(2)     = uninit_int ! base seconds of base day
      integer :: nbdate(2)     = uninit_int ! base date
      integer :: nbsec(2)      = uninit_int ! base seconds of base date
      integer :: mdt(2)        = uninit_int ! time step
!
! Loop through files, get variable info (if available), and print if diffs are found
!
      do n=1,numcases
         call wrap_inq_varid (nfid(n), 'ntrm', ntrmid(n), -1)
         call wrap_inq_varid (nfid(n), 'ntrn', ntrnid(n), -1)
         call wrap_inq_varid (nfid(n), 'ntrk', ntrkid(n), -1)
         call wrap_inq_varid (nfid(n), 'ndbase', ndbaseid(n), -1)
         call wrap_inq_varid (nfid(n), 'nsbase', nsbaseid(n), -1)
         call wrap_inq_varid (nfid(n), 'nbdate', nbdateid(n), -1)
         call wrap_inq_varid (nfid(n), 'nbsec', nbsecid(n), -1)
         call wrap_inq_varid (nfid(n), 'mdt', mdtid(n), -1)
!
! Cannot use optional arg to get_var_int or it complains about scalar/array problems
!
         if (ntrmid(n) > 0)   call wrap_get_var_int (nfid(n), ntrmid(n), ntrm(n))
         if (ntrnid(n) > 0)   call wrap_get_var_int (nfid(n), ntrnid(n), ntrn(n))
         if (ntrkid(n) > 0)   call wrap_get_var_int (nfid(n), ntrkid(n), ntrk(n))
         if (ndbaseid(n) > 0) call wrap_get_var_int (nfid(n), ndbaseid(n), ndbase(n))
         if (nsbaseid(n) > 0) call wrap_get_var_int (nfid(n), nsbaseid(n), nsbase(n))
         if (nbdateid(n) > 0) call wrap_get_var_int (nfid(n), nbdateid(n), nbdate(n))
         if (nbsecid(n) > 0)  call wrap_get_var_int (nfid(n), nbsecid(n), nbsec(n))
         if (mdtid(n) > 0)    call wrap_get_var_int (nfid(n), mdtid(n), mdt(n))
      end do

      write(6,*)'SUMMARY OF TIME-INDEPENDENT HEADER DIFFERENCES:'
      if (ntrm(1)   /= ntrm(2)  ) write(6,100)'ntrm  :', ntrm(1)  , ntrm(2)  
      if (ntrn(1)   /= ntrn(2)  ) write(6,100)'ntrn  :', ntrn(1)  , ntrn(2)  
      if (ntrk(1)   /= ntrk(2)  ) write(6,100)'ntrk  :', ntrk(1)  , ntrk(2)  
      if (ndbase(1) /= ndbase(2)) write(6,100)'ndbase:', ndbase(1), ndbase(2)
      if (nsbase(1) /= nsbase(2)) write(6,100)'nsbase:', nsbase(1), nsbase(2)
      if (nbdate(1) /= nbdate(2)) write(6,100)'nbdate:', nbdate(1), nbdate(2)
      if (nbsec(1)  /= nbsec(2) ) write(6,100)'nbsec :', nbsec(1) , nbsec(2) 
      if (mdt(1)    /= mdt(2)   ) write(6,100)'mdt   :', mdt(1)   , mdt(2)   

      write(6,300)
      return

100   format(a,2i9)
300   format(132('*'))
   end subroutine prhddiff

end module printers
