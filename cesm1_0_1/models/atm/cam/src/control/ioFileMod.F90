module ioFileMod
!---------------------------------------------------------------------
!
! Purpose:
!
!	Input/Output file manipulations. Mind file on archival system, or local
!	disk etc.
!
! Author: Mariana Vertenstein
!
!---------------------------------------------------------------------
 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils,   only: endrun
   use spmd_utils,   only: masterproc
   use cam_logfile,  only: iulog

   implicit none

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

   private

   public getfil      ! Get file from archive
   public opnfil      ! Open file

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

!=======================================================================
   contains
!=======================================================================
 
   subroutine getfil (fulpath, locfn, iflag)
 
! --------------------------------------------------------------------
! obtain local copy of file
! o first check current working directory
! o next check full pathname[fulpath] on disk
! o by default, abort if file not found.  Optional iflag arg overrides this behavior.
! --------------------------------------------------------------------
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in)  :: fulpath !Archival permanent disk full pathname
   character(len=*), intent(out) :: locfn   !output local file name
   integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer i               !loop index
   integer klen            !length of fulpath character string
   integer maxlen          ! length of locfn input variable
   integer ierr            !error status
   logical lexist          !true if local file exists
   logical abort_on_failure
! --------------------------------------------------------------------
 
   abort_on_failure = .true.
   if (present(iflag)) then
      if (iflag==1) abort_on_failure = .false.
   end if
   maxlen = len(locfn)

! get local file name from full name: start at end. look for first "/"
 
   klen = len_trim(fulpath)
   do i = klen, 1, -1
      if (fulpath(i:i).eq.'/') go to 100
   end do
   i = 0

100 if((klen-i)>maxlen) then
      if(abort_on_failure) then
         call endrun('(GETFIL): local filename variable is too short for path length')
      else
         if(masterproc) write(iulog,*) '(GETFIL): local filename variable is too short for path length',klen-i,maxlen
         RETURN
      end if
   end if

   locfn = fulpath(i+1:klen)
   if (len_trim(locfn) == 0) then
      call endrun ('(GETFIL): local filename has zero length')
   else if(masterproc) then
      write(iulog,*)'(GETFIL): attempting to find local file ', trim(locfn)
   endif
 
! first check if file is in current working directory.
 
   inquire (file=locfn,exist=lexist)
   if (lexist) then
      if(masterproc) write(iulog,*) '(GETFIL): using ',trim(locfn), ' in current working directory'
      RETURN
   endif
 
! second check for full pathname on disk
 
   if(klen>maxlen) then
      if(abort_on_failure) then
         call endrun('(GETFIL): local filename variable is too short for path length')
      else
         if(masterproc) write(iulog,*) '(GETFIL): local filename variable is too short for path length',klen,maxlen
         RETURN
      end if
   end if
   inquire(file=fulpath,exist=lexist)
   if (lexist) then
      locfn = trim(fulpath)
      if(masterproc) write(iulog,*)'(GETFIL): using ',trim(fulpath)
      return
   else
      if(masterproc) write(iulog,*)'(GETFIL): all tries to get file have been unsuccessful: ',trim(fulpath)
      if (abort_on_failure) then
         call endrun ('GETFIL: FAILED to get '//trim(fulpath))
      else
         RETURN
      endif
   endif

   end subroutine getfil
 
!=======================================================================
 
 
   subroutine opnfil (locfn, iun, form, status)
 
!-----------------------------------------------------------------------
! open file locfn in unformatted or formatted form on unit iun
!-----------------------------------------------------------------------
 
! ------------------------ input variables ---------------------------
   character(len=*), intent(in):: locfn  !file name
   integer, intent(in):: iun             !fortran unit number
   character(len=1), intent(in):: form   !file format: u = unformatted. f = formatted
   character(len=*), optional, intent(in):: status !file status
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer ioe             !error return from fortran open
   character(len=11) ft    !format type: formatted. unformatted
   character(len=11) st    !file status: old or unknown
! --------------------------------------------------------------------
 
   if (len_trim(locfn) == 0) then
      call endrun ('(OPNFIL): local filename has zero length')
   endif
   if (form=='u' .or. form=='U') then
      ft = 'unformatted'
   else
      ft = 'formatted  '
   end if
   if ( present(status) ) then
      st = status
   else
      st = "unknown"
   end if
   open (unit=iun,file=locfn,status=st, form=ft,iostat=ioe)
   if (ioe /= 0) then
      if(masterproc) write(iulog,*)'(OPNFIL): failed to open file ',trim(locfn), ' on unit ',iun,' ierr=',ioe
      call endrun ('opnfil') 
   else
      if(masterproc) write(iulog,*)'(OPNFIL): Successfully opened file ',trim(locfn), ' on unit= ',iun
   end if
 
   return
   end subroutine opnfil
 
!=======================================================================
 
 
end module ioFileMod
