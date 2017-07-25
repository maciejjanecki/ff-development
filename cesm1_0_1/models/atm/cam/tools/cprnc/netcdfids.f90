module netcdfids
!
! Netcdf dimension ids and variable ids.  Initialize to invalid
!
! Dimension ids
!
   integer :: unlimdimid(2)  = -1
!
! Variable ids
!   
   integer :: dateid(2)         = -1
   integer :: datesecid(2)      = -1
   integer :: nstephid(2)       = -1
   integer :: date_writtenid(2) = -1
   integer :: time_writtenid(2) = -1
   integer :: timevid(2)        = -1
   integer :: timedid(2)        = -1
   character(len=*),parameter :: timestr = 'time'

end module netcdfids
