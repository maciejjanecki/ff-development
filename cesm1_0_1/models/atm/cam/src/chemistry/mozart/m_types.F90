
      module m_types
!---------------------------------------------------------------------
! 	... Derived types definition and related parameters
!---------------------------------------------------------------------

      implicit none

!!$      integer, parameter :: hstdim = 200
!!$      type hst_list
!!$	 character(len=32) :: list(hstdim)
!!$      end type hst_list
!!$
!!$      type hst_usr_lst
!!$	 character(len=32) :: name
!!$	 character(len=8)  :: levs
!!$	 character(len=8)  :: units
!!$      end type hst_usr_lst

      type filespec
         character(len=168) :: &
            local_path, &                    ! local file path info
            remote_path, &                   ! remote path info (only used if NCAR is defined)
            nl_filename                      ! filename
         character(len=30) :: &
            hor_res                          ! horizontal resolution
      end type filespec

!!$      type pdiag                             ! routine specific print diagnostics
!!$	 logical :: adv
!!$	 logical :: physlic
!!$	 logical :: imp_slv
!!$	 logical :: negtrc
!!$      end type pdiag

      type time_ramp
	 character(len=8) :: type
	 integer          :: date
	 integer          :: yr_offset
      end type time_ramp

      end module m_types
