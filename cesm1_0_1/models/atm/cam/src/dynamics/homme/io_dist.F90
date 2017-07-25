module io_dist
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use ppgrid,         only: pcols
  use phys_grid, only : get_gcol_all_p, get_ncols_p, gather_chunk_to_field, get_chunk_indices_p
  implicit none
  interface dyn_gather
     module procedure dyn_gather_r4
     module procedure dyn_gather_r8
  end interface

  integer, parameter :: phys_decomp = 1     ! flag indicates physics decomposition
  integer, parameter :: dyn_decomp  = 2     ! flag indicates dynamics decomposition

contains
  subroutine dyn_gather_r4( decomp_type, buff_in, buff_out)
    integer,intent(in) :: decomp_type
    real(r4), intent(in) :: buff_in(:,:,:)
    real(r4), intent(out) :: buff_out(:,:,:)
    integer :: gcols(pcols)
    integer :: lchnk, ncol, i, c, begchunk, endchunk


#ifdef SPMD
#else
    call get_chunk_indices_p(begchunk,endchunk)
    if(decomp_type == phys_decomp) then
      do lchnk=begchunk,endchunk
         c=lchnk-begchunk+1
         ncol = get_ncols_p(lchnk)
         call get_gcol_all_p(lchnk,pcols ,gcols)
         do i=1,ncol
            buff_out(gcols(i),:,1)=buff_in(i,:,c)
         end do
      end do
    else
    end if
#endif
  end subroutine dyn_gather_r4
  subroutine dyn_gather_r8( decomp_type, buff_in, buff_out)
    integer,intent(in) :: decomp_type
    real(r8), intent(in) :: buff_in(:,:,:)
    real(r8), intent(out) :: buff_out(:,:,:)
    integer :: gcols(pcols)
    integer :: c, lchnk, ncol, i, begchunk, endchunk
#ifdef SPMD
    if(decomp_type == phys_decomp) then
       call gather_chunk_to_field (1,size(buff_in,2),1, &
            size(buff_out,1),buff_in,buff_out)
    end if
#else
    call get_chunk_indices_p(begchunk,endchunk)
    if(decomp_type == phys_decomp) then
      do lchnk=begchunk,endchunk
         c=lchnk-begchunk+1
         ncol = get_ncols_p(lchnk)
         call get_gcol_all_p(lchnk,pcols ,gcols)
         do i=1,ncol
            buff_out(gcols(i),:,1)=buff_in(i,:,c)
         end do
      end do
    else
!      do ie=1,nelemd
!         do i=1,elem(ie)%NumUniqueP
!            buff_out(elem(ie)%UniqueOffsetP+i,:,1)=buff_in(



    end if
#endif
  end subroutine dyn_gather_r8

end module io_dist
