module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use kinds, only : real_kind, iulog
use dimensions_mod, only : nv, np, nlev,qsize
use hybrid_mod, only : hybrid_t
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit
use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, edgevunpackmin, initEdgeBuffer, FreeEdgeBuffer
use bndry_mod, only : bndry_exchangev

implicit none
public :: biharmonic_wk
#ifdef _PRIM
public :: biharmonic_wk_scalar
#endif
public :: compute_zeta_C0
public :: compute_zeta_C0_2d
public :: test_ibyp

contains

#ifdef _PRIM
subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
#else
subroutine biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(nv,nv,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
integer :: nt,nets,nete
#ifdef _PRIM
real (kind=real_kind), dimension(nv,nv,nets:nete) :: pstens
#endif

! local
integer :: k,kptr,i,j,ie,ic
real (kind=real_kind), dimension(:,:), pointer :: spheremv,rspheremv
real (kind=real_kind), dimension(nv,nv) :: lap_ps
real (kind=real_kind), dimension(nv,nv,nlev) :: T
real (kind=real_kind), dimension(nv,nv,2) :: v


   do ie=nets,nete
      
      spheremv     => elem(ie)%spheremv(:,:)

#ifdef _PRIM
      ! should filter lnps + PHI_s/RT?
      pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie))
#endif
      
      do k=1,nlev
         do j=1,nv
         do i=1,nv
#ifdef _PRIM
            T(i,j,k)=elem(ie)%state%T(i,j,k,nt) 
#else            
            ! filter surface height, not thickness
            T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%ps(i,j)
#endif
         enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie))
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

#ifdef _PRIM
      kptr=3*nlev
      call edgeVpack(edge3, pstens(1,1,ie),1,kptr,elem(ie)%desc)
#endif
   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremv(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      
      ! apply inverse mass matrix, then apply laplace again
      do k=1,nlev
         do j=1,nv
         do i=1,nv
            T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
            v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
            v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
         enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie))
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie))
      enddo
         
#ifdef _PRIM
      kptr=3*nlev
      call edgeVunpack(edge3, pstens(1,1,ie), 1, kptr, elem(ie)%desc)
      ! apply inverse mass matrix, then apply laplace again
      lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
      pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie))
#endif

   enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


#ifdef _PRIM


subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nt,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  Q (stored in elem()%state
!    output: qtens  with weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv
integer :: nt,nets,nete

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(:,:), pointer :: spheremv,rspheremv
real (kind=real_kind), dimension(nv,nv) :: lap_p

   do ie=nets,nete
      spheremv     => elem(ie)%spheremv(:,:)
      do q=1,qsize      
      do k=1,nlev
!         qtens(:,:,k,q,ie)=laplace_sphere_wk(elem(ie)%state%Q(:,:,k,q,nt),deriv,elem(ie))
         qtens(:,:,k,q,ie)=laplace_sphere_wk(qtens(:,:,k,q,ie),deriv,elem(ie))
      enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
   enddo
   
   call bndry_exchangeV(hybrid,edgeq)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremv(:,:)
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
      do q=1,qsize      
      do k=1,nlev
         do j=1,nv
         do i=1,nv
            lap_p(i,j)=rspheremv(i,j)*qtens(i,j,k,q,ie)
         enddo
         enddo
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie))
      enddo
      enddo
   enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

#endif





subroutine make_C0_2d(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,nets:nete) :: zeta
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge1


call initEdgeBuffer(edge1,1)

do ie=nets,nete
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%spheremv(:,:)
   kptr=0
   call edgeVpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%rspheremv(:,:)
enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

call FreeEdgeBuffer(edge1) 
end subroutine


subroutine make_C0(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,nlev,nets:nete) :: zeta
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge1


call initEdgeBuffer(edge1,nlev)

do ie=nets,nete
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremv(:,:)
   enddo
   kptr=0
   call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremv(:,:)
   enddo
enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

call FreeEdgeBuffer(edge1) 
end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,2,nlev,nets:nete) :: v
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2


call initEdgeBuffer(edge2,2*nlev)

do ie=nets,nete
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremv(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremv(:,:)
   enddo
   kptr=0
   call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge2)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremv(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremv(:,:)
   enddo
enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

call FreeEdgeBuffer(edge2) 
end subroutine





subroutine compute_zeta_C0_2d(zeta,elem,hybrid,nets,nete,nt,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,nets:nete) :: zeta
integer :: nt,nets,nete,k

! local
integer :: i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo

call make_C0_2d(zeta,elem,hybrid,nets,nete)

end subroutine



subroutine compute_zeta_C0(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(nv,nv,nlev,nets:nete) :: zeta
integer :: nt,nets,nete

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine



subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,nt,min_neigh,max_neigh)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete,nt
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in) :: elem(:)
type (EdgeBuffer_t)  , intent(in) :: edgeMinMax
#ifdef _PRIM
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: Qmin(nv,nv,nlev,qsize)
real (kind=real_kind) :: Qmax(nv,nv,nlev,qsize)
#else
real (kind=real_kind) :: min_neigh(nlev,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,nets:nete)
real (kind=real_kind) :: Qmin(nv,nv,nlev)
real (kind=real_kind) :: Qmax(nv,nv,nlev)
#endif


! local
integer :: ie,k,q

#ifdef _PRIM
    ! compute Qmin, Qmax
    do ie=nets,nete
       do q=1,qsize	
          do k=1,nlev
             Qmin(:,:,k,q)=minval(elem(ie)%state%Q(:,:,k,q,nt))
             Qmax(:,:,k,q)=-maxval(elem(ie)%state%Q(:,:,k,q,nt))
          enddo
       end do
       call edgeVpack(edgeMinMax,Qmax,nlev*qsize,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmin,nlev*qsize,nlev*qsize,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edgeMinMax)
       
    do ie=nets,nete
       do q=1,qsize	
          do k=1,nlev
             Qmin(:,:,k,q)=minval(elem(ie)%state%Q(:,:,k,q,nt))
             Qmax(:,:,k,q)=-maxval(elem(ie)%state%Q(:,:,k,q,nt))
          enddo
       end do
       call edgeVunpackMin(edgeMinMax,Qmax,nlev*qsize,0,elem(ie)%desc)
       call edgeVunpackMin(edgeMinMax,Qmin,nlev*qsize,nlev*qsize,elem(ie)%desc)
       do q=1,qsize
          do k=1,nlev
             max_neigh(k,q,ie)=-min( minval(Qmax(1,:,k,q)), minval(Qmax(nv,:,k,q)),&
                   minval(Qmax(:,1,k,q)),minval(Qmax(:,1,k,q)) )
             min_neigh(k,q,ie)=min(minval(Qmin(1,:,k,q)),minval(Qmin(nv,:,k,q)),&
                  minval(Qmin(:,1,k,q)),minval(Qmin(:,nv,k,q)) )
             min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
          enddo
       end do
    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
#else
    ! compute p min, max
    do ie=nets,nete
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=-maxval(elem(ie)%state%p(:,:,k,nt))
       enddo
       call edgeVpack(edgeMinMax,Qmax,nlev,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmin,nlev,nlev,elem(ie)%desc)
    enddo
    
    call bndry_exchangeV(hybrid,edgeMinMax)
       
    do ie=nets,nete
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=-maxval(elem(ie)%state%p(:,:,k,nt))
       enddo
       call edgeVunpackMin(edgeMinMax,Qmax,nlev,0,elem(ie)%desc)
       call edgeVunpackMin(edgeMinMax,Qmin,nlev,nlev,elem(ie)%desc)
       do k=1,nlev
          max_neigh(k,ie)=-min( minval(Qmax(1,:,k)), minval(Qmax(nv,:,k)),&
               minval(Qmax(:,1,k)),minval(Qmax(:,1,k)) )
          min_neigh(k,ie)=min(minval(Qmin(1,:,k)),minval(Qmin(nv,:,k)),&
               minval(Qmin(:,1,k)),minval(Qmin(:,nv,k)) )
          min_neigh(k,ie) = max(min_neigh(k,ie),0d0)
       enddo
    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
#endif

end subroutine





subroutine local_minmax(elem,hybrid,edgeMinMax,nets,nete,nt,Qmin,Qmax,elem_only)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete,nt
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in) :: elem(:)
type (EdgeBuffer_t)  , intent(in) :: edgeMinMax
logical :: elem_only
#ifdef _PRIM
real (kind=real_kind) :: Qmin(nv,nv,nlev,qsize,nets:nete)
real (kind=real_kind) :: Qmax(nv,nv,nlev,qsize,nets:nete)
#else
real (kind=real_kind) :: Qmin(nv,nv,nlev,nets:nete)
real (kind=real_kind) :: Qmax(nv,nv,nlev,nets:nete)
#endif


! local
integer :: ie,k,q,im,ip,jm,jp,i,j

#ifdef _PRIM
    ! compute Qmin, Qmax
    do ie=nets,nete
       do q=1,qsize	
          do k=1,nlev
             do i=1,nv
             do j=1,nv                
                ip=i+1
                if (i==nv) ip=nv
                im=i-1
                if (i==1) im=1
                jp=j+1
                if (j==nv) jp=nv
                jm=j-1
                if (j==1) jm=1

                Qmin(i,j,k,q,ie)= minval(elem(ie)%state%Q(im:ip,jm:jp,k,q,nt))
                Qmax(i,j,k,q,ie)=-maxval(elem(ie)%state%Q(im:ip,jm:jp,k,q,nt))
             enddo
             enddo
          enddo
       end do
       call edgeVpack(edgeMinMax,Qmax(1,1,1,1,ie),nlev*qsize,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmin(1,1,1,1,ie),nlev*qsize,nlev*qsize,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edgeMinMax)
       
    do ie=nets,nete
       call edgeVunpackMin(edgeMinMax,Qmax(1,1,1,1,ie),nlev*qsize,0,elem(ie)%desc)
       Qmax(:,:,:,:,ie)=-Qmax(:,:,:,:,ie)
       call edgeVunpackMin(edgeMinMax,Qmin(1,1,1,1,ie),nlev*qsize,nlev*qsize,elem(ie)%desc)
    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
#else
    ! compute Qmin, Qmax
    do ie=nets,nete
       do k=1,nlev
          do i=1,nv
             do j=1,nv                
                ip=i+1
                if (i==nv) ip=nv
                im=i-1
                if (i==1) im=1
                jp=j+1
                if (j==nv) jp=nv
                jm=j-1
                if (j==1) jm=1

                Qmin(i,j,k,ie)= minval(elem(ie)%state%p(im:ip,jm:jp,k,nt))
                Qmax(i,j,k,ie)=-maxval(elem(ie)%state%p(im:ip,jm:jp,k,nt))
             enddo
          enddo
       end do
       call edgeVpack(edgeMinMax,Qmax(1,1,1,ie),nlev,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmin(1,1,1,ie),nlev,nlev,elem(ie)%desc)
    enddo
    if (elem_only) then
       Qmax(:,:,:,:)=-Qmax(:,:,:,:)
    else
       call bndry_exchangeV(hybrid,edgeMinMax)
       do ie=nets,nete
          call edgeVunpackMin(edgeMinMax,Qmax(1,1,1,ie),nlev,0,elem(ie)%desc)
          Qmax(:,:,:,ie)=-Qmax(:,:,:,ie)
          call edgeVunpackMin(edgeMinMax,Qmin(1,1,1,ie),nlev,nlev,elem(ie)%desc)
       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    endif
#endif

end subroutine





  subroutine test_ibyp(elem, hybrid,  nets,   nete)

    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use physical_constants, only : rearth 
    ! ---------------------
    use dimensions_mod, only : nv, np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
                               divergence_sphere_wk, curl_sphere
    use global_norms_mod

    implicit none

    type (element_t)     , intent(inout), target :: elem(:)

    type (hybrid_t)      , intent(in) :: hybrid

    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

#if 0
#undef CURLGRAD_TEST
#define IBYP_TEST
    ! =================
    ! Local
    ! =================
    ! pointer ...
    real (kind=real_kind), dimension(:,:), pointer :: rspheremv,spheremv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens2
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens3

    real (kind=real_kind), dimension(nv,nv,2,nets:nete)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(nv,nv,nets:nete)            :: E          ! kinetic energy term
    real (kind=real_kind), dimension(nv,nv,nets:nete)  :: divbig
    real (kind=real_kind), dimension(nv,nv,nets:nete)  :: wdivbig
    real (kind=real_kind), dimension(nv,nv,2,nets:nete)    :: gradbig

    real (kind=real_kind), dimension(nv,nv,2)      :: ulatlon
    real (kind=real_kind), dimension(nv,nv,2)      :: grade
    real (kind=real_kind), dimension(nv,nv,2)      :: grade2
    real (kind=real_kind), dimension(nv,nv)      :: vor
    real (kind=real_kind), dimension(nv,nv)      :: div  
    real (kind=real_kind), dimension(nv,nv)      :: wdiv  

    real (kind=real_kind) ::  v1,v2,v3

    real*8                :: st,et, time_adv
    integer    :: i,j,k,ie,iie,ii,jj
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep

    type (derivative_t)          :: deriv
    call derivinit(deriv)



    ! ===================================
    ! construct test function
    ! ===================================

    do iie=nets,nete
       do jj=1,nv
          do ii=1,nv
             ! test for cardinal function at iie,jj,ii

    write(iulog,'(a,3i4,2e20.10)') 'carinal function:  ',iie,ii,jj

    ! construct two global cardinal functions  pv and E    
    do ie=nets,nete
       do j=1,nv
          do i=1,nv


             E(i,j,ie)=0
#ifdef CURLGRAD_TEST
             if (ie==iie .and. i==ii .and. j==jj) E(i,j,ie)=1
#else
             if (ie==1 .and. i==1 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==1 .and. j==2) E(i,j,ie)=1
             if (ie==1 .and. i==2 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==3 .and. j==3) E(i,j,ie)=1
#endif
             
             ! delta function in contra coordinates
             v1     = 0
             v2     = 0
             if (ie==iie .and. i==ii .and. j==jj) then
                !v1=rearth
                v2=rearth
             endif
             
             ulatlon(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
             ulatlon(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             pv(i,j,1,ie) = ulatlon(i,j,1)
             pv(i,j,2,ie) = ulatlon(i,j,2)
             
          end do
       end do
    enddo
    call make_C0(E,elem,hybrid,nets,nete)
    call make_C0_vector(pv,elem,hybrid,nets,nete) 


#ifdef CURLGRAD_TEST
    ! check curl(grad(E)) 
    do ie=nets,nete
       if ( maxval(abs(E(:,:,ie))) > 0 ) then
       !write(iulog,'(a,i4,2e20.10)') 'maxval: E =',ie,maxval(E(:,:,ie))
       grade=curl_sphere(E(:,:,ie),deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       vor=vorticity_sphere(grade,deriv,elem(ie))
       if (maxval(abs(div))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: div(curl),  vor(curl)=',ie,maxval(abs(div))*rearth**2,maxval(abs(vor))*rearth**2
       endif

       grade=gradient_sphere(E(:,:,ie),deriv,elem(ie))
       vor=vorticity_sphere(grade,deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       if (maxval(abs(vor))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: curl(grad), div(grad)=',ie,maxval(abs(vor))*rearth**2,maxval(abs(div))*rearth**2
       endif
       endif
    enddo

    ! check div(curl(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=curl_sphere(E(:,:,ie),deriv,elem(ie))
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=divergence_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [div([curl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo


    ! check curl(grad(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=gradient_sphere(E(:,:,ie),deriv,elem(ie))
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=vorticity_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [curl([gradl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo
#endif


#ifdef IBYP_TEST
    ! compare <grad(E) dot pv> and <E div(pv)>  < E weak_div(pv) >
    v2=0
    do ie=nets,nete
       spheremv     => elem(ie)%spheremv(:,:)

       div = divergence_sphere(pv(1,1,1,ie),deriv,elem(ie))      ! latlon vector -> scalar 
       grade = gradient_sphere(E(1,1,ie),deriv,elem(ie))
       wdiv = divergence_sphere_wk(pv(1,1,1,ie),deriv,elem(ie)) 



       do j=1,nv
          do i=1,nv
!             write(iulog,'(3i3,3e22.14)') ie,i,j,pv(i,j,1,ie),pv(i,j,2,ie),div(i,j)

             ! (grad(E) dot pv )
             ptens3(i,j,ie) = grade(i,j,1)*pv(i,j,1,ie) + grade(i,j,2)*pv(i,j,2,ie)
             v2=v2+wdiv(i,j)*E(i,j,ie)
          end do
       end do
       ptens(:,:,ie)=div(:,:)*E(:,:,ie)   ! < E div(pv) >
       ptens2(:,:,ie)=wdiv(:,:)*E(:,:,ie)/spheremv(:,:)   ! < wdiv E >
       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       divbig(:,:,ie)=div(:,:)*spheremv(:,:)
       wdivbig(:,:,ie)=wdiv(:,:)
    end do
    call make_C0(divbig,elem,hybrid,nets,nete)
    call make_C0(wdivbig,elem,hybrid,nets,nete)

    
!!$    v1=global_integral(elem,ptens,hybrid,nv,nets,nete)
!!$    v3=global_integral(elem,ptens3,hybrid,nv,nets,nete)
!!$    print *,'< E div(pv) >   =',v1
!!$    print *,'< E div_wk(pv) >=',v2/(4*4*atan(1d0))
!!$    v2=global_integral(elem,ptens2,hybrid,nv,nets,nete)
!!$    print *,'< E div_wk(pv) >=',v2
!!$    print *,'-<grad(E),pv >  =',-v3
!!$!    print *,'sum1-sum2/max',(v1-v2)/max(abs(v1),abs(v2))
!!$


    do ie=nets,nete
       div(:,:)=divbig(:,:,ie)
       wdiv(:,:)=wdivbig(:,:,ie)
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
          do j=1,nv
             do i=1,nv
                ! < div(pv) >   vs < div_wk(pv) >
                if ( abs(div(i,j)-wdiv(i,j)) > .15e-17) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif
                if ( E(i,j,ie)/=0 ) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif
                
             end do
          end do


    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    write(iulog,'(a,3i4,2e20.10)') 'max diff div-wdiv: ',iie,ii,jj,maxval(abs(divbig(:,:,:)-wdivbig(:,:,:))),maxval(divbig(:,:,:))
#endif


    enddo
    enddo
    enddo
    stop
#endif
  end subroutine test_ibyp










end module
