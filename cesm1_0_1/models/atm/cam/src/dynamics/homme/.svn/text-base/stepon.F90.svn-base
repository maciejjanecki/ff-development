!BOP
!
! !MODULE: stepon -- FV Dynamics specific time-stepping
!
! !INTERFACE:
module stepon

! !USES:
! from cam
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_sys_mod,    only: shr_sys_flush
   use pmgrid,         only: plev, plevp, plat
   use spmd_utils,     only: iam, masterproc
   use constituents,   only: pcnst, cnst_name, cnst_longname
   use abortutils,     only: endrun
   use ppgrid,         only: begchunk, endchunk
   use physconst,      only: zvir, cappa
   use physics_types,  only: physics_state, physics_tend
   use dyn_comp,       only: dyn_import_t, dyn_export_t
#if defined ( SPMD )
   use mpishorthand,   only: mpicom
#endif
   use perf_mod
! from homme
   use derivative_mod, only : derivinit, deriv_print, derivative_t
   use quadrature_mod, only : gauss, gausslobatto, quadrature_t
   use edge_mod,       only: EdgeBuffer_t, initEdgeBuffer, FreeEdgeBuffer, &
                             edgeVpack, edgeVunpack
 
   implicit none

   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
  public stepon_init   ! Initialization
  public stepon_run1    ! run method phase 1
  public stepon_run2    ! run method phase 2
  public stepon_run3    ! run method phase 3
  public stepon_final  ! Finalization

!----------------------------------------------------------------------
!
! !DESCRIPTION: Module for dynamics time-stepping.
!
! !REVISION HISTORY:
!
! 2006.05.31  JPE    Created
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
  type (derivative_t)   :: deriv           ! derivative struct
  type (quadrature_t)   :: gv,gp           ! quadratures on velocity and pressure grids
  integer :: nets, nete
  type (EdgeBuffer_t) :: edgebuf              ! edge buffer
!-----------------------------------------------------------------------


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_init --- Time stepping initialization
!
! !INTERFACE:
subroutine stepon_init( gw, etamid, dyn_in, dyn_out )
! !USES:
  use dimensions_mod, only: nlev
  use hycoef,             only: hyam, hybm
  use cam_history,      only: phys_decomp, addfld, add_default


! !OUTPUT PARAMETERS
!
  real(r8), intent(out) :: gw(plat)           ! Gaussian weights
  real(r8), intent(out) :: etamid(plev)       ! vertical coords at midpoints
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container

  integer :: m
! !DESCRIPTION:
!
! Allocate data, initialize values, setup grid locations and other
! work needed to prepare the dynamics to run. Return weights and 
! vertical coords to atmosphere calling program.
!
!EOP
!-----------------------------------------------------------------------
!BOC

  call initEdgeBuffer(edgebuf, (3+pcnst)*nlev)

  etamid(:) = hyam(:) + hybm(:)

  ! fields that are written by the dycore
  ! these calls cant be done in dyn_init() because physics grid
  ! is not initialized at that point if making a restart runs
  !
  ! Forcing from physics
  ! FU, FV, other dycores, doc, says "m/s" but I think that is m/s^2
  call addfld ('FU      ','m/s2    ',nlev, 'A','Zonal wind forcing term',phys_decomp)
  call addfld ('FV      ','m/s2    ',nlev, 'A','Meridional wind forcing term',phys_decomp)
  call addfld ('CONVU   ','m/s2    ',nlev, 'A','Zonal component IE->KE conversion term',phys_decomp)
  call addfld ('CONVV   ','m/s2    ',nlev, 'A','Meridional component IE->KE conversion term',phys_decomp)
  call addfld ('DIFFU   ','m/s2    ',nlev, 'A','U horizontal diffusion',phys_decomp)
  call addfld ('DIFFV   ','m/s2    ',nlev, 'A','V horizontal diffusion',phys_decomp)
  
  call addfld ('ETADOT  ','1/s ',plevp,'A','Vertical (eta) velocity',phys_decomp)
  call addfld ('U&IC    ','m/s ',plev, 'I','Zonal wind'        ,phys_decomp )
  call addfld ('V&IC    ','m/s ',plev, 'I','Meridional wind'    ,phys_decomp )
  call add_default ('U&IC       ',0, 'I')
  call add_default ('V&IC       ',0, 'I')

  call addfld ('PS&IC      ','Pa      ',1,    'I','Surface pressure'    ,phys_decomp)
  call addfld ('T&IC       ','K       ',plev, 'I','Temperature'         ,phys_decomp)
  call add_default ('PS&IC      ',0, 'I')
  call add_default ('T&IC       ',0, 'I')
  do m = 1,pcnst
     call addfld (trim(cnst_name(m))//'&IC','kg/kg   ',plev, 'I',cnst_longname(m), phys_decomp)
  end do
  do m = 1,pcnst
     call add_default(trim(cnst_name(m))//'&IC',0, 'I')
  end do
  call addfld ('PHIS&IC      ','m2/s2      ',1,    'I','Surface pressure'    ,phys_decomp)
  call add_default ('PHIS&IC       ',0, 'I')


end subroutine stepon_init

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run1 -- Phase 1 of dynamics run method.
!
! !INTERFACE:
subroutine stepon_run1( dtime_out, phys_state, phys_tend,               &
                        dyn_in, dyn_out )
  use phys_buffer,    only: pbuf
  use dp_coupling, only: d_p_coupling
  use time_mod,    only: tstep          ! dynamics timestep
  use time_manager, only: dtime
  implicit none
!
! !OUTPUT PARAMETERS:
!

   real(r8), intent(out) :: dtime_out   ! Time-step
   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   type(physics_tend), intent(out):: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container

#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
!-----------------------------------------------------------------------

   ! in the case of no subcycling of any kind (tstep=dtime)
   ! we are using leapfrog, so physics dt = 2*tstep
   dtime_out = max(dtime,2*int(tstep))
   if(tstep <= 0)  stop 'bad tstep'
   if(dtime <= 0)  stop 'bad dtime'
   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling (phys_state, phys_tend, pbuf, dyn_out )
   call t_stopf('d_p_coupling')
   
end subroutine stepon_run1

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out )
   use bndry_mod,   only: bndry_exchangeV
   use dimensions_mod, only: nlev, nelemd, nv, np
   use dp_coupling,    only: p_d_coupling
   use parallel_mod,   only: par
   use dyn_comp,       only: TimeLevel
   use phys_buffer,    only: pbuf
   use time_manager,     only: dtime   ! namelist timestep
   use time_mod,        only: tstep,homme_nsplit   !  dynamics typestep
   use control_mod,     only: tracer_advection_formulation,TRACERADV_TOTAL_DIVERGENCE,&
        ftype, qsplit
   use hycoef,             only: hyam, hybm, hyai, hybi, ps0

   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout):: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
   integer :: kptr, ie, ic, i, j, k, tl_f, tl_2
   real(r8) :: rec2dt, dyn_ps0
   real(r8) :: dp(nv,nv,nlev),dp_tmp,fq,fq0,qn0
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
 
   ! copy from phys structures -> dynamics structures
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   call p_d_coupling(phys_state, phys_tend, pbuf, dyn_in)
   call t_stopf('p_d_coupling')

   call t_startf('bndry_exchange')
   ! do boundary exchange
   do ie=1,nelemd
      kptr=0
      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:,1),2*nlev,kptr,dyn_in%elem(ie)%desc)
      kptr=kptr+2*nlev

      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:,1),nlev,kptr,dyn_in%elem(ie)%desc)
      kptr=kptr+nlev
      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FQ(:,:,:,:,1),nlev*pcnst,kptr,dyn_in%elem(ie)%desc)
   end do

   call bndry_exchangeV(par, edgebuf)
   rec2dt = 1./max(dtime,2*int(tstep))  ! see dtime_out above

   do ie=1,nelemd
      kptr=0

      call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:,1),2*nlev,kptr,dyn_in%elem(ie)%desc)
      kptr=kptr+2*nlev

      call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:,1),nlev,kptr,dyn_in%elem(ie)%desc)
      kptr=kptr+nlev

      call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FQ(:,:,:,:,1),nlev*pcnst,kptr,dyn_in%elem(ie)%desc)

      if (qsplit>1) then
         tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics
         tl_2 = TimeLevel%nm1    ! other timelevel to also apply forcing
      else
         tl_f = TimeLevel%nm1   ! timelevel which was adjusted by physics
         tl_2 = TimeLevel%n0    ! other timelevel to also apply forcing
      endif
      dyn_ps0=ps0/100.D0     



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=2:  apply forcing to Q,ps.  Return dynamics tendencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype==2) then
         ! apply forcing to states tl_f 
         if (homme_nsplit==1 .and. qsplit==1) call endrun("ftype==2 and nsplit==1 not allowed")
!$omp parallel do private(k)
         do k=1,nlev
            dp(:,:,k) = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                 ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(:,:,tl_f)
         enddo
         do k=1,nlev
            do j=1,nv
               do i=1,nv

                  do ic=1,pcnst
                     ! back out tendency: Qdp*dtime 
                     fq = dp(i,j,k)*(  dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f))
                     
                     ! apply forcing to Qdp
!                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f) = &
!                          dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f) + fq 
                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f) = &
                          dp(i,j,k)*dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) 

! BEWARE critical region if using OpenMP over k (AAM)
                     if (ic==1) then
                        ! force ps_v to conserve mass:  
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_f)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_f) + fq
                     endif
                  enddo
               end do
            end do
         end do

!$omp parallel do private(k, j, i, ic, dp_tmp)
         do k=1,nlev
          do ic=1,pcnst
            do j=1,nv
               do i=1,nv
                  ! make Q consistent now that we have updated ps_v above
                  ! recompute dp, since ps_v was changed above
                  dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                       ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                  dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f)= &
                       dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f)/dp_tmp

               end do
            end do
          end do
         end do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=1:  apply all forcings as an adjustment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype==1) then
         ! apply forcing to states tl_f and tl_2
         ! when running forward-in-time, we do not need to apply forcings to tl_2
         if (homme_nsplit==1 .and. qsplit==1) call endrun("ftype==1 and nsplit==1 not allowed")
!$omp parallel do private(k)
         do k=1,nlev
            dp(:,:,k) = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                 ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(:,:,tl_f)
         enddo
         do k=1,nlev
            do j=1,nv
               do i=1,nv

                  if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
                  do ic=1,pcnst
                     ! back out tendency: Qdp*dtime 
                     fq = dp(i,j,k)*(  dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f))
                     
                     ! apply forcing to Qdp
                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f) = &
                          dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f) + fq 

                     ! apply forcing to other leapfrog timelevel 
                     fq0=fq
                     qn0=dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_2)
                     ! forcing would drive q negative
                     if (fq0<0 .and. qn0+fq0 < 0) fq0=-qn0
                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_2) = qn0 + fq0

! BEWARE critical region if using OpenMP over k (AAM)
                     if (ic==1 .and. tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
                        ! force ps_v to conserve mass:  
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_f)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_f) + fq
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_2)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_2) + fq0
                     endif
                  enddo
                  else
                  do ic=1,pcnst
                     ! back out tendency: Qdp*dtime 
                     fq = (  dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f))
                     
                     ! apply forcing to Q
                     dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f) = &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f) + fq 

                     fq0=fq
                     qn0=dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_2)
                     ! forcing would drive q negative
                     if (fq0<0 .and. qn0+fq0 < 0) fq0=-qn0
                     dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_2) = qn0 + fq0

#if 1
! BEWARE critical region if using OpenMP over k (AAM)
                     if (ic==1) then
                        ! force ps_v to conserve mass:  
                        fq=fq*dp(i,j,k)
                        fq0=fq0*dp(i,j,k)
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_f)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_f) + fq
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_2)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_2) + fq0
                     endif
#endif
                  enddo
                  endif


                  ! force V, T, both timelevels
                  dyn_in%elem(ie)%state%v(i,j,:,k,tl_f)= &
                       dyn_in%elem(ie)%state%v(i,j,:,k,tl_f) +  &
                       dtime*dyn_in%elem(ie)%derived%FM(i,j,:,k,1)
                  
                  dyn_in%elem(ie)%state%T(i,j,k,tl_f)= &
                       dyn_in%elem(ie)%state%T(i,j,k,tl_f) + &
                       dtime*dyn_in%elem(ie)%derived%FT(i,j,k,1)
                  
                  dyn_in%elem(ie)%state%v(i,j,:,k,tl_2)= &
                       dyn_in%elem(ie)%state%v(i,j,:,k,tl_2) +  &
                       dtime*dyn_in%elem(ie)%derived%FM(i,j,:,k,1)
                  
                  dyn_in%elem(ie)%state%T(i,j,k,tl_2)= &
                       dyn_in%elem(ie)%state%T(i,j,k,tl_2) + &
                       dtime*dyn_in%elem(ie)%derived%FT(i,j,k,1)
               end do
            end do
         end do

         if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
!$omp parallel do private(k, j, i, ic, dp_tmp)
          do k=1,nlev
           do ic=1,pcnst
            do j=1,nv
               do i=1,nv
                  ! make Q consistent now that we have updated ps_v above
                  dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                       ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                  dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f)= &
                       dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_f)/dp_tmp

                  dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                       ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_2)
                  dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_2)= &
                       dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_2)/dp_tmp
               end do
            end do
           end do
          end do
         endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=0 and ftype<0 (debugging options):  just return tendencies to dynamics
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype<=0) then

         do ic=1,pcnst
            if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
               ! Q  =  data used for forcing, at timelevel nm1   t
               ! FQ =  adjusted Q returned by forcing,  at time  t+dt
               ! tendency = (FQ*dp - Q*dp) / dt 
               ! Convert this to a tendency on Qdp:  CAM physics does not change ps
               ! so use ps_v at t.  (or, if physics worked with Qdp
               ! and returned FQdp, tendency = (FQdp-Qdp)/2dt and since physics
               ! did not change dp, dp would be the same in both terms)
!$omp parallel do private(k, j, i, dp_tmp)
               do k=1,nlev
                  do j=1,nv
                     do i=1,nv
                        dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                             ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                        
                        dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1)=&
                             (  dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) - &
                             dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f))*rec2dt*dp_tmp
                     end do
                  end do
               end do
               
            else
               
               ! Q  =  data used for forcing, at timelevel nm1   t-dt
               ! FQ =  adjusted Q returned by forcing,  at time  t+dt
               ! tendency = (FQ - Q) / 2dt 
!$omp parallel do private(k, j, i)
               do k=1,nlev
                  do j=1,nv
                     do i=1,nv
                        dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1)=&
                             (  dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1) - &
                             dyn_in%elem(ie)%state%Q(i,j,k,ic,tl_f))*rec2dt
                     end do
                  end do
               end do
            end if
         end do
      endif
   end do
   call t_stopf('bndry_exchange')
	

end subroutine stepon_run2


subroutine stepon_run3( dtime, etamid, cam_out, phys_state, dyn_in, dyn_out )
   use camsrfexch_types, only: cam_out_t     
   use dyn_comp,    only: dyn_run
   use time_mod,    only: tstep
   real(r8), intent(in) :: dtime   ! Time-step
   real(r8), intent(in)  :: etamid(plev)
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(out)   :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out)   :: dyn_out ! Dynamics export container
   integer :: rc
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(dyn_out,rc)	
   call t_stopf  ('dyn_run')

end subroutine stepon_run3


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)

! !PARAMETERS:
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC


!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------


end module stepon
