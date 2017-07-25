    module modal_aer_opt

#ifdef MODAL_AERO
!     parameterizes aerosol coefficients using chebychev polynomial
!     parameterize aerosol radiative properties in terms of
!     surface mode wet radius and wet refractive index

!     Ghan and Zaveri, JGR 2007.

!     uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pcols, pver, pverp, begchunk, endchunk
use abortutils,        only: endrun
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use error_messages,    only: handle_ncerr
use modal_aero_data,   only: maxd_amode, maxd_aspectype, ntot_amode, ntot_aspectype, &
                             nspec_amode, lspectype_amode, lmassptr_amode, numptr_amode, alnsg_amode, &
                             specrefndxsw, specrefndxlw, specdens_amode,   &
!++xl added 11/24/2009, dust AOD
                             lptr_dust_a_amode, specdens_dust_amode, &
!--dust AOD
		             sigmag_amode, species_class, spec_class_aerosol
   use physconst,      only: rhoh2o
   use cam_history,    only: phys_decomp, addfld, add_default
   use perf_mod
   use radconstants,   only: nswbands, nlwbands
   use cam_logfile,    only: iulog

implicit none
private
save
      public :: modal_aer_opt_init, modal_size_parameters, modal_aero_sw, modal_aero_lw


      integer nrefr,nrefi,nr,ni
      integer, parameter :: prefr=7,prefi=10

      integer, parameter, public :: ncoef=5

!     coefficients for parameterizing aerosol radiative properties
!     in terms of refractive index and wet radius
      real(r8) extpsw(ncoef,prefr,prefi,maxd_amode,nswbands) ! specific extinction
      real(r8) abspsw(ncoef,prefr,prefi,maxd_amode,nswbands) ! specific absorption
      real(r8) asmpsw(ncoef,prefr,prefi,maxd_amode,nswbands) ! asymmetry factor
      real(r8) extplw(ncoef,prefr,prefi,maxd_amode,nlwbands) ! specific extinction
      real(r8) absplw(ncoef,prefr,prefi,maxd_amode,nlwbands) ! specific absorption
      real(r8) asmplw(ncoef,prefr,prefi,maxd_amode,nlwbands) ! asymmetry factor
#ifdef AEROCOM
      logical saveaerocom
      complex crefw_aerocom(naerocom)
      real(r8) extp_aerocom(ncoef,prefr,prefi,maxd_amode,naerocom)
      real(r8) absp_aerocom(ncoef,prefr,prefi,maxd_amode,naerocom)
      real(r8) asmp_aerocom(ncoef,prefr,prefi,maxd_amode,naerocom)
      real(r8) refrtab_aerocom(prefr,naerocom),refitab_aerocom(prefi,naerocom)
#endif
      integer, parameter :: pnangle=1

!     specabs = absorption coeff / unit dry mass
!     specscat = scattering coeff / unit dry mass
      complex crefin
      complex, public :: crefwsw(nswbands) ! complex refractive index for water visible
      complex, public :: crefwlw(nlwbands) ! complex refractive index for water infrared
      real(r8), pointer :: refrwsw(:),refiwsw(:) ! real, imaginary ref index for water visible
      real(r8), pointer :: refrwlw(:),refiwlw(:) ! real, imaginary ref index for water infrared
!      data refrwsw/19*1.33/
!      data refiwsw/8*1.e-9,1.5e-8,7*1.e-5,6.e-2,2*1.e-2/
      real(r8) pie

      real(r8) rmmin,rmmax ! min, max aerosol surface mode radius treated (m)
      real(r8) xrmin,xrmax

      real(r8) refrmin ! minimum of real part of refractive index
      real(r8) refrmax ! maximum of real part of refractive index
      real(r8) refimin ! minimum of imag part of refractive index
      real(r8) refimax ! maximum of imag part of refractive index
      real(r8) refrtabsw(prefr,nswbands) ! table of real refractive indices for aerosols visible
      real(r8) refitabsw(prefi,nswbands) ! table of imag refractive indices for aerosols visible
      real(r8) refrtablw(prefr,nlwbands) ! table of real refractive indices for aerosols infrared
      real(r8) refitablw(prefi,nlwbands) ! table of imag refractive indices for aerosols infrared
#if (defined BACKSCAT)
      real(r8)  angmin, angmax=3.14 ! minimum, maximum scattering angle
      real(r8)  ang(pnangle), xmu(pnangle)
#endif
!===============================================================================
CONTAINS
!===============================================================================
   subroutine modal_aer_opt_init(wavminsw,wavmaxsw,wavenumber1_longwave,wavenumber2_longwave)

  use rad_constituents, only: rad_cnst_get_clim_aer_props
   use ioFileMod,    only: getfil
   use filenames,    only: modal_optics
#if ( defined SPMD )
   use mpishorthand,   only: mpicom, mpiint, mpir8, mpichar
#endif

implicit none

      real(r8), intent(in):: wavmaxsw(nswbands)
      real(r8), intent(in):: wavminsw(nswbands)
      real(r8), intent(in):: wavenumber1_longwave(nlwbands) ! (cm-1)
      real(r8), intent(in):: wavenumber2_longwave(nlwbands) ! (cm-1)
      real(r8) wavmidsw(nswbands) ! middle wavelength of shortwave band (m)
      real(r8) wavmidlw(nlwbands) ! middle wavelength of longwave band (m)
      real(r8) refr,refi
      integer lw_band_len
      integer  :: nc_id ! index to netcdf file
      character*256 outfilename,infilename
      logical, parameter :: calc_optics = .false.
      integer :: ierr ! error codes from pio


      integer ns, m, l, ltype

      call rad_cnst_get_clim_aer_props(1, &
                            refindex_real_water_sw=refrwsw,  &
                            refindex_im_water_sw=refiwsw,  &
			    refindex_real_water_lw=refrwlw,  &
			    refindex_im_water_lw=refiwlw)

      pie=4._r8*atan(1._r8)

      nrefr=prefr
      nrefi=prefi
      rmmin=0.01e-6
      rmmax=25.e-6
      xrmin=log(rmmin)
      xrmax=log(rmmax)

#if (defined BACKSCAT)
      angmin=1.1 ! minimum scattering angle between satellite and sun
      angmax=3.14 ! backscattering
      ang(1)=angmax ! backscatter angle only
      do i=1,pnangle
         xmu(i)=cos(ang(i))
      enddo
#endif

      if(calc_optics)then

!        wavelength loop

         do 200 ns=1,nswbands

!           wavelength (m)

            wavmidsw(ns) = 0.5e-6*(wavminsw(ns) + wavmaxsw(ns))
!            write(iulog,*)'wavelength=',wavmidsw(ns)

!           first find min,max of real and imaginary parts of refractive index

!           real and imaginary parts of water refractive index

            crefwsw(ns)=cmplx(refrwsw(ns),refiwsw(ns))

            refrmin=real(crefwsw(ns))
            refrmax=real(crefwsw(ns))
            refimin=aimag(crefwsw(ns))
            refimax=aimag(crefwsw(ns))

!           aerosol mode loop

            do m=1,ntot_amode

!              aerosol species loop

               do l=1,nspec_amode(m)

                  ltype=lspectype_amode(l,m)

!                 real and imaginary parts of aerosol refractive index

                  refr=real(specrefndxsw(ns,ltype))
                  refi=aimag(specrefndxsw(ns,ltype))

                  refrmin=min(refrmin,refr)
                  refrmax=max(refrmax,refr)
                  refimin=min(refimin,refi)
                  refimax=max(refimax,refi)

               enddo
            enddo

!	    write(iulog,*)'sw band ',ns
!            write(iulog,*)'refrmin=',refrmin
!            write(iulog,*)'refrmax=',refrmax
!            write(iulog,*)'refimin=',refimin
!            write(iulog,*)'refimax=',refimax
            if(refimin<-1.e-3)then
               write(iulog,*)'negative refractive indices not allowed'
               call endrun('error in modal_aer_opt_init')
            endif

	    call miefit(wavmidsw(ns),refrtabsw(1,ns),refitabsw(1,ns),   &
	       extpsw(1,1,1,1,ns),abspsw(1,1,1,1,ns),asmpsw(1,1,1,1,ns) )
  200    continue

         do 300 ns=1,nlwbands

!            wavmidir=10.e-6 ! wavelength (m) at 10 microns
            wavmidlw(ns) = 0.5e-2*(1._r8/wavenumber1_longwave(ns) + 1._r8/wavenumber1_longwave(ns))

!           first find min,max of real and imaginary parts of infrared (10 micron) refractive index

!           real and imaginary parts of water refractive index

!            refrwir=1.179 ! Kent et al., Applied Optics, 22, 1655-1665, 1983
!	    refiwir=0.6777 ! Kent et al., Applied Optics, 22, 1655-1665, 1983
            crefwlw(ns)=cmplx(refrwlw(ns),refiwlw(ns))

            refrmin=real(crefwlw(ns))
            refrmax=real(crefwlw(ns))
            refimin=aimag(crefwlw(ns))
            refimax=aimag(crefwlw(ns))

!           aerosol mode loop

            do m=1,ntot_amode

!              aerosol species loop

               do l=1,nspec_amode(m)

                  ltype=lspectype_amode(l,m)

!                 real and imaginary parts of aerosol refractive index

		  refr=real(specrefndxlw(ns,ltype))
		  refi=aimag(specrefndxlw(ns,ltype))
                  refrmin=min(refrmin,refr)
                  refrmax=max(refrmax,refr)
                  refimin=min(refimin,refi)
                  refimax=max(refimax,refi)

               enddo
            enddo

!	    write(iulog,*)'longwave band ',ns
!            write(iulog,*)'refrmin=',refrmin
!            write(iulog,*)'refrmax=',refrmax
!            write(iulog,*)'refimin=',refimin
!            write(iulog,*)'refimax=',refimax

	    call miefit(wavmidlw(ns),refrtablw(1,ns),refitablw(1,ns), &
	        extplw(1,1,1,1,ns),absplw(1,1,1,1,ns),asmplw(1,1,1,1,ns) )
  300    continue

            outfilename='modal_optics.nc'

            call write_modal_optics(outfilename,  maxd_amode, nlwbands, nswbands,  &
                 refrtabsw, refitabsw, refrtablw, refitablw, &
                 prefr, prefi, sigmag_amode, ncoef, extpsw, abspsw, asmpsw, absplw )

        
     else


            infilename='modal_optics.nc'
            call getfil(modal_optics, infilename)
	    write(iulog,*)'modal_optics filename=',infilename
            call read_modal_optics(infilename,  maxd_amode, nlwbands, nswbands,  &
	       refrtabsw, refitabsw, refrtablw, refitablw, &
	       prefr, prefi, ncoef, extpsw, abspsw, asmpsw, absplw )

	 do ns=1,nswbands

            wavmidsw(ns) = 0.5e-6*(wavminsw(ns) + wavmaxsw(ns))

            crefwsw(ns)=cmplx(refrwsw(ns),refiwsw(ns))

         end do

	 do ns=1,nlwbands
!            wavmidir=10.e-6 ! wavelength (m) at 10 microns
!            wavmidlw(ns) = 0.5e-6*(wavminlw(ns) + wavmaxlw(ns))
            wavmidlw(ns) = 0.5e-2*(1._r8/wavenumber1_longwave(ns) + 1._r8/wavenumber1_longwave(ns))

!           real and imaginary parts of water refractive index

           crefwlw(ns)=cmplx(refrwlw(ns),refiwlw(ns))
	 end do

    endif

     call addfld ('EXTINCT','/m  ',pver,    'A','Aerosol extinction',phys_decomp, flag_xyfill=.true.)
     call addfld ('ABSORB','/m  ',pver,    'A','Aerosol absorption',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODVIS','  ',1,    'A','Aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODABS','  ',1,    'A','Aerosol absorption optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE1','  ',1,    'A','Aerosol optical depth 550 nm mode 1',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE2','  ',1,    'A','Aerosol optical depth 550 nm mode 2',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE3','  ',1,    'A','Aerosol optical depth 550 nm mode 3',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE4','  ',1,    'A','Aerosol optical depth 550 nm mode 4',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE5','  ',1,    'A','Aerosol optical depth 550 nm mode 5',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE6','  ',1,    'A','Aerosol optical depth 550 nm mode 6',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODMODE7','  ',1,    'A','Aerosol optical depth 550 nm mode 7',phys_decomp, flag_xyfill=.true.)
!++xl added 11/24/2009, dust AOD
     call addfld ('AODDUST1','  ',1,    'A','Aerosol optical depth 550 nm model 1 from dust',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODDUST2','  ',1,    'A','Aerosol optical depth 550 nm model 2 from dust',phys_decomp, flag_xyfill=.true.)
     call addfld ('AODDUST3','  ',1,    'A','Aerosol optical depth 550 nm model 3 from dust',phys_decomp, flag_xyfill=.true.)
!--dust AOD
     call addfld ('BURDEN1','kg/m2',1,    'A','Aerosol burden mode 1',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN2','kg/m2',1,    'A','Aerosol burden mode 2',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN3','kg/m2',1,    'A','Aerosol burden mode 3',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN4','kg/m2',1,    'A','Aerosol burden mode 4',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN5','kg/m2',1,    'A','Aerosol burden mode 5',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN6','kg/m2',1,    'A','Aerosol burden mode 6',phys_decomp, flag_xyfill=.true.)
     call addfld ('BURDEN7','kg/m2',1,    'A','Aerosol burden mode 7',phys_decomp, flag_xyfill=.true.)
     call addfld ('SSAVIS','  ',1,    'A','Aerosol singel-scatter albedo',phys_decomp, flag_xyfill=.true.)
     call add_default ('EXTINCT', 1, ' ')
     call add_default ('ABSORB', 1, ' ')
     call add_default ('AODVIS', 1, ' ')
     call add_default ('AODABS', 1, ' ')
     call add_default ('AODMODE1', 1, ' ')
     call add_default ('AODMODE2', 1, ' ')
     call add_default ('AODMODE3', 1, ' ')
#if (defined MODAL_AERO_7MODE)
     call add_default ('AODMODE4', 1, ' ')
     call add_default ('AODMODE5', 1, ' ')
     call add_default ('AODMODE6', 1, ' ')
     call add_default ('AODMODE7', 1, ' ')
#endif
!++xl added 11/24/2009, dust AOD
     call add_default ('AODDUST1', 1, ' ')
     call add_default ('AODDUST2', 1, ' ')
     call add_default ('AODDUST3', 1, ' ')
!--dust AOD
     call add_default ('BURDEN1', 1, ' ')
     call add_default ('BURDEN2', 1, ' ')
     call add_default ('BURDEN3', 1, ' ')
#if (defined MODAL_AERO_7MODE)
     call add_default ('BURDEN4', 1, ' ')
     call add_default ('BURDEN5', 1, ' ')
     call add_default ('BURDEN6', 1, ' ')
     call add_default ('BURDEN7', 1, ' ')
#endif
     call add_default ('SSAVIS', 1, ' ')


   end subroutine modal_aer_opt_init

   subroutine modal_size_parameters(ncol, dgnumwet, radsurf, logradsurf, cheb)

   use phys_buffer,      only: pbuf_size_max, pbuf_fld,  pbuf_get_fld_idx

   implicit none

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dgnumwet(pcols,pver,maxd_amode) ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(pcols,pver,maxd_amode) ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(pcols,pver,maxd_amode) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(ncoef,maxd_amode,pcols,pver)

   integer m,l,lmass,lnum,i,k,nc,ifld,lwater
   real(r8) explnsigma
   real(r8) xrad(pcols) ! normalized aerosol radius

   do m=1,ntot_amode

      explnsigma=exp(2.0_r8*alnsg_amode(m)*alnsg_amode(m))

      do k=1,pver
      do i=1,ncol
!        convert from number mode diameter to surface area
         radsurf(i,k,m)=0.5*dgnumwet(i,k,m)*explnsigma
         logradsurf(i,k,m)=log(radsurf(i,k,m))
!        normalize size parameter
         xrad(i)=max(logradsurf(i,k,m),xrmin)
         xrad(i)=min(xrad(i),xrmax)
         xrad(i)=(2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
!        chebyshev polynomials
         cheb(1,m,i,k)=1._r8
         cheb(2,m,i,k)=xrad(i)
         do nc=3,ncoef
            cheb(nc,m,i,k)=2._r8*xrad(i)*cheb(nc-1,m,i,k)-cheb(nc-2,m,i,k)
         enddo
      end do
      end do

   end do

   end subroutine modal_size_parameters

      subroutine miefit(wavmid,refrtab,refitab, extparam,absparam,asmparam)

   implicit none

      integer nsiz ! number of wet particle sizes
      integer nlog ! number of log-normals modes to product mie fit
      parameter (nsiz=200,nlog=30)

      real(r8) rad(nsiz)
      real(r8) wavmid
      real(r8) extparam(ncoef,prefr,prefi,maxd_amode)
      real(r8) absparam(ncoef,prefr,prefi,maxd_amode)
      real(r8) asmparam(ncoef,prefr,prefi,maxd_amode)
      real(r8) refrtab(prefr),refitab(prefi)
      real(r8) size ! 2 pi radpart / waveleng = size parameter
      real(r8) qext(nsiz) ! array of extinction efficiencies
      real(r8) qsca(nsiz) ! array of scattering efficiencies
      real(r8) qabs(nsiz) ! array of absorption efficiencies
      real(r8) gqsc(nsiz) ! array of asymmetry factor * scattering efficiency
      real(r8) asymm(nsiz) ! array of asymmetry factor
      complex*16 sforw,sback,tforw(2),tback(2),refindx
      integer :: nmom,ipolzn,momdim,numang
      real(r8) pmom(0:1,1)
      logical :: perfct,prnt(2)
      complex*16 crefd
      real(r8) :: mimcut
      logical :: anyang

      real(r8) xmu
      complex*16 s1(1),s2(1)

      real(r8) dsdlogr(nsiz),dlogr
      real(r8) rmin,rmax ! min, max aerosol size bin
      real(r8) xr
      real(r8) drefr ! increment in real part of refractive index
      real(r8) drefi ! increment in imag part of refractive index
      integer m,n,nr,ni,nl

      real(r8) exparg
      real(r8) volwet ! sum of volume of wet aerosols
      real(r8) sumabs ! sum of specific absorption coefficients
      real(r8) sumsca ! sum of specific scattering coefficients
      real(r8) sumg   ! sum of asymmetry factors
      real(r8) rs(nlog) ! surface mode radius (cm)
      real(r8) specscat(nlog) ! specific scattering (m2/g)
      real(r8) specabs(nlog)  ! specific absorption (m2/g)
      real(r8) specext(nlog)  ! specfiic extinction (m2/g)
      real(r8) abs(nlog)      ! single scattering albedo
      real(r8) asym(nlog)     ! asymmetry factor
      real(r8) logr(nsiz)
      real(r8) bma,bpa
      real(r8) sq2pi

      sq2pi=2.5066283

      nmom=0
      ipolzn=0
      momdim=1
      perfct=.false.
      prnt(1)=.false.
      prnt(2)=.false.
      mimcut=0.0_r8

      anyang=.false.
      numang=0


        drefr=(refrmax-refrmin)/(nrefr-1)

!       size bins (m)

        rmin=0.001e-6
        rmax=100.e-6

        dlogr=log(rmax/rmin)/(nsiz-1)
        logr(1)=log(rmin)
        do n=2,nsiz
           logr(n)=logr(n-1)+dlogr
        enddo

        do n=1,nsiz
           rad(n)=exp(logr(n))
        enddo

!           calibrate parameterization with range of refractive indices

            do 120 nr=1,nrefr
            do 120 ni=1,nrefi

               refrtab(nr)=refrmin+(nr-1)*drefr
               refitab(ni)=refimax*0.3**(nrefi-ni)
	       if(ni.eq.1)refitab(ni)=0.
               crefd=dcmplx(refrtab(nr),refitab(ni))
!               write(iulog,*)'nr,ni,refr,refi=',nr,ni,refrtab(nr),refitab(ni)

!              mie calculations of optical efficiencies

               do n=1,nsiz

!                 size parameter and weighted refractive index

                  size=2.*pie*rad(n)/wavmid
                  size=min(size,400.d0)

                  call miev0(size,crefd,perfct,mimcut,anyang,          &
                     numang,xmu,nmom,ipolzn,momdim,prnt,               &
                     qext(n),qsca(n),gqsc(n),pmom,sforw,sback,s1,      &
                     s2,tforw,tback )

                  qsca(n)=min(qsca(n),qext(n))
                  qabs(n)=qext(n)-qsca(n)
                  asymm(n)=gqsc(n)/qsca(n)
!                 if(nr.eq.1.and.ni.eq.1)then
!                  write(iulog,*)'rad=',rad(n),' qsca,qabs=',qsca(n),qabs(n)
!                endif


               enddo

!            now consider a variety of lognormal size distributions
!            with surface mode radius rs (m) and log standard deviation alnsg_amode

!            aerosol mode loop

             do 110 m=1,ntot_amode

!              size units are m

!              bounds of surface mode radius

               bma=0.5*log(rmmax/rmmin)
               bpa=0.5*log(rmmax*rmmin)

               do 100 nl=1,nlog

!                 pick rs to match zeroes of chebychev

                  xr=cos(pie*(nl-0.5)/nlog)
                  rs(nl)=exp(xr*bma+bpa)

!                 define wet size distribution

!                 integrate over all size bins

                  volwet=0
                  sumsca=0
                  sumabs=0
                  sumg=0

                  do n=1,nsiz
                     exparg=log(rad(n)/rs(nl))/alnsg_amode(m)
                     dsdlogr(n)=exp(-0.5*exparg*exparg)
                     volwet=volwet+4./3.*rad(n)*dsdlogr(n)*dlogr
                     sumabs=sumabs+qabs(n)*dsdlogr(n)*dlogr
                     sumsca=sumsca+qsca(n)*dsdlogr(n)*dlogr
                     sumg=sumg+asymm(n)*qsca(n)*dsdlogr(n)*dlogr
                  enddo


!                 coefficients wrt wet mass, using water density
                  specscat(nl)=sumsca/(volwet*rhoh2o)
                  specabs(nl)=sumabs/(volwet*rhoh2o)
!                 coefficients wrt wet mass
                  specext(nl)=specscat(nl)+specabs(nl)
                  asym(nl)=sumg/sumsca

                 !if(nr.eq.1.and.ni.eq.1)then
                  if(m.eq.3.and.wavmid.gt.0.45e-6.and.wavmid.lt.0.6e-6)then
!                    write(98,*)'rs=',rs(nl),' cref=',refrtab(nr),refitab(ni),' specext=',specext(nl)
                  endif
                ! endif

  100          continue

               call fitcurv(rs,specext,extparam(1,nr,ni,m),ncoef,nlog)
               call fitcurvlin(rs,specabs,absparam(1,nr,ni,m),ncoef,nlog)
               call fitcurvlin(rs,asym,asmparam(1,nr,ni,m),ncoef,nlog)

  110         continue
  120       continue
!       write(iulog,*)'refitab=',(refitab(ni),ni=1,nrefi)
      return
      end subroutine miefit

      subroutine modal_aero_sw(ncol,lchnk,iswband, nnite, idxnite, tauxar,wa,ga,fa, raer, state, &
                               qaerwat, radsurf, logradsurf, cheb )

!  calculates aerosol radiative properties

   use ppgrid
   use physconst,        only: rga, rair
   use physics_types,    only: physics_state
   use radconstants,     only: idx_sw_diag
   use cam_history,      only: outfld, fillvalue

   implicit none

   integer, intent(in) :: ncol  ! number of columns
   integer, intent(in) :: lchnk  ! chunk index
   integer, intent(in) :: iswband  ! spectral interval
   integer, intent(in) :: nnite          ! number of night columns
   integer, intent(in) :: idxnite(nnite) ! local column indices of night columns
   real(r8), intent(in) :: raer(pcols, pver, pcnst)
   type(physics_state), intent(in) :: state
   real(r8), intent(in) :: qaerwat(pcols,pver,maxd_amode) ! aerosol water (g/g)
   real(r8), intent(in) :: radsurf(pcols,pver,maxd_amode) ! aerosol surface mode radius
   real(r8), intent(in) :: logradsurf(pcols,pver,maxd_amode) ! log(aerosol surface mode radius)
   real(r8), intent(in) :: cheb(ncoef,maxd_amode,pcols,pver)

         real(r8), intent(out) :: tauxar(pcols,0:pver) ! layer extinction optical depth
	 real(r8), intent(out) :: wa(pcols,0:pver) ! layer single-scatter albedo
	 real(r8), intent(out) :: ga(pcols,0:pver) ! asymmetry factor
	 real(r8), intent(out) :: fa(pcols,0:pver) ! forward scattered fraction

         integer i,k
         real(r8) pext(pcols) ! parameterized specific extinction (m2/kg)
         real(r8) specpext(pcols) ! specific extinction (m2/kg)
         real(r8) pabs(pcols) ! parameterized specific absorption (m2/kg)
         real(r8) palb(pcols) ! parameterized single scattering albedo
         real(r8) pasm(pcols) ! parameterized asymmetry factor
      real(r8) :: mass(pcols,pver) ! layer mass
      real(r8) totmass(pcols)
      real(r8) vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
      real(r8) volsol(pcols)   ! volume concentration of soluble aerosol specie (m3/kg)
      real(r8) volinsol(pcols) ! volume concentration of insoluble aerosol specie (m3/kg)
      real(r8) dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
      real(r8) dryvolsol(pcols)   ! volume concentration of soluble aerosol mode (m3/kg)
      real(r8) dryvolinsol(pcols) ! volume concentration of insoluble aerosol mode (m3/kg)
      real(r8) wetvol(pcols)   ! volume concentration of wet mode (m3/kg)
      real(r8) watervol(pcols) ! volume concentration of water in each mode (m3/kg)
      real(r8) volf            ! volume fraction of insoluble aerosol
      real(r8) refr(pcols)     ! real part of refractive index
      real(r8) refi(pcols)     ! imaginary part of refractive index
      real(r8) :: aodvis(pcols) ! extinction optical depth
      real(r8) :: aodabs(pcols) ! absorption optical depth
      real(r8) :: aodmode(pcols,maxd_amode)
      real(r8) :: burden(pcols,maxd_amode)
!++xl added 11/24/2009, dust AOD
      real(r8) :: dustvol(pcols)  ! volume concentration of dust in aerosol mode (m3/kg)
      real(r8) :: dustaodmode(pcols,maxd_amode)  ! dust aod in aerosol mode
!--dust AOD
      real(r8) :: ssavis(pcols)
      real(r8) :: extinct(pcols,pver)
      real(r8) :: absorb(pcols,pver)
   logical :: savaervis ! true if visible wavelength (0.55 micron)

      complex crefin(pcols), crefsol(pcols), crefinsol(pcols) ! complex refractive index
      complex crefinsq, crefsolsq, crefinsolsq
      complex crefa, crefb
      integer lmass,ltype
      integer itab(pcols),jtab(pcols)
      real(r8) ttab(pcols),utab(pcols)
      real(r8) cext(pcols,ncoef),cabs(pcols,ncoef),casm(pcols,ncoef)
      real(r8)  air_density(pcols,pver) ! (kg/m3)
      integer, parameter :: nerrmax_dopaer=1000
      integer nerr_dopaer
      save nerr_dopaer
      data nerr_dopaer/0/
      real(r8) dopaer(pcols) ! aerosol optical depth in layer
      real(r8) third
      real(r8) colext(pcols,maxd_amode)
	 integer m,l,nc,ns
	
                do i=1,ncol
                   extinct(i,:)=0.0_r8
                   absorb(i,:)=0.0_r8
                   aodvis(i)=0.0_r8
                   aodmode(i,:)=0.0_r8
                   aodabs(i)=0.0_r8
!++xl added 11/24/2009, dust AOD
                   dustaodmode(i,:)=0.0_r8
!--dust AOD
                   ssavis(i)=0.0_r8
		end do

      do k=1,pver
         tauxar(:ncol,k)=0._r8
         wa(:ncol,k)=0._r8
         ga(:ncol,k)=0._r8
	 fa(:ncol,k)=0._r8
	 mass(:ncol,k)=state%pdeldry(:ncol,k)*rga
	 air_density(:ncol,k)=state%pmid(:ncol,k)/(rair*state%t(:ncol,k))
      enddo
      burden(:ncol,:ntot_amode)=0._r8
      third=1./3.
      savaervis=iswband.eq.idx_sw_diag

! zero'th layer does not contain aerosol
!
   tauxar(1:ncol,0)  = 0._r8
   wa(1:ncol,0)      = 0.925_r8
   ga(1:ncol,0)      = 0.850_r8
   fa(1:ncol,0)      = 0.7225_r8

   do k=1,pver

!     loop over all aerosol modes

      do m=1,ntot_amode

!        form bulk refractive index
         crefin(:ncol)=0._r8
         dryvol(:ncol)=0._r8
!++xl added 11/24/2009, dust AOD
         dustvol(:ncol)=0._r8
!--dust AOD

!        aerosol species loop

         do l=1,nspec_amode(m)
            ltype=lspectype_amode(l,m)
            lmass=lmassptr_amode(l,m)
	    do i=1,ncol
	       burden(i,m)=burden(i,m) + raer(i,k,lmass)*mass(i,k)
               vol(i)=raer(i,k,lmass)/specdens_amode(ltype)
               dryvol(i)=dryvol(i)+vol(i)
               crefin(i)=crefin(i)+vol(i)*specrefndxsw(iswband,ltype)
	    end do
         enddo
!++xl added 11/24/2009, dust AOD
         do i=1,ncol
           if (lptr_dust_a_amode(m)>0) then
             dustvol(i)=raer(i,k,lptr_dust_a_amode(m))/specdens_dust_amode
           else
             dustvol(i)=0.0_r8
           endif
         enddo
!--dust AOD
         do i=1,ncol
       	    watervol(i)=qaerwat(i,k,m)/rhoh2o
	    wetvol(i)=watervol(i)+dryvol(i)
            if(watervol(i)<0._r8)then
               if(abs(watervol(i)).gt.1.e-1*wetvol(i))then
                  write(iulog,'(a,4e10.2,a)')'watervol,wetvol,dryvolsol,dryvolinsol=', &
		                   watervol(i),wetvol(i),dryvolsol(i),dryvolinsol(i) &
                         ,' in modal_aero_sw'
!               !  call endrun()
               endif
                 watervol(i)=0._r8
                 wetvol(i)=dryvol(i)
            endif
! volume mixing
            crefin(i)=crefin(i)+watervol(i)*crefwsw(iswband)
            crefin(i)=crefin(i)/max(wetvol(i),1.e-60_r8)
            refr(i)=real(crefin(i))
            refi(i)=abs(aimag(crefin(i)))
          end do ! ncol

!        call t_startf('binterp')

!        interpolate coefficients linear in refractive index
!        first call calcs itab,jtab,ttab,utab
         itab(:ncol)=0
         call binterp(extpsw(1,1,1,m,iswband),ncol,ncoef,nrefr,nrefi,   &
                          refr,refi,refrtabsw(1,iswband),refitabsw(1,iswband),itab,jtab,  &
                          ttab,utab,cext)

         call binterp(abspsw(1,1,1,m,iswband),ncol,ncoef,nrefr,nrefi,   &
                          refr,refi,refrtabsw(1,iswband),refitabsw(1,iswband),itab,jtab,  &
                          ttab,utab,cabs)
         call binterp(asmpsw(1,1,1,m,iswband),ncol,ncoef,nrefr,nrefi,   &
                          refr,refi,refrtabsw(1,iswband),refitabsw(1,iswband),itab,jtab,  &
                          ttab,utab,casm)

!       call t_stopf('binterp')

!        parameterized optical properties
         do i=1,ncol
            if(logradsurf(i,k,m).le.xrmax)then
               pext(i)=0.5_r8*cext(i,1)
               do nc=2,ncoef
                  pext(i)=pext(i)+cheb(nc,m,i,k)*cext(i,nc)
               enddo
               pext(i)=exp(pext(i))
            else
               pext(i)=1.5_r8/(radsurf(i,k,m)*rhoh2o) ! geometric optics
            endif
!           convert from m2/kg water to m2/kg aerosol
            specpext(i)=pext(i)
            pext(i)=pext(i)*wetvol(i)*rhoh2o
!           write(iulog,*)'pext=',pext(i)
            pabs(i)=0.5_r8*cabs(i,1)
            pasm(i)=0.5_r8*casm(i,1)
            do nc=2,ncoef
               pabs(i)=pabs(i)+cheb(nc,m,i,k)*cabs(i,nc)
               pasm(i)=pasm(i)+cheb(nc,m,i,k)*casm(i,nc)
            enddo
            pabs(i)=pabs(i)*wetvol(i)*rhoh2o
            pabs(i)=max(0._r8,pabs(i))
	    pabs(i)=min(pext(i),pabs(i))

	    palb(i)=1._r8-pabs(i)/max(pext(i),1.e-40_r8)
	    palb(i)=1._r8-pabs(i)/max(pext(i),1.e-40_r8)

            dopaer(i)=pext(i)*mass(i,k)
	 end do

!           Save aerosol optical depth at longest visible wavelength
!           sum over layers
            if(savaervis)then
!               aerosol extinction (/m)
                do i=1,ncol
                   extinct(i,k)=extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                   absorb(i,k)=absorb(i,k)+pabs(i)*air_density(i,k)
                   aodvis(i)=aodvis(i)+dopaer(i)
                   aodabs(i)=aodabs(i)+pabs(i)*mass(i,k)
                   aodmode(i,m)=aodmode(i,m)+dopaer(i)
!++xl added 11/24/2009, dust AOD
                   if (wetvol(i)>1.e-40_r8) then
                     dustaodmode(i,m)=dustaodmode(i,m)+dopaer(i)*dustvol(i)/wetvol(i)
 !                  else
 !                    dustaodmode(i,m)=0.0_r8
                   endif
!--dust AOD
                   ssavis(i)=ssavis(i)+dopaer(i)*palb(i)
	end do
             endif
             do i=1,ncol
                if ((dopaer(i)>-1.e-10) .and. (dopaer(i)<30.))goto 125

                    write(iulog,*)'dopaer(',i,',',k,',',m,',',lchnk,')=',dopaer(i)
!                    write(iulog,*)'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                    write(iulog,*)'k=',k,' pext=',pext(i),' specext=',specpext(i)
                    write(iulog,*)'wetvol=',wetvol(i),' dryvol=',dryvol(i),' watervol=',watervol(i)
!                    write(iulog,*)'cext=',(cext(i,l),l=1,ncoef)
!                    write(iulog,*)'crefin=',crefin(i)
                    write(iulog,*)'nspec_amode(m)=',nspec_amode(m)
!                    write(iulog,*)'cheb=',(cheb(nc,m,i,k),nc=2,ncoef)
                     do l=1,nspec_amode(m)
                        lmass=lmassptr_amode(l,m)
                        ltype=lspectype_amode(l,m)
                        volf=raer(i,k,lmass)/specdens_amode(ltype)
                       write(iulog,*)'l=',l,'vol(l)=',volf
                       write(iulog,*)'specrefndxsw(l)=',specrefndxsw(iswband,ltype)
                       write(iulog,*)'specdens_amode(ltype)=',specdens_amode(ltype)
                     enddo

                     l=numptr_amode(m)
                     if ((l > 0) .and. (l .le. pcnst)) then
!                       write(iulog,*)'number=',raer(i,k,l)
                     else
!                       write(iulog,*)'number= ????, numptr=', l
                     endif

 !                    call endrun('exit from modal_aero_sw')
                     nerr_dopaer = nerr_dopaer + 1
                     if (nerr_dopaer >= nerrmax_dopaer) then
!                       write(iulog,*)'*** halting in modal_aero_sw after nerr_dopaer =', nerr_dopaer
!                       call endrun('exit from modal_aero_sw')
                     end if

125               continue
               end do
               do i=1,ncol
                  tauxar(i,k)=tauxar(i,k)+dopaer(i)
                  wa(i,k)=wa(i,k)+dopaer(i)*palb(i)
                  ga(i,k)=ga(i,k)+dopaer(i)*palb(i)*pasm(i)
                  fa(i,k)=fa(i,k)+dopaer(i)*palb(i)*pasm(i)*pasm(i)
	       end do
            end do ! ntot_amode
  end do ! pver
         if(savaervis)then
            do i=1,ncol
               if(aodvis(i)>1.e-10)then
                  ssavis(i)=ssavis(i)/aodvis(i)
               else
                  ssavis(i)=0.925_r8
               endif
			   do m=1,ntot_amode
			      colext(i,m)=0.001*aodmode(i,m)/burden(i,m)
!				  if(aodmode(i,m).gt.1.)write(iulog,*)'m,burden,tau,ext=',m, &
!                                                 burden(i,m),aodmode(i,m),colext(i,m)
			   end do
	    end do
            do i = 1, nnite
               extinct(idxnite(i),:) = fillvalue
               absorb(idxnite(i),:) = fillvalue
               aodvis(idxnite(i)) = fillvalue
               aodabs(idxnite(i)) = fillvalue
               aodmode(idxnite(i),:) = fillvalue
!++xl added 11/24/2009, dust AOD
               dustaodmode(idxnite(i),:) = fillvalue
!--xl
               burden(idxnite(i),:) = fillvalue
               ssavis(idxnite(i)) = fillvalue
            end do
            call outfld('EXTINCT',extinct ,pcols,lchnk)
            call outfld('ABSORB',absorb ,pcols,lchnk)
            call outfld('AODVIS',aodvis ,pcols,lchnk)
            call outfld('AODABS',aodabs ,pcols,lchnk)
            call outfld('AODMODE1',aodmode(1,1) ,pcols,lchnk)
            call outfld('AODMODE2',aodmode(1,2) ,pcols,lchnk)
            call outfld('AODMODE3',aodmode(1,3) ,pcols,lchnk)
#if (defined MODAL_AERO_7MODE)
            call outfld('AODMODE4',aodmode(1,4) ,pcols,lchnk)
            call outfld('AODMODE5',aodmode(1,5) ,pcols,lchnk)
            call outfld('AODMODE6',aodmode(1,6) ,pcols,lchnk)
            call outfld('AODMODE7',aodmode(1,7) ,pcols,lchnk)
#endif
!++xl added 11/24/2009, dust AOD
            call outfld('AODDUST1',dustaodmode(1,1) ,pcols,lchnk)
            call outfld('AODDUST2',dustaodmode(1,2) ,pcols,lchnk)
            call outfld('AODDUST3',dustaodmode(1,3) ,pcols,lchnk)
!--xl
            call outfld('BURDEN1',burden(1,1) ,pcols,lchnk)
            call outfld('BURDEN2',burden(1,2) ,pcols,lchnk)
            call outfld('BURDEN3',burden(1,3) ,pcols,lchnk)
#if (defined MODAL_AERO_7MODE)
            call outfld('BURDEN4',burden(1,4) ,pcols,lchnk)
            call outfld('BURDEN5',burden(1,5) ,pcols,lchnk)
            call outfld('BURDEN6',burden(1,6) ,pcols,lchnk)
            call outfld('BURDEN7',burden(1,7) ,pcols,lchnk)
#endif
            call outfld('SSAVIS',ssavis ,pcols,lchnk)
         endif

      return
      end subroutine modal_aero_sw

      subroutine modal_aero_lw(ncol,lchnk, pbuf, raer, tauxar, pdeldry)

!  calculates aerosol radiative properties

   use ppgrid
   use constituents,  only: pcnst
   use physconst, only: rga
   use phys_buffer,      only: pbuf_size_max, pbuf_fld,  pbuf_get_fld_idx

   implicit none

   integer, intent(in) :: ncol  ! number of columns
   integer, intent(in) :: lchnk  ! chunk index
   real(r8), intent(in) :: raer(pcols, pver, pcnst)
      real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth
      real(r8), intent(in) :: pdeldry(pcols,pver) ! pressure thickness (Pa)
    type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf


      integer i,k
      real(r8) pabs(pcols) ! parameterized specific absorption (m2/kg)
      real(r8) :: mass(pcols,pver) ! layer mass
      real(r8) totmass(pcols)
      real(r8) vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
      real(r8) dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
      real(r8) wetvol(pcols)   ! volume concentration of wet mode (m3/kg)
      real(r8) watervol(pcols) ! volume concentration of water in each mode (m3/kg)
      real(r8) volf            ! volume fraction of insoluble aerosol
      real(r8) refr(pcols)     ! real part of refractive index
      real(r8) refi(pcols)     ! imaginary part of refractive index
      complex crefin(pcols)    ! complex refractive index
      integer lmass,ltype
      integer itab(pcols),jtab(pcols)
      real(r8) ttab(pcols),utab(pcols)
      real(r8) cext(pcols,ncoef),cabs(pcols,ncoef),casm(pcols,ncoef)
      real(r8) xrad(pcols)
      real(r8) cheby(ncoef,pcols,pver,maxd_amode) ! chebychef polynomials
      integer, parameter ::  nerrmax_dopaer=1000
      integer nerr_dopaer
      save nerr_dopaer
      data nerr_dopaer/0/
      real(r8) dopaer(pcols) ! aerosol optical depth in layer
      real(r8) third
	 integer m,l,nc,ns,ifld
      real(r8), pointer, dimension(:,:,:) ::  dgnumwet  ! wet number mode diameter (m)
      real(r8), pointer, dimension(:,:,:) ::  qaerwat  ! aerosol water (g/g)
      integer ilwband

      ifld = pbuf_get_fld_idx('DGNUMWET')
      dgnumwet  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1:maxd_amode)
      ifld  = pbuf_get_fld_idx( 'QAERWAT' )
      if ( associated(pbuf(ifld)%fld_ptr) ) then
         qaerwat => pbuf(ifld)%fld_ptr( 1, 1:pcols, 1:pver, lchnk, 1:maxd_amode )
      else
         call endrun( 'pbuf for QAERWAT not allocated in modal_aer_day' )
      end if


      do k=1,pver
	 mass(:ncol,k)=pdeldry(:ncol,k)*rga !
      enddo
      third=1./3.

!     calc size parameter for all columns
      do k=1,pver
      do m=1,ntot_amode
         do i=1,ncol
!           convert from number diameter to surface area
            xrad(i)=log(0.5*dgnumwet(i,k,m))+2.0_r8*alnsg_amode(m)*alnsg_amode(m)
!           normalize size parameter
            xrad(i)=max(xrad(i),xrmin)
            xrad(i)=min(xrad(i),xrmax)
            xrad(i)=(2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
!           chebyshev polynomials
            cheby(1,i,k,m)=1.0_r8
            cheby(2,i,k,m)=xrad(i)
            do nc=3,ncoef
               cheby(nc,i,k,m)=2.0_r8*xrad(i)*cheby(nc-1,i,k,m)-cheby(nc-2,i,k,m)
            enddo
         enddo
      enddo
      enddo


   do 1000 ilwband=1,nlwbands


   do k=1,pver

      tauxar(:ncol,k,ilwband)=0._r8

!     loop over all aerosol modes

      do m=1,ntot_amode

!        form bulk refractive index. Use volume mixing for infrared

         do i=1,ncol
            crefin(i)=0._r8
            dryvol(i)=0._r8
	 end do

!        aerosol species loop

         do l=1,nspec_amode(m)
            ltype=lspectype_amode(l,m)
            lmass=lmassptr_amode(l,m)
	    do i=1,ncol
               vol(i)=raer(i,k,lmass)/specdens_amode(ltype)
               crefin(i)=crefin(i)+vol(i)*specrefndxlw(ilwband,ltype)
               dryvol(i)=dryvol(i)+vol(i)
	    end do
         enddo
         do i=1,ncol
       	    watervol(i)=qaerwat(i,k,m)/rhoh2o
	    wetvol(i)=watervol(i)+dryvol(i)
            if(watervol(i)<0.0_r8)then
               if(abs(watervol(i)).gt.1.e-1*wetvol(i))then
                  write(iulog,*)'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in modal_aero_lw'
              !   call endrun()
               endif
               watervol(i)=0._r8
               wetvol(i)=dryvol(i)
            endif
            crefin(i)=crefin(i)+watervol(i)*crefwlw(ilwband)
            if(wetvol(i)>1.e-40)crefin(i)=crefin(i)/wetvol(i)
            refr(i)=real(crefin(i))
            refi(i)=aimag(crefin(i))

         end do ! ncol

!        interpolate coefficients linear in refractive index
!        first call calcs itab,jtab,ttab,utab
         itab(:ncol)=0
         call binterp(absplw(1,1,1,m,ilwband),ncol,ncoef,nrefr,nrefi,   &
            refr,refi,refrtablw(1,ilwband),refitablw(1,ilwband),itab,jtab,  &
            ttab,utab,cabs)

!        parameterized optical properties
         do i=1,ncol
            pabs(i)=0.5*cabs(i,1)
            do nc=2,ncoef
               pabs(i)=pabs(i)+cheby(nc,i,k,m)*cabs(i,nc)
            enddo
            pabs(i)=pabs(i)*wetvol(i)*rhoh2o
            pabs(i)=max(0._r8,pabs(i))
            dopaer(i)=pabs(i)*mass(i,k)
	 end do

             do i=1,ncol
                  if ((dopaer(i)>-1.e-10) .and. (dopaer(i)<20.))goto 225

                     write(iulog,*)'dopaer(',i,',',k,',',m,',',lchnk,')=',dopaer(i)
                     write(iulog,*)'k=',k,' pabs=',pabs(i)
                     write(iulog,*)'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                             ' watervol=',watervol(i)
                     write(iulog,*)'cabs=',(cabs(i,l),l=1,ncoef)
                     write(iulog,*)'crefin=',crefin(i)
                     write(iulog,*)'nspec_amode(m)=',nspec_amode(m)
                     do l=1,nspec_amode(m)
                        lmass=lmassptr_amode(l,m)
                        ltype=lspectype_amode(l,m)
                        vol=raer(i,k,lmass)/specdens_amode(ltype)
                        write(iulog,*)'l=',l,'vol(l)=',vol
                        write(iulog,*)'specrefndxlw(l)=',specrefndxlw(ilwband,ltype)
                        write(iulog,*)'specdens_amode(ltype)=',specdens_amode(ltype)
                     enddo

                     l=numptr_amode(m)
                     if ((l > 0) .and. (l .le. pcnst)) then
                        write(iulog,*)'number=',raer(i,k,l)
                     else
                        write(iulog,*)'number= ????, numptr=', l
                     endif

!                    call endrun()
                     nerr_dopaer = nerr_dopaer + 1
                     if (nerr_dopaer >= nerrmax_dopaer) then
                        write(iulog,*)'*** halting in modal_aero_lw after nerr_dopaer =', nerr_dopaer
                        call endrun()
                     end if

225               continue
               end do
               do i=1,ncol
                  tauxar(i,k,ilwband)=tauxar(i,k,ilwband)+dopaer(i)
	       end do
            end do ! ntot_amode
  end do ! pver

  1000 continue

      return
      end subroutine modal_aero_lw

      subroutine fitcurv(rs,yin,coef,ncoef,maxm)

!     fit y(x) using Chebyshev polynomials

      implicit none
      integer nmodes,maxm,m,ncoef
      parameter (nmodes=300)

      real(r8) rs(maxm) ! surface mode radius (cm)
      real(r8) yin(maxm),coef(ncoef)
      real(r8) x(nmodes),y(nmodes)
      real(r8) xmin,xmax

      if(maxm.gt.nmodes)then
        write(iulog,*)'nmodes too small in fitcurv',maxm
        call endrun()
      endif

      xmin=1.e20
      xmax=-1.e20
      do 100 m=1,maxm
        x(m)=log(rs(m))
        xmin=min(xmin,x(m))
        xmax=max(xmax,x(m))
        y(m)=log(yin(m))
  100 continue

      do 110 m=1,maxm
      x(m)=(2*x(m)-xmax-xmin)/(xmax-xmin)
  110 continue

      call chebft(coef,ncoef,y,maxm)

      return
      end subroutine fitcurv

      subroutine fitcurvlin(rs,yin,coef,ncoef,maxm)

!     fit y(x) using Chebychev polynomials
!     calculates ncoef coefficients coef

      implicit none
      integer nmodes,m,ncoef
      parameter (nmodes=300)

      integer,intent(in) :: maxm
      real(r8),intent(in) :: rs(maxm),yin(maxm)
      real(r8) x(nmodes),y(nmodes)
      real(r8),intent(out) :: coef(ncoef)
      real(r8) xmin,xmax

      if(maxm.gt.nmodes)then
        write(iulog,*)'nmodes too small in fitcurv',maxm
        call endrun
      endif

      do 100 m=1,maxm
      x(m)=log(rs(m))
      y(m)=(yin(m))
  100 continue

      call chebft(coef,ncoef,y,maxm)

      return
      end subroutine fitcurvlin

      subroutine chebft(c,m,f,n)
!     given a function f with values at zeroes x_k of Chebychef polynomial
!     T_n(x), calculate coefficients c_j such that
!     f(x)=sum(k=1,n) c_k t_(k-1)(y) - 0.5*c_1
!     where y=(x-0.5*(xmax+xmin))/(0.5*(xmax-xmin))
!     See Numerical Recipes, pp. 148-150.

      implicit none
      real(r8) pi
      parameter (pi=3.14159265)
      real(r8) c(m),f(n),fac,sum
      integer n,j,k,m

      fac=2./n
      do j=1,m
         sum=0
         do k=1,n
            sum=sum+f(k)*cos((pi*(j-1))*((k-0.5)/n))
         enddo
         c(j)=fac*sum
      enddo
      return
      end subroutine chebft

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*)'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp

      subroutine  write_modal_optics(outfilename,  modes, lw_band_len, sw_band_len,  &
	 refindex_real_sw,refindex_im_sw, refindex_real_lw, refindex_im_lw, &
	 nrefr, nrefi, sigma_logr_aer, ncoef, extpsw, abspsw, asmpsw, absplw )
        use pio, only : pio_def_var, pio_def_dim, pio_put_var, pio_clobber, pio_put_att, pio_double, &
             file_desc_T, var_desc_t, pio_global, pio_enddef, pio_closefile
        use cam_pio_utils, only : cam_pio_createfile
        implicit none

      character*(*) outfilename
! error status return
      integer  ierr
! netCDF id
      type(file_desc_t) ::  ncid
! dimension ids
      integer  lw_band_dim
      integer  sw_band_dim
      integer  refindex_real_dim
      integer  refindex_im_dim
      integer  mode_dim
      integer  ncoef_dim
! dimension lengths
      integer, intent(in) :: nrefr,nrefi ! number of refractive indices
      integer, intent(in) :: ncoef ! number of coefficients in polynomial expressions
      integer, intent(in) :: sw_band_len ! number of solar wavelengths
      integer, intent(in) :: lw_band_len ! number of infrared wavelengths
      integer, intent(in) :: modes ! number of aerosol modes

! variable ids
      type(var_desc_t) ::   lw_abs_id
      type(var_desc_t) ::   sw_asm_id
      type(var_desc_t) ::   sw_ext_id
      type(var_desc_t) ::   sw_abs_id
      type(var_desc_t) ::   sigma_logr_aer_id
      type(var_desc_t) ::   refindex_real_lw_id
      type(var_desc_t) ::   refindex_im_lw_id
      type(var_desc_t) ::   refindex_real_sw_id
      type(var_desc_t) ::   refindex_im_sw_id
! rank (number of dimensions) for each variable
      integer  lw_abs_rank
      integer  sw_asm_rank
      integer  sw_ext_rank
      integer  sw_abs_rank
      integer  sigma_logr_aer_rank
      integer refindex_rank
      parameter (lw_abs_rank = 5)
      parameter (sw_asm_rank = 5)
      parameter (sw_ext_rank = 5)
      parameter (sw_abs_rank = 5)
      parameter (refindex_rank=2)
      parameter (sigma_logr_aer_rank = 1)
! variable shapes
      integer  lw_abs_dims(lw_abs_rank)
      integer  sw_asm_dims(sw_asm_rank)
      integer  sw_ext_dims(sw_ext_rank)
      integer  sw_abs_dims(sw_abs_rank)
      integer refindex_dims(refindex_rank)
      integer coef_dim
      integer sigma_logr_aer_dims(sigma_logr_aer_rank)
! data variables

!     coefficients for parameterizing aerosol radiative properties
!     in terms of refractive index and wet radius
      real(r8), intent(in) :: extpsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! specific extinction
      real(r8), intent(in) :: abspsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! specific absorption
      real(r8), intent(in) :: asmpsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! asymmetry factor
      real(r8), intent(in) :: absplw(ncoef,nrefr,nrefi,modes,lw_band_len) ! specific absorption

      real(r8), intent(in) ::  sigma_logr_aer(maxd_amode) ! geometric standard deviation of size distribution

!     tables of real refractive indices for aerosols
      real(r8), intent(in) :: refindex_real_sw(nrefr,sw_band_len) !
      real(r8), intent(in) :: refindex_im_sw(nrefi,sw_band_len) !
      real(r8), intent(in) :: refindex_real_lw(nrefr,lw_band_len) !
      real(r8), intent(in) :: refindex_im_lw(nrefi,lw_band_len) !

      integer lenchr


      integer len
! attribute vectors
! enter define mode
      call cam_pio_createfile(ncid, outfilename, PIO_CLOBBER)
! define dimensions
      ierr = pio_def_dim(ncid, 'lw_band', lw_band_len, lw_band_dim) 
      ierr = pio_def_dim(ncid, 'sw_band', sw_band_len, sw_band_dim)
      ierr = pio_def_dim(ncid, 'refindex_real', nrefr, refindex_real_dim)
      ierr = pio_def_dim(ncid, 'refindex_im', nrefi, refindex_im_dim)
      ierr = pio_def_dim(ncid, 'mode', modes, mode_dim)
      ierr = pio_def_dim(ncid, 'coef_number', ncoef, coef_dim)

! define variables

      lw_abs_dims(5) = lw_band_dim
      lw_abs_dims(4) = mode_dim
      lw_abs_dims(3) = refindex_im_dim
      lw_abs_dims(2) = refindex_real_dim
      lw_abs_dims(1) = coef_dim
      ierr = pio_def_var(ncid, 'absplw', PIO_DOUBLE, lw_abs_dims, lw_abs_id )
             
      sw_asm_dims(5) = sw_band_dim
      sw_asm_dims(4) = mode_dim
      sw_asm_dims(3) = refindex_im_dim
      sw_asm_dims(2) = refindex_real_dim
      sw_asm_dims(1) = coef_dim
      ierr = pio_def_var(ncid, 'asmpsw', PIO_DOUBLE, sw_asm_dims, sw_asm_id )

      sw_ext_dims(5) = sw_band_dim
      sw_ext_dims(4) = mode_dim
      sw_ext_dims(3) = refindex_im_dim
      sw_ext_dims(2) = refindex_real_dim
      sw_ext_dims(1) = coef_dim
      ierr = pio_def_var(ncid, 'extpsw', PIO_DOUBLE, sw_ext_dims, sw_ext_id )

      sw_abs_dims(5) = sw_band_dim
      sw_abs_dims(4) = mode_dim
      sw_abs_dims(3) = refindex_im_dim
      sw_abs_dims(2) = refindex_real_dim
      sw_abs_dims(1) = coef_dim
      ierr = pio_def_var(ncid, 'abspsw', PIO_DOUBLE, sw_abs_dims, sw_abs_id )

      refindex_dims(2) = sw_band_dim
      refindex_dims(1) = refindex_real_dim
      ierr = pio_def_var(ncid, 'refindex_real_sw', PIO_DOUBLE,  refindex_dims, &
          refindex_real_sw_id)
      refindex_dims(2) = sw_band_dim
      refindex_dims(1) = refindex_im_dim
      ierr = pio_def_var(ncid, 'refindex_im_sw', PIO_DOUBLE, refindex_dims, &
          refindex_im_sw_id)
      refindex_dims(2) = lw_band_dim
      refindex_dims(1) = refindex_real_dim
      ierr = pio_def_var(ncid, 'refindex_real_lw', PIO_DOUBLE, refindex_dims, &
          refindex_real_lw_id)
      refindex_dims(2) = lw_band_dim
      refindex_dims(1) = refindex_im_dim
      ierr = pio_def_var(ncid, 'refindex_im_lw', PIO_DOUBLE,  refindex_dims, &
          refindex_im_lw_id)
      sigma_logr_aer_dims(1) = mode_dim
      ierr = pio_def_var(ncid, 'sigma_logr_aer', PIO_DOUBLE, sigma_logr_aer_dims, &
          sigma_logr_aer_id)
! assign attributes
      ierr = pio_put_att(ncid, lw_abs_id, 'long_name',   &
         'coefficients of polynomial expression for longwave absorption')
      ierr = pio_put_att(ncid, lw_abs_id, 'units',  'meter^2 kilogram^-1')
      ierr = pio_put_att(ncid, sw_ext_id, 'long_name',  &
         'coefficients of polynomial expression for shortwave extinction') 
      ierr = pio_put_att(ncid, sw_ext_id, 'units',  'meter^2 kilogram^-1')
      ierr = pio_put_att(ncid, sw_asm_id, 'long_name',  &
         'coefficients of polynomial expression for shortwave asymmetry parameter')
      ierr = pio_put_att(ncid, sw_asm_id, 'units',  'fraction')
      ierr = pio_put_att(ncid, sw_abs_id, 'long_name',  &
         'coefficients of polynomial expression for shortwave absorption')
      ierr = pio_put_att(ncid, sw_abs_id, 'units',  'meter^2 kilogram^-1')
      ierr = pio_put_att(ncid, refindex_real_lw_id, 'long_name',  &
            'real refractive index of aerosol - longwave')
      ierr = pio_put_att(ncid, refindex_im_lw_id,   'long_name',  &
            'imaginary refractive index of aerosol - longwave')
      ierr = pio_put_att(ncid, refindex_real_sw_id, 'long_name',  &
            'real refractive index of aerosol - shortwave')
      ierr = pio_put_att(ncid, refindex_im_sw_id,   'long_name',  &
            'imaginary refractive index of aerosol - shortwave')
      ierr = pio_put_att(ncid, sigma_logr_aer_id, 'long_name',  &
            'geometric standard deviation of aerosol')
      ierr = pio_put_att(ncid, sigma_logr_aer_id, 'units',  '-')

      ierr = pio_put_att(ncid, PIO_GLOBAL, 'source',  'OPAC + miev0')
      ierr = pio_put_att(ncid, PIO_GLOBAL, 'history',  'Ghan and Zaveri, JGR 2007')

! leave define mode

      ierr = pio_enddef(ncid)

      ierr = pio_put_var(ncid,  sigma_logr_aer_id,  sigma_logr_aer)

      ierr = pio_put_var(ncid, refindex_real_lw_id, refindex_real_lw)

      ierr = pio_put_var(ncid, refindex_im_lw_id, refindex_im_lw)

      ierr = pio_put_var(ncid, refindex_real_sw_id, refindex_real_sw)

      ierr = pio_put_var(ncid, refindex_im_sw_id, refindex_im_sw)

      ierr = pio_put_var(ncid, lw_abs_id, absplw)

       ierr = pio_put_var(ncid, sw_ext_id, extpsw)

       ierr = pio_put_var(ncid, sw_abs_id, abspsw)

      ierr = pio_put_var(ncid, sw_asm_id, asmpsw)

      call pio_closefile(ncid)

      end subroutine write_modal_optics

      subroutine  read_modal_optics(infilename,  modes, lw_band_len, sw_band_len,  &
	 refindex_real_sw,refindex_im_sw, refindex_real_lw, refindex_im_lw, &
	 nrefr, nrefi, ncoef, extpsw, abspsw, asmpsw, absplw )

      use pio, only : pio_inq_dimlen, pio_inq_dimid, pio_get_var, pio_inq_varid, file_desc_t, var_desc_t, &
           pio_nowrite, pio_closefile
      use cam_pio_utils, only : cam_pio_openfile

      implicit none

      character*(*), intent(in) :: infilename

! netCDF id
      type(file_desc_t) ::  ncid
! dimension ids
      integer  lw_band_dim
      integer  sw_band_dim
      integer  refindex_real_dim
      integer  refindex_im_dim
      integer  mode_dim
      integer  ncoef_dim
! dimension lengths
      integer, intent(in) :: nrefr,nrefi ! number of refractive indices
      integer, intent(in) :: ncoef ! number of coefficients in polynomial expressions
      integer, intent(in) :: sw_band_len ! number of solar wavelengths
      integer, intent(in) :: lw_band_len ! number of infrared wavelengths
      integer, intent(in) :: modes ! number of aerosol modes

! variable ids
      type(var_desc_t) ::  lw_abs_id
      type(var_desc_t) ::  sw_asm_id
      type(var_desc_t) ::  sw_ext_id
      type(var_desc_t) ::  sw_abs_id
      type(var_desc_t) ::  sigma_logr_aer_id
      type(var_desc_t) ::  refindex_real_lw_id
      type(var_desc_t) ::  refindex_im_lw_id
      type(var_desc_t) ::  refindex_real_sw_id
      type(var_desc_t) ::  refindex_im_sw_id
! rank (number of dimensions) for each variable
      integer  lw_abs_rank
      integer  sw_asm_rank
      integer  sw_ext_rank
      integer  sw_abs_rank
      integer  sigma_logr_aer_rank
      integer refindex_rank
      parameter (lw_abs_rank = 5)
      parameter (sw_asm_rank = 5)
      parameter (sw_ext_rank = 5)
      parameter (sw_abs_rank = 5)
      parameter (refindex_rank=2)
      parameter (sigma_logr_aer_rank = 1)
! variable shapes
      integer  lw_abs_dims(lw_abs_rank)
      integer  sw_asm_dims(sw_asm_rank)
      integer  sw_ext_dims(sw_ext_rank)
      integer  sw_abs_dims(sw_abs_rank)
      integer refindex_dims(refindex_rank)
      integer coef_dim
      integer sigma_logr_aer_dims(sigma_logr_aer_rank)
! data variables

!     coefficients for parameterizing aerosol radiative properties
!     in terms of refractive index and wet radius
      real(r8), intent(out) :: extpsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! specific extinction
      real(r8), intent(out) :: abspsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! specific absorption
      real(r8), intent(out) :: asmpsw(ncoef,nrefr,nrefi,modes,sw_band_len) ! asymmetry factor
      real(r8), intent(out) :: absplw(ncoef,nrefr,nrefi,modes,lw_band_len) ! specific absorption

      real(r8) ::  sigma_logr_aer(maxd_amode) ! geometric standard deviation of size distribution

!     tables of real refractive indices for aerosols
      real(r8), intent(out) :: refindex_real_sw(nrefr,sw_band_len) !
      real(r8), intent(out) :: refindex_im_sw(nrefi,sw_band_len) !
      real(r8), intent(out) :: refindex_real_lw(nrefr,lw_band_len) !
      real(r8), intent(out) :: refindex_im_lw(nrefi,lw_band_len) !

      integer lenchr
      integer dimread
      integer start(5),count(5)
      integer m, ierr

      integer len
! open file
       call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)
! inquire dimensions
      ierr = pio_inq_dimid(ncid,'lw_band',lw_band_dim)
      ierr = pio_inq_dimlen(ncid, lw_band_dim, dimread)
      if(dimread.ne.lw_band_len)then
         write(iulog,*)'lw_band_len=',dimread,' from ',infilename,' ne lw_band_len=',lw_band_len
	 call endrun
      endif
      ierr = pio_inq_dimid(ncid,'sw_band',sw_band_dim)
      ierr = pio_inq_dimlen(ncid, sw_band_dim, dimread)
      if(dimread.ne.sw_band_len)then
         write(iulog,*)'sw_band_len=',dimread,' from ',infilename,' ne sw_band_len=',sw_band_len
	 call endrun
      endif
      ierr = pio_inq_dimid(ncid,'refindex_real',refindex_real_dim)
      ierr = pio_inq_dimlen(ncid, refindex_real_dim, dimread)
      if(dimread.ne.nrefr)then
         write(iulog,*)'refindex_real=',dimread,' from ',infilename,' ne nrefr=',nrefr
	 call endrun
      endif
      ierr = pio_inq_dimid(ncid,'refindex_im',refindex_im_dim)
      ierr = pio_inq_dimlen(ncid, refindex_im_dim, dimread)
      if(dimread.ne.nrefi)then
         write(iulog,*)'refindex_im=',dimread,' from ',infilename,' ne nrefi=',nrefi
	 call endrun
      endif
      ierr = pio_inq_dimid(ncid,'mode',mode_dim)
      ierr = pio_inq_dimlen(ncid, mode_dim, dimread)
      if(dimread.ne.modes)then
         write(iulog,*)'mode=',dimread,' from ',infilename,' ne modes=',modes
	 call endrun
      endif
      ierr = pio_inq_dimid(ncid,'coef_number',coef_dim)
      ierr = pio_inq_dimlen(ncid, coef_dim, dimread)
      if(dimread.ne.ncoef)then
         write(iulog,*)'coef_number=',dimread,' from ',infilename,' ne ncoef=',ncoef
	 call endrun
      endif

! read variables
      ierr = pio_inq_varid(ncid,'absplw',lw_abs_id)
      ierr = pio_get_var(ncid, lw_abs_id, absplw)

      ierr = pio_inq_varid(ncid,'extpsw',sw_ext_id)
      ierr = pio_get_var(ncid, sw_ext_id, extpsw)

      ierr = pio_inq_varid(ncid,'abspsw',sw_abs_id)
      ierr = pio_get_var(ncid, sw_abs_id, abspsw)

      ierr = pio_inq_varid(ncid,'asmpsw',sw_asm_id)
      ierr = pio_get_var(ncid, sw_asm_id, asmpsw)

      ierr = pio_inq_varid(ncid, 'refindex_real_sw', refindex_real_sw_id)
      ierr = pio_get_var(ncid, refindex_real_sw_id, refindex_real_sw)

      ierr = pio_inq_varid(ncid, 'refindex_im_sw', refindex_im_sw_id)
      ierr = pio_get_var(ncid, refindex_im_sw_id, refindex_im_sw)

      ierr = pio_inq_varid(ncid, 'refindex_real_lw', refindex_real_lw_id)
      ierr = pio_get_var(ncid, refindex_real_lw_id, refindex_real_lw)

      ierr = pio_inq_varid(ncid, 'refindex_im_lw', refindex_im_lw_id)
      ierr = pio_get_var(ncid, refindex_im_lw_id, refindex_im_lw)

      ierr = pio_inq_varid(ncid, 'sigma_logr_aer', sigma_logr_aer_id)
      ierr = pio_get_var(ncid, sigma_logr_aer_id, sigma_logr_aer)
      do m=1,modes
         if(sigma_logr_aer(m)/sigmag_amode(m).lt.0.99.or. &
            sigma_logr_aer(m)/sigmag_amode(m).gt.1.01)then
	    write(iulog,*)'sigma_logr_aer,sigmag_amode=',sigma_logr_aer(m),sigmag_amode(m)
	    call endrun
	 end if
      end do

      call pio_closefile(ncid)

      end subroutine read_modal_optics
#endif
!MODAL_AERO
end module modal_aer_opt
