! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Feb 2008 - A. Bodas-Salcedo - Initial version
!

#include "cosp_defs.h"
PROGRAM COSPTEST
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  USE MOD_COSP
  USE MOD_COSP_IO
  IMPLICIT NONE

  ! Local variables
  character(len=64)  :: cosp_input_nl='cosp_input_nl.txt'
  character(len=64)  :: cosp_output_nl='cosp_output_nl.txt'
!   character(len=64)  :: cosp_output_nl='./cfmip2/cosp_output_cfmip2_short_offline.txt'
  character(len=512) :: finput ! Input file name
  character(len=512) :: cmor_nl
  character(len=8)  :: wmode ! Writing mode 'replace' or 'append'
  integer :: overlap   !  overlap type: 1=max, 2=rand, 3=max/rand
  integer :: isccp_topheight,isccp_topheight_direction
  integer :: Ncolumns     ! Number of subcolumns in SCOPS
  integer :: Npoints      ! Number of gridpoints
  integer :: Nlevels      ! Number of levels
  integer :: Nlr          ! Number of levels in statistical outputs
  integer :: Npoints_it   ! Max number of gridpoints to be processed in one iteration
  integer,parameter :: ntsteps=5 
  integer :: i
  type(cosp_config) :: cfg   ! Configuration options
  type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
  type(cosp_subgrid) :: sgx     ! Subgrid outputs
  type(cosp_sgradar) :: sgradar ! Output from radar simulator
  type(cosp_sglidar) :: sglidar ! Output from lidar simulator
  type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
  type(cosp_modis)   :: modis   ! Output from MODIS simulator
  type(cosp_misr)    :: misr    ! Output from MISR simulator
#ifdef RTTOV 
  type(cosp_rttov)   :: rttov   ! Output from RTTOV 
#endif 
  type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  real,dimension(:),allocatable :: lon,lat
  real,dimension(:,:),allocatable,target :: p,ph,zlev,zlev_half,T,sh,rh,tca,cca, &
                    mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
                    fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
  real,dimension(:,:,:),allocatable :: Reff
  real,dimension(:),allocatable :: skt,landmask,sfc_height,u_wind,v_wind,sunlit
  integer :: t0,t1,count_rate,count_max
  integer :: Nlon,Nlat,geomode
  real :: radar_freq,k2,ZenAng,co2,ch4,n2o,co,emsfc_lw
  integer,dimension(RTTOV_MAX_CHANNELS) :: Channels
  real,dimension(RTTOV_MAX_CHANNELS) :: Surfem
  integer :: surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay
  integer :: Nprmts_max_hydro,Naero,Nprmts_max_aero,lidar_ice_type
  integer :: platform,satellite,Instrument,Nchannels
  logical :: use_vgrid,csat_vgrid,use_precipitation_fluxes,use_reff
  namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight,isccp_topheight_direction, &
              npoints,npoints_it,ncolumns,nlevels,use_vgrid,nlr,csat_vgrid,finput, &
              radar_freq,surface_radar,use_mie_tables, &
              use_gas_abs,do_ray,melt_lay,k2,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
              lidar_ice_type,use_precipitation_fluxes,use_reff, &
              platform,satellite,Instrument,Nchannels, &
              Channels,Surfem,ZenAng,co2,ch4,n2o,co

  double precision :: time(ntsteps)=(/1.D0,2.D0,3.D0,4.D0,5.D0/)
  double precision :: time_bnds(2,ntsteps)

  !---------------- End of declaration of variables --------------
   

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read COSP namelists
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(10,file=cosp_input_nl,status='old')
  read(10,nml=cosp_input)
  close(10)
  call read_cosp_output_nl(cosp_output_nl,cfg)
  print *, 'Processing file: ', finput
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Allocate local arrays
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(lon(Npoints),lat(Npoints), &
           p(Npoints,Nlevels),ph(Npoints,Nlevels), &
           zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels),T(Npoints,Nlevels), &
           sh(Npoints,Nlevels),rh(Npoints,Nlevels), &
           tca(Npoints,Nlevels),cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels), &
           mr_lsice(Npoints,Nlevels),mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels), &
           fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels), &
           fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels), &
           dem_s(Npoints,Nlevels),dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints),sfc_height(Npoints), &
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints))
  
  ! Example that processes ntsteps. It always uses the same input data
  wmode = 'replace' ! Only for first iteration
  do i=1,1
        time_bnds(:,i) = (/time(i),time(i)/) ! This is just for example purposes
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Read input geophysical variables from NetCDF file
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! input : surface to top
        call nc_read_input_file(finput,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, &
                mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff, &
                dtau_s,dtau_c,dem_s,dem_c,skt,landmask,sfc_height,mr_ozone,u_wind,v_wind,sunlit, &
                emsfc_lw,geomode,Nlon,Nlat)
                ! geomode = 2 for (lon,lat) mode.
                ! geomode = 3 for (lat,lon) mode.
                ! In those modes it returns Nlon and Nlat with the correct values
        call system_clock(t0,count_rate,count_max) !!! Only for testing purposes
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Allocate memory for gridbox type
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Allocating memory for gridbox type...'
        call construct_cosp_gridbox(time(i),time_bnds(:,i),radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
                                    do_ray,melt_lay,k2, &
                                    Npoints,Nlevels,Ncolumns,N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                                    lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
                                    use_precipitation_fluxes,use_reff, &
                                    Platform,Satellite,Instrument,Nchannels,ZenAng, &
                                    channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,gbx)
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Here code to populate input structure
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Populating input structure...'
        gbx%longitude = lon
        gbx%latitude = lat
        gbx%p = p
        gbx%ph = ph
        gbx%zlev = zlev
        gbx%zlev_half = zlev_half
        gbx%T = T
        gbx%q = rh
        gbx%sh = sh
        gbx%cca = cca
        gbx%tca = tca
        gbx%psfc = ph(:,1)
        gbx%skt  = skt
        gbx%land = landmask
!         gbx%sfc_height  = sfc_height
        gbx%mr_ozone  = mr_ozone
        gbx%u_wind  = u_wind
        gbx%v_wind  = v_wind
        gbx%sunlit  = sunlit
        
        gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
        gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice
        gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
        gbx%rain_ls = fl_lsrain
        gbx%snow_ls = fl_lssnow
        gbx%grpl_ls = fl_lsgrpl
        gbx%rain_cv = fl_ccrain
        gbx%snow_cv = fl_ccsnow
        
        gbx%Reff = Reff

        ! ISCCP simulator
        gbx%dtau_s   = dtau_s
        gbx%dtau_c   = dtau_c
        gbx%dem_s    = dem_s
        gbx%dem_c    = dem_c

               
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Define new vertical grid
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Defining new vertical grid...'
        call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Allocate memory for other types
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Allocating memory for other types...'
        call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
        call construct_cosp_sgradar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
        call construct_cosp_radarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)
        call construct_cosp_sglidar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
        call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
        call construct_cosp_isccp(cfg,Npoints,Ncolumns,Nlevels,isccp)
        call construct_cosp_modis(cfg,Npoints,modis)
        call construct_cosp_misr(cfg,Npoints,misr)
#ifdef RTTOV 
        call construct_cosp_rttov(Npoints,Nchannels,rttov) 
#endif
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Call simulator
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Calling simulator...'
#ifdef RTTOV 
        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
#else
        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
#endif
        

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Write outputs to CMOR-compliant NetCDF
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (i /= 1) wmode = 'append'
        gbx%time = time(i)
        if (cfg%Lwrite_output) then
            print *, 'Writing outputs...'
            if (geomode == 1) then 
#ifdef RTTOV 
               call nc_write_cosp_1d(cmor_nl,wmode,cfg,vgrid,gbx,sgx,sgradar,sglidar, &
                                     isccp,misr,modis,rttov,stradar,stlidar)
#else
               call nc_write_cosp_1d(cmor_nl,wmode,cfg,vgrid,gbx,sgx,sgradar,sglidar, &
                                     isccp,misr,modis,stradar,stlidar)
#endif
            elseif (geomode >  1) then
#ifdef RTTOV 
               call nc_write_cosp_2d(cmor_nl,wmode,cfg,vgrid,gbx,sgx,sgradar,sglidar, &
                                     isccp,misr,modis,rttov,stradar,stlidar,geomode,Nlon,Nlat)
#else
               call nc_write_cosp_2d(cmor_nl,wmode,cfg,vgrid,gbx,sgx,sgradar,sglidar, &
                                     isccp,misr,modis,stradar,stlidar,geomode,Nlon,Nlat)
#endif
            endif
        endif

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Deallocate memory in derived types
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Deallocating memory...'
        call free_cosp_gridbox(gbx)
        call free_cosp_subgrid(sgx)
        call free_cosp_sgradar(sgradar)
        call free_cosp_radarstats(stradar)
        call free_cosp_sglidar(sglidar)
        call free_cosp_lidarstats(stlidar)
        call free_cosp_isccp(isccp)
        call free_cosp_misr(misr)
        call free_cosp_modis(modis)
#ifdef RTTOV 
        call free_cosp_rttov(rttov) 
#endif
        call free_cosp_vgrid(vgrid)  
  enddo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Deallocate memory in local arrays
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  deallocate(lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, mr_lsliq,mr_lsice,mr_ccliq,mr_ccice, &
           fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,dem_c,skt, &
           landmask,sfc_height,mr_ozone,u_wind,v_wind,sunlit)

  ! Time in s. Only for testing purposes
  call system_clock(t1,count_rate,count_max)
  print *,(t1-t0)*1.0/count_rate
    
END PROGRAM COSPTEST
