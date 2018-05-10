  program read_grid

  implicit none
  integer, parameter :: char_len  = 256, &
                            int_kind  = kind(1),  &
                            log_kind  = kind(.true.), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13)

  integer (kind=int_kind), parameter :: imt = 1000, &
!       jmt = 640, kmt = 66
        jmt = 640, kmt = 33

  integer(selected_int_kind(12)) :: polei, &
      polej, idist, jdist,i,j,k,in3,nrecl,ii,jj,kk,nrec
  real (dbl_kind),parameter :: c1=1, c4=4

  real (dbl_kind) ::  tlat2(0:imt+2,0:jmt+1),tlat(0:imt+2,0:jmt+1), &
       tlat_m(0:imt+2,0:jmt+1), tlon_m(0:imt+2,0:jmt+1), tlon(0:imt+2,0:jmt+1), &
       ulat(0:imt+1,0:jmt), &
       ulon(0:imt+1,0:jmt)
  real (dbl_kind) :: ulato(imt,jmt), tlat_model(imt,jmt), tlon_model(imt,jmt), &
      ulono(imt,jmt),htn(imt,jmt),hte(imt,jmt), &
      huw(imt,jmt),hus(imt,jmt),angle(imt,jmt),full_grid(imt,jmt,7), &
      !data for restart
      rstData(imt,jmt)
  real (dbl_kind) :: radius, radian, radr, deg_km, &
       radrad, dxdeg, dydeg, dx, dy, dvx, ridist, rjdist, pi, &
       tlonij, tlonpj, tlonip, tlonpp,  &
       ridist_rad, rjdist_rad
  real(dbl_kind) :: rtmp
  logical,parameter :: writeParameters = .true.
        

        character(char_len) :: in_file,outFile

        in_file = 'init_ts_bs01v1_20171214.ieeer8'
        in_file = 'simple_ts.ieeer8'
!      in_file = '../../../../../../cesm_input_data/ocn/pop/bs2v3/ic/init_ts_bs2v3_jan_20130618.ieeer8'
        in_file = 'temp_in'
 2000 format(32x,i2)

      in3=22
 2002 format('temp_WP_rec_',i2.2,'.bin')
      nrecl=imt*jmt*2*4
!      write(*,*) trim(in_file)
!      stop
      open(in3,file=trim(in_file), &
       form='unformatted',access='direct',recl=nrecl, convert='big_endian')
      write(*,*) 'infile:',trim(in_file)
        pi=atan(c1)*c4
        radius=6370.0
     nrec = 1
!     do j = 1,24
     do i = 1,kmt
       read(in3,rec=nrec) rstData
       if (writeParameters) then  
!           if (j==1) then  
                write(outFile,2002) i
                open(55,file=trim(outFile),&
                          form='unformatted',access='direct',recl=nrecl,status='replace')!, convert='big_endian')
                write(55,rec=1) rstData
                close(55)
!           endif
       endif
       write(*,'(3i6,A,2f14.6)') i,j,nrec,'   max/min = ',maxval(rstData),minval(rstData)
       nrec = nrec + 1
      
     enddo
!     enddo
     close(in3)
!     open(in3,file=trim(in_file), &
!       form='unformatted',access='direct',recl=8*4, convert='big_endian')
!      write(*,*) 'infile:',trim(in_file)
!      do i=1,640*1000
!         read(in3,rec=i) rtmp
!         write(*,'(i10,f16.6)') i,rtmp
!      enddo
    close(in3) 
    
end program
