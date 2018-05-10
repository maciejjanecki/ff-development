  program read_grid

  implicit none
  integer, parameter :: char_len  = 256, &
                            int_kind  = kind(1),  &
                            log_kind  = kind(.true.), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13)

  integer (kind=int_kind), parameter :: imt = 1000, &
       jmt = 640

  integer(selected_int_kind(12)) :: polei, &
      polej, idist, jdist,i,j,k,in3,nrecl,ii,jj,kk
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
       ridist_rad, rjdist_rad,wrk,wrk1
  logical,parameter :: writeParameters = .true.

        character(char_len) :: in_file,outFile

!        in_file = 'waterPuck.pop.r.2010-01-01-00000'
      in_file = '115_zeros.pop.r.2016-07-01-00000'       
 2000 format(32x,i2)

      in3=22
 2002 format('grids/grid_bs01v1_rec_',i2.2,'.bin')
      nrecl=imt*jmt*2*4
      open(in3,file=trim(in_file), &
       form='unformatted',access='direct',recl=nrecl, convert='big_endian')
      write(*,*) 'infile:',trim(in_file)
        pi=atan(c1)*c4
        radius=6370.0
     do i = 1,344
       write(outFile,2002) i
!       open(55,file=trim(outFile),&
!       form='unformatted',access='direct',recl=nrecl)!, convert='big_endian')
       read(in3,rec=i) rstData
       write(*,'(i4,A,2f14.6)') i,'   max/min = ',maxval(rstData),minval(rstData)
!       if (i <= 2 .or. i == 7) then
!          full_grid(:,:,i)=full_grid(:,:,i)*180./pi
!       endif
!       write(55,rec=1) full_grid(:,:,i)
!       close(55)
     enddo

!     write(*,*) 'checking for the NaN'
!     do kk = 1,7
!       do ii = 1,imt
!       do jj = 1,jmt
!        if (full_grid(ii,jj,kk) .ne. full_grid(ii,jj,kk)) then
!          write(*,*) kk,ii,jj,full_grid(ii,jj,kk)
!        endif
!       enddo
!       enddo
!     enddo


      close (in3)
! create init_ts file with 33 additional records
      open(in3,file=trim(in_file), &
          form='unformatted',access='direct',recl=nrecl, convert='big_endian')
      open(33,file='init_ts_bs01v1_20160701.ieeer8', &
          form='unformatted',access='direct',recl=nrecl, convert='big_endian',status = 'replace')
      !write temperature
      k=1
      wrk1 = 7.     
      do i=147,179
           read(in3,rec=i) rstData
           wrk = sum(rstData)/imt/jmt
           where (rstData<=0.) rstData=wrk
!           rstData = wrk1
           write(33,rec=k) rstData
           write(*,'(2i6,A,3f14.6)') i,k,'   max/min = ',maxval(rstData),minval(rstData),sum(rstData)/imt/jmt
           k = k + 1
           wrk1 = wrk1 - 0.1
      enddo
      wrk1 = 7.
      do i=213,245
           read(in3,rec=i) rstData
           rstData = rstData * 1000._dbl_kind
!           rstData = wrk
           wrk = sum(rstData)/imt/jmt
           where (rstData <= 0.) rstData=wrk 
!           rstData = wrk1
           write(33,rec=k) rstData
           write(*,'(2i6,A,3f14.6)') i,k,'   max/min = ',maxval(rstData),minval(rstData),sum(rstData)/imt/jmt
           k = k + 1
           wrk1 = wrk1 + 0.1
      enddo
      do i=213,245
           rstData = 0._dbl_kind
           write(33,rec=k) rstData
           k = k + 1
      enddo
      close(in3)
      close(33)
      
      

end

