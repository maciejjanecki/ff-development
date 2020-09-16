 program TS_boundaries_data

 integer,parameter :: ni=1000,nj=640,nk=26,mnths=12, &
                      r8 = selected_real_kind(13)
!ustawienia*
 integer, parameter :: nhours = 4,ntracers=2,yearOffset = 0
 
!***********
 real*8 :: t(nk),s(nk),odata(ni,nj,nk),factors(ntracers)
 integer :: ii,jj,kk,ll,lvls,daysInMonths(mnths),year,month,day,hour, &
            sDate(4),eDate(4),hours(nhours),sHH,eHH,seconds(nhours)
 character*256 :: pathIn,pathOut
 character*4 :: tracer(ntracers)
 integer :: yy,mm,dd,hh,dayOfYear
 integer :: iidx,jidx
 character(len=8),dimension(4) :: outFNames
 character*256 :: inFile,outFile
 character*4   :: sDateC,eDateC
 logical ::         lbath = .true., &
            writeTestData = .false.

!ustawienia*
 data sDate /2012, 1, 1, 1/ !data poczatkowa 
 data eDate /2019,12,31, 4/ !data koncowa
 data seconds /3600,25200,46800,68400/
! data hours /0,6,12,18/
 data hours /3,9,15,21/
 data outFnames /"SST","SSS","TEMP3D","SALT3D"/

 write(sDateC,'(i4.4)') sDate(1)
 write(eDateC,'(i4.4)') eDate(1)

!!! >>>>>>  pathIn  = '../../tmp_data/ARTUR/lbc_FF/'//sDateC
! pathIn  
! 2000 format('../../../tmp_data/ARTUR/lbc_FF/',i4.4)
 2000 format('/users/work/ffish/boundary_575m/',i4.4) 
 pathOut = '/users/work/ffish/cesm_input_data/ocn/pop/bs05v1/forcing' 
 1000 format(i4.4,'-',i2.2,'-',i2.2,'-',i5.5,'_',a4,'_1000_0640_0026_0001.ieeer8') 
 1010 format(A,'.',i4.4,'.',i3.3,'.',i2.2)
!***********

 data tracer /'TEMP','SALT'/
 data factors /1.,1./
 data daysInMonths /31,28,31,30,31,30,31,31,30,31,30,31/
 

 do yy=sDate(1),eDate(1)
   dayOfYear = 1
   do mm=sDate(2),eDate(2)
      do dd=1,daysInMonths(mm) !sDate(3),eDate(3)
         do hh=sDate(4),eDate(4)
           do ii = 1,ntracers
              !read in data
 !!! =====>>>>>>             !!!!write(inFile,1000) yy,mm,dd,seconds(hh),tracer(ii)
              write(pathIn,2000) yy
              write(inFile,1000) yy,mm,dd,seconds(hh),tracer(ii)
              inFile = trim(pathIn)//'/'//trim(inFile)
              write(*,'(A)') trim(inFile)
              open(10,file=trim(inFile),status='old',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              do kk=1,nk
!                 write(*,*) 'record number ',kk
                 read(10,rec=kk) odata(:,:,kk)
                 odata(:,:,kk) = odata(:,:,kk)*factors(ii)
!                 write(*,'(i4,2f14.6,"   ",A)') kk,maxval(odata(:,:,kk)),minval(odata(:,:,kk)),tracer(ii) 
                 write(outFile,"(a4,i2.2,'.bin')") tracer(ii),kk
!                 write(*,'(A)') outFile
                 
                 if (lbath) then
                    !! correction for changet bathymetry
                           ! set edges
                           !    bath(1:3,:)=ic0
                           !    bath(imt-2:imt,:)=ic0
                           !    bath(:,1:3)=ic0
                           !    bath(:,jmt-2:jmt)=ic0
                           !! set boundaries
                           !    bath(1:50,316:510) = 5
                           !    bath(51:80,330:510) = 5
                           !    bath(327:900,590:640) = 5
                    do lvls=3,5
                       if (lvls==4) then
                             odata(493,628,lvls) = ((odata(492,628,lvls)+odata(494,628,lvls))/2 + &
                                                    (odata(493,627,lvls)+odata(493,629,lvls))/2)/2
                             odata(494,627,lvls) = ((odata(493,627,lvls)+odata(495,627,lvls))/2 + &
                                                    (odata(494,626,lvls)+odata(494,628,lvls))/2)/2
                       endif
                       if (lvls==5) then
                             do jidx=630,631
                                odata(492:495,jidx,lvls) = (odata(491,jidx,lvls)+odata(496,jidx,lvls))/2
                             enddo
                             jidx = 629
                             odata(491:496,jidx,lvls) = (odata(490,jidx,lvls)+odata(497,jidx,lvls))/2
                             odata(500:503,jidx,lvls) = (odata(499,jidx,lvls)+odata(504,jidx,lvls))/2
                             do jidx=627,628
                                odata(491:503,jidx,lvls) = (odata(490,jidx,lvls)+odata(504,jidx,lvls))/2
                             enddo
                             jidx = 626
                             odata(492:500,jidx,lvls) = (odata(491,jidx,lvls)+odata(501,jidx,lvls))/2
                             do jidx=624,625
                                odata(488:500,jidx,lvls) = (odata(487,jidx,lvls)+odata(501,jidx,lvls))/2
                             enddo
                             jidx = 623
                             odata(489:497,jidx,lvls) = (odata(488,jidx,lvls)+odata(498,jidx,lvls))/2
                             jidx = 622
                             odata(486:495,jidx,lvls) = (odata(485,jidx,lvls)+odata(496,jidx,lvls))/2
                             jidx = 621
                             odata(487:494,jidx,lvls) = (odata(486,jidx,lvls)+odata(495,jidx,lvls))/2
                             jidx = 620
                             odata(487:493,jidx,lvls) = (odata(486,jidx,lvls)+odata(494,jidx,lvls))/2
                             jidx = 619
                             odata(486:493,jidx,lvls) = (odata(485,jidx,lvls)+odata(494,jidx,lvls))/2
                             jidx = 618
                             odata(485:493,jidx,lvls) = (odata(484,jidx,lvls)+odata(494,jidx,lvls))/2
                             jidx = 618
                             odata(485:493,jidx,lvls) = (odata(484,jidx,lvls)+odata(494,jidx,lvls))/2
                             jidx = 617
                             odata(485:492,jidx,lvls) = (odata(484,jidx,lvls)+odata(493,jidx,lvls))/2
                             jidx = 616
                             odata(485:490,jidx,lvls) = (odata(484,jidx,lvls)+odata(491,jidx,lvls))/2
                             jidx = 615
                             odata(486:491,jidx,lvls) = (odata(485,jidx,lvls)+odata(492,jidx,lvls))/2
                             jidx = 614
                             odata(486:490,jidx,lvls) = (odata(485,jidx,lvls)+odata(491,jidx,lvls))/2
                             jidx = 613
                             odata(487:490,jidx,lvls) = (odata(486,jidx,lvls)+odata(491,jidx,lvls))/2
                             jidx = 612
                             odata(487:488,jidx,lvls) = (odata(486,jidx,lvls)+odata(489,jidx,lvls))/2
                       endif
                       do jidx = 590,640
                          odata(857:900,jidx,lvls) =  odata(857,jidx,lvls)
                          !odata(400:410,jidx,lvls) =  odata(410,jidx,lvls)
                          if (jidx<=597) then
                             odata(327:370,jidx,lvls) =  (odata(326,jidx,lvls)+odata(370,jidx,lvls))/2
                          else
                             odata(327:370,jidx,lvls) =  odata(370,jidx,lvls)
                          endif
                          odata(390:393,jidx,lvls) =  odata(390,jidx,lvls)
                          odata(766:793,jidx,lvls) = (odata(766,jidx,lvls)+odata(793,jidx,lvls))/2
                          odata(484:494,jidx,lvls) = (odata(484,jidx,lvls)+odata(494,jidx,lvls))/2
                          odata(387:410,jidx,lvls) = (odata(387,jidx,lvls)+odata(410,jidx,lvls))/2
                          odata(387:410,jidx,lvls) = (odata(387,jidx,lvls)+odata(410,jidx,lvls))/2
                       enddo
                       do iidx=1,80
                          odata(iidx,469:510,lvls)=odata(iidx,469,lvls)
                          odata(iidx,310:367,lvls)=odata(iidx,367,lvls)
                       enddo
                       if (lvls == 5) then
                          do iidx=3,28
                             odata(iidx,422:435,lvls) = (odata(iidx,421,lvls)+odata(iidx,436,lvls))/2
                          enddo
                       endif
                    enddo
                 endif
                 if (tracer(ii)=='SALT') then 
                 !    write(*,*) 'max/min tracer(',ii,') : ',tracer(ii),maxval(odata(:,:,kk)),minval(odata(:,:,kk))
                     where(odata(:,:,kk) <= 0._r8)
                          odata(:,:,kk) = 0._r8
                     endwhere
                  !   write(*,*) 'max/min tracer(',ii,') : ',tracer(ii),maxval(odata(:,:,kk)),minval(odata(:,:,kk))
                 endif
                 if (writeTestData) then  
                    open(22,file=trim(outFile),status='replace',access='direct', &
                                 form='unformatted',recl=ni*nj*8)
                    write(22,rec=1) odata(:,:,kk)
                    close(22)
                 endif
              enddo
              close(10)

              
              !save data in POP's flux correction formats
              !SST and SSS
              write(outFile,1010) trim(outFnames(ii)),yy+yearOffset,dayOfYear,hours(hh)
              outFile = trim(pathOut)//'/'//trim(outFile)
              write(*,'(A)') trim(outFile)
              open(10,file=trim(outFile),status='replace',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              write(10,rec=1) odata(:,:,1) !write surface data into the separate file
              close(10)
              !TEMP and SALT 3D
              write(outFile,1010) trim(outFnames(ii+2)),yy+yearOffset,dayOfYear,hours(hh)
              outFile = trim(pathOut)//'/'//trim(outFile)
              write(*,'(A)') trim(outFile)
              open(10,file=trim(outFile),status='replace',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              do kk=1,nk
                 write(10,rec=kk) odata(:,:,kk) !write 3D data 
              enddo                
           enddo
!           stop
        enddo
        dayOfYear = dayOfYear + 1
        
     enddo !dd
   enddo
 enddo

 stop

!old part of the script
 open(10,file='init_ts_bs01v1_09052018.ieeer8',status='replace',access='direct', &
         form='unformatted',recl=ni*nj*8)
! open(5,file='ic_bs2v3_v3.txt',status='old')
 open(11,file='temp_in',access='direct',status='old',form='unformatted',recl=ni*nj*8)
 jj=1
 do ii=1,nk
    read(11,rec=ii) odata
    write(10,rec=jj) odata
    write(*,'(2i4,2f14.6)') ii,jj,maxval(odata),minval(odata)
    jj=jj+1
 enddo
 close(11)
 open(11,file='salt_in',access='direct',status='old',form='unformatted',recl=ni*nj*8)
 do ii=1,nk
    read(11,rec=ii) odata
    odata=odata*1000.
    write(10,rec=jj) odata
    write(*,'(2i4,2f14.6)') ii,jj,maxval(odata),minval(odata)
    jj=jj+1
 enddo
 close(11)
 odata=0.
 do ii=1,nk
    write(10,rec=jj) odata
    write(*,'(2i4,2f14.6)') ii,jj,maxval(odata),minval(odata)
    jj=jj+1
 enddo
 stop
 do ii=1,nk
    odata=t(ii)
    write(10,rec=ii) odata
    write(*,*) ii,maxval(odata),minval(odata)
 enddo
 do ii=1,nk
    odata=s(ii)
    write(10,rec=ii+nk) odata
    write(*,*) ii+nk,maxval(odata),minval(odata)
 enddo
 odata = 0.
 do ii=1,nk
    write(10,rec=ii+nk+nk) odata
    write(*,*) ii+nk+nk,maxval(odata),minval(odata)
 enddo
 close(10)
 close(5)
 end program
    
