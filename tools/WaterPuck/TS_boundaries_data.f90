 program TS_boundaries_data

 integer,parameter :: ni=1000,nj=640,nk=33,mnths=12
!ustawienia*
 integer, parameter :: nhours = 4,ntracers=2,yearOffset = 0
 
!***********
 real*8 :: t(nk),s(nk),odata(ni,nj,nk),factors(ntracers)
 integer :: ii,jj,kk,ll,daysInMonths(mnths),year,month,day,hour, &
            sDate(4),eDate(4),hours(nhours),sHH,eHH,seconds(nhours)
 character*256 :: pathIn,pathOut
 character*4 :: tracer(ntracers)
 integer :: yy,mm,dd,hh,dayOfYear
 character(len=8),dimension(4) :: outFNames
 character*256 :: inFile,outFile
 character*4   ::sDateC,eDateC

!ustawienia*
 data sDate /2019, 1, 1, 1/ !data poczatkowa 
 data eDate /2019,12,31, 4/ !data koncowa
 data seconds /3600,25200,46800,68400/
! data hours /0,6,12,18/
 data hours /3,9,15,21/
 data outFnames /"SST","SSS","TEMP3D","SALT3D"/

 write(sDateC,'(i4.4)') sDate(1)
 write(eDateC,'(i4.4)') eDate(1)

 pathIn  = '../../../tmp_data/ARTUR/lbc/'//sDateC
 pathOut = '/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/forcing' 
 1000 format(i4.4,'-',i2.2,'-',i2.2,'-',i5.5,'_',a4,'_1000_0640_0033_0001.ieeer8') 
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
              write(inFile,1000) yy,mm,dd,seconds(hh),tracer(ii)
              inFile = trim(pathIn)//'/'//trim(inFile)
              write(*,'(A)') trim(inFile)
              open(10,file=trim(inFile),status='old',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              do kk=1,nk
                 read(10,rec=kk) odata(:,:,kk)
                 odata(:,:,kk) = odata(:,:,kk)*factors(ii)
!                 write(*,'(i4,2f14.6,"   ",A)') kk,maxval(odata(:,:,kk)),minval(odata(:,:,kk)),tracer(ii) 
                 write(outFile,"(a4,i2.2,'.bin')") tracer(ii),kk
!                 write(*,'(A)') outFile
                 open(22,file=trim(outFile),status='replace',access='direct', &
                           form='unformatted',recl=ni*nj*8)
                 write(22,rec=1) odata(:,:,kk)
                 close(22)
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
    
