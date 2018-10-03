 program create_init_ts_bs01v1

 integer,parameter :: ni=1000,nj=640,mnths=12,r8=selected_real_kind(13)
!ustawienia*
 integer, parameter :: nhours = 4,ntracers=1,nk=1
! real(r8), parameter :: mean_sea_level = 29.024161 ![cm] average from one point (x=500,y=600)
 real(r8), parameter :: mean_sea_level = 29.5161514642064 ![cm] average from whole domain of Water Puck
!***********
 real*8 :: t(nk),s(nk),odata(ni,nj,nk),factors(ntracers),mxdata,lmask(ni,nj)
 integer :: ii,jj,kk,ll,daysInMonths(mnths),year,month,day,hour, &
            sDate(4),eDate(4),hours(nhours),sHH,eHH,seconds(nhours)
 character*256 :: pathIn,pathOut
 character*4 :: tracer(ntracers)
 integer :: yy,mm,dd,hh,dayOfYear,midx,kmt(ni,nj)
 character(len=8),dimension(1) :: outFNames
 character*256 :: inFile,outFile, kmtFile

 kmtFile = '/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/grid/kmt.bs01v1.ocn.20180432.ieeer4'

!ustawienia*
 data sDate /2013, 1, 1, 1/ !data poczatkowa 
 data eDate /2013,12,31, 4/ !data koncowa
 data seconds /3600,25200,46800,68400/
! data hours /0,6,12,18/
 data hours /3,9,15,21/
 data outFnames /"SSH"/

! pathIn  = '../../tmp_data/ARTUR/lbc/2011'
 2000 format('../../tmp_data/ARTUR/lbc/',i4.4)
 write(pathIn,2000) sDate(1)
 pathIn = trim(pathIn)
 pathOut = '/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/forcing'
! for testing only pathOut = '.' 
 1000 format(i4.4,'-',i2.2,'-',i2.2,'-',i5.5,'_',a3,'_1000_0640_0001_0001.ieeer8') 
!2011/2011-01-01-03600_SSH_1000_0640_0001_0001.ieeer8 
 1010 format(A,'.',i4.4,'.',i3.3,'.',i2.2)
!***********
!read in kmt file
 open(10,file=trim(kmtFile),status='old',access='direct', &
                           form='unformatted',recl=ni*nj*4) 
 read(10,rec=1) kmt
 close(10)
 

 lmask=0._r8
 where (kmt >0) lmask=1._r8 
 data tracer /'SSH'/
 data factors /0.01/
 data daysInMonths /31,28,31,30,31,30,31,31,30,31,30,31/
 
 midx=0
 mxdata = 0.
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
!              stop
              open(10,file=trim(inFile),status='old',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              do kk=1,nk
                 read(10,rec=kk) odata(:,:,kk)
                 odata(:,:,kk) = (odata(:,:,kk)-mean_sea_level)*factors(ii)
                 midx=midx+1
                 mxdata=mxdata+sum(odata(:,:,kk)*lmask) 
                 write(*,'(i4,3f14.6,"   ",A)') kk,maxval(odata(:,:,kk)),minval(odata(:,:,kk)),& 
                          mxdata/sum(lmask)/real(midx,r8),tracer(ii) 
!                 write(outFile,"(a4,i2.2,'.bin')") tracer(ii),kk
!                 write(*,'(A)') outFile
!                 open(22,file=trim(outFile),status='replace',access='direct', &
!                           form='unformatted',recl=ni*nj*8)
!                 write(22,rec=1) odata(:,:,kk)
!                 close(22)
              enddo
              close(10)
              
              !save data in POP's flux correction formats
              !SST and SSS
              write(outFile,1010) trim(outFnames(ii)),yy,dayOfYear,hours(hh)
              outFile = trim(pathOut)//'/'//trim(outFile)
              write(*,'(A)') trim(outFile)
              open(10,file=trim(outFile),status='replace',access='direct', &
                           form='unformatted',recl=ni*nj*8)
              write(10,rec=1) odata(:,:,1) !write surface data into the separate file
              close(10)
!for SSH is not needed              !TEMP and SALT 3D
!              write(outFile,1010) trim(outFnames(ii+2)),yy,dayOfYear,hours(hh)
!              outFile = trim(pathOut)//'/'//trim(outFile)
!              write(*,'(A)') trim(outFile)
!              open(10,file=trim(outFile),status='replace',access='direct', &
!                           form='unformatted',recl=ni*nj*8)
!              do kk=1,nk
!                 write(10,rec=kk) odata(:,:,kk) !write 3D data 
!              enddo                
           enddo
!           stop
        enddo
        dayOfYear = dayOfYear + 1
     enddo !dd
   enddo
 enddo
 write(*,*) 'mean slevel = ',mxdata/sum(lmask)/real(midx,r8)

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
    
