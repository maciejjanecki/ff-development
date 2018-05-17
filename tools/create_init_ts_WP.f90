 program create_init_ts_bs01v1

 integer,parameter :: ni=1000,nj=640,nk=33,r8=selected_real_kind(13)
 real*8 :: t(nk),s(nk),odata(ni,nj),rmask(ni,nj),srm,otmp(ni,nj),dN
 integer :: ii,jj,kmt(ni,nj),imask(ni,nj)
 logical :: res
 open(10,file='init_ts_bs01v1_09052018.ieeer8',status='replace',access='direct', &
         form='unformatted',recl=ni*nj*8)
! open(5,file='ic_bs2v3_v3.txt',status='old')
 open(11,file='temp_in',access='direct',status='old',form='unformatted',recl=ni*nj*8)
 
!read in kmt file and fill all land values
 open(100,file='/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/grid/kmt.bs01v1.ocn.20180432.ieeer4', &
         access='direct',status='old',form='unformatted',recl=ni*nj*4)
 read(100,rec=1) kmt
 write(*,*) 'max/min kmt = ',maxval(kmt),minval(kmt)
 close(100)
 imask = 0
 where (kmt > 0) imask = 1
 rmask=real(imask,r8)
 srm = sum(rmask)
!check that NaN test works
 dN = 1/0.
 res=isnan(dN)
 write(*,*) res
 if (dN .ne.dN) then 
   write(*,*) 'NaN works'
 else
   
   stop ('NaN does not work')
 endif
 jj=1
 do ii=1,nk
    read(11,rec=ii) odata
    otmp = 0._r8
    where (imask == 1) otmp = odata
    tmp=sum(odata)/srm
    odata = otmp
    where (imask == 0) odata = tmp 
    write(10,rec=jj) odata
    write(*,'(2i4,2f14.6)') ii,jj,maxval(odata),minval(odata)
    jj=jj+1
 enddo
 close(11)
 open(11,file='salt_in',access='direct',status='old',form='unformatted',recl=ni*nj*8)
 do ii=1,nk
    read(11,rec=ii) odata
    odata=odata*1000.
    otmp = 0._r8
    where (imask == 1) otmp = odata
    tmp=sum(odata)/srm
    odata = otmp
    where (imask == 0) odata = tmp
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
    
