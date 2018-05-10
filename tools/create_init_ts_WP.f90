 program create_init_ts_bs01v1

 integer,parameter :: ni=1000,nj=640,nk=33
 real*8 :: t(nk),s(nk),odata(ni,nj)
 integer :: ii,jj
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
    
