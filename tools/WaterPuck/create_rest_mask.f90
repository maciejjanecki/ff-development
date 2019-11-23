 program create_init_ts_bs01v1

 integer,parameter :: ni=1000,nj=640,nk=33,r8=selected_real_kind(13)
 real*8 :: t(nk),s(nk),odata(ni,nj),mask(ni,nj),srm,otmp(ni,nj),dN,dp
 integer :: ii,jj,kmt(ni,nj),imask(ni,nj),maskPoints,ib,jb,npoints
 logical :: res
 jb = 530
 ib = 1
 mask = 0._r8
 maskPoints=30
 mask(ib:maskPoints,jb:nj) = 1._r8
 mask(:,nj-maskPoints+1:nj) = 1._r8
 npoints=20
 dp=real(npoints,r8)
 do ii=1,npoints
    mask(maskPoints+ii,jb:nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,'|',jb,':',640-maskPoints-ii,'|||',1._r8-real((ii-1),r8)*dp
 enddo
 write(*,*) '======================='
 do ii=1,npoints
    mask(maskPoints+ii:ni,nj-maskPoints-1-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,':',1000,'|',jb,':',640-maskPoints-1-ii,'|||',1._r8-real((ii-1),r8)*0.1_r8
 enddo

 open(10,file='tsrst_mask.bin',access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10)

 
  
     



 end program
    
