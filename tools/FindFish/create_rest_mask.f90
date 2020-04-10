 program create_restoring_mask_bs05v1

 integer,parameter :: ni=1000,nj=640,nk=33,r8=selected_real_kind(13)
 real*8 :: t(nk),s(nk),odata(ni,nj),mask(ni,nj),srm,otmp(ni,nj),dN,dp,real_kmt_mask(ni,nj)
 integer :: ii,jj,kmt(ni,nj),imask(ni,nj),maskPoints,ib,jb,npoints
 logical :: res
 character*256 :: kmtFile

  kmtFile = 'kmt.bs05v1.ocn.20170627.ieeei4'
!read in kmt file
  open(22,file=trim('kmt_file'),access='direct',form='unformatted', &
               status='old',recl=ni*nj*4)
  read(22,rec=1) imask
  close(22)
  real_kmt_mask = 0._r8
  where (imask > 0) real_kmt_mask = 1._r8



 jb = 250
 ib = 1
 mask = 0._r8
 maskPoints=50
 mask(ib:maskPoints,jb:nj) = 1._r8
 mask(:,nj-maskPoints+1:nj) = 1._r8
 npoints=20
 dp=real(npoints,r8)
 dp=1._r8/dp
 do ii=1,npoints
    mask(maskPoints+ii,jb:nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,'|',jb,':',640-maskPoints-ii,'|||',1._r8-real((ii-1),r8)*dp
 enddo
 write(*,*) '======================='
 do ii=1,npoints
    mask(maskPoints+ii:ni,nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,':',1000,'|',jb,':',640-maskPoints-1-ii,'|||',1._r8-real((ii-1),r8)*0.1_r8
 enddo

 mask = mask * real_kmt_mask
 mask(50:100,300:325)=0._r8
 mask(150:327,550:640)=0._r8
 mask(1:150,220:300)=0._r8

 open(10,file='tsrst_mask_ff_15032020.bin',access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10)

 
  
     



 end program
    
