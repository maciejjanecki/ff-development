 program create_sea_level_mask !find fish

 integer,parameter :: ni=1000,nj=640,nk=33,r8=selected_real_kind(13)
 real*8 :: t(nk),s(nk),odata(ni,nj),mask(ni,nj),srm,otmp(ni,nj),dN,dp,real_kmt_mask(ni,nj)
 integer :: ii,jj,kmt(ni,nj),kmtNew(ni,nj),imask(ni,nj),maskPoints,ib,jb,npoints
 logical :: res
 character*256 :: slvlFile,rstFile,dataPath, kmtFile
 integer :: slvlMaskPoints,rstMaskPoints,slvlNPoints,rstNPoints 
 slvlFile = 'slvl_mask_ff_06032020.bin'
 rstFile  = 'tsrst_mask_ff_06032020.bin'
 dataPath = '/users/work/ffish/cesm_input_data/ocn/pop/bs05v1/grid'
! kmt file /scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs05v1/grid/kmt.bs05v1.ocn.20170627.ieeei4 
 kmtFile = 'kmt.bs05v1.ocn.20170627.ieeei4'
!read in kmt file
  open(22,file=trim('kmt_file'),access='direct',form='unformatted', &
               status='old',recl=ni*nj*4)
  read(22,rec=1) imask
  close(22)
  real_kmt_mask = 0._r8
  where (imask > 0) real_kmt_mask = 1._r8
 ! slvlMaskPoints=34  !1st version
 slvlMaskPoints=60
 rstMaskPoints =60
 ! slvlNPoints = 44
 slvlNPoints = 60
 rstNPoints = 30 
 jb = 300
 ib = 3
 mask = 0._r8
 maskPoints=slvlMaskPoints
 npoints=slvlNPoints
 mask(ib:maskPoints,jb:nj) = 1._r8
 mask(:,nj-maskPoints-1:nj) = 1._r8
 dp=1._r8/real(npoints,r8)
 do ii=1,npoints
    mask(maskPoints+ii,jb:nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,'|',jb,':',nj-maskPoints-ii,'| |',1._r8-real((ii-1),r8)*dp
 enddo
 write(*,*) '======================='
 do ii=1,npoints
    mask(maskPoints+ii:ni,nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,':',ni,'|',jb,':',nj-maskPoints-1-ii,'|||',1._r8-real((ii-1),r8)*dp
 enddo
 !extend mask into whole domain
 where (mask <= 0.25_r8) mask=0.25_r8
 mask = mask * real_kmt_mask
 mask(50:100,300:325)=0._r8
 mask(150:327,550:640)=0._r8
 write(*,*) trim(slvlFile)
 open(10,file=trim(slvlFile),access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10)
 slvlFile = trim(dataPath)//'/'//trim(slvlFile)
 open(10,file=trim(slvlFile),access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10)
 write(*,*) trim(slvlFile)
 stop
!===== restoring
!clear mask
 mask = 0._r8
 maskPoints=rstMaskPoints
 npoints=rstNPoints
 mask(ib:maskPoints,jb:nj) = 1._r8
 mask(:,nj-maskPoints-1:nj) = 1._r8
 dp=1._r8/real(npoints,r8)
 do ii=1,npoints
    mask(maskPoints+ii,jb:nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,'|',jb,':',nj-maskPoints-ii,'| |',1._r8-real((ii-1),r8)*dp
 enddo
 write(*,*) '======================='
 do ii=1,npoints
    mask(maskPoints+ii:ni,nj-maskPoints-ii) = 1._r8-real((ii-1),r8)*dp
    write(*,*) maskPoints+ii,':',ni,'|',jb,':',nj-maskPoints-1-ii,'|||',1._r8-real((ii-1),r8)*dp
 enddo
 write(*,*) trim(rstFile)
 open(10,file=trim(rstFile),access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10)
 rstFile = trim(dataPath)//'/'//trim(rstFile)
 open(10,file=trim(rstFile),access='direct',form='unformatted', &
               status='replace',recl=ni*nj*8)
 write(10,rec=1) mask
 close(10) 
 write(*,*) trim(rstFile)

 write(*,*) '============== kmt =============='
 kmtFile=trim(dataPath)//'/'//trim(kmtFile)
 open(10,file = trim(kmtFile),status = 'old',access = 'direct', &
      form = 'unformatted',recl=ni*nj*4)
 read(10,rec=1) kmt
 kmtNew = kmt
 imask = 0
 where (kmt > 0) imask = 1
 write(*,*) 'kmt max/min = ',maxval(kmt),minval(kmt)
 close(10)
 write(*,'(A)') trim(kmtFile)
 where(mask >= 0.99_r8) kmtNew=15
 kmtNew = kmtNew*imask
 kmtFile = 'kmt.bs01v1.ocn.20180602.ieeer4'
 kmtFile=trim(dataPath)//'/'//trim(kmtFile)
  open(10,file = trim(kmtFile),status = 'replace',access = 'direct', &
      form = 'unformatted',recl=ni*nj*4)
 write(10,rec=1) kmtNew
 close(10)
 write(*,'(A)') trim(kmtFile)
 kmtFile = 'kmt.bs01v1.ocn.20180602.ieeer4'
 open(10,file = trim(kmtFile),status = 'replace',access = 'direct', &
      form = 'unformatted',recl=ni*nj*4)
 write(10,rec=1) kmt!  New
 close(10)
 write(*,'(A)') trim(kmtFile)

 

 
  
     



 end program
    
