 program create_init_ts_bs01v1

 integer,parameter :: ni=1000,nj=640,nk=33,r8=selected_real_kind(13)
 real*8 :: t(nk),s(nk),odata(ni,nj),mask(ni,nj),srm,otmp(ni,nj),dN,dp
 integer :: ii,jj,kmt(ni,nj),kmtNew(ni,nj),imask(ni,nj),maskPoints,ib,jb,npoints
 logical :: res
 character*256 :: slvlFile,rstFile,dataPath, kmtFile
 integer :: slvlMaskPoints,rstMaskPoints,slvlNPoints,rstNPoints 
 slvlFile = 'slvl_mask.bin'
 rstFile  = 'tsrst_mask.bin'
 dataPath = '/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/grid'
! kmt file /scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs01v1/grid/kmt.bs01v1.ocn.20180432.ieeer4 
 kmtFile = 'kmt.bs01v1.ocn.20180432.ieeer4'
 slvlMaskPoints=20
 rstMaskPoints =40
 slvlNPoints = 40
 rstNPoints = 30 
 jb = 510
 ib = 1
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
    
