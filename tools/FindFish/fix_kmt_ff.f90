  program read_grid

  implicit none
  integer, parameter :: char_len  = 256, &
                            int_kind  = kind(1),  &
                            log_kind  = kind(.true.), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13)

  integer (kind=int_kind), parameter :: imt = 1000, &
!       jmt = 640, kmt = 66
        jmt = 640, kmt = 26, &
        ic0 = 0

  integer(selected_int_kind(12)) :: polei, &
      polej, idist, jdist,i,j,k,in3,nrecl,ii,jj,kk,nrec
  integer(int_kind) :: bath(imt,jmt),bath0(imt,jmt)
  real (dbl_kind),parameter :: c1=1, c4=4
  character(char_len) :: inFile,outFile

        inFile = 'kmt_file'
        outFile= 'FF_kmt_out.bin'
 2000 format(32x,i2)

      in3=22
      nrecl=imt*jmt*4
      open(in3,file=trim(inFile), &
       form='unformatted',access='direct',recl=nrecl, convert='big_endian', status='old')
      write(*,*) 'infile:',trim(inFile)
      read(in3,rec=1) bath0
    close(in3) 
    bath = bath0
    write(*,*) 'max/min bath = ',maxval(bath),minval(bath)
! set edges
    bath(1:3,:)=ic0
    bath(imt-2:imt,:)=ic0
    bath(:,1:3)=ic0
    bath(:,jmt-2:jmt)=ic0 
! set boundaries
    bath(1:50,316:510) = 5
    bath(51:80,330:510) = 5
    bath(327:900,590:640) = 5
    where (bath0 <=0) bath = 0
         open(in3,file=trim(outFile), &
       form='unformatted',access='direct',recl=nrecl, convert='big_endian', status='replace')
      write(*,*) 'outfile:',trim(outFile)
      write(in3,rec=1) bath
    close(in3)

end program
