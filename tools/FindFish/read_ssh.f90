  program read_ssh

  implicit none
  integer, parameter :: char_len  = 256, &
                            int_kind  = kind(1),  &
                            log_kind  = kind(.true.), &
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13)

  integer (kind=int_kind), parameter :: imt = 1000, &
!       jmt = 640, kmt = 66
        jmt = 640, kmt = 33

  integer(selected_int_kind(12)) :: polei, &
      polej, idist, jdist,i,j,k,in3,nrecl,ii,jj,kk,nrec
  real (dbl_kind),parameter :: c1=1, c4=4
  real(dbl_kind) :: rtmp,rstData(imt,jmt)
  logical,parameter :: writeParameters = .true.
        

        character(char_len) :: in_file,outFile

      in_file = '/scratch/lustre/plgjjakacki/LD/cesm_input_data/ocn/pop/bs05v1/forcing/SSH.2011.002.03'
      in3=22
      nrecl = imt*jmt*8
      open(in3,file=trim(in_file), &
       form='unformatted',access='direct',recl=nrecl, convert='big_endian', status='old')
      write(*,'(A,A)') 'infile:',trim(in_file)
      nrec=1
      read(in3,rec=nrec) rstData
      write(*,'(3i6,A,2f14.6)') i,j,nrec,'   max/min = ',maxval(rstData),minval(rstData)
    close(in3) 
    
end program
