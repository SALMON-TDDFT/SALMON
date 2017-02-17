subroutine gceed(nprocs,procid)

  implicit none
  integer :: nprocs,procid
  character(30) :: cfunction2

  call read_input_gceed(procid,cfunction2)

  if(cfunction2=="scf")then
    call real_space_dft(nprocs,procid)
  else if(cfunction2=="rt")then
    call real_time_dft(nprocs,procid)
  end if

end subroutine gceed
