program main

  implicit none
  integer :: nprocs,myrank
  character(30) :: cfunction

  call setup_parallel(nprocs,myrank)
  call read_input(myrank,cfunction)

  if(cfunction=="nanostructure") then
    call gceed(nprocs,myrank)
  end if

  call end_parallel

end program main
