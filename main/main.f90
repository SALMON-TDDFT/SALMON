program main

  implicit none
  integer :: nprocs,procid
  character(30) :: cfunction

  call setup_parallel(nprocs,procid)
  call read_input(procid,cfunction)

  if(cfunction=="nanostructure") then
    call gceed(nprocs,procid)
  end if

  call end_parallel

end program main
