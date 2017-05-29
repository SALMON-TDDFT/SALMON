program main
  use input
  implicit none
  integer :: nprocs,myrank
  character(30) :: cfunction

  call setup_parallel(nprocs,myrank)
  call read_stdin(myrank,cfunction)

  select case(cfunction)
  case("nanostructure")
    call gceed(nprocs,myrank)
  case("singlecell", "multiscale")
    call arted(nprocs, myrank, cfunction)
  end select

  call end_parallel

end program main
