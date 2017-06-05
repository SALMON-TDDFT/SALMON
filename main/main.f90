program main
  use salmon_global
  use inputoutput
  implicit none
  integer :: nprocs,myrank

  call setup_parallel(nprocs,myrank)
  call read_stdin(myrank)
  call read_input_common(myrank)
  call dump_input_common(myrank)
! read_stdin, read_input_common and dump_input_common will be wrapped by 
! a single routine routine after implimentation for MPI wrapper

  select case(iperiodic)
  case(0)
    call gceed(nprocs,myrank)
  case(3)
    call arted(nprocs, myrank)
  case default
    stop 'invalid iperiodic'
  end select

  call end_parallel

end program main
