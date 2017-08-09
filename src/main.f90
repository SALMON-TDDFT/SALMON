program main
  use salmon_global
  use salmon_parallel
  use inputoutput
  implicit none

  call setup_parallel
  call read_input

  select case(iperiodic)
  case(0)
    call gceed(nproc_size_global, nproc_id_global)
  case(3)
    call arted
  case default
    stop 'invalid iperiodic'
  end select

  call end_parallel

end program main
