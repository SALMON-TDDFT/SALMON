program main
  use salmon_global
  use salmon_parallel
  use inputoutput
  implicit none

  call setup_parallel
  call read_stdin(nproc_id_global)
  call read_input_common(nproc_id_global)
  call dump_input_common(nproc_id_global)
! read_stdin, read_input_common and dump_input_common will be wrapped by 
! a single routine routine after implimentation for MPI wrapper

  select case(iperiodic)
  case(0)
    call gceed(nproc_size_global, nproc_id_global)
  case(3)
    call arted(nproc_size_global, nproc_id_global)
  case default
    stop 'invalid iperiodic'
  end select

  call end_parallel

end program main
