subroutine setup_parallel(nprocs,myrank)
subroutine setup_parallel
  use salmon_parallel
  use salmon_communication
  implicit none

  call comm_init
  call comm_get_globalinfo(nproc_group_global, nproc_id_global, nproc_size_global)
end subroutine setup_parallel
