subroutine setup_parallel(nprocs,myrank)
  use mpi

  implicit none
  integer,intent(out) :: nprocs,myrank
  integer :: ierr

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  call mPI_comm_rank(mpi_comm_world,myrank,ierr)

end subroutine setup_parallel
