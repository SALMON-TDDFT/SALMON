subroutine setup_parallel(nprocs,procid)

  implicit none
  include 'mpif.h'
  integer,intent(out) :: nprocs,procid
  integer :: ierr

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  call mPI_comm_rank(mpi_comm_world,procid,ierr)

end subroutine setup_parallel
