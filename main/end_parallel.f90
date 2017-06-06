subroutine end_parallel
  use mpi

  implicit none
  integer :: ierr

  call mpi_finalize(ierr)

end subroutine end_parallel
