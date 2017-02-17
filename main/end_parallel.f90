subroutine end_parallel

  implicit none
  include 'mpif.h'
  integer :: ierr

  call mpi_finalize(ierr)

end subroutine end_parallel