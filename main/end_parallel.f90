subroutine end_parallel
  use salmon_communication
  implicit none

  call comm_finalize
end subroutine end_parallel
