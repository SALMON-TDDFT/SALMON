subroutine read_input(myrank,cfunction)

  implicit none
  include 'mpif.h'
  character(30),intent(out) :: cfunction
  integer :: ierr
  integer :: myrank
  namelist / group_function / cfunction

  if(myrank==0)then
    read(*,nml=group_function)
  end if
  call mpi_bcast(cfunction,30,mpi_character,0,mpi_comm_world,ierr)

end subroutine read_input
