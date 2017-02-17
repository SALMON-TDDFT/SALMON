subroutine read_input_gceed(procid,cfunction2)

  implicit none
  include 'mpif.h'
  character(30),intent(out) :: cfunction2
  integer :: procid
  integer :: ierr
  namelist / group_function2 / cfunction2

  if(procid==0)then
    read(*,nml=group_function2)
  end if
  call mpi_bcast(cfunction2,30,mpi_character,0,mpi_comm_world,ierr)

end subroutine read_input_gceed
