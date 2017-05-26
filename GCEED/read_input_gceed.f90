subroutine read_input_gceed(procid,cfunction2)
  use inputoutput
  implicit none
  include 'mpif.h'
  character(30),intent(out) :: cfunction2
  integer :: procid
  integer :: ierr
  namelist / group_function2 / cfunction2

  if(procid==0)then
    open(fh_namelist, file='.namelist.tmp', status='old')
    read(fh_namelist,nml=group_function2)
    close(fh_namelist)
  end if
  call mpi_bcast(cfunction2,30,mpi_character,0,mpi_comm_world,ierr)

end subroutine read_input_gceed
