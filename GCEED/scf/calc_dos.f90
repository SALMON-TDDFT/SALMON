subroutine calc_dos
use scf_data
use allocate_psl_sub
use new_world_sub
implicit none
integer :: iob,iobmax,iob_allob,iene
real(8) :: rbox_dos(-300:300)
real(8) :: dos(-300:300)
real(8),parameter :: sigma_gd=0.01d0

call calc_pmax(iobmax)

rbox_dos=0.d0

do iob=1,iobmax
  call calc_allob(iob,iob_allob)
  do iene=-300,300
    rbox_dos(iene)=rbox_dos(iene)  &
              +exp(-(dble(iene)/10d0/2.d0/Ry-esp(iob_allob,1))**2/(2.d0*sigma_gd**2))/sqrt(2.d0*Pi*sigma_gd**2)
  end do
end do
call MPI_Allreduce(rbox_dos,dos,601,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_grid,ierr) 

if(myrank==0)then
  open(101,file="dos.data")
  do iene=-300,300
    write(101,'(f10.5,f14.8)') dble(iene)/10.d0,dos(iene)
  end do
  close(101)
end if

end subroutine calc_dos
