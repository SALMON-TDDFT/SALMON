subroutine add_in_broyden(trho,trho_in,tfb,trho_temp)
  use inputoutput
  use salmon_parallel, only: nproc_group_orbital
  use salmon_communication, only: comm_summation
  use scf_data
  implicit none
  integer :: ix,iy,iz
  integer :: nn,nn2
  real(8) :: trho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: trho_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: tfb(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: trho_temp(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  
  trho(:,:,:)=trho_in(:,:,:)+alpha_mb*tfb(:,:,:)-trho_temp(:,:,:)
  
  nn=0
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    if(trho(ix,iy,iz) < 0.d0) then
      nn=nn+1
    end if
  end do
  end do
  end do
  call comm_summation(nn,nn2,nproc_group_orbital)
  if(nn2 > 0) then
    trho(:,:,:)=trho_in(:,:,:)+alpha_mb*tfb(:,:,:)
  end if
  
  return

end subroutine add_in_broyden
