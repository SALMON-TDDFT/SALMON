subroutine inner_product5(matbox1,matbox2,cbox)
  use salmon_parallel, only: nproc_group_korbital
  use salmon_communication, only: comm_summation
  use misc_routines, only: get_wtime
  use scf_data
  use new_world_sub
  !$ use omp_lib
  implicit none
  integer :: iob,iob_allob
  integer :: ix,iy,iz
  complex(8) :: matbox1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum)
  complex(8) :: matbox2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum)
  complex(8) :: cbox(1:itotMST),cbox2(1:itotMST)
  complex(8) :: sum0
  
  cbox2(:)=0.d0
  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
    sum0=0.d0
    !$OMP parallel do collapse(2) reduction(+ : sum0)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      sum0=sum0+conjg(matbox1(ix,iy,iz,iob))*matbox2(ix,iy,iz,iob)
    end do
    end do
    end do
    cbox2(iob_allob)=sum0*Hvol
  end do
  call comm_summation(cbox2,cbox,itotMST,nproc_group_korbital)

end subroutine inner_product5
