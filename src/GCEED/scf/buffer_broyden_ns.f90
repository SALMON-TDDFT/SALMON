subroutine buffer_broyden_ns(iter)
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use scf_data, only: ng_sta,ng_end,ng_num,num_rho_stock,rho,rho_in,rho_out,  &
                      ilsda,rho_s,rho_s_in,rho_s_out
  integer,intent(in) :: iter
  real(8) :: vecr(1:ng_num(1)*ng_num(2)*ng_num(3))
  real(8) :: vecr_in(1:ng_num(1)*ng_num(2)*ng_num(3),num_rho_stock+1)
  real(8) :: vecr_out(1:ng_num(1)*ng_num(2)*ng_num(3),num_rho_stock+1)
  integer :: i
  integer :: iix

  if(ilsda==0)then

    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
      vecr(iix)=rho(ix,iy,iz)
    end do
    end do
    end do

    do i=1,num_rho_stock
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
        vecr_in(iix,i)=rho_in(ix,iy,iz,i)
        vecr_out(iix,i)=rho_out(ix,iy,iz,i)
      end do
      end do
      end do
    end do

    call broyden(vecr,vecr_in,vecr_out,ng_num(1)*ng_num(2)*ng_num(3),  &
                 iter,num_rho_stock,num_rho_stock)
  
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
      rho(ix,iy,iz)= vecr(iix)
    end do
    end do
    end do

    do i=1,num_rho_stock+1
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
        rho_in(ix,iy,iz,i)=vecr_in(iix,i)
        rho_out(ix,iy,iz,i)=vecr_out(iix,i)
      end do
      end do
      end do
    end do

  else if(ilsda==1)then
    
    do is=1,2
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
        vecr(iix)=rho_s(ix,iy,iz,is)
      end do
      end do
      end do

      do i=1,num_rho_stock
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
          vecr_in(iix,i)=rho_s_in(ix,iy,iz,is,i)
          vecr_out(iix,i)=rho_s_out(ix,iy,iz,is,i)
        end do
        end do
        end do
      end do
 
      call broyden(vecr,vecr_in,vecr_out,ng_num(1)*ng_num(2)*ng_num(3),   &
                   iter,num_rho_stock,num_rho_stock)
  
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
        rho_s(ix,iy,iz,is)= vecr(iix)
      end do
      end do
      end do
      
      do i=num_rho_stock,num_rho_stock+1
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          iix=(ix-ng_sta(1)+1)+(iy-ng_sta(2))*ng_num(1)+(iz-ng_sta(3))*ng_num(1)*ng_num(2)
          rho_s_in(ix,iy,iz,is,i)=vecr_in(iix,i)
          rho_s_out(ix,iy,iz,is,i)=vecr_out(iix,i)
        end do
        end do
        end do
      end do
    end do

  end if

end subroutine buffer_broyden_ns
