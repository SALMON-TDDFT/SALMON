subroutine broyden(iter)
  use inputoutput
  use salmon_parallel, only: nproc_group_global,nproc_group_orbital
  use salmon_communication, only: comm_summation
  use scf_data
  implicit none
  integer,parameter :: iter_mb=0
  real(8),parameter :: romega0=0.01d0
  integer :: iter_s,iter_e
  integer :: iter,i,j,ix,iy,iz
  integer :: is,is_sta,is_end
  integer :: ibox
  real(8),allocatable :: fb(:,:,:,:)
  real(8),allocatable :: del_fb(:,:,:,:),del_x(:,:,:,:)
  real(8),allocatable :: romega_mb(:)
  real(8) :: s,rho_temp(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8),allocatable :: aa(:,:),beta(:,:)
  real(8) :: rbox,rbox2
  
  if(ilsda==0)then
    is_sta=1
    is_end=1
  else if(ilsda==1)then
    is_sta=1
    is_end=2
  end if
  
  if (iter <= iter_mb+1) then
    if(ilsda==0)then
      rho(:,:,:)=rho_in(:,:,:,num_rho_stock)+alpha_mb*(rho_out(:,:,:,num_rho_stock)-rho_in(:,:,:,num_rho_stock))
    else
      rho_s(:,:,:,1:2)=rho_s_in(:,:,:,1:2,num_rho_stock)  &
                         +alpha_mb*(rho_s_out(:,:,:,1:2,num_rho_stock)-rho_s_in(:,:,:,1:2,num_rho_stock))
    end if
  else
    iter_s=max(iter_MB+1,iter-nmemory_mb)
    iter_e=iter-1
    allocate(fb(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iter_s:iter))
    allocate(del_fb(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iter_s:iter_e))
    allocate(del_x(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iter_s:iter_e))
    allocate(romega_mb(iter_s:iter_e))
    allocate(beta(iter_s:iter_e,iter_s:iter_e))
    allocate(aa(1:iter_e-iter_s+1,1:iter_e-iter_s+1))
    do is=is_sta,is_end
      if(ilsda==0)then
        do i=iter_s,iter
          ibox=num_rho_stock-(iter-i)
          fb(:,:,:,i)=rho_out(:,:,:,ibox)-rho_in(:,:,:,ibox)
        end do
      else if(ilsda==1)then
        do i=iter_s,iter
          ibox=num_rho_stock-(iter-i)
          fb(:,:,:,i)=rho_s_out(:,:,:,is,ibox)-rho_s_in(:,:,:,is,ibox)
        end do
      end if
      romega_mb(iter_s:iter_e)=1.d0
  !     &        /sqrt(sum((rho_out_b(:,iter_mb+1:iter-1)-rho_in_b(:,1:iter_mb+1:iter-1))**2))
      if(ilsda==0)then
        do i=iter_s,iter_e
          ibox=num_rho_stock-(iter-i)
          del_x(:,:,:,i)=rho_in(:,:,:,ibox+1)-rho_in(:,:,:,ibox)
          del_fb(:,:,:,i)=fb(:,:,:,i+1)-fb(:,:,:,i)
        end do
      else if(ilsda==1)then
        do i=iter_s,iter_e
          ibox=num_rho_stock-(iter-i)
          del_x(:,:,:,i)=rho_s_in(:,:,:,is,ibox+1)-rho_s_in(:,:,:,is,ibox)
          del_fb(:,:,:,i)=fb(:,:,:,i+1)-fb(:,:,:,i)
        end do
      end if
      do i=iter_s,iter_e
        rbox=0.d0
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox=rbox+del_fb(ix,iy,iz,i)**2
        end do
        end do
        end do
    
        call comm_summation(rbox,rbox2,nproc_group_global)
    
        s=sqrt(rbox2)
    
        del_x(:,:,:,i)=del_x(:,:,:,i)/s
        del_fb(:,:,:,i)=del_fb(:,:,:,i)/s
      end do
      do i=1,iter_e-iter_s+1
        do j=1,iter_e-iter_s+1
          rbox=0.d0
          do iz=ng_sta(3),ng_end(3)
          do iy=ng_sta(2),ng_end(2)
          do ix=ng_sta(1),ng_end(1)
            rbox=rbox+del_fb(ix,iy,iz,iter_s-1+i)*del_fb(ix,iy,iz,iter_s-1+j)
          end do
          end do
          end do
          call comm_summation(rbox,rbox2,nproc_group_global)
          aa(i,j)=romega_mb(iter_s-1+i)*romega_mb(iter_s-1+j)*rbox2
          if (i == j) then
            aa(i,j)=aa(i,j)+romega0**2
          end if
        end do
      end do
      call matrix_inverse_ns(aa,iter_e-iter_s+1)
      beta(iter_s:iter_e,iter_s:iter_e)=aa(1:iter_e-iter_s+1,1:iter_e-iter_s+1)
      rho_temp(:,:,:)=0.d0
      do i=iter_s,iter_e
        do j=iter_s,iter_e
          rbox=0.d0
          do iz=ng_sta(3),ng_end(3)
          do iy=ng_sta(2),ng_end(2)
          do ix=ng_sta(1),ng_end(1)
            rbox=rbox+del_fb(ix,iy,iz,i)*fb(ix,iy,iz,iter)
          end do
          end do
          end do
          call comm_summation(rbox,rbox2,nproc_group_orbital)
          rho_temp(:,:,:)=rho_temp(:,:,:)&
              &+romega_mb(i)*romega_mb(j)*beta(i,j)*rbox2*(alpha_mb*del_fb(:,:,:,j)+del_x(:,:,:,j))
        end do
      end do
  
    
      if(ilsda==0)then
        call add_in_broyden(rho(:,:,:),rho_in(:,:,:,num_rho_stock),fb(:,:,:,iter),rho_temp(:,:,:))
      else if(ilsda==1)then
        call add_in_broyden(rho_s(:,:,:,is),rho_s_in(:,:,:,is,num_rho_stock),fb(:,:,:,iter),rho_temp(:,:,:))
      end if
    
    end do
    
    deallocate(fb,del_fb,del_x,romega_mb,beta,aa)
  end if
  
  if(ilsda==1)then
    rho(:,:,:)=rho_s(:,:,:,1)+rho_s(:,:,:,2)
  end if
  
  return
end subroutine broyden
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120----------
!Reference
!D.D. Johnson, Phys. Rev. B 38 12807 (1988)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120----------
subroutine matrix_inverse_ns(aa,nn)
  implicit none
  integer :: nn
  real(8) :: aa(nn,nn)
  integer :: ii,jj,kk,ll
  real(8) :: bb(nn,nn),work(nn),xx(nn,nn)
  
  bb(1:nn,1:nn)=0.d0
  do ii=1,nn
    bb(ii,ii)=1.d0
  end do
  
  !Transforming to Upper triangular matrix
  do kk=1,nn-1
    do ii=kk+1,nn
      do jj=kk+1,nn
        aa(ii,jj)=aa(ii,jj)/aa(ii,kk)-aa(kk,jj)/aa(kk,kk)
      end do
      bb(ii,:)=bb(ii,:)/aa(ii,kk)-bb(kk,:)/aa(kk,kk)
    end do
  end do
  !Back substitution
  do kk=nn,1,-1
    work(1:nn)=0.d0
    do ll=kk+1,nn
      work(1:nn)=work(1:nn)+aa(kk,ll)*xx(ll,1:nn)
    end do
    xx(kk,1:nn)=(bb(kk,1:nn)-work(1:nn))/aa(kk,kk)
  end do
  
  aa(1:nn,1:nn)=xx(1:nn,1:nn)
  
  return
end Subroutine matrix_inverse_ns
