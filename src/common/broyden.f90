!
!  Copyright 2017 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine broyden(vecr,vecr_in,vecr_out,nl,iter,iter_mod,nstock)
  use inputoutput, only: iperiodic, alpha_mb, nmemory_mb
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_summation
  use salmon_math
  implicit none
  integer, intent(in) :: nl
  integer, intent(in) :: iter
  integer, intent(in) :: iter_mod
  integer, intent(in) :: nstock
  real(8), intent(inout) :: vecr(1:nl)
  real(8), intent(inout) :: vecr_in(1:nl,1:nstock+1)
  real(8), intent(inout) :: vecr_out(1:nl,1:nstock+1)
  integer,parameter :: iter_mb=0
  real(8),parameter :: omega0=0.01d0
  integer :: iter_s,iter_e
  integer :: i,j
  real(8),allocatable :: vecf(:,:)
  real(8) :: vecr_tmp(1:nl)
  real(8),allocatable :: del_vecf(:,:),del_vecx(:,:)
  real(8),allocatable :: omega_mb(:)
  real(8),allocatable :: aa(:,:),aa_tmp1(:,:)
  real(8),allocatable :: beta(:,:)
  real(8),allocatable :: ss_tmp1(:), ss_tmp2(:)
  real(8),allocatable :: ss_tmp3(:,:), ss_tmp4(:,:)
  real(8) :: ss
  integer :: nnegative,nnegative_tmp

  vecr_out(1:nl,iter_mod)=vecr(1:nl)
  if (iter <= iter_mb+1) then
    allocate(vecf(1:nl,iter:iter))
    vecf(1:nl,iter)=vecr_out(1:nl,iter_mod)-vecr_in(1:nl,iter_mod)
    vecr_in(1:nl,iter_mod+1)=vecr_in(1:nl,iter_mod)+alpha_mb*vecf(1:nl,iter)
    deallocate(vecf)
  else
    iter_s=max(iter_mb+1+(iter_mod-iter),iter_mod-nmemory_mb)
    iter_e=iter_mod-1
    allocate(vecf(1:nl,iter_s:iter_e+1))
    allocate(del_vecf(1:nl,iter_s:iter_e))
    allocate(del_vecx(1:nl,iter_s:iter_e))
    allocate(omega_mb(iter_s:iter_e))
    allocate(beta(iter_s:iter_e,iter_s:iter_e))
    allocate(aa(1:iter_e-iter_s+1,1:iter_e-iter_s+1))
    allocate(aa_tmp1(1:iter_e-iter_s+1,1:iter_e-iter_s+1))
    if(iperiodic==0)then
      allocate(ss_tmp1(iter_s:iter_e))
      allocate(ss_tmp2(iter_s:iter_e))
    end if

    vecf(1:nl,iter_s:iter_mod)=vecr_out(1:nl,iter_s:iter_mod)  &
                             -vecr_in(1:nl,iter_s:iter_mod)
    omega_mb(iter_s:iter_e)=1.d0
    do i=iter_s,iter_e
      del_vecx(1:nl,i)=vecr_in(1:nl,i+1)-vecr_in(1:nl,i)
      del_vecf(1:nl,i)=vecf(1:nl,i+1)-vecf(1:nl,i)
    end do

    if(iperiodic==0)then
      do i=iter_s,iter_e
        ss_tmp1(i)=sum(del_vecf(1:nl,i)**2)
      end do
      call comm_summation(ss_tmp1,ss_tmp2,iter_e-iter_s+1,nproc_group_global)
      do i=iter_s,iter_e
        ss=sqrt(ss_tmp2(i))
        del_vecx(1:nl,i)=del_vecx(1:nl,i)/ss
        del_vecf(1:nl,i)=del_vecf(1:nl,i)/ss
      end do
    else if(iperiodic==3)then
      do i=iter_s,iter_e
        ss=sum(del_vecf(1:nl,i)**2)
        del_vecx(1:nl,i)=del_vecx(1:nl,i)/ss
        del_vecf(1:nl,i)=del_vecf(1:nl,i)/ss
      end do
    end if

    if(iperiodic==0)then
      do i=1,iter_e-iter_s+1
        do j=1,iter_e-iter_s+1
          aa_tmp1(i,j)=sum(del_vecf(1:nl,iter_s-1+i)*del_vecf(1:nl,iter_s-1+j))
        end do
      end do
      call comm_summation(aa_tmp1,aa,(iter_e-iter_s+1)**2,nproc_group_global)
      do i=1,iter_e-iter_s+1
        do j=1,iter_e-iter_s+1
          aa(i,j)=omega_mb(iter_s-1+i)*omega_mb(iter_s-1+j)*aa(i,j)
          if(i==j)then
            aa(i,j)=aa(i,j)+omega0**2
          end if
        end do
      end do
    else if(iperiodic==3)then
      do i=1,iter_e-iter_s+1
        do j=1,iter_e-iter_s+1
          aa(i,j)=omega_mb(iter_s-1+i)*omega_mb(iter_s-1+j)*sum(del_vecf(1:nl,iter_s-1+i)*del_vecf(1:nl,iter_s-1+j))
          if(i==j)then
            aa(i,j)=aa(i,j)+omega0**2
          end if
        end do
      end do
    end if

    call matrix_inverse(aa)

    beta(iter_s:iter_e,iter_s:iter_e)=aa(1:iter_e-iter_s+1,1:iter_e-iter_s+1)
    if(iperiodic==0)then
      do i=iter_s,iter_e
        ss_tmp1(i)=sum(del_vecf(:,i)*vecf(:,iter_mod))
      end do
      call comm_summation(ss_tmp1,ss_tmp2,iter_e-iter_s+1,nproc_group_global)
      vecr_tmp(1:nl)=0.d0
      do i=iter_s,iter_e
        do j=iter_s,iter_e
          vecr_tmp(1:nl)=vecr_tmp(1:nl)&
              &+omega_mb(i)*omega_mb(j)*beta(i,j)*ss_tmp2(i)*(alpha_mb*del_vecf(1:nl,j)+del_vecx(1:nl,j))
        end do
      end do
    else if(iperiodic==3)then
      vecr_tmp(1:nl)=0.d0
      do i=iter_s,iter_e
        do j=iter_s,iter_e
          vecr_tmp(1:nl)=vecr_tmp(1:nl)&
              &+omega_mb(i)*omega_mb(j)*beta(i,j)*sum(del_vecf(:,i)*vecf(:,iter_mod))*(alpha_mb*del_vecf(1:nl,j)+del_vecx(1:nl,j))
        end do
      end do
    end if
  
    vecr_in(1:nl,iter_mod+1)=vecr_in(1:nl,iter_mod)+alpha_mb*vecf(1:nl,iter_mod)-vecr_tmp(1:nl)

    if(iperiodic==0)then
      nnegative_tmp=0
      do i=1,nl
        if(vecr_in(i,iter_mod+1) < 0.d0) then
          nnegative_tmp=nnegative_tmp+1
        end if
      end do
      call comm_summation(nnegative_tmp,nnegative,nproc_group_global)
    else
      nnegative=0
      do i=1,nl
        if(vecr_in(i,iter_mod+1) < 0.d0) then
          nnegative=nnegative+1
        end if
      end do
    end if
    if(nnegative > 0) then
      vecr_in(1:nl,iter_mod+1)=vecr_in(1:nl,iter_mod)+alpha_mb*vecf(1:nl,iter_mod)
    end if

    deallocate(vecf,del_vecf,del_vecx,omega_mb,beta,aa,aa_tmp1)
    if(iperiodic==0)then
      deallocate(ss_tmp1,ss_tmp2)
    end if
  end if
  vecr(1:nl)=vecr_in(1:nl,iter_mod+1)

  return
end subroutine broyden
