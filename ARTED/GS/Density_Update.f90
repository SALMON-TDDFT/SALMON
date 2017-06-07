!
!  Copyright 2016 ARTED developers
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
!This file is "Density_Update.f"
!This file contain two subroutines.
!Subroutine Density_Update(iter)
!Subroutine Matrix_Inverse(a,n)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine Density_Update(iter)
  use Global_Variables
  implicit none
  integer,parameter :: iter_MB=0
  real(8),parameter :: omega0=0.01d0
  integer :: iter_s,iter_e
  integer :: iter,i,j
  real(8) :: F(1:NL,1:iter)
  real(8),allocatable :: del_F(:,:),del_X(:,:)
  real(8),allocatable :: omega_MB(:)
  real(8) :: s,rho_temp(1:NL)
  real(8),allocatable :: A(:,:),beta(:,:)
  integer :: Nnegative
      
  rho_out(1:NL,iter)=rho(1:NL)
  if (iter <= iter_MB+1) then
    F(1:NL,iter)=rho_out(1:NL,iter)-rho_in(1:NL,iter)
    rho_in(1:NL,iter+1)=rho_in(1:NL,iter)+alpha_MB*F(1:NL,iter)
  else
    iter_s=max(iter_MB+1,iter-Nmemory_MB)
    iter_e=iter-1
    allocate(del_F(1:NL,iter_s:iter_e))
    allocate(del_X(1:NL,iter_s:iter_e))
    allocate(omega_MB(iter_s:iter_e))
    allocate(beta(iter_s:iter_e,iter_s:iter_e))
    allocate(A(1:iter_e-iter_s+1,1:iter_e-iter_s+1))

    F(1:NL,iter_s:iter)=rho_out(1:NL,iter_s:iter)-rho_in(1:NL,iter_s:iter)
    omega_MB(iter_s:iter_e)=1.d0
!     &        /sqrt(sum((rho_out(:,iter_MB+1:iter-1)-rho_in(:,1:iter_MB+1:iter-1))**2))
    do i=iter_s,iter_e
      del_X(1:NL,i)=rho_in(1:NL,i+1)-rho_in(1:NL,i)
      del_F(1:NL,i)=F(1:NL,i+1)-F(1:NL,i)
    end do
    do i=iter_s,iter_e
      s=sqrt(sum(del_F(:,i)**2))
      del_X(1:NL,i)=del_X(1:NL,i)/s
      del_F(1:NL,i)=del_F(1:NL,i)/s
    end do
    do i=1,iter_e-iter_s+1
      do j=1,iter_e-iter_s+1
        A(i,j)=omega_MB(iter_s-1+i)*omega_MB(iter_s-1+j)*sum(del_F(:,iter_s-1+i)*del_F(:,iter_s-1+j))
        if (i == j) then
          A(i,j)=A(i,j)+omega0**2
        end if
      end do
    end do
    call Matrix_Inverse(a,iter_e-iter_s+1)
    beta(iter_s:iter_e,iter_s:iter_e)=A(1:iter_e-iter_s+1,1:iter_e-iter_s+1)
    rho_temp(1:NL)=0.d0
    do i=iter_s,iter_e
      do j=iter_s,iter_e
        rho_temp(1:NL)=rho_temp(1:NL)&
            &+omega_MB(i)*omega_MB(j)*beta(i,j)*sum(del_F(:,i)*F(:,iter))*(alpha_MB*del_F(1:NL,j)+del_X(1:NL,j))
      end do
    end do

    rho_in(1:NL,iter+1)=rho_in(1:NL,iter)+alpha_MB*F(1:NL,iter)-rho_temp(1:NL)
    Nnegative=0
    do i=1,NL
      if(rho_in(i,iter+1) < 0.d0) then
        Nnegative=Nnegative+1
      end if
    end do
    if(Nnegative > 0) then
      rho_in(1:NL,iter+1)=rho_in(1:NL,iter)+alpha_MB*F(1:NL,iter)
    end if

    deallocate(del_F,del_X,omega_MB,beta,A)
  end if
  rho(1:NL)=rho_in(1:NL,iter+1)

  return
End Subroutine Density_Update
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!Reference
!D.D. Johnson, Phys. Rev. B 38 12807 (1988)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine Matrix_Inverse(a,n)
  implicit none
  integer :: n
  real(8) :: a(n,n)
  integer :: i,j,k,l
  real(8) :: b(n,n),work(n),x(n,n)

  b(1:N,1:N)=0.d0
  do i=1,N
    b(i,i)=1.d0
  end do

!Transforming to Upper triangular matrix
  do k=1,n-1
    do i=k+1,n
      do j=k+1,n
        a(i,j)=a(i,j)/a(i,k)-a(k,j)/a(k,k)
      end do
      b(i,:)=b(i,:)/a(i,k)-b(k,:)/a(k,k)
    end do
  end do
!Back substitution
  do k=n,1,-1
    work(1:n)=0.d0
    do l=k+1,n
      work(1:n)=work(1:n)+a(k,l)*x(l,1:n)
    end do
    x(k,1:n)=(b(k,1:n)-work(1:n))/a(k,k)
  end do

  a(1:n,1:n)=x(1:n,1:n)

  return
End Subroutine Matrix_Inverse
