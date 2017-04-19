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
!This file is "psi_rho.f90"
!This file contain one subroutine.
!Subroutine psi_rho
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef __INTEL_COMPILER
# define OMP_SIMD simd
#else
# define OMP_SIMD
#endif

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

subroutine psi_rho_GS
  use global_variables, only: zu_GS,NB
  implicit none
  call psi_rho_impl(zu_GS,NB)
end subroutine

subroutine psi_rho_RT(zu)
  use global_variables, only: NL,NBoccmax,NK_s,NK_e
  implicit none
  complex(8),intent(in) :: zu(NL,NBoccmax,NK_s:NK_e)
  call psi_rho_impl(zu,NBoccmax)
end subroutine

subroutine psi_rho_impl(zutmp,zu_NB)
  use global_variables
  use timer
  use opt_variables
  use communication
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  implicit none
  integer,intent(in)    :: zu_NB
  complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

  call timer_begin(LOG_PSI_RHO)
  ! write(*,*) "Sym:", Sym
  ! stop
  select case(Sym)
  case(1)
    call sym1(zutmp,zu_NB,rho_l)
  case(4)
    if(crystal_structure == 'diamond')then
      call sym4(zutmp,zu_NB,rho_l,rho_tmp1)
    else
      call err_finalize('Bad crystal structure')
    end if
  case(8)
    if(crystal_structure == 'diamond')then
       call sym8(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else if(crystal_structure == 'diamond2')then
       call sym8_diamond2(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else if(crystal_structure == 'tetragonal')then
       call sym8_tetragonal(zutmp,zu_NB,rho_l,rho_tmp1,rho_tmp2)
    else
       call err_finalize('Bad crystal structure')
    end if
  case default
    call err_finalize('Bad Symmetry')
  end select
  call timer_end(LOG_PSI_RHO)

  call timer_begin(LOG_ALLREDUCE)
  call comm_summation(rho_l,rho,NL,proc_group(2))
  call timer_end(LOG_ALLREDUCE)


contains
  subroutine reduce(tid,zfac,zutmp,zu_NB)
    use global_variables
    use opt_variables, only: zrhotmp
    use misc_routines, only: ceiling_pow2
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: tid
    real(8),intent(in)    :: zfac
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

    integer :: ib,ik,i,mytid

    mytid = tid

#ifdef ARTED_REDUCE_FOR_MANYCORE
    zrhotmp(:,mytid)=0.d0

!$omp do private(ik,ib,i) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
    do i=0,NL-1
      zrhotmp(i,mytid)=zrhotmp(i,mytid)+(zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
    end do
    end do
    end do
!$omp end do

    i = ceiling_pow2(NUMBER_THREADS)/2
    do while(i > 0)
      if(mytid < i) then
        zrhotmp(0:NL-1,mytid) = zrhotmp(0:NL-1,mytid) + zrhotmp(0:NL-1,mytid + i)
      end if
      i = i/2
!$omp barrier
    end do
#else
    mytid = 0

!$omp single
    zrhotmp(:,mytid) = 0.d0
!$omp end single

    do ik=NK_s,NK_e
    do ib=1,NBoccmax
!$omp do private(i)
    do i=0,NL-1
      zrhotmp(i,mytid)=zrhotmp(i,mytid)+(zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
    end do
!$omp end do
    end do
    end do
#endif
  end subroutine

  subroutine reduce_acc(zfac, zutmp, zu_NB, zrho)
    use global_variables
    implicit none
    real(8),intent(in)    :: zfac
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1, zu_NB, NK_s:NK_e)
    real(8),intent(out)   :: zrho(0:NL-1)

    integer :: ib,ik,i
    real(8) :: my_zrho

    NVTX_BEG("reduce_acc()",1)
!$acc kernels pcopy(zrho) pcopyin(zutmp,occ)
    zrho(:)=0.d0

!$acc loop gang vector(1)
    do ik=NK_s,NK_e
!$acc loop gang vector(128) private(my_zrho)
      do i=0,NL-1
        my_zrho = 0.d0
!$acc loop seq
        do ib=1,NBoccmax
          my_zrho = my_zrho + (zfac*occ(ib,ik))*abs(zutmp(i,ib,ik))**2
        end do
!$acc atomic update
        zrho(i) = zrho(i) + my_zrho
!$acc end atomic
      end do
    end do
!$acc end kernels
    NVTX_END()
  end subroutine

  subroutine sym1(zutmp,zu_NB,zrho_l)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    integer :: tid

    NVTX_BEG("sym1()",1)

#ifdef _OPENACC
    call reduce_acc(1.0d0,zutmp,zu_NB,zrho_l)
#else
!$omp parallel private(tid)
!$  tid=omp_get_thread_num()
    call reduce(tid,1.0d0,zutmp,zu_NB)
!$omp end parallel

    zrho_l(:) = zrhotmp(0:NL-1,0)
#endif

    NVTX_END()
  end subroutine

  !====== diamond(4) structure =========================!
  subroutine sym4(zutmp,zu_NB,zrho_l,zrhotmp1)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    NVTX_BEG("sym4()",2)

    zfac=1.0d0/4d0

#ifdef _OPENACC
    call reduce_acc(zfac,zutmp,zu_NB,zrhotmp(:,0))
#endif

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

#ifndef _OPENACC
    call reduce(tid,zfac,zutmp,zu_NB)
#endif

! 1.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp(i,0)+zrhotmp(itable_sym(1,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel

    NVTX_END()
  end subroutine

  !====== diamond(8) structure =========================!
  subroutine sym8(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    NVTX_BEG("sym8()",3)

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/32d0

#ifdef _OPENACC
    call reduce_acc(zfac,zutmp,zu_NB,zrhotmp(:,0))
#endif

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

#ifndef _OPENACC
    call reduce(tid,zfac,zutmp,zu_NB)
#endif

! 1.T_4
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp(i,0)+zrhotmp(itable_sym(4,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(5,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(3,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_1
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 2.T_2
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

   !====== diamond2(8) structure =========================!
  subroutine sym8_diamond2(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    NVTX_BEG("sym8_diamond2()",4)

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/16d0

#ifdef _OPENACC
    call reduce_acc(zfac,zutmp,zu_NB,zrhotmp(:,0))
#endif

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

#ifndef _OPENACC
    call reduce(tid,zfac,zutmp,zu_NB)
#endif

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp(i,0)+zrhotmp(itable_sym(4,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(3,i+1)-1)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp1(i)+zrhotmp1(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel
  end subroutine

  !====== tetragonal(8) structure =========================!
  subroutine sym8_tetragonal(zutmp,zu_NB,zrho_l,zrhotmp1,zrhotmp2)
    use global_variables
    use opt_variables, only: zrhotmp
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(in)    :: zu_NB
    complex(8),intent(in) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)
    real(8),intent(out)   :: zrho_l(0:NL-1)

    real(8) :: zrhotmp1(0:NL-1)
    real(8) :: zrhotmp2(0:NL-1)
    integer :: i,tid
    real(8) :: zfac

    ! wk(ik)=8.0,(ikx==iky >. wk(ik)=4.0)
    zfac=1.0d0/16d0

#ifdef _OPENACC
    call reduce_acc(zfac,zutmp,zu_NB,zrhotmp(:,0))
#endif

!$omp parallel private(tid)
!$  tid=omp_get_thread_num()

#ifndef _OPENACC
    call reduce(tid,zfac,zutmp,zu_NB)
#endif

! 1.T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp(i,0)+zrhotmp(itable_sym(3,i+1)-1,0)
    end do
!$omp end do OMP_SIMD

! 2.T_3*T_3
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp1(i)=zrhotmp2(i)+zrhotmp2(itable_sym(4,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 3.T_1
!$omp do OMP_SIMD
    do i=0,NL-1
      zrhotmp2(i)=zrhotmp1(i)+zrhotmp1(itable_sym(1,i+1)-1)
    end do
!$omp end do OMP_SIMD

! 4.T_2
!$omp do OMP_SIMD
    do i=0,NL-1
      zrho_l(i)=zrhotmp2(i)+zrhotmp2(itable_sym(2,i+1)-1)
    end do
!$omp end do OMP_SIMD
!$omp end parallel

    NVTX_END()
  end subroutine
end subroutine psi_rho_impl

