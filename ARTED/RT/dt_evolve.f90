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
!This file is "dt_evolve.f90"
!This file contain one subroutine.
!Subroutine dt_evolve
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

subroutine dt_evolve_KB(iter)
  use global_variables, only: propagator,kAc,kAc0,Ac_tot
  implicit none
  integer, intent(in) :: iter

  select case(propagator)
    case('default')
      call default_propagator
    case('etrs')
      call etrs_propagator
    case default
      call err_finalize('invalid propagator')
  end select

contains
  subroutine default_propagator
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+0.5d0*(Ac_tot(iter,ixyz) + Ac_tot(iter+1,ixyz) )
    enddo
!$acc update device(kAc)
    call dt_evolve_omp_KB
  end subroutine

  subroutine etrs_propagator
    use global_variables, only: kAc_new
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter,ixyz)
      kAc_new(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter+1,ixyz)
    enddo
!$acc update device(kAc,kAc_new)
    call dt_evolve_etrs_omp_KB
  end subroutine
end subroutine

subroutine dt_evolve_KB_MS(ix_m,iy_m)
  use global_variables, only: propagator,kAc,kAc0,Ac_new_m,Ac_m
  implicit none
  integer, intent(in) :: ix_m, iy_m

  select case(propagator)
    case('default')
      call default_propagator
    case('etrs')
      call etrs_propagator
    case default
      call err_finalize('invalid propagator')
  end select

contains
  subroutine default_propagator
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+(Ac_new_m(ixyz,ix_m,iy_m)+Ac_m(ixyz,ix_m,iy_m))/2d0
    enddo
!$acc update device(kAc)
    call dt_evolve_omp_KB_MS
  end subroutine

  subroutine etrs_propagator
    use global_variables, only: kAc_new
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_m(ixyz,ix_m,iy_m)
      kAc_new(:,ixyz)=kAc0(:,ixyz)+Ac_new_m(ixyz,ix_m,iy_m)
    enddo
!$acc update device(kAc,kAc_new)
    call dt_evolve_etrs_omp_KB
  end subroutine
end subroutine

! ---------------------------------------------

Subroutine dt_evolve_omp_KB
  use Global_Variables
  use timer
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  use opt_variables
  implicit none
  integer    :: ik,ib
  integer    :: ia,j,i,ix,iy,iz
  real(8)    :: kr
  integer    :: thr_id,omp_get_thread_num,ikb

  NVTX_BEG('dt_evolve_omp_KB()',1)
  call timer_begin(LOG_DT_EVOLVE)

!$acc data pcopy(zu, vloc) pcopyout(ekr_omp)

!Constructing nonlocal part
  NVTX_BEG('dt_evolve_omp_KB(): nonlocal part',2)
#ifdef _OPENACC
!$acc kernels pcopy(ekr_omp) pcopyin(Mps, Jxyz,Jxx,Jyy,Jzz, kAc, Lx,Ly,Lz)
!$acc loop collapse(2) independent gang
#else
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
#endif
  do ik=NK_s,NK_e
  do ia=1,NI
!$acc loop independent vector(128)
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
#ifdef _OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif
  NVTX_END()

! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ')
!$acc update self(zu, ekr_omp, vloc)

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do
    
  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call hamiltonian(.false.)

  call psi_rho_RT
  call Hartree
  call Exc_Cor('RT')

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if


!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

!$acc update device(zu, vloc)
  end select
! yabana

  NVTX_BEG('dt_evolve_omp_KB(): hamiltonian',3)
  call hamiltonian(.true.)
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB(): psi_rho_RT',4)
  call psi_rho_RT
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB(): Hartree',5)
  call Hartree
  NVTX_END()

! yabana
  NVTX_BEG('dt_evolve_omp_KB(): Exc_Cor',6)
  call Exc_Cor('RT')
  NVTX_END()
! yabana


#ifdef _OPENACC
!$acc kernels pcopy(Vloc) pcopyin(Vh,Vpsl,Vexc)
#else
!$omp parallel do
#endif
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do
!$acc end kernels

!$acc end data

  call timer_end(LOG_DT_EVOLVE)
  NVTX_END()

  return
End Subroutine dt_evolve_omp_KB
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine dt_evolve_etrs_omp_KB
  use Global_Variables
  use timer
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  use opt_variables
  implicit none
  integer    :: ik,ib
  integer    :: ia,j,i,ix,iy,iz
  real(8)    :: kr,dt_t
  integer    :: thr_id,omp_get_thread_num,ikb

  NVTX_BEG('dt_evolve_omp_KB()',1)
  call timer_begin(LOG_DT_EVOLVE)

  dt_t = dt; dt = 0.5d0*dt

!$acc data pcopy(zu, vloc) pcopyout(ekr_omp)

!Constructing nonlocal part
  NVTX_BEG('dt_evolve_omp_KB(): nonlocal part',2)
#ifdef _OPENACC
!$acc kernels pcopy(ekr_omp) pcopyin(Mps, Jxyz,Jxx,Jyy,Jzz, kAc, Lx,Ly,Lz)
!$acc loop collapse(2) independent gang
#else
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
#endif
  do ik=NK_s,NK_e
  do ia=1,NI
!$acc loop independent vector(128)
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
#ifdef _OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif
  NVTX_END()


!$acc update self(zu, ekr_omp, vloc)

  NVTX_BEG('dt_evolve_omp_KB(): hamiltonian',3)
  call hamiltonian(.false.)
  NVTX_END()


  Vloc_t=Vloc
  Vloc_new(:) = 3d0*Vloc(:) - 3d0*Vloc_old(:,1) + Vloc_old(:,2)
  Vloc_old(:,2) = Vloc_old(:,1)
  Vloc_old(:,1) = Vloc(:)
  Vloc(:) = Vloc_new(:)

  kAc=kAc_new
!$acc update device(kAc)

!Constructing nonlocal part
  NVTX_BEG('dt_evolve_omp_KB(): nonlocal part',2)
#ifdef _OPENACC
!$acc kernels pcopy(ekr_omp) pcopyin(Mps, Jxyz,Jxx,Jyy,Jzz, kAc, Lx,Ly,Lz)
!$acc loop collapse(2) independent gang
#else
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
#endif
  do ik=NK_s,NK_e
  do ia=1,NI
!$acc loop independent vector(128)
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
#ifdef _OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif
  NVTX_END()

!== predictor-corrector ==
  select case(functional)
  case('VS98','TPSS','TBmBJ')
!$acc update self(zu, ekr_omp, vloc)

!$omp parallel do private(ik,ib)
     do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)
        zu_GS(:,ib,ik)=zu(:,ib,ik)
     end do

     NVTX_BEG('dt_evolve_omp_KB(): hamiltonian',3)
     call hamiltonian(.false.)
     NVTX_END()

     NVTX_BEG('dt_evolve_omp_KB(): psi_rho_RT',4)
     call psi_rho_RT
     NVTX_END()

     NVTX_BEG('dt_evolve_omp_KB(): Hartree',5)
     call Hartree
     NVTX_END()

     NVTX_BEG('dt_evolve_omp_KB(): Exc_Cor',6)
     call Exc_Cor('RT')
     NVTX_END()

#ifdef _OPENACC
!$acc kernels pcopy(Vloc) pcopyin(Vh,Vpsl,Vexc)
#else
!$omp parallel do
#endif
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do
!$acc end kernels


!$omp parallel do private(ik,ib)
     do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)
        zu(:,ib,ik)=zu_GS(:,ib,ik)
     end do

  end select


  NVTX_BEG('dt_evolve_omp_KB(): hamiltonian',3)
  call hamiltonian(.true.)
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB(): psi_rho_RT',4)
  call psi_rho_RT
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB(): Hartree',5)
  call Hartree
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB(): Exc_Cor',6)
  call Exc_Cor('RT')
  NVTX_END()


#ifdef _OPENACC
!$acc kernels pcopy(Vloc) pcopyin(Vh,Vpsl,Vexc)
#else
!$omp parallel do
#endif
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do
!$acc end kernels

!$acc end data

  dt = dt_t
  call timer_end(LOG_DT_EVOLVE)
  NVTX_END()

  return
End Subroutine dt_evolve_etrs_omp_KB
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine dt_evolve_omp_KB_MS
  use Global_Variables
  use timer
  use nvtx
  use opt_variables
  implicit none
  integer    :: ik,ib
  integer    :: ia,j,i,ix,iy,iz
  real(8)    :: kr
  integer    :: thr_id,omp_get_thread_num,ikb

  NVTX_BEG('dt_evolve_omp_KB_MS()',1)
  call timer_begin(LOG_DT_EVOLVE)

!$acc data pcopy(zu, vloc) pcopyout(ekr_omp)

!Constructing nonlocal part ! sato
  NVTX_BEG('dt_evolve_omp_KB_MS(): nonlocal part',2)
#ifdef _OPENACC
!$acc kernels pcopy(ekr_omp) pcopyin(Mps, Jxyz,Jxx,Jyy,Jzz, kAc, Lx,Ly,Lz)
!$acc loop collapse(2) independent gang
#else
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
#endif
  do ik=NK_s,NK_e
  do ia=1,NI
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
#ifdef _OPENACC
!$acc end kernels
#else
!$omp end parallel
#endif

! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ')
!$acc update self(zu, ekr_omp, vloc)


!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do

  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call hamiltonian(.false.)

  call psi_rho_RT
  call Hartree
  call Exc_Cor('RT')

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

  end select
! yabana

  NVTX_BEG('dt_evolve_omp_KB_MS(): hamiltonian',3)
  call hamiltonian(.true.)
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB_MS(): psi_rho_RT',4)
  call psi_rho_RT
  NVTX_END()

  NVTX_BEG('dt_evolve_omp_KB_MS(): Hartree',5)
  call Hartree
  NVTX_END()

! yabana
  NVTX_BEG('dt_evolve_omp_KB_MS(): Hartree',5)
  call Exc_Cor('RT')
  NVTX_END()
! yabana

#ifdef _OPENACC
!$acc kernels pcopy(Vloc) pcopyin(Vh,Vpsl,Vexc)
#else
!$omp parallel do
#endif
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do
!$acc end kernels

!$acc end data

  call timer_end(LOG_DT_EVOLVE)
  NVTX_END()

  return
End Subroutine dt_evolve_omp_KB_MS
