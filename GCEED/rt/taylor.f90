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
subroutine taylor(tzpsi_in,tzpsi_out,htpsi)
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
implicit none

integer :: nn,ix,iy,iz
complex(8) :: tzpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                    mg_sta(2)-Nd:mg_end(2)+Nd,    &
                    mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
complex(8) :: tzpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)

!-------------------------------------------------------------------------------------
real(8),allocatable :: V_wrk(:,:,:,:)
!-------------------------------------------------------------------------------------

iwk_size=2
call make_iwksta_iwkend

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox(ix,iy,iz) = 0.d0
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox_s(ix,iy,iz,1) = 0.d0
    rhobox_s(ix,iy,iz,2) = 0.d0
  end do
  end do
  end do
end if

if(ihpsieff==1)then
  if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then 
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
      Vlocal2(ix,iy,iz,2)=Vlocal(ix,iy,iz,2)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  end if
end if

!-------------------------------------------------------------------------------------
allocate(V_wrk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin))
if(ihpsieff==1)then
  V_wrk = Vlocal2
else
  V_wrk = Vlocal
end if
do nn=1,N_hamil
  if(mod(nn,2)==1)then
    call hpsi_test1(tzpsi_in,htpsi,V_wrk)
    call mode_add_polynomial(tzpsi_in,htpsi,tzpsi_out,nn)
  else
    call hpsi_test1(htpsi,tzpsi_in,V_wrk)
    call mode_add_polynomial(htpsi,tzpsi_in,tzpsi_out,nn)
  end if
end do
deallocate(V_wrk)
!-------------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------------
subroutine hpsi_test1(tpsi0,htpsi0,V)
  use hpsi_sub
  implicit none
  complex(8) :: tpsi0(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
  complex(8) :: htpsi0(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                      mg_sta(2)-Nd:mg_end(2)+Nd,    &
                      mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
  real(8) :: V(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin)
  !
  integer :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end
  integer :: ispin,i_all,Norb,i,iobmax,Nspin,Nk
  complex(8),allocatable :: tpsi(:,:,:,:,:,:),htpsi(:,:,:,:,:,:)

  call calc_pmax(iobmax)
  Norb = iobnum
  Nspin = numspin
  Nk = 1

  ipx_sta = mg_sta(1)-Nd
  ipx_end = mg_end(1)+Nd+1
  ipy_sta = mg_sta(2)-Nd
  ipy_end = mg_end(2)+Nd
  ipz_sta = mg_sta(3)-Nd
  ipz_end = mg_end(3)+Nd

  allocate( tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Nspin,1:Norb,1:Nk) &
          ,htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Nspin,1:Norb,1:Nk))
  tpsi = 0d0
  htpsi = 0d0

  do i=1,iobmax
    call calc_allob(i,i_all)
    call set_ispin(i_all,ispin)
    tpsi(:,:,:,ispin,i,1) = tpsi0(:,:,:,i,1)
  end do

  call hpsi(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Nspin,Norb,Nk,Nd &
            ,V,mg_sta(1),mg_end(1),mg_sta(2),mg_end(2),mg_sta(3),mg_end(3),nproc_Mxin_mul)

  do i=1,iobmax
    call calc_allob(i,i_all)
    call set_ispin(i_all,ispin)
    htpsi0(:,:,:,i,1) = htpsi(:,:,:,ispin,i,1)
  end do

  deallocate(tpsi,htpsi)
  return
end subroutine hpsi_test1
!-------------------------------------------------------------------------------------

! hpsi_groupob.f90 --> taylor.f90
!-------------------------------------------------------------------------------------
subroutine mode_add_polynomial(tpsi,htpsi,tpsi_out,nn)
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
  complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                      mg_sta(2)-Nd:mg_end(2)+Nd,    &
                      mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
  complex(8) :: tpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
  integer :: nn
  !
  integer :: iobmax
  call calc_pmax(iobmax)

  if(N_hamil==1)then
    if(ikind_eext==0)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,5)
    else
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,4)
    end if
  else if(N_hamil==4)then
    if(nn==1)then
      if(ikind_eext==0)then
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,1)
      else
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,0)
      end if
    else if(nn==3)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,-1)
    else if(nn==4)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,6)
    end if
  else
    if(nn==1)then
      if(ikind_eext==0)then
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,1)
      else
        call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,0)
      end if
    else if(nn==N_hamil)then
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,6)
    else
      call add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,2)
    end if
  end if

  return
end subroutine
!-------------------------------------------------------------------------------------

end subroutine taylor

