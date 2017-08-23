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
do nn=1,N_hamil
  if(mod(nn,2)==1) then
    if(ihpsieff==1) then
      call hpsi_test1(tzpsi_in,htpsi,Vlocal2)
    else
      call hpsi_test1(tzpsi_in,htpsi,Vlocal)
    end if
    call mode_add_polynomial(tzpsi_in,htpsi,tzpsi_out,nn)
  else
    if(ihpsieff==1) then
      call hpsi_test1(htpsi,tzpsi_in,Vlocal2)
    else
      call hpsi_test1(htpsi,tzpsi_in,Vlocal)
    end if
    call mode_add_polynomial(htpsi,tzpsi_in,tzpsi_out,nn)
  end if
end do
!-------------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------------
subroutine hpsi_test1(tpsi,htpsi,V)
  use hpsi_sub
  use init_sendrecv_sub
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                     mg_sta(2)-Nd:mg_end(2)+Nd,    &
                     mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum)
  complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                      mg_sta(2)-Nd:mg_end(2)+Nd,    &
                      mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum)
  real(8) :: V(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin)
  !
  integer :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
            ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
            ,ispin,i_all,Norb,i,iobmax,Nspin,Nk,ind,j,irank_overlap(6)
  real(8) :: lap0,lapt(4,3)
  integer, allocatable :: is_table(:),idx(:),idy(:),idz(:)

  irank_overlap(1) = iup_array(1)
  irank_overlap(2) = idw_array(1)
  irank_overlap(3) = jup_array(1)
  irank_overlap(4) = jdw_array(1)
  irank_overlap(5) = kup_array(1)
  irank_overlap(6) = kdw_array(1)

  lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      lapt(ind,j) = cNmat(ind,Nd)/Hgs(j)**2
    end do
  end do

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

  ix_sta = mg_sta(1)
  ix_end = mg_end(1)
  iy_sta = mg_sta(2)
  iy_end = mg_end(2)
  iz_sta = mg_sta(3)
  iz_end = mg_end(3)

  allocate(is_table(Norb))
  do i=1,iobmax
    call calc_allob(i,i_all)
    call set_ispin(i_all,ispin)
    is_table(i) = ispin
  end do

  allocate(idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4))
  do j=ix_sta-4,ix_end+4
    idx(j) = j
  end do
  do j=iy_sta-4,iy_end+4
    idy(j) = j
  end do
  do j=iz_sta-4,iz_end+4
    idz(j) = j
  end do

  call hpsi_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                 ,V,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                 ,idx,idy,idz,lap0,lapt,is_table,Nk,nproc_Mxin_mul,irank_overlap)

  deallocate(is_table,idx,idy,idz)
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

