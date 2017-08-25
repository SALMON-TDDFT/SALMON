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
MODULE hpsi2_sub

use scf_data
use laplacian2_sub
use new_world_sub
use gradient2_sub
use allocate_mat_sub
implicit none

INTERFACE hpsi2
!   MODULE PROCEDURE R_hpsi2,C_hpsi2
  MODULE PROCEDURE hpsi_test2_R,hpsi_test2_C
END INTERFACE

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine hpsi_test2_R(tpsi0,htpsi0,iob,nn,isub)
  use salmon_parallel, only: nproc_group_orbital
  use init_sendrecv_sub
  use hpsi_sub
  implicit none
  real(8) :: tpsi0(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
  real(8) :: htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
  integer :: iob,nn,isub

  integer :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
            ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
            ,ispin,Norb,Nspin,is_table(1),ind,j,irank_overlap(6),icomm_overlap,icomm_pseudo &
            ,NI,Nps,Nlma
  real(8) :: lap0,lapt(4,3)
  integer,allocatable :: idx(:),idy(:),idz(:)
  real(8),allocatable :: htpsi(:,:,:,:)

  integer,allocatable :: ia_table(:),Mps(:),Jxyz_wrk(:,:,:)
  real(8),allocatable :: uV_wrk(:,:),uVu_wrk(:)

  irank_overlap(1) = iup_array(1)
  irank_overlap(2) = idw_array(1)
  irank_overlap(3) = jup_array(1)
  irank_overlap(4) = jdw_array(1)
  irank_overlap(5) = kup_array(1)
  irank_overlap(6) = kdw_array(1)

  lap0 = -0.5d0*cNmat(0,4)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      lapt(ind,j) = cNmat(ind,4)/Hgs(j)**2
    end do
  end do

  Norb = 1
  Nspin = numspin

  call set_ispin(iob,ispin)
  is_table(1) = ispin

  ipx_sta = iwksta(1)
  ipx_end = iwkend(1)
  ipy_sta = iwksta(2)
  ipy_end = iwkend(2)
  ipz_sta = iwksta(3)
  ipz_end = iwkend(3)

  ix_sta = mg_sta(1)
  ix_end = mg_end(1)
  iy_sta = mg_sta(2)
  iy_end = mg_end(2)
  iz_sta = mg_sta(3)
  iz_end = mg_end(3)

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

  icomm_overlap = nproc_group_orbital
  call convert_pseudo_GCEED(NI,Nps,Nlma,ia_table,Mps,uV_wrk,uVu_wrk,Jxyz_wrk,icomm_pseudo)

  allocate(htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb))
  htpsi = 0d0

  call hpsi_R(tpsi0,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                   ,Vlocal,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                   ,idx,idy,idz,lap0,lapt,is_table &
                   ,NI,Nps,Nlma,ia_table,Mps,Jxyz_wrk,uV_wrk,uVu_wrk &
                   ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo)

  htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3)) &
    = htpsi(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3),1)

  deallocate(idx,idy,idz,htpsi,uV_wrk,Jxyz_wrk,Mps,ia_table,uVu_wrk)
end subroutine hpsi_test2_R

subroutine hpsi_test2_C(tpsi0,htpsi0,iob,nn,isub)
  use salmon_parallel, only: nproc_group_orbital
  use init_sendrecv_sub
  use hpsi_sub
  implicit none
  complex(8) :: tpsi0(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
  complex(8) :: htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
  integer :: iob,nn,isub

  integer :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
            ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
            ,ispin,Norb,Nspin,Nk,is_table(1),ind,j,irank_overlap(6),icomm_overlap,icomm_pseudo &
            ,NI,Nps,Nlma
  real(8) :: lap0,lapt(4,3)
  integer,allocatable :: idx(:),idy(:),idz(:)
  complex(8),allocatable :: htpsi(:,:,:,:)

  integer,allocatable :: ia_table(:),Mps(:),Jxyz_wrk(:,:,:)
  real(8),allocatable :: uV_wrk(:,:),uVu_wrk(:)

  irank_overlap(1) = iup_array(1)
  irank_overlap(2) = idw_array(1)
  irank_overlap(3) = jup_array(1)
  irank_overlap(4) = jdw_array(1)
  irank_overlap(5) = kup_array(1)
  irank_overlap(6) = kdw_array(1)

  lap0 = -0.5d0*cNmat(0,4)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  do j=1,3
    do ind=1,4
      lapt(ind,j) = cNmat(ind,4)/Hgs(j)**2
    end do
  end do

  Norb = 1
  Nspin = numspin
  Nk = 1

  call set_ispin(iob,ispin)
  is_table(1) = ispin

  ipx_sta = iwksta(1)
  ipx_end = iwkend(1)
  ipy_sta = iwksta(2)
  ipy_end = iwkend(2)
  ipz_sta = iwksta(3)
  ipz_end = iwkend(3)

  ix_sta = mg_sta(1)
  ix_end = mg_end(1)
  iy_sta = mg_sta(2)
  iy_end = mg_end(2)
  iz_sta = mg_sta(3)
  iz_end = mg_end(3)

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

  icomm_overlap = nproc_group_orbital
  call convert_pseudo_GCEED(NI,Nps,Nlma,ia_table,Mps,uV_wrk,uVu_wrk,Jxyz_wrk,icomm_pseudo)

  allocate(htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb))
  htpsi = 0d0

  call hpsi_C(tpsi0,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                   ,Vlocal,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                   ,idx,idy,idz,lap0,lapt,is_table,Nk &
                   ,NI,Nps,Nlma,ia_table,Mps,Jxyz_wrk,uV_wrk,uVu_wrk &
                   ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo)

  htpsi0(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3)) &
    = htpsi(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3),1)

  deallocate(idx,idy,idz,htpsi,uV_wrk,uVu_wrk,Jxyz_wrk,Mps,ia_table)
end subroutine hpsi_test2_C

subroutine convert_pseudo_GCEED(NI,Nps,Nlma,ia_table,Mps,uV_new,uVu_new,Jxyz_new,icomm_pseudo)
  use scf_data, only: MI,Kion,Mlps,uVu,iwk_size,max_jMps_l,uV,Jxyz,jMps_l,max_jMps_l_s,jMps_l_s,uVu,Hvol ! GCEED
  use salmon_parallel, only: nproc_group_orbital, nproc_group_h
  integer :: NI,Nps,Nlma,icomm_pseudo
  integer,allocatable :: ia_table(:),Mps(:),Jxyz_new(:,:,:)
  real(8),allocatable :: uV_new(:,:),uVu_new(:)
  !
  integer :: jj,iatom,ik,lm,ilma,jm
  integer,allocatable :: jMps(:,:)

  NI = MI

  allocate(Mps(NI))

  if(iwk_size>=1.and.iwk_size<=2)then
    Mps = max_jMps_l
    icomm_pseudo = nproc_group_orbital
  else if(iwk_size>=11.and.iwk_size<=12)then
    Mps = max_jMps_l_s
    icomm_pseudo = nproc_group_h
  end if

  allocate(jMps(maxval(Mps),NI))

  if(iwk_size>=1.and.iwk_size<=2)then
    jMps = jMps_l
  else if(iwk_size>=11.and.iwk_size<=12)then
    jMps = jMps_l_s
  end if

  Nps = maxval(Mps)
  allocate(Jxyz_new(3,Nps,NI))

  do iatom=1,MI
    do jj=1,Mps(iatom)
      jm = jMps(jj,iatom)

      Jxyz_new(:,jj,iatom) = Jxyz(:,jm,iatom)
    end do
  end do

  ilma = 0
  do iatom=1,MI
    ik=Kion(iatom)
    loop_lm : do lm=1,(Mlps(ik)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm
      ilma = ilma + 1
    end do loop_lm
  end do
  Nlma = ilma

  allocate(ia_table(Nlma),uV_new(Nps,Nlma),uVu_new(Nlma))
  uV_new = 0d0

  ilma = 0
  do iatom=1,MI
    ik=Kion(iatom)
    loop_lm2 : do lm=1,(Mlps(ik)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
      ilma = ilma + 1

      ia_table(ilma) = iatom
      do jj=1,Mps(iatom)
        jm = jMps(jj,iatom)

        uV_new(jj,ilma) = uV(jm,lm,iatom)
      end do
      uVu_new(ilma) = Hvol/uVu(lm,iatom)

    end do loop_lm2
  end do

  deallocate(jMps)
  return
end subroutine convert_pseudo_GCEED
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE hpsi2_sub

