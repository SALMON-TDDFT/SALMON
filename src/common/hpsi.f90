!
!  Copyright 2017-2019 SALMON developers
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
module hpsi_sub

interface hpsi
  module procedure hpsi_R,hpsi_C
end interface hpsi

contains

!===================================================================================================================================

SUBROUTINE hpsi_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                 ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                 ,idx,idy,idz,lap0,lapt,is_table &
                 ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
                 ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo)
  use update_overlap_sub, only: update_overlap_R
  use stencil_sub, only: stencil_R
  implicit none
  integer,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                        ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                        ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4) &
                        ,is_table(Norb) &
                        ,NI,Nps,Nlma,ia_table(Nlma),Mps(NI),Jxyz(3,Nps,NI) &
                        ,nproc_Mxin_mul,irank_overlap(6),icomm_overlap,icomm_pseudo
  real(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) &
                        ,V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end,1:Nspin) &
                        ,lap0,lapt(4,3),uV(Nps,Nlma),uVu(Nlma)
  real(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: is,iorb

  if(nproc_Mxin_mul.ne.1) then
    call update_overlap_R(tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) &
                         ,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,4 &
                         ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap,icomm_overlap)
  end if

! stencil

  do iorb=1,Norb
    is = is_table(iorb)
    call stencil_R(tpsi(:,:,:,iorb),htpsi(:,:,:,iorb),ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                  ,V_local(:,:,:,is),ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                  ,idx,idy,idz,lap0,lapt)

  end do

! pseudopotential

  call pseudo_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
               ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
               ,nproc_Mxin_mul,icomm_pseudo)

  return
end subroutine hpsi_R

!===================================================================================================================================

SUBROUTINE hpsi_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb & ! wavefunctions
                 ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &         ! local potential
                 ,idx,idy,idz,lap0,lapt &                                           ! stencil
                 ,is_table,Nk &                                                     ! spin table & # of k points
                 ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &                            ! pseudo potential
                 ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo &         ! MPI
                 ,ik_table,nabt,kAc,exp_ikr,ttpsi)                                  ! k vector & ttpsi
  use update_overlap_sub, only: update_overlap_C
  use stencil_sub, only: stencil_C
  implicit none
  integer   ,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                           ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                           ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4) &
                           ,is_table(Norb),Nk &
                           ,NI,Nps,Nlma,ia_table(Nlma),Mps(NI),Jxyz(3,Nps,NI) &
                           ,nproc_Mxin_mul,irank_overlap(6),icomm_overlap,icomm_pseudo
  complex(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  real(8)   ,intent(in)  :: V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end,1:Nspin) &
                           ,lap0,lapt(4,3),uV(Nps,Nlma),uVu(Nlma)
  complex(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer   ,optional,intent(in)  :: ik_table(Norb)
  real(8)   ,optional,intent(in)  :: nabt(4,3),kAc(Nk,3)
  complex(8),optional,intent(in)  :: exp_ikr(Nps,NI,Nk)
  complex(8),optional,intent(out) :: ttpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: is,iorb,ik,ix,iy,iz
  real(8) :: k_nabt(4,3),k2
  logical :: if_kAc

  if_kAc = present(kAc)

  if(nproc_Mxin_mul.ne.1) then
    call update_overlap_C(tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) &
                       ,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,4 &
                       ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap,icomm_overlap)
  end if

! stencil

  do iorb=1,Norb
    is = is_table(iorb)

    if(if_kAc) then
      ik = ik_table(iorb)
      k2 = 0.5d0* sum(kAc(ik,:)**2)
      k_nabt(:,1) = kAc(ik,1) * nabt(:,1)
      k_nabt(:,2) = kAc(ik,2) * nabt(:,2)
      k_nabt(:,3) = kAc(ik,3) * nabt(:,3)
    else
      k2 = 0d0
      k_nabt = 0d0
    end if

    call stencil_C(tpsi(:,:,:,iorb),htpsi(:,:,:,iorb),ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                  ,V_local(:,:,:,is),ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                  ,idx,idy,idz,lap0+k2,lapt,k_nabt)

  end do

! subtraction

  if(present(ttpsi)) then
    do iorb=1,Norb
      is = is_table(iorb)
      do iz=iz_sta,iz_end
        do iy=iy_sta,iy_end
          do ix=ix_sta,ix_end
            ttpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) - V_local(ix,iy,iz,is) * tpsi(ix,iy,iz,iorb)
          end do
        end do
      end do
    end do
  end if

! pseudopotential

  call pseudo_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
               ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
               ,Nk,nproc_Mxin_mul,icomm_pseudo &
               ,ik_table,exp_ikr)

  return
end subroutine hpsi_C

!===================================================================================================================================

subroutine pseudo_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                   ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
                   ,nproc_Mxin_mul,icomm_pseudo)
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                        ,NI,Nps,Nlma,ia_table(Nlma),Mps(NI),Jxyz(3,Nps,NI) &
                        ,nproc_Mxin_mul,icomm_pseudo
  real(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  real(8),intent(in)  :: uV(Nps,Nlma),uVu(Nlma)
  !
  real(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: ilma,ia,j,ix,iy,iz,iorb
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)


  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))

    do iorb=1,Norb
      do ilma=1,Nlma
        ia = ia_table(ilma)
        uVpsi=0.d0
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
        end do
        uVpsi = uVpsi * uVu(ilma)
        uVpsibox(ilma,iorb) = uVpsi
      end do
    end do
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
    do iorb=1,Norb
      do ilma=1,Nlma
        ia = ia_table(ilma)
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          wrk = uVpsibox2(ilma,iorb) * uV(j,ilma)
          htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
        end do
      end do
    end do

    deallocate(uVpsibox,uVpsibox2)

  else

    do iorb=1,Norb
      do ilma=1,Nlma
        ia = ia_table(ilma)
        uVpsi=0.d0
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
        end do
        uVpsi = uVpsi * uVu(ilma)
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          wrk = uVpsi * uV(j,ilma)
          htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
        end do
      end do
    end do

  end if

  return
end subroutine pseudo_R

!===================================================================================================================================

subroutine pseudo_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                   ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
                   ,Nk,nproc_Mxin_mul,icomm_pseudo &
                   ,ik_table,exp_ikr)
  use salmon_communication, only: comm_summation
  implicit none
  integer   ,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                          ,NI,Nps,Nlma,ia_table(Nlma),Mps(NI),Jxyz(3,Nps,NI) &
                          ,Nk,nproc_Mxin_mul,icomm_pseudo
  complex(8),intent(in) :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  real(8)   ,intent(in) :: uV(Nps,Nlma),uVu(Nlma)
  !
  integer   ,optional,intent(in) :: ik_table(Norb)
  complex(8),optional,intent(in) :: exp_ikr(Nps,NI,Nk)
  !
  complex(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer    :: ilma,ia,j,ix,iy,iz,ik,iorb
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)


  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))
    if(present(exp_ikr)) then

      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * exp_ikr(j,ia,ik) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = conjg(exp_ikr(j,ia,ik)) * uVpsibox2(ilma,iorb) * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    else

      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = uVpsibox2(ilma,iorb) * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    end if
    deallocate(uVpsibox,uVpsibox2)

  else

    if(present(exp_ikr)) then

      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * exp_ikr(j,ia,ik) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = conjg(exp_ikr(j,ia,ik)) * uVpsi * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    else

      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = uVpsi * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    end if

  end if

  return
end subroutine pseudo_C

!===================================================================================================================================

end module hpsi_sub
