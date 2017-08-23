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
module hpsi_sub
complex(8), parameter :: zI=(0.d0,1.d0)

interface hpsi
  module procedure hpsi_R,hpsi_C
end interface hpsi

contains

!==================================================================================================

SUBROUTINE hpsi_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                 ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                 ,idx,idy,idz,lap0,lapt,is_table,nproc_Mxin_mul)
  use scf_data, only: iflag_ps ! GCEED
  use update_overlap_sub, only: update_overlap_R
  implicit none
  integer,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                        ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                        ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4) &
                        ,is_table(Norb),nproc_Mxin_mul
  real(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) &
                        ,V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end,1:Nspin) &
                        ,lap0,lapt(4,3)
  real(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: is,iorb

  if(nproc_Mxin_mul.ne.1) then
    call update_overlap_R(tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) &
                         ,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,4 &
                         ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end)
  end if

! stencil

  do iorb=1,Norb
    is = is_table(iorb)
    call stencil_R(tpsi(:,:,:,iorb),htpsi(:,:,:,iorb),ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                  ,V_local(:,:,:,is),ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                  ,idx,idy,idz,lap0,lapt)

  end do

! pseudopotential

  if(iflag_ps.eq.1)then
    call pseudo_GCEED_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul)
  end if

  return
end subroutine hpsi_R

!==================================================================================================

SUBROUTINE hpsi_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                 ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                 ,idx,idy,idz,lap0,lapt,is_table,Nk,nproc_Mxin_mul &
                 ,ik_table,nabt,kAc,exp_ikr,ttpsi)
  use scf_data, only: iflag_ps ! GCEED
  use Global_Variables, only: Nps,NI ! ARTED
  use update_overlap_sub, only: update_overlap_C
  implicit none
  integer   ,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                           ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,Nspin &
                           ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4) &
                           ,is_table(Norb),Nk,nproc_Mxin_mul
  complex(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  real(8)   ,intent(in)  :: V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end,1:Nspin) &
                           ,lap0,lapt(4,3)
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
                       ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end)
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
            ttpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) - V_local(ix,iy,iz,is) * tpsi(iz,iy,iz,iorb)
          end do
        end do
      end do
    end do
  end if

! pseudopotential

  if(if_kAc) then
    do iorb=1,Norb
      ik = ik_table(iorb)
      call pseudo_ARTED(tpsi(:,:,:,iorb),htpsi(:,:,:,iorb) &
                       ,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                       ,exp_ikr(:,:,ik),Nps,NI)
    end do
  else
    if(iflag_ps.eq.1)then
      call pseudo_GCEED_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul)
    end if
  end if

  return
end subroutine hpsi_C

!==================================================================================================

# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

subroutine stencil_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                    ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                    ,idx,idy,idz,lap0,lapt)
  implicit none
  integer,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                        ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                        ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4)
  real(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end) &
                        ,V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end),lap0,lapt(4,3)
  real(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end)
  !
  integer :: iz,iy,ix
  real(8) :: v

!$OMP parallel
!$OMP do private(iz,iy,ix,v)
  do iz=iz_sta,iz_end
  do iy=iy_sta,iy_end
  do ix=ix_sta,ix_end

    v =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
        +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
        +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
        +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
        +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
        +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
        +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4))) + v

    v =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
        +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
        +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
        +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4))) + v

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_R

subroutine stencil_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                    ,V_local,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                    ,idx,idy,idz,lap0,lapt,nabt)
  implicit none
  integer   ,intent(in)  :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
                           ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
                           ,idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4)
  complex(8),intent(in)  :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end)
  real(8)   ,intent(in)  :: V_local(ix_sta:ix_end,iy_sta:iy_end,iz_sta:iz_end),lap0,lapt(4,3),nabt(4,3)
  complex(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end)
  !
  integer :: iz,iy,ix
  complex(8) :: v,w

!$OMP parallel
!$OMP do private(iz,iy,ix,v,w)
  do iz=iz_sta,iz_end
  do iy=iy_sta,iy_end
  do ix=ix_sta,ix_end

    v =  lapt(1,1)*(tpsi(DX(1)) + tpsi(DX(-1))) &
        +lapt(2,1)*(tpsi(DX(2)) + tpsi(DX(-2))) &
        +lapt(3,1)*(tpsi(DX(3)) + tpsi(DX(-3))) &
        +lapt(4,1)*(tpsi(DX(4)) + tpsi(DX(-4)))

    v =  lapt(1,2)*(tpsi(DY(1)) + tpsi(DY(-1))) &
        +lapt(2,2)*(tpsi(DY(2)) + tpsi(DY(-2))) &
        +lapt(3,2)*(tpsi(DY(3)) + tpsi(DY(-3))) &
        +lapt(4,2)*(tpsi(DY(4)) + tpsi(DY(-4))) + v

    v =  lapt(1,3)*(tpsi(DZ(1)) + tpsi(DZ(-1))) &
        +lapt(2,3)*(tpsi(DZ(2)) + tpsi(DZ(-2))) &
        +lapt(3,3)*(tpsi(DZ(3)) + tpsi(DZ(-3))) &
        +lapt(4,3)*(tpsi(DZ(4)) + tpsi(DZ(-4))) + v

    w =  nabt(1,1)*(tpsi(DX(1)) - tpsi(DX(-1))) &
        +nabt(2,1)*(tpsi(DX(2)) - tpsi(DX(-2))) &
        +nabt(3,1)*(tpsi(DX(3)) - tpsi(DX(-3))) &
        +nabt(4,1)*(tpsi(DX(4)) - tpsi(DX(-4)))

    w =  nabt(1,2)*(tpsi(DY(1)) - tpsi(DY(-1))) &
        +nabt(2,2)*(tpsi(DY(2)) - tpsi(DY(-2))) &
        +nabt(3,2)*(tpsi(DY(3)) - tpsi(DY(-3))) &
        +nabt(4,2)*(tpsi(DY(4)) - tpsi(DY(-4))) + w

    w =  nabt(1,3)*(tpsi(DZ(1)) - tpsi(DZ(-1))) &
        +nabt(2,3)*(tpsi(DZ(2)) - tpsi(DZ(-2))) &
        +nabt(3,3)*(tpsi(DZ(3)) - tpsi(DZ(-3))) &
        +nabt(4,3)*(tpsi(DZ(4)) - tpsi(DZ(-4))) + w

    htpsi(ix,iy,iz) = ( V_local(ix,iy,iz) + lap0 )*tpsi(ix,iy,iz) - 0.5d0 * v - zI * w
  end do
  end do
  end do
!$OMP end do
!$OMP end parallel

  return
end subroutine stencil_C

!==================================================================================================

subroutine pseudo_ARTED(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,exp_ikr,Nps,NI)
  use Global_Variables, only: Mps,uV,iuV,Hxyz,Nlma,a_tbl,Jxyz,NLy,NLz ! ARTED
  implicit none
  integer,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Nps,NI
  complex(8),dimension(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end) :: tpsi,htpsi
  complex(8),intent(in) :: exp_ikr(Nps,NI)
  !
  integer    :: ilma,ia,j,i,ix,iy,iz ! i = ix*NLy*NLz + iy*NLz + iz + 1, ix=0:NLx-1, ...
  complex(8) :: uVpsi,wrk
  do ilma=1,Nlma
    ia=a_tbl(ilma)
    uVpsi=0.d0
    do j=1,Mps(ia)
      i = Jxyz(j,ia)
      iz = mod(i-1,NLz)
      iy = mod((i-1-iz)/NLz,NLy)
      ix = (i-1-iz-iy*NLy)/(NLy*NLz)
      uVpsi = uVpsi + uV(j,ilma) * exp_ikr(j,ia) * tpsi(iz,iy,ix) ! tpsi(0:NLz-1,0:NLy-1,0:NLx-1) ! x <--> z
    end do
    uVpsi=uVpsi*Hxyz*iuV(ilma)
    do j=1,Mps(ia)
      i = Jxyz(j,ia)
      iz = mod(i-1,NLz)
      iy = mod((i-1-iz)/NLz,NLy)
      ix = (i-1-iz-iy*NLy)/(NLy*NLz)
      wrk = conjg(exp_ikr(j,ia)) * uVpsi * uV(j,ilma)
      htpsi(iz,iy,ix) = htpsi(iz,iy,ix) + wrk ! x <--> z
    end do
  end do
end subroutine pseudo_ARTED

!==================================================================================================

subroutine pseudo_GCEED_R(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul)
  use scf_data, only: MI,maxlm,Kion,Mlps,uVu,iwk_size,max_jMps_l,uV,Jxyz,jMps_l,max_jMps_l_s,jMps_l_s,Hvol,uVu ! GCEED
  use salmon_parallel, only: nproc_group_orbital, nproc_group_h
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul
  real(8),dimension(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) :: tpsi,htpsi

  integer :: iatom,jj,lm,ik
  integer :: iorb

  real(8) :: sumbox

  real(8), allocatable :: uVpsibox(:,:,:)
  real(8), allocatable :: uVpsibox2(:,:,:)

  allocate (uVpsibox(1:maxlm,1:MI,1:Norb))
  allocate (uVpsibox2(1:maxlm,1:MI,1:Norb))

  if(nproc_Mxin_mul==1)then
!$OMP parallel do private(iatom,lm)
    do iorb=1,Norb
    do iatom=1,MI
      do lm=1,maxlm
        uVpsibox2(lm,iatom,iorb)=0.d0
      end do
    end do
    end do
  else
!$OMP parallel do private(iatom,lm)
    do iorb=1,Norb
    do iatom=1,MI
      do lm=1,maxlm
        uVpsibox(lm,iatom,iorb)=0.d0
      end do
    end do
    end do
  end if

  if(nproc_Mxin_mul==1)then
    do iorb=1,Norb
    do iatom=1,MI
      ik=Kion(iatom)
      loop_lm : do lm=1,(Mlps(ik)+1)**2
        if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm
        sumbox=0.d0
        if(iwk_size>=1.and.iwk_size<=2)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l(iatom)
            sumbox=sumbox+uV(jMps_l(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l(jj,iatom),iatom),iorb)
          end do
        else if(iwk_size>=11.and.iwk_size<=12)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l_s(iatom)
            sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l_s(jj,iatom),iatom),iorb)
          end do
        end if
        uVpsibox2(lm,iatom,iorb)=sumbox*Hvol/uVu(lm,iatom)
      end do loop_lm
    end do
    end do
  else
    do iorb=1,Norb
    do iatom=1,MI
      ik=Kion(iatom)
      loop_lm2 : do lm=1,(Mlps(ik)+1)**2
        if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
        sumbox=0.d0
        if(iwk_size>=1.and.iwk_size<=2)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l(iatom)
            sumbox=sumbox+uV(jMps_l(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l(jj,iatom),iatom),iorb)
          end do
        else if(iwk_size>=11.and.iwk_size<=12)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l_s(iatom)
            sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l_s(jj,iatom),iatom),iorb)
          end do
        end if
        uVpsibox(lm,iatom,iorb)=sumbox*Hvol/uVu(lm,iatom)
      end do loop_lm2
    end do
    end do
  end if

  if(nproc_Mxin_mul==1)then
    continue
  else
    if(iwk_size>=1.and.iwk_size<=2)then
      call comm_summation(uVpsibox,uVpsibox2,maxlm*MI*Norb,nproc_group_orbital)
    else if(iwk_size>=11.and.iwk_size<=12)then
      call comm_summation(uVpsibox,uVpsibox2,maxlm*MI*Norb,nproc_group_h)
    end if
  end if


  if(iwk_size>=1.and.iwk_size<=2)then
    do iorb=1,Norb
    do iatom=1,MI
      ik=Kion(iatom)
      do jj=1,max_jMps_l(iatom)
        do lm=1,(Mlps(ik)+1)**2
          htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom),iorb)= &
            htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom),iorb) + &
            uVpsibox2(lm,iatom,iorb)*uV(jMps_l(jj,iatom),lm,iatom)
        end do
      end do
    end do
    end do
  else if(iwk_size>=11.and.iwk_size<=12)then
    do iorb=1,Norb
    do iatom=1,MI
      ik=Kion(iatom)
      do jj=1,max_jMps_l_s(iatom)
        do lm=1,(Mlps(ik)+1)**2
          htpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),Jxyz(3,jMps_l_s(jj,iatom),iatom),iorb)= &
            htpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),Jxyz(3,jMps_l_s(jj,iatom),iatom),iorb) + &
            uVpsibox2(lm,iatom,iorb)*uV(jMps_l_s(jj,iatom),lm,iatom)
        end do
      end do
    end do
    end do
  end if
  deallocate(uVpsibox2,uVpsibox)
  return
end subroutine pseudo_GCEED_R

!==================================================================================================

subroutine pseudo_GCEED_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul)

! GCEED
!-------------------------------------------------------------------------------------
  use scf_data, only: MI,maxlm,Kion,Mlps,uVu,max_jMps_l_s,uV,uV2nd,Jxyz,Jxyz2nd,Hvol,jMps_l,jMps_l_s &
                     ,Mps3rd_2nd,iatomnum_ps_2nd,ikind_eext,icalcforce,iflag_md,iwk_size,max_jMps_l &
                     ,maxMps_all,numatom_ps_2nd,jja,Jxyz_all
!-------------------------------------------------------------------------------------

  use salmon_parallel, only: nproc_group_orbital, nproc_group_h
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,nproc_Mxin_mul
  complex(8),dimension(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb) :: tpsi,htpsi
  integer :: ja,jj,iatom,lm,ikoa,iix,iiy,iiz
  integer :: iorb
  complex(8) :: sumbox
  complex(8),allocatable :: uVpsibox3(:,:,:,:)
  complex(8),allocatable :: uVpsibox4(:,:,:,:)

  allocate (uVpsibox3(1:maxlm,1:MI,1:Norb,1))
  allocate (uVpsibox4(1:maxlm,1:MI,1:Norb,1))

!-------------------------------------------------------------------------------------
  do iorb=1,Norb
!$OMP parallel
!$OMP do private(iatom,jj)
    do iatom=1,MI
    do jj=1,maxlm
      uVpsibox3(jj,iatom,iorb,1)=0.d0
      uVpsibox4(jj,iatom,iorb,1)=0.d0
    end do
    end do
!$OMP end parallel
  end do
!-------------------------------------------------------------------------------------

!$OMP parallel
  if(nproc_Mxin_mul==1)then

! uV2nd, Jxyz2nd <--> hpsi2.f90
!-------------------------------------------------------------------------------------
    do iorb=1,Norb
!$OMP do private(iatom,ikoa,lm,jj,sumbox) schedule(static, 1)
      do iatom=1,MI
        ikoa=Kion(iatom)
        loop_lm : do lm=1,(Mlps(ikoa)+1)**2
          if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm
          sumbox=0.d0
          do jj=1,max_jMps_l(iatom)
            sumbox=sumbox+uV2nd(jj,lm,iatom)*tpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iorb)
          end do
          uVpsibox4(lm,iatom,iorb,1)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm
      end do
!$OMP end do nowait
    end do
!-------------------------------------------------------------------------------------

  else

! common uV,Jxyz,etc. @ hpsi2.f90 (Norb==1)
!-------------------------------------------------------------------------------------
    do iorb=1,Norb
!$OMP do private(iatom,ikoa,lm,jj,sumbox) schedule(static, 1)
      do iatom=1,MI
        ikoa=Kion(iatom)
        loop_lm2 : do lm=1,(Mlps(ikoa)+1)**2
          if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
          sumbox=0.d0
          if(iwk_size>=1.and.iwk_size<=2)then
            do jj=1,max_jMps_l(iatom)
              sumbox=sumbox+uV(jMps_l(jj,iatom),lm,iatom)*  &
                       tpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),  &
                            Jxyz(3,jMps_l(jj,iatom),iatom),iorb)
            end do
          else if(iwk_size>=11.and.iwk_size<=12)then
            do jj=1,max_jMps_l_s(iatom)
              sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                       tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                            Jxyz(3,jMps_l_s(jj,iatom),iatom),iorb)
            end do
          end if
          uVpsibox3(lm,iatom,iorb,1)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm2
      end do
!$OMP end do nowait
    end do
    if(iwk_size>=1.and.iwk_size<=2)then
!      elp3(705)=get_wtime()
      call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*Norb,nproc_group_orbital)
!      elp3(706)=get_wtime()
!      elp3(744)=elp3(744)+elp3(706)-elp3(705)
    else if(iwk_size>=11.and.iwk_size<=12)then
      call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*Norb,nproc_group_h)
    end if
!-------------------------------------------------------------------------------------

  end if
!$OMP end parallel

! Pseudopotential 2 (non-local)
  if(ikind_eext==0.and.icalcforce==0.and.iflag_md==0)then

!-------------------------------------------------------------------------------------
!$OMP parallel
    do iorb=1,Norb
!$OMP do private(jj,iatom,ja,iix,iiy,iiz,ikoa,lm)
      do jj=1,maxMps_all
        iix=Jxyz_all(1,jj)
        iiy=Jxyz_all(2,jj)
        iiz=Jxyz_all(3,jj)
        do iatom=1,numatom_ps_2nd(jj)
          ja=jja(iatom,jj)
          ikoa=Kion(iatomnum_ps_2nd(ja))
          do lm=1,(Mlps(ikoa)+1)**2
          htpsi(iix,iiy,iiz,iorb)= &
            htpsi(iix,iiy,iiz,iorb) + &
            uVpsibox4(lm,iatomnum_ps_2nd(ja),iorb,1)*uV(Mps3rd_2nd(ja),lm,iatomnum_ps_2nd(ja))
          end do
        end do
      end do
!$OMP end do nowait
    end do
!$OMP end parallel
!-------------------------------------------------------------------------------------

  else

!-------------------------------------------------------------------------------------
    do iorb=1,Norb
      do iatom=1,MI
      ikoa=Kion(iatom)
        do jj=1,max_jMps_l(iatom)
          do lm=1,(Mlps(ikoa)+1)**2
            htpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iorb)= &
              htpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iorb) + &
                uVpsibox4(lm,iatom,iorb,1)*uV2nd(jj,lm,iatom)
          end do
        end do
      end do
    end do
!-------------------------------------------------------------------------------------

  end if

  deallocate(uVpsibox3,uVpsibox4)
  return
end subroutine pseudo_GCEED_C

end module hpsi_sub
