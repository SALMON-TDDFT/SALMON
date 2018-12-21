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
!=======================================================================
!========================= Hamiltonian Operation (for complex funcitons)

SUBROUTINE hpsi_groupob(tpsi,htpsi,tpsi_out,tVlocal,nn,isub)
use salmon_parallel, only: nproc_group_korbital, nproc_group_h
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use gradient2_sub
use allocate_mat_sub
use sendrecv_groupob_sub
use sendrecv_groupob_tmp_sub
use init_sendrecv_sub
!use sendrecv_groupob_ngp_sub
implicit none
complex(8) :: tpsi(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
complex(8) :: htpsi(iwk2sta(1):iwk2end(1)+1,  &
                    iwk2sta(2):iwk2end(2),      &
                    iwk2sta(3):iwk2end(3),     &
                   1:iobnum,k_sta:k_end)
complex(8) :: tpsi_out(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,k_sta:k_end)
real(8) :: tVlocal(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin)

integer :: ist,ix,iy,iz,jj,iatom,lm,ikoa,iob,j,ind
integer :: nn,isub,iik
integer :: jspin

complex(8) :: sumbox

complex(8), parameter :: zi=(0.d0,1.d0)
complex(8),allocatable :: uVpsibox3(:,:,:,:)
complex(8),allocatable :: uVpsibox4(:,:,:,:)

complex(8) :: ekr(maxMps,MI,k_sta:k_end)
real(8) :: x,y,z

complex(8) :: fdN0
complex(8) :: fdN1(0:12,3)
complex(8) :: fdN2(0:12,3)

real(8) :: f0
integer :: iobmax
integer :: iob_allob

f0=(1.d0/Hgs(1)**2   &
   +1.d0/Hgs(2)**2   &
   +1.d0/Hgs(3)**2)

call calc_pmax(iobmax)

allocate (uVpsibox3(1:maxlm,1:MI,1:iobmax,k_sta:k_end))
allocate (uVpsibox4(1:maxlm,1:MI,1:iobmax,k_sta:k_end))

if(iperiodic==3)then
  do iatom=1,MI
    ikoa=Kion(iatom)
    do jj=1,Mps(iatom)
      x=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
      y=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
      z=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
      do iik=k_sta,k_end
        ekr(jj,iatom,iik)=exp(zi*(k_rd(1,iik)*x+k_rd(2,iik)*y+k_rd(3,iik)*z))
      end do
    end do
  end do
end if

elp3(701)=get_wtime()

select case(iperiodic)
case(0)
  if(nproc_Mxin_mul/=1)then
    call sendrecv_groupob(tpsi)
  end if
case(3)
  call sendrecv_groupob_tmp(tpsi)
end select

elp3(704)=get_wtime()
elp3(743)=elp3(743)+elp3(704)-elp3(701)

select case(iperiodic)
case(0)
  if(Nd==4)then
    fdN0=-0.5d0*cNmat(0,Nd)*f0
    do j=1,3
      do ind=1,4
        fdN1(ind,j)=-0.5d0*cNmat(ind,Nd)/Hgs(j)**2
      end do
    end do
    do iik=k_sta,k_end
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      call set_ispin(iob_allob,jspin)
!$OMP parallel private(iz)
      do iz=iwk3sta(3),iwk3end(3)
!$OMP do private(iy,ix)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          htpsi(ix,iy,iz,iob,iik) =   &
            ( tVlocal(ix,iy,iz,jspin)+fdN0)*tpsi(ix,iy,iz,iob,iik)  &
            +fdN1(1,1)*(tpsi(ix+1,iy,iz,iob,iik) + tpsi(ix-1,iy,iz,iob,iik))  &
            +fdN1(2,1)*(tpsi(ix+2,iy,iz,iob,iik) + tpsi(ix-2,iy,iz,iob,iik))  &
            +fdN1(3,1)*(tpsi(ix+3,iy,iz,iob,iik) + tpsi(ix-3,iy,iz,iob,iik))  &
            +fdN1(4,1)*(tpsi(ix+4,iy,iz,iob,iik) + tpsi(ix-4,iy,iz,iob,iik))  &
            +fdN1(1,2)*(tpsi(ix,iy+1,iz,iob,iik) + tpsi(ix,iy-1,iz,iob,iik))  &
            +fdN1(2,2)*(tpsi(ix,iy+2,iz,iob,iik) + tpsi(ix,iy-2,iz,iob,iik))  &
            +fdN1(3,2)*(tpsi(ix,iy+3,iz,iob,iik) + tpsi(ix,iy-3,iz,iob,iik))  &
            +fdN1(4,2)*(tpsi(ix,iy+4,iz,iob,iik) + tpsi(ix,iy-4,iz,iob,iik))  
        end do
        end do
!$OMP end do nowait
!$OMP do private(iy,ix)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          htpsi(ix,iy,iz,iob,iik) = htpsi(ix,iy,iz,iob,iik)  &
            +fdN1(1,3)*(tpsi(ix,iy,iz+1,iob,iik) + tpsi(ix,iy,iz-1,iob,iik))  &
            +fdN1(2,3)*(tpsi(ix,iy,iz+2,iob,iik) + tpsi(ix,iy,iz-2,iob,iik))  &
            +fdN1(3,3)*(tpsi(ix,iy,iz+3,iob,iik) + tpsi(ix,iy,iz-3,iob,iik))  &
            +fdN1(4,3)*(tpsi(ix,iy,iz+4,iob,iik) + tpsi(ix,iy,iz-4,iob,iik))
        end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel
    end do
    end do
  else
    do iik=k_sta,k_end
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      call set_ispin(iob_allob,jspin)
      fdN0=-0.5d0*cNmat(0,Nd)*f0
      do j=1,3
        do ind=1,Nd
          fdN1(ind,j)=-0.5d0*cNmat(ind,Nd)/Hgs(j)**2
        end do
      end do
!$OMP parallel private(iz,ist)
      do iz=iwk3sta(3),iwk3end(3)
!$OMP do private(iy,ix)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          htpsi(ix,iy,iz,iob,iik) = (tVlocal(ix,iy,iz,jspin)+fdN0) *tpsi(ix,iy,iz,iob,iik)  
          do ist=1,Nd  
            htpsi(ix,iy,iz,iob,iik) = htpsi(ix,iy,iz,iob,iik)         &
              +fdN1(ist,1)* (tpsi(ix+ist,iy,iz,iob,iik) + tpsi(ix-ist,iy,iz,iob,iik))     &
              +fdN1(ist,2)* (tpsi(ix,iy+ist,iz,iob,iik) + tpsi(ix,iy-ist,iz,iob,iik))
          end do 
        end do
        end do
!$OMP end do nowait
!$OMP do private(iy,ix)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          do ist=1,Nd  
            htpsi(ix,iy,iz,iob,iik) = htpsi(ix,iy,iz,iob,iik)  &
              +fdN1(ist,3)* (tpsi(ix,iy,iz+ist,iob,iik) + tpsi(ix,iy,iz-ist,iob,iik))
          end do 
        end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel
    end do
    end do
  end if
case(3)
  do iik=k_sta,k_end
    fdN0=0.5d0*ksquare(iik)-0.5d0*cNmat(0,Nd)*f0
    do jj=1,3
      do ind=1,4
        fdN1(ind,jj)=-0.5d0*cNmat(ind,Nd)/Hgs(jj)**2-zi*k_rd(jj,iik)*bNmat(ind,Nd)/Hgs(jj)
        fdN2(ind,jj)=-0.5d0*cNmat(ind,Nd)/Hgs(jj)**2+zi*k_rd(jj,iik)*bNmat(ind,Nd)/Hgs(jj)
      end do
    end do
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      call set_ispin(iob_allob,jspin)
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=iwk3sta(3),iwk3end(3)
      do iy=iwk3sta(2),iwk3end(2)
      do ix=iwk3sta(1),iwk3end(1)
        htpsi(ix,iy,iz,iob,iik) =                                        &
        ( tVlocal(ix,iy,iz,jspin)+fdN0)*tpsi(ix,iy,iz,iob,iik)  &
          +fdN1(1,1)* tpsi(ix+1,iy,iz,iob,iik) + fdN2(1,1)* tpsi(ix-1,iy,iz,iob,iik)  &
          +fdN1(2,1)* tpsi(ix+2,iy,iz,iob,iik) + fdN2(2,1)* tpsi(ix-2,iy,iz,iob,iik)  &
          +fdN1(3,1)* tpsi(ix+3,iy,iz,iob,iik) + fdN2(3,1)* tpsi(ix-3,iy,iz,iob,iik)  &
          +fdN1(4,1)* tpsi(ix+4,iy,iz,iob,iik) + fdN2(4,1)* tpsi(ix-4,iy,iz,iob,iik)  &
          +fdN1(1,2)* tpsi(ix,iy+1,iz,iob,iik) + fdN2(1,2)* tpsi(ix,iy-1,iz,iob,iik)  &
          +fdN1(2,2)* tpsi(ix,iy+2,iz,iob,iik) + fdN2(2,2)* tpsi(ix,iy-2,iz,iob,iik)  &
          +fdN1(3,2)* tpsi(ix,iy+3,iz,iob,iik) + fdN2(3,2)* tpsi(ix,iy-3,iz,iob,iik)  &
          +fdN1(4,2)* tpsi(ix,iy+4,iz,iob,iik) + fdN2(4,2)* tpsi(ix,iy-4,iz,iob,iik)  &
          +fdN1(1,3)* tpsi(ix,iy,iz+1,iob,iik) + fdN2(1,3)* tpsi(ix,iy,iz-1,iob,iik)  &
          +fdN1(2,3)* tpsi(ix,iy,iz+2,iob,iik) + fdN2(2,3)* tpsi(ix,iy,iz-2,iob,iik)  &
          +fdN1(3,3)* tpsi(ix,iy,iz+3,iob,iik) + fdN2(3,3)* tpsi(ix,iy,iz-3,iob,iik)  &
          +fdN1(4,3)* tpsi(ix,iy,iz+4,iob,iik) + fdN2(4,3)* tpsi(ix,iy,iz-4,iob,iik)
      end do
      end do
      end do
    end do
  end do
end select



! Pseudopotential 1 (non-local)
if(iflag_ps.eq.1)then

  do iik=k_sta,k_end
  do iob=1,iobmax
!$OMP parallel
!$OMP do private(iatom,jj)
    do iatom=1,MI
    do jj=1,maxlm
      uVpsibox3(jj,iatom,iob,iik)=0.d0
      uVpsibox4(jj,iatom,iob,iik)=0.d0
    end do
    end do
!$OMP end parallel
  end do
  end do

  select case(iperiodic)
  case(0)
    do iik=k_sta,k_end
    do iob=1,iobmax
!$OMP parallel
!$OMP do private(iatom,ikoa,lm,jj,sumbox) schedule(static, 1)
      do iatom=1,MI
        ikoa=Kion(iatom)
        loop_lm2 : do lm=1,(Mlps(ikoa)+1)**2
          if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
          sumbox=0.d0
          do jj=1,Mps(iatom)
            sumbox=sumbox+uV(jj,lm,iatom)*  &
                     tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
          end do
          uVpsibox3(lm,iatom,iob,iik)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm2
      end do
!$OMP end do nowait
!$OMP end parallel
    end do
    end do
  case(3)
    do iik=k_sta,k_end
    do iob=1,iobmax
!$OMP parallel
!$OMP do private(iatom,ikoa,lm,jj,sumbox) schedule(static, 1)
      do iatom=1,MI
        ikoa=Kion(iatom)
        loop_lm4 : do lm=1,(Mlps(ikoa)+1)**2
          if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm4
          sumbox=0.d0
          do jj=1,Mps(iatom)
            sumbox=sumbox+uV(jj,lm,iatom)*tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)*ekr(jj,iatom,iik)
          end do
          uVpsibox3(lm,iatom,iob,iik)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm4
      end do
!$OMP end do nowait
!$OMP end parallel
    end do
    end do
  end select

  elp3(705)=get_wtime()
  call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*iobmax*k_num,nproc_group_korbital)
  elp3(706)=get_wtime()
  elp3(744)=elp3(744)+elp3(706)-elp3(705)

end if

! Pseudopotential 2 (non-local)
if(iflag_ps==1) then
  select case(iperiodic)
  case(0)
    do iik=k_sta,k_end
    do iob=1,iobmax
      do iatom=1,MI
        ikoa=Kion(iatom)
!$OMP parallel do private(jj,lm)
        do jj=1,Mps(iatom)
          do lm=1,(Mlps(ikoa)+1)**2
            htpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)= &
              htpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik) + &
                uVpsibox4(lm,iatom,iob,iik)*uV(jj,lm,iatom)
          end do
        end do
      end do
    end do
    end do
  case(3)
    do iik=k_sta,k_end
    do iob=1,iobmax
      do iatom=1,MI
        ikoa=Kion(iatom)
!$OMP parallel do private(jj,lm)
        do jj=1,Mps(iatom)
          do lm=1,(Mlps(ikoa)+1)**2
            htpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)= &
              htpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik) + &
                uVpsibox4(lm,iatom,iob,iik)*uV(jj,lm,iatom)*conjg(ekr(jj,iatom,iik))
          end do
        end do
      end do
    end do
    end do
  end select
end if

if(isub==1)then
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
end if

return

END SUBROUTINE hpsi_groupob
