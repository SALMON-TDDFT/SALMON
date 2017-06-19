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
use salmon_parallel, only: nproc_group_orbital, nproc_group_h
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use gradient2_sub
use allocate_mat_sub
use allocate_sendrecv_groupob_sub
use sendrecv_groupob_sub
use init_sendrecv_sub
!use sendrecv_groupob_ngp_sub
implicit none
complex(8) :: tpsi(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,1)
complex(8) :: htpsi(iwk2sta(1):iwk2end(1)+1,  &
                    iwk2sta(2):iwk2end(2),      &
                    iwk2sta(3):iwk2end(3),     &
                   1:iobnum,1)
complex(8) :: tpsi_out(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,1)
real(8) :: tVlocal(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),numspin)

integer :: ist,ix,iy,iz,jj,iatom,lm,ikoa,iob,j,ind
integer :: ja
integer :: nn,isub
integer :: jspin
integer :: iix,iiy,iiz

complex(8) :: sumbox

complex(8), parameter :: zi=(0.d0,1.d0)
complex(8),allocatable :: uVpsibox3(:,:,:,:)
complex(8),allocatable :: uVpsibox4(:,:,:,:)

real(8) :: fdN0
real(8) :: fdN1(0:12,3)

real(8) :: f0
integer :: iobmax
integer :: iob_allob

f0=(1.d0/Hgs(1)**2   &
   +1.d0/Hgs(2)**2   &
   +1.d0/Hgs(3)**2)

call calc_pmax(iobmax)

allocate (uVpsibox3(1:maxlm,1:MI,1:iobmax,1))
allocate (uVpsibox4(1:maxlm,1:MI,1:iobmax,1))


elp3(701)=get_wtime()

if(nproc_Mxin_mul/=1)then
  call sendrecv_groupob(tpsi)
end if

elp3(704)=get_wtime()
elp3(743)=elp3(743)+elp3(704)-elp3(701)

if(Nd==4)then
  fdN0=-0.5d0*cNmat(0,Nd)*f0
  do j=1,3
    do ind=1,4
      fdN1(ind,j)=-0.5d0*cNmat(ind,Nd)/Hgs(j)**2
    end do
  end do
  do iob=1,iobmax
    call calc_allob(iob,iob_allob)
    call set_ispin(iob_allob,jspin)
!$OMP parallel private(iz)
    do iz=iwk3sta(3),iwk3end(3)
!$OMP do
      do iy=iwk3sta(2),iwk3end(2)
      do ix=iwk3sta(1),iwk3end(1)
        htpsi(ix,iy,iz,iob,1) =   &
          ( tVlocal(ix,iy,iz,jspin)+fdN0)*tpsi(ix,iy,iz,iob,1)  &
          +fdN1(1,1)*(tpsi(ix+1,iy,iz,iob,1) + tpsi(ix-1,iy,iz,iob,1))  &
          +fdN1(2,1)*(tpsi(ix+2,iy,iz,iob,1) + tpsi(ix-2,iy,iz,iob,1))  &
          +fdN1(3,1)*(tpsi(ix+3,iy,iz,iob,1) + tpsi(ix-3,iy,iz,iob,1))  &
          +fdN1(4,1)*(tpsi(ix+4,iy,iz,iob,1) + tpsi(ix-4,iy,iz,iob,1))  &
          +fdN1(1,2)*(tpsi(ix,iy+1,iz,iob,1) + tpsi(ix,iy-1,iz,iob,1))  &
          +fdN1(2,2)*(tpsi(ix,iy+2,iz,iob,1) + tpsi(ix,iy-2,iz,iob,1))  &
          +fdN1(3,2)*(tpsi(ix,iy+3,iz,iob,1) + tpsi(ix,iy-3,iz,iob,1))  &
          +fdN1(4,2)*(tpsi(ix,iy+4,iz,iob,1) + tpsi(ix,iy-4,iz,iob,1))  
      end do
      end do
!$OMP end do nowait
!$OMP do
      do iy=iwk3sta(2),iwk3end(2)
      do ix=iwk3sta(1),iwk3end(1)
        htpsi(ix,iy,iz,iob,1) = htpsi(ix,iy,iz,iob,1)  &
          +fdN1(1,3)*(tpsi(ix,iy,iz+1,iob,1) + tpsi(ix,iy,iz-1,iob,1))  &
          +fdN1(2,3)*(tpsi(ix,iy,iz+2,iob,1) + tpsi(ix,iy,iz-2,iob,1))  &
          +fdN1(3,3)*(tpsi(ix,iy,iz+3,iob,1) + tpsi(ix,iy,iz-3,iob,1))  &
          +fdN1(4,3)*(tpsi(ix,iy,iz+4,iob,1) + tpsi(ix,iy,iz-4,iob,1))
      end do
      end do
!$OMP end do nowait
    end do
!$OMP end parallel
  end do
else
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
!$OMP do
      do iy=iwk3sta(2),iwk3end(2)
      do ix=iwk3sta(1),iwk3end(1)
        htpsi(ix,iy,iz,iob,1) = (tVlocal(ix,iy,iz,jspin)+fdN0) *tpsi(ix,iy,iz,iob,1)  
        do ist=1,Nd  
          htpsi(ix,iy,iz,iob,1) = htpsi(ix,iy,iz,iob,1)         &
            +fdN1(ist,1)* (tpsi(ix+ist,iy,iz,iob,1) + tpsi(ix-ist,iy,iz,iob,1))     &
            +fdN1(ist,2)* (tpsi(ix,iy+ist,iz,iob,1) + tpsi(ix,iy-ist,iz,iob,1))
        end do 
      end do
      end do
!$OMP end do nowait
!$OMP do
      do iy=iwk3sta(2),iwk3end(2)
      do ix=iwk3sta(1),iwk3end(1)
        do ist=1,Nd  
          htpsi(ix,iy,iz,iob,1) = htpsi(ix,iy,iz,iob,1)  &
            +fdN1(ist,3)* (tpsi(ix,iy,iz+ist,iob,1) + tpsi(ix,iy,iz-ist,iob,1))
        end do 
      end do
      end do
!$OMP end do nowait
    end do
!$OMP end parallel
  end do
end if

! Pseudopotential 1 (non-local)
if(iflag_ps.eq.1)then

  do iob=1,iobmax
!$OMP parallel
!$OMP do private(iatom,jj)
    do iatom=1,MI
    do jj=1,maxlm
      uVpsibox3(jj,iatom,iob,1)=0.d0
      uVpsibox4(jj,iatom,iob,1)=0.d0
    end do
    end do
!$OMP end parallel
  end do

  do iob=1,iobmax
!$OMP parallel
    if(nproc_Mxin_mul==1)then
!$OMP do private(iatom,ikoa,lm,jj,sumbox) schedule(static, 1)
      do iatom=1,MI
        ikoa=Kion(iatom)
        loop_lm : do lm=1,(Mlps(ikoa)+1)**2
          if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm
          sumbox=0.d0
          do jj=1,max_jMps_l(iatom)
            sumbox=sumbox+uV2nd(jj,lm,iatom)*tpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iob,1)
          end do
          uVpsibox4(lm,iatom,iob,1)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm
      end do
!$OMP end do nowait
    else
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
                            Jxyz(3,jMps_l(jj,iatom),iatom),iob,1)
            end do
          else if(iwk_size>=11.and.iwk_size<=12)then
            do jj=1,max_jMps_l_s(iatom)
              sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                       tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                            Jxyz(3,jMps_l_s(jj,iatom),iatom),iob,1)
            end do
          end if
          uVpsibox3(lm,iatom,iob,1)=sumbox*Hvol/uVu(lm,iatom)
        end do loop_lm2
      end do
!$OMP end do nowait
    end if
!$OMP end parallel
  end do

  if(nproc_Mxin_mul/=1)then
    if(iwk_size>=1.and.iwk_size<=2)then
      elp3(705)=get_wtime()
      call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*iobmax,nproc_group_orbital)
      elp3(706)=get_wtime()
      elp3(744)=elp3(744)+elp3(706)-elp3(705)
    else if(iwk_size>=11.and.iwk_size<=12)then
      call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*iobmax,nproc_group_h)
    end if
  end if

end if

! Pseudopotential 2 (non-local)
if(iflag_ps==1) then
  if(ikind_eext==0.and.icalcforce==0.and.iflag_md==0)then
!$OMP parallel
    do iob=1,iobmax
!$OMP do private(jj,iatom,ja,iix,iiy,iiz,ikoa,lm)
      do jj=1,maxMps_all
        iix=Jxyz_all(1,jj)
        iiy=Jxyz_all(2,jj)
        iiz=Jxyz_all(3,jj)
        do iatom=1,numatom_ps_2nd(jj)
          ja=jja(iatom,jj)
          ikoa=Kion(iatomnum_ps_2nd(ja))
          do lm=1,(Mlps(ikoa)+1)**2
          htpsi(iix,iiy,iiz,iob,1)= &
            htpsi(iix,iiy,iiz,iob,1) + &
            uVpsibox4(lm,iatomnum_ps_2nd(ja),iob,1)*uV(Mps3rd_2nd(ja),lm,iatomnum_ps_2nd(ja))
          end do
        end do
      end do
!$OMP end do nowait
    end do
!$OMP end parallel
  else
    do iob=1,iobmax
      do iatom=1,MI
      ikoa=Kion(iatom)
        do jj=1,max_jMps_l(iatom)
          do lm=1,(Mlps(ikoa)+1)**2
            htpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iob,1)= &
              htpsi(Jxyz2nd(1,jj,iatom),Jxyz2nd(2,jj,iatom),Jxyz2nd(3,jj,iatom),iob,1) + &
                uVpsibox4(lm,iatom,iob,1)*uV2nd(jj,lm,iatom)
          end do
        end do
      end do
    end do
  end if
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
