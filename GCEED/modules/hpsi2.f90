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

   MODULE PROCEDURE R_hpsi2,C_hpsi2

END INTERFACE

CONTAINS

!=======================================================================
!========================== Hamiltonian Operation ( for real functions )

SUBROUTINE R_hpsi2(tpsi,htpsi,iob,nn,isub)
!$ use omp_lib
implicit none
real(8) :: tpsi(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
real(8) :: htpsi(iwk3sta(1):iwk3end(1),  &
                               iwk3sta(2):iwk3end(2),      &
                               iwk3sta(3):iwk3end(3))

integer :: iatom,ix,iy,iz,jj,lm,ik
integer :: iob,nn,isub

real(8) :: sumbox

real(8), allocatable :: uVpsibox(:,:)
real(8), allocatable :: uVpsibox2(:,:)
real(8) :: rlap_wk(iwk3sta(1):iwk3end(1),  &
                  iwk3sta(2):iwk3end(2),      &
                  iwk3sta(3):iwk3end(3))
real(8) :: wk2(iwk2sta(1):iwk2end(1),  &
               iwk2sta(2):iwk2end(2),      &
               iwk2sta(3):iwk2end(3))

integer :: ispin

call set_ispin(iob,ispin)

if(iflag_ps==1) then
  allocate (uVpsibox(1:maxlm,1:MI))
  allocate (uVpsibox2(1:maxlm,1:MI))
end if

! Pseudopotential 1 (non-local)
if(iflag_ps.eq.1)then
  if(nproc_Mxin_mul==1)then
!$OMP parallel do
    do iatom=1,MI
      do lm=1,maxlm
        uVpsibox2(lm,iatom)=0.d0
      end do
    end do
  else
!$OMP parallel do
    do iatom=1,MI
      do lm=1,maxlm
        uVpsibox(lm,iatom)=0.d0
      end do
    end do
  end if

  if(nproc_Mxin_mul==1)then
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
                          Jxyz(3,jMps_l(jj,iatom),iatom))
          end do
        else if(iwk_size>=11.and.iwk_size<=12)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l_s(iatom)
            sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l_s(jj,iatom),iatom))
          end do
        end if
        uVpsibox2(lm,iatom)=sumbox*Hvol/uVu(lm,iatom)
      end do loop_lm
    end do
  else
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
                          Jxyz(3,jMps_l(jj,iatom),iatom))
          end do
        else if(iwk_size>=11.and.iwk_size<=12)then
!$OMP parallel do reduction( + : sumbox )
          do jj=1,max_jMps_l_s(iatom)
            sumbox=sumbox+uV(jMps_l_s(jj,iatom),lm,iatom)*  &
                     tpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),  &
                          Jxyz(3,jMps_l_s(jj,iatom),iatom))
          end do
        end if
        uVpsibox(lm,iatom)=sumbox*Hvol/uVu(lm,iatom)
      end do loop_lm2
    end do
  end if

  if(nproc_Mxin_mul==1)then
    continue
  else
    if(iwk_size>=1.and.iwk_size<=2)then
      call MPI_allreduce(uVpsibox(1,1),uVpsibox2(1,1),maxlm*MI,      &
                     MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_orbital,ierr)
    else if(iwk_size>=11.and.iwk_size<=12)then
      call MPI_allreduce(uVpsibox(1,1),uVpsibox2(1,1),maxlm*MI,      &
                     MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)
    end if
  end if

end if

! Kinetic energy

! Pseudo(local), Hartree, Exchange-Correlation potential
if(ihpsieff==0)then
!$OMP parallel do
  do iz=iwk3sta(3),iwk3end(3)
  do iy=iwk3sta(2),iwk3end(2)
  do ix=iwk3sta(1),iwk3end(1)
    htpsi(ix,iy,iz) = Vlocal(ix,iy,iz,ispin) *tpsi(ix,iy,iz) 
  end do
  end do
  end do
else if(ihpsieff==1)then
!$OMP parallel do
  do iz=iwk3sta(3),iwk3end(3)
  do iy=iwk3sta(2),iwk3end(2)
  do ix=iwk3sta(1),iwk3end(1)
    htpsi(ix,iy,iz) = (Vlocal(ix,iy,iz,ispin)+Vbox(ix,iy,iz)) *tpsi(ix,iy,iz) 
  end do
  end do
  end do
end if

! Pseudopotential 2 (non-local)
if(iflag_ps==1) then
  if(iwk_size>=1.and.iwk_size<=2)then
    do iatom=1,MI
      ik=Kion(iatom)
      do jj=1,max_jMps_l(iatom)
        do lm=1,(Mlps(ik)+1)**2
          htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom))= &
            htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom)) + &
            uVpsibox2(lm,iatom)*uV(jMps_l(jj,iatom),lm,iatom)
        end do
      end do
    end do
  else if(iwk_size>=11.and.iwk_size<=12)then
    do iatom=1,MI
      ik=Kion(iatom)
      do jj=1,max_jMps_l_s(iatom)
        do lm=1,(Mlps(ik)+1)**2
          htpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),Jxyz(3,jMps_l_s(jj,iatom),iatom))= &
            htpsi(Jxyz(1,jMps_l_s(jj,iatom),iatom),Jxyz(2,jMps_l_s(jj,iatom),iatom),Jxyz(3,jMps_l_s(jj,iatom),iatom)) + &
            uVpsibox2(lm,iatom)*uV(jMps_l_s(jj,iatom),lm,iatom)
        end do
      end do
    end do
  end if
end if

wk2=0.d0
!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  wk2(ix,iy,iz)=tpsi(ix,iy,iz)
end do
end do
end do
call sendrecv(wk2)
call calc_laplacian2(wk2,rlap_wk)

!$OMP parallel do              
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  htpsi(ix,iy,iz)=htpsi(ix,iy,iz)-0.5d0*rlap_wk(ix,iy,iz)
end do                         
end do
end do

if(iflag_ps==1) deallocate(uVpsibox,uVpsibox2)

return

END SUBROUTINE R_hpsi2

!=======================================================================
!========================= Hamiltonian Operation (for complex funcitons)

SUBROUTINE C_hpsi2(tpsi,htpsi,iob,nn,isub)
!$ use omp_lib
implicit none
complex(8) :: tpsi(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
complex(8) :: htpsi(iwk3sta(1):iwk3end(1),  &
                    iwk3sta(2):iwk3end(2),      &
                    iwk3sta(3):iwk3end(3))

integer :: iatom,ix,iy,iz,jj,lm,ik
integer :: iob,nn,isub

complex(8) :: sumbox

complex(8) :: clap_wk( iwk3sta(1):iwk3end(1),  &
                      iwk3sta(2):iwk3end(2),      &
                      iwk3sta(3):iwk3end(3))

complex(8) :: grad_wk( 3,iwk3sta(1):iwk3end(1),  &   
                         iwk3sta(2):iwk3end(2),      &
                         iwk3sta(3):iwk3end(3))
complex(8), allocatable :: uVpsibox(:,:)
complex(8), allocatable :: uVpsibox2(:,:)

complex(8), parameter :: zi=(0.d0,1.d0)

real(8) :: f0

integer :: ispin

real(8) :: fdN0
real(8) :: fdN1(0:12,3)

integer :: j,ind

call set_ispin(iob,ispin)

if(nn<=1)then
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    htpsi(ix,iy,iz) = 0.d0
  end do
  end do
  end do
end if

elp3(751)=MPI_Wtime()

! Pseudopotential 1 (non-local)

if(iflag_ps==1) then
  allocate (uVpsibox(1:maxlm,1:MI))
  allocate (uVpsibox2(1:maxlm,1:MI))
end if

if(iflag_ps.eq.1)then
!$OMP parallel do
  do iatom=1,MI
    do lm=1,maxlm
      uVpsibox(lm,iatom)=0.d0
    end do
  end do

  do iatom=1,MI
    ik=Kion(iatom)
    loop_lm2 : do lm=1,(Mlps(ik)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
      sumbox=0.d0
!$OMP parallel do reduction( + : sumbox )
      do jj=1,max_jMps_l(iatom)
        sumbox=sumbox+uV(jMps_l(jj,iatom),lm,iatom)*  &
                 tpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),  &
                      Jxyz(3,jMps_l(jj,iatom),iatom))
      end do
      uVpsibox(lm,iatom)=sumbox*Hvol/uVu(lm,iatom)
    end do loop_lm2
  end do

  call MPI_allreduce(uVpsibox(1,1),uVpsibox2(1,1),maxlm*MI,      &
                 MPI_DOUBLE_COMPLEX,MPI_SUM,newworld_comm_orbital,ierr)
end if

! Pseudopotential 2 (non-local)
if(iflag_ps==1) then
  do iatom=1,MI
    ik=Kion(iatom)
    do jj=1,max_jMps_l(iatom)
      do lm=1,(Mlps(ik)+1)**2
        htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom))= &
          htpsi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom)) + &
          uVpsibox2(lm,iatom)*uV(jMps_l(jj,iatom),lm,iatom)
      end do
    end do
  end do
end if
 
elp3(755)=MPI_Wtime()
elp3(775)=elp3(775)+elp3(755)-elp3(754)
elp3(785)=elp3(785)+elp3(755)-elp3(751)

! Kinetic energy


if(isub==0)then

  if(nproc_Mxin_mul>1) call sendrecv(tpsi)

  f0=(1.d0/Hgs(1)**2   &
     +1.d0/Hgs(2)**2   &
     +1.d0/Hgs(3)**2)
  fdN0=-0.5d0*cNmat(0,Nd)*f0
  do j=1,3
    do ind=1,4
      fdN1(ind,j)=-0.5d0*cNmat(ind,Nd)/Hgs(j)**2
    end do
  end do

  elp3(757)=MPI_Wtime()
!$OMP parallel do
  do iz=iwk3sta(3),iwk3end(3)
  do iy=iwk3sta(2),iwk3end(2)
  do ix=iwk3sta(1),iwk3end(1)
    htpsi(ix,iy,iz)=(Vlocal(ix,iy,iz,ispin)+fdN0)*tpsi(ix,iy,iz)     &
      +fdN1(1,1)*(tpsi(ix+1,iy,iz) + tpsi(ix-1,iy,iz))      &
      +fdN1(1,2)*(tpsi(ix,iy+1,iz) + tpsi(ix,iy-1,iz))      &
      +fdN1(1,3)*(tpsi(ix,iy,iz+1) + tpsi(ix,iy,iz-1))      &
      +fdN1(2,1)*(tpsi(ix+2,iy,iz) + tpsi(ix-2,iy,iz))      &
      +fdN1(2,2)*(tpsi(ix,iy+2,iz) + tpsi(ix,iy-2,iz))      &
      +fdN1(2,3)*(tpsi(ix,iy,iz+2) + tpsi(ix,iy,iz-2))      &
      +fdN1(3,1)*(tpsi(ix+3,iy,iz) + tpsi(ix-3,iy,iz))      &
      +fdN1(3,2)*(tpsi(ix,iy+3,iz) + tpsi(ix,iy-3,iz))      &
      +fdN1(3,3)*(tpsi(ix,iy,iz+3) + tpsi(ix,iy,iz-3))      &
      +fdN1(4,1)*(tpsi(ix+4,iy,iz) + tpsi(ix-4,iy,iz))      &
      +fdN1(4,2)*(tpsi(ix,iy+4,iz) + tpsi(ix,iy-4,iz))      &
      +fdN1(4,3)*(tpsi(ix,iy,iz+4) + tpsi(ix,iy,iz-4)) 
  end do
  end do
  end do

  elp3(758)=MPI_Wtime()
  elp3(778)=elp3(778)+elp3(758)-elp3(757)

else if(isub==1)then
  if(nproc_Mxin_mul>1) call sendrecv(tpsi)
  elp3(757)=MPI_Wtime()
  call calc_laplacian2(tpsi,clap_wk)
  call calc_gradient2(tpsi,grad_wk)

  if(nn==N_hamil)then
!$OMP parallel do
    do iz=iwk3sta(3),iwk3end(3)
    do iy=iwk3sta(2),iwk3end(2)
    do ix=iwk3sta(1),iwk3end(1)
!ocl prefetch_write(htpsi(ix+16,iy,iz))
!ocl prefetch_write(zpsi(ix+16,iy,iz))
      htpsi(ix,iy,iz) = htpsi(ix,iy,iz) + Vlocal(ix,iy,iz,ispin) *tpsi(ix,iy,iz) -1d0/2d0*clap_wk(ix,iy,iz)
      zpsi(ix,iy,iz,iob,1)=zpsi(ix,iy,iz,iob,1)+zc(nn)*htpsi(ix,iy,iz)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+abs(zpsi(ix,iy,iz,iob,1))**2*rocc(iob,1)*wtk(1)
    end do
    end do
    end do
  else
!$OMP parallel do
    do iz=iwk3sta(3),iwk3end(3)
    do iy=iwk3sta(2),iwk3end(2)
    do ix=iwk3sta(1),iwk3end(1)
!ocl prefetch_write(htpsi(ix+16,iy,iz))
!ocl prefetch_write(zpsi(ix+16,iy,iz))
      htpsi(ix,iy,iz) = htpsi(ix,iy,iz) + Vlocal(ix,iy,iz,ispin) *tpsi(ix,iy,iz) -1d0/2d0*clap_wk(ix,iy,iz)
      zpsi(ix,iy,iz,iob,1)=zpsi(ix,iy,iz,iob,1)+zc(nn)*htpsi(ix,iy,iz)
      tpsi(ix,iy,iz)=htpsi(ix,iy,iz)
      htpsi(ix,iy,iz)=0.d0
    end do
    end do
    end do
  end if
  elp3(758)=MPI_Wtime()
  elp3(778)=elp3(778)+elp3(758)-elp3(757)
end if


elp3(759)=MPI_Wtime()
elp3(780)=elp3(780)+elp3(759)-elp3(751)

return

END SUBROUTINE C_hpsi2

END MODULE hpsi2_sub

