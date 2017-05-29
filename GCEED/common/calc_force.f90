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
subroutine calc_force
use scf_data
use allocate_mat_sub
use read_pslfile_sub
use new_world_sub
implicit none
integer :: ix,iy,iz,iob,ikoa,jj,j2,iatom,ia,ib,lm,ikoa2
real(8) :: rbox1,rbox2
real(8),allocatable :: uVpsibox(:,:,:,:),uVpsibox2(:,:,:,:)
real(8) :: rforce1(3,MI),rforce2(3,MI),rforce3(3,MI)
real(8) :: rab
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)

do iatom=1,MI
do j2=1,3
  rforce(j2,iatom)=0.d0
  rforce1(j2,iatom)=0.d0
  rforce2(j2,iatom)=0.d0
  rforce3(j2,iatom)=0.d0
end do
end do

! ion-ion
do ia=1,MI
  ikoa=Kion(ia)
  do ib=1,MI
    if(ia/=ib)then
      ikoa2=Kion(ib)
      rab=sqrt((Rion(1,ia)-Rion(1,ib))**2+   &
               (Rion(2,ia)-Rion(2,ib))**2+   &
               (Rion(3,ia)-Rion(3,ib))**2)
      do j2=1,3
        rforce(j2,ia)=rforce(j2,ia)+Zps(ikoa)*Zps(ikoa2)*(Rion(j2,ia)-Rion(j2,ib))/rab**3
        rforce1(j2,ia)=rforce1(j2,ia)+Zps(ikoa)*Zps(ikoa2)*(Rion(j2,ia)-Rion(j2,ib))/rab**3
      end do
    end if
  end do
end do

tpsi=0.d0
do iob=1,iobnum
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)
  end do
  end do
  end do
end do

call calc_gradient_fast(tpsi,rgrad_wk)

! local part of force
do iatom=1,MI
do j2=1,3
  rbox1=0.d0
  do iob=1,iobnum
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rbox1=rbox1-2.d0*rocc(iob,1)*rgrad_wk(ix,iy,iz,iob,1,j2)*Vpsl_atom(ix,iy,iz,iatom)*psi(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
  call MPI_Allreduce(rbox1,rbox2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rforce(j2,iatom)=rforce(j2,iatom)+rbox2*Hvol
  rforce2(j2,iatom)=rbox2*Hvol
end do
end do

! nonlocal part of force

allocate (uVpsibox(1:iobnum,1,1:maxlm,1:MI))
allocate (uVpsibox2(1:iobnum,1,1:maxlm,1:MI))

do iatom=1,MI
  do lm=1,maxlm
    do iob=1,iobnum
      uVpsibox(iob,1,lm,iatom)=0.d0
    end do
  end do
end do

do iatom=1,MI
  ikoa=Kion(iatom)
  do iob=1,iobnum
    loop_lm2 : do lm=1,(Mlps(ikoa)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm2
      rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 )
      do jj=1,max_jMps_l(iatom)
        rbox1=rbox1+uV(jMps_l(jj,iatom),lm,iatom)*  &
                       psi(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),  &
                           Jxyz(3,jMps_l(jj,iatom),iatom),iob,1)
      end do
      uVpsibox(iob,1,lm,iatom)=rbox1*Hvol/uVu(lm,iatom)
    end do loop_lm2
  end do
end do

call MPI_allreduce(uVpsibox,uVpsibox2,iobnum*maxlm*MI,      &
                     MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_orbital,ierr)

do iatom=1,MI
  ikoa=Kion(iatom)
  do j2=1,3
    rbox1=0.d0
    do iob=1,iobnum
      do jj=1,max_jMps_l(iatom)
        do lm=1,(Mlps(ikoa)+1)**2
          rbox1=rbox1-2.d0*rocc(iob,1)*uV(jMps_l(jj,iatom),lm,iatom)*   &
                  rgrad_wk(Jxyz(1,jMps_l(jj,iatom),iatom),Jxyz(2,jMps_l(jj,iatom),iatom),Jxyz(3,jMps_l(jj,iatom),iatom),iob,1,j2)* &
                  uVpsibox2(iob,1,lm,iatom)
        end do
      end do
    end do
    call MPI_allreduce(rbox1,rbox2,1,      &
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    rforce(j2,iatom)=rforce(j2,iatom)+rbox2*Hvol
    rforce3(j2,iatom)=rbox2*Hvol
  end do
end do

if(myrank==0)then
  write(*,*) "===== force ====="
  do iatom=1,MI
    write(*,'(i6,3e16.8)') iatom,(rforce(j2,iatom)*2.d0*Ry/a_B,j2=1,3)
  end do
end if

deallocate(uVpsibox,uVpsibox2)

end subroutine calc_force
