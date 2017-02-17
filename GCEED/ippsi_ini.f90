! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.


MODULE ippsi_ini_sub

use scf_data
use new_world_sub
implicit none
INTERFACE ippsi_ini

   MODULE PROCEDURE R_ippsi_ini,C_ippsi_ini

END INTERFACE

CONTAINS

!=============================================================================================
subroutine R_ippsi_ini(matbox,ista_ini,iend_ini,file_OUT_data_ini)
real(8) :: box
real(8) :: matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))
real(8),allocatable :: matbox3(:,:,:)
real(8) :: sum0,sum1
real(8) :: epsilon=1.d-10
integer :: ista_ini(3),iend_ini(3)
real(8) :: coo_bound_sta(3),coo_bound_end(3),coo_ini(3)
real(8) :: t_ini(3)
integer :: ibox(3)
integer :: inmol,i1,i2,i3
integer :: ix,iy,iz
character(100) :: file_OUT_data_ini

allocate( matbox3(lg_sta_ini(1):lg_end_ini(1),lg_sta_ini(2):lg_end_ini(2),lg_sta_ini(3):lg_end_ini(3)) )
  
do inmol=1,num_mol
  if(imesh_oddeven==1)then
    coo_bound_sta(:)=Hgs(:)*lg_sta(:)+coo_mol_ini(:,inmol)
    coo_bound_end(:)=Hgs(:)*(lg_end(:)-1)+coo_mol_ini(:,inmol)
  else if(imesh_oddeven==2)then
    coo_bound_sta(:)=Hgs(:)*lg_sta(:)+coo_mol_ini(:,inmol)-0.5d0*H
    coo_bound_end(:)=Hgs(:)*(lg_end(:)-1)+coo_mol_ini(:,inmol)-0.5d0*H
  end if
  matbox3=0.d0
  sum0=0.d0
  if(((num_datafiles_OUT==1.or.num_datafiles_OUT>nproc).and.myrank==0).or.   &
     ((num_datafiles_OUT>1.and.num_datafiles_OUT<=nproc).and.myrank<num_datafiles_OUT))then
    do i1=ista_ini(1),iend_ini(1)
    do i2=ista_ini(2),iend_ini(2)
    do i3=ista_ini(3),iend_ini(3)
      if(imesh_oddeven==1)then
        coo_ini(1)=H_ini(1)*i1
        coo_ini(2)=H_ini(2)*i2
        coo_ini(3)=H_ini(3)*i3
      else if(imesh_oddeven==2)then
        coo_ini(1)=H_ini(1)*i1-0.5d0*H_ini(1)
        coo_ini(2)=H_ini(2)*i2-0.5d0*H_ini(2)
        coo_ini(3)=H_ini(3)*i3-0.5d0*H_ini(3)
      end if
      if(coo_ini(1)>coo_bound_sta(1)+epsilon.and.  &
         coo_ini(1)<coo_bound_end(1)-epsilon.and.  &
         coo_ini(2)>coo_bound_sta(2)+epsilon.and.  &
         coo_ini(2)<coo_bound_end(2)-epsilon.and.  &
         coo_ini(3)>coo_bound_sta(3)+epsilon.and.  &
         coo_ini(3)<coo_bound_end(3)-epsilon) then
        ibox(:)=int((coo_ini(:)-coo_bound_sta(:))/Hgs(:))+lg_sta(:)
        t_ini(:)=(coo_ini(:)-coo_bound_sta(:)-dble(ibox(:)-lg_sta(:))*Hgs(:))/Hgs(:)
        matbox3(i1,i2,i3)=(1.d0-t_ini(1))*(1.d0-t_ini(2))*(1.d0-t_ini(3))*matbox(ibox(1),ibox(2),ibox(3)) &
           +(1.d0-t_ini(1))*(1.d0-t_ini(2))*t_ini(3)*matbox(ibox(1),ibox(2),ibox(3)+1) &
           +(1.d0-t_ini(1))*t_ini(2)*(1.d0-t_ini(3))*matbox(ibox(1),ibox(2)+1,ibox(3)) &
           +(1.d0-t_ini(1))*t_ini(2)*t_ini(3)*matbox(ibox(1),ibox(2)+1,ibox(3)+1) &
           +t_ini(1)*(1.d0-t_ini(2))*(1.d0-t_ini(3))*matbox(ibox(1)+1,ibox(2),ibox(3)) &
           +t_ini(1)*(1.d0-t_ini(2))*t_ini(3)*matbox(ibox(1)+1,ibox(2),ibox(3)+1) &
           +t_ini(1)*t_ini(2)*(1.d0-t_ini(3))*matbox(ibox(1)+1,ibox(2)+1,ibox(3)) &
           +t_ini(1)*t_ini(2)*t_ini(3)*matbox(ibox(1)+1,ibox(2)+1,ibox(3)+1)
        sum0=sum0+matbox3(i1,i2,i3)**2
      else
        matbox3(i1,i2,i3)=0.d0
      end if
    end do
    end do
    end do
  end if
  call MPI_Allreduce(sum0,sum1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  sum1=sqrt(sum1*H_ini(1)*H_ini(2)*H_ini(3))
  if(((num_datafiles_OUT==1.or.num_datafiles_OUT>nproc).and.myrank==0).or.   &
     ((num_datafiles_OUT>1.and.num_datafiles_OUT<=nproc).and.myrank<num_datafiles_OUT))then
    write(67) ((( matbox3(ix,iy,iz)/sum1,ix=ista_ini(1),iend_ini(1)),iy=ista_ini(2),iend_ini(2)),iz=ista_ini(3),iend_ini(3))
  end if
end do

deallocate(matbox3)

end subroutine R_ippsi_ini

!=============================================================================================
subroutine C_ippsi_ini(matbox,ista_ini,iend_ini,file_OUT_data_ini)
complex(8) :: box
complex(8) :: matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))
complex(8),allocatable :: matbox3(:,:,:)
real(8) :: sum0,sum1
real(8) :: epsilon=1.d-10
integer :: ista_ini(3),iend_ini(3)
real(8) :: coo_bound_sta(3),coo_bound_end(3),coo_ini(3)
real(8) :: t_ini(3)
integer :: ibox(3)
integer :: inmol,i1,i2,i3
integer :: ix,iy,iz
character(100) :: file_OUT_data_ini

allocate( matbox3(lg_sta_ini(1):lg_end_ini(1),lg_sta_ini(2):lg_end_ini(2),lg_sta_ini(3):lg_end_ini(3)) )

do inmol=1,num_mol
  if(imesh_oddeven==1)then
    coo_bound_sta(:)=Hgs(:)*lg_sta(:)+coo_mol_ini(:,inmol)
    coo_bound_end(:)=Hgs(:)*(lg_end(:)-1)+coo_mol_ini(:,inmol)
  else if(imesh_oddeven==2)then
    coo_bound_sta(:)=Hgs(:)*lg_sta(:)+coo_mol_ini(:,inmol)-0.5d0*Hgs(:)
    coo_bound_end(:)=Hgs(:)*(lg_end(:)-1)+coo_mol_ini(:,inmol)-0.5d0*Hgs(:)
  end if
  matbox3=0.d0
  sum0=0.d0
  if(((num_datafiles_OUT==1.or.num_datafiles_OUT>nproc).and.myrank==0).or.   &
     ((num_datafiles_OUT>1.and.num_datafiles_OUT<=nproc).and.myrank<num_datafiles_OUT))then
    do i1=ista_ini(1),iend_ini(1)
    do i2=ista_ini(2),iend_ini(2)
    do i3=ista_ini(3),iend_ini(3)
      if(imesh_oddeven==1)then
        coo_ini(1)=H_ini(1)*i1
        coo_ini(2)=H_ini(2)*i2
        coo_ini(3)=H_ini(3)*i3
      else if(imesh_oddeven==2)then
        coo_ini(1)=H_ini(1)*i1-0.5d0*H_ini(1)
        coo_ini(2)=H_ini(2)*i2-0.5d0*H_ini(2)
        coo_ini(3)=H_ini(3)*i3-0.5d0*H_ini(3)
      end if
      if(coo_ini(1)>coo_bound_sta(1)+epsilon.and.  &
         coo_ini(1)<coo_bound_end(1)-epsilon.and.  &
         coo_ini(2)>coo_bound_sta(2)+epsilon.and.  &
         coo_ini(2)<coo_bound_end(2)-epsilon.and.  &
         coo_ini(3)>coo_bound_sta(3)+epsilon.and.  &
         coo_ini(3)<coo_bound_end(3)-epsilon) then
        ibox(:)=int((coo_ini(:)-coo_bound_sta(:))/Hgs(:))+lg_sta(:)
        t_ini(:)=(coo_ini(:)-coo_bound_sta(:)-dble(ibox(:)-lg_sta(:))*Hgs(:))/Hgs(:)
        matbox3(i1,i2,i3)=(1.d0-t_ini(1))*(1.d0-t_ini(2))*(1.d0-t_ini(3))*matbox(ibox(1),ibox(2),ibox(3)) &
           +(1.d0-t_ini(1))*(1.d0-t_ini(2))*t_ini(3)*matbox(ibox(1),ibox(2),ibox(3)+1) &
           +(1.d0-t_ini(1))*t_ini(2)*(1.d0-t_ini(3))*matbox(ibox(1),ibox(2)+1,ibox(3)) &
           +(1.d0-t_ini(1))*t_ini(2)*t_ini(3)*matbox(ibox(1),ibox(2)+1,ibox(3)+1) &
           +t_ini(1)*(1.d0-t_ini(2))*(1.d0-t_ini(3))*matbox(ibox(1)+1,ibox(2),ibox(3)) &
           +t_ini(1)*(1.d0-t_ini(2))*t_ini(3)*matbox(ibox(1)+1,ibox(2),ibox(3)+1) &
           +t_ini(1)*t_ini(2)*(1.d0-t_ini(3))*matbox(ibox(1)+1,ibox(2)+1,ibox(3)) &
           +t_ini(1)*t_ini(2)*t_ini(3)*matbox(ibox(1)+1,ibox(2)+1,ibox(3)+1)
        sum0=sum0+abs(matbox3(i1,i2,i3))**2
      else
        matbox3(i1,i2,i3)=0.d0
      end if
    end do
    end do
    end do
  end if
  call MPI_Allreduce(sum0,sum1,1,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)
  sum1=sqrt(sum1*H_ini(1)*H_ini(2)*H_ini(3))
  if(((num_datafiles_OUT==1.or.num_datafiles_OUT>nproc).and.myrank==0).or.   &
     ((num_datafiles_OUT>1.and.num_datafiles_OUT<=nproc).and.myrank<num_datafiles_OUT))then
    write(67) ((( matbox3(ix,iy,iz)/sum1,ix=ista_ini(1),iend_ini(1)),iy=ista_ini(2),iend_ini(2)),iz=ista_ini(3),iend_ini(3))
  end if
end do

deallocate(matbox3)

end subroutine C_ippsi_ini

END MODULE ippsi_ini_sub
