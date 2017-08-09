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
!======================================================================
!======================================================================
subroutine writepsi
  use inputoutput, only: au_length_aa
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use scf_data
  use allocate_mat_sub
  implicit none
  integer :: iob,ix,iy,iz
  integer :: p0,iob_myob,icheck_corrkob
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit
 
  if(iSCFRT==1)then
    do p0=1,itotMST
      call conv_p0(p0,iob)
      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,icheck_corrkob)
  !OMP parallel do private(iz,iy,ix)
      do iz=lg_sta(3),lg_end(3)
      do iy=lg_sta(2),lg_end(2)
      do ix=lg_sta(1),lg_end(1)
        matbox_l(ix,iy,iz)=0.d0
      end do
      end do
      end do
      if(icheck_corrkob==1)then
  !OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_l(ix,iy,iz)=psi(ix,iy,iz,iob_myob,1)
        end do
        end do
        end do
      end if

      if(format3d=='avs')then
  !OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_l(ix,iy,iz)=matbox_l(ix,iy,iz)/sqrt(au_length_aa)**3
        end do
        end do
        end do
      end if

      call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

      write(filenum, '(i5)') p0
      suffix = "psi"//trim(adjustl(filenum))
      phys_quantity = "psi"
      if(format3d=='avs')then
        header_unit = "A**(-3/2)"
        call writeavs(103,suffix,header_unit,matbox_l2)
      else if(format3d=='cube')then
        call writecube(103,suffix,phys_quantity,matbox_l2)
      end if
    end do
  end if
  
end subroutine writepsi

