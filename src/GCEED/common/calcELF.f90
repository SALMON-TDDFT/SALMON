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
subroutine calcELF(filename_ELF,ttmp)
use salmon_parallel, only: nproc_group_global,nproc_group_kgrid
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use gradient_sub
use allocate_mat_sub
use new_world_sub
implicit none

integer :: iob,ix,iy,iz
integer :: p_allob
integer :: ttmp

real(8) :: elftau(mg_sta(1):mg_end(1),   &
                  mg_sta(2):mg_end(2),   &
                  mg_sta(3):mg_end(3))
real(8) :: mrelftau(mg_sta(1):mg_end(1),   &
                    mg_sta(2):mg_end(2),   &
                    mg_sta(3):mg_end(3))
real(8) :: curden(mg_sta(1):mg_end(1),   &
                  mg_sta(2):mg_end(2),   &
                  mg_sta(3):mg_end(3))
real(8) :: mrcurden(mg_sta(1):mg_end(1),   &
                    mg_sta(2):mg_end(2),   &
                    mg_sta(3):mg_end(3))
real(8) :: gradpsi(3,mg_sta(1):mg_end(1),   &
                     mg_sta(2):mg_end(2),   &
                     mg_sta(3):mg_end(3))
complex(8) :: gradzpsi(3,mg_sta(1):mg_end(1),   &
                         mg_sta(2):mg_end(2),   &
                         mg_sta(3):mg_end(3))
real(8) :: gradrho(3,mg_sta(1):mg_end(1),   &
                     mg_sta(2):mg_end(2),   &
                     mg_sta(3):mg_end(3))
real(8) :: gradrho2(mg_sta(1):mg_end(1),   &
                    mg_sta(2):mg_end(2),   &
                    mg_sta(3):mg_end(3))
real(8) :: elfc(mg_sta(1):mg_end(1),   &
                mg_sta(2):mg_end(2),   &
                mg_sta(3):mg_end(3))
real(8) :: elfcuni(mg_sta(1):mg_end(1),   &
                   mg_sta(2):mg_end(2),   &
                   mg_sta(3):mg_end(3))
real(8) :: rho_half(mg_sta(1):mg_end(1),   &
                    mg_sta(2):mg_end(2),   &
                    mg_sta(3):mg_end(3))
character(50) :: filename_ELF

elp3(801)=get_wtime()

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho_half(ix,iy,iz)=rho(ix,iy,iz)/2.d0
end do
end do
end do
mrelftau=0.d0
mrcurden=0.d0

elp3(802)=get_wtime()
elp3(832)=elp3(832)+elp3(802)-elp3(801)

iwk_size=1
call make_iwksta_iwkend

if(iSCFRT==1)then
  if(iperiodic==0)then

    do iob=1,iobnum
      call calc_allob(iob,p_allob)
      if((ilsda==0.and.p_allob<=ifMST(1)).or.  &
         (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then
  
        call calc_gradient(psi(:,:,:,iob,1),gradpsi(:,:,:,:))
  
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradpsi(1,ix,iy,iz))**2      &
                             +abs(gradpsi(2,ix,iy,iz))**2      &
                             +abs(gradpsi(3,ix,iy,iz))**2
        end do
        end do
        end do
      end if
    end do

  else if(iperiodic==3)then

    do iob=1,iobnum
      call calc_allob(iob,p_allob)
      if((ilsda==0.and.p_allob<=ifMST(1)).or.  &
         (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then
  
        call calc_gradient(zpsi(:,:,:,iob,1),gradzpsi(:,:,:,:))
  
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                             +abs(gradzpsi(2,ix,iy,iz))**2      &
                             +abs(gradzpsi(3,ix,iy,iz))**2
        end do
        end do
        end do
      end if
    end do

  end if

  call comm_summation(mrelftau,elftau,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_kgrid)


  call calc_gradient(rho_half(:,:,:),gradrho(:,:,:,:))
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
          +gradrho(2,ix,iy,iz)**2      &
          +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0
  end do
  end do
  end do

else 

  do iob=1,iobnum

    call calc_allob(iob,p_allob)
    if((ilsda==0.and.p_allob<=ifMST(1)).or.   &
       (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then

      cmatbox_m=0.d0
      if(iSCFRT==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox_m(ix,iy,iz)=zpsi(ix,iy,iz,iob,1)
        end do
        end do
        end do
      else if(iSCFRT==2)then
        if(mod(ttmp,2)==1)then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
            cmatbox_m(ix,iy,iz)=zpsi_out(ix,iy,iz,iob,1)
          end do
          end do
          end do
        else if(mod(ttmp,2)==0)then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
            cmatbox_m(ix,iy,iz)=zpsi_in(ix,iy,iz,iob,1)
          end do
          end do
          end do
        end if
      end if
 
      call calc_gradient(cmatbox_m(:,:,:),gradzpsi(:,:,:,:))

      elp3(807)=get_wtime()

!$OMP parallel do private(iz,iy,ix)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
  
        mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                           +abs(gradzpsi(2,ix,iy,iz))**2      &
                           +abs(gradzpsi(3,ix,iy,iz))**2
  
        mrcurden(ix,iy,iz)=mrcurden(ix,iy,iz)      &
             +( abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(1,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(1,ix,iy,iz)))**2      &
               +abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(2,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(2,ix,iy,iz)))**2      &
               +abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(3,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(3,ix,iy,iz)))**2 )/2.d0
  
      end do
      end do
      end do
      
      elp3(808)=get_wtime()
      elp3(838)=elp3(838)+elp3(808)-elp3(807)

    end if
  end do

  elp3(809)=get_wtime()
  elp3(839)=elp3(839)+elp3(809)-elp3(808)
  
  call comm_summation(mrelftau,elftau,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_kgrid)
  call comm_summation(mrcurden,curden,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_kgrid)

  
  elp3(810)=get_wtime()
  elp3(840)=elp3(840)+elp3(810)-elp3(809)
  
  call calc_gradient(rho_half(:,:,:),gradrho(:,:,:,:))

  elp3(815)=get_wtime()

  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
          +gradrho(2,ix,iy,iz)**2      &
          +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0  &
                                   -curden(ix,iy,iz)/rho_half(ix,iy,iz)
  end do
  end do
  end do

elp3(816)=get_wtime()
elp3(846)=elp3(846)+elp3(816)-elp3(815)

end if

! matbox_l stores ELF
matbox_l=0.d0
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  elfcuni(ix,iy,iz)=3.d0/5.d0*(6.d0*Pi**2)**(2.d0/3.d0)      &
            *rho_half(ix,iy,iz)**(5.d0/3.d0)
  matbox_l(ix,iy,iz)=1.d0/(1.d0+elfc(ix,iy,iz)**2/elfcuni(ix,iy,iz)**2)
end do
end do
end do

call comm_summation(matbox_l,elf,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

elp3(817)=get_wtime()
elp3(847)=elp3(847)+elp3(817)-elp3(816)

end subroutine calcELF
