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
!SUBROUTINE Hartree_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
subroutine Hartree_FFTE(trho,tVh)
  use salmon_parallel, only: nproc_group_global
  use salmon_parallel, only: nproc_id_icommy, nproc_group_icommy
  use salmon_parallel, only: nproc_id_icommz, nproc_group_icommz
  use salmon_parallel, only: nproc_group_icommw
  use salmon_communication, only: comm_summation
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  use allocate_psl_sub
!$  use omp_lib
  implicit none
  integer :: ix,iy,iz
  integer :: iix,iiy,iiz
  integer :: iz_sta,iz_end,iy_sta,iy_end
  real(8) :: trho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: tVh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: inv_lgnum3
  complex(8),parameter :: zI=(0.d0,1.d0)
  integer :: n
  real(8) :: bLx,bLy,bLz
  complex(8) :: A_FFTE_tmp(1:lg_num(1),1:lg_num(2)/NPUY,1:lg_num(3)/NPUZ)
  integer :: ng_sta_2(3),ng_end_2(3),ng_num_2(3)
  integer :: lg_sta_2(3),lg_end_2(3),lg_num_2(3)

  lg_sta_2(1:3)=lg_sta(1:3)
  lg_end_2(1:3)=lg_end(1:3)
  lg_num_2(1:3)=lg_num(1:3)
  
  ng_sta_2(1:3)=ng_sta(1:3)
  ng_end_2(1:3)=ng_end(1:3)
  ng_num_2(1:3)=ng_num(1:3)

  bLx=2.d0*Pi/(Hgs(1)*dble(lg_num_2(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg_num_2(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg_num_2(3)))

  inv_lgnum3=1.d0/(lg_num_2(1)*lg_num_2(2)*lg_num_2(3))

  iz_sta=1
  iz_end=lg_num_2(3)/NPUZ
  iy_sta=1
  iy_end=lg_num_2(2)/NPUY
  
!  rhoe_G_tmp=0.d0

  if(icheck_ascorder==1)then
    if(NPUW==1)then
!$OMP parallel do private(iiz,iiy)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
          A_FFTE(1:lg_end(1),iy,iz)=trho(1:lg_end(1),iiy,iiz)
        end do
      end do
    else
      A_FFTE_tmp=0.d0
!$OMP parallel do private(iiz,iiy,ix)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
          do iix=ng_sta(1),ng_end(1)
            ix=iix-lg_sta(1)+1
            A_FFTE_tmp(ix,iy,iz)=trho(iix,iiy,iiz)
          end do
        end do
      end do
      call comm_summation(A_FFTE_tmp,A_FFTE,lg_num(1)*lg_num(2)/NPUY*lg_num(3)/NPUZ,nproc_group_icommw)
    end if
  else
!$OMP parallel do
    do iz = lg_sta(3),lg_end(3)
    do iy = lg_sta(2),lg_end(2)
    do ix = lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
!$OMP parallel do
    do iz = ng_sta(3),ng_end(3)
    do iy = ng_sta(2),ng_end(2)
    do ix = ng_sta(1),ng_end(1)
      matbox_l(ix,iy,iz)=trho(ix,iy,iz)
    end do
    end do
    end do

    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

!!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
        A_FFTE(1:lg_end(1),iy,iz)=matbox_l2(1:lg_end(1),iiy,iiz)
      end do
    end do
  end if

  CALL PZFFT3DV_MOD(A_FFTE,B_FFTE,lg_num_2(1),lg_num_2(2),lg_num_2(3),NPUY,NPUZ,0) 
  CALL PZFFT3DV_MOD(A_FFTE,B_FFTE,lg_num_2(1),lg_num_2(2),lg_num_2(3),NPUY,NPUZ,-1) 

!$OMP parallel do private(n)
  do iz=iz_sta,iz_end
    do iy=iy_sta,iy_end
      do ix=1,lg_num(1)
        n=(iz-1)*lg_num_2(2)/NPUY*lg_num_2(1)+(iy-1)*lg_num_2(1)+ix
        rhoe_G(n)=B_FFTE(ix,iy,iz)*inv_lgnum3
        B_FFTE(ix,iy,iz)=B_FFTE(ix,iy,iz)*coef_poisson(ix,iy,iz)
      end do
    end do
  end do
  if(nproc_id_icommz==0.and.nproc_id_icommy==0)then
    rhoe_G(1)=0.d0
  end if

  CALL PZFFT3DV_MOD(B_FFTE,A_FFTE,lg_num_2(1),lg_num_2(2),lg_num_2(3),NPUY,NPUZ,1)

  if(icheck_ascorder==1)then
    if(NPUW==1)then
!$OMP parallel do private(iiz,iiy)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
          tVh(1:lg_end(1),iiy,iiz)=A_FFTE(1:lg_end(1),iy,iz)
        end do
      end do
    else
!$OMP parallel do private(iiz,iiy,ix)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
          do iix=ng_sta(1),ng_end(1)
            ix=iix-lg_sta(1)+1
            tVh(iix,iiy,iiz)=A_FFTE(ix,iy,iz)
          end do
        end do
      end do
    end if
  else
!$OMP parallel do
    do iz = lg_sta(3),lg_end(3)
    do iy = lg_sta(2),lg_end(2)
    do ix = lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
!!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg_num_2(3)/NPUZ
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg_num_2(2)/NPUY
        matbox_l(1:lg_end(1),iiy,iiz)=A_FFTE(1:lg_end(1),iy,iz)
      end do
    end do
    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
!$OMP parallel do
    do iz = mg_sta(3),mg_end(3)
    do iy = mg_sta(2),mg_end(2)
    do ix = mg_sta(1),mg_end(1)
      tVh(ix,iy,iz)=matbox_l2(ix,iy,iz)
    end do
    end do
    end do
  end if 

  return
end subroutine Hartree_FFTE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
