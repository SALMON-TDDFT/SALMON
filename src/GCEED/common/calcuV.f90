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
SUBROUTINE calcuV
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use scf_data
use allocate_psl_sub
implicit none
integer :: iatom,jj,lm

  integer :: a,i,ik,ix,iy,iz,l,m
  integer :: nl

  real(8) :: alx,aly,alz
  real(8) :: hx,hy,hz
  integer :: lx(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: ly(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: lz(lg_num(1)*lg_num(2)*lg_num(3))

  integer :: nlma
  integer :: lma
  character(17) :: property
  
  real(8),allocatable :: save_udVtbl_a(:,:,:)
  real(8),allocatable :: save_udVtbl_b(:,:,:)
  real(8),allocatable :: save_udVtbl_c(:,:,:)
  real(8),allocatable :: save_udVtbl_d(:,:,:)

  integer,allocatable :: lma_tbl(:,:)

  logical :: flag_use_grad_wf_on_force
  

  property='initial'
  flag_use_grad_wf_on_force=.false. 

  nl=lg_num(1)*lg_num(2)*lg_num(3)

  hx=Hgs(1)
  hy=Hgs(2)
  hz=Hgs(3)
  alx=Hgs(1)*dble(lg_num(1))
  aly=Hgs(2)*dble(lg_num(2))
  alz=Hgs(3)*dble(lg_num(3))

  if(iperiodic==0)then
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      i=(iz-lg_sta(3))*lg_num(1)*lg_num(2)+(iy-lg_sta(2))*lg_num(1)+ix-lg_sta(1)+1
      lx(i)=ix
      ly(i)=iy
      lz(i)=iz
    end do
    end do
    end do
  else if(iperiodic==3)then
    do iz=1,lg_num(3)
    do iy=1,lg_num(2)
    do ix=1,lg_num(1)
      i=(iz-1)*lg_num(1)*lg_num(2)+(iy-1)*lg_num(1)+ix
      lx(i)=ix-1
      ly(i)=iy-1
      lz(i)=iz-1
    end do
    end do
    end do
  end if

  lma=0
  do a=1,MI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  nlma=lma

  allocate(lma_tbl((pp%lmax+1)**2,MI))

  lma=0
  do a=1,MI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo


  allocate( save_udVtbl_a(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_b(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_c(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_d(pp%nrmax,0:pp%lmax,natom) )
     
  call calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
               nlma,lx,ly,lz,nl,hx,hy,hz,alx,aly,alz,  &
               lma_tbl,flag_use_grad_wf_on_force,property)

  lma = 0
  do iatom=1,MI
    ik=Kion(iatom)
    do l=0,Mlps(ik)
      if ( pp%inorm(l,ik)==0) then
        do lm=l**2+1,(l+1)**2
          do jj=1,ppg%mps(iatom)
            uV_all(jj,lm,iatom) = 0.d0
          end do
          uVu(lm,iatom)=1.d-10
        end do 
      else
        do lm=l**2+1,(l+1)**2
          lma = lma + 1
          do jj=1,ppg%mps(iatom)
            uV_all(jj,lm,iatom) = ppg%uv(jj,lma)
          end do
          uVu(lm,iatom)=pp%inorm(l,ik)
        end do 
      end if
    end do
  end do

return

END SUBROUTINE calcuV
