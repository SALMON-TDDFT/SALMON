!
!  Copyright 2017-2019 SALMON developers
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
!================================================= Initial wave function

SUBROUTINE init_wf_ns(ifunc)
use scf_data
implicit none

integer :: ik,iob,iseed,a,ix,iy,iz
integer :: is,iss,pstart(2),pend(2)
real(8) :: xx,yy,zz,x1,y1,z1,rr,rnd,Xmax,Ymax,Zmax
integer :: iob_myob
integer :: icorr_p
integer :: ifunc
integer :: icheck

if(ilsda == 0)then
  iss=1
  pstart(1)=1
  pend(1)=itotMST
else if(ilsda == 1)then
  iss=2
  pstart(1)=1
  pend(1)=MST(1)
  pstart(2)=MST(1)+1
  pend(2)=itotMST
end if

select case(iperiodic)
case(0)
  Xmax=0.d0 ; Ymax=0.d0 ; Zmax=0.d0
  if(iflag_ps.eq.1)then
    do a=1,MI
      if ( abs(Rion(1,a)) > Xmax ) Xmax=abs(Rion(1,a))
      if ( abs(Rion(2,a)) > Ymax ) Ymax=abs(Rion(2,a))
      if ( abs(Rion(3,a)) > Zmax ) Zmax=abs(Rion(3,a))
    end do
  end if
  
  Xmax=Xmax+1.d0/a_B ; Ymax=Ymax+1.d0/a_B ; Zmax=Zmax+1.d0/a_B
  
  iseed=123
  do is=1,iss
  do iob=pstart(is),pend(is)
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,1,icorr_p)
    call quickrnd_ns ; x1=Xmax*(2.d0*rnd-1.d0)
    call quickrnd_ns ; y1=Ymax*(2.d0*rnd-1.d0)
    call quickrnd_ns ; z1=Zmax*(2.d0*rnd-1.d0)
    call check_init_wf(icheck)
    if(icheck==1.and.icorr_p==1)then
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        xx=gridcoo(ix,1) ; yy=gridcoo(iy,2) ; zz=gridcoo(iz,3)
        rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
        psi(ix,iy,iz,iob_myob,1)=exp(-0.5d0*(rr*a_B)**2)*(a_B)**(3/2)
      end do
      end do
      end do
    end if
  end do
  end do
case(3)
  Xmax=rLsize(1,1)
  Ymax=rLsize(2,1)
  Zmax=rLsize(3,1)

  iseed=123
  do is=1,iss
  do ik=1,num_kpoints_rd
  do iob=pstart(is),pend(is)
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,ik,icorr_p)
    call quickrnd_ns ; x1=Xmax*rnd
    call quickrnd_ns ; y1=Ymax*rnd
    call quickrnd_ns ; z1=Zmax*rnd
    call check_init_wf(icheck)
    if(icheck==1.and.icorr_p==1)then
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        xx=gridcoo(ix,1) ; yy=gridcoo(iy,2) ; zz=gridcoo(iz,3)
        rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
        zpsi(ix,iy,iz,iob_myob,ik)=exp(-0.5d0*rr**2)
      end do
      end do
      end do
    end if
  end do
  end do
  end do
end select
  
return

CONTAINS

  subroutine quickrnd_ns
  implicit none
  integer,parameter :: im=6075,ia=106,ic=1283
  iseed=mod(iseed*ia+ic,im) ; rnd=real(iseed,8)/real(im,8)
  end subroutine quickrnd_ns

  subroutine check_init_wf(icheck)
  implicit none
  integer :: icheck
    if(ifunc==1)then
      icheck=1
    else if(ifunc==2)then
      if(ilsda==0)then
        if(iob>MST0(1))then
          icheck=1
        else
          icheck=0 
        end if
      else if(ilsda==1)then
        if((iob<=MST(1).and.iob>MST0(1)).or.(iob>MST(1).and.iob>=MST(1)+MST0(2)))then
          icheck=1
        else
          icheck=0
        end if
      end if
    end if
        
  end subroutine check_init_wf

END SUBROUTINE init_wf_ns


