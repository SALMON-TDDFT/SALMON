!
!  Copyright 2016 ARTED developers
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
!This file is "init_wf.f90"
!This file contain two subroutines
!SUBROUTINE init_wf
!SUBROUTINE quickrnd(iseed,rnd)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine init_wf
  use Global_Variables
  implicit none
  integer :: iseed,ib,ik,i
  real(8) :: r2,x1,y1,z1,rnd

  zu_GS=0.d0
  do ik=NK_s,NK_e
    iseed= 123 + ik
    do ib=1,NB
      call quickrnd(iseed,rnd)
      x1=aLx*rnd
      call quickrnd(iseed,rnd)
      y1=aLy*rnd
      call quickrnd(iseed,rnd)
      z1=aLz*rnd
      do i=1,NL
        r2=(Lx(i)*Hx-x1)**2+(Ly(i)*Hy-y1)**2+(Lz(i)*Hz-z1)**2
        zu_GS(i,ib,ik)=exp(-0.5d0*r2)
      enddo
    enddo
  enddo

  return
End Subroutine init_wf
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine quickrnd(iseed,rnd)
  implicit none
  integer,parameter :: im=6075,ia=106,ic=1283
  integer :: iseed
  real(8) :: rnd

  iseed=mod(iseed*ia+ic,im)
  rnd=dfloat(iseed)/dfloat(im)

  return
End Subroutine quickrnd
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130

