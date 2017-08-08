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
!This file is "Ylm_dYlm.f"
!This file contain two functions.
!Function Ylm
!Function dYlm
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! This function Ylm is related to the real spherical harmonics Ylm0 by
!c Ylm=sqrt(4*pi/(2l+1))*r**l*Ylm0 and is a monomial of x,y,z
Function Ylm(x,y,z,il,im)
  implicit none
  real(8),intent(IN) :: x,y,z
  integer,intent(IN) :: il,im
  integer :: ilm
  real(8) :: Ylm,r2

  ilm=il*il+il+1+im
  r2=x*x+y*y+z*z
  select case( ilm )
  case(1)  ; Ylm=1.d0                                     ! ilm=1  (0  0)
  case(2)  ; Ylm=-y                                       ! ilm=2  (1 -1)
  case(3)  ; Ylm=z                                        ! ilm=3  (1  0)
  case(4)  ; Ylm=-x                                       ! ilm=4  (1  1)
  case(5)  ; Ylm=sqrt(3.d0)*x*y                           ! ilm=5  (2 -2)
  case(6)  ; Ylm=-sqrt(3.d0)*y*z                          ! ilm=6  (2 -1)
  case(7)  ; Ylm=(3.d0*z*z-r2)/2.d0                       ! ilm=7  (2  0)
  case(8)  ; Ylm=-sqrt(3.d0)*x*z                          ! ilm=8  (2  1)
  case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)                ! ilm=9  (2  2)
  case(10) ; Ylm=-sqrt(5.d0/8.d0)*y*(3*x*x-y*y)           ! ilm=10 (3 -3)
  case(11) ; Ylm=sqrt(15.d0)*x*y*z                        ! ilm=11 (3 -2)
  case(12) ; Ylm=-sqrt(3.d0/8.d0)*y*(5.d0*z*z-r2)         ! ilm=12 (3 -1)
  case(13) ; Ylm=z*(5.d0*z*z-3*r2)/2.d0                   ! ilm=13 (3  0)
  case(14) ; Ylm=-sqrt(3.d0/8.d0)*x*(5.d0*z*z-r2)         ! ilm=14 (3  1)
  case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)             ! ilm=15 (3  2)
  case(16) ; Ylm=-sqrt(5.d0/8.d0)*x*(x*x-3.d0*y*y)        ! ilm=16 (3  3)
  case(17) ; Ylm=sqrt(35.d0)/2.d0*x*y*(x*x-y*y)           ! ilm=17 (4 -4)
  case(18) ; Ylm=-sqrt(35.d0/8.d0)*y*z*(3.d0*x*x-y*y)     ! ilm=18 (4 -3)
  case(19) ; Ylm=sqrt(5.d0)/2.d0*x*y*(7.d0*z*z-r2)        ! ilm=19 (4 -2)
  case(20) ; Ylm=-sqrt(5.d0/8.d0)*y*z*(7.d0*z*z-3.d0*r2)  ! ilm=20 (4 -1)
  case(21) ; Ylm=(35.d0*z**4-30.d0*z*z*r2+3.d0*r2*r2)/8.d0! ilm=21 (4  0)
  case(22) ; Ylm=-sqrt(5.d0/8.d0)*x*z*(7.d0*z**2-3.d0*r2) ! ilm=22 (4  1)
  case(23) ; Ylm=sqrt(5.d0)/4.d0*(7.d0*z*z-r2)*(x*x-y*y)  ! ilm=23 (4  2)
  case(24) ; Ylm=-sqrt(35.d0/8.d0)*x*z*(x*x-3.d0*y*y)     ! ilm=24 (4  3)
  case(25) ; Ylm=sqrt(35.d0)/8.d0*(x**4+y**4-6.d0*x*x*y*y)! ilm=25 (4  4)
  end select

  return
End Function Ylm
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Function dYlm(x,y,z,il,im,idir)
  implicit none
  real(8),intent(IN) :: x,y,z
  integer,intent(IN) :: il,im,idir
  integer :: ilm
  real(8) :: dYlm

  ilm=il*il+il+1+im
  if ((ilm > 16) .or. (il >= 4)) then
    write(*,*) 'dYlm routine not prepared for il>=4'
    stop
  endif
  if((ilm ==  1 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=1  (0  0)
  if((ilm ==  1 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=1  (0  0)
  if((ilm ==  1 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=1  (0  0)
  if((ilm ==  2 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=2  (1 -1)
  if((ilm ==  2 ).and.(idir == 2)) dYlm=-1.d0                           ! ilm=2  (1 -1)
  if((ilm ==  2 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=2  (1 -1)
  if((ilm ==  3 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=3  (1  0)
  if((ilm ==  3 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=3  (1  0)
  if((ilm ==  3 ).and.(idir == 3)) dYlm=1.d0                            ! ilm=3  (1  0)
  if((ilm ==  4 ).and.(idir == 1)) dYlm=-1.d0                           ! ilm=4  (1  1)
  if((ilm ==  4 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=4  (1  1)
  if((ilm ==  4 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=4  (1  1)
  if((ilm ==  5 ).and.(idir == 1)) dYlm=sqrt(3.d0)*y                    ! ilm=5  (2 -2)
  if((ilm ==  5 ).and.(idir == 2)) dYlm=sqrt(3.d0)*x                    ! ilm=5  (2 -2)
  if((ilm ==  5 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=5  (2 -2)
  if((ilm ==  6 ).and.(idir == 1)) dYlm=0.d0                            ! ilm=6  (2 -1)
  if((ilm ==  6 ).and.(idir == 2)) dYlm=-sqrt(3.d0)*z                   ! ilm=6  (2 -1)
  if((ilm ==  6 ).and.(idir == 3)) dYlm=-sqrt(3.d0)*y                   ! ilm=6  (2 -1)
  if((ilm ==  7 ).and.(idir == 1)) dYlm=-x                              ! ilm=7  (2  0)
  if((ilm ==  7 ).and.(idir == 2)) dYlm=-y                              ! ilm=7  (2  0)
  if((ilm ==  7 ).and.(idir == 3)) dYlm=2.d0*z                          ! ilm=7  (2  0)
  if((ilm ==  8 ).and.(idir == 1)) dYlm=-sqrt(3.d0)*z                   ! ilm=8  (2  1)
  if((ilm ==  8 ).and.(idir == 2)) dYlm=0.d0                            ! ilm=8  (2  1)
  if((ilm ==  8 ).and.(idir == 3)) dYlm=-sqrt(3.d0)*x                   ! ilm=8  (2  1)
  if((ilm ==  9 ).and.(idir == 1)) dYlm=sqrt(3.d0)*x                    ! ilm=9  (2  2)
  if((ilm ==  9 ).and.(idir == 2)) dYlm=-sqrt(3.d0)*y                   ! ilm=9  (2  2)
  if((ilm ==  9 ).and.(idir == 3)) dYlm=0.d0                            ! ilm=9  (2  2)
  if((ilm == 10 ).and.(idir == 1)) dYlm=-3.d0*sqrt(5.d0/2.d0)*x*y                         ! ilm=10 (3 -3)
  if((ilm == 10 ).and.(idir == 2)) dYlm=-3.d0/2.d0*sqrt(5.d0/2.d0)*(x*x-y*y)               ! ilm=10 (3 -3)
  if((ilm == 10 ).and.(idir == 3)) dYlm=0.d0                                              ! ilm=10 (3 -3)
  if((ilm == 11 ).and.(idir == 1)) dYlm=sqrt(15.d0)*y*z                                   ! ilm=11 (3 -2)
  if((ilm == 11 ).and.(idir == 2)) dYlm=sqrt(15.d0)*z*x                                   ! ilm=11 (3 -2)
  if((ilm == 11 ).and.(idir == 3)) dYlm=sqrt(15.d0)*x*y                                   ! ilm=11 (3 -2)
  if((ilm == 12 ).and.(idir == 1)) dYlm=sqrt(3.d0/2.d0)*x*y                               ! ilm=12 (3 -1)
  if((ilm == 12 ).and.(idir == 2)) dYlm=-1.d0/2.d0*sqrt(3.d0/2.d0)*(4.d0*z*z-x*x-3.d0*y*y)! ilm=12 (3 -1)
  if((ilm == 12 ).and.(idir == 3)) dYlm=-4.d0*sqrt(3.d0/2.d0)*y*z                         ! ilm=12 (3 -1)
  if((ilm == 13 ).and.(idir == 1)) dYlm=-3.d0*z*x                                         ! ilm=13 (3  0)
  if((ilm == 13 ).and.(idir == 2)) dYlm=-3.d0*y*z                                         ! ilm=13 (3  0)
  if((ilm == 13 ).and.(idir == 3)) dYlm=3.d0/2.d0*(2.d0*z*z-x*x-y*y)                      ! ilm=13 (3  0)
  if((ilm == 14 ).and.(idir == 1)) dYlm=-1.d0/2.d0*sqrt(3.d0/2.d0)*(4.d0*z*z-3.d0*x*x-y*y)! ilm=14 (3  1)
  if((ilm == 14 ).and.(idir == 2)) dYlm=sqrt(3.d0/2.d0)*x*y                               ! ilm=14 (3  1)
  if((ilm == 14 ).and.(idir == 3)) dYlm=-4.d0*sqrt(3.d0/2.d0)*z*x                         ! ilm=14 (3  1)
  if((ilm == 15 ).and.(idir == 1)) dYlm=sqrt(15.d0)*z*x                                   ! ilm=15 (3  2)
  if((ilm == 15 ).and.(idir == 2)) dYlm=-sqrt(15.d0)*y*z                                  ! ilm=15 (3  2)
  if((ilm == 15 ).and.(idir == 3)) dYlm=1.d0/2.d0*sqrt(15.d0)*(x*x-y*y)                   ! ilm=15 (3  2)
  if((ilm == 16 ).and.(idir == 1)) dYlm=-3.d0/2.d0*sqrt(5.d0/2.d0)*(x*x-y*y)              ! ilm=16 (3  3)
  if((ilm == 16 ).and.(idir == 2)) dYlm=3.d0*sqrt(5.d0/2.d0)*x*y                          ! ilm=16 (3  3)
  if((ilm == 16 ).and.(idir == 3)) dYlm=0.d0                                              ! ilm=16 (3  3)

  return
End Function dYlm
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

