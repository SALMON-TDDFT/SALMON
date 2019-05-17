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
!=========================================== Spherical harmonic function

SUBROUTINE Ylm_sub(xx,yy,zz,lm,Ylm)
!$ use omp_lib

real(8),intent(IN) :: xx,yy,zz
integer,intent(IN) :: lm
real(8) :: rrrr
real(8) :: Ylm

select case( lm )

case(1)  ; Ylm=1.d0
case(2)  ; Ylm=yy
case(3)  ; Ylm=zz
case(4)  ; Ylm=xx
case(5)  ; Ylm=sqrt(3.d0)*xx*yy                         ! lm=5  (2 -2)
case(6)  ; Ylm=sqrt(3.d0)*yy*zz                         ! lm=6  (2 -1)
case(7)  ; Ylm=(2*zz*zz-xx*xx-yy*yy)/2.d0               ! lm=7  (2 0)
case(8)  ; Ylm=sqrt(3.d0)*xx*zz                         ! lm=8  (2 1)
case(9)  ; Ylm=sqrt(3.d0/4.d0)*(xx*xx-yy*yy)            ! lm=9  (2 2)
case(10) ; Ylm=sqrt(5.d0/8.d0)*yy*(3*xx*xx-yy*yy)       ! lm=10 (3 -3)
case(11) ; Ylm=sqrt(15.d0)*xx*yy*zz                     ! lm=11 (3 -2)
case(12) ; Ylm=sqrt(3.d0/8.d0)*yy*(4*zz*zz-xx*xx-yy*yy) ! lm=12 (3 -1)
case(13) ; Ylm=zz*(2*zz*zz-3*xx*xx-3*yy*yy)/2.d0        ! lm=13 (3 0)
case(14) ; Ylm=sqrt(3.d0/8.d0)*xx*(4*zz*zz-xx*xx-yy*yy) ! lm=14 (3 1)
case(15) ; Ylm=sqrt(15.d0/4.d0)*zz*(xx*xx-yy*yy)        ! lm=15 (3 2)
case(16) ; Ylm=sqrt(5.d0/8.d0)*xx*(xx*xx-3*yy*yy)       ! lm=16 (3 3)

! for Hartree routine   

case(17)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(35.d0)/2.d0*xx*yy*(xx**2-yy**2)
case(18)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(35.d0/8.d0)*yy*zz*(3*xx**2-yy**2)
case(19)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(5.d0)/2.d0*xx*yy*(7*zz**2-1)
case(20)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(5.d0/8.d0)*yy*(7*zz**3-3*zz)
case(21)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*(35*zz**4-30*zz**2+3.d0)/8.d0
case(22)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(5.d0/8.d0)*xx*(7*zz**3-3*zz)
case(23)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(5.d0)/4.d0*(7*zz**2-1)*(xx**2-yy**2)
case(24)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(35.d0/8.d0)*zz*xx*(xx**2-3*yy**2)
case(25)
   rrrr=(xx**2+yy**2+zz**2)**2
   Ylm=rrrr*sqrt(35.d0)/8.d0*(xx**4+yy**4-6*xx**2*yy**2)

end select

END SUBROUTINE Ylm_sub

