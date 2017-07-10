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
!-----------------------------------------------------------------------------------------
module salmon_math
  implicit none
  private

  real(8),parameter :: Pi=3.141592653589793d0

  public :: erf_salmon, &
            erfc_salmon, &
            fbessel_j1
contains
!--------------------------------------------------------------------------------
!! Error function and its complement are implemented based on the reference
!! W. J. Cody, Math. Comp. 23 631 (1969)
  function erf_salmon(x) result(y)
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: xa

    xa = abs(x)
    if(xa < 0.46875d0)then
       y = sign( erf_salmon_short(xa) , x )
    else if(xa < 4d0)then
       y = erfc_salmon_mid(xa)
       y = sign(1d0 - y,x)
    else
       y = erfc_salmon_long(xa)
       y = sign(1d0 - y,x)
    end if

  end function erf_salmon
!--------------------------------------------------------------------------------
  function erfc_salmon(x) result(y)
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: xa

    xa = abs(x)
    if(xa < 0.46875d0)then
       y = 1d0 - erf_salmon(xa)
    else if(xa < 4d0)then
       y = erfc_salmon_mid(xa)
    else
       y = erfc_salmon_long(xa)
    end if

    if(x<0d0) y = 2d0 - y

  end function erfc_salmon
!--------------------------------------------------------------------------------
  function erf_salmon_short(x) result(y) ! 0d0 <=x<0.5d0
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2
    real(8),parameter :: &
       a(0:4) = (/  &
       & 3.209377589138469472562d3, &
       & 3.774852376853020208137d2, &
       & 1.138641541510501556495d2, &
       & 3.161123743870565596947d0, &
       & 1.857777061846031526730d-1 &
       &/), &
       b(0:4) = (/  &
       & 2.844236833439170622273d3, &
       & 1.282616526077372275645d3, &
       & 2.440246379344441733056d2, &
       & 2.360129095234412093499d1, &
       & 1d0 &
       &/)

    x2 = x**2
    y = x*(a(0) + a(1)*x2 + a(2)*x2**2 + a(3)*x2**3 + a(4)*x2**4 ) &
         /(b(0) + b(1)*x2 + b(2)*x2**2 + b(3)*x2**3 + b(4)*x2**4 )


  end function erf_salmon_short
!--------------------------------------------------------------------------------
  function erfc_salmon_mid(x) result(y) ! 0.46875d0 <x< 4d0
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2
    real(8),parameter :: &
       p(0:7) = (/  &
       & 3.004592610201616005d2, &
       & 4.519189537118729422d2, &
       & 3.393208167343436870d2, &
       & 1.529892850469404039d2, &
       & 4.316222722205673530d1, &
       & 7.211758250883093659d0, &
       & 5.641955174789739711d-1, &
       &-1.368648573827167067d-7 & 
       &/), &
       q(0:7) = (/ &
       & 3.004592609569832933d2, &
       & 7.909509253278980272d2, &
       & 9.313540948506096211d2, &
       & 6.389802644656311665d2, &
       & 2.775854447439876434d2, &
       & 7.700015293522947295d1, &
       & 1.278272731962942351d1, &
       & 1d0 &
       &/)

    x2 = x**2

    y = (p(0) + p(1)*x + p(2)*x**2 + p(3)*x**3 + p(4)*x**4 & 
         + p(5)*x**5 + p(6)*x**6 + p(7)*x**7)/ &
         (q(0) + q(1)*x + q(2)*x**2 + q(3)*x**3 + q(4)*x**4 & 
         + q(5)*x**5 + q(6)*x**6 + q(7)*x**7)
    y = exp(-x2)*y


  end function erfc_salmon_mid
!--------------------------------------------------------------------------------
  function erfc_salmon_long(x) result(y) ! 4d0 <x
    implicit none
    real(8),intent(in) :: x
    real(8) :: y
    real(8) :: x2,x2_i
    real(8),parameter :: &
       r(0:4) = (/ &
       & -2.99610707703542174d-3, &
       & -4.94730910623250734d-2, &
       & -2.26956593539686930d-1, &
       & -2.78661308609647788d-1, &
       & -2.23192459734184686d-2  &
       &/), &
       s(0:4) = (/ &
       & 1.06209230528467918d-2, &
       & 1.91308926107829841d-1, &
       & 1.05167510706793207d0,  &
       & 1.98733201817135256d0,  &
       & 1d0 &
       &/)

    x2 = x**2

    x2_i = 1d0/x2
    y = 1d0/sqrt(pi) &
         +x2_i*(r(0) + r(1)*x2_i + r(2)*x2_i**2 + r(3)*x2_i**3 + r(4)*x2_i**4) &
         /(s(0) + s(1)*x2_i + s(2)*x2_i**2 + s(3)*x2_i**3 + s(4)*x2_i**4) 
    y = exp(-x2)/x*y

  end function erfc_salmon_long
  
  
!--------------------------------------------------------------------------------
  real(8) function fbessel_j1(x)
    implicit none
    real(8), intent(in) :: x
    integer, parameter :: order = 30
    real(8) :: c, s
    integer :: m

    c = 0.5 * x
    s = c
    do m = 1, 30
      c = -0.25d0 * x * x / (m * (m + 1)) * c
      s = s + c
    end do
    fbessel_j1 = s
    return
  end function fbessel_j1
  
end module salmon_math
!--------------------------------------------------------------------------------
