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
module math_constants
  implicit none

  public :: set_math_constants
  public :: is_zero
  public :: is_nonzero

private
  real    :: eps_sp
  real(8) :: eps_dp

  interface is_zero
    module procedure is_zero_sp
    module procedure is_zero_dp
  end interface

  interface is_nonzero
    module procedure is_nonzero_sp
    module procedure is_nonzero_dp
  end interface

  interface init_epsilon
    module procedure init_epsilon_sp
    module procedure init_epsilon_dp
  end interface

contains

  subroutine set_math_constants
    implicit none
    call init_epsilon(eps_sp)
    call init_epsilon(eps_dp)
  end subroutine

  function is_zero_sp(x) result(r)
    implicit none
    real, intent(in) :: x
    logical :: r
    r = (abs(x) < eps_sp)
  end function

  function is_zero_dp(x) result(r)
    implicit none
    real(8), intent(in) :: x
    logical :: r
    r = (abs(x) < eps_dp)
  end function

  function is_nonzero_sp(x) result(r)
    implicit none
    real, intent(in) :: x
    logical :: r
    r = (.NOT. is_zero(x))
  end function

  function is_nonzero_dp(x) result(r)
    implicit none
    real(8), intent(in) :: x
    logical :: r
    r = (.NOT. is_zero(x))
  end function

  subroutine init_epsilon_sp(eps)
    implicit none
    real, intent(out) :: eps
    real :: x
    eps = 1.0
    x   = 2.0
    do while (x > 1.0)
      eps = eps * 0.5
      x   = eps + 1.0
    end do
  end subroutine

  subroutine init_epsilon_dp(eps)
    implicit none
    real(8), intent(out) :: eps
    real(8) :: x
    eps = 1.0d0
    x   = 2.0d0
    do while (x > 1.0d0)
      eps = eps * 0.5d0
      x   = eps + 1.0d0
    end do
  end subroutine

end module math_constants
!--------------------------------------------------------------------------------
  
