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

!===============================================================
!! TODO: Create new beam source function

real(8) function sin2cos(t, tw, omega, cep)
  use Global_Variables, only: pi
  implicit none
  real(8), intent(in) :: t, tw, omega, cep
  real(8) :: theta1, theta2

  if ((0 <= t) .and. (t <= tw)) then
    theta1 = pi / tw * t
    theta2 = omega * t + 2 * pi * cep
    sin2cos = sin(theta1) ** 2 * cos(theta2)
  else
    sin2cos = 0.0d0
  endif
  return
end function sin2cos
!===============================================================
subroutine incident_bessel_beam()
  use Global_Variables, only: amplitude1, rlaser_int_wcm2_1, omega1, pulse_tw1, &
                            & Epdir_re1, phi_CEP1, &
                            & nx_m, ny_m, &
                            & nx1_m, nx2_m, ny1_m, ny2_m, nz_origin_m, &
                            & my1_m, my2_m, &
                            & HX_m, HY_m, Ac_ms, Ac_new_ms, dt, &
                            & pi, c_light
  use salmon_math, only: bessel_j1_salmon
  implicit none
  real(8) :: f0_1, wpulse_1
  integer :: ix_m, iy_m, iz_m
  real(8) :: lx, ly, x, y, kx, ky, k, vx
  real(8) :: f(3), j, tau
  real(8) sin2cos
  
  iz_m = nz_origin_m
  
  ! First pulse
  if(rlaser_int_wcm2_1 < 0d0)then
    f0_1 = amplitude1
  else
    f0_1=5.338d-9*sqrt(rlaser_int_wcm2_1)      ! electric field in a.u.
  end if
!  omega_1 = omegaeV_1 / (2d0*13.6058d0)  ! frequency in a.u.
!  tpulse_1 = tpulsefs_1 / 0.02418d0 ! pulse_duration in a.u.
  wpulse_1 = 2*pi/pulse_tw1

  lx = NX_m * HX_m
  ly = NY_m * HY_m

  ky = 3.8317d0 / ly
  k = omega1 / c_light
  
  if (ky > k) then
    call err_finalize('Error: NY,HY is not enough')
  endif
  
  kx = SQRT(k*k-ky*ky)
  vx = omega1 / kx

  Ac_ms = 0.0
  Ac_new_ms = 0.0
  do iy_m = ny1_m, ny2_m
     y = HY_m * iy_m
     j = bessel_j1_salmon(ky * y)
     f = j * f0_1 / 0.58186d0 * Epdir_re1
     do ix_m = nx1_m, nx2_m
       x = ix_m * HX_m
       ! Ac_new_ms
       tau = - x / vx
       Ac_new_ms(:,ix_m,iy_m,iz_m) = f * sin2cos(tau, pulse_tw1, omega1, phi_CEP1)
       ! Ac_ms (previous time-step)
       tau = - x / vx - dt
       Ac_ms(:,ix_m,iy_m,iz_m) =  f * sin2cos(tau, pulse_tw1, omega1, phi_CEP1)
     end do
  end do
  
  ! impose boundary condition
  Ac_new_ms(:,:,my2_m,iz_m) = 0.0d0
  Ac_new_ms(2:3,:,my1_m,iz_m) = 0.0d0
  Ac_new_ms(1,:,my1_m,iz_m) = Ac_new_ms(1,:,ny1_m,iz_m)
  Ac_ms(:,:,my2_m,iz_m) = 0.0d0
  Ac_ms(2:3,:,my1_m,iz_m) = 0.0d0
  Ac_ms(1,:,my1_m,iz_m) = Ac_ms(1,:,ny1_m,iz_m)
  return
end subroutine incident_bessel_beam
