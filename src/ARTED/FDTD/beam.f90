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

!===============================================================
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
                            & NXvacL_m, NXvacR_m, NYvacB_m, NYvacT_m, &
                            & HX_m, HY_m, Ac_m, Ac_new_m, dt, &
                            & pi, c_light
  use salmon_math, only: bessel_j1_salmon
  implicit none
  real(8) :: f0_1, wpulse_1
  integer :: ix_m, iy_m
  real(8) :: lx, ly, x, y, kx, ky, k, vx
  real(8) :: f(3), j, tau
  
  real(8) sin2cos
  
  ! First pulse
  if(rlaser_int_wcm2_1 < 0d0)then
    f0_1 = amplitude1
  else
    f0_1=5.338d-9*sqrt(rlaser_int_wcm2_1)      ! electric field in a.u.
  end if
!  omega_1 = omegaeV_1 / (2d0*13.6058d0)  ! frequency in a.u.
!  tpulse_1 = tpulsefs_1 / 0.02418d0 ! pulse_duration in a.u.
  wpulse_1 = 2*pi/pulse_tw1

  lx = (NXvacR_m-NXvacL_m) * HX_m
  ly = (NYvacT_m+1) * HY_m

  ky = 3.8317d0 / ly
  k = omega1 / c_light
  
  if (ky > k) then
    call err_finalize('Error: NY,HY is not enough')
  endif
  
  kx = SQRT(k*k-ky*ky)
  vx = omega1 / kx

  Ac_m = 0.0
  Ac_new_m = 0.0
  do iy_m = NYvacB_m, NYvacT_m
     y = HY_m * iy_m
     j = bessel_j1_salmon(ky * y)
     f = j * f0_1 / 0.58186d0 * Epdir_re1
     do ix_m = NXvacL_m-1, 0
       x = ix_m * HX_m
       ! Ac_new_m
       tau = - x / vx
       Ac_new_m(:,ix_m,iy_m) = f * sin2cos(tau, pulse_tw1, omega1, phi_CEP1)
       ! Ac_m (previous time-step)
       tau = - x / vx - dt
       Ac_m(:,ix_m,iy_m) =  f * sin2cos(tau, pulse_tw1, omega1, phi_CEP1)
     end do
  end do
  
  ! impose boundary condition
  Ac_new_m(:,:,NYvacT_m+1) = 0.0d0
  Ac_new_m(2:3,:,NYvacB_m-1) = 0.0d0
  Ac_new_m(1,:,NYvacB_m-1) = Ac_new_m(1,:,NYvacB_m)
  Ac_m(:,:,NYvacT_m+1) = 0.0d0
  Ac_m(2:3,:,NYvacB_m-1) = 0.0d0
  Ac_m(1,:,NYvacB_m-1) = Ac_m(1,:,NYvacB_m)
  return
end subroutine incident_bessel_beam


subroutine read_initial_ac_from_file()
  use Global_Variables, only: SYSName, directory,file_ac_init, &
                            & Ac_m, Ac_new_m
  use salmon_parallel
  use salmon_communication, only: comm_is_root, comm_bcast
  implicit none
  integer :: ix_m, iy_m
  integer :: nx1_m, nx2_m, ny1_m, ny2_m
  
  write(file_ac_init, "(A,A,'_Ac_init.dat')") trim(directory), trim(SYSname)
  if (comm_is_root(nproc_id_global)) then
    Ac_m = 0.0
    Ac_new_m = 0.0
    open(944, file=trim(file_ac_init))
    read(944, *) nx1_m, nx2_m
    read(944, *) ny1_m, ny2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        read(944, *) Ac_m(:, ix_m, iy_m), Ac_new_m(:, ix_m, iy_m)
      end do
    end do
    close(944)
  end if
  call comm_bcast(Ac_m, nproc_group_global)
  call comm_bcast(Ac_new_m, nproc_group_global)
  return
end subroutine 
