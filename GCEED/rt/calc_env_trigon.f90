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
!=======================================================================
!=======================================================================

subroutine calc_env_trigon(ipulse,tenv_trigon)
  !$ use omp_lib
  use inputoutput
  use scf_data
  
  implicit none
  integer,intent(in)  :: ipulse  ! 1: first pulse, 2: second pulse
  real(8),intent(out) :: tenv_trigon

  real(8) :: alpha1,alpha2
  real(8) :: theta1,theta2
  character(16) :: tae_shape

  if(ipulse==1)then
    tae_shape=ae_shape1
    ! sin(alpha1*theta1)**2
    alpha1=Pi/pulse_tw1
    theta1=dble(itt)*dt
    ! sin(alpha2*theta2) or cos(alhpa2*theta2)
    alpha2=omega1
    theta2=dble(itt)*dt-0.5d0*pulse_tw1+phi_cep1*2d0*pi
  else if(ipulse==2)then
    tae_shape=ae_shape2
    ! sin(alpha1*theta1)**2
    alpha1=Pi/pulse_tw2
    theta1=dble(itt)*dt-t1_t2
    ! sin(alpha2*theta2) or cos(alhpa2*theta2)
    alpha2=omega2
    theta2=dble(itt)*dt-t1_t2-0.5d0*pulse_tw2+phi_cep2*2d0*pi
  end if

  if(tae_shape=='esin2sin')then
    tenv_trigon=sin(alpha1*theta1)**2*sin(alpha2*theta2)
!  else if(tae_shape=='esin2cos')then
!    tenv_trigon=sin(alpha1*theta1)**2*cos(alpha2*theta2)
!  else if(tae_shape=='asin2sin')then
!    tenv_trigon=alpha1*sin(2.d0*alpha1*theta1)*sin(alpha2*theta2)   &
!               +alpha2*sin(alpha1*theta1)**2*cos(alpha2*theta2)
  else if(tae_shape=='asin2cos')then
    tenv_trigon=alpha1*sin(2.d0*alpha1*theta1)*cos(alpha2*theta2)   &
               -alpha2*sin(alpha1*theta1)**2*sin(alpha2*theta2)
  end if
    
  return
  
end subroutine calc_env_trigon
