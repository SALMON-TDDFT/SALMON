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

  real(8) :: alpha,beta
  real(8) :: theta1,theta2
  character(16) :: tae_shape

  if(ipulse==1)then
    tae_shape=ae_shape1
    ! cos(theta1)**2
    theta1=Pi/pulse_tw1*(dble(itt)*dt-0.5d0*pulse_tw1)
    alpha=Pi/pulse_tw1
    ! cos(theta2)
    theta2=omega1*(dble(itt)*dt-0.5d0*pulse_tw1)+phi_cep1*2d0*pi
    beta=omega1
  else if(ipulse==2)then
    tae_shape=ae_shape2
    ! cos(theta1)**2
    theta1=Pi/pulse_tw2*(dble(itt)*dt-t1_t2-0.5d0*pulse_tw2)
    alpha=Pi/pulse_tw2
    ! cos(theta2)
    theta2=omega2*(dble(itt)*dt-t1_t2-0.5d0*pulse_tw2)+phi_cep2*2d0*pi
    beta=omega2
  end if

  if(tae_shape=='Ecos2')then
    tenv_trigon=cos(theta1)**2*cos(theta2)
  else if(tae_shape=='Acos2')then
    tenv_trigon=-(-alpha*sin(2.d0*theta1)*cos(theta2)   &
                  -beta*cos(theta1)**2*sin(theta2))
  end if
    
  return
  
end subroutine calc_env_trigon
