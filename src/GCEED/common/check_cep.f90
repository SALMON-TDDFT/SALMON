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
subroutine check_cep
  use inputoutput
  implicit none
  integer :: round_phi
  real(8) :: udp_phi  ! udp: under dicimal point

  round_phi=int((phi_cep1-0.25d0)*2.d0)
  udp_phi=(phi_cep1-0.25d0)*2.d0-round_phi  
  if(ae_shape1=="Ecos2".and.abs(udp_phi)>=1.d-12)then
    stop "phi_cep1 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape1."
  end if
  
  round_phi=int((phi_cep2-0.25d0)*2.d0)
  udp_phi=(phi_cep2-0.25d0)*2.d0-round_phi
  if(ae_shape2=="Ecos2".and.abs(udp_phi)>=1.d-12)then
    stop "phi_cep2 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape2."
  end if

  return

end subroutine check_cep
