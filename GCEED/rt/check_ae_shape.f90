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
subroutine check_ae_shape
  use inputoutput
  implicit none

  select case(ae_shape1)
  case("impulse","esin2cos","asin2cos")
    continue
  case default
    stop 'set ae_shape1 to "impulse", "esin2cos", or "asin2cos"'
  end select 

  select case(ae_shape2)
  case("none","impulse","esin2cos","asin2cos")
    continue
  case default 
    stop 'set ae_shape2 to "none", "impulse", "esin2cos", or "asin2cos"'
  end select 

end subroutine check_ae_shape
