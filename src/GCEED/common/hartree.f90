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
!============================ Hartree potential (Solve Poisson equation)
SUBROUTINE Hartree_ns
use scf_data

implicit none

if(iSCFRT==1)then
  select case(iperiodic)
  case(0)
    call Hartree_cg(rho,Vh)
  case(3)
    select case(iflag_hartree)
    case(2)
      call Hartree_periodic(rho,Vh)
    case(4)
      call Hartree_FFTE(rho,Vh)
    end select
  end select
else if(iSCFRT==2)then
  select case(iperiodic)
  case(0)
    if(mod(itt,2)==1)then
      call Hartree_cg(rho,Vh_stock2)
    else
      call Hartree_cg(rho,Vh_stock1)
    end if
  case(3)
    if(mod(itt,2)==1)then
      select case(iflag_hartree)
      case(2)
        call Hartree_periodic(rho,Vh_stock2)
      case(4)
        call Hartree_FFTE(rho,Vh_stock2)
      end select
    else
      select case(iflag_hartree)
      case(2)
        call Hartree_periodic(rho,Vh_stock1)
      case(4)
        call Hartree_FFTE(rho,Vh_stock1)
      end select
    end if
  end select
end if

return

END SUBROUTINE Hartree_ns
