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
subroutine check_dos_pdos
  use inputoutput, only: ispin, out_dos, out_pdos, out_dos_method
  implicit none

  if(out_dos=='y'.or.out_pdos=='y')then
    if(ispin==1)then
      stop "Sorry, out_dos is not implemented for ispin=1."
    end if
    select case(out_dos_method)
    case("gaussian","lorentzian")
      continue
    case default
      stop 'set out_dos_meshotd to "gaussian" or "lorentzian"'
    end select 
  end if

end subroutine check_dos_pdos

