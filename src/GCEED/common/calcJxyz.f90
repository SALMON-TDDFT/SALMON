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
SUBROUTINE calcJxyz
use scf_data
use allocate_psl_sub
implicit none
integer :: iatom,jj,ibox,ikoa,lm

uV=0.d0
do iatom=1,MI
  ikoa=Kion(iatom)
  ibox=0
  do jj=1,Mps_all(iatom)
    if(Jxyz_all(1,jj,iatom)>=mg_sta(1).and.Jxyz_all(1,jj,iatom)<=mg_end(1).and.  &
       Jxyz_all(2,jj,iatom)>=mg_sta(2).and.Jxyz_all(2,jj,iatom)<=mg_end(2).and.  &
       Jxyz_all(3,jj,iatom)>=mg_sta(3).and.Jxyz_all(3,jj,iatom)<=mg_end(3)) then
      ibox=ibox+1
      do lm=1,(Mlps(ikoa)+1)**2
        uV(ibox,lm,iatom)=uV_all(jj,lm,iatom)
      end do
    end if
  end do
end do

return

END SUBROUTINE calcJxyz
