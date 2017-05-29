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
SUBROUTINE calcJxyz2nd
use scf_data
use allocate_psl_sub
implicit none
integer :: iatom,jj,ibox,ikoa,lm

Jxyz2nd=0.d0
uV2nd=0.d0
do iatom=1,MI
  ikoa=Kion(iatom)
  ibox=0
  do jj=1,Mps(iatom)
    if(Jxyz(1,jj,iatom)>=mg_sta(1).and.Jxyz(1,jj,iatom)<=mg_end(1).and.  &
       Jxyz(2,jj,iatom)>=mg_sta(2).and.Jxyz(2,jj,iatom)<=mg_end(2).and.  &
       Jxyz(3,jj,iatom)>=mg_sta(3).and.Jxyz(3,jj,iatom)<=mg_end(3)) then
      ibox=ibox+1
      jMps_l(ibox,iatom)=jj
      Jxyz2nd(1:3,ibox,iatom)=Jxyz(1:3,jj,iatom)
      Jxxyyzz2nd(1:3,ibox,iatom)=Jxxyyzz(1:3,jj,iatom)
      do lm=1,(Mlps(ikoa)+1)**2
        uV2nd(ibox,lm,iatom)=uV(jj,lm,iatom)
      end do
    end if
  end do
  max_jMps_l(iatom)=ibox
end do

if(iSCFRT==1)then
  do iatom=1,MI
    ibox=0
    do jj=1,Mps(iatom)
      if(Jxyz(1,jj,iatom)>=ng_sta(1).and.Jxyz(1,jj,iatom)<=ng_end(1).and.  &
         Jxyz(2,jj,iatom)>=ng_sta(2).and.Jxyz(2,jj,iatom)<=ng_end(2).and.  &
         Jxyz(3,jj,iatom)>=ng_sta(3).and.Jxyz(3,jj,iatom)<=ng_end(3)) then
        ibox=ibox+1
        jMps_l_s(ibox,iatom)=jj
      end if
    end do
    max_jMps_l_s(iatom)=ibox
  end do
end if

return

END SUBROUTINE calcJxyz2nd
