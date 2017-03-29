! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine setlg(lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
                 Hgs,Nd,rLsize1,imesh_oddeven)
implicit none
integer :: lg_sta(3),lg_end(3),lg_num(3)
integer :: ista_Mx_ori(3),iend_Mx_ori(3),inum_Mx_ori(3)
integer :: img
integer :: Nd
integer :: imesh_oddeven
real(8) :: Hgs(3)
real(8) :: rLsize1(3)
real(8),parameter :: epsilon=1.d-10

iend_Mx_ori(:)=int((rLsize1(:)+epsilon)/2.d0/Hgs(:))+Nd
lg_end(:)=int((rLsize1(:)+epsilon)/2.d0/Hgs(:))
if(imesh_oddeven==1)then
  ista_Mx_ori(:)=-(int((rLsize1(:)+epsilon)/2.d0/Hgs(:))+Nd)
  lg_sta(:)=-(int((rLsize1(:)+epsilon)/2.d0/Hgs(:)))
else if(imesh_oddeven==2)then
  ista_Mx_ori(:)=-(int((rLsize1(:)+epsilon)/2.d0/Hgs(:))+Nd)+1
  lg_sta(:)=-(int((rLsize1(:)+epsilon)/2.d0/Hgs(:)))+1
end if

inum_Mx_ori(:)=iend_Mx_ori(:)-ista_Mx_ori(:)+1
lg_num(:)=lg_end(:)-lg_sta(:)+1

end subroutine setlg
