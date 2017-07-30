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
subroutine setmg(mg_sta,mg_end,mg_num,mg_sta_all,mg_end_all,mg_num_all,  &
                 lg_sta,lg_num,nproc,nproc_id_global,nproc_Mxin,nproc_ob,isequential)
implicit none
integer :: i1,j1,j2,j3
integer :: ibox
integer :: isequential
integer :: nproc,nproc_id_global
integer :: nproc_Mxin(3)
integer :: nproc_ob
integer :: mg_sta(3),mg_end(3),mg_num(3)
integer :: lg_sta(3),lg_num(3)
integer :: mg_sta_all(3,0:nproc-1),mg_end_all(3,0:nproc-1),mg_num_all(3,0:nproc-1)
integer :: nproc_Mxin_mul

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
if(isequential==1)then
  do j3=0,nproc_Mxin(3)-1
  do j2=0,nproc_Mxin(2)-1
  do j1=0,nproc_Mxin(1)-1
    do i1=0,nproc_ob-1
      ibox = i1 + nproc_ob*( j1 + nproc_Mxin(1)*j2 + nproc_Mxin(1)*nproc_Mxin(2)*j3 )
      mg_sta_all(1,ibox)=j1*lg_num(1)/nproc_Mxin(1)+lg_sta(1)
      mg_end_all(1,ibox)=(j1+1)*lg_num(1)/nproc_Mxin(1)+lg_sta(1)-1
      mg_sta_all(2,ibox)=j2*lg_num(2)/nproc_Mxin(2)+lg_sta(2)
      mg_end_all(2,ibox)=(j2+1)*lg_num(2)/nproc_Mxin(2)+lg_sta(2)-1
      mg_sta_all(3,ibox)=j3*lg_num(3)/nproc_Mxin(3)+lg_sta(3)
      mg_end_all(3,ibox)=(j3+1)*lg_num(3)/nproc_Mxin(3)+lg_sta(3)-1
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i1=0,nproc_ob-1
    do j3=0,nproc_Mxin(3)-1
    do j2=0,nproc_Mxin(2)-1
    do j1=0,nproc_Mxin(1)-1
      ibox = j1 + nproc_Mxin(1)*j2 + nproc_Mxin(1)*nproc_Mxin(2)*j3 + i1*nproc_Mxin_mul
      mg_sta_all(1,ibox)=j1*lg_num(1)/nproc_Mxin(1)+lg_sta(1)
      mg_end_all(1,ibox)=(j1+1)*lg_num(1)/nproc_Mxin(1)+lg_sta(1)-1
      mg_sta_all(2,ibox)=j2*lg_num(2)/nproc_Mxin(2)+lg_sta(2)
      mg_end_all(2,ibox)=(j2+1)*lg_num(2)/nproc_Mxin(2)+lg_sta(2)-1
      mg_sta_all(3,ibox)=j3*lg_num(3)/nproc_Mxin(3)+lg_sta(3)
      mg_end_all(3,ibox)=(j3+1)*lg_num(3)/nproc_Mxin(3)+lg_sta(3)-1
    end do
    end do
    end do
  end do
end if
mg_num_all(:,:)=mg_end_all(:,:)-mg_sta_all(:,:)+1

mg_sta(:)=mg_sta_all(:,nproc_id_global)
mg_end(:)=mg_end_all(:,nproc_id_global)
mg_num(:)=mg_num_all(:,nproc_id_global)

end subroutine setmg
