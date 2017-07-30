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
subroutine setng(ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s,   &
                 nproc,nproc_id_global,nproc_Mxin,nproc_Mxin_s_dm,ista_Mxin,iend_Mxin,isequential)
implicit none
integer :: ii,i1,i2,i3,i4
integer :: ibox
integer :: nproc,nproc_id_global
integer :: ng_sta(3),ng_end(3),ng_num(3)
integer :: ista_Mxin_s(3,0:nproc-1),iend_Mxin_s(3,0:nproc-1),inum_Mxin_s(3,0:nproc-1)
integer :: ista_Mxin(3,0:nproc-1),iend_Mxin(3,0:nproc-1)
integer :: nproc_Mxin_mul,nproc_Mxin_mul_s_dm
integer :: isequential
integer :: nproc_Mxin(3),nproc_Mxin_s_dm(3)

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

if(isequential==1)then
  do ii=0,nproc_Mxin_mul-1
    do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox= i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)       &
              +i4*nproc_Mxin_mul_s_dm   &
              +ii*nproc/nproc_Mxin_mul
      ista_Mxin_s(1,ibox)=i1*(iend_Mxin(1,ibox)-ista_Mxin(1,ibox)+1)/nproc_Mxin_s_dm(1)+ista_Mxin(1,ibox)
      iend_Mxin_s(1,ibox)=(i1+1)*(iend_Mxin(1,ibox)-ista_Mxin(1,ibox)+1)/nproc_Mxin_s_dm(1)+ista_Mxin(1,ibox)-1
      ista_Mxin_s(2,ibox)=i2*(iend_Mxin(2,ibox)-ista_Mxin(2,ibox)+1)/nproc_Mxin_s_dm(2)+ista_Mxin(2,ibox)
      iend_Mxin_s(2,ibox)=(i2+1)*(iend_Mxin(2,ibox)-ista_Mxin(2,ibox)+1)/nproc_Mxin_s_dm(2)+ista_Mxin(2,ibox)-1
      ista_Mxin_s(3,ibox)=i3*(iend_Mxin(3,ibox)-ista_Mxin(3,ibox)+1)/nproc_Mxin_s_dm(3)+ista_Mxin(3,ibox)
      iend_Mxin_s(3,ibox)=(i3+1)*(iend_Mxin(3,ibox)-ista_Mxin(3,ibox)+1)/nproc_Mxin_s_dm(3)+ista_Mxin(3,ibox)-1
    end do
    end do
    end do
    end do
  end do
else if(isequential==2)then
  do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do ii=0,nproc_Mxin_mul-1
      ibox=ii+(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))*nproc_Mxin_mul   &
            +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      ista_Mxin_s(1,ibox)=i1*(iend_Mxin(1,ii)-ista_Mxin(1,ii)+1)/nproc_Mxin_s_dm(1)+ista_Mxin(1,ii)
      iend_Mxin_s(1,ibox)=(i1+1)*(iend_Mxin(1,ii)-ista_Mxin(1,ii)+1)/nproc_Mxin_s_dm(1)+ista_Mxin(1,ii)-1
      ista_Mxin_s(2,ibox)=i2*(iend_Mxin(2,ii)-ista_Mxin(2,ii)+1)/nproc_Mxin_s_dm(2)+ista_Mxin(2,ii)
      iend_Mxin_s(2,ibox)=(i2+1)*(iend_Mxin(2,ii)-ista_Mxin(2,ii)+1)/nproc_Mxin_s_dm(2)+ista_Mxin(2,ii)-1
      ista_Mxin_s(3,ibox)=i3*(iend_Mxin(3,ii)-ista_Mxin(3,ii)+1)/nproc_Mxin_s_dm(3)+ista_Mxin(3,ii)
      iend_Mxin_s(3,ibox)=(i3+1)*(iend_Mxin(3,ii)-ista_Mxin(3,ii)+1)/nproc_Mxin_s_dm(3)+ista_Mxin(3,ii)-1
    end do
  end do
  end do
  end do
  end do
end if
inum_Mxin_s(1:3,0:nproc-1)=iend_Mxin_s(1:3,0:nproc-1)-ista_Mxin_s(1:3,0:nproc-1)+1

ng_sta(:)=ista_Mxin_s(:,nproc_id_global)
ng_end(:)=iend_Mxin_s(:,nproc_id_global)
ng_num(:)=inum_Mxin_s(:,nproc_id_global)

end subroutine setng
