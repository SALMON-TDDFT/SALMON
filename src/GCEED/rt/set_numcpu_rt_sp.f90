!
!  Copyright 2018 SALMON developers
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
subroutine set_numcpu_rt_sp
  use salmon_parallel, only: nproc_size_global
  use scf_data
  use new_world_sub
  implicit none
  integer :: ii
  integer :: nproc_size_global_tmp
  integer :: nproc_mxin_tmp(3)
  
  integer :: num_factor2
  integer :: num_factor3
  integer :: num_factor5

  integer :: icount
 
  nproc_size_global_tmp=nproc_size_global
  
  ! this code treats the situation that nproc_size_global is less than or equal to 48,828,125
  
  num_factor2=0
  do ii=1,26
    if(mod(nproc_size_global_tmp,2)==0)then
      num_factor2=num_factor2+1
      nproc_size_global_tmp=nproc_size_global_tmp/2
    end if
  end do
  
  num_factor3=0
  do ii=1,17
    if(mod(nproc_size_global_tmp,3)==0)then
      num_factor3=num_factor3+1
      nproc_size_global_tmp=nproc_size_global_tmp/3
    end if
  end do
  
  num_factor5=0
  do ii=1,11
    if(mod(nproc_size_global_tmp,5)==0)then
      num_factor5=num_factor5+1
      nproc_size_global_tmp=nproc_size_global_tmp/5
    end if
  end do
  
  if(nproc_size_global_tmp/=1)then
    stop "In automatic process assignment, prime factors for number of processes must be combination of 2, 3 or 5."
  end if

  nproc_mxin_tmp(1:3)=1
 
  icount=0

  do ii=1,num_factor5
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_mxin_tmp(3)=nproc_mxin_tmp(3)*5
    else if(mod(icount,3)==2)then
      nproc_mxin_tmp(2)=nproc_mxin_tmp(2)*5
    else
      nproc_mxin_tmp(1)=nproc_mxin_tmp(1)*5
    end if
  end do

  do ii=1,num_factor3
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_mxin_tmp(3)=nproc_mxin_tmp(3)*3
    else if(mod(icount,3)==2)then
      nproc_mxin_tmp(2)=nproc_mxin_tmp(2)*3
    else
      nproc_mxin_tmp(1)=nproc_mxin_tmp(1)*3
    end if
  end do

  do ii=1,num_factor2
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_mxin_tmp(3)=nproc_mxin_tmp(3)*2
    else if(mod(icount,3)==2)then
      nproc_mxin_tmp(2)=nproc_mxin_tmp(2)*2
    else
      nproc_mxin_tmp(1)=nproc_mxin_tmp(1)*2
    end if
  end do

  nproc_k=1
  nproc_ob=1
  nproc_mxin(1:3)=nproc_mxin_tmp(1:3)
  nproc_mxin_s(1:3)=nproc_mxin_tmp(1:3)
  nproc_mxin_s_dm(1:3)=nproc_mxin_s(1:3)/nproc_mxin(1:3)
  
end subroutine set_numcpu_rt_sp
