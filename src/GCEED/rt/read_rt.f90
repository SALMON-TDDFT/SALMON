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
SUBROUTINE OUT_data_rt
use salmon_parallel, only: nproc_id_global, nproc_group_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use scf_data
use writebox_rt_sub
use allocate_mat_sub
implicit none
integer       :: i1,i2,i3,jj,iob,is,it2,kk
integer       :: ix,iy,iz
integer :: ibox
character(100) :: file_OUT_rt_data
integer :: ii,j1,j2,j3
integer :: nproc_id_global_datafiles
integer :: ista_Mxin_datafile(3)
integer :: iend_Mxin_datafile(3)
integer :: inum_Mxin_datafile(3)
integer :: nproc_xyz_datafile(3)
character(8) :: fileNumber_data
integer :: iob_myob
integer :: icorr_p

if(comm_is_root(nproc_id_global)) open(99,file=file_OUT_rt,form='unformatted')

if(comm_is_root(nproc_id_global)) write(99) itotNtime

if(num_datafiles_OUT>=2.and.num_datafiles_OUT<=nproc_size_global)then
  if(nproc_id_global<num_datafiles_OUT)then
    nproc_id_global_datafiles=nproc_id_global

    ibox=1
    nproc_xyz_datafile=1
    do ii=1,19
      do jj=3,1,-1
        if(ibox<num_datafiles_OUT)then
          nproc_xyz_datafile(jj)=nproc_xyz_datafile(jj)*2
          ibox=ibox*2
        end if
      end do
    end do

    do j3=0,nproc_xyz_datafile(3)-1
    do j2=0,nproc_xyz_datafile(2)-1
    do j1=0,nproc_xyz_datafile(1)-1
      ibox = j1 + nproc_xyz_datafile(1)*j2 + nproc_xyz_datafile(1)*nproc_xyz_datafile(2)*j3 
      if(ibox==nproc_id_global_datafiles)then
        ista_Mxin_datafile(1)=j1*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)
        iend_Mxin_datafile(1)=(j1+1)*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)-1
        ista_Mxin_datafile(2)=j2*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)
        iend_Mxin_datafile(2)=(j2+1)*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)-1
        ista_Mxin_datafile(3)=j3*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)
        iend_Mxin_datafile(3)=(j3+1)*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)-1
      end if
    end do
    end do
    end do
    inum_Mxin_datafile(:)=iend_Mxin_datafile(:)-ista_Mxin_datafile(:)+1

    write(fileNumber_data, '(i8)') nproc_id_global_datafiles
    file_Out_rt_data = trim(file_Out_rt)//"."//adjustl(fileNumber_data)
    open(89,file=file_Out_rt_data,form='unformatted')

  end if
end if

do iob=1,itotMST
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,icorr_p)

  cmatbox_l2=0.d0
    if(mod(itotNtime,2)==1)then
      if(icorr_p==1)then
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox_l2(ix,iy,iz)=zpsi_out(ix,iy,iz,iob_myob,1)
        end do
        end do
        end do
      end if
    else
      if(icorr_p==1)then
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox_l2(ix,iy,iz)=zpsi_in(ix,iy,iz,iob_myob,1)
        end do
        end do
        end do
      end if
    end if

  call comm_summation(cmatbox_l2,cmatbox_l,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

  if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc_size_global)then
    if(comm_is_root(nproc_id_global))then
      write(99) ((( cmatbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
  else
    if(nproc_id_global<num_datafiles_OUT)then
      write(89) ((( cmatbox_l(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                        iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                        iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
    end if
  end if
  
end do

if(num_datafiles_OUT>=2.and.num_datafiles_OUT<=nproc_size_global)then
  if(nproc_id_global<num_datafiles_OUT)then
    close(89)
  end if
end if

matbox_l=0.d0
do i3=ng_sta(3),ng_end(3)
do i2=ng_sta(2),ng_end(2)
do i1=ng_sta(1),ng_end(1)
   matbox_l(i1,i2,i3)=rho(i1,i2,i3)
end do
end do
end do

call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
if(comm_is_root(nproc_id_global))then
  write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(ilsda == 1)then
  do is=1,2

    matbox_l=0.d0
    do i3=ng_sta(3),ng_end(3)
    do i2=ng_sta(2),ng_end(2)
    do i1=ng_sta(1),ng_end(1)
      matbox_l(i1,i2,i3)=rho_s(i1,i2,i3,is)
    end do
    end do
    end do

    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
    if(comm_is_root(nproc_id_global))then
      write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
  end do
end if


matbox_l=0.d0
do i3=ng_sta(3),ng_end(3)
do i2=ng_sta(2),ng_end(2)
do i1=ng_sta(1),ng_end(1)
   matbox_l(i1,i2,i3)=Vh(i1,i2,i3)
end do
end do
end do
call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
if(comm_is_root(nproc_id_global))then
  write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(ilsda == 0)then
  matbox_l=0.d0
  do i3=ng_sta(3),ng_end(3)
  do i2=ng_sta(2),ng_end(2)
  do i1=ng_sta(1),ng_end(1)
    matbox_l(i1,i2,i3)=Vxc(i1,i2,i3)
  end do
  end do
  end do
  call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
  if(comm_is_root(nproc_id_global))then
    write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
  end if
else if(ilsda == 1)then
  do is=1,2

    matbox_l=0.d0
    do i3=ng_sta(3),ng_end(3)
    do i2=ng_sta(2),ng_end(2)
    do i1=ng_sta(1),ng_end(1)
      matbox_l(i1,i2,i3)=Vxc_s(i1,i2,i3,is)
    end do
    end do
    end do

    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
    if(comm_is_root(nproc_id_global))then
      write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
  end do
end if


matbox_l=0.d0
do i3=ng_sta(3),ng_end(3)
do i2=ng_sta(2),ng_end(2)
do i1=ng_sta(1),ng_end(1)
   matbox_l(i1,i2,i3)=Vh_stock1(i1,i2,i3)
end do
end do
end do
call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
if(comm_is_root(nproc_id_global))then
  write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if
matbox_l=0.d0
do i3=ng_sta(3),ng_end(3)
do i2=ng_sta(2),ng_end(2)
do i1=ng_sta(1),ng_end(1)
   matbox_l(i1,i2,i3)=Vh_stock2(i1,i2,i3)
end do
end do
end do
call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
if(comm_is_root(nproc_id_global))then
  write(99) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(comm_is_root(nproc_id_global)) then

  write(99) (vecDs(jj),jj=1,3)
  do it2=1,itotNtime
     write(99) (Dp(jj,it2),jj=1,3)
  end do
  if(iflag_dip2==1)then
    write(99) ((vecDs2(jj,kk),jj=1,3),kk=1,num_dip2)
    do it2=1,itotNtime
      write(99) ((Dp2(jj,it2,kk),jj=1,3),kk=1,num_dip2)
    end do
  end if
  
  close(99)
  
end if


END SUBROUTINE OUT_data_rt

!---------------------------------------------------------------------------
SUBROUTINE IN_data_rt(Ntime)
use salmon_parallel, only: nproc_id_global, nproc_group_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use scf_data
use new_world_sub
use readbox_rt_sub
use allocate_mat_sub
implicit none
integer       :: i1,i2,i3,jj,iob,is,it2,kk
integer       :: ix,iy,iz
integer       :: Ntime
character(100) :: file_IN_rt_data
integer :: ibox
integer :: ii,j1,j2,j3
integer :: nproc_id_global_datafiles
integer :: ista_Mxin_datafile(3)
integer :: iend_Mxin_datafile(3)
integer :: inum_Mxin_datafile(3)
integer :: nproc_xyz_datafile(3)
character(8) :: fileNumber_data
integer :: iob_myob
integer :: icorr_p

if(comm_is_root(nproc_id_global))then

   open(98,file=file_IN_rt,form='unformatted')

end if

if(comm_is_root(nproc_id_global)) read(98) Miter_rt
call comm_bcast(Miter_rt,nproc_group_global)
itotNtime=Ntime+Miter_rt

allocate(rIe(0:itotNtime))
allocate(Dp(3,0:itotNtime))
if(iflag_dip2==1) then
    allocate(rIe2(0:itotNtime,1:num_dip2))
    allocate(Dp2(3,0:itotNtime,1:num_dip2))
    allocate(Qp2(3,3,0:itotNtime,1:num_dip2))
end if
allocate(Qp(3,3,0:itotNtime))

if(num_datafiles_IN>=2.and.num_datafiles_IN<=nproc_size_global)then
  if(nproc_id_global<num_datafiles_IN)then
    nproc_id_global_datafiles=nproc_id_global

    ibox=1
    nproc_xyz_datafile=1
    do ii=1,19
      do jj=3,1,-1
        if(ibox<num_datafiles_IN)then
          nproc_xyz_datafile(jj)=nproc_xyz_datafile(jj)*2
          ibox=ibox*2
        end if
      end do
    end do

    do j3=0,nproc_xyz_datafile(3)-1
    do j2=0,nproc_xyz_datafile(2)-1
    do j1=0,nproc_xyz_datafile(1)-1
      ibox = j1 + nproc_xyz_datafile(1)*j2 + nproc_xyz_datafile(1)*nproc_xyz_datafile(2)*j3 
      if(ibox==nproc_id_global_datafiles)then
        ista_Mxin_datafile(1)=j1*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)
        iend_Mxin_datafile(1)=(j1+1)*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)-1
        ista_Mxin_datafile(2)=j2*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)
        iend_Mxin_datafile(2)=(j2+1)*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)-1
        ista_Mxin_datafile(3)=j3*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)
        iend_Mxin_datafile(3)=(j3+1)*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)-1
      end if
    end do
    end do
    end do
    inum_Mxin_datafile(:)=iend_Mxin_datafile(:)-ista_Mxin_datafile(:)+1

    write(fileNumber_data, '(i8)') nproc_id_global_datafiles
    file_IN_rt_data = trim(file_IN_rt)//"."//adjustl(fileNumber_data)
    open(88,file=file_IN_rt_data,form='unformatted')

  end if
end if

do iob=1,itotMST
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,icorr_p)
  
  if(num_datafiles_IN==1.or.num_datafiles_IN>nproc_size_global)then
    if(comm_is_root(nproc_id_global))then
      read(98) ((( cmatbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
    call comm_bcast(cmatbox_l(lg_sta(1):lg_sta(1)+lg_num(1),lg_sta(2):lg_sta(2)+lg_num(2),lg_sta(3):lg_sta(3)+lg_num(3)),  &
                    nproc_group_global)
  else
    cmatbox_l2=0.d0
    if(nproc_id_global<num_datafiles_IN)then
      read(88) ((( cmatbox_l2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                        iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                        iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
    end if
    call comm_summation(cmatbox_l2,cmatbox_l,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
  end if
  if(mod(Miter_rt,2)==1)then
    if(icorr_p==1)then
      zpsi_out(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
          mg_sta(3):mg_end(3),iob_myob,1)=  &
      cmatbox_l(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
             mg_sta(3):mg_end(3))
    end if
  else
    if(icorr_p==1)then
      zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
          mg_sta(3):mg_end(3),iob_myob,1)=  &
      cmatbox_l(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
             mg_sta(3):mg_end(3))
    end if
  end if
  
end do

if(num_datafiles_IN>=2.and.num_datafiles_IN<=nproc_size_global)then
  if(nproc_id_global<num_datafiles_IN)then
    close(88)
  end if
end if

if(comm_is_root(nproc_id_global))then
  read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if
call comm_bcast(matbox_l,nproc_group_global)
do i3=mg_sta(3),mg_end(3)
do i2=mg_sta(2),mg_end(2)
do i1=mg_sta(1),mg_end(1)
  rho(i1,i2,i3)=matbox_l(i1,i2,i3)
end do
end do
end do
if(ilsda == 1)then
  do is=1,2
    if(comm_is_root(nproc_id_global))then
      read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
    call comm_bcast(matbox_l,nproc_group_global)
    do i3=lg_sta(3),lg_end(3)
    do i2=lg_sta(2),lg_end(2)
    do i1=lg_sta(1),lg_end(1)
      rho_s(i1,i2,i3,is)=matbox_l(i1,i2,i3)
    end do
    end do
    end do
  end do
end if

if(comm_is_root(nproc_id_global))then
  read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if
call comm_bcast(matbox_l,nproc_group_global)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Vh(ix,iy,iz)=matbox_l(ix,iy,iz)
end do
end do
end do

if(ilsda == 0)then
  if(comm_is_root(nproc_id_global))then
    read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
  end if
  call comm_bcast(matbox_l,nproc_group_global)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    Vxc(ix,iy,iz)=matbox_l(ix,iy,iz)
  end do
  end do
  end do
else if(ilsda == 1)then
  do is=1,2
    if(comm_is_root(nproc_id_global))then
      read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
    call comm_bcast(matbox_l,nproc_group_global)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vxc_s(ix,iy,iz,is)=matbox_l(ix,iy,iz)
    end do
    end do
    end do

  end do
end if

if(comm_is_root(nproc_id_global))then
  read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if
call comm_bcast(matbox_l,nproc_group_global)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Vh_stock1(ix,iy,iz)=matbox_l(ix,iy,iz)
end do
end do
end do


if(comm_is_root(nproc_id_global))then
  read(98) ((( matbox_l(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if
call comm_bcast(matbox_l,nproc_group_global)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Vh_stock2(ix,iy,iz)=matbox_l(ix,iy,iz)
end do
end do
end do

itt=Miter_rt
call allgatherv_vlocal

if(comm_is_root(nproc_id_global))then
  read(98) (vecDs(jj),jj=1,3)
  do it2=1,Miter_rt
    read(98) (Dp(jj,it2),jj=1,3)
  end do
end if

call comm_bcast(vecDs,nproc_group_global)

do jj=1,3
do it2=1,Miter_rt
   call comm_bcast(Dp(jj,it2),nproc_group_global)
end do
end do

if(iflag_dip2==1)then
  if(comm_is_root(nproc_id_global))then
    read(98) ((vecDs2(jj,kk),jj=1,3),kk=1,num_dip2)
    do it2=1,Miter_rt
      read(98) ((Dp2(jj,it2,kk),jj=1,3),kk=1,num_dip2)
    end do
  end if
  call comm_bcast(vecDs2,nproc_group_global)
  do kk=1,num_dip2
    call comm_bcast(Dp2(:,:,kk),nproc_group_global)
  end do
end if

if(comm_is_root(nproc_id_global)) close(98)

END SUBROUTINE IN_data_rt

