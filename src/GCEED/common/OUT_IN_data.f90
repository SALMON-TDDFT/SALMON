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
!=======================================================================

SUBROUTINE OUT_data
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global, nproc_group_h
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use scf_data
use new_world_sub
use read_pslfile_sub
use allocate_psl_sub
use allocate_mat_sub
implicit none
integer :: is,iob,jj,ik
integer :: ix,iy,iz
real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
character(100) :: file_OUT_data
character(100) :: file_OUT_data_ini
integer :: ibox
integer :: ii,j1,j2,j3
integer :: myrank_datafiles
integer :: ista_Mxin_datafile(3)
integer :: iend_Mxin_datafile(3)
integer :: inum_Mxin_datafile(3)
integer :: nproc_xyz_datafile(3)
character(8) :: fileNumber_data
integer :: iob_myob
integer :: icorr_p

if(comm_is_root(nproc_id_global))then

  open(97,file=file_OUT,form='unformatted')
  
!version number
  version_num(1)=40
  version_num(2)=1
  write(97) version_num(1),version_num(2)
  write(97) Nd
  write(97) ilsda
  write(97) iflag_ps
  write(97) iend_Mx_ori(:3)
  write(97) lg_end(:3)
  if(ilsda == 0)then
    write(97) MST(1)
    write(97) ifMST(1)
  else if(ilsda == 1)then
    write(97) (MST(is),is=1,2)
    write(97) (ifMST(is),is=1,2)
  end if
  if(iflag_ps.eq.1)then
    write(97) MI,MKI,maxMps,Mlmps
  end if
  write(97) (Hgs(jj),jj=1,3)
  write(97) (rLsize(jj,ntmg),jj=1,3)
  write(97) Miter
  write(97) MEO
  
  if(iflag_ps.eq.1)then
    write(97) Jxyz_all(1:3,1:maxMps,1:MI),Mps_all(1:MI)
  end if
  
  if(iflag_ps.eq.1)then
    write(97) Kion(:MI)
    write(97) Rion(:,:MI)
    write(97) iZatom(:MKI)
    write(97) pseudo_file(:MKI) !ipsfileform(:MKI)
    write(97) Zps(:MKI),Rps(:MKI)
    write(97) AtomName(:MI) 
    write(97) iAtomicNumber(:MI) 
  end if
  
end if

if(comm_is_root(nproc_id_global))then
  if(iflag_ps.eq.1)then
    write(97) uV_all(:maxMps,:Mlmps,:MI),uVu(:Mlmps,:MI)
    write(97) Mlps(:MKI),Lref(:MKI)
  end if
end if

allocate(matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(matbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(cmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))

if(OC<=2)then
  if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc_size_global)then
    file_OUT_data_ini = file_OUT_ini
  else
    if(nproc_id_global<num_datafiles_OUT)then
      myrank_datafiles=nproc_id_global

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
        if(ibox==myrank_datafiles)then
          ista_Mxin_datafile(1)=j1*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)
          iend_Mxin_datafile(1)=(j1+1)*lg_num(1)/nproc_xyz_datafile(1)+lg_sta(1)-1
          ista_Mxin_datafile(2)=j2*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)
          iend_Mxin_datafile(2)=(j2+1)*lg_num(2)/nproc_xyz_datafile(2)+lg_sta(2)-1
          ista_Mxin_datafile(3)=j3*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)
          iend_Mxin_datafile(3)=(j3+1)*lg_num(3)/nproc_xyz_datafile(3)+lg_sta(3)-1
          if(OC==2)then
            mg_sta_ini(1)=j1*lg_num_ini(1)/nproc_xyz_datafile(1)+lg_sta_ini(1)
            mg_end_ini(1)=(j1+1)*lg_num_ini(1)/nproc_xyz_datafile(1)+lg_sta_ini(1)-1
            mg_sta_ini(2)=j2*lg_num_ini(2)/nproc_xyz_datafile(2)+lg_sta_ini(2)
            mg_end_ini(2)=(j2+1)*lg_num_ini(2)/nproc_xyz_datafile(2)+lg_sta_ini(2)-1
            mg_sta_ini(3)=j3*lg_num_ini(3)/nproc_xyz_datafile(3)+lg_sta_ini(3)
            mg_end_ini(3)=(j3+1)*lg_num_ini(3)/nproc_xyz_datafile(3)+lg_sta_ini(3)-1
          end if
        end if
      end do
      end do
      end do
      inum_Mxin_datafile(:)=iend_Mxin_datafile(:)-ista_Mxin_datafile(:)+1

      write(fileNumber_data, '(i6.6)') myrank_datafiles
      file_OUT_data = trim(adjustl(sysname))//"_gs_"//trim(adjustl(fileNumber_data))//".bin"
      open(87,file=file_OUT_data,form='unformatted')
    end if
  end if
end if

select case(iperiodic)
case(0)
  if(OC<=2)then
  
    do ik=1,num_kpoints_rd
    do iob=1,itotMST
      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,ik,icorr_p)
  
      matbox_l=0.d0
      if(icorr_p==1)then
        matbox_l(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))   &
          = psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iob_myob,ik)
      end if
  
      call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
  
  
      if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc_size_global)then
        if(comm_is_root(nproc_id_global))then
          write(97) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
        end if
      else
        if(nproc_id_global<num_datafiles_OUT)then
          write(87) ((( matbox_l2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                            iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                            iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
        end if
      end if
    end do
    end do
  
  else if(OC==3)then
    do iob=1,iobnum
      write(87,rec=iob) ((( psi(ix,iy,iz,iob,1),ix=mg_sta(1),mg_end(1)),   &
                                                iy=mg_sta(2),mg_end(2)),   &
                                                iz=mg_sta(3),mg_end(3))
    end do
  end if
case(3)
  if(OC<=2)then

    do ik=1,num_kpoints_rd
    do iob=1,itotMST

      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,ik,icorr_p)
      cmatbox_l=0.d0
      if(icorr_p==1)then
        cmatbox_l(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))   &
            = zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iob_myob,ik)
      end if

      call comm_summation(cmatbox_l,cmatbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
      if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc_size_global)then
      if(comm_is_root(nproc_id_global))then
        write(97) ((( cmatbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
      end if
      else
        if(nproc_id_global<num_datafiles_OUT)then
          write(87) ((( cmatbox_l2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                             iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                             iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
        end if
      end if
    end do
    end do

  else if(OC==3)then
    do ik=k_sta,k_end
    do iob=1,iobnum
      write(87,rec=iob) ((( zpsi(ix,iy,iz,iob,ik),ix=mg_sta(1),mg_end(1)),   &
                                                  iy=mg_sta(2),mg_end(2)),   &
                                                  iz=mg_sta(3),mg_end(3))
    end do
    end do
  end if
end select



if(OC<=2)then
  if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc_size_global)then
    if(comm_is_root(nproc_id_global).and.OC==2) close(67)
  else
    if(nproc_id_global<num_datafiles_OUT)then
      close(87)
      if(OC==2) close(67)
    end if
  end if
else if(OC==3)then
  close(87)
end if

matbox2=0.d0
matbox2(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))   &
   = rho(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))

call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

if(comm_is_root(nproc_id_global))then
  write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

do ii=1,num_rho_stock+1
  matbox2=0.d0
  matbox2(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))   &
     = rho_in(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3),ii)

  call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

  if(comm_is_root(nproc_id_global))then
    write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
  end if
end do

do ii=1,num_rho_stock
  matbox2=0.d0
  matbox2(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))   &
     = rho_out(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3),ii)

  call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)
  if(comm_is_root(nproc_id_global))then
    write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
  end if
end do

if(ilsda == 1)then
  do is=1,2
    matbox2=0.d0
    matbox2(ng_sta(1):ng_end(1),   &
            ng_sta(2):ng_end(2),   &
            ng_sta(3):ng_end(3))   &
      = rho_s(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3),is)

    call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

    if(comm_is_root(nproc_id_global))then
      write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if

    do ii=1,num_rho_stock+1
      matbox2=0.d0
      matbox2(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3))   &
        = rho_s_in(ng_sta(1):ng_end(1),   &
                ng_sta(2):ng_end(2),   &
                ng_sta(3):ng_end(3),ii,is)

      call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

      if(comm_is_root(nproc_id_global))then
        write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
      end if
    end do

    do ii=1,num_rho_stock
      matbox2=0.d0
      matbox2(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3))   &
        = rho_s_out(ng_sta(1):ng_end(1),   &
                ng_sta(2):ng_end(2),   &
                ng_sta(3):ng_end(3),ii,is)

      call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

      if(comm_is_root(nproc_id_global))then
        write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
      end if
    end do
    
  end do

end if

if(comm_is_root(nproc_id_global))then
  write(97) esp(:itotMST,:num_kpoints_rd),rocc(:itotMST,:num_kpoints_rd)
end if

matbox2=0.d0
matbox2(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))   &
   = Vh(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))

call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

if(comm_is_root(nproc_id_global))then
  write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(ilsda == 0)then
  matbox2=0.d0
  matbox2(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))   &
     = Vxc(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))

  call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

  if(comm_is_root(nproc_id_global))then
    write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
  end if
else if(ilsda == 1) then
  do is=1,2
    matbox2=0.d0
    matbox2(ng_sta(1):ng_end(1),   &
            ng_sta(2):ng_end(2),   &
            ng_sta(3):ng_end(3))   &
     = Vxc_s(ng_sta(1):ng_end(1),   &
            ng_sta(2):ng_end(2),   &
            ng_sta(3):ng_end(3),is)

    call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

    if(comm_is_root(nproc_id_global))then
      write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if
  end do
end if


matbox2=0.d0
matbox2(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))   &
   = Vpsl(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))

  call comm_summation(matbox2,matbox,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

if(comm_is_root(nproc_id_global))then
  write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(comm_is_root(nproc_id_global))then
  close(97)
end if

deallocate(matbox,matbox2)
deallocate(cmatbox,cmatbox2)

END SUBROUTINE OUT_data

!=======================================================================

SUBROUTINE IN_data
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global, nproc_id_spin
use salmon_parallel, only: nproc_id_orbitalgrid, nproc_id_kgrid
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
integer :: NI0,Ndv0,Nps0,Nd0
integer :: ii,is,iob,jj,ibox,j1,j2,j3,ik
integer :: ix,iy,iz
real(8),allocatable :: matbox(:,:,:)
real(8),allocatable :: matbox2(:,:,:)
real(8),allocatable :: matbox3(:,:,:)
real(8),allocatable :: esp0(:,:),rocc0(:,:)
complex(8),allocatable :: cmatbox(:,:,:)
complex(8),allocatable :: cmatbox2(:,:,:)
character(100) :: file_IN_data
character(8) :: cha_version_num(2)
integer :: version_num_box(2)
integer :: myrank_datafiles
integer :: ista_Mxin_datafile(3)
integer :: iend_Mxin_datafile(3)
integer :: inum_Mxin_datafile(3)
integer :: nproc_xyz_datafile(3)
character(8) :: fileNumber_data
integer :: maxMdvbox
integer :: iob_myob
integer :: icheck_corrkob
integer :: pstart(2),pend(2)
integer :: is_sta,is_end
integer :: p0
integer :: iobnum0
integer :: icount
complex(8),parameter :: zi=(0.d0,1.d0)

integer :: ig_sta(3),ig_end(3),ig_num(3)
real(8),allocatable :: matbox_read(:,:,:)
real(8),allocatable :: matbox_read2(:,:,:)
complex(8),allocatable :: cmatbox_read(:,:,:)
complex(8),allocatable :: cmatbox_read2(:,:,:)
real(8),allocatable :: matbox_read3(:,:,:)
complex(8),allocatable :: cmatbox_read3(:,:,:)
integer :: icheck_read
integer :: ifilenum_data
integer :: icomm
integer :: ifMST0(2)
integer :: imesh_oddeven0
integer :: itmg

if(comm_is_root(nproc_id_global))then
  write(*,*) file_IN
  open(96,file=file_IN,form='unformatted')

  read(96) version_num_box(1),version_num_box(2)
end if

call comm_bcast(version_num_box,nproc_group_global)

if(version_num_box(1)>=40)then
  continue
else if((version_num_box(1)==17.and.version_num_box(2)>=13).or.version_num_box(1)>=18) then
  if(comm_is_root(nproc_id_global))then
    read(96) imesh_oddeven0
  end if
  call comm_bcast(imesh_oddeven0,nproc_group_global)
else
  continue
end if

if(comm_is_root(nproc_id_global)) then
   read(96) Nd0
   read(96) ilsda
   if(version_num_box(1)<=36)then
     read(96) iflag_ps,ibox
   else
     read(96) iflag_ps
   end if
   if(version_num_box(1)==17.and.version_num_box(2)<=10)then
     read(96) NI0,Ndv0,Nps0,Nd0
   end if

   write(cha_version_num(1), '(i8)') version_num_box(1)
   write(cha_version_num(2), '(i8)') version_num_box(2)
   if((version_num_box(1)==17.and.version_num_box(2)==22).or.  &
      (version_num_box(1)==18.and.version_num_box(2)==17).or.  &
      (version_num_box(1)==23.and.version_num_box(2)==62).or.  &
      (version_num_box(1)==25.and.version_num_box(2)==17).or.  &
      (version_num_box(1)==26.and.version_num_box(2)==3).or.  &
      (version_num_box(1)==27.and.version_num_box(2)==9).or.  &
      (version_num_box(1)==28.and.version_num_box(2)==1).or.  &
      (version_num_box(1)==29.and.version_num_box(2)==1))then
     write(*,'(a,a)') "A version of input data file is ", &
      "1."//trim(adjustl(cha_version_num(1)))
   else if(version_num_box(1)==30.and.version_num_box(2)>=18)then
     write(*,'(a,a)') "A version of input data file is ", &
      trim(adjustl(cha_version_num(2)))
   else
     write(*,'(a,a)') "A version of input data file is ", &
      "1."//trim(adjustl(cha_version_num(1)))//"."//trim(adjustl(cha_version_num(2)))
   end if
end if

call comm_bcast(ilsda,nproc_group_global)
call comm_bcast(iflag_ps,nproc_group_global)

if(comm_is_root(nproc_id_global))then
  read(96) iend_Mx_ori(:3)
  read(96) lg_end(:3)
  if(ilsda == 0) then
    read(96) MST0(1)
!    read(96) ifMST(1)
    if(iSCFRT==2)then
      read(96) ifMST(1)
    else
      read(96) ifMST0(1)
    endif
  else if(ilsda == 1)then
    read(96) (MST0(is),is=1,2)
!    read(96) (ifMST(is),is=1,2)
    if(iSCFRT==2)then
      read(96) (ifMST(is),is=1,2)
    else
      read(96) (ifMST0(is),is=1,2)
    endif
  end if
  if(version_num_box(1)<=31)then
    if(iflag_ps.eq.1)then
      read(96) MI_read,MKI,maxMdvbox,maxMps,Mlmps
    end if
  else
    if(iflag_ps.eq.1)then
      read(96) MI_read,MKI,maxMps,Mlmps
    end if
  end if
  if(version_num_box(1)>=35)then
    read(96) (Hgs(jj),jj=1,3)
  else
    read(96) Hgs(1)
    Hgs(2)=Hgs(1)
    Hgs(3)=Hgs(1)
  end if
  Hvol=Hgs(1)*Hgs(2)*Hgs(3)
  read(96) (rLsize(jj,1),jj=1,3)
  read(96) Miter
  read(96) ibox
end if

call comm_bcast(iend_Mx_ori,nproc_group_global)
call comm_bcast(lg_end,nproc_group_global)
call comm_bcast(MST0,nproc_group_global)
call comm_bcast(ifMST,nproc_group_global)
call comm_bcast(Hgs,nproc_group_global)
call comm_bcast(Hvol,nproc_group_global)
call comm_bcast(rLsize,nproc_group_global)
call comm_bcast(Miter,nproc_group_global)

itmg=1
call set_imesh_oddeven(itmg)

if(version_num_box(1)>=40)then
  continue
else if((version_num_box(1)==17.and.version_num_box(2)>=13).or.version_num_box(1)>=18) then
  if(imesh_oddeven0==1.and.imesh_oddeven(1)==1.and.imesh_oddeven(2)==1.and.imesh_oddeven(3)==1)then
    continue
  else if(imesh_oddeven0==2.and.imesh_oddeven(1)==2.and.imesh_oddeven(2)==2.and.imesh_oddeven(3)==2)then
    continue
  else
    stop "You cannot use data files of this version because imesh_oddeven and values of Lsize/Hgs are a mixture of odd and even."
  end if
else
  if(imesh_oddeven(1)==2.and.imesh_oddeven(2)==2.and.imesh_oddeven(3)==2)then
    continue
  else
    stop "You cannot use data files of this version because imesh_oddeven is not 2."
  end if
end if

if(iSCFRT==2) then
  if(ilsda == 0) then
    MST(1)=ifMST(1)
  else if(ilsda == 1) then
    MST(1:2)=ifMST(1:2)
  end if
end if

select case(iperiodic)
case(0)
  do jj=1,3
    select case(imesh_oddeven(jj))
      case(1)
        ista_Mx_ori(jj)=-iend_Mx_ori(jj)
        lg_sta(jj)=-lg_end(jj)
      case(2)
        ista_Mx_ori(jj)=-iend_Mx_ori(jj)+1
        lg_sta(jj)=-lg_end(jj)+1
    end select
  end do
case(3)
  ista_Mx_ori(:)=1-Nd
  lg_sta(:)=1
end select
inum_Mx_ori(:)=iend_Mx_ori(:)-ista_Mx_ori(:)+1

lg_num(:)=lg_end(:)-lg_sta(:)+1

call check_fourier

call set_gridcoo

allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1))
allocate(inum_Mxin(3,0:nproc_size_global-1))

call setmg(mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential)

if(ilsda == 0) then
  itotMST0=MST0(1)
  itotMST=MST(1)
  itotfMST=ifMST(1)
else if(ilsda == 1) then
  itotMST0=MST0(1)+MST0(2)
  itotMST=MST(1)+MST(2)
  itotfMST=ifMST(1)+ifMST(2)
end if

call init_mesh_s
call check_ng

if(iflag_ps.eq.1)then
  call comm_bcast(MI_read,nproc_group_global)
  call comm_bcast(MKI,nproc_group_global)
  call comm_bcast(maxMps,nproc_group_global)
  call comm_bcast(Mlmps,nproc_group_global)
  MI=MI_read
end if



if(iflag_ps.eq.1)then
   if(comm_is_root(nproc_id_global))then
     if(version_num_box(1)<=31)then
       read(96) 
       read(96) 
     else
       read(96) 
     end if
   end if

  if(iSCFRT==2) then
!    allocate( Kion(MI),Rion(3,MI) )
  end if
  if(iSCFRT==2) allocate( AtomName(MI), iAtomicNumber(MI) )
  if(comm_is_root(nproc_id_global))then
    read(96) Kion(:MI_read)
    read(96) Rion(:,:MI_read)
    read(96) iZatom(:MKI)
    if(version_num_box(1)>=34)then
      read(96) pseudo_file(:MKI) !ipsfileform(:MKI)
    else
      stop "This version is already invalid."
    end if
    read(96) 
    read(96) AtomName(:MI_read)
    read(96) iAtomicNumber(:MI_read)
  end if
  
  call comm_bcast(Kion,nproc_group_global)
  call comm_bcast(Rion,nproc_group_global)
  call comm_bcast(iZatom,nproc_group_global)
  call comm_bcast(pseudo_file,nproc_group_global)
  call comm_bcast(AtomName,nproc_group_global)
  call comm_bcast(iAtomicNumber,nproc_group_global)

end if

if(ilsda==1)then
  nproc_ob_spin(1)=(nproc_ob+1)/2
  nproc_ob_spin(2)=nproc_ob/2
end if

if(iSCFRT==2) call make_new_world

call setk(k_sta, k_end, k_num, num_kpoints_rd, nproc_k, nproc_id_orbitalgrid)
if(ilsda==0)then
  call calc_iobnum(itotMST,nproc_ob,nproc_id_kgrid,iobnum,nproc_ob,iparaway_ob)
else if(ilsda==1)then
  if(nproc_ob==1)then
    iobnum=itotMST
  else
    if(nproc_id_spin<nproc_ob_spin(1))then
      call calc_iobnum(MST(1),nproc_ob_spin(1),nproc_id_kgrid,iobnum,nproc_ob_spin(1),iparaway_ob)
    else
      call calc_iobnum(MST(2),nproc_ob_spin(2),nproc_id_kgrid,iobnum,nproc_ob_spin(2),iparaway_ob)
    end if
  end if
end if

if(iSCFRT==2)then
  call allocate_mat
  call set_icoo1d
end if

allocate(k_rd0(3,num_kpoints_rd),ksquare0(num_kpoints_rd))
if(iperiodic==3)then
  call init_k_rd(k_rd0,ksquare0,3)
end if

allocate( matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
allocate( cmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
allocate( matbox3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )

if(iSCFRT==1)then
  select case(iperiodic)
  case(0)
    if(iobnum.ge.1)then
      allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                    mg_sta(3):mg_end(3), &
&                   1:iobnum,k_sta:k_end) )
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      allocate( psi_mesh(ng_sta(1):ng_end(1),  &
                       ng_sta(2):ng_end(2),   &
                       ng_sta(3):ng_end(3), &
                       1:itotMST,num_kpoints_rd) )
    end if
  case(3)
    allocate( ttpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    if(iobnum.ge.1)then
      allocate( zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3), &
      &                    1:iobnum,k_sta:k_end) )
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      allocate( zpsi_mesh(ng_sta(1):ng_end(1),  &
                          ng_sta(2):ng_end(2),   &
                          ng_sta(3):ng_end(3), &
                          1:itotMST,num_kpoints_rd) )
    end if
  end select
else if(iSCFRT==2)then
  if(iobnum.ge.1)then
    allocate( zpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum,k_sta:k_end) )
    allocate( zpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum,k_sta:k_end) )
    zpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                    1:iobnum,k_sta:k_end) = 0.d0
    zpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                    1:iobnum,k_sta:k_end) = 0.d0
  end if
  if(iwrite_projection==1)then
    if(ilsda==0)then
      call calc_iobnum(itotMST0,nproc_ob,nproc_id_kgrid,iobnum0,nproc_ob,iparaway_ob)
    else if(ilsda==1)then
      if(nproc_ob==1)then
        iobnum0=itotMST0
      else
        if(nproc_id_spin<nproc_ob_spin(1))then
          call calc_iobnum(MST0(1),nproc_ob_spin(1),nproc_id_kgrid,iobnum0,nproc_ob_spin(1),iparaway_ob)
        else
          call calc_iobnum(MST0(2),nproc_ob_spin(2),nproc_id_kgrid,iobnum0,nproc_ob_spin(2),iparaway_ob)
        end if
      end if
    end if
    if(iobnum0.ge.1)then
      allocate( zpsi_t0(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum0,k_sta:k_end) )
      zpsi_t0(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum0,k_sta:k_end) = 0.d0
    end if
  end if
end if

allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate( rho0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate( rho_diff(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
if(iSCFRT==1)then
  allocate( rho_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1))
  allocate( rho_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1))
  rho_in=0.d0
  rho_out=0.d0
end if

if(ilsda == 0) then
  continue
else if(ilsda == 1) then
  allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2))
  allocate( rho_s_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2))
  allocate( rho_s_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2))
  rho_s_in=0.d0
  rho_s_out=0.d0
end if

if(iSCFRT==1)then
  allocate( esp0(itotMST0,num_kpoints_rd))
  allocate( esp(itotMST,num_kpoints_rd))
  allocate( rocc0(itotMST0,num_kpoints_rd))
else if(iSCFRT==2)then
  allocate( esp(itotMST,num_kpoints_rd),rocc(itotMST,num_kpoints_rd))
  allocate( esp0(itotMST0,num_kpoints_rd),rocc0(itotMST0,num_kpoints_rd))
  allocate( esp2(itotMST,num_kpoints_rd))
  esp2=0.d0
end if
allocate( Vh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
if(ilsda == 0) then
  allocate( Vxc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
else if(ilsda == 1) then
  allocate( Vxc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
end if
if(ilsda==0)then
  allocate( Vlocal(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1) )
  allocate( Vlocal2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1) )
else
  allocate( Vlocal(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
  allocate( Vlocal2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
end if
allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
if(icalcforce==1) allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )

if(comm_is_root(nproc_id_global))then
  if(version_num_box(1)>=32)then
    if(iflag_ps.eq.1)then
      read(96) 
      read(96) Mlps(:MKI),Lref(:MKI)
    end if
  end if
end if
if(version_num_box(1)>=32)then
  if(iflag_ps.eq.1)then
    call comm_bcast(Mlps,nproc_group_global)
    call comm_bcast(Lref,nproc_group_global)
  end if
end if 

allocate( cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )

if(num_datafiles_IN==1)then
  ifilenum_data=96
else
  ifilenum_data=86
end if

!set ista_Mxin_datafile etc.
if(IC<=2)then
  if(num_datafiles_IN<=nproc_size_global)then
    if(nproc_id_global<num_datafiles_IN)then
      myrank_datafiles=nproc_id_global

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
        if(ibox==myrank_datafiles)then
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

      if(num_datafiles_IN>=2.and.nproc_id_global<num_datafiles_IN)then
        write(fileNumber_data, '(i6.6)') myrank_datafiles
        file_IN_data = trim(adjustl(sysname))//"_gs_"//trim(adjustl(fileNumber_data))//".bin"
        open(86,file=file_IN_data,form='unformatted')
      end if
    end if
  end if
end if

if(ilsda == 0)then
  is_sta=1
  is_end=1
  pstart(1)=1
  pend(1)=itotMST0
else if(ilsda == 1)then
  is_sta=1
  is_end=2
  pstart(1)=1
  pend(1)=MST0(1)
  pstart(2)=MST0(1)+1
  pend(2)=itotMST0
end if

allocate( matbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )

ig_sta(:)=lg_sta(:)
ig_end(:)=lg_end(:)
ig_num(:)=lg_num(:)

allocate( matbox_read(ig_sta(1):ig_end(1),ig_sta(2):ig_end(2),ig_sta(3):ig_end(3)) )
allocate( matbox_read2(ig_sta(1):ig_end(1),ig_sta(2):ig_end(2),ig_sta(3):ig_end(3)) )
allocate( cmatbox_read(ig_sta(1):ig_end(1),ig_sta(2):ig_end(2),ig_sta(3):ig_end(3)) )
allocate( cmatbox_read2(ig_sta(1):ig_end(1),ig_sta(2):ig_end(2),ig_sta(3):ig_end(3)) )

allocate( matbox_read3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
allocate( cmatbox_read3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )

!$OMP parallel do private(iz,iy,ix) 
do iz=ig_sta(3),ig_end(3)
do iy=ig_sta(2),ig_end(2)
do ix=ig_sta(1),ig_end(1)
  matbox_read2(ix,iy,iz)=0.d0
  cmatbox_read2(ix,iy,iz)=0.d0
end do
end do
end do

icount=0

do ik=1,num_kpoints_rd
do is=is_sta,is_end
do p0=pstart(is),pend(is)

! read file
  call conv_p0(p0,iob)
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,ik,icheck_corrkob)

  if(IC<=2)then
    if(nproc_id_global<num_datafiles_IN)then
      icheck_read=1
    else
      icheck_read=0
    end if
  else if(IC==3.or.IC==4)then
    if(icheck_corrkob==1)then
      icheck_read=1
    else
      icheck_read=0
    end if
  end if
  
  if(icheck_read==1)then
    icount=icount+1
    if(IC<=2)then
      select case(iperiodic)
      case(0)
        read(ifilenum_data) ((( matbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                       iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                       iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
      case(3)
        read(ifilenum_data) ((( cmatbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                        iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                        iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
      end select
    else if(IC==3.or.IC==4)then
      select case(iperiodic)
      case(0)
        read(ifilenum_data,rec=icount) ((( matbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                     iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                     iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
      case(3)
        read(ifilenum_data,rec=icount) ((( cmatbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                      iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                      iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
      end select
    end if
  end if
  
  icomm=nproc_group_global

  select case(iperiodic)
  case(0)
    call comm_summation(matbox_read2,matbox_read,ig_num(1)*ig_num(2)*ig_num(3),icomm)
  case(3)
    call comm_summation(cmatbox_read2,cmatbox_read,ig_num(1)*ig_num(2)*ig_num(3),icomm)
  end select

  if(icheck_corrkob==1)then
    if(iSCFRT==1)then
      select case(iperiodic)
      case(0)
        psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
            mg_sta(3):mg_end(3),iob_myob,ik)=  &
        matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
               mg_sta(3):mg_end(3))
      case(3)
        zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3),iob_myob,ik)=  &
        cmatbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                mg_sta(3):mg_end(3))
      end select
    else if(iSCFRT==2)then
      select case(iperiodic)
      case(0)
        if(iwrite_projection==1)then
          zpsi_t0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,ik)=  &
          matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                  mg_sta(3):mg_end(3))
        else
          if((ilsda==0.and.p0<=MST(1)).or.  &
             (ilsda==1.and.(p0<=MST0(1).and.p0<=MST(1)).or.(p0>MST0(1).and.p0<=MST0(1)+MST(2))))then
            zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                    mg_sta(3):mg_end(3),iob_myob,ik)=  &
            matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                        mg_sta(3):mg_end(3))
          end if
        end if
      case(3)
        if(iwrite_projection==1)then
          zpsi_t0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,ik)=  &
             cmatbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                          mg_sta(3):mg_end(3))
        else
          if((ilsda==0.and.p0<=MST(1)).or.  &
             (ilsda==1.and.(p0<=MST0(1).and.p0<=MST(1)).or.(p0>MST0(1).and.p0<=MST0(1)+MST(2))))then
            zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                    mg_sta(3):mg_end(3),iob_myob,ik)=  &
              cmatbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                           mg_sta(3):mg_end(3))
          end if
        end if
      end select
      if(iwrite_projection==1)then
        if((ilsda==0.and.p0<=MST(1)).or.  &
           (ilsda==1.and.(p0<=MST0(1).and.p0<=MST(1)).or.(p0>MST0(1).and.p0<=MST0(1)+MST(2))))then
          zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,ik)=  &
          zpsi_t0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,ik)
        end if
      end if
    end if
  end if
  
end do
end do
end do

if(iSCFRT==1.and.itotMST>itotMST0) call init_wf_ns(2)

if(IC<=2)then
  call read_copy_pot(rho,matbox_read,ig_sta,ig_end)
 
  if(version_num_box(1)<=29.or.(version_num_box(1)==30.and.version_num_box(2)<=18))then
    if(comm_is_root(nproc_id_global))then
      read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
    end if
    if(iSCFRT==1)then
      call comm_bcast(matbox_read,nproc_group_global)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        rho_in(ix,iy,iz,num_rho_stock+1)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end if
  
    if(comm_is_root(nproc_id_global))then
      read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
    end if
    if(iSCFRT==1)then
      call comm_bcast(matbox_read,nproc_group_global)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        rho_out(ix,iy,iz,num_rho_stock)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end if
  
    if(ilsda == 1)then
      do is=1,2
        if(comm_is_root(nproc_id_global))then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        if(comm_is_root(nproc_id_global))then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rho_s_in(ix,iy,iz,num_rho_stock,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        if(comm_is_root(nproc_id_global))then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rho_s_out(ix,iy,iz,num_rho_stock,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
      end do
    end if
  else
    do ii=1,num_rho_stock+1
      if(comm_is_root(nproc_id_global))then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      if(iSCFRT==1)then
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rho_in(ix,iy,iz,ii)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  
    do ii=1,num_rho_stock
      if(comm_is_root(nproc_id_global))then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      if(iSCFRT==1)then
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rho_out(ix,iy,iz,ii)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  
    if(ilsda == 1)then
      do is=1,2
        if(comm_is_root(nproc_id_global))then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call comm_bcast(matbox_read,nproc_group_global)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        do ii=1,num_rho_stock+1
          if(comm_is_root(nproc_id_global))then
            read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
          end if
          if(iSCFRT==1)then
            call comm_bcast(matbox_read,nproc_group_global)
            do iz=ng_sta(3),ng_end(3)
            do iy=ng_sta(2),ng_end(2)
            do ix=ng_sta(1),ng_end(1)
              rho_s_in(ix,iy,iz,ii,is)=matbox_read(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
    
        do ii=1,num_rho_stock
          if(comm_is_root(nproc_id_global))then
            read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
          end if
          if(iSCFRT==1)then
            call comm_bcast(matbox_read,nproc_group_global)
            do iz=ng_sta(3),ng_end(3)
            do iy=ng_sta(2),ng_end(2)
            do ix=ng_sta(1),ng_end(1)
              rho_s_out(ix,iy,iz,ii,is)=matbox_read(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
      end do
    end if
  end if
end if

if(comm_is_root(nproc_id_global))then
  read(96) esp0(:itotMST0,:num_kpoints_rd),rocc0(:itotMST0,:num_kpoints_rd)
  if(itotMST0>=itotMST)then
    if(ilsda == 0)then
      is_sta=1
      is_end=1
      pstart(1)=1
      pend(1)=itotMST
    else if(ilsda == 1)then
      is_sta=1
      is_end=2
      pstart(1)=1
      pend(1)=MST(1)
      pstart(2)=MST(1)+1
      pend(2)=itotMST
    end if
    do ik=1,num_kpoints_rd
    do is=is_sta,is_end
      do iob=pstart(is),pend(is)
        call conv_p(iob,p0)
        esp(iob,ik)=esp0(p0,ik)
        rocc(iob,ik)=rocc0(p0,ik)
      end do
    end do
    end do
  else
    if(ilsda == 0)then
      is_sta=1
      is_end=1
      pstart(1)=1
      pend(1)=itotMST0
    else if(ilsda == 1)then
      is_sta=1
      is_end=2
      pstart(1)=1
      pend(1)=MST0(1)
      pstart(2)=MST0(1)+1
      pend(2)=itotMST0
    end if
    esp(:,:)=0.d0
    rocc(:,:)=0.d0
    do ik=1,num_kpoints_rd
    do is=is_sta,is_end
      do iob=pstart(is),pend(is)
        call conv_p0(p0,iob)
        esp(p0,ik)=esp0(iob,ik)
        rocc(p0,ik)=rocc0(iob,ik)
      end do
    end do
    end do
  end if
end if

if(IC<=2)then
  call read_copy_pot(Vh,matbox_read,ig_sta,ig_end)
  
  if(ilsda == 0)then
    call read_copy_pot(Vxc,matbox_read,ig_sta,ig_end)
  else if(ilsda == 1)then
    do is=1,2
      if(comm_is_root(nproc_id_global))then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      call comm_bcast(matbox_read,nproc_group_global)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        Vxc_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end do
  end if
  
  call read_copy_pot(Vpsl,matbox_read,ig_sta,ig_end)
end if

if(comm_is_root(nproc_id_global))then
  if(version_num_box(1)<=31)then
    if(iflag_ps.eq.1)then
      read(96) 
      read(96) Mlps(:MKI),Lref(:MKI)
    end if
  end if

close(96)

end if

call comm_bcast(rocc,nproc_group_global)
call comm_bcast(esp,nproc_group_global)

if(version_num_box(1)<=31)then
  if(iflag_ps.eq.1)then
    call comm_bcast(Mlps,nproc_group_global)
    call comm_bcast(Lref,nproc_group_global)
  end if
end if

if(iSCFRT==2)then
  allocate(Vh_stock1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  allocate(Vh_stock2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
!$OMP parallel do private(iz,iy,ix) 
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    Vh_stock1(ix,iy,iz) = Vh(ix,iy,iz)
    Vh_stock2(ix,iy,iz) = Vh(ix,iy,iz)
  end do
  end do
  end do
end if

call allgatherv_vlocal

deallocate( esp0,rocc0 )

deallocate(matbox,matbox2,matbox3)
deallocate(cmatbox,cmatbox2)

END SUBROUTINE IN_data

