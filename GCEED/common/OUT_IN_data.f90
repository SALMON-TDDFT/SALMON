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

SUBROUTINE OUT_data
use scf_data
use new_world_sub
use read_pslfile_sub
use allocate_psl_sub
use allocate_mat_sub
implicit none
integer :: is,iob,jj
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

if(myrank.eq.0)then

  open(97,file=file_OUT,form='unformatted')
  
!version number
  version_num(1)=39
  version_num(2)=1
  write(97) version_num(1),version_num(2)
  write(97) imesh_oddeven
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
    write(97) Jxyz(1:3,1:maxMps,1:MI),Mps(1:MI)
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

if(myrank==0)then
  if(iflag_ps.eq.1)then
    write(97) uV(:maxMps,:Mlmps,:MI),uVu(:Mlmps,:MI)
    write(97) Mlps(:MKI),Lref(:MKI)
  end if
end if

allocate(matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(matbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(cmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))

if(OC<=2)then
  if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc)then
    file_OUT_data_ini = file_OUT_ini
  else
    if(myrank<num_datafiles_OUT)then
      myrank_datafiles=myrank

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

      write(fileNumber_data, '(i8)') myrank_datafiles
      file_OUT_data = trim(file_OUT)//"."//adjustl(fileNumber_data)
      open(87,file=file_OUT_data,form='unformatted')
    end if
  end if
end if

if(OC<=2)then

  do iob=1,itotMST
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,icorr_p)

    matbox_l=0.d0
    if(icorr_p==1)then
      matbox_l(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))   &
        = psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),iob_myob,1)
    end if

    call MPI_Allreduce(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3), &
&               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


    if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc)then
      if(myrank.eq.0)then
        write(97) ((( matbox_l2(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
      end if
    else
      if(myrank<num_datafiles_OUT)then
        write(87) ((( matbox_l2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                          iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                          iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
      end if
    end if
  end do

else if(OC==3)then
  do iob=1,iobnum
    write(87,rec=iob) ((( psi(ix,iy,iz,iob,1),ix=mg_sta(1),mg_end(1)),   &
                                              iy=mg_sta(2),mg_end(2)),   &
                                              iz=mg_sta(3),mg_end(3))
  end do
end if

if(OC<=2)then
  if(num_datafiles_OUT==1.or.num_datafiles_OUT>nproc)then
    if(myrank==0.and.OC==2) close(67)
  else
    if(myrank<num_datafiles_OUT)then
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

call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

if(myrank.eq.0)then
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

  call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

  if(myrank.eq.0)then
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

  call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)
  if(myrank.eq.0)then
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

    call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

    if(myrank==0)then
      write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
    end if

    do ii=1,num_rho_stock+1
      matbox2=0.d0
      matbox2(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3))   &
        = rho_s_in(ng_sta(1):ng_end(1),   &
                ng_sta(2):ng_end(2),   &
                ng_sta(3):ng_end(3),is,ii)

      call MPI_Allreduce(matbox2,matbox, &
&               lg_num(1)*lg_num(2)*lg_num(3), &
&               MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

      if(myrank==0)then
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
                ng_sta(3):ng_end(3),is,ii)

      call MPI_Allreduce(matbox2,matbox, &
&               lg_num(1)*lg_num(2)*lg_num(3), &
&               MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

      if(myrank==0)then
        write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
      end if
    end do
    
  end do

end if

if(myrank==0)then
  write(97) esp(:itotMST,1),rocc(:itotMST,1)
end if

matbox2=0.d0
matbox2(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))   &
   = Vh(ng_sta(1):ng_end(1),   &
        ng_sta(2):ng_end(2),   &
        ng_sta(3):ng_end(3))

call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

if(myrank.eq.0)then
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

  call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

  if(myrank.eq.0)then
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

    call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

    if(myrank.eq.0)then
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

  call MPI_Allreduce(matbox2,matbox, &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)

if(myrank==0)then
  write(97) ((( matbox(ix,iy,iz),ix=lg_sta(1),lg_end(1)),iy=lg_sta(2),lg_end(2)),iz=lg_sta(3),lg_end(3))
end if

if(myrank==0)then
  close(97)
end if

deallocate(matbox,matbox2)
deallocate(cmatbox,cmatbox2)

END SUBROUTINE OUT_data

!=======================================================================

SUBROUTINE IN_data
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
integer :: NI0,Ndv0,Nps0,Nd0
integer :: ii,is,iob,jj,ibox,j1,j2,j3
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

if(myrank.eq.0)then

   write(*,*) file_IN
   open(96,file=file_IN,form='unformatted')

   read(96) version_num_box(1),version_num_box(2)
   if((version_num_box(1)==17.and.version_num_box(2)>=13).or.version_num_box(1)>=18) then
     read(96) imesh_oddeven
   else
     imesh_oddeven=2
   end if
end if
call MPI_Bcast(version_num_box,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(imesh_oddeven,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(myrank==0) then
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

call MPI_Bcast(ilsda,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_ps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(myrank.eq.0)then
  read(96) iend_Mx_ori(:3)
  read(96) lg_end(:3)
  if(ilsda == 0) then
    read(96) MST0(1)
    read(96) ifMST(1)
  else if(ilsda == 1)then
    read(96) (MST0(is),is=1,2)
    read(96) (ifMST(is),is=1,2)
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

call MPI_Bcast(iend_Mx_ori,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lg_end,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(ilsda == 0) then
   call MPI_Bcast(MST0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ifMST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
else if(ilsda == 1)then
   call MPI_Bcast(MST0,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ifMST,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end if

if(iSCFRT==2) then
  if(ilsda == 0) then
    MST(1)=ifMST(1)
  else if(ilsda == 1) then
    MST(1:2)=ifMST(1:2)
  end if
end if

if(imesh_oddeven==1)then
  ista_Mx_ori(:)=-iend_Mx_ori(:)
  lg_sta(:)=-lg_end(:)
else if(imesh_oddeven==2)then
  ista_Mx_ori(:)=-iend_Mx_ori(:)+1
  lg_sta(:)=-lg_end(:)+1
end if
inum_Mx_ori(:)=iend_Mx_ori(:)-ista_Mx_ori(:)+1

lg_num(:)=lg_end(:)-lg_sta(:)+1

call MPI_Bcast(Hgs,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Hvol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(rLsize,3*ntmg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Miter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call set_gridcoo

allocate(ista_Mxin(3,0:nproc-1),iend_Mxin(3,0:nproc-1))
allocate(inum_Mxin(3,0:nproc-1))

call setmg(mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_end,lg_num,nproc,myrank,nproc_Mxin,nproc_ob,isequential)

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

if(iflag_ps.eq.1)then
  call MPI_Bcast(MI_read,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(MKI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(maxMps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Mlmps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  MI=MI_read
end if



if(iflag_ps.eq.1)then
   if(myrank.eq.0)then
     if(version_num_box(1)<=31)then
       read(96) 
       read(96) 
     else
       read(96) 
     end if
   end if

  if(iSCFRT==2) then
    allocate( Kion(MI),Rion(3,MI) )
  end if
  if(iSCFRT==2) allocate( AtomName(MI), iAtomicNumber(MI) )
  if(myrank.eq.0)then
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
  
  call MPI_Bcast(Kion(1),MI_read,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Rion(1,1),MI_read*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iZatom,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(pseudo_file,256*MKI,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(AtomName(1),8*MI_read,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iAtomicNumber(1),MI_read,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end if

if(ilsda==1)then
  nproc_ob_spin(1)=(nproc_ob+1)/2
  nproc_ob_spin(2)=nproc_ob/2
end if

if(iSCFRT==2) call make_new_world

if(ilsda==0)then
  call calc_iobnum(itotMST,nproc_ob,newrank_comm_grid,iobnum,nproc_ob,iparaway_ob)
else if(ilsda==1)then
  if(nproc_ob==1)then
    iobnum=itotMST
  else
    if(newrank_comm_spin<nproc_ob_spin(1))then
      call calc_iobnum(MST(1),nproc_ob_spin(1),newrank_comm_grid,iobnum,nproc_ob_spin(1),iparaway_ob)
    else
      call calc_iobnum(MST(2),nproc_ob_spin(2),newrank_comm_grid,iobnum,nproc_ob_spin(2),iparaway_ob)
    end if
  end if
end if

if(iSCFRT==2) call allocate_mat

allocate( matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
allocate( cmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
allocate( matbox3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )

if(iSCFRT==1)then
  if(iobnum.ge.1)then
    allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                  mg_sta(3):mg_end(3), &
&                 1:iobnum,1) )
  end if
  if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
    allocate( psi_mesh(ng_sta(1):ng_end(1),  &
                     ng_sta(2):ng_end(2),   &
                     ng_sta(3):ng_end(3), &
                     1:itotMST,1) )
  end if
else if(iSCFRT==2)then
  if(iobnum.ge.1)then
    allocate( zpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum,1) )
    allocate( zpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum,1) )
    zpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                    1:iobnum,1) = 0.d0
    zpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                    1:iobnum,1) = 0.d0
  end if
  if(iwrite_projection==1)then
    if(ilsda==0)then
      call calc_iobnum(itotMST0,nproc_ob,newrank_comm_grid,iobnum0,nproc_ob,iparaway_ob)
    else if(ilsda==1)then
      if(nproc_ob==1)then
        iobnum0=itotMST0
      else
        if(newrank_comm_spin<nproc_ob_spin(1))then
          call calc_iobnum(MST0(1),nproc_ob_spin(1),newrank_comm_grid,iobnum0,nproc_ob_spin(1),iparaway_ob)
        else
          call calc_iobnum(MST0(2),nproc_ob_spin(2),newrank_comm_grid,iobnum0,nproc_ob_spin(2),iparaway_ob)
        end if
      end if
    end if
    if(iobnum0.ge.1)then
      allocate( zpsi_t0(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum0,1) )
      zpsi_t0(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
&                  1:iobnum0,1) = 0.d0
    end if
  end if
end if

allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate( rho0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate( rho_diff(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
if(iSCFRT==1)then
  allocate( rho_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:num_rho_stock+1))
  allocate( rho_out(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:num_rho_stock))
  rho_in=0.d0
  rho_out=0.d0
end if

if(ilsda == 0) then
  continue
else if(ilsda == 1) then
  allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2))
  allocate( rho_s_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2,1:num_rho_stock+1))
  allocate( rho_s_out(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2,1:num_rho_stock))
  rho_s_in=0.d0
  rho_s_out=0.d0
end if

if(iSCFRT==1)then
  allocate( esp0(itotMST0,1))
  allocate( esp(itotMST,1))
  allocate( rocc0(itotMST0,1))
else if(iSCFRT==2)then
  allocate( esp(itotMST,1),rocc(itotMST,1))
  allocate( esp0(itotMST0,1),rocc0(itotMST0,1))
  allocate( esp2(itotMST,1))
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

if(myrank==0)then
  if(version_num_box(1)>=32)then
    if(iflag_ps.eq.1)then
      read(96) 
      read(96) Mlps(:MKI),Lref(:MKI)
    end if
  end if
end if
if(version_num_box(1)>=32)then
  if(iflag_ps.eq.1)then
    call MPI_Bcast(Mlps,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(Lref,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
  if(num_datafiles_IN<=nproc)then
    if(myrank<num_datafiles_IN)then
      myrank_datafiles=myrank

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

      if(num_datafiles_IN>=2.and.myrank<num_datafiles_IN)then
        write(fileNumber_data, '(i8)') myrank_datafiles
        file_IN_data = trim(file_IN)//"."//adjustl(fileNumber_data)
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

!$OMP parallel do
do iz=ig_sta(3),ig_end(3)
do iy=ig_sta(2),ig_end(2)
do ix=ig_sta(1),ig_end(1)
  matbox_read2(ix,iy,iz)=0.d0
  cmatbox_read2(ix,iy,iz)=0.d0
end do
end do
end do

icount=0

do is=is_sta,is_end
do p0=pstart(is),pend(is)

! read file
  call conv_p0(p0,iob)
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,icheck_corrkob)

  if(IC<=2)then
    if(myrank<num_datafiles_IN)then
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
      read(ifilenum_data) ((( matbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                     iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                     iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
    else if(IC==3.or.IC==4)then
      read(ifilenum_data,rec=icount) ((( matbox_read2(ix,iy,iz),ix=ista_Mxin_datafile(1),iend_Mxin_datafile(1)),   &
                                                     iy=ista_Mxin_datafile(2),iend_Mxin_datafile(2)),   &
                                                     iz=ista_Mxin_datafile(3),iend_Mxin_datafile(3))
    end if
  end if
  
  icomm=MPI_COMM_WORLD

  call MPI_Allreduce(matbox_read2,matbox_read,  &
&           ig_num(1)*ig_num(2)*ig_num(3), &
&           MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)

  if(icheck_corrkob==1)then
    if(iSCFRT==1)then
      psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
          mg_sta(3):mg_end(3),iob_myob,1)=  &
      matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
             mg_sta(3):mg_end(3))
    else if(iSCFRT==2)then
      if(iwrite_projection==1)then
        zpsi_t0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3),iob_myob,1)=  &
        matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                mg_sta(3):mg_end(3))
      else
        if((ilsda==0.and.p0<=MST(1)).or.  &
           (ilsda==1.and.(p0<=MST0(1).and.p0<=MST(1)).or.(p0>MST0(1).and.p0<=MST0(1)+MST(2))))then
          zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,1)=  &
          matbox_read(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),   &
                      mg_sta(3):mg_end(3))
        end if
      end if
      if(iwrite_projection==1)then
        if((ilsda==0.and.p0<=MST(1)).or.  &
           (ilsda==1.and.(p0<=MST0(1).and.p0<=MST(1)).or.(p0>MST0(1).and.p0<=MST0(1)+MST(2))))then
          zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,1)=  &
          zpsi_t0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                  mg_sta(3):mg_end(3),iob_myob,1)
        end if
      end if
    end if
  end if
  
end do
end do

if(iSCFRT==1.and.itotMST>itotMST0) call init_wf_ns(2)

if(IC<=2)then
  call read_copy_pot(rho,matbox_read,ig_sta,ig_end,ig_num)
 
  if(version_num_box(1)<=29.or.(version_num_box(1)==30.and.version_num_box(2)<=18))then
    if(myrank.eq.0)then
      read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
    end if
    if(iSCFRT==1)then
      call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rho_in(ix,iy,iz,num_rho_stock+1)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end if
  
    if(myrank.eq.0)then
      read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
    end if
    if(iSCFRT==1)then
      call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rho_out(ix,iy,iz,num_rho_stock)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end if
  
    if(ilsda == 1)then
      do is=1,2
        if(myrank.eq.0)then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        if(myrank.eq.0)then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s_in(ix,iy,iz,is,num_rho_stock)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        if(myrank.eq.0)then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s_out(ix,iy,iz,is,num_rho_stock)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
      end do
    end if
  else
    do ii=1,num_rho_stock+1
      if(myrank.eq.0)then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      if(iSCFRT==1)then
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_in(ix,iy,iz,ii)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  
    do ii=1,num_rho_stock
      if(myrank.eq.0)then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      if(iSCFRT==1)then
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_out(ix,iy,iz,ii)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  
    if(ilsda == 1)then
      do is=1,2
        if(myrank.eq.0)then
          read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
        end if
        call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rho_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
        end do
        end do
        end do
  
        do ii=1,num_rho_stock+1
          if(myrank.eq.0)then
            read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
          end if
          if(iSCFRT==1)then
            call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            do iz=mg_sta(3),mg_end(3)
            do iy=mg_sta(2),mg_end(2)
            do ix=mg_sta(1),mg_end(1)
              rho_s_in(ix,iy,iz,is,ii)=matbox_read(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
    
        do ii=1,num_rho_stock
          if(myrank.eq.0)then
            read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
          end if
          if(iSCFRT==1)then
            call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            do iz=mg_sta(3),mg_end(3)
            do iy=mg_sta(2),mg_end(2)
            do ix=mg_sta(1),mg_end(1)
              rho_s_out(ix,iy,iz,is,ii)=matbox_read(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
      end do
    end if
  end if
end if

if(myrank==0)then
  read(96) esp0(:itotMST0,1),rocc0(:itotMST0,1)
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
    do is=is_sta,is_end
      do iob=pstart(is),pend(is)
        call conv_p(iob,p0)
        esp(iob,1)=esp0(p0,1)
        rocc(iob,1)=rocc0(p0,1)
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
    do is=is_sta,is_end
      do iob=pstart(is),pend(is)
        call conv_p0(p0,iob)
        esp(p0,1)=esp0(iob,1)
        rocc(p0,1)=rocc0(iob,1)
      end do
    end do
  end if
end if

if(IC<=2)then
  call read_copy_pot(Vh,matbox_read,ig_sta,ig_end,ig_num)
  
  if(ilsda == 0)then
    call read_copy_pot(Vxc,matbox_read,ig_sta,ig_end,ig_num)
  else if(ilsda == 1)then
    do is=1,2
      if(myrank==0)then
        read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1)),iy=ig_sta(2),ig_end(2)),iz=ig_sta(3),ig_end(3))
      end if
      call MPI_Bcast(matbox_read,lg_num(1)*lg_num(2)*lg_num(3),&
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        Vxc_s(ix,iy,iz,is)=matbox_read(ix,iy,iz)
      end do
      end do
      end do
    end do
  end if
  
  call read_copy_pot(Vpsl,matbox_read,ig_sta,ig_end,ig_num)
end if

if(myrank==0)then
  if(version_num_box(1)<=31)then
    if(iflag_ps.eq.1)then
      read(96) 
      read(96) Mlps(:MKI),Lref(:MKI)
    end if
  end if

close(96)

end if

call MPI_Bcast(rocc,itotMST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(esp,itotMST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

if(version_num_box(1)<=31)then
  if(iflag_ps.eq.1)then
    call MPI_Bcast(Mlps,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(Lref,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  end if
end if

if(iSCFRT==2)then
  allocate(Vh_stock1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  allocate(Vh_stock2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    Vh_stock1(ix,iy,iz) = Vh(ix,iy,iz)
    Vh_stock2(ix,iy,iz) = Vh(ix,iy,iz)
  end do
  end do
  end do
end if

allocate(icoo1d(3,lg_num(1)*lg_num(2)*lg_num(3)))
icount=0
do iz=lg_sta(3),lg_end(3),1
do iy=lg_sta(2),lg_end(2),1
do ix=lg_sta(1),lg_end(1),1
  icount=icount+1
  icoo1d(1,icount)=ix
  icoo1d(2,icount)=iy
  icoo1d(3,icount)=iz
end do
end do
end do

call mpi_allgatherv_vlocal

deallocate( esp0,rocc0 )

deallocate(matbox,matbox2,matbox3)
deallocate(cmatbox,cmatbox2)

END SUBROUTINE IN_data

