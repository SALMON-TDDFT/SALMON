! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi, Shunsuke A. Sato
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

subroutine read_input_scf(file_IN,file_OUT,LDA_Info,file_ini,iDiterYBCG,file_atoms_coo)
use input
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ii,iatom
integer :: ibox2
integer :: icheck1,icheck2
character(100) :: file_OUT,file_IN,LDA_Info,file_ini
character(100) :: file_atoms_coo
integer :: iDiterYBCG
namelist / group_fundamental / imesh_oddeven,iflag_stopt,iter_stopt,Ncg,iDiterYBCG,   &
                               iflag_subspace_diag,iDiter_nosubspace_diag, &
                               Harray, rLsize, ntmg, iDiter, ilsda,  &
                               MST, ifMST, iflag_convergence, ithresholdVh, threshold_norm_diff_rho, &
                               threshold_square_norm_diff_Vlocal,      &
                               mixrate, Nmemory_MB, icalcforce
namelist / group_parallel/ nproc_ob,nproc_Mxin,nproc_Mxin_s,  &
                           isequential,num_datafiles_IN,num_datafiles_OUT,imesh_s_all
namelist / group_hartree / Hconv, MEO, num_pole_xyz, lmax_MEO
namelist / group_file / IC,OC,file_IN,file_OUT,LDA_Info
namelist / group_atom / MI,MKI,iZatom,ipsfileform,file_atoms_coo, Mlps, Lref
namelist / group_scf_analysis / iflag_writepsi, iflag_ELF, iflag_dos, iflag_pdos
namelist / group_others / iparaway_ob,iscf_order,iswitch_orbital_mesh,iflag_psicube,  &
                          lambda1_diis, lambda2_diis, file_ini

! 'Hartree' parameter

iterVh = 0         ! Iteration counter


if(myrank ==0)then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if
!===== namelist for group_fundamental =====
imesh_oddeven=2
ik_oddeven=2
iflag_stopt=0
iter_stopt=1
minroutine=1
Ncg=4
iDiterYBCG=20
iflag_subspace_diag=0
iDiter_nosubspace_diag=10
ntmg=1
iflag_convergence=2
ithresholdVh(:)=1.d-11
threshold_norm_diff_rho(:)=1.d-11
threshold_square_norm_diff_Vlocal(:)=1.d-11
mixrate=0.1d0
Nmemory_MB=8
icalcforce=0

Harray(1:3,1:maxntmg)=0.d0
rLsize(1:3,1:maxntmg)=0.d0
iDiter(1:maxntmg)=1000

ilsda=0

MST(1:2)=0
ifMST(1:2)=0

if(myrank==0)then
  read(fh_namelist,NML=group_fundamental) 
  rewind(fh_namelist)
  if(iflag_stopt==0) iter_stopt=1 ! overwrite iter_stopt
  if(imesh_oddeven<=0.or.imesh_oddeven>=3)then
    write(*,*) "imesh_oddeven must be equal to 1 or 2."
    stop
  end if
  if(minroutine<=0.or.minroutine==2.or.minroutine==3.or.minroutine>=5)then
    write(*,*) "minroutine must be equal to 1 or 4."
    stop
  end if
  if(iflag_subspace_diag<=-1.or.iflag_subspace_diag>=2)then
    write(*,*) "iflag_subspace_diag must be equal to 0 or 1."
    stop
  end if
end if

call MPI_Bcast(imesh_oddeven,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_stopt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iter_stopt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(minroutine,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Ncg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iDiterYBCG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_subspace_diag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iDiter_nosubspace_diag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

allocate(wtk(1))
wtk(:)=1.d0

call MPI_Bcast(iflag_convergence,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_Bcast(ithresholdVh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(threshold_norm_diff_rho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
threshold_norm_diff_rho=threshold_norm_diff_rho*a_B**3
call MPI_Bcast(threshold_square_norm_diff_Vlocal,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
threshold_square_norm_diff_Vlocal=threshold_square_norm_diff_Vlocal/(2.d0*Ry)**2/a_B**3

call MPI_Bcast(mixrate,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Nmemory_MB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_Bcast(icalcforce,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(myrank==0)then
  if(iflag_stopt==1.and.icalcforce==0)then
    write(*,*) "icalcforce should be set to 1 when iflag_stopt = 1"
    stop
  end if
end if

call MPI_Bcast(ntmg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Harray,3*10,MPI_DOUBLE_PRECISION,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(rLsize,3*10,MPI_DOUBLE_PRECISION,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iDiter,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(ilsda,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_Bcast(MST,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(ifMST,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

Harray=Harray/a_B
rLsize=rLsize/a_B

if(ilsda==1)then
  nproc_ob_spin(1)=(nproc_ob+1)/2
  nproc_ob_spin(2)=nproc_ob/2
end if

if(ilsda == 0) then
  itotMST=MST(1)
  itotfMST=ifMST(1)
else if(ilsda == 1) then
  itotMST=MST(1)+MST(2)
  itotfMST=ifMST(1)+ifMST(2)
end if

allocate( rocc(itotMST,1) )

rocc(:,:)=0.d0             
if(ilsda == 0) then
  rocc(:ifMST(1),1) = 2.d0   ! Occupation number
else if(ilsda == 1) then
  rocc(:ifMST(1),1) = 1.d0   ! Occupation number
  rocc(MST(1)+1:MST(1)+ifMST(2),1) = 1.d0   ! Occupation number
end if

!===== namelist for group_parallel =====
nproc_ob=0
nproc_Mxin(1:3)=0
nproc_Mxin_s(1:3)=0
isequential=2
num_datafiles_IN=1
num_datafiles_OUT=1
imesh_s_all=1
if(myrank==0)then
  read(fh_namelist,NML=group_parallel)
  rewind(fh_namelist)
  ibox2=1
  icheck1=0
  icheck2=0
  do ii=1,19
    if(num_datafiles_IN==ibox2) icheck1=1
    if(num_datafiles_OUT==ibox2) icheck2=1
    ibox2=ibox2*2 
  end do
  if(icheck1/=1.or.icheck2/=1)then
    write(*,*) "num_datafiles_IN and num_datafiles_OUT must be equal to nth power of 2. (n: positive integer)"
    stop
  end if
end if
call MPI_Bcast(nproc_ob,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(nproc_Mxin,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(nproc_Mxin_s,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(isequential,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_datafiles_IN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_datafiles_OUT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(imesh_s_all,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myrank==0.and.nproc_ob==0)then
  write(*,*) "set nproc_ob."
  stop
else if(myrank==0.and.nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)==0)then
  write(*,*) "set nproc_Mxin."
  stop
else if(myrank==0.and.nproc_Mxin_s(1)*nproc_Mxin_s(2)*nproc_Mxin_s(3)==0)then
  write(*,*) "set nproc_Mxin_s."
  stop
end if

call check_numcpu

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

!===== namelist for group_hartree =====
Hconv=1.d-8
num_pole_xyz(1:3)=-1
MEO=3
lmax_MEO=4
if(myrank==0)then
  read(fh_namelist,NML=group_hartree) 
  rewind(fh_namelist)
  if(MEO<=0.or.MEO>=4)then
    write(*,*) "MEO must be equal to 1 or 2 or 3."
    stop
  end if
  if(MEO==3.and.(num_pole_xyz(1)==-1.or.num_pole_xyz(2)==-1.or.num_pole_xyz(3)==-1))then
    write(*,*) "num_pole_xyz must be set when MEO=3."
    stop
  end if
end if
call MPI_Bcast(Hconv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(MEO,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_pole_xyz,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lmax_MEO,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

num_pole=num_pole_xyz(1)*num_pole_xyz(2)*num_pole_xyz(3)
Hconv  = Hconv/(2d0*Ry)**2d0/a_B**3     ! Convergence criterion

!===== namelist for group_file =====
IC=0
OC=1
file_IN='file_IN'
file_OUT='file_OUT'
LDA_Info='LDA_Info'
if(myrank==0)then
  read(fh_namelist,NML=group_file) 
  rewind(fh_namelist)
  if(IC<0.or.IC>=2)then
    write(*,*) "IC must be equal to 0 or 1."
    stop
  end if
end if
call MPI_Bcast(IC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(OC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_IN,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_OUT,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!===== namelist for group_atom =====
iflag_ps=1
MI=0
MKI=0 
iZatom(:)=0
ipsfileform(:)=1
file_atoms_coo='file_atoms_coo'
Mlps(:)=-1
Lref(:)=-1
if(myrank==0)then
  read(fh_namelist,NML=group_atom) 
  rewind(fh_namelist)
end if

call MPI_Bcast(MI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(MKI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myrank==0) write(*,*) "MI =",MI

call MPI_Bcast(iZatom,MKI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(ipsfileform,MKI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Mlps,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Lref,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

allocate( Rion(3,MI),Kion(MI),istopt_a(MI) )
allocate( AtomName(MI),iAtomicNumber(MI) )

if(myrank.eq.0) then
  open(91,file=file_atoms_coo)
  do iatom=1,MI
    if(iflag_stopt==1)then
      read(91,*) AtomName(iatom), Rion(1,iatom),Rion(2,iatom),Rion(3,iatom),      &
           Kion(iatom),istopt_a(iatom)
    else
      read(91,*) AtomName(iatom), Rion(1,iatom),Rion(2,iatom),Rion(3,iatom),      &
           Kion(iatom)
      istopt_a(iatom)=0
    end if
    if(AtomName(iatom)=="Na".or. AtomName(iatom)=="na")then
      iAtomicNumber(iatom) = 11
    elseif(AtomName(iatom)=="H" .or. AtomName(iatom)=="h") then
      iAtomicNumber(iatom) = 1
    elseif(AtomName(iatom)=="C" .or. AtomName(iatom)=="c") then
      iAtomicNumber(iatom) = 6
    elseif(AtomName(iatom)=="Ag".or. AtomName(iatom)=="ag")then
      iAtomicNumber(iatom) = 47
    elseif(AtomName(iatom)=="Au".or. AtomName(iatom)=="au")then
      iAtomicNumber(iatom) = 79
    elseif(AtomName(iatom)=="O" .or. AtomName(iatom)=="o") then
      iAtomicNumber(iatom) = 8
    elseif(AtomName(iatom)=="N" .or. AtomName(iatom)=="n") then
      iAtomicNumber(iatom) = 7
    endif
  end do
  close(91)
end if

call MPI_Bcast(Rion,3*MI,MPI_DOUBLE_PRECISION,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Kion,MI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(istopt_a,MI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(AtomName,8*MI,MPI_CHARACTER,      &
          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iAtomicNumber,MI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)


iflag_writepsi=0
iflag_ELF=0
iflag_dos=0
iflag_pdos=0

if(myrank==0)then
  read(fh_namelist,NML=group_scf_analysis)
  rewind(fh_namelist)
  if(iflag_writepsi<0.or.iflag_writepsi>1)then
    write(*,*) "iflag_writepsi must be equal to 0 or 1."
    stop
  end if
  if(iflag_ELF<0.or.iflag_ELF>1)then
    write(*,*) "iflag_ELF must be equal to 0 or 1."
    stop
  end if
  if(iflag_dos<0.or.iflag_dos>1)then
    write(*,*) "iflag_dos must be equal to 0 or 1."
    stop
  end if
  if(iflag_pdos<0.or.iflag_pdos>1)then
    write(*,*) "iflag_pdos must be equal to 0 or 1."
    stop
  end if
end if 
call MPI_Bcast(iflag_writepsi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_ELF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_dos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_pdos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!===== namelist for group_others =====
iparaway_ob=2
iscf_order=1
iswitch_orbital_mesh=0
iflag_psicube=0
lambda1_diis=0.5d0
lambda2_diis=0.3d0
if(myrank==0)then
  read(fh_namelist,NML=group_others) 
  rewind(fh_namelist)
  if(iparaway_ob<=0.or.iparaway_ob>=3)then
    write(*,*) "iparaway_ob must be equal to 1 or 2."
    stop
  end if
  if(iflag_dos<=-1.or.iflag_dos>=2)then
    write(*,*) "iflag_dos must be equal to 0 or 1."
    stop
  end if
  if(iflag_pdos<=-1.or.iflag_pdos>=2)then
    write(*,*) "iflag_pdos must be equal to 0 or 1."
    stop
  end if
end if

call MPI_Bcast(iparaway_ob,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iscf_order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iswitch_orbital_mesh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_psicube,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lambda1_diis,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lambda2_diis,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_ini,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_dos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_pdos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(iflag_ps==1)then
  do ii=1,3
    do iatom=1,MI
      Rion(ii,iatom)=Rion(ii,iatom)/a_B
    end do
  end do
end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
call make_new_world

if(myrank ==0)close(fh_namelist)

return

end subroutine read_input_scf
