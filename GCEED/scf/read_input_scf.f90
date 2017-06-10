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
subroutine read_input_scf(iDiterYBCG,file_atoms_coo)
use salmon_global
use inputoutput
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ii,iatom
integer :: ibox2
integer :: icheck1,icheck2
character(100) :: file_atoms_coo
integer :: iDiterYBCG
integer :: inml_group_fundamental, &
         & inml_group_parallel, &
         & inml_group_hartree, &
         & inml_group_file, &
         & inml_group_others
namelist / group_fundamental / imesh_oddeven,iDiterYBCG,   &
                               iDiter_nosubspace_diag, &
                               ntmg, MST, &
                               ifMST, ithresholdVh, threshold_norm_diff_rho, &
                               threshold_square_norm_diff_Vlocal

namelist / group_parallel/ isequential,imesh_s_all

namelist / group_hartree / Hconv, lmax_MEO
namelist / group_file / IC,OC
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
!Ncg=4
iDiterYBCG=20
iflag_subspace_diag=0
iDiter_nosubspace_diag=10
ntmg=1
iflag_convergence=2
ithresholdVh(:)=1

! Convergence criterion, ||rho(i)-rho(i-1)||**2/(# of grids), 1.d-17 a.u. = 6.75d-17 AA**(-3)
threshold_norm_diff_rho(:)=1.d-17/ulength_from_au**3

! Convergence criterion, ||Vlocal(i)-Vlocal(i-1)||**2/(# of grids), 1.d-17 a.u. = 1.10d-15 eV**2*AA**3
threshold_square_norm_diff_Vlocal(:)=1.d-17*uenergy_from_au**2*ulength_from_au**3
mixrate=0.1d0
icalcforce=0

Harray(1:3,1:maxntmg)=0.d0
rLsize(1:3,1:maxntmg)=0.d0
iDiter(1:maxntmg)=1000

ilsda=0

MST(1:2)=0
ifMST(1:2)=0

if(myrank==0)then
  read(fh_namelist,NML=group_fundamental, iostat=inml_group_fundamental) 
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
end if


select case(use_geometry_opt)
case('y')
  iflag_stopt = 1
case('n')
  iflag_stopt = 0
case default
  stop 'invalid iflag_stopt'
end select


iter_stopt = ngeometry_opt

select case(subspace_diagonalization)
case('y')
  iflag_subspace_diag = 1
case('n')
  iflag_subspace_diag = 0
end select

Harray(1:3,1)=dl(1:3)
rLsize(1:3,1)=al(1:3)

iDiter = nscf
ilsda = ispin

select case(convergence)
case('rho')
  iflag_convergence = 1
case('vh')
  iflag_convergence = 2
case default
  stop 'invalid   iflag_convergence'
end select

mixrate = rmixrate

!Nmemory_MB = Nmemory_MB
select case(use_force)
case('y')
  icalcforce = 1
case('n')
  icalcforce = 0
end select

call MPI_Bcast(imesh_oddeven,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_stopt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iter_stopt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(minroutine,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!call MPI_Bcast(Ncg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
!nproc_ob=0
nproc_Mxin(1:3)=0
nproc_Mxin_s(1:3)=0
isequential=2
imesh_s_all=1
if(myrank==0)then
  read(fh_namelist,NML=group_parallel, iostat=inml_group_parallel)
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

nproc_Mxin = nproc_domain
nproc_Mxin_s = nproc_domain_s

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
! Convergence criterion, ||Vh(i)-Vh(i-1)||**2/(# of grids), 1.d-15 a.u. = 1.10d-13 eV**2*AA**3
Hconv=1.d-15*uenergy_from_au**2*ulength_from_au**3

!num_pole_xyz(1:3)=-1
!MEO=3
lmax_MEO=4
if(myrank==0)then
  read(fh_namelist,NML=group_hartree, iostat=inml_group_hartree) 
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
Hconv  = Hconv*uenergy_to_au**2*ulength_to_au**3

!===== namelist for group_file =====
IC=0
OC=1
if(myrank==0)then
  read(fh_namelist,NML=group_file, iostat=inml_group_file ) 
  rewind(fh_namelist)
  if(IC<0.or.IC>=2)then
    write(*,*) "IC must be equal to 0 or 1."
    stop
  end if
end if
call MPI_Bcast(IC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(OC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!===== namelist for group_atom =====
iflag_ps=1
MI=natom
MKI=nelem
!iZatom(:)=0
!ipsfileform(:)=1
!ps_format = 'default'
file_atoms_coo=trim(file_atom_coor)
!Lmax_ps(:)=-1
!Lloc_ps(:)=-1
if(myrank==0)then
!  read(fh_namelist,NML=group_atom) 
!  rewind(fh_namelist)

  Mlps(:) = Lmax_ps(:)
  Lref(:) = Lloc_ps(:)

!ps format conversion
!  do iatom = 1,MKI
!    select case(ps_format(iatom))
!    case('default')
!    case('KY')        ; ipsfileform(iatom)=n_Yabana_Bertsch_psformat
!    case('ABINIT')    ; ipsfileform(iatom)=n_ABINIT_psformat
!    case('FHI')       ; ipsfileform(iatom)=n_FHI_psformat
!    case('ABINITFHI') !; ipsfileform(iatom)=n_ABINITFHI_psformat
!      write(*,"(A)") "Invalid ps_format. ABINITFHI format is not supported for isolated systems."
!    stop
!    case default
!      write(*,"(A)") "Invalid ps_format."
!    stop
!    end select
!  end do
end if


call MPI_Bcast(MI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(MKI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myrank==0) write(*,*) "MI =",MI

call MPI_Bcast(iZatom,MKI,MPI_INTEGER,      &
          0,MPI_COMM_WORLD,ierr)
!call MPI_Bcast(ipsfileform,MKI,MPI_INTEGER,      &
!          0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Mlps,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Lref,MKI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

allocate(istopt_a(MI) ); istopt_a = 0
allocate( AtomName(MI),iAtomicNumber(MI) )

if(myrank.eq.0) then
  do iatom=1,MI
    if(flag_geo_opt_atom(iatom) == 'y')istopt_a(iatom)=1
  end do
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


select case(out_psi)
case('y')
  iflag_writepsi=1
case('n')
  iflag_writepsi=0
case default
  stop 'invalid iflag_writepsi'
end select

select case(out_elf_rt)
case('y')
iflag_ELF=1
case('n')
iflag_ELF=0
case default
  stop 'invalid iflag_ELF'
end select

select case(out_dos)
case('y')
  iflag_dos=1
case('n')
  iflag_dos=0
case default
  stop 'invalid iflag_dos'
end select

select case(out_pdos)
case('y')
  iflag_pdos=1
case('n')
  iflag_pdos=0
case default
  stop 'invalid iflag_pdos'
end select


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
  read(fh_namelist,NML=group_others, iostat=inml_group_others) 
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

!if(iflag_ps==1)then
!  do ii=1,3
!    do iatom=1,MI
!      Rion(ii,iatom)=Rion(ii,iatom)*ulength_to_au
!    end do
!  end do
!end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
call make_new_world

if(myrank ==0)close(fh_namelist)

return

end subroutine read_input_scf
