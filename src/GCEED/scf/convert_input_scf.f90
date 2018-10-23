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
subroutine convert_input_scf(file_atoms_coo)
use salmon_global
use salmon_parallel, only: nproc_group_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_bcast
use inputoutput
use scf_data
use new_world_sub
implicit none
integer :: ii,iatom
integer :: ibox2
integer :: icheck1,icheck2
character(100) :: file_atoms_coo
real(8) :: dip_spacing

ik_oddeven=2
iterVh = 0         ! Iteration counter
ilsda = ispin

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if
!===== namelist for group_fundamental =====

select case(amin_routine)
  case('cg','diis','cg-diis')
    continue
  case default
    stop 'Specify either "cg", "diis", or "cg-diis" for amin_routine.'
end select

select case(amixing)
  case ('simple','broyden')
    continue
  case default
    stop 'Specify either "simple" or "broyden" for amixing.'
end select
    
Harray(1:3,1:maxntmg)=0.d0
rLsize(1:3,1:maxntmg)=0.d0
iDiter(1:maxntmg)=1000

Harray(1:3,1)=dl(1:3)
rLsize(1:3,1)=al(1:3)
iDiter(1) = nscf

if(ispin == 0)then
  MST(1)=nstate
  if(temperature_k>=0.d0)then
    ifMST(1)=nstate
    rNetot=dble(nelec)
  else
    ifMST(1)=nelec/2
  end if
  MST(2)=0
  ifMST(2)=0
else if(ispin == 1)then
  if(nstate /= 0 .and. sum(nstate_spin) ==0)then
    MST(1:2)=nstate
    if(temperature_k>=0.d0)then
      ifMST(1:2)=nstate
      rNetot=dble(nelec)
    else
      ifMST(1)=nelec - nelec/2
      ifMST(2)=nelec/2
    end if
  else if(nstate == 0 .and. sum(nstate_spin) /=0)then
    MST(1:2)=nstate_spin(1:2)
    if(temperature_k>=0.d0)then
      ifMST(1:2)=nstate_spin(1:2)
      rNetot=dble(nelec_spin(1))+dble(nelec_spin(2))
    else
      ifMST(1:2)=nelec_spin(1:2)
    end if
  else
    write(*,*)"'nstate' or 'nstate_spin' should be spacified in input. "
  end if
else 
  write(*,*)"'ispin' should be 0 or 1. "
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

select case(use_force)
case('y')
  icalcforce = 1
case('n')
  icalcforce = 0
end select

num_kpoints_3d(1:3)=num_kgrid(1:3)
num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)
allocate(wtk(num_kpoints_rd))
wtk(:)=1.d0/dble(num_kpoints_rd)

if(comm_is_root(nproc_id_global))then
  if(iflag_stopt==1.and.icalcforce==0)then
    write(*,*) "use_force should be set to 'y' when use_geometry_opt is 'y'"
    stop
  end if
end if

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

allocate( rocc(itotMST,num_kpoints_rd) )

rocc(:,:)=0.d0             
if(ilsda == 0) then
  rocc(:ifMST(1),:num_kpoints_rd) = 2.d0   ! Occupation number
else if(ilsda == 1) then
  rocc(:ifMST(1),:num_kpoints_rd) = 1.d0   ! Occupation number
  rocc(MST(1)+1:MST(1)+ifMST(2),:num_kpoints_rd) = 1.d0   ! Occupation number
end if

!===== namelist for group_parallel =====
isequential=2
imesh_s_all=1
if(comm_is_root(nproc_id_global))then
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

if(nproc_ob==0.and.nproc_mxin(1)==0.and.nproc_mxin(2)==0.and.nproc_mxin(3)==0.and.  &
                   nproc_mxin_s(1)==0.and.nproc_mxin_s(2)==0.and.nproc_mxin_s(3)==0) then
  call set_numcpu_scf
else
  call check_numcpu
end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

!===== namelist for group_hartree =====
if(meo<=0.or.meo>=4)then
  stop "meo must be equal to 1 or 2 or 3."
else if(meo==3)then
  if(num_pole_xyz(1)==0.and.num_pole_xyz(2)==0.and.num_pole_xyz(3)==0)then
    continue
  else if(num_pole_xyz(1)<=0.or.num_pole_xyz(2)<=0.or.num_pole_xyz(3)<=0)then
    stop "num_pole_xyz must be largar than 0 when they are not default values."
  end if
  if(num_pole_xyz(1)==0.and.num_pole_xyz(2)==0.and.num_pole_xyz(3)==0)then
    dip_spacing = 8.d0/au_length_aa  ! approximate spacing of multipoles 
    num_pole_xyz(:)=int((al(:)+dip_spacing)/dip_spacing-1.d-8)
  end if
end if

num_pole=num_pole_xyz(1)*num_pole_xyz(2)*num_pole_xyz(3)

!===== namelist for group_file =====
if(comm_is_root(nproc_id_global))then
  if(IC<0.or.IC>=2)then
    write(*,*) "IC must be equal to 0 or 1."
    stop
  end if
end if

!===== namelist for group_atom =====
iflag_ps=1
MI=natom
MKI=nelem
!iZatom(:)=0
!ipsfileform(:)=1
!ps_format = 'default'
if(file_atom_coor/='none')then
   file_atoms_coo=trim(file_atom_coor)
else
   file_atoms_coo='.atomic_coor.tmp'
end if

!Lmax_ps(:)=-1
!Lloc_ps(:)=-1
if(comm_is_root(nproc_id_global))then
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

call comm_bcast(Mlps,nproc_group_global)
call comm_bcast(Lref,nproc_group_global)

if(comm_is_root(nproc_id_global)) write(*,*) "MI =",MI

allocate(istopt_a(MI) ); istopt_a = 0
allocate( AtomName(MI),iAtomicNumber(MI) )

if(comm_is_root(nproc_id_global)) then
  do iatom=1,MI
    if(flag_geo_opt_atom(iatom) == 'y')istopt_a(iatom)=1
  end do
end if

call comm_bcast(istopt_a,      nproc_group_global)

Atomname(:)='dm' ! dummy
iAtomicNumber(:)=0 ! dummy
!call comm_bcast(AtomName,      nproc_group_global)
!call comm_bcast(iAtomicNumber, nproc_group_global)


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

call comm_bcast(iflag_writepsi, nproc_group_global)
call comm_bcast(iflag_ELF,      nproc_group_global)
call comm_bcast(iflag_dos,      nproc_group_global)
call comm_bcast(iflag_pdos,     nproc_group_global)

!===== namelist for group_others =====
if(comm_is_root(nproc_id_global))then
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

!if(iflag_ps==1)then
!  do ii=1,3
!    do iatom=1,MI
!      Rion(ii,iatom)=Rion(ii,iatom)*ulength_to_au
!    end do
!  end do
!end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

select case(fourier)
case('ft','FT')
  iflag_hartree=2
case('ffte','FFTE')
  iflag_hartree=4
end select

if(temperature>=0.d0)then
  write(*,*) "At the moment, temperature must be given in a variable temperature_k"
  stop 
end if

call make_new_world

if(comm_is_root(nproc_id_global))close(fh_namelist)

return

end subroutine convert_input_scf
