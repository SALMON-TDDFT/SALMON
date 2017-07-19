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
subroutine read_input_rt(IC_rt,OC_rt,Ntime)
use salmon_parallel, only: nproc_id_global, nproc_group_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_bcast
use inputoutput
use scf_data
use new_world_sub
implicit none
integer :: ii
integer :: IC_rt,OC_rt
integer :: Ntime
integer :: inml_group_fundamental, &
         & inml_group_parallel, &
         & inml_group_hartree, &
         & inml_group_file, &
         & inml_group_others
real(8) :: dip_spacing


namelist / group_fundamental / idisnum, iwrite_projection, &
                               itwproj, iwrite_projnum, itcalc_ene
namelist / group_parallel /  isequential, imesh_s_all, iflag_comm_rho
namelist / group_hartree / Hconv, lmax_MEO
namelist / group_file / IC,IC_rt,OC_rt
namelist / group_others / iparaway_ob,num_projection,iwrite_projection_ob,iwrite_projection_k,  &
                          filename_pot, &
    & iwrite_external,iflag_dip2,iflag_intelectron,num_dip2, dip2boundary, dip2center,& 
    & iflag_fourier_omega, num_fourier_omega, fourier_omega, itotNtime2, &
    & iwdenoption,iwdenstep, iflag_Estatic

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if

!===== namelist for group_fundamental =====
icalcforce=0
iflag_md=0
idisnum(1)=1
idisnum(2)=2
iwrite_projection=0
iwrite_projnum=0
itwproj=-1
itcalc_ene=1
if(comm_is_root(nproc_id_global))then
  read(fh_namelist,NML=group_fundamental, iostat=inml_group_fundamental)
  rewind(fh_namelist)
end if

select case(use_force)
case('y')
  icalcforce = 1
case('n')
  icalcforce = 0
case default
  stop 'invald icalcforce'
end select

select case(use_ehrenfest_md)
case('y')
  iflag_md = 1
case('n')
  iflag_md = 0
case default
  stop 'invald iflag_md'
end select

call comm_bcast(Nenergy,           nproc_group_global)
call comm_bcast(N_hamil,           nproc_group_global)
call comm_bcast(icalcforce,        nproc_group_global)
call comm_bcast(iflag_md,          nproc_group_global)
call comm_bcast(idisnum,           nproc_group_global)
call comm_bcast(iwrite_projection, nproc_group_global)
call comm_bcast(iwrite_projnum,    nproc_group_global)
call comm_bcast(itwproj,           nproc_group_global)
call comm_bcast(itcalc_ene,        nproc_group_global)

allocate(wtk(1))
wtk(:)=1.d0

if(iwrite_projection==1.and.itwproj==-1)then
  write(*,*) "Please specify itwproj when iwrite_projection=1."
  stop
end if

!===== namelist for group_parallel =====
isequential=2
imesh_s_all=1
iflag_comm_rho=1
if(comm_is_root(nproc_id_global))then
  read(fh_namelist,NML=group_parallel, iostat=inml_group_parallel)
  rewind(fh_namelist)
  if(isequential<=0.or.isequential>=3)then
    write(*,*) "isequential must be equal to 1 or 2."
    stop
  end if
end if

call comm_bcast(isequential,       nproc_group_global)
call comm_bcast(num_datafiles_IN,  nproc_group_global)
call comm_bcast(num_datafiles_OUT, nproc_group_global)
call comm_bcast(imesh_s_all,       nproc_group_global)
call comm_bcast(iflag_comm_rho,    nproc_group_global)

nproc_Mxin = nproc_domain
nproc_Mxin_s = nproc_domain_s

if(nproc_ob==0.and.nproc_mxin(1)==0.and.nproc_mxin(2)==0.and.nproc_mxin(3)==0.and.  &
                   nproc_mxin_s(1)==0.and.nproc_mxin_s(2)==0.and.nproc_mxin_s(3)==0) then
  call set_numcpu_rt
else
  call check_numcpu
end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

!===== namelist for group_propagation =====
!dt=0.d0
Ntime=nt
!if(comm_is_root(nproc_id_global))then
!  read(fh_namelist,NML=group_propagation)
!  rewind(fh_namelist)
!end if
call comm_bcast(dt,    nproc_group_global)
call comm_bcast(Ntime, nproc_group_global)
if(dt<=1.d-10)then
  write(*,*) "please set dt."
  stop
end if
if(Ntime==0)then
  write(*,*) "please set Ntime."
  stop
end if

!===== namelist for group_hartree =====
! Convergence criterion, ||Vh(i)-Vh(i-1)||**2/(# of grids), 1.d-15 a.u. = 1.10d-13 eV**2*AA**3
Hconv=1.d-15*uenergy_from_au**2*ulength_from_au**3
lmax_MEO=4
if(comm_is_root(nproc_id_global))then
  read(fh_namelist,NML=group_hartree, iostat=inml_group_hartree ) 
  rewind(fh_namelist)
end if
call comm_bcast(Hconv,        nproc_group_global)
Hconv  = Hconv*uenergy_to_au**2*ulength_to_au**3
call comm_bcast(lmax_MEO,     nproc_group_global)

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
IC=1
IC_rt=0
OC_rt=0
if(comm_is_root(nproc_id_global))then
  read(fh_namelist,NML=group_file, iostat=inml_group_file)
  rewind(fh_namelist)
end if
call comm_bcast(IC,    nproc_group_global)
call comm_bcast(IC_rt, nproc_group_global)
call comm_bcast(OC_rt, nproc_group_global)

if(IC==3.and.num_datafiles_IN/=nproc_size_global)then
  if(comm_is_root(nproc_id_global))then
    write(*,*) "num_datafiles_IN is set to nproc."
  end if
  num_datafiles_IN=nproc_size_global
end if

!===== namelist for group_extfield =====
if(ae_shape1 == 'impulse')then
  ikind_eext = 0
else
  ikind_eext = 1
end if

Fst = e_impulse
romega = omega1
pulse_T = pulse_tw1
rlaser_I = rlaser_int1

if(epdir_im1(1)**2+epdir_im1(2)**2+epdir_im1(3)**2+ &
   epdir_im2(1)**2+epdir_im2(2)**2+epdir_im2(3)**2>=1.d-12)then
  circular='y'
else
  circular='n'
end if

!===== namelist for group_others =====

iparaway_ob=2
num_projection=1
do ii=1,200
  iwrite_projection_ob(ii)=ii
  iwrite_projection_k(ii)=1
end do
filename_pot='pot'
iwrite_external=0
iflag_dip2=0
iflag_intelectron=0
num_dip2=1
dip2boundary(:)=0.d0
dip2center(:)=0.d0
iflag_fourier_omega=0
num_fourier_omega=1
fourier_omega(:)=0.d0
itotNtime2=Ntime
iwdenoption=0
iwdenstep=0
iflag_Estatic=0
if(comm_is_root(nproc_id_global))then
  read(fh_namelist,NML=group_others, iostat=inml_group_others)
  rewind(fh_namelist)
end if
call comm_bcast(iparaway_ob,          nproc_group_global)
call comm_bcast(num_projection,       nproc_group_global)
call comm_bcast(iwrite_projection_ob, nproc_group_global)
call comm_bcast(iwrite_projection_k,  nproc_group_global)
call comm_bcast(filename_pot,         nproc_group_global)
call comm_bcast(iwrite_external,      nproc_group_global)
call comm_bcast(iflag_dip2,           nproc_group_global)
call comm_bcast(iflag_intelectron,    nproc_group_global)
call comm_bcast(num_dip2,             nproc_group_global)
call comm_bcast(dip2boundary,         nproc_group_global)
dip2boundary = dip2boundary*ulength_to_au
call comm_bcast(dip2center,           nproc_group_global)
dip2center = dip2center*ulength_to_au
call comm_bcast(iflag_fourier_omega,  nproc_group_global)
call comm_bcast(num_fourier_omega,    nproc_group_global)
call comm_bcast(fourier_omega,        nproc_group_global)
fourier_omega = fourier_omega*uenergy_to_au
call comm_bcast(itotNtime2,           nproc_group_global)
call comm_bcast(iwdenoption,          nproc_group_global)
call comm_bcast(iwdenstep,            nproc_group_global)
call comm_bcast(iflag_Estatic,        nproc_group_global)

if(iflag_dip2==0)then
  allocate(rto(1))
  allocate(idip2int(1))
else if(iflag_dip2==1)then
  allocate(rto(1:num_dip2-1))
  allocate(idip2int(1:num_dip2-1))

!  dip2center(:)=dip2center(:)/a_B
end if

if(comm_is_root(nproc_id_global))then
  if(iwdenoption/=0.and.iwdenoption/=1)then
    write(*,*)  'iwdenoption must be equal to 0 or 1.'
    stop
  end if
end if
if(iwdenoption==0)then
  iwdenstep=0
end if

if(comm_is_root(nproc_id_global))close(fh_namelist)

end subroutine read_input_rt
