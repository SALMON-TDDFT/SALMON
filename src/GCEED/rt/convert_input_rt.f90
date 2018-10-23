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
subroutine convert_input_rt(Ntime)
use salmon_parallel, only: nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_bcast
use inputoutput
use scf_data
use new_world_sub
implicit none
integer :: Ntime
real(8) :: dip_spacing

ik_oddeven=2
ilsda=ispin

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if

!===== namelist for group_fundamental =====
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

num_kpoints_3d(1:3)=num_kgrid(1:3)
num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)

allocate(wtk(num_kpoints_rd))
wtk(:)=1.d0/dble(num_kpoints_rd)

if(iwrite_projection==1.and.itwproj==-1)then
  write(*,*) "Please specify itwproj when iwrite_projection=1."
  stop
end if

!===== namelist for group_parallel =====
if(isequential<=0.or.isequential>=3)then
  write(*,*) "isequential must be equal to 1 or 2."
  stop
end if

nproc_Mxin = nproc_domain
nproc_Mxin_s = nproc_domain_s

if(nproc_ob==0.and.nproc_mxin(1)==0.and.nproc_mxin(2)==0.and.nproc_mxin(3)==0.and.  &
                   nproc_mxin_s(1)==0.and.nproc_mxin_s(2)==0.and.nproc_mxin_s(3)==0) then
  if(ilsda==0)then
    call set_numcpu_rt
  else if(ilsda==1)then
    call set_numcpu_rt_sp
  end if
else
  call check_numcpu
end if

nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

!===== namelist for group_propagation =====
Ntime=nt
if(dt<=1.d-10)then
  write(*,*) "please set dt."
  stop
end if
if(Ntime==0)then
  write(*,*) "please set nt."
  stop
end if

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
if(ic==0)then
  ic=1
end if

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
rlaser_I = rlaser_int_wcm2_1

if(epdir_im1(1)**2+epdir_im1(2)**2+epdir_im1(3)**2+ &
   epdir_im2(1)**2+epdir_im2(2)**2+epdir_im2(3)**2>=1.d-12)then
  circular='y'
else
  circular='n'
end if

!===== namelist for group_others =====

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

select case(trans_longi)
case('tr')
  iflag_indA=0
case('lo')
  iflag_indA=1
end select

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

if(comm_is_root(nproc_id_global))close(fh_namelist)

end subroutine convert_input_rt
