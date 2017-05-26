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

subroutine read_input_rt(IC_rt,OC_rt,Ntime,Nenergy,dE,file_IN,file_RT,file_alpha,file_RT_q,file_alpha_q,file_RT_e, &
    & file_RT_dip2,file_alpha_dip2,file_RT_dip2_q,file_alpha_dip2_q,file_RT_dip2_e,file_external, &
    & file_IN_rt,file_OUT_rt)
use inputoutput
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ii
integer :: IC_rt,OC_rt
integer :: Ntime, Nenergy
real(8) :: dE
character(LEN=100) :: file_IN
character(LEN=100) :: file_OUT_rt, file_IN_rt
character(100) :: file_RT
character(100) :: file_alpha
character(100) :: file_RT_q
character(100) :: file_alpha_q
character(100) :: file_RT_e
character(100) :: file_RT_dip2
character(100) :: file_alpha_dip2
character(100) :: file_RT_dip2_q
character(100) :: file_alpha_dip2_q
character(100) :: file_RT_dip2_e
character(100) :: file_external
namelist / group_fundamental / Nenergy,dE, N_hamil, &
                               icalcforce, iflag_md, idisnum, iwrite_projection, &
                               itwproj, iwrite_projnum, itcalc_ene
namelist / group_parallel / nproc_ob,nproc_Mxin,nproc_Mxin_s,  &
                            isequential, num_datafiles_IN, num_datafiles_OUT, imesh_s_all, iflag_comm_rho
namelist / group_hartree / Hconv, MEO, num_pole_xyz, lmax_MEO
namelist / group_file / IC,IC_rt,OC_rt,file_IN,file_RT,file_alpha,file_RT_q,file_alpha_q,  &
                        file_RT_e,file_RT_dip2,file_alpha_dip2,file_RT_dip2_q,file_alpha_dip2_q, &
    & file_RT_dip2_e,file_IN_rt,file_OUT_rt,fileTmp, fileTmp2, file_Projection
namelist / group_extfield / ikind_eext, Fst, dir, dir2, romega, pulse_T,rlaser_I,tau, &
                            romega2, pulse_T2, rlaser_I2, tau2, delay, rcycle 
namelist / group_propagation / dt, Ntime
namelist / group_others / iparaway_ob,num_projection,iwrite_projection_ob,iwrite_projection_k,  &
                          filename_pot,lasbound_sta,lasbound_end, &
    & iwrite_external,iflag_dip2,iflag_quadrupole,iflag_intelectron,num_dip2, dip2boundary, dip2center,& 
    & iflag_fourier_omega, num_fourier_omega, fourier_omega, itotNtime2, &
    & iwdenoption,iwdenstep, numfile_movie, iflag_Estatic

if(myrank ==0)then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if

!===== namelist for group_fundamental =====
Nenergy=2000
dE=0.01d0
N_hamil=4
icalcforce=0
iflag_md=0
idisnum(1)=1
idisnum(2)=2
iwrite_projection=0
iwrite_projnum=0
itwproj=-1
itcalc_ene=1
if(myrank==0)then
  read(fh_namelist,NML=group_fundamental)
  rewind(fh_namelist)
end if
call MPI_Bcast(Nenergy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(dE,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr); dE = dE*uenergy_to_au
call MPI_Bcast(N_hamil,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(icalcforce,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_md,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(idisnum,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwrite_projection,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwrite_projnum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(itwproj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(itcalc_ene,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

allocate(wtk(1))
wtk(:)=1.d0

if(iwrite_projection==1.and.itwproj==-1)then
  write(*,*) "Please specify itwproj when iwrite_projection=1."
  stop
end if

!===== namelist for group_parallel =====
nproc_ob=0
nproc_Mxin(1:3)=0
nproc_Mxin_s(1:3)=0
isequential=2
num_datafiles_IN=1
num_datafiles_OUT=1
imesh_s_all=1
iflag_comm_rho=1
if(myrank==0)then
  read(fh_namelist,NML=group_parallel)
  rewind(fh_namelist)
  if(isequential<=0.or.isequential>=3)then
    write(*,*) "isequential must be equal to 1 or 2."
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
call MPI_Bcast(iflag_comm_rho,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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

!===== namelist for group_propagation =====
dt=0.d0
Ntime=0
if(myrank==0)then
  read(fh_namelist,NML=group_propagation)
  rewind(fh_namelist)
end if
call MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr); dt=dt*utime_to_au
call MPI_Bcast(Ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(dt<=1.d-10)then
  write(*,*) "please set dt."
  stop
end if
if(Ntime==0)then
  write(*,*) "please set Ntime."
  stop
end if

!===== namelist for group_hartree =====
Hconv=1.d-12
num_pole_xyz(1:3)=1
MEO=2
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
IC=1
IC_rt=0
OC_rt=0
file_IN='file_IN'
file_RT='file_RT'
file_alpha='file_alpha'
file_RT_q='file_RT_q'
file_alpha_q='file_alpha_q'
file_RT_e='file_RT_e'
file_RT_dip2='dip2.data'
file_alpha_dip2='sf2.data'
file_RT_dip2_q='qp2.data'
file_alpha_dip2_q='sfq2.data'
file_RT_dip2_q='ie2.data'
file_external='ext.data'
file_IN_rt='file_IN_rt'
file_OUT_rt='file_OUT_rt'
file_Projection='projection.data'
fileTmp='progress'
fileTmp2='diff'
if(myrank==0)then
  read(fh_namelist,NML=group_file)
  rewind(fh_namelist)
end if
call MPI_Bcast(IC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(IC_rt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(OC_rt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_IN,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_IN_rt,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_OUT_rt,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(file_Projection,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(fileTmp,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(fileTmp2,20,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

if(IC==3.and.num_datafiles_IN/=nproc)then
  if(myrank==0)then
    write(*,*) "num_datafiles_IN is set to nproc."
  end if
  num_datafiles_IN=nproc
end if

!===== namelist for group_extfield =====
ikind_eext=-1
Fst=0.25d0
dir='w'
dir2='w+'
romega=0.d0
pulse_T=0.d0
rlaser_I=0.d0
tau=0.d0
romega2(1:2)=0.d0
pulse_T2(1:2)=0.d0
rlaser_I2(1:2)=0.d0
tau2(1:2)=0.d0
delay=0.d0
rcycle=0.d0
if(myrank==0)then
  read(fh_namelist,NML=group_extfield)
  rewind(fh_namelist)
  if(ikind_eext==-1)then
    write(*,*) "please set ikind_eext."
    stop
  end if
  if(ikind_eext==0)then
    if(dir=='w')then
      write(*,*) "please set dir."
      stop
    end if
  end if
  if(ikind_eext==1)then
    if(romega<=1.d-12)then
      write(*,*) "please set romega."
      stop
    end if
    if(pulse_T<=1.d-12)then
      write(*,*) "please set pulse_T."
      stop
    end if
    if(rlaser_I<=1.d-12)then
      write(*,*) "please set rlaser_I."
      stop
    end if
    if(tau<=1.d-12)then
      write(*,*) "please set tau."
      stop
    end if
    if(dir=='w')then
      write(*,*) "please set dir."
      stop
    end if
  end if
  if(ikind_eext==4)then
    if(romega2(1)<=1.d-12.or.romega2(2)<=1.d-12)then
      write(*,*) "please set romega2."
      stop
    end if
    if(pulse_T2(1)<=1.d-12.or.pulse_T2(2)<=1.d-12)then
      write(*,*) "please set pulse_T2."
      stop
    end if
    if(rlaser_I2(1)<=1.d-12.or.rlaser_I2(2)<=1.d-12)then
      write(*,*) "please set rlaser_I2."
      stop
    end if
    if(tau2(1)<=1.d-12.or.tau2(2)<=1.d-12)then
      write(*,*) "please set tau2."
      stop
    end if
    if(dir2=='w+')then
      write(*,*) "please set dir2."
      stop
    end if
  end if
end if
call MPI_Bcast(ikind_eext,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Fst,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  Fst=Fst*(uenergy_to_au/ulength_to_au)
call MPI_Bcast(dir,3,MPI_Character,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(dir2,2,MPI_Character,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(romega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
romega = romega*uenergy_to_au
call MPI_Bcast(pulse_T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
pulse_T=pulse_T*utime_to_au
call MPI_Bcast(rlaser_I,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
tau_1=tau_1*utime_to_au
call MPI_Bcast(romega2,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
romega2 = romega2*uenergy_to_au
call MPI_Bcast(pulse_T2,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
pulse_T2=pulse_T2*utime_to_au
call MPI_Bcast(rlaser_I2,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(tau2,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
tau_2=tau_2*utime_to_au

call MPI_Bcast(delay,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
delay=delay*utime_to_au
call MPI_Bcast(rcycle,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!===== namelist for group_others =====

iparaway_ob=2
lasbound_sta(:)=-1.d7
lasbound_end(:)=1.d7
num_projection=1
do ii=1,200
  iwrite_projection_ob(ii)=ii
  iwrite_projection_k(ii)=1
end do
filename_pot='pot'
iwrite_external=0
iflag_dip2=0
iflag_quadrupole=0
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
numfile_movie=1
iflag_Estatic=0
if(myrank==0)then
  read(fh_namelist,NML=group_others)
  rewind(fh_namelist)
end if
call MPI_Bcast(iparaway_ob,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_projection,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwrite_projection_ob,200,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwrite_projection_k,200,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(filename_pot,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lasbound_sta,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(lasbound_end,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwrite_external,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_dip2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_quadrupole,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_intelectron,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_dip2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(dip2boundary,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(dip2center,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_fourier_omega,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(num_fourier_omega,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(fourier_omega,200,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(itotNtime2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwdenoption,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iwdenstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(numfile_movie,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(iflag_Estatic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(iflag_dip2==0)then
  allocate(rto(1))
  allocate(idip2int(1))
else if(iflag_dip2==1)then
  allocate(rto(1:num_dip2-1))
  allocate(idip2int(1:num_dip2-1))

  dip2center(:)=dip2center(:)/a_B
end if

if(myrank==0)then
  if(iwdenoption/=0.and.iwdenoption/=1)then
    write(*,*)  'iwdenoption must be equal to 0 or 1.'
    stop
  end if
end if
if(iwdenoption==0)then
  iwdenstep=0
end if

if(myrank ==0)close(fh_namelist)

end subroutine read_input_rt
