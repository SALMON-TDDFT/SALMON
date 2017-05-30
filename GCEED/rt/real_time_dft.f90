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
MODULE global_variables_rt
use inputoutput
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
use init_sendrecv_sub
use new_world_sub
use Total_Energy_sub
use read_pslfile_sub
use allocate_psl_sub

implicit none

integer       :: Ntime
integer       :: Nenergy
real(8)       :: dE
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

real(8)       :: debye2au   ! [D]  -> [a.u.] 
integer       :: iii

real(8), allocatable :: alpha2(:,:,:,:)

END MODULE global_variables_rt

!=======================================================================

subroutine Real_Time_DFT(nprocs,nprocid)
use global_variables_rt
use allocate_sendrecv_groupob_sub
implicit none

INTEGER :: IC_rt, OC_rt
character(LEN=100) :: file_IN
character(LEN=100) :: file_OUT_rt, file_IN_rt
real(8),allocatable :: alpha_R(:,:),alpha_I(:,:) 
real(8),allocatable :: alphaq_R(:,:,:),alphaq_I(:,:,:) 
real(8),allocatable :: alpha2_R(:,:,:),alpha2_I(:,:,:) 
real(8),allocatable :: alpha2q_R(:,:,:,:),alpha2q_I(:,:,:,:) 
real(8),allocatable :: Dp_box(:,:),alpha_R_box(:,:),alpha_I_box(:,:) 
real(8),allocatable :: Qp_box(:,:,:),alpha_Rq_box(:,:,:),alpha_Iq_box(:,:,:) 
real(8),allocatable :: Sf(:),Sf2(:,:),Sq2(:,:,:)
real(8) :: plz !polarizability
integer :: jj
integer :: iene,nntime,ix,iy,iz
character(100):: timeFile
character(100):: alpha2OutFile
integer :: ia,ib
real(8) :: rab
integer :: nprocs,nprocid

nproc=nprocs
myrank=nprocid

elp3(:)=0.d0
elp5(:)=0.d0
elp3(401)=MPI_Wtime()

iSCFRT=2
OC=0
img=1

iwdenstep=30 
denplane='xy'
idensum=0
posplane=0.d0

inumcpu_check=0

call setbN
call setcN

call read_input_rt(IC_rt,OC_rt,Ntime,Nenergy,dE,file_IN,file_RT,file_alpha,file_RT_q,file_alpha_q,file_RT_e, &
    & file_RT_dip2,file_alpha_dip2,file_RT_dip2_q,file_alpha_dip2_q,file_RT_dip2_e,file_external, &
    & file_IN_rt,file_OUT_rt)

if(myrank.eq.0)then
  write(*,*)
  write(*,*) "Total time step      =",Ntime
  write(*,*) "Time step[fs]        =",dt*au_time_fs
  write(*,*) "Field strength[?]    =",Fst
  if(ikind_eext <= 1)then
    write(*,*) "      direction      =  ",dir
  end if
  write(*,*) "Energy range         =",Nenergy
  write(*,*) "Energy resolution[eV]=",dE*au_energy_ev
  write(*,*) "ikind_eext is           ", ikind_eext
  write(*,*) "Step for writing dens=", iwdenstep
  write(*,*) "Plane showing density=", denplane
  write(*,*) "idensum              =", idensum 
  if(idensum==0) write(*,*) "Position of the plane=", posplane
  select case (ikind_eext)
    case(1,6,7,8,15)
      write(*,'(a21,f5.2,a4)') "Laser frequency     =",       &
                           romega*au_energy_ev, "[eV]"
      write(*,'(a21,f16.8,a4)') "Pulse width of laser=",      &
                           pulse_T*au_time_fs,"[fs]"
      write(*,'(a21,e16.8,a8)') "Laser intensity      =",      &
                           rlaser_I, "[W/cm^2]"
      write(*,'(a21,e16.8,a8)') "tau                  =",      &
                           tau*au_time_fs, "[fs]"
    case(4,12)
      write(*,'(a21,2f5.2,a4)') "Laser frequency     =",       &
                          romega2(1)*au_energy_ev &
                          ,romega2(2)*au_energy_ev, "[eV]"
      write(*,'(a21,2f16.8,a4)') "Pulse width of laser=",      &
                          pulse_T2(1)*au_time_fs&
                          ,pulse_T2(2)*au_time_fs,"[fs]"
      write(*,'(a21,2e16.8,a8)') "Laser intensity      =",      &
                          rlaser_I2(1),rlaser_I2(2), "[W/cm^2]"
      write(*,'(a21,f16.8,a4)') "delay time           =",      &
                          delay*au_time_fs, "[fs]"
      write(*,'(a21,f16.8)') "rcycle                =",rcycle
  end select
  
  if(iflag_dip2 == 1) then
    write(*,'(a21)',advance="no") "dipole boundary      ="
    do jj=1,num_dip2-2
      write(*,'(1e16.8,a8)',advance="no") dip2boundary(jj)*au_length_aa, "[A],"
    end do
    write(*,'(1e16.8,a8)',advance="yes") dip2boundary(num_dip2-1)*au_length_aa, "[A]"
  end if
  
  if(iflag_fourier_omega == 1) then
    write(*,'(a61)') "===== List of frequencies for fourier transform (in eV) ====="
    do jj=1,num_fourier_omega  
      write(*,'(f16.8)') fourier_omega(jj)*au_energy_ev
    end do
    write(*,'(a61)') "============================================================="
  end if

end if

debye2au = 0.393428d0

select case (ikind_eext)
  case(0,10)
    Fst=Fst !/5.14223d1
end select
dE=dE !/2d0/Ry 
dt=dt !*fs2eVinv*2.d0*Ry!a.u. ! 1[fs] = 1.51925 [1/eV]  !2.d0*Ry*1.51925d0

if(idensum==0) posplane=posplane/a_B 

select case (ikind_eext)
  case(1,6:8,11,15)
    romega=romega !/2.d0/Ry 
    pulse_T=pulse_T !*fs2eVinv*2.d0*Ry 
    Fst=sqrt(rlaser_I)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    tau=tau !*fs2eVinv*2.d0*Ry 
    lasbound_sta(1:3)=lasbound_sta(1:3) !/a_B
    lasbound_end(1:3)=lasbound_end(1:3) !/a_B
  case(4,12)
    romega2=romega2 !/2.d0/Ry 
    pulse_T2=pulse_T2 !*fs2eVinv*2.d0*Ry 
    Fst2=sqrt(rlaser_I2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    tau2=tau2 !*fs2eVinv*2.d0*Ry 
    delay  = delay !*fs2eVinv*2.d0*Ry 
end select

if(iflag_fourier_omega==1)then
   fourier_omega(1:num_fourier_omega)=fourier_omega(1:num_fourier_omega) !/2.d0/Ry 
end if

elp3(402)=MPI_Wtime()

! Read SCF data
call IN_data(file_IN)

if(myrank==0)then
  if(icalcforce==1.and.iflag_md==1)then
    do jj=1,2
      if(idisnum(jj)>MI) then
        write(*,*) "idisnum is larger than MI"
        stop
      end if
    end do
  end if
end if

call read_pslfile
call allocate_psl
call init_ps

if(ikind_eext==0.and.icalcforce==0.and.iflag_md==0) call calc_Mps3rd

call init_updown
call init_itype
call init_sendrecv_matrix

call allocate_sendrecv_groupob

if(ilsda==0)then
  numspin=1
else if(ilsda==1)then
  numspin=2
end if

if(MEO==2.or.MEO==3) call make_corr_pole
call make_icoobox_bound
elp3(403)=MPI_Wtime()


if(iflag_dip2==1) then
  if(imesh_oddeven==1)then
    dip2boundary(1:num_dip2-1)=dip2boundary(1:num_dip2-1) !/a_B
    idip2int(1:num_dip2-1)=nint(dip2boundary(1:num_dip2-1)/Hgs(1))
    rto(1:num_dip2-1)=(dip2boundary(1:num_dip2-1)-((dble(idip2int(1:num_dip2-1))-0.5d0)*Hgs(1)))/Hgs(1)
  else if(imesh_oddeven==2)then
    dip2boundary(1:num_dip2-1)=dip2boundary(1:num_dip2-1) !/a_B
    idip2int(1:num_dip2-1)=nint(dip2boundary(1:num_dip2-1)/Hgs(1)+0.5d0)
    rto(1:num_dip2-1)=(dip2boundary(1:num_dip2-1)-((dble(idip2int(1:num_dip2-1))-1.0d0)*Hgs(1)))/Hgs(1)
  end if
end if

if(iflag_fourier_omega==1) then
   allocate(alpha2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3),num_fourier_omega))
end if

Eion=0.d0
do ia=1,MI
do ib=1,ia-1
  rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
           +(Rion(2,ia)-Rion(2,ib))**2      &
           +(Rion(3,ia)-Rion(3,ib))**2)
  Eion=Eion+Zps(Kion(ia))*Zps(Kion(ib))/rab
end do
end do

elp3(404)=MPI_Wtime()

allocate(Ex_fast(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(Ec_fast(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho0(ix,iy,iz) = rho(ix,iy,iz)
end do
end do
end do

allocate( Vh0(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate( Ex_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 
allocate( Ey_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 
allocate( Ez_static(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))) 

!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Ex_static(ix,iy,iz)=0.d0; Ey_static(ix,iy,iz)=0.d0; Ez_static(ix,iy,iz)=0.d0
end do
end do
end do

if(IC_rt==0) then
  allocate( rIe(0:Ntime) )
  allocate( Dp(3,0:Ntime) )
  allocate( Qp(3,3,0:Ntime) )

  if(iflag_dip2==1) then
    allocate( rIe2(0:Ntime,1:num_dip2) ) 
    allocate( Dp2(3,0:Ntime,1:num_dip2) ) 
    allocate( Qp2(3,3,0:Ntime,1:num_dip2) )
  end if
  itotNtime=Ntime
  Miter_rt=0
else if(IC_rt==1) then
  call IN_data_rt(file_IN_rt,IC_rt,Ntime)
end if

elp3(405)=MPI_Wtime()

allocate( alpha_R(3,0:Nenergy), & 
                    alpha_I(3,0:Nenergy), Sf(3) )
allocate( alphaq_R(3,3,0:Nenergy), & 
                    alphaq_I(3,3,0:Nenergy) )

if(iflag_dip2==1)then
  allocate( alpha2_R(3,0:Nenergy,1:num_dip2), & 
                    alpha2_I(3,0:Nenergy,1:num_dip2), Sf2(3,1:num_dip2) )
  allocate( alpha_R_box(3,0:Nenergy), alpha_I_box(3,0:Nenergy) )
  allocate( Dp_box(3,0:Ntime) )

  allocate( alpha2q_R(3,3,0:Nenergy,1:num_dip2), alpha2q_I(3,3,0:Nenergy,1:num_dip2), Sq2(3,3,1:num_dip2) )
  allocate( alpha_Rq_box(3,3,0:Nenergy), alpha_Iq_box(3,3,0:Nenergy) )
  allocate( Qp_box(3,3,0:Ntime) )
end if

ntmg=1
! 'Hartree' parameter

Hconv  = Hconv !/(2d0*Ry)**2d0/a_B**3   ! Convergence criterion
iterVh = 0        ! Iteration counter


if(myrank==0)then
  write(*, *) 
  write(*, *) "dip2boundary", dip2boundary(1), dip2boundary(2)
  write(*, *) "dip2center", dip2center(1), dip2center(2)
  write(*, *) "dip2boundary[A]", dip2boundary(1)*a_B, dip2boundary(2)*a_B
  write(*, *) "dip2center[A]", dip2center(1)*a_B, dip2center(2)*a_B
  write(*, *) 
end if

call Time_Evolution(IC_rt)

elp3(409)=MPI_Wtime()

if(OC_rt==1) call OUT_data_rt(file_OUT_rt)
elp3(410)=MPI_Wtime()


! Output

if(iwrite_external==1)then
  if(myrank==0)then
    open(1,file=file_external)
    if(ikind_eext==1)then
      do nntime=0,itotNtime
        write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
        if(dt*itt <= tau)then
          write(1,'(e16.8)',advance="yes") Fst*sin(romega*dble(nntime)*dt)*sin(Pi*dble(nntime)*dt/pulse_T)**2
        else
          write(1,'(e16.8)',advance="yes") 0.d0
        end if
      end do
    end if
    close(1)
  end if
end if

call Fourier3D(Dp,alpha_R,alpha_I) 
if(iflag_quadrupole==1)then
  do iii=1,3
    call Fourier3D(Qp(iii,:,:),alphaq_R(iii,:,:),alphaq_I(iii,:,:)) 
  end do
end if
if(myrank.eq.0)then
  open(1,file=file_RT)
  write(1,*) "# time[fs],    dipoleMoment(x,y,z)[A?]" 
   do nntime=0,itotNtime
      write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
      write(1,'(3e16.8)',advance="yes") (Dp(iii,nntime)*a_B, iii=1,3)
   end do
  close(1)

  if(iflag_quadrupole==1)then
    open(1,file=file_RT_q)
    write(1,*) "# time[fs],    quadrupoleMoment(xx,yy,zz,xy,yz,zx)[a.u.]" 
    do nntime=0,itotNtime
       write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
       write(1,'(6e16.8)',advance="yes") (Qp(iii,iii,nntime)*a_B**2, iii=1,3), &
           & Qp(1,2,nntime)*a_B**2,Qp(2,3,nntime)*a_B**2,Qp(3,1,nntime)*a_B**2
    end do
    close(1)
  end if

  if(iflag_intelectron==1)then
    open(1,file=file_RT_e)
    write(1,*) "# time[fs],    integrated electron density" 
     do nntime=0,itotNtime
        write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
        write(1,'(e16.8)',advance="yes") rIe(nntime)
     end do
    close(1)
  end if

  if(iflag_dip2==1)then
    open(1,file=file_RT_dip2)
    write(1,*) "# time[fs],    dipoleMoment(x,y,z)[A?]" 
      do nntime=0,itotNtime
        write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
        do jj=1,num_dip2-1
          write(1,'(3e16.8)',advance="no") (Dp2(iii,nntime,jj)*a_B, iii=1,3)
        end do
        write(1,'(3e16.8)',advance="yes") (Dp2(iii,nntime,num_dip2)*a_B, iii=1,3)
      end do
    close(1)

    if(iflag_quadrupole==1)then
      open(1,file=file_RT_dip2_q)
      write(1,*) "# time[fs],    quadrupoleMoment(xx,yy,zz,xy,yz,zx)[A?]" 
        do nntime=0,itotNtime
          write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
          do jj=1,num_dip2-1
            write(1,'(6e16.8)',advance="no") (Qp2(iii,iii,nntime,jj)*a_B**2, iii=1,3), &
                & Qp2(1,2,nntime,jj)*a_B**2,Qp2(2,3,nntime,jj)*a_B**2,Qp2(3,1,nntime,jj)*a_B**2  
          end do
          write(1,'(6e16.8)',advance="yes") (Qp2(iii,iii,nntime,num_dip2)*a_B**2, iii=1,3), &
              & Qp2(1,2,nntime,num_dip2)*a_B**2,Qp2(2,3,nntime,num_dip2)*a_B**2,Qp2(3,1,nntime,num_dip2)*a_B**2
        end do
      close(1)
    end if

    if(iflag_intelectron==1)then
      open(1,file=file_RT_dip2_e)
      write(1,*) "# time[fs],    integrated electron density" 
        do nntime=0,itotNtime
          write(1,'(e13.5)',advance="no") nntime*dt/2.d0/Ry/fs2eVinv
          do jj=1,num_dip2-1
            write(1,'(e16.8)',advance="no") rIe2(nntime,jj)
          end do
          write(1,'(e16.8)',advance="yes") rIe2(nntime,num_dip2)
        end do
      close(1)
    end if
  end if

! Alpha
  open(1,file=file_alpha)
  plz=0.d0
  write(1,*) "# energy[eV], Re[alpha](x,y,z), Im[alpha](x,y,z), S(x,y,z)" 
   do iene=0,Nenergy
      Sf(:)=2*iene*dE/(Pi)*alpha_I(:,iene)
      write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
      write(1,'(3e16.8)',advance="no") (alpha_R(iii,iene)*(a_B)**3, iii=1,3)
      write(1,'(3e16.8)',advance="no") (alpha_I(iii,iene)*(a_B)**3, iii=1,3)
      write(1,'(3e16.8)',advance="yes") (Sf(iii)/2d0/Ry, iii=1,3)
      if(iene.ge.1) plz=plz+Sf(3)/(dble(iene)*dE)**2*dE
   end do
   write(*,'("===== polarizability in z direction =====")')
   write(*,*) "in a.u.:", plz
   write(*,*) "in A^3 :", plz*a_B**3
   write(*,'("=========================================")')
  close(1)

  if(iflag_quadrupole==1)then
    open(1,file=file_alpha_q)
    write(1,*) "# energy[eV], Re[alpha](xx,yy,zz,xy,yz,zx), Im[alpha](xx,yy,zz,xy,yz,zx)" 
     do iene=0,Nenergy
       Sf(:)=2*iene*dE/(Pi)*alpha_I(:,iene)
       write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
       write(1,'(6e16.8)',advance="no") (alphaq_R(iii,iii,iene), iii=1,3),alphaq_R(1,2,iene),alphaq_R(2,3,iene),alphaq_R(3,1,iene)
       write(1,'(6e16.8)',advance="yes") (alphaq_I(iii,iii,iene), iii=1,3),alphaq_I(1,2,iene),alphaq_I(2,3,iene),alphaq_I(3,1,iene)
     end do
    close(1)
  end if

  if(iflag_dip2==1)then
    open(1,file=file_alpha_dip2)
    write(1,*) "# energy[eV], Re[alpha1](x,y,z), Im[alpha1](x,y,z), S1(x,y,z), Re[alpha2](x,y,z), ..."
    do jj=1,num_dip2
      Dp_box(:,:)=Dp2(:,:,jj)
      call Fourier3D(Dp_box,alpha_R_box,alpha_I_box)
      alpha2_R(:,:,jj)=alpha_R_box(:,:)
      alpha2_I(:,:,jj)=alpha_I_box(:,:)
    end do
    do iene=0,Nenergy
      Sf2(1:3,1:num_dip2)=2*iene*dE/(Pi)*alpha2_I(1:3,iene,1:num_dip2)
      write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
      do jj=1,num_dip2-1
        write(1,'(3e16.8)',advance="no") (alpha2_R(iii,iene,jj)*(a_B)**3, iii=1,3)
        write(1,'(3e16.8)',advance="no") (alpha2_I(iii,iene,jj)*(a_B)**3, iii=1,3)
        write(1,'(3e16.8)',advance="no") (Sf2(iii,jj)/2d0/Ry, iii=1,3)
      end do
      write(1,'(3e16.8)',advance="no") (alpha2_R(iii,iene,num_dip2)*(a_B)**3, iii=1,3)
      write(1,'(3e16.8)',advance="no") (alpha2_I(iii,iene,num_dip2)*(a_B)**3, iii=1,3)
      write(1,'(3e16.8)',advance="yes") (Sf2(iii,num_dip2)/2d0/Ry, iii=1,3)
    end do
    close(1)

    if(iflag_quadrupole==1)then
      open(1,file=file_alpha_dip2_q)
      write(1,*) "# energy[eV], Im[alpha1](x,y,z), Im[alpha2](x,y,z), ..."
      do jj=1,num_dip2
        Qp_box(:,:,:)=Qp2(:,:,:,jj)
        do iii=1,3
          call Fourier3D(Qp_box(iii,:,:),alpha_Rq_box(iii,:,:),alpha_Iq_box(iii,:,:)) 
        end do
        alpha2q_R(:,:,:,jj)=alpha_Rq_box(:,:,:)
        alpha2q_I(:,:,:,jj)=alpha_Iq_box(:,:,:)
      end do
      do iene=0,Nenergy
        write(1,'(e13.5)',advance="no") iene*dE*2d0*Ry
        do jj=1,num_dip2-1
          write(1,'(6e16.8)',advance="no") (alpha2q_I(iii,iii,iene,jj), iii=1,3),  &
                                            alpha2q_I(1,2,iene,jj),alpha2q_I(2,3,iene,jj),alpha2q_I(3,1,iene,jj)
        end do
        write(1,'(6e16.8)',advance="yes") (alpha2q_I(iii,iii,iene,num_dip2), iii=1,3), &
            & alpha2q_I(1,2,iene,num_dip2),alpha2q_I(2,3,iene,num_dip2),alpha2q_I(3,1,iene,num_dip2)
      end do
      close(1)
    end if
  end if
end if

if(iflag_fourier_omega==1)then

  call MPI_Allreduce(zalpha2,zalpha3,   &
                     lg_num(1)*lg_num(2)*lg_num(3)*num_fourier_omega,    &
                     MPI_DOUBLE_COMPLEX,MPI_SUM,newworld_comm_h,ierr)

  if(myrank.eq.0)then
    alpha2=real(zalpha3,8)*dt/a_B**3/fs2eVinv/2.d0/Ry
    do jj=1,num_fourier_omega
      write(fileNumber, '(i8)') jj
      alpha2OutFile = trim("fourier3d.")//adjustl(fileNumber)
      open(1,file=alpha2OutFile)
      do iz=lg_sta(3),lg_end(3),1
      do iy=lg_sta(2),lg_end(2),1
      do ix=lg_sta(1),lg_end(1),1
        if(abs(alpha2(ix,iy,iz,jj))>=1.0d-6) then
          write(1,'(e20.8)') alpha2(ix,iy,iz,jj)
        else
          write(1,'(a1)') "0"
        end if
      end do
      end do
      end do
      close(1)
    end do
  end if
end if


elp3(411)=MPI_Wtime()

write(fileNumber, '(i8)') myrank
timeFile = "cputime"//adjustl(fileNumber)
open(79,file=timeFile)

   write(79,'(a)') "==================== elapsed time ===================="
   write(79,'(a,f16.8)') "elapsed time bef. reading lda data [s] = ", elp3(402)-elp3(401)
   write(79,'(a,f16.8)') "elapsed time for reading lda data [s]  = ", elp3(403)-elp3(402)
   write(79,'(a,f16.8)') "elapsed time for allocating zpsi [s]   = ", elp3(404)-elp3(403)
   write(79,'(a,f16.8)') "elapsed time for reading rt data [s]   = ", elp3(405)-elp3(404)
   write(79,'(a,f16.8)') "elapsed time for rt initialization [s] = ", elp3(406)-elp3(405)
   write(79,'(a,f16.8)') "elapsed time for prep. time prop. [s]  = ", elp3(407)-elp3(406)
   write(79,'(a,f16.8)') "elapsed time for prev iterations [s]   = ", elp3(413)-elp3(412)
   write(79,'(a,f16.8)') "elapsed time for rt iterations [s]     = ", elp3(414)-elp3(413)
   write(79,'(a,f16.8)') "elapsed time for aft iterations [s]    = ", elp3(415)-elp3(414)
   write(79,'(a,f16.8)') "elapsed time aft. rt iterations [s]    = ", elp3(409)-elp3(408)
   write(79,'(a,f16.8)') "elapsed time for writing rt data [s]   = ", elp3(410)-elp3(409)
   write(79,'(a,f16.8)') "elapsed time aft. writing rt data [s]  = ", elp3(411)-elp3(410)
   write(79,'(a,f16.8)') "total time [s]                         = ", elp3(411)-elp3(401)
   write(79,'(a)') "======================================================"
   write(79,'(a)') "=========== elapsed time for rt iterations ==========="
   write(79,'(a,f16.8)') "elapsed time for Vbox [s]              = ", elp5(532)
   write(79,'(a,f16.8)') "elapsed time for time propagation [s]  = ", elp5(533)
   write(79,'(a,f16.8)') "elapsed time for calculating rho [s]   = ", elp5(534)
   write(79,'(a,f16.8)') "elapsed time for Allreduce rho [s]     = ", elp5(535)
   write(79,'(a,f16.8)') "elapsed time for Hartree routine [s]   = ", elp5(536)
   write(79,'(a,f16.8)') "elapsed time for Exc_Cor routine [s]   = ", elp5(537)
   write(79,'(a,f16.8)') "elapsed time for bcast Vhxc [s]        = ", elp5(538)
   write(79,'(a,f16.8)') "elapsed time for calculating Dp [s]    = ", elp5(539)
   write(79,'(a,f16.8)') "elapsed time for calculating Etot [s]  = ", elp5(540)
   write(79,'(a,f16.8)') "elapsed time for writing info etc. [s] = ", elp5(541)
   write(79,'(a,f16.8)') "total time for rt iterations [s]       = ", elp5(542)
   write(79,'(a)') "======================================================"
   write(79,'(a)') "======================================================"
   write(79,'(a)') "=========== communication time ======================="
   write(79,'(a,f16.8)') "copy (1) [s]                           = ", elp5(731)
   write(79,'(a,f16.8)') "copy (2) [s]                           = ", elp5(732)
   write(79,'(a,f16.8)') "copy (3) [s]                           = ", elp5(733)
   write(79,'(a,f16.8)') "copy (4) [s]                           = ", elp5(734)
   write(79,'(a,f16.8)') "copy (5) [s]                           = ", elp5(735)
   write(79,'(a,f16.8)') "copy (6) [s]                           = ", elp5(736)
   write(79,'(a,f16.8)') "copy (7) [s]                           = ", elp5(737)
   write(79,'(a,f16.8)') "copy (8) [s]                           = ", elp5(738)
   write(79,'(a,f16.8)') "copy (9) [s]                           = ", elp5(739)
   write(79,'(a,f16.8)') "copy (10) [s]                          = ", elp5(740)
   write(79,'(a,f16.8)') "copy (11) [s]                          = ", elp5(741)
   write(79,'(a,f16.8)') "copy (12) [s]                          = ", elp5(742)
   write(79,'(a,f16.8)') "time for sendrecv [s]                  = ", elp5(743)
   write(79,'(a,f16.8)') "Allreduce for nonlocal [s]             = ", elp5(744)
   write(79,'(a,f16.8)') "Allreduce for rhobox [s]               = ", elp5(760)
   write(79,'(a,f16.8)') "Allreduce in Hartree (1) [s]           = ", elp5(251)
   write(79,'(a,f16.8)') "Allreduce in Hartree (2) [s]           = ", elp5(252)
   write(79,'(a,f16.8)') "Allreduce in Hartree (3) [s]           = ", elp5(253)
   write(79,'(a,f16.8)') "Allreduce in Hartree (4) [s]           = ", elp5(254)
   write(79,'(a,f16.8)') "Allreduce in Hartree (5) [s]           = ", elp5(255)
   write(79,'(a,f16.8)') "Allreduce in Hartree (6) [s]           = ", elp5(256)
   write(79,'(a,f16.8)') "Allgatherv [s]                         = ", elp5(781)
   write(79,'(a,f16.8)') "Allreduce in total_energy_ex (1) [s]   = ", elp5(782)
   write(79,'(a,f16.8)') "Allreduce in total_energy_ex (2) [s]   = ", elp5(783)
   write(79,'(a,f16.8)') "Allreduce in dipole calc. [s]          = ", elp5(784)
   write(79,'(a)') "=========== analysis ================================="
   write(79,'(a,f16.8)') "isend, irecv, wait [s]                 = ", elp5(743)-sum(elp5(731:742))
   write(79,'(a,f16.8)') "Allgatherv [s]                         = ", elp5(781)
   write(79,'(a,f16.8)') "Allreduce [s]                          = ", & 
     elp5(744)+elp5(760)+sum(elp5(251:256))+sum(elp5(782:783))+elp5(784)
   write(79,'(a,f16.8)') "Allreduce (related to num of nodes )   = ", elp5(744)
   write(79,'(a,f16.8)') "Allreduce (not related to num of nodes)= ", &
     elp5(760)+sum(elp5(251:256))+sum(elp5(782:783))+elp5(784)

call deallocate_mat

END subroutine Real_Time_DFT

!=======================================================================

SUBROUTINE Time_Evolution(IC_rt)
use global_variables_rt

implicit none

complex(8),parameter :: zi=(0.d0,1.d0)
integer :: ii,iob,i1,i2,i3,ix,iy,iz,jj,mm
real(8),allocatable :: R1(:,:,:),R2(:,:,:,:)
character(10):: fileLaser
integer:: idensity, idiffDensity, ielf
integer :: IC_rt
integer :: iob_allob
real(8) :: absr2

real(8)    :: rbox_array(10)
real(8)    :: rbox_array2(10)
real(8)    :: rbox_arrayq(3,3)
real(8)    :: rbox_arrayq2(3,3)
real(8)    :: rbox1q,rbox1q12,rbox1q23,rbox1q31

complex(8), allocatable :: shtpsi(:,:,:,:,:)

if(myrank==0.and.iflag_md==1)then
  open(15,file="distance.data")
  if(MI<=9)then
    wmaxMI=MI
  else
    wmaxMI=9
  end if
  do ii=1,wmaxMI
    write(fileNumber, '(i8)') ii
    rtOutFile = "coo"//trim(adjustl(fileNumber))//".data"
    open(20+ii,file=rtOutFile)
    rtOutFile = "force"//trim(adjustl(fileNumber))//".data"
    open(30+ii,file=rtOutFile)
  end do
end if

cumnum=0.d0

idensity=0
idiffDensity=1
ielf=2
fileELF ="ELF"
fileLaser= "laser.out"

allocate (R1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3))) 
if(ikind_eext.ne.0)then
  allocate( Veff(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3))) 
  allocate( Vbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3)))
endif
if(ikind_eext == 4)then
  allocate( R2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3),2)) 
  allocate( Veff2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3),2)) 
endif

allocate( elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2), & 
                              lg_sta(3):lg_end(3))) 

allocate(rhobox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
if(ilsda==1)then
  allocate(rhobox_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2))
end if
if(ilsda==0)then
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox(ix,iy,iz) = 0.d0
  end do
  end do
  end do
  
  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+zpsi_in(ix,iy,iz,iob,1)*conjg(zpsi_in(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
    end do
    end do
    end do
  end do
  call MPI_allreduce(rhobox,rho,      &
                     mg_num(1)*mg_num(2)*mg_num(3),      &
                     MPI_DOUBLE_PRECISION,MPI_SUM,      &
                     newworld_comm_grid,ierr)
else if(ilsda==1)then
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox(ix,iy,iz) = 0.d0
    rhobox_s(ix,iy,iz,1:2) = 0.d0
  end do
  end do
  end do
  
  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
    if(iob_allob<=MST(1))then
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+zpsi_in(ix,iy,iz,iob,1)*conjg(zpsi_in(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
    else
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+zpsi_in(ix,iy,iz,iob,1)*conjg(zpsi_in(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
    end if
  end do
  call MPI_allreduce(rhobox_s,rho_s,      &
                     mg_num(1)*mg_num(2)*mg_num(3)*2,      &
                     MPI_DOUBLE_PRECISION,MPI_SUM,      &
                     newworld_comm_grid,ierr)
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz)=rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
end if

!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho0(ix,iy,iz)=rho(ix,iy,iz)
end do
end do
end do

allocate(zc(N_hamil))

! External Field Direction
select case (ikind_eext)
  case(0:5,7,15)
    lascenter(1:3)=0.d0
  case(6,8)
    lascenter(1:3)=(lasbound_sta(1:3)+lasbound_end(1:3))/2.d0
end select 

if(imesh_oddeven==1)then
  do i1=lg_sta(1),lg_end(1)
    vecR(1,i1,:,:)=dble(i1)-lascenter(1)/Hgs(1)
  end do
  do i2=lg_sta(2),lg_end(2)
    vecR(2,:,i2,:)=dble(i2)-lascenter(2)/Hgs(2)
  end do
  do i3=lg_sta(3),lg_end(3)
    vecR(3,:,:,i3)=dble(i3)-lascenter(3)/Hgs(3)
  end do
else if(imesh_oddeven==2)then
  do i1=lg_sta(1),lg_end(1)
    vecR(1,i1,:,:)=dble(i1)-0.5d0-lascenter(1)/Hgs(1)
  end do
  do i2=lg_sta(2),lg_end(2)
    vecR(2,:,i2,:)=dble(i2)-0.5d0-lascenter(2)/Hgs(2)
  end do
  do i3=lg_sta(3),lg_end(3)
    vecR(3,:,:,i3)=dble(i3)-0.5d0-lascenter(3)/Hgs(3)
  end do
end if


select case (ikind_eext)
  case(0:3,6:8)
    select case(dir)
      case('x') 
        R1(:,:,:) = vecR(1,:,:,:)*Hgs(1)
      case('y') 
        R1(:,:,:) = vecR(2,:,:,:)*Hgs(2)
      case('z') 
        R1(:,:,:) = vecR(3,:,:,:)*Hgs(3)
      case('xx') 
        R1(:,:,:) = vecR(1,:,:,:)*Hgs(1)*vecR(1,:,:,:)*Hgs(1)
      case('yy') 
        R1(:,:,:) = vecR(2,:,:,:)*Hgs(2)*vecR(2,:,:,:)*Hgs(2)
      case('zz') 
        R1(:,:,:) = vecR(3,:,:,:)*Hgs(3)*vecR(3,:,:,:)*Hgs(3)
      case('xy','yx') 
        R1(:,:,:) = vecR(1,:,:,:)*Hgs(1)*vecR(2,:,:,:)*Hgs(2)
      case('yz','zy') 
        R1(:,:,:) = vecR(2,:,:,:)*Hgs(2)*vecR(3,:,:,:)*Hgs(3)
      case('zx','xz') 
        R1(:,:,:) = vecR(3,:,:,:)*Hgs(3)*vecR(1,:,:,:)*Hgs(1)
      case('x+y','y+x') 
        R1(:,:,:) = vecR(1,:,:,:)*Hgs(1)+vecR(2,:,:,:)*Hgs(2)
      case('y+z','z+y') 
        R1(:,:,:) = vecR(2,:,:,:)*Hgs(2)+vecR(3,:,:,:)*Hgs(3)
      case('z+x','x+z') 
        R1(:,:,:) = vecR(3,:,:,:)*Hgs(3)+vecR(1,:,:,:)*Hgs(1)
    end select
  case(4,12)
    select case(dir2)
      case('x+')
        R2(:,:,:,1) = vecR(2,:,:,:)*Hgs(2)
        R2(:,:,:,2) = vecR(3,:,:,:)*Hgs(3)
      case('x-')
        R2(:,:,:,1) = vecR(3,:,:,:)*Hgs(3)
        R2(:,:,:,2) = vecR(2,:,:,:)*Hgs(2)
      case('y+')
        R2(:,:,:,1) = vecR(3,:,:,:)*Hgs(3)
        R2(:,:,:,2) = vecR(1,:,:,:)*Hgs(1)
      case('y-')
        R2(:,:,:,1) = vecR(1,:,:,:)*Hgs(1)
        R2(:,:,:,2) = vecR(3,:,:,:)*Hgs(3)
      case('z+')
        R2(:,:,:,1) = vecR(1,:,:,:)*Hgs(1)
        R2(:,:,:,2) = vecR(2,:,:,:)*Hgs(2)
      case('z-')
        R2(:,:,:,1) = vecR(2,:,:,:)*Hgs(2)
        R2(:,:,:,2) = vecR(1,:,:,:)*Hgs(1)
    end select
end select

if(iflag_dip2==1) then
  allocate(rbox_array_dip2(4,num_dip2),rbox_array2_dip2(4,num_dip2))
  allocate(rbox_array_dip2q(3,3,num_dip2),rbox_array2_dip2q(3,3,num_dip2))
  allocate(rbox_array_dip2e(num_dip2),rbox_array2_dip2e(num_dip2))
  allocate(rto_ix(lg_sta(1):lg_end(1),num_dip2))
  allocate(vecDs2(1:3,1:num_dip2))
  allocate(vecQs2(1:3,1:3,1:num_dip2))

  rto_ix(:,1:num_dip2-1)=0.d0
  rto_ix(:,num_dip2)=1.d0
  do jj=1,num_dip2-1
    do ix=lg_sta(1),lg_end(1)
      if(ix<idip2int(jj))then
        rto_ix(ix,jj)=rto_ix(ix,jj)+1.d0
        rto_ix(ix,jj+1)=rto_ix(ix,jj+1)-1.d0
      else if(ix==idip2int(jj))then
        rto_ix(ix,jj)=rto_ix(ix,jj)+rto(jj)
        rto_ix(ix,jj+1)=rto_ix(ix,jj+1)-rto(jj)
      end if
    end do
  end do
end if

if(IC_rt==0)then
  rbox_array=0.d0
  do i1=1,3
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rbox_array(i1)=rbox_array(i1)+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)
    end do
    end do
    end do
  end do
  
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rbox_array(4)=rbox_array(4)+rho(ix,iy,iz)
  end do
  end do
  end do

  call MPI_allreduce(rbox_array,rbox_array2,4,MPI_DOUBLE_PRECISION,MPI_SUM,      &
           newworld_comm_h,ierr)
  vecDs(1:3)=rbox_array2(1:3)*Hgs(1:3)*Hvol

  if(iflag_quadrupole==1)then
    do i1=1,3
      rbox1q=0.d0
 !$OMP parallel do reduction( + : rbox1q ) private(absr2)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        absr2=vecR(1,ix,iy,iz)**2+vecR(2,ix,iy,iz)**2+vecR(3,ix,iy,iz)**2
        rbox1q=rbox1q+(3.d0*vecR(i1,ix,iy,iz)*vecR(i1,ix,iy,iz)-absr2)*rho(ix,iy,iz)
      end do
      end do
      end do
      rbox_arrayq(i1,i1)=rbox1q
    end do

    rbox1q12=0.d0
    rbox1q23=0.d0
    rbox1q31=0.d0
 !$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 )
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rbox1q12=rbox1q12+3.d0*vecR(1,ix,iy,iz)*vecR(2,ix,iy,iz)*rho(ix,iy,iz)
      rbox1q23=rbox1q23+3.d0*vecR(2,ix,iy,iz)*vecR(3,ix,iy,iz)*rho(ix,iy,iz)
      rbox1q31=rbox1q31+3.d0*vecR(3,ix,iy,iz)*vecR(1,ix,iy,iz)*rho(ix,iy,iz)
    end do
    end do
    end do

    rbox_arrayq(1,2)=rbox1q12 ; rbox_arrayq(2,1)=rbox1q12
    rbox_arrayq(2,3)=rbox1q23 ; rbox_arrayq(3,2)=rbox1q23
    rbox_arrayq(3,1)=rbox1q31 ; rbox_arrayq(1,3)=rbox1q31

    call MPI_allreduce(rbox_arrayq,rbox_arrayq2,9,MPI_DOUBLE_PRECISION,MPI_SUM,      &
             newworld_comm_h,ierr)
    do i1=1,3
      vecQs(1:3,i1)=rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
    end do
  end if

  if(iflag_dip2==1)then
    rbox_array_dip2=0.d0
    do jj=1,num_dip2
      do i1=1,3
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox_array_dip2(i1,jj)=rbox_array_dip2(i1,jj)+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
      end do
    end do

    do jj=1,num_dip2
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        rbox_array_dip2(4,jj)=rbox_array_dip2(4,jj)+rho(ix,iy,iz)*rto_ix(ix,jj)
      end do
      end do
      end do
    end do

    call MPI_allreduce(rbox_array_dip2,rbox_array2_dip2,4*num_dip2,MPI_DOUBLE_PRECISION,MPI_SUM,      &
             newworld_comm_h,ierr)
    do ii=1,num_dip2
      vecDs2(1:3,ii)=rbox_array2_dip2(1:3,ii)*Hgs(1:3)*Hvol
    end do

    if(iflag_quadrupole==1)then
      do jj=1,num_dip2
        vecR_tmp(:,:,:,:)=vecR(:,:,:,:)
        vecR_tmp(1,:,:,:)=vecR_tmp(1,:,:,:)-dip2center(jj)
        do i1=1,3
          rbox1q=0.d0
 !$OMP parallel do reduction( + : rbox1q ) private(absr2)
          do iz=ng_sta(3),ng_end(3)
          do iy=ng_sta(2),ng_end(2)
          do ix=ng_sta(1),ng_end(1)
            absr2=vecR_tmp(1,ix,iy,iz)**2+vecR_tmp(2,ix,iy,iz)**2+vecR_tmp(3,ix,iy,iz)**2
            rbox1q=rbox1q+(3.d0*vecR_tmp(i1,ix,iy,iz)*vecR_tmp(i1,ix,iy,iz)-absr2)*rho(ix,iy,iz)*rto_ix(ix,jj)
          end do
          end do
          end do
          rbox_array_dip2q(i1,i1,jj)=rbox1q
        end do
      end do
        
      do jj=1,num_dip2
        rbox1q12=0.d0
        rbox1q23=0.d0
        rbox1q31=0.d0
 !$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 )
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1q12=rbox1q12+3.d0*vecR_tmp(1,ix,iy,iz)*vecR_tmp(2,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q23=rbox1q23+3.d0*vecR_tmp(2,ix,iy,iz)*vecR_tmp(3,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q31=rbox1q31+3.d0*vecR_tmp(3,ix,iy,iz)*vecR_tmp(1,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2q(1,2,jj)=rbox1q12 ; rbox_array_dip2q(2,1,jj)=rbox1q12
        rbox_array_dip2q(2,3,jj)=rbox1q23 ; rbox_array_dip2q(3,2,jj)=rbox1q23
        rbox_array_dip2q(3,1,jj)=rbox1q31 ; rbox_array_dip2q(1,3,jj)=rbox1q31
      end do

      call MPI_allreduce(rbox_array_dip2q,rbox_array2_dip2q,9*num_dip2,MPI_DOUBLE_PRECISION,MPI_SUM,      &
               newworld_comm_h,ierr)

      do jj=1,num_dip2
        do i1=1,3
          vecQs2(1:3,i1,jj)=rbox_array2_dip2q(1:3,i1,jj)*Hgs(1:3)*Hvol
        end do
      end do
      if (myrank==0)then
        write(*, *) "dip2center maxx", dip2center(2), vecR(1,ng_end(1),ng_end(2),ng_end(3))
        write(*, *) "initial vecQs2", vecQs2(1,1,2)
      end if
    end if

  end if

end if
if(myrank==0)then
  write(*,'(a30)', advance="no") "Static dipole moment(xyz) ="
  write(*,'(3e15.8)') (vecDs(i1)*a_B, i1=1,3)
  write(*,*)
  if(iflag_quadrupole==1)then
    write(*,'(a30)', advance="no") "Static quadrupole moment ="
    write(*,'(6e15.8)') (vecQs(i1,i1), i1=1,3),vecQs(1,2),vecQs(2,3),vecQs(3,1)
    write(*,*)
  end if
endif

! Static dipole moment

select case (ikind_eext)
  case(1,7)
    Veff(:,:,:)=R1(:,:,:)*Fst
  case(4)
    Veff2(:,:,:,1)=R2(:,:,:,1)*Fst2(1)  
    Veff2(:,:,:,2)=R2(:,:,:,2)*Fst2(2)  
  case(6,8)
    Veff(:,:,:)=0.d0
    if(imesh_oddeven==1)then
      ilasbound_sta(1:3)=nint(lasbound_sta(1:3)/Hgs(:))
      ilasbound_end(1:3)=nint(lasbound_end(1:3)/Hgs(:))
    else if(imesh_oddeven==2)then
      ilasbound_sta(1:3)=nint(lasbound_sta(1:3)/Hgs(:)+0.5d0)
      ilasbound_end(1:3)=nint(lasbound_end(1:3)/Hgs(:)+0.5d0)
    end if
    do jj=1,3
      if(ilasbound_sta(jj)<lg_sta(jj))then
        ilasbound_sta(jj)=lg_sta(jj)
      end if
      if(ilasbound_end(jj)>lg_end(jj))then
        ilasbound_end(jj)=lg_end(jj)
      end if
    end do
    do iz=ilasbound_sta(3),ilasbound_end(3)
    do iy=ilasbound_sta(2),ilasbound_end(2)
    do ix=ilasbound_sta(1),ilasbound_end(1)
      Veff(ix,iy,iz)=R1(ix,iy,iz)*Fst  
    end do
    end do
    end do
end select

! Initial wave function
if(IC_rt==0)then
if(iobnum.ge.1)then
  do iob=1,iobnum
    select case (ikind_eext)
      case(0)
        zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
           mg_sta(3):mg_end(3),iob,1)  &
        = exp(zi*Fst*R1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
           mg_sta(3):mg_end(3)))   &
           *  zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                   mg_sta(3):mg_end(3),iob,1) 
      case(3)
        zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  & 
           mg_sta(3):mg_end(3),iob,1)  &
        = exp(zi*Veff(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
           mg_sta(3):mg_end(3)))  &
           *  zpsi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                   mg_sta(3):mg_end(3),iob,1) 
    end select 
  end do
end if
end if

rIe(0)=rbox_array2(4)*Hvol
Dp(:,0)=0.d0
Qp(:,:,0)=0.d0
if(iflag_dip2==1)then
  rIe2(0,:)=rbox_array2_dip2(4,:)*Hvol
  Dp2(:,0,:)=0.d0 
  Qp2(:,:,0,:)=0.d0 
end if


!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Vh0(ix,iy,iz)=Vh(ix,iy,iz)
end do
end do
end do

! calculate Estatic
if(iflag_Estatic==1)then
  call calcEstatic
end if

! write DFT data
! WriteDensity writes rho(xyz) if it recieves "idensity", elseif it
! recieves "idiffDensity", it writes rho-rho0(xyz) to, for example, diff30
  if (iwdenstep /= 0) then
    do itt=0,0
      iSCFRT=2     
!    call calcELF
      write(fileNumber, '(i8)') 0
      rtOutFile = trim(fileTmp)//adjustl(fileNumber)
      rtDiffOutFile = trim(fileTmp2)//adjustl(fileNumber)
      call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                 denplane,idensum,posplane,idensity)
      call WriteDensity(rtDiffOutFile,iwdenoption,iwdenstep,      &
                 denplane,idensum,posplane,idiffDensity)
      if(iflag_Estatic==1)then
        fileTmp3="Exsta"
        rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
        call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                   denplane,idensum,posplane,10)
        fileTmp3="Eysta"
        rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
        call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                 denplane,idensum,posplane,11)
        fileTmp3="Ezsta"
        rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
        call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                 denplane,idensum,posplane,12)
      end if
    end do
  end if

allocate (shtpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,   &
                 1:iobnum,1))

do iob=1,iobnum
!$OMP parallel do
  do iz=mg_sta(3)-Nd,mg_end(3)+Nd
  do iy=mg_sta(2)-Nd,mg_end(2)+Nd
  do ix=mg_sta(1)-Nd,mg_end(1)+Nd+1
    shtpsi(ix,iy,iz,iob,1)=0.d0
  end do
  end do
  end do
end do

if(iflag_comm_rho==2)then
  allocate(rhobox1_all(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))) 
  allocate(rhobox2_all(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))) 
!$OMP parallel do
  do iz=lg_sta(3),lg_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=lg_sta(1),lg_end(1)
    rhobox1_all(ix,iy,iz) = 0.d0
  end do
  end do
  end do
end if


if(iflag_fourier_omega==1)then
  do mm=1,num_fourier_omega
!$OMP parallel do
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      zalpha2(ix,iy,iz,mm)=0.d0
    end do
    end do
    end do
  end do
end if

!-------------------------------------------------- Time evolution
if(iflag_md==1)then
  call calc_force_c(zpsi_in)
end if

if(myrank.eq.0)then
  write(*,'(1x,a10,a10,a25,a15,a25,a10)') " timestep ","time[fs]",      &
                           " Dipole moment(xyz)[A]"      &
        ,"      electrons","      Total energy[eV]","   iterVh"
  write(*,*) "-------------------------------------"      &
     ,"------------------"
  if(iflag_md==1)then
    write(15,'(2a16,a24)') "        Time[fs]",      &
                   "    Distance [A]",     &
                   " Total energy [eV]   "  
    write(15,*) "-------------------------------------"      &
     ,"------------------"
    do ii=1,wmaxMI
      write(20+ii,'(a16,a28)') "        Time[fs]",      &
                   "    Cooordinate (xyz) [A]   "
      write(20+ii,*) "-------------------------------------"      &
       ,"------------------"
      write(30+ii,'(a16,a28)') "        Time[fs]",      &
                   "       Force (xyz) [eV/A]   "
      write(30+ii,*) "-------------------------------------"      &
       ,"------------------"
    end do
    write(15,'(3f16.8)') dble(0)*dt*0.0241889d0,  &
                sqrt((Rion(1,idisnum(1))-Rion(1,idisnum(2)))**2   &
                    +(Rion(2,idisnum(1))-Rion(2,idisnum(2)))**2   &
                    +(Rion(3,idisnum(1))-Rion(3,idisnum(2)))**2)*a_B, Etot*2.d0*Ry
    do ii=1,wmaxMI
      write(20+ii,'(4f16.8)') dble(0)*dt*0.0241889d0, (Rion(jj,ii)*a_B,jj=1,3)
      write(30+ii,'(4f16.8)') dble(0)*dt*0.0241889d0, (rforce(jj,ii)*2.d0*Ry/a_B,jj=1,3)
    end do
  end if
  if(iwrite_projection==1)then
    open(41,file=file_Projection)
    write(41,'("#",a13,a56)') "time[fs]", "    projection    projection    projection    projection" 
    write(41,'("#",13x,a9,i5,a9,i5,a9,i5,a9,i5)') " orbital",iwrite_projection_ob(1),&
                                              " orbital",iwrite_projection_ob(2),&
                                              " orbital",iwrite_projection_ob(3),&
                                              " orbital",iwrite_projection_ob(4)
    write(41,'("#",13x,a9,i5,a9,i5,a9,i5,a9,i5)') "k",iwrite_projection_k(1),&
                                              "k",iwrite_projection_k(2),&
                                              "k",iwrite_projection_k(3),&
                                              "k",iwrite_projection_k(4)
    write(41,'("#",a)') "---------------------------------------------------------------------"
  end if
end if


elp3(406)=MPI_Wtime()

call taylor_coe

elp3(407)=MPI_Wtime()

if(itotNtime-Miter_rt<=10000)then

  elp3(412)=MPI_Wtime()
  elp3(413)=MPI_Wtime()
  TE : do itt=Miter_rt+1-1,itotNtime
    
    if(iwrite_projection==1.and.itt==Miter_rt+1-1) then
      if(mod(itt,2)==1)then 
        call projection(zpsi_out)
      else
        call projection(zpsi_in)
      end if
    end if

    if(itt>=Miter_rt+1) call time_evolution_step(shtpsi)
  end do TE
  elp3(414)=MPI_Wtime()
  elp3(415)=MPI_Wtime()
  elp5(1:400)=elp3(1:400)
  elp5(431:3000)=elp3(431:3000)


else

  elp3(412)=MPI_Wtime()
  TE1 : do itt=Miter_rt+1-1,Miter_rt+10
    if(iwrite_projection==1.and.itt==Miter_rt+1-1) then
      if(mod(itt,2)==1)then 
        call projection(zpsi_out)
      else
        call projection(zpsi_in)
      end if
    end if

    if(itt>=Miter_rt+1) call time_evolution_step(shtpsi)
  end do TE1
  elp3(413)=MPI_Wtime()

  elp3(1:400)=0.d0
  elp3(431:3000)=0.d0

  TE2 : do itt=Miter_rt+11,itotNtime-5
    call time_evolution_step(shtpsi)
  end do TE2

  elp5(1:400)=elp3(1:400)
  elp5(431:3000)=elp3(431:3000)

  elp3(414)=MPI_Wtime()

  TE3 : do itt=itotNtime-4,itotNtime
    call time_evolution_step(shtpsi)
  end do TE3
  elp3(415)=MPI_Wtime()

end if

elp3(408)=MPI_Wtime()

close(030) ! laser

deallocate (R1)
deallocate (Vlocal)
if(ikind_eext.ne.0)then
  deallocate (Veff, Vbox)
endif
END SUBROUTINE Time_Evolution

!=======================================================================
! Fourier transform for 3D

SUBROUTINE Fourier3D(Dp_t,alpha_R,alpha_I)
use global_variables_rt
implicit none

real(8),intent(IN) :: Dp_t(3,0:Ntime)
real(8),intent(OUT) :: alpha_R(3,0:Nenergy),alpha_I(3,0:Nenergy)
complex(8),parameter   :: zi=(0.d0,1.d0)
complex(8),allocatable :: zalpha(:)
integer :: iene,nntime
real(8) :: t2,hw,TT
allocate(zalpha(3))

! Fourier Transform

TT = dt*itotNtime ! [a.u.]

do iene=0,Nenergy
  hw=iene*dE ; zalpha=(0.d0,0.d0)  ! [a.u.]
  do nntime=1,itotNtime
     t2=nntime*dt ; zalpha(:)=zalpha(:)+exp(zi*hw*t2)*Dp_t(:,nntime) & !hw*t is unitless      
                       *(1-3*(t2/TT)**2+2*(t2/TT)**3)
  end do
  if(ikind_eext==0.or.ikind_eext==10)then
    zalpha=zalpha/Fst*dt
  else
    zalpha=zalpha*dt 
  end if
  alpha_R(:,iene)=real(zalpha(:),8)    ! Real part
  alpha_I(:,iene)=imag(zalpha(:))      ! Imaginary part
end do

deallocate(zalpha)
END SUBROUTINE Fourier3D

