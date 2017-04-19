!
!  Copyright 2016 ARTED developers
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

!===============================================================
subroutine write_result(index)
  use Global_Variables
  use communication, only: procid, nprocs
  implicit none
  integer, intent(in) :: index
  integer :: iter,ix_m,iy_m
  ! export results of each calculation step
  ! (written by M.Uemoto on 2016-11-22)
  iter = Nstep_write * (nprocs(1) * index + procid(1))
  write(file_ac, "(A,A,'_Ac_',I6.6,'.out')") trim(process_directory), trim(SYSname), iter
  open(902, file=file_ac)
  select case(FDTDdim)
  case("1D")
    write(902,*) "# X Acx Acy Acz Ex Ey Ez Bx By Bz Jx Jy Jz E(EM) E(J) E(Mat) E(Tot)"
    do ix_m=NXvacL_m,NXvacR_m
      write(902,'(17e26.16E3)') ix_m*HX_m, data_out(1:16,ix_m,1,index)
    end do
  case("2D", "2DC")
    write(902,*) "# X Y Acx Acy Acz Ex Ey Ez Bx By Bz Jx Jy Jz E_EM E_J E_Mat E_Tot"
    do iy_m=NYvacB_m,NYvacT_m
      do ix_m=NXvacL_m,NXvacR_m
        write(902,'(18e26.16E3)') ix_m*HX_m,iy_m*HY_m,data_out(1:16,ix_m,iy_m,index)
      end do
    end do
  end select
  close(902)
  return
end subroutine write_result
!===============================================================
subroutine write_result_all()
  use Global_Variables
  use communication
  implicit none
  integer :: index, n
  
  do index = 0, Ndata_out_per_proc
    n = nprocs(1) * index + procid(1)
    if (n <= Ndata_out) then
      call write_result(index)
    end if
  end do
  return
end subroutine write_result_all
!===============================================================
subroutine write_energy(iter)
  use Global_Variables
  implicit none
  integer iter,ix_m,iy_m
  character(30) wf,citer
  
  write(citer,'(I6.6)')iter
  wf=trim(directory)//trim(SYSname) &
       &//'_energy_'//trim(citer)//'.out'
  open(903,file=wf)

  if(NY_m==1) then
    do ix_m=NXvacL_m,NXvacR_m
      write(903,'(4e26.16E3)') ix_m*HX_m,energy_elemag(ix_m,1),energy_elec(ix_m,1) &
        &,energy_total(ix_m,1)
    end do
  else
    do iy_m=1,NY_m
      do ix_m=NXvacL_m,NXvacR_m
        write(903,'(5e26.16E3)') ix_m*HX_m,iy_m*HY_m,energy_elemag(ix_m,iy_m),energy_elec(ix_m,iy_m) &
          &,energy_total(ix_m,iy_m)
      end do
    end do
  end if
  close(903)
! 1000 format(1x,2(1pe10.3,1x))
  return
end subroutine write_energy
!===============================================================
subroutine write_excited_electron(iter)
  use Global_Variables
  implicit none
  integer iter,ix_m,iy_m
  character(30) wf,citer
  
  write(citer,'(I6.6)')iter
  wf=trim(directory)//trim(SYSname) &
       &//'_exc_elec_'//trim(citer)//'.out'
  open(903,file=wf)

  if(NY_m==1) then
    do ix_m=1,NX_m
      write(903,'(2e26.16E3)') ix_m*HX_m,excited_electron(ix_m,1)
    end do
  else
    do iy_m=1,NY_m
      do ix_m=1,NX_m
        write(903,'(3e26.16E3)') ix_m*HX_m,iy_m*Hy_m,excited_electron(ix_m,iy_m)
      end do
    end do
  end if
  close(903)
! 1000 format(1x,2(1pe10.3,1x))
  return
end subroutine write_excited_electron
!===============================================================
subroutine init_Ac_ms_2dc()
  use Global_variables
  use communication
  implicit none
  ! initialization for 2D cylinder mode
  ! (written by M.Uemoto on 2016-11-22)
  Elec=0d0
  Bmag=0d0
  j_m=0d0
  jmatter_m=0d0
  jmatter_m_l=0d0
  energy_joule=0d0
  select case(AE_shape)
  case('Asin2cos')
    call incident_bessel_beam()
  case('input')
    call read_initial_ac_from_file()
  end select 
  call comm_sync_all
  return
end subroutine init_Ac_ms_2dc
!===============================================================
subroutine init_Ac_ms
  use Global_variables
  use communication
  implicit none
  real(8) x,y
  integer ix_m,iy_m
  real(8) Xstart
  real(8) wpulse_1
  real(8) wpulse_2
! 2D parameter  
  real(8) angle,kabs,kx,ky
  real(8) length_y

  
  call comm_sync_all
  if(comm_is_root())write(*,*)'ok8-1-1'

!  BC_my='isolated'
!  BC_my='periodic'
  f0_1=5.338d-9*sqrt(IWcm2_1)  ! electric field in a.u.
  omega_1=omegaeV_1/(2d0*13.6058d0)  ! frequency in a.u.
  tpulse_1=tpulsefs_1/0.02418d0 ! pulse_duration in a.u.
  Xstart=5*HX_m
  wpulse_1=2*pi/tpulse_1
  f0_2=5.338d-9*sqrt(IWcm2_2)  ! electric field in a.u.
  omega_2=omegaeV_2/(2d0*13.6058d0)  ! frequency in a.u.
  tpulse_2=tpulsefs_2/0.02418d0 ! pulse_duration in a.u.
  wpulse_2=2*pi/tpulse_2
  T1_T2=T1_T2fs/0.02418d0 ! pulse_duration in a.u.

!  aY=40000d0

  
  Elec=0d0
  Bmag=0d0
  g=0d0
  j_m=0d0
  jmatter_m=0d0
  jmatter_m_l=0d0      
  energy_joule=0.0
      
  Ac_m=0d0
  Ac_old_m=0d0
  Ac_new_m=0d0



  call comm_sync_all

  select case(FDTDdim)
  case('1D')
     if(AE_shape == 'Esin2sin') then
     do iy_m=1,NY_m
        y=iy_m*HY_m
        do ix_m=NXvacL_m,0
           x=(ix_m-1)*HX_m
           if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
              Ac_m(3,ix_m,iy_m)=Epdir_1(3)*(-f0_1/(2*omega_1)*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1+wpulse_1))*cos((omega_1+wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1-wpulse_1))*cos((omega_1-wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light))

              
              g(3,ix_m,iy_m)=-Epdir_1(3)*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                   &*sin(omega_1*(x+Xstart+tpulse_1*c_light)/c_light)

              Ac_m(2,ix_m,iy_m)=Epdir_1(2)*(-f0_1/(2*omega_1)*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1+wpulse_1))*cos((omega_1+wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light) &
                   &+f0_1/(4*(omega_1-wpulse_1))*cos((omega_1-wpulse_1)*(x+Xstart+tpulse_1*c_light)/c_light)) 
              
              g(2,ix_m,iy_m)=-Epdir_1(2)*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                   &*sin(omega_1*(x+Xstart+tpulse_1*c_light)/c_light)
           endif
 
          if(x > -Xstart-(tpulse_1+T1_T2)*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light) then
      Ac_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)+Epdir_2(3)*(-f0_2/(2*omega_2)*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2+wpulse_2))*cos((omega_2+wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2-wpulse_2))*cos((omega_2-wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light))

              
              g(3,ix_m,iy_m)=g(3,ix_m,iy_m)-Epdir_2(3)*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                   &*sin(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)

       Ac_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)+Epdir_2(2)*(-f0_2/(2*omega_2)*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2+wpulse_2))*cos((omega_2+wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light) &
                   &+f0_2/(4*(omega_2-wpulse_2))*cos((omega_2-wpulse_2)*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)) 

             
              g(2,ix_m,iy_m)=g(2,ix_m,iy_m)-Epdir_2(2)*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                   &*sin(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light)
           endif
        enddo
     enddo

     do iy_m=1,NY_m
       do ix_m=NXvacL_m+1,NXvacR_m-1
         Ac_new_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)+dt*g(2,ix_m,iy_m)+0.5d0*(c_light*dt/HX_m)**2 &
           *(Ac_m(2,ix_m+1,iy_m)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m-1,iy_m))+0.5d0*(c_light*dt/HY_m)**2 &
           *(Ac_m(2,ix_m,iy_m+1)-2*Ac_m(2,ix_m,iy_m)+Ac_m(2,ix_m,iy_m-1))-4*pi*dt*dt*j_m(2,ix_m,iy_m) 
         
         Ac_new_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)+dt*g(3,ix_m,iy_m)+0.5d0*(c_light*dt/HX_m)**2 &
           *(Ac_m(3,ix_m+1,iy_m)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m-1,iy_m))+0.5d0*(c_light*dt/HY_m)**2 &
           *(Ac_m(3,ix_m,iy_m+1)-2*Ac_m(3,ix_m,iy_m)+Ac_m(3,ix_m,iy_m-1))-4*pi*dt*dt*j_m(3,ix_m,iy_m) 
       end do
     enddo
     
   else if(AE_shape == 'Asin2cos') then
     do iy_m=1,NY_m
       y=iy_m*HY_m
       do ix_m=NXvacL_m,0
         x=(ix_m-1)*HX_m
         if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
           Ac_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
           Ac_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
         endif
         
         if(x > -Xstart-(tpulse_1+T1_T2)*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light ) then
           
    Ac_m(3,ix_m,iy_m)=Ac_m(3,ix_m,iy_m)-Epdir_2(3)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
       &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
          
    Ac_m(2,ix_m,iy_m)=Ac_m(2,ix_m,iy_m)-Epdir_2(2)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
                  &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
         endif
         
         x=x-dt*c_light
         if(x > -Xstart-tpulse_1*c_light+dt*c_light .and. x < -Xstart+dt*c_light) then
           Ac_new_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
           
           Ac_new_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
             &*cos(omega_1*(x+Xstart+tpulse_1*c_light)/c_light+phi_CEP_1*2d0*pi)
              
         endif
         
         if(x > -Xstart-(tpulse_1+T1_T2)*c_light+dt*c_light .and. x < -Xstart-(tpulse_1+T1_T2-tpulse_2)*c_light+dt*c_light ) then
    Ac_new_m(3,ix_m,iy_m)=Ac_new_m(3,ix_m,iy_m)&
             &-Epdir_2(3)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
             &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
           
    Ac_new_m(2,ix_m,iy_m)=Ac_new_m(2,ix_m,iy_m)&
             &-Epdir_2(2)/omega_2*f0_2*sin(pi*(x+Xstart+(tpulse_1+T1_T2)*c_light)/(tpulse_2*c_light))**2 &
             &*cos(omega_2*(x+Xstart+(tpulse_1+T1_T2)*c_light)/c_light+phi_CEP_2*2d0*pi)
         endif
       enddo
     enddo
   end if

 case('2D')

   angle=45d0*pi/180d0
   kabs=omega_1/c_light
   kx=kabs*cos(angle)
   ky=kabs*sin(angle)
   length_y=2d0*pi/ky
   HY_m=length_y/dble(NY_m)
   do iy_m=1,NY_m
     y=iy_m*HY_m
     do ix_m=NXvacL_m,0
       x=(ix_m-1)*HX_m
       if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
         Ac_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y)
         
         Ac_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y)
       endif
       if(x > -Xstart-tpulse_1*c_light .and. x < -Xstart) then
         Ac_new_m(3,ix_m,iy_m)=-Epdir_1(3)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
            &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y-omega_1*dt)            
         
         Ac_new_m(2,ix_m,iy_m)=-Epdir_1(2)/omega_1*f0_1*sin(pi*(x+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
           &*cos(kx*(x+Xstart+tpulse_1*c_light)+ky*y-omega_1*dt)            
       endif
     end do
   end do
 end select




  
  if(comm_is_root()) write(*,*)'ok fdtd 02'
  call comm_sync_all

  select case(TwoD_shape)
  case('periodic')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,NY_m)
    Ac_new_m(:,:,NY_m+1)=Ac_new_m(:,:,1)
    Ac_m(:,:,0)=Ac_m(:,:,NY_m)
    Ac_m(:,:,NY_m+1)=Ac_m(:,:,1)
    g(:,:,0)=g(:,:,NY_m)
    g(:,:,NY_m+1)=g(:,:,1)
  case('isolated')
    Ac_new_m(:,:,0)=Ac_new_m(:,:,1)
    Ac_new_m(:,:,NY_m+1)=0d0
    Ac_m(:,:,0)=Ac_m(:,:,1)
    Ac_m(:,:,NY_m+1)=0d0
    g(:,:,0)=g(:,:,1)
    g(:,:,NY_m+1)=0d0
  case default
    stop 'boundary condition is not good'
  end select
  return
end subroutine init_Ac_ms
!===========================================================
subroutine dt_evolve_Ac_1d
  use Global_variables, only: Ac_old_m, Ac_m, Ac_new_m, &
                            & NXvacL_m,NXvacR_m, HX_m, dt, &
                            & j_m, c_light, pi, c_light  
  implicit none
  integer :: ix_m
  real(8) :: RR(3) ! rot rot Ac

!$omp parallel do default(none) &
!$    private(ix_m) &
!$    shared(NXvacL_m,NXvacR_m,Ac_m,Ac_old_m,Ac_new_m)
  do ix_m=NXvacL_m-1,NXvacR_m+1
    Ac_old_m(:,ix_m,1) = Ac_m    (:,ix_m,1)
    Ac_m    (:,ix_m,1) = Ac_new_m(:,ix_m,1)
  end do
!$omp end parallel do

!$omp parallel do default(none) &
!$    private(ix_m,RR) &
!$    shared(NXvacL_m,NXvacR_m,Ac_m,HX_m,j_m,Ac_old_m,Ac_new_m) &
!$    firstprivate(dt)
  do ix_m=NXvacL_m,NXvacR_m
    RR(1) = 0.0d0
    RR(2) = -(Ac_m(2,ix_m+1,1) - 2*Ac_m(2,ix_m,1) + Ac_m(2,ix_m-1,1)) * (1.0 / HX_m**2) 
    RR(3) = -(Ac_m(3,ix_m+1,1) - 2*Ac_m(3,ix_m,1) + Ac_m(3,ix_m-1,1)) * (1.0 / HX_m**2) 
    Ac_new_m(:,ix_m,1) = (2*Ac_m(:,ix_m,1) - Ac_old_m(:,ix_m,1) &
      & - j_m(:,ix_m,1)*4.0d0*pi*(dt**2) - RR(:)*(c_light*dt)**2 )
  end do
!$omp end parallel do
  return
end subroutine dt_evolve_Ac_1d
!===========================================================
subroutine dt_evolve_Ac_2d
  use Global_variables, only: Ac_old_m, Ac_m, Ac_new_m, &
                            & NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m, TwoD_shape, &
                            & HX_m, HY_m, dt, j_m, c_light, pi, c_light
  implicit none
  integer :: ix_m,iy_m
  real(8) :: RR(3) ! rot rot Ac
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,Ac_m,Ac_old_m,Ac_new_m)
  do iy_m=NYvacB_m-1,NYvacT_m+1
    do ix_m=NXvacL_m-1,NXvacR_m+1
      Ac_old_m(:,ix_m,iy_m) = Ac_m    (:,ix_m,iy_m)
      Ac_m    (:,ix_m,iy_m) = Ac_new_m(:,ix_m,iy_m)
    end do
  end do
!$omp end parallel do

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m,RR) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,HX_m,HY_m,Ac_m,j_m,Ac_old_m,Ac_new_m) &
!$    firstprivate(dt)
  do iy_m=NYvacB_m,NYvacT_m
    do ix_m=NXvacL_m,NXvacR_m
      RR(1) = +(-1.00d0/HY_m**2) * Ac_m(1, ix_m, iy_m-1) &
            & +(+2.00d0/HY_m**2) * Ac_m(1, ix_m, iy_m) &
            & +(-1.00d0/HY_m**2) * Ac_m(1, ix_m, iy_m+1) &
            & +(+0.25d0/HX_m/HY_m) * Ac_m(2, ix_m-1, iy_m-1) &
            & +(-0.25d0/HX_m/HY_m) * Ac_m(2, ix_m-1, iy_m+1) &
            & +(-0.25d0/HX_m/HY_m) * Ac_m(2, ix_m+1, iy_m-1) &
            & +(+0.25d0/HX_m/HY_m) * Ac_m(2, ix_m+1, iy_m+1)
      RR(2) = +(+0.25d0/HX_m/HY_m) * Ac_m(1, ix_m-1, iy_m-1) &
            & +(-0.25d0/HX_m/HY_m) * Ac_m(1, ix_m-1, iy_m+1) &
            & +(-0.25d0/HX_m/HY_m) * Ac_m(1, ix_m+1, iy_m-1) &
            & +(+0.25d0/HX_m/HY_m) * Ac_m(1, ix_m+1, iy_m+1) &
            & +(-1.00d0/HX_m**2) * Ac_m(2, ix_m-1, iy_m) &
            & +(+2.00d0/HX_m**2) * Ac_m(2, ix_m, iy_m) &
            & +(-1.00d0/HX_m**2) * Ac_m(2, ix_m+1, iy_m)
      RR(3) = +(-1.00d0/HX_m**2) * Ac_m(3, ix_m-1, iy_m) &
            & +(-1.00d0/HY_m**2) * Ac_m(3, ix_m, iy_m-1) &
            & +(+2.00d0/HY_m**2 +2.00d0/HX_m**2) * Ac_m(3, ix_m, iy_m) &
            & +(-1.00d0/HY_m**2) * Ac_m(3, ix_m, iy_m+1) &
            & +(-1.00d0/HX_m**2) * Ac_m(3, ix_m+1, iy_m)
      Ac_new_m(:,ix_m,iy_m) = (2 * Ac_m(:,ix_m,iy_m) - Ac_old_m(:,ix_m,iy_m) &
        & -j_m(:,ix_m,iy_m) * 4.0*pi*(dt**2) - RR(:)*(c_light*dt)**2 )
    end do
  end do
!$omp end parallel do

  ! Boundary Condition
  select case(TwoD_shape)
  case('periodic')
!$omp parallel do default(none) &
!$    private(ix_m) &
!$    shared(Ac_new_m,NXvacL_m,NXvacR_m,NYvacB_m,NYvacT_m)
    do ix_m=NXvacL_m-1,NXvacR_m+1
      Ac_new_m(:,ix_m,NYvacB_m-1)=Ac_new_m(:,ix_m,NYvacT_m)
      Ac_new_m(:,ix_m,NYvacT_m+1)=Ac_new_m(:,ix_m,NYvacB_m)
    enddo
!$omp end parallel do
  case('isolated')
!$omp parallel do default(none) &
!$    private(ix_m) &
!$    shared(Ac_new_m,NXvacL_m,NXvacR_m,NYvacB_m,NYvacT_m)
    do ix_m=NXvacL_m-1,NXvacR_m+1
      Ac_new_m(:,ix_m,NYvacB_m-1)=Ac_new_m(:,ix_m,NYvacB_m)
      Ac_new_m(:,ix_m,NYvacT_m+1)=0.0d0
    enddo
!$omp end parallel do
  end select
  return
end subroutine dt_evolve_Ac_2d
!===========================================================
subroutine dt_evolve_Ac_2dc()
  use Global_variables, only: Ac_old_m, Ac_m, Ac_new_m, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, &
                            & HX_m, HY_m, dt, j_m, c_light, pi, c_light
  implicit none
  integer :: ix_m, iy_m
  real(8) :: Y, RR(3) ! rot rot Ac

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,Ac_m,Ac_old_m,Ac_new_m)
  do iy_m=NYvacB_m-1,NYvacT_m+1
    do ix_m=NXvacL_m-1,NXvacR_m+1
      Ac_old_m(:,ix_m,iy_m) = Ac_m    (:,ix_m,iy_m)
      Ac_m    (:,ix_m,iy_m) = Ac_new_m(:,ix_m,iy_m)
    end do
  end do
!$omp end parallel do

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m,RR,Y) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,HX_m,HY_m,Ac_m,j_m,Ac_old_m,Ac_new_m) &
!$    firstprivate(dt)
  do iy_m=NYvacB_m,NYvacT_m
    do ix_m=NXvacL_m,NXvacR_m
      Y = iy_m * HY_m
      RR(1) = +(+0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_m(1,ix_m+0,iy_m-1) &
            & +2.00d0*(1.00d0/HY_m**2)*Ac_m(1,ix_m+0,iy_m+0) &
            & +(-0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_m(1,ix_m+0,iy_m+1) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(2,ix_m-1,iy_m-1) &
            & -0.50d0*(1.00d0/HX_m)*(1.00d0/Y)*Ac_m(2,ix_m-1,iy_m+0) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(2,ix_m-1,iy_m+1) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(2,ix_m+1,iy_m-1) &
            & +0.50d0*(1.00d0/HX_m)*(1.00d0/Y)*Ac_m(2,ix_m+1,iy_m+0) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(2,ix_m+1,iy_m+1)
      RR(2) = +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(1,ix_m-1,iy_m-1) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(1,ix_m-1,iy_m+1) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(1,ix_m+1,iy_m-1) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_m(1,ix_m+1,iy_m+1) &
            & -(1.00d0/HX_m**2)*Ac_m(2,ix_m-1,iy_m+0) &
            & +2.00d0*(1.00d0/HX_m**2)*Ac_m(2,ix_m+0,iy_m+0) &
            & -(1.00d0/HX_m**2)*Ac_m(2,ix_m+1,iy_m+0)
      RR(3) = -(1.00d0/HX_m**2)*Ac_m(3,ix_m-1,iy_m+0) &
            & +(+0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_m(3,ix_m+0,iy_m-1) &
            & +(+(1.00d0/Y**2)+2.00d0*(1.00d0/HX_m**2)+2.00d0*(1.00d0/HY_m**2))*Ac_m(3,ix_m+0,iy_m+0) &
            & +(-0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_m(3,ix_m+0,iy_m+1) &
            & -(1.00d0/HX_m**2)*Ac_m(3,ix_m+1,iy_m+0)
      Ac_new_m(:,ix_m,iy_m) = (2 * Ac_m(:,ix_m,iy_m) - Ac_old_m(:,ix_m,iy_m) &
        & -J_m(:,ix_m,iy_m) * 4.0*pi*(dt**2) - RR(:)*(c_light*dt)**2 )
    end do
  end do
!$omp end parallel do

  ! Boundary condition
!$omp parallel do default(none) &
!$    private(ix_m) &
!$    shared(Ac_new_m,NXvacL_m,NXvacR_m,NYvacB_m,NYvacT_m)
  do ix_m=NXvacL_m-1,NXvacR_m+1
    Ac_new_m(1,ix_m,NYvacB_m-1)=Ac_new_m(1,ix_m,NYvacB_m)
    !!Following BCs are automatically satisfied by adequate initial state.
    !Ac_new_m(2:3,ix_m,NYvacB_m-1)=0.0d0
    !Ac_new_m(:,ix_m,NYvacT_m+1)=0.0d0
  enddo
!$omp end parallel do
  return
end subroutine dt_evolve_Ac_2dc
!===========================================================
subroutine dt_evolve_Ac
  use Global_variables
  use timer
  implicit none
  ! (written by M.Uemoto on 2016-11-22)
  call timer_begin(LOG_DT_EVOLVE_AC)
  
  select case(FDTDdim)
  case('1D')
    call dt_evolve_Ac_1d()
  case('2D')
    call dt_evolve_Ac_2d()
  case('2DC')
    call dt_evolve_Ac_2dc()
  end select
  
  call timer_end(LOG_DT_EVOLVE_AC)
      
  return
end subroutine dt_evolve_Ac
!===========================================================
subroutine calc_elec_field()
  use Global_variables, only: Ac_old_m, Ac_new_m, Elec, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, &
                            & dt
  implicit none
  integer ix_m,iy_m

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,Elec,Ac_new_m,Ac_old_m) &
!$    firstprivate(dt)
  do iy_m=NYvacB_m, NYvacT_m
    do ix_m=NXvacL_m, NXvacR_m
      Elec(:,ix_m,iy_m)=-(Ac_new_m(:,ix_m,iy_m)-Ac_old_m(:,ix_m,iy_m))/(2d0*dt)
    end do
  end do
!$omp end parallel do
end subroutine calc_elec_field
!===========================================================
subroutine calc_bmag_field_1d()
  use Global_variables, only: Ac_m, Bmag, NXvacL_m, NXvacR_m, HX_m, c_light
  implicit none
  integer :: ix_m
  real(8) :: Rc(3)  ! rot Ac
  ! calculate the magnetic field from the vector potential (1D case)
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do default(none) &
!$    private(ix_m,Rc) &
!$    shared(NXvacL_m,NXvacR_m,HX_m,Bmag,Ac_m)
  do ix_m=NXvacL_m, NXvacR_m
    Rc(1) = 0.0d0
    Rc(2) = - (Ac_m(3,ix_m+1,1) - Ac_m(3,ix_m-1,1)) / (2 * HX_m)
    Rc(3) = + (Ac_m(2,ix_m+1,1) - Ac_m(2,ix_m-1,1)) / (2 * HX_m)
    Bmag(:,ix_m,1) = Rc(:) * c_light
  end do
!$omp end parallel do
  return
end subroutine calc_bmag_field_1d
!===========================================================
subroutine calc_bmag_field_2d()
  use Global_variables, only: Ac_m, Bmag, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, &
                            & HX_m, HY_m, c_light
  implicit none
  integer :: ix_m,iy_m
  real(8) :: rc(3)  ! rot Ac
  ! calculate the magnetic field from the vector potential (2D case)
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m,Rc) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,HY_m,HX_m,Bmag,Ac_m)
  do iy_m=NYvacB_m, NYvacT_m
    do ix_m=NXvacL_m, NXvacR_m
      Rc(1) = + (Ac_m(3, ix_m, iy_m+1) - Ac_m(3, ix_m, iy_m-1)) / (2*HY_m)
      Rc(2) = - (Ac_m(3, ix_m+1, iy_m) - Ac_m(3, ix_m-1, iy_m)) / (2*HX_m)
      Rc(3) = + (Ac_m(2, ix_m+1, iy_m) - Ac_m(2, ix_m-1, iy_m)) / (2*HX_m) &
            & - (Ac_m(1, ix_m, iy_m+1) - Ac_m(1, ix_m, iy_m-1)) / (2*HY_m)
      Bmag(:,ix_m,iy_m) = Rc(:) * c_light
    end do
  end do
!$omp end parallel do
  return
end subroutine calc_bmag_field_2d
!===========================================================
subroutine calc_bmag_field_2dc()
  use Global_variables, only: Ac_m, Bmag, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, &
                            & HX_m, HY_m, c_light
  implicit none
  integer :: ix_m,iy_m
  real(8) :: Y, Rc(3)  ! rot Ac
  ! calculate the magnetic field from the vector potential (2D cylindal case)
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m,Rc,Y) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,HY_m,HX_m,Bmag,Ac_m)
  do iy_m=NYvacB_m, NYvacT_m
    do ix_m=NXvacL_m, NXvacR_m
      Y = iy_m * HY_m
      Rc(1) = -0.50d0*(1.00d0/HY_m)*Ac_m(3,ix_m+0,iy_m-1) &
            & +1.00d0*(1.00d0/Y)*Ac_m(3,ix_m+0,iy_m+0) &
            & +0.50d0*(1.00d0/HY_m)*Ac_m(3,ix_m+0,iy_m+1)
      Rc(2) = +0.50d0*(1.00d0/HX_m)*Ac_m(3,ix_m-1,iy_m+0) &
            & -0.50d0*(1.00d0/HX_m)*Ac_m(3,ix_m+1,iy_m+0)
      Rc(3) = +0.50d0*(1.00d0/HY_m)*Ac_m(1,ix_m+0,iy_m-1) &
            & -0.50d0*(1.00d0/HY_m)*Ac_m(1,ix_m+0,iy_m+1) &
            & -0.50d0*(1.00d0/HX_m)*Ac_m(2,ix_m-1,iy_m+0) &
            & +0.50d0*(1.00d0/HX_m)*Ac_m(2,ix_m+1,iy_m+0)
      Bmag(:,ix_m,iy_m) = Rc(:) * c_light
    end do
  end do
!$omp end parallel do

  return
end subroutine calc_bmag_field_2dc
!===========================================================
subroutine calc_bmag_field()
  use Global_variables
  implicit none
  select case(FDTDdim)
  case('1D')
    call calc_bmag_field_1d()
  case('2D')
    call calc_bmag_field_2d()
  case('2DC')
    call calc_bmag_field_2dc()
  end select
  return
end subroutine calc_bmag_field
!===========================================================
subroutine calc_energy_joule()
  use Global_variables, only: energy_joule, Elec, j_m, dt, aLxyz, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, pi
  implicit none
  integer :: ix_m,iy_m
  ! calculate the Ohmic losses in the media
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,energy_joule,j_m,Elec,aLxyz) &
!$    firstprivate(dt)
  do iy_m=NYvacB_m,NYvacT_m
    do ix_m=NXvacL_m,NXvacR_m
      energy_joule(ix_m, iy_m) = energy_joule(ix_m, iy_m) &
          & + sum(-j_m(:,ix_m,iy_m) * Elec(:,ix_m,iy_m)) * dt * aLxyz
    end do
  end do
!$omp end parallel do

  return
end subroutine calc_energy_joule
!===========================================================
subroutine calc_energy_elemag()
  use Global_variables, only: Elec, Bmag, energy_elemag, aLxyz, &
                            & NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, pi
  implicit none
  integer :: ix_m,iy_m
  real(8) :: e2, b2
  ! calculate the total electromagnetic energy
  ! (written by M.Uemoto on 2016-11-22)

!$omp parallel do collapse(2) default(none) &
!$    private(iy_m,ix_m,e2,b2) &
!$    shared(NYvacB_m,NYvacT_m,NXvacL_m,NXvacR_m,Elec,Bmag,energy_elemag,aLxyz)
  do iy_m=NYvacB_m,NYvacT_m
    do ix_m=NXvacL_m,NXvacR_m
      e2 = sum(Elec(:, ix_m, iy_m) ** 2)
      b2 = sum(Bmag(:, ix_m, iy_m) ** 2)
      energy_elemag(ix_m, iy_m) = (1.0 / (8.0 * pi)) * aLxyz * (e2 + b2)
    end do
  end do
!$omp end parallel do

  return
end subroutine calc_energy_elemag


real(8) function calc_pulse_xcenter()
  use Global_variables, only: energy_elemag, NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, HX_m
  implicit none
  integer :: ix_m, iy_m
  real(8) :: x, ex_tot, e_tot

  ex_tot = 0.0
  e_tot = 0.0
!$omp parallel do collapse(2) default(none) &
!$    private(iy_m, ix_m, x) & 
!$    shared(NYvacB_m, NYvacT_m, NXvacL_m, NXvacR_m, HX_m, energy_elemag) &
!$    reduction(+: ex_tot, e_tot)
do iy_m = NYvacB_m, NYvacT_m
  do ix_m = NXvacL_m, NXvacR_m
    x = ix_m * HX_m
    e_tot = e_tot + energy_elemag(ix_m, iy_m)
    ex_tot = ex_tot + x * energy_elemag(ix_m, iy_m)
  end do
end do
!$omp end parallel do
  calc_pulse_xcenter = ex_tot / e_tot
  return
end function calc_pulse_xcenter
