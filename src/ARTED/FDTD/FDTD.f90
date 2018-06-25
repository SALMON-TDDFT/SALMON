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

!===============================================================
subroutine write_excited_electron(iter)
  use Global_Variables
  use salmon_file
  implicit none
  integer iter, imacro, ix_m, iy_m, iz_m, fh
  character(100) file_exc_elec
  
  write(file_exc_elec,'(A, A,"_exc_elec_",I6.6,".out")') trim(directory), trim(SYSname), iter
  fh = open_filehandle(file_exc_elec)
  write(fh, '("#",1X,A)') "imacro ix_m iy_m iz_m excited_electron"
  do imacro = 1, nmacro
    ix_m = macropoint(1, imacro)
    iy_m = macropoint(2, imacro)
    iz_m = macropoint(3, imacro)
    write(fh, '(4I6, e26.16E3)') imacro, ix_m, iy_m, iz_m, excited_electron_new_m(imacro)
  end do
  close(fh)
  return
end subroutine write_excited_electron
!===============================================================
subroutine init_ac_ms_2dc()
  use Global_variables
  use salmon_communication, only: comm_sync_all
  implicit none
  select case(AE_shape1)
  case('Asin2cos')
    call incident_bessel_beam()
  end select 
  call comm_sync_all
  return
end subroutine init_ac_ms_2dc
!===============================================================
subroutine init_ac_ms
  use Global_variables
  use salmon_communication, only: comm_sync_all, comm_is_root
  implicit none
  ! real(8) x,y
  real(8) x
  integer ix_m,iy_m,iz_m
  integer :: icount
  real(8) Xstart
  real(8) wpulse_1
  real(8) wpulse_2
  integer :: npower
! 2D parameter  
  ! real(8) angle,kabs,kx,ky
  ! real(8) length_y
  
  call comm_sync_all
  if(rlaser_int_wcm2_1 < 0d0)then
    f0_1 = amplitude1
  else
    f0_1=5.338d-9*sqrt(rlaser_int_wcm2_1)      ! electric field in a.u.
  end if
!  omega_1=omegaeV_1/(2d0*13.6058d0)  ! frequency in a.u.
!  tpulse_1=tpulsefs_1/0.02418d0 ! pulse_duration in a.u.
  Xstart=5*HX_m
  wpulse_1=2*pi/pulse_tw1
  if(rlaser_int_wcm2_1 < 0d0)then
    f0_2 = amplitude2
  else
    f0_2=5.338d-9*sqrt(rlaser_int_wcm2_2)      ! electric field in a.u.
  end if
!  omega_2=omegaeV_2/(2d0*13.6058d0)  ! frequency in a.u.
!  tpulse_2=tpulsefs_2/0.02418d0 ! pulse_duration in a.u.
  wpulse_2=2*pi/pulse_tw2
!  T1_T2=T1_T2fs/0.02418d0 ! pulse_duration in a.u.

!  aY=40000d0

  call comm_sync_all

  select case(FDTDdim)
  case('1D','1d','2D','2d','3D','3d')
!Pump
     select case(ae_shape1)
     case('Acos2','Acos3','Acos4','Acos6','Acos8')
       select case(ae_shape1)
       case('Acos2'); npower = 2
       case('Acos3'); npower = 3
       case('Acos4'); npower = 4
       case('Acos6'); npower = 6
       case('Acos8'); npower = 8
       case default
         stop 'Error in init_Ac.f90'
       end select

       do iy_m = ny1_m, ny2_m
       do iz_m = nz1_m, nz2_m
       
           do ix_m = nx1_m, nx2_m
              x=(ix_m-1)*HX_m + Xstart + 0.5d0*pulse_tw1*c_light

              if(abs(x) < 0.5d0*pulse_tw1*c_light) then
                 Ac_ms(3,ix_m, iy_m, iz_m)=f0_1/omega1*cos(pi*x/(pulse_tw1*c_light))**npower &
                      *aimag( (epdir_re1(3) + zI*epdir_im1(3)) &
                      *exp(zI*(omega1*x/c_light+phi_CEP1*2d0*pi)))
                 
                 Ac_ms(2,ix_m, iy_m, iz_m)=f0_1/omega1*cos(pi*x/(pulse_tw1*c_light))**npower &
                      *aimag( (epdir_re1(2) + zI*epdir_im1(2)) &
                      *exp(zI*(omega1*x/c_light+phi_CEP1*2d0*pi)))
              endif


              x=x-dt*c_light

              if(abs(x) < 0.5d0*pulse_tw1*c_light) then
                 Ac_new_ms(3,ix_m, iy_m, iz_m)=f0_1/omega1*cos(pi*x/(pulse_tw1*c_light))**npower &
                      *aimag( (epdir_re1(3) + zI*epdir_im1(3)) &
                      *exp(zI*(omega1*x/c_light+phi_CEP1*2d0*pi)))
                 
                 Ac_new_ms(2,ix_m, iy_m, iz_m)=f0_1/omega1*cos(pi*x/(pulse_tw1*c_light))**npower &
                      *aimag( (epdir_re1(2) + zI*epdir_im1(2)) &
                      *exp(zI*(omega1*x/c_light+phi_CEP1*2d0*pi)))
              endif

           end do
       end do
       end do

     case('Asin2cos')
    
       do iy_m = ny1_m, ny2_m
       do iz_m = nz1_m, nz2_m
    
           do ix_m = nx1_m, nx2_m
              x=(ix_m-1)*HX_m

              if(x > -Xstart-pulse_tw1*c_light .and. x < -Xstart) then
                 Ac_ms(3,ix_m, iy_m, iz_m)=-Epdir_Re1(3)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
                      &*cos(omega1*(x+Xstart+pulse_tw1*c_light)/c_light+phi_CEP1*2d0*pi)
                 
                 Ac_ms(2,ix_m, iy_m, iz_m)=-Epdir_Re1(2)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
                      &*cos(omega1*(x+Xstart+pulse_tw1*c_light)/c_light+phi_CEP1*2d0*pi)
                 
              endif


              x=x-dt*c_light
              if(x > -Xstart-pulse_tw1*c_light+dt*c_light .and. x < -Xstart+dt*c_light) then
                 Ac_new_ms(3,ix_m, iy_m, iz_m) = &
                  & -Epdir_Re1(3)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light) &
                  & /(pulse_tw1*c_light))**2 &
                  & *cos(omega1*(x+Xstart+pulse_tw1*c_light)/c_light+phi_CEP1*2d0*pi)
                 
                 Ac_new_ms(2,ix_m, iy_m, iz_m) = &
                  & -Epdir_Re1(2)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light) &
                  & /(pulse_tw1*c_light))**2 &
                  & *cos(omega1*(x+Xstart+pulse_tw1*c_light)/c_light+phi_CEP1*2d0*pi)
              endif
              
           end do
       end do
       end do
         
     case('none')
     case default
        call Err_finalize("Invalid pulse_shape_1 parameter!")
     end select
!Probe         
     select case(ae_shape2)
     case('Acos2','Acos3','Acos4','Acos6','Acos8')
       select case(ae_shape2)
       case('Acos2'); npower = 2
       case('Acos3'); npower = 3
       case('Acos4'); npower = 4
       case('Acos6'); npower = 6
       case('Acos8'); npower = 8
       case default
         stop 'Error in init_Ac.f90'
       end select

       do iy_m = ny1_m, ny2_m
       do iz_m = nz1_m, nz2_m

           do ix_m=nx1_m, nx2_m
              x=(ix_m-1)*HX_m + Xstart + (0.5d0*pulse_tw1 + T1_T2)*c_light

              if(abs(x) < 0.5d0*pulse_tw2*c_light) then
                 Ac_ms(3,ix_m, iy_m, iz_m)=Ac_ms(3,ix_m, iy_m, iz_m) &
                      +f0_2/omega2*cos(pi*x/(pulse_tw2*c_light))**npower &
                      *aimag( (epdir_re2(3) + zI*epdir_im2(3)) &
                      *exp(zI*(omega2*x/c_light+phi_CEP2*2d0*pi)))
                 Ac_ms(2,ix_m, iy_m, iz_m)=Ac_ms(2,ix_m, iy_m, iz_m) &
                      +f0_2/omega2*cos(pi*x/(pulse_tw2*c_light))**npower &
                      *aimag( (epdir_re2(2) + zI*epdir_im2(2)) &
                      *exp(zI*(omega2*x/c_light+phi_CEP2*2d0*pi)))                 
              endif

              x=x-dt*c_light

              if(abs(x) < 0.5d0*pulse_tw2*c_light) then
                 Ac_new_ms(3,ix_m, iy_m, iz_m)=Ac_new_ms(3,ix_m, iy_m, iz_m) &
                      +f0_2/omega2*cos(pi*x/(pulse_tw2*c_light))**npower &
                      *aimag( (epdir_re2(3) + zI*epdir_im2(3)) &
                      *exp(zI*(omega2*x/c_light+phi_CEP2*2d0*pi)))
                 Ac_new_ms(2,ix_m, iy_m, iz_m)=Ac_new_ms(2,ix_m, iy_m, iz_m) &
                      +f0_2/omega2*cos(pi*x/(pulse_tw2*c_light))**npower &
                      *aimag( (epdir_re2(2) + zI*epdir_im2(2)) &
                      *exp(zI*(omega2*x/c_light+phi_CEP2*2d0*pi)))                 
              endif

           end do
           end do
           end do

     case('Asin2cos')
       do iy_m = ny1_m, ny2_m
       do iz_m = nz1_m, nz2_m

           do ix_m = nx1_m, nx2_m
              x=(ix_m-1)*HX_m

              if(x > -Xstart-(pulse_tw1+T1_T2)*c_light .and. x < -Xstart-(pulse_tw1+T1_T2-pulse_tw2)*c_light ) then
                 Ac_ms(3,ix_m, iy_m, iz_m)=Ac_ms(3,ix_m, iy_m, iz_m)-Epdir_Re2(3)/omega2*f0_2 &
                      *sin(pi*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/(pulse_tw2*c_light))**2 &
                      &*cos(omega2*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/c_light+phi_CEP2*2d0*pi)
          
                 Ac_ms(2,ix_m, iy_m, iz_m)=Ac_ms(2,ix_m, iy_m, iz_m)-Epdir_Re2(2)/omega2*f0_2 &
                      &*sin(pi*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/(pulse_tw2*c_light))**2 &
                      &*cos(omega2*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/c_light+phi_CEP2*2d0*pi)
              endif
         

              x=x-dt*c_light         
              if(x > -Xstart-(pulse_tw1+T1_T2)*c_light+dt*c_light &
                   .and. x < -Xstart-(pulse_tw1+T1_T2-pulse_tw2)*c_light+dt*c_light ) then
                 Ac_new_ms(3,ix_m, iy_m, iz_m)=Ac_new_ms(3,ix_m, iy_m, iz_m)&
                      &-Epdir_Re2(3)/omega2*f0_2*sin(pi*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/(pulse_tw2*c_light))**2 &
                      &*cos(omega2*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/c_light+phi_CEP2*2d0*pi)
           
                 Ac_new_ms(2,ix_m, iy_m, iz_m)=Ac_new_ms(2,ix_m, iy_m, iz_m)&
                      &-Epdir_Re2(2)/omega2*f0_2*sin(pi*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/(pulse_tw2*c_light))**2 &
                      &*cos(omega2*(x+Xstart+(pulse_tw1+T1_T2)*c_light)/c_light+phi_CEP2*2d0*pi)
              endif
           enddo
           enddo
           enddo

     case('none')
     case default
        call Err_finalize("Invalid pulse_shape_1 parameter!")
     end select

 ! case('2D', '2d', '3D', '3d')
 ! 
 !   angle=0d0 !45d0*pi/180d0
 !   kabs=omega1/c_light
 !   kx=kabs*cos(angle)
 !   ky=kabs*sin(angle)
 !   length_y=2d0*pi/ky
 !   HY_m=length_y/dble(NY_m)
 !   do iz_m=nz1_m, nz2_m
 !     do iy_m=ny1_m, ny2_m
 !       y=iy_m*HY_m
 !       do ix_m=nx1_m, nx2_m
 !         x=(ix_m-1)*HX_m
 !         if(x > -Xstart-pulse_tw1*c_light .and. x < -Xstart) then
 !           Ac_ms(3,ix_m, iy_m, iz_m)=-Epdir_Re1(3)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
 !             &*cos(kx*(x+Xstart+pulse_tw1*c_light)+ky*y)
 !           
 !           Ac_ms(2,ix_m, iy_m, iz_m)=-Epdir_Re1(2)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
 !             &*cos(kx*(x+Xstart+pulse_tw1*c_light)+ky*y)
 !         endif
 !         if(x > -Xstart-pulse_tw1*c_light .and. x < -Xstart) then
 !           Ac_new_ms(3,ix_m, iy_m, iz_m)=-Epdir_Re1(3)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
 !              &*cos(kx*(x+Xstart+pulse_tw1*c_light)+ky*y-omega1*dt)            
 !           
 !           Ac_new_ms(2,ix_m, iy_m, iz_m)=-Epdir_Re1(2)/omega1*f0_1*sin(pi*(x+Xstart+pulse_tw1*c_light)/(pulse_tw1*c_light))**2 &
 !             &*cos(kx*(x+Xstart+pulse_tw1*c_light)+ky*y-omega1*dt)            
 !         endif
 !       end do
 !     end do
 !  end do
 end select
 
 do icount = 1, ninit_acfield
   ix_m = init_acfield_point(1, icount)
   iy_m = init_acfield_point(2, icount)
   iz_m = init_acfield_point(3, icount)
   Ac_ms(1:3, ix_m, iy_m, iz_m) = init_acfield_val(1:3, icount)
   Ac_new_ms(1:3, ix_m, iy_m, iz_m) = init_acfield_val(4:6, icount)
 end do


  call comm_sync_all

  select case(TwoD_shape)
  case('periodic')
!$omp parallel do collapse(2) default(shared) private(ix_m, iy_m)
  do ix_m = mx1_m, mx2_m
    do iy_m = my1_m, my2_m
      Ac_new_ms(1:3, ix_m, iy_m, mz1_m) =Ac_new_ms(1:3, ix_m, iy_m, nz2_m)
      Ac_new_ms(1:3, ix_m, iy_m, mz2_m) =Ac_new_ms(1:3, ix_m, iy_m, nz1_m)
      Ac_ms(1:3, ix_m, iy_m, mz1_m) =Ac_ms(1:3, ix_m, iy_m, nz2_m)
      Ac_ms(1:3, ix_m, iy_m, mz2_m) =Ac_new_ms(1:3, ix_m, iy_m, nz1_m)
    end do
  end do
!end omp parallel do        
!$omp parallel do collapse(2) default(shared) private(ix_m, iz_m)
  do ix_m = mx1_m, mx2_m
    do iz_m = mz1_m, mz2_m
      Ac_new_ms(1:3, ix_m, my1_m, iz_m) = Ac_new_ms(1:3, ix_m, ny2_m, iz_m) 
      Ac_new_ms(1:3, ix_m, my2_m, iz_m) = Ac_new_ms(1:3, ix_m, ny1_m, iz_m) 
      Ac_ms(1:3, ix_m, my1_m, iz_m) = Ac_ms(1:3, ix_m, ny2_m, iz_m) 
      Ac_ms(1:3, ix_m, my2_m, iz_m) = Ac_ms(1:3, ix_m, ny1_m, iz_m) 
    end do
  end do
!end omp parallel do
!$omp parallel do collapse(2) default(shared) private(iy_m, iz_m)
  do iy_m = my1_m, my2_m
    do iz_m = mz1_m, mz2_m
      Ac_new_ms(1:3, mx1_m, iy_m, iz_m) = Ac_new_ms(1:3, nx2_m, iy_m, iz_m)
      Ac_new_ms(1:3, mx2_m, iy_m, iz_m) = Ac_new_ms(1:3, nx1_m, iy_m, iz_m)
      Ac_ms(1:3, mx1_m, iy_m, iz_m) = Ac_ms(1:3, nx2_m, iy_m, iz_m)
      Ac_ms(1:3, mx2_m, iy_m, iz_m) = Ac_ms(1:3, nx1_m, iy_m, iz_m)
    end do
  end do
!end omp parallel do
    
  case('isolated')
!$omp parallel do default(shared) private(ix_m)
    do ix_m = mx1_m, mx2_m
      Ac_new_ms(:, ix_m, my1_m, iz_m) = Ac_new_ms(:, ix_m, ny1_m, iz_m)
      Ac_new_ms(:, ix_m, my2_m, iz_m) = 0d0
      Ac_ms(:, ix_m, my1_m, iz_m) = Ac_ms(:, ix_m, ny1_m, iz_m)
      Ac_ms(:, ix_m, my2_m, iz_m) = 0d0
    enddo
!$omp end parallel do
  case default
    stop 'boundary condition is not good'
  end select
  return
end subroutine init_ac_ms
!

subroutine add_init_ac_ms
  ! work only for 1D with read_rt_wfn_k_ms=y option
  ! this is just connection from pumpin simulation to probing simulation, 
  ! where vector potential on material at each macro-grid is necessary
  ! to be consistent restarting
  use Global_variables, only:Ac_new_ms,Ac_ms,add_Ac_new_ms,add_Ac_ms
  implicit none

  Ac_new_ms(:,:,:,:) = Ac_new_ms(:,:,:,:) + add_Ac_new_ms(:,:,:,:)
  Ac_ms(:,:,:,:)     = Ac_ms(:,:,:,:)     + add_Ac_ms(:,:,:,:)

  return
end subroutine add_init_ac_ms


!===========================================================


!===========================================================
subroutine dt_evolve_Ac_1d
  use Global_variables
  implicit none
  integer :: ix_m
  integer :: iy_m
  integer :: iz_m
  real(8) :: rr(3) ! rot rot Ac

  iz_m = nz_origin_m
  iy_m = ny_origin_m
!$omp parallel do default(shared) private(ix_m)
  do ix_m = nx1_m, nx2_m
    rr(1) = 0d0
    rr(2:3) = -( &
            &      + Ac_ms(2:3,ix_m+1, iy_m, iz_m) &
            & -2d0 * Ac_ms(2:3,ix_m,   iy_m, iz_m) &
            &      + Ac_ms(2:3,ix_m-1, iy_m, iz_m) &
            & ) * (1d0 / HX_m ** 2)
    Ac_new_ms(:,ix_m, iy_m, iz_m) = (2 * Ac_ms(:,ix_m, iy_m, iz_m) - Ac_old_ms(:,ix_m, iy_m, iz_m) &
      & -Jm_ms(:,ix_m, iy_m, iz_m) * 4.0*pi*(dt**2) - rr(:)*(c_light*dt)**2 )
  end do
!$omp end parallel do


!!(add ion current)
  if(use_ehrenfest_md=='y')then
!$omp parallel do default(shared) private(ix_m)
    do ix_m = nx1_m, nx2_m
       Ac_new_ms(:,ix_m,iy_m,iz_m) = Ac_new_ms(:,ix_m,iy_m,iz_m) &
                                   & - Jm_ion_ms(:,ix_m,iy_m,iz_m)* 4d0*pi*(dt**2)
    end do
!$omp end parallel do
  endif

  return
end subroutine dt_evolve_Ac_1d
!===========================================================
subroutine dt_evolve_Ac_2d
  use Global_variables
  implicit none
  integer :: ix_m,iy_m
  integer :: iz_m
  real(8) :: rr(3) ! rot rot Ac
  iz_m = nz1_m
!$omp parallel do collapse(2) default(shared) private(ix_m, iy_m, rr)
  do iy_m = ny1_m, ny2_m
    do ix_m = nx1_m, nx2_m
      rr(1) = +(-1.00d0/HY_m**2) * Ac_ms(1, ix_m, iy_m-1, iz_m) &
            & +(+2.00d0/HY_m**2) * Ac_ms(1, ix_m, iy_m, iz_m) &
            & +(-1.00d0/HY_m**2) * Ac_ms(1, ix_m, iy_m+1, iz_m) &
            & +(+0.25d0/HX_m/HY_m) * Ac_ms(2, ix_m-1, iy_m-1, iz_m) &
            & +(-0.25d0/HX_m/HY_m) * Ac_ms(2, ix_m-1, iy_m+1, iz_m) &
            & +(-0.25d0/HX_m/HY_m) * Ac_ms(2, ix_m+1, iy_m-1, iz_m) &
            & +(+0.25d0/HX_m/HY_m) * Ac_ms(2, ix_m+1, iy_m+1, iz_m)
      rr(2) = +(+0.25d0/HX_m/HY_m) * Ac_ms(1, ix_m-1, iy_m-1, iz_m) &
            & +(-0.25d0/HX_m/HY_m) * Ac_ms(1, ix_m-1, iy_m+1, iz_m) &
            & +(-0.25d0/HX_m/HY_m) * Ac_ms(1, ix_m+1, iy_m-1, iz_m) &
            & +(+0.25d0/HX_m/HY_m) * Ac_ms(1, ix_m+1, iy_m+1, iz_m) &
            & +(-1.00d0/HX_m**2) * Ac_ms(2, ix_m-1, iy_m, iz_m) &
            & +(+2.00d0/HX_m**2) * Ac_ms(2, ix_m, iy_m, iz_m) &
            & +(-1.00d0/HX_m**2) * Ac_ms(2, ix_m+1, iy_m, iz_m)
      rr(3) = +(-1.00d0/HX_m**2) * Ac_ms(3, ix_m-1, iy_m, iz_m) &
            & +(-1.00d0/HY_m**2) * Ac_ms(3, ix_m, iy_m-1, iz_m) &
            & +(+2.00d0/HY_m**2 +2.00d0/HX_m**2) * Ac_ms(3, ix_m, iy_m, iz_m) &
            & +(-1.00d0/HY_m**2) * Ac_ms(3, ix_m, iy_m+1, iz_m) &
            & +(-1.00d0/HX_m**2) * Ac_ms(3, ix_m+1, iy_m, iz_m)
      Ac_new_ms(:,ix_m, iy_m, iz_m) = (2 * Ac_ms(:,ix_m, iy_m, iz_m) - Ac_old_ms(:,ix_m, iy_m, iz_m) &
        & -Jm_ms(:,ix_m, iy_m, iz_m) * 4.0*pi*(dt**2) - rr(:)*(c_light*dt)**2 )
    end do
  end do
!$omp end parallel do

  ! Boundary Condition
  select case(TwoD_shape)
  case('periodic')
!$omp parallel do default(shared) private(ix_m)
    do ix_m = mx1_m, mx2_m
      Ac_new_ms(:, ix_m, my1_m, iz_m) = Ac_new_ms(:, ix_m, ny2_m, iz_m)
      Ac_new_ms(:, ix_m, my2_m, iz_m) = Ac_new_ms(:, ix_m, ny1_m, iz_m)
    enddo
!$omp end parallel do
  case('isolated')
!$omp parallel do default(shared) private(ix_m)
    do ix_m = mx1_m, mx2_m
      Ac_new_ms(:, ix_m, my1_m, iz_m) = Ac_new_ms(:, ix_m, ny1_m, iz_m)
      Ac_new_ms(:, ix_m, my2_m, iz_m) = 0d0
    enddo
!$omp end parallel do
  case default
    stop 'boundary condition is not good'
  end select
  return
end subroutine dt_evolve_Ac_2d
!===========================================================
subroutine dt_evolve_Ac_2dc()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m
  integer :: iz_m
  real(8) :: Y, rr(3) ! rot rot Ac
  iz_m = nz_origin_m
  
  !! Propagator
!$omp parallel do collapse(2) default(shared) private(ix_m, iy_m, rr,Y)
  do iy_m = ny1_m, ny2_m
    do ix_m = nx1_m, nx2_m
      Y = iy_m * HY_m
      rr(1) = +(+0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_ms(1,ix_m+0,iy_m-1,iz_m) &
            & +2.00d0*(1.00d0/HY_m**2)*Ac_ms(1,ix_m+0,iy_m+0,iz_m) &
            & +(-0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_ms(1,ix_m+0,iy_m+1,iz_m) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(2,ix_m-1,iy_m-1,iz_m) &
            & -0.50d0*(1.00d0/HX_m)*(1.00d0/Y)*Ac_ms(2,ix_m-1,iy_m+0,iz_m) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(2,ix_m-1,iy_m+1,iz_m) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(2,ix_m+1,iy_m-1,iz_m) &
            & +0.50d0*(1.00d0/HX_m)*(1.00d0/Y)*Ac_ms(2,ix_m+1,iy_m+0,iz_m) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(2,ix_m+1,iy_m+1,iz_m)
      rr(2) = +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(1,ix_m-1,iy_m-1,iz_m) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(1,ix_m-1,iy_m+1,iz_m) &
            & -0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(1,ix_m+1,iy_m-1,iz_m) &
            & +0.25d0*(1.00d0/(HX_m*HY_m))*Ac_ms(1,ix_m+1,iy_m+1,iz_m) &
            & -(1.00d0/HX_m**2)*Ac_ms(2,ix_m-1,iy_m+0,iz_m) &
            & +2.00d0*(1.00d0/HX_m**2)*Ac_ms(2,ix_m+0,iy_m+0,iz_m) &
            & -(1.00d0/HX_m**2)*Ac_ms(2,ix_m+1,iy_m+0,iz_m)
      rr(3) = -(1.00d0/HX_m**2)*Ac_ms(3,ix_m-1,iy_m+0,iz_m) &
            & +(+0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_ms(3,ix_m+0,iy_m-1,iz_m) &
            & +(+(1.00d0/Y**2)+2.00d0*(1.00d0/HX_m**2)+2.00d0*(1.00d0/HY_m**2))*Ac_ms(3,ix_m+0,iy_m+0,iz_m) &
            & +(-0.50d0*(1.00d0/HY_m)*(1.00d0/Y)-(1.00d0/HY_m**2))*Ac_ms(3,ix_m+0,iy_m+1,iz_m) &
            & -(1.00d0/HX_m**2)*Ac_ms(3,ix_m+1,iy_m+0,iz_m)
      Ac_new_ms(:,ix_m, iy_m, iz_m) = (2 * Ac_ms(:,ix_m, iy_m, iz_m) - Ac_old_ms(:,ix_m, iy_m, iz_m) &
        & -Jm_ms(:,ix_m, iy_m, iz_m) * 4.0*pi*(dt**2) - rr(:)*(c_light*dt)**2)
    end do
  end do
!$omp end parallel do

  !! Boundary condition
  select case(TwoD_shape)
  case('isolated')
!$omp parallel do default(shared) private(ix_m)
  do ix_m = nx1_m, nx2_m
    Ac_new_ms(:, ix_m, my1_m, iz_m) = 0d0
    Ac_new_ms(:, ix_m, my2_m, iz_m) = 0d0
  enddo
!$omp end parallel do
  end select

  return
end subroutine dt_evolve_Ac_2dc
!===========================================================
subroutine dt_evolve_Ac_3d
  use Global_variables
  implicit none
  integer :: ix_m,iy_m
  integer :: iz_m
  real(8) :: rr(3) ! rot rot Ac
  real(8) :: rinv_dx, rinv_dy, rinv_dz
  
  rinv_dx = 1.00 / HX_m
  rinv_dy = 1.00 / HY_m
  rinv_dz = 1.00 / HZ_m
  
!$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, iz_m,rr)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        ! Calculate Rot Rot A
        rr(1) = - (rinv_dy**2) * Ac_ms(1, ix_m+0, iy_m-1, iz_m+0) &
              & - (rinv_dz**2) * Ac_ms(1, ix_m+0, iy_m+0, iz_m-1) &
              & + (2d0*(rinv_dy**2 + rinv_dz**2)) * Ac_ms(1, ix_m+0, iy_m+0, iz_m+0) &
              & - (rinv_dz**2) * Ac_ms(1, ix_m+0, iy_m+0, iz_m+1) &
              & - (rinv_dy**2) * Ac_ms(1, ix_m+0, iy_m+1, iz_m+0) &
              & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m-1, iy_m-1, iz_m+0) &
              & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m-1, iy_m+1, iz_m+0) &
              & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m+1, iy_m-1, iz_m+0) &
              & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(2, ix_m+1, iy_m+1, iz_m+0) &
              & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m-1, iy_m+0, iz_m-1) &
              & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m-1, iy_m+0, iz_m+1) &
              & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m+1, iy_m+0, iz_m-1) &
              & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(3, ix_m+1, iy_m+0, iz_m+1)
        rr(2) = + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m-1, iy_m-1, iz_m+0) &
              & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m-1, iy_m+1, iz_m+0) &
              & - (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m+1, iy_m-1, iz_m+0) &
              & + (rinv_dx*rinv_dy*0.25d0) * Ac_ms(1, ix_m+1, iy_m+1, iz_m+0) &
              & - (rinv_dx**2) * Ac_ms(2, ix_m-1, iy_m+0, iz_m+0) &
              & - (rinv_dz**2) * Ac_ms(2, ix_m+0, iy_m+0, iz_m-1) &
              & + (2d0*(rinv_dx**2 + rinv_dz**2)) * Ac_ms(2, ix_m+0, iy_m+0, iz_m+0) &
              & - (rinv_dz**2) * Ac_ms(2, ix_m+0, iy_m+0, iz_m+1) &
              & - (rinv_dx**2) * Ac_ms(2, ix_m+1, iy_m+0, iz_m+0) &
              & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m-1, iz_m-1) &
              & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m-1, iz_m+1) &
              & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m+1, iz_m-1) &
              & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(3, ix_m+0, iy_m+1, iz_m+1)
        rr(3) = + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m-1, iy_m+0, iz_m-1) &
              & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m-1, iy_m+0, iz_m+1) &
              & - (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m+1, iy_m+0, iz_m-1) &
              & + (rinv_dx*rinv_dz*0.25d0) * Ac_ms(1, ix_m+1, iy_m+0, iz_m+1) &
              & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m-1, iz_m-1) &
              & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m-1, iz_m+1) &
              & - (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m+1, iz_m-1) &
              & + (rinv_dy*rinv_dz*0.25d0) * Ac_ms(2, ix_m+0, iy_m+1, iz_m+1) &
              & - (rinv_dx**2) * Ac_ms(3, ix_m-1, iy_m+0, iz_m+0) &
              & - (rinv_dy**2) * Ac_ms(3, ix_m+0, iy_m-1, iz_m+0) &
              & + (2d0*(rinv_dx**2 + rinv_dy**2)) * Ac_ms(3, ix_m+0, iy_m+0, iz_m+0) &
              & - (rinv_dy**2) * Ac_ms(3, ix_m+0, iy_m+1, iz_m+0) &
              & - (rinv_dx**2) * Ac_ms(3, ix_m+1, iy_m+0, iz_m+0)
        Ac_new_ms(:,ix_m, iy_m, iz_m) = &
        & + (2 * Ac_ms(:,ix_m, iy_m, iz_m) - Ac_old_ms(:,ix_m, iy_m, iz_m) &
        & - Jm_ms(:,ix_m, iy_m, iz_m) * 4.0*pi*(dt**2) - rr(:)*(c_light*dt)**2 )
      end do
    end do
  end do
  !$omp end parallel do
  
  ! select case(TwoD_shape)
  ! case('periodic')
  !! Only periodic boundary condition is available
  
!$omp parallel do collapse(2) default(shared) private(ix_m, iy_m)
  do ix_m = mx1_m, mx2_m
    do iy_m = my1_m, my2_m
      Ac_new_ms(1:3, ix_m, iy_m, mz1_m) =Ac_new_ms(1:3, ix_m, iy_m, nz2_m)
      Ac_new_ms(1:3, ix_m, iy_m, mz2_m) =Ac_new_ms(1:3, ix_m, iy_m, nz1_m)
    end do
  end do
!end omp parallel do        
!$omp parallel do collapse(2) default(shared) private(ix_m, iz_m)
  do ix_m = mx1_m, mx2_m
    do iz_m = mz1_m, mz2_m
      Ac_new_ms(1:3, ix_m, my1_m, iz_m) = Ac_new_ms(1:3, ix_m, ny2_m, iz_m) 
      Ac_new_ms(1:3, ix_m, my2_m, iz_m) = Ac_new_ms(1:3, ix_m, ny1_m, iz_m) 
    end do
  end do
!end omp parallel do
!$omp parallel do collapse(2) default(shared) private(iy_m, iz_m)
  do iy_m = my1_m, my2_m
    do iz_m = mz1_m, mz2_m
      Ac_new_ms(1:3, mx1_m, iy_m, iz_m) = Ac_new_ms(1:3, nx2_m, iy_m, iz_m)
      Ac_new_ms(1:3, mx2_m, iy_m, iz_m) = Ac_new_ms(1:3, nx1_m, iy_m, iz_m)
    end do
  end do
!end omp parallel do

  ! case default
  !   stop 'boundary condition is not good'
  ! end select
  return
end subroutine dt_evolve_Ac_3d
!===========================================================
subroutine dt_evolve_Ac
  use Global_variables
  use timer
  implicit none
  select case(FDTDdim)
  case('1d', '1D')
    call dt_evolve_Ac_1d()
  case('2d', '2D')
    call dt_evolve_Ac_2d()
  case('2dc', '2DC')
    call dt_evolve_Ac_2dc()
  case('3d', '3D')
    call dt_evolve_Ac_3d()
  end select
  return
end subroutine dt_evolve_Ac
!===========================================================
subroutine calc_elec_field()
  use Global_variables
  implicit none
  integer ix_m, iy_m, iz_m
  !$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, iz_m)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        elec_ms(:, ix_m, iy_m, iz_m) = - (Ac_new_ms(:, ix_m, iy_m, iz_m) - Ac_old_ms(:, ix_m, iy_m, iz_m)) / (2d0 * dt)
      end do
    end do
  end do
  !$omp end parallel do
end subroutine calc_elec_field
!===========================================================
subroutine calc_bmag_field_1d()
  use Global_variables
  implicit none
  integer :: ix_m
  integer :: iy_m
  integer :: iz_m
  real(8) :: Rc(3) 
  iy_m = ny1_m
  iz_m = nz1_m
  !$omp parallel do default(shared) private(ix_m, Rc)
  do ix_m = nx1_m, nx2_m
    Rc(1) = 0.0d0
    Rc(2) = - (Ac_ms(3, ix_m+1, iy_m, iz_m) - Ac_ms(3, ix_m-1, iy_m, iz_m)) / (2 * HX_m)
    Rc(3) = + (Ac_ms(2, ix_m+1, iy_m, iz_m) - Ac_ms(2, ix_m-1, iy_m, iz_m)) / (2 * HX_m)
    bmag_ms(:, ix_m, iy_m, iz_m) = Rc(:) * c_light
  end do
  !$omp end parallel do
  return
end subroutine calc_bmag_field_1d
!===========================================================
subroutine calc_bmag_field_2d()
  use Global_variables
  implicit none
  integer :: ix_m,iy_m
  integer :: iz_m
  real(8) :: rAc(3)  ! rot Ac
  real(8) :: rinv_dx, rinv_dy
  
  rinv_dx = 1.00 / HX_m
  rinv_dy = 1.00 / HY_m
  iz_m = nz1_m
  !$omp parallel do collapse(2) default(shared) private(ix_m, iy_m, rAc)
  do iy_m = ny1_m, ny2_m
    do ix_m = nx1_m, nx2_m
      rAc(1) = +(+Ac_ms(3, ix_m+0, iy_m+1, iz_m+0)-Ac_ms(3, ix_m+0, iy_m-1, iz_m+0)) * (0.5d0 * rinv_dy)
      rAc(2) = -(+Ac_ms(3, ix_m+1, iy_m+0, iz_m+0)-Ac_ms(3, ix_m-1, iy_m+0, iz_m+0)) * (0.5d0 * rinv_dx) 
      rAc(3) = +(+Ac_ms(2, ix_m+1, iy_m+0, iz_m+0)-Ac_ms(2, ix_m-1, iy_m+0, iz_m+0)) * (0.5d0 * rinv_dx) &
             & -(+Ac_ms(1, ix_m+0, iy_m+1, iz_m+0)-Ac_ms(1, ix_m+0, iy_m-1, iz_m+0)) * (0.5d0 * rinv_dy) 
      bmag_ms(:, ix_m, iy_m, iz_m) = rAc(:) * c_light
    end do
  end do
  !$omp end parallel do
  return
end subroutine calc_bmag_field_2d
!===========================================================
subroutine calc_bmag_field_2dc()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: Y, Rc(3)  ! rot Ac
  iz_m = nz_origin_m
  !$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, Y, Rc)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        Y = iy_m * HY_m
        Rc(1) = -0.50d0*(1.00d0/HY_m)*Ac_ms(3,ix_m+0,iy_m-1,iz_m) &
              & +1.00d0*(1.00d0/Y)*Ac_ms(3,ix_m+0,iy_m+0,iz_m) &
              & +0.50d0*(1.00d0/HY_m)*Ac_ms(3,ix_m+0,iy_m+1,iz_m)
        Rc(2) = +0.50d0*(1.00d0/HX_m)*Ac_ms(3,ix_m-1,iy_m+0,iz_m) &
              & -0.50d0*(1.00d0/HX_m)*Ac_ms(3,ix_m+1,iy_m+0,iz_m)
        Rc(3) = +0.50d0*(1.00d0/HY_m)*Ac_ms(1,ix_m+0,iy_m-1,iz_m) &
              & -0.50d0*(1.00d0/HY_m)*Ac_ms(1,ix_m+0,iy_m+1,iz_m) &
              & -0.50d0*(1.00d0/HX_m)*Ac_ms(2,ix_m-1,iy_m+0,iz_m) &
              & +0.50d0*(1.00d0/HX_m)*Ac_ms(2,ix_m+1,iy_m+0,iz_m)
        bmag_ms(:, ix_m, iy_m, iz_m) = Rc(:) * c_light
      end do
    end do
  end do
  !$omp end parallel do
  return
end subroutine calc_bmag_field_2dc
!===========================================================
subroutine calc_bmag_field_3d()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: rac(3)  ! rot Ac
  real(8) :: rinv_dx, rinv_dy, rinv_dz
  
  rinv_dx = 1.00 / HX_m
  rinv_dy = 1.00 / HY_m
  rinv_dz = 1.00 / HZ_m
  !$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, iz_m, rAc)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        rAc(1) = +(+Ac_ms(3, ix_m+0, iy_m+1, iz_m+0)-Ac_ms(3, ix_m+0, iy_m-1, iz_m+0)) * (0.5d0 * rinv_dy) &
               & -(+Ac_ms(2, ix_m+0, iy_m+0, iz_m+1)-Ac_ms(2, ix_m+0, iy_m+0, iz_m-1)) * (0.5d0 * rinv_dz) 
        rAc(2) = +(+Ac_ms(1, ix_m+0, iy_m+0, iz_m+1)-Ac_ms(1, ix_m+0, iy_m+0, iz_m-1)) * (0.5d0 * rinv_dz) &
               & -(+Ac_ms(3, ix_m+1, iy_m+0, iz_m+0)-Ac_ms(3, ix_m-1, iy_m+0, iz_m+0)) * (0.5d0 * rinv_dx) 
        rAc(3) = +(+Ac_ms(2, ix_m+1, iy_m+0, iz_m+0)-Ac_ms(2, ix_m-1, iy_m+0, iz_m+0)) * (0.5d0 * rinv_dx) &
               & -(+Ac_ms(1, ix_m+0, iy_m+1, iz_m+0)-Ac_ms(1, ix_m+0, iy_m-1, iz_m+0)) * (0.5d0 * rinv_dy) 
        bmag_ms(:, ix_m, iy_m, iz_m) = rAc * c_light
      end do
    end do
  end do
  !$omp end parallel do
  return
end subroutine calc_bmag_field_3d
!===========================================================
subroutine calc_bmag_field()
  use Global_variables
  implicit none
  select case(FDTDdim)
  case('1d', '1D')
    call calc_bmag_field_1d()
  case('2d', '2D')
    call calc_bmag_field_2d()
  case('2dc', '2DC')
    call calc_bmag_field_2dc()
  case('3d', '3D')
    call calc_bmag_field_3d()
  end select
  return
end subroutine calc_bmag_field
!===========================================================
subroutine calc_energy_joule()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: elec_mid_old(3), jm_mid_old(3), ohm_mid_old
  !$omp parallel do collapse(3) default(shared) private( &
  !$omp     ix_m, iy_m, iz_m, elec_mid_old, jm_mid_old, ohm_mid_old &
  !$omp     )
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        jm_mid_old = (Jm_ms(:, ix_m, iy_m, iz_m) + Jm_old_ms(:, ix_m, iy_m, iz_m)) * 0.5
        if(use_ehrenfest_md=='y') then
           jm_mid_old = jm_mid_old + (Jm_ion_ms(:,ix_m,iy_m,iz_m) + Jm_ion_old_ms(:,ix_m,iy_m,iz_m))*0.5
        endif
        elec_mid_old = - (Ac_ms(:, ix_m, iy_m, iz_m) - Ac_old_ms(:, ix_m, iy_m, iz_m)) / dt
        ohm_mid_old = sum(-jm_mid_old * elec_mid_old)
        energy_joule_ms(ix_m, iy_m, iz_m) = energy_joule_ms(ix_m, iy_m, iz_m) &
                                          & + ohm_mid_old * aLxyz * dt 
      end do
    end do
  end do
  !$omp end parallel do
  return
end subroutine calc_energy_joule
!===========================================================
subroutine calc_energy_elemag()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: e2, b2
  !$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, iz_m, e2, b2)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        e2 = sum(elec_ms(:, ix_m, iy_m, iz_m) ** 2)
        b2 = sum(bmag_ms(:, ix_m, iy_m, iz_m) ** 2)
        energy_elemag_ms(ix_m, iy_m, iz_m) = (1.0 / (8.0 * pi)) * (e2 + b2) * aLxyz
      end do
    end do
  end do
  !$omp end parallel do
  return
end subroutine calc_energy_elemag
!===========================================================
subroutine calc_total_energy()
  use Global_variables
  implicit none
  integer :: ix_m, iy_m, iz_m
  real(8) :: e_em, e_ej, e_ex
  real(8) :: cell_scalling 
  
  e_em=0d0; e_ej=0d0; e_ex=0d0
  
  select case (FDTDdim)
  case ("1D", "1d")
    cell_scalling = HX_m / aLx
  case ("2D", "2d")
    cell_scalling = HX_m / aLx * HY_m / aLy
  case ("3D", "3d")
    cell_scalling = HX_m / aLx * HY_m / aLy * HZ_m / aLz
  end select
  
!$omp parallel do collapse(3) default(shared) private(ix_m, iy_m, iz_m) reduction(+: e_em, e_ej, e_ex)
  do iz_m = nz1_m, nz2_m
    do iy_m = ny1_m, ny2_m
      do ix_m = nx1_m, nx2_m
        e_em = e_em + energy_elemag_ms(ix_m, iy_m, iz_m) * cell_scalling
        e_ej = e_ej + energy_joule_ms(ix_m, iy_m, iz_m) * cell_scalling
        e_ex = e_ex + energy_elec_ms(ix_m, iy_m, iz_m) * cell_scalling
      end do
    end do
  end do
!$end omp parallel do
  total_energy_elemag_old = total_energy_elemag
  total_energy_absorb_old = total_energy_absorb
  total_energy_em_old = total_energy_em
  total_energy_elec_old = total_energy_elec
  total_energy_elemag = e_em
  total_energy_absorb = e_ej
  total_energy_em = e_em + e_ej
  total_energy_elec = e_ex
  return
end subroutine calc_total_energy

  
