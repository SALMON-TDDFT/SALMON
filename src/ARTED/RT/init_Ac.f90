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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine init_Ac
  use Global_Variables
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_is_root
  use Ac_alocal_laser
  implicit none
  integer :: iter, npower
  real(8) :: tt

  javt = 0d0
  Ac_ext = 0d0
  Ac_ind = 0d0


  if(rlaser_int_wcm2_1 < 0d0)then
    f0_1 = amplitude1
  else
    f0_1=5.338d-9*sqrt(rlaser_int_wcm2_1)      ! electric field in a.u.
  end if
  if(rlaser_int_wcm2_2 < 0d0)then
    f0_2 = amplitude2
  else
    f0_2=5.338d-9*sqrt(rlaser_int_wcm2_2)      ! electric field in a.u.
  end if


!  f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
!  f0_2=5.338d-9*sqrt(IWcm2_2)      ! electric field in a.u.
!  omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
!  tpulse_1=tpulsefs_1/0.02418d0 ! pulse duration in a.u.

!  omega_2=omegaev_2/(2d0*Ry)  ! frequency in a.u.
!  tpulse_2=tpulsefs_2/0.02418d0 ! pulse duration in a.u.
!  T1_T2=T1_T2fs/0.02418d0 ! pulse duration in a.u.
  javt=0.d0
  Ac_ext=0.d0



  select case(AE_shape1)
  case('impulse')
    Ac_ext(:,1)=epdir_re1(1)*e_impulse
    Ac_ext(:,2)=epdir_re1(2)*e_impulse
    Ac_ext(:,3)=epdir_re1(3)*e_impulse
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

   !do iter=0,Nt+1
     !tt=iter*dt - 0.5d0*pulse_tw1
    do iter=-1,Nt+1
      if(iter==-1 .and. t1_delay.ge.0d0) cycle !only for restart of field using rt_wfn_k
      tt=iter*dt - 0.5d0*pulse_tw1 - t1_delay
      if (abs(tt)<0.5d0*pulse_tw1) then
        Ac_ext(iter,:)=-f0_1/omega1*(cos(pi*tt/pulse_tw1))**npower &
          *aimag( (epdir_re1(:) + zI*epdir_im1(:)) &
          *exp(zI*(omega1*tt+phi_CEP1*2d0*pi))  &
          )
      end if
    enddo
    T1_T2 = T1_T2 + t1_delay

  case('Ecos2')
    if(phi_CEP1 /= 0.75d0)then
      call Err_finalize("Error: phi_cep1 should be 0.75 when ae_shape1 is 'Ecos2'.")
    end if
    do iter=0,Nt+1
      tt=iter*dt - 0.5d0*pulse_tw1
      if (abs(tt)<0.5d0*pulse_tw1) then
        Ac_ext(iter,:)=-epdir_re1(:)*f0_1/(8d0*pi**2*omega1 - 2d0*pulse_tw1**2*omega1**3) &
          *( &
          (-4d0*pi**2+pulse_tw1**2*omega1**2 + pulse_tw1**2*omega1**2*cos(2d0*pi*tt/pulse_tw1))*cos(omega1*tt) &
          +2d0*pi*(2d0*pi*cos(pulse_tw1*omega1/2d0) &
          +pulse_tw1*omega1*sin(2d0*pi*tt/pulse_tw1)*sin(omega1*tt)))
      end if
    enddo

  case('Esin2sin')
    do iter=0,Nt+1
      tt=iter*dt
      if (tt<pulse_tw1) then
        Ac_ext(iter,:)=epdir_re1(:)*f0_1*(&
          &-(cos(omega1*tt+phi_CEP1*2d0*pi)-cos(phi_CEP1*2d0*pi))/(2*omega1)&
          &+(cos((omega1+2*Pi/pulse_tw1)*tt+phi_CEP1*2d0*pi)-cos(phi_CEP1*2d0*pi))/(4*(omega1+2*Pi/pulse_tw1))&
          &+(cos((omega1-2*Pi/pulse_tw1)*tt+phi_CEP1*2d0*pi)-cos(phi_CEP1*2d0*pi))/(4*(omega1-2*Pi/pulse_tw1))&
          &-0.5*(cos(omega1*pulse_tw1+phi_CEP1*2d0*pi)-cos(phi_CEP1*2d0*pi))&
          &/(omega1*((omega1*pulse_tw1/(2*pi))**2-1.d0))&
          &*(3.d0*(tt/pulse_tw1)**2-2.d0*(tt/pulse_tw1)**3)&
          &)
      else
        Ac_ext(iter,:)=Ac_ext(iter-1,:)
      endif
    enddo
    
  case('Asin2cos')
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do iter=0,Nt+1
      tt=iter*dt
      if (tt<pulse_tw1) then
        Ac_ext(iter,:)=-epdir_re1(:)*f0_1/omega1*(sin(pi*tt/pulse_tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    enddo
    
  case('input')
    Ac_ext=0d0
    if(comm_is_root(nproc_id_global))then
      open(899,file='input_Ac.dat')
      do iter=0,Nt+1
        read(899,*)Ac_ext(iter,1),Ac_ext(iter,2),Ac_ext(iter,3)
      end do
      close(899)
    end if
    call comm_bcast(Ac_ext,nproc_group_global)
    
  case('Asin2_cw')
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do iter=0,Nt+1
      tt=iter*dt
      if (tt<pulse_tw1*0.5d0) then
        Ac_ext(iter,:)=-Epdir_re1(:)*f0_1/omega1*(sin(pi*tt/pulse_tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      else
        Ac_ext(iter,:)=-Epdir_re1(:)*f0_1/omega1*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    enddo
  case('none')
    !there is no pump
  case default
    call Err_finalize("Invalid pulse_shape_1 parameter!")
  end select
  

! Probe
  select case(ae_shape2)
  case('impulse')
    do iter=0,Nt+1
      tt=iter*dt
      if(tt > T1_T2)then
        Ac_ext(iter,1)=Ac_ext(iter,1) + epdir_re2(1)*e_impulse
        Ac_ext(iter,2)=Ac_ext(iter,2) + epdir_re2(2)*e_impulse
        Ac_ext(iter,3)=Ac_ext(iter,3) + epdir_re2(3)*e_impulse
      end if
    end do
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

    do iter=0,Nt+1
      tt=iter*dt - 0.5d0*pulse_tw1 - T1_T2
      if (abs(tt)<0.5d0*pulse_tw2) then
        Ac_ext(iter,:)=Ac_ext(iter,:) &
          -f0_2/omega2*(cos(pi*tt/pulse_tw2))**npower &
          *aimag( (epdir_re2(:) + zI*epdir_im2(:)) &
          *exp(zI*(omega2*tt+phi_CEP2*2d0*pi))  &
          )
      end if
    end do

  case('Ecos2')
    if(phi_CEP2 /= 0.75d0)then
      call Err_finalize("Error: phi_cep2 should be 0.75 when ae_shape2 is 'Ecos2'.")
    end if
    do iter=0,Nt+1
      tt=iter*dt - 0.5d0*pulse_tw1 - T1_T2
      if (abs(tt)<0.5d0*pulse_tw2) then
        Ac_ext(iter,:)=Ac_ext(iter,:) &
          -epdir_re2(:)*f0_2/(8d0*pi**2*omega2 - 2d0*pulse_tw2**2*omega2**3) &
          *( &
          (-4d0*pi**2+pulse_tw2**2*omega2**2 + pulse_tw2**2*omega2**2*cos(2d0*pi*tt/pulse_tw2))*cos(omega2*tt) &
          +2d0*pi*(2d0*pi*cos(pulse_tw2*omega2/2d0) &
          +pulse_tw2*omega2*sin(2d0*pi*tt/pulse_tw2)*sin(omega2*tt)))
      end if
    enddo

  case('Esin2sin')
      ! probe laser
    do iter=0,Nt+1
      tt=iter*dt
      if(tt-T1_T2 <= 0d0)then
        
      else if ( (tt-T1_T2>0d0) .and. (tt-T1_T2<pulse_tw2) ) then
        Ac_ext(iter,:)=Ac_ext(iter,:)+Epdir_re2(:)*f0_2*(&
          &-(cos(omega2*(tt-T1_T2)+phi_CEP2*2d0*pi)-cos(phi_CEP2*2d0*pi))/(2*omega2)&
          &+(cos((omega2+2*Pi/pulse_tw2)*(tt-T1_T2)+phi_CEP2*2d0*pi)-cos(phi_CEP2*2d0*pi))&
          &/(4*(omega2+2*Pi/pulse_tw2))&
          &+(cos((omega2-2*Pi/pulse_tw2)*(tt-T1_T2)+phi_CEP2*2d0*pi)-cos(phi_CEP2*2d0*pi))&
          /(4*(omega2-2*Pi/pulse_tw2))&
          &-0.5*(cos(omega2*pulse_tw2+phi_CEP2*2d0*pi)-cos(phi_CEP2*2d0*pi))/(omega2*((omega2*pulse_tw2/(2*pi))**2-1.d0))&
          &*(3.d0*((tt-T1_T2)/pulse_tw2)**2-2.d0*((tt-T1_T2)/pulse_tw2)**3)&
          &)
      else
        Ac_ext(iter,:)=Ac_ext(iter-1,:)
      endif
    enddo
    
  case('Asin2cos')
      ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! probe laser
    do iter=0,Nt+1
      tt=iter*dt
      if ( (tt-T1_T2>0d0) .and. (tt-T1_T2<pulse_tw2) ) then
        Ac_ext(iter,:)=Ac_ext(iter,:) &
          &-Epdir_re2(:)*f0_2/omega2*(sin(pi*(tt-T1_T2)/pulse_tw2))**2*cos(omega2*(tt-T1_T2)+phi_CEP2*2d0*pi)
      endif
    enddo
  case('input')
    !There is no probe
  case('Asin2_cw')
    !There is no probe
  case('none')
  case default
    call Err_finalize("Invalid pulse_shape_2 parameter!")
  end select
  Ac_ind=0.d0

  Ac_tot=Ac_ind+Ac_ext

  !hidden option (AY)
  if(alocal_laser=='y') call init_Ac_alocal

  return
End Subroutine init_Ac
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
