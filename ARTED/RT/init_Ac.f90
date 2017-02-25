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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine init_Ac
  use Global_Variables
  use communication
  implicit none
  integer :: iter
  real(8) :: tt

  select case(ext_field)
  case('LR')
    Ac_ext(:,1)=Epdir_1(1)*dAc
    Ac_ext(:,2)=Epdir_1(2)*dAc
    Ac_ext(:,3)=Epdir_1(3)*dAc
    Ac_ind=0.d0
    javt=0.d0
  case('LF')
    f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
    omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
    tpulse_1=tpulsefs_1/0.02418d0 ! pulse duration in a.u.
    f0_2=5.338d-9*sqrt(IWcm2_2)      ! electric field in a.u.
    omega_2=omegaev_2/(2d0*Ry)  ! frequency in a.u.
    tpulse_2=tpulsefs_2/0.02418d0 ! pulse duration in a.u.
    T1_T2=T1_T2fs/0.02418d0 ! pulse duration in a.u.
    javt=0.d0
    Ac_ext=0.d0

    select case(AE_shape)
    case('Esin2sin')
! pulse shape : E(t)=f0*sin(Pi t/T)**2 *sin (omega t+phi_CEP*2d0*pi) 
! pump laser
      do iter=0,Nt+1
        tt=iter*dt
        if (tt<tpulse_1) then
          Ac_ext(iter,:)=Epdir_1(:)*f0_1*(&
            &-(cos(omega_1*tt+phi_CEP_1*2d0*pi)-cos(phi_CEP_1*2d0*pi))/(2*omega_1)&
            &+(cos((omega_1+2*Pi/tpulse_1)*tt+phi_CEP_1*2d0*pi)-cos(phi_CEP_1*2d0*pi))/(4*(omega_1+2*Pi/tpulse_1))&
            &+(cos((omega_1-2*Pi/tpulse_1)*tt+phi_CEP_1*2d0*pi)-cos(phi_CEP_1*2d0*pi))/(4*(omega_1-2*Pi/tpulse_1))&
            &-0.5*(cos(omega_1*tpulse_1+phi_CEP_1*2d0*pi)-cos(phi_CEP_1*2d0*pi))/(omega_1*((omega_1*tpulse_1/(2*pi))**2-1.d0))&
            &*(3.d0*(tt/tpulse_1)**2-2.d0*(tt/tpulse_1)**3)&
            &)
        else
          Ac_ext(iter,:)=Ac_ext(iter-1,:)
        endif
      enddo

! probe laser
      do iter=0,Nt+1
        tt=iter*dt
        if(tt-T1_T2 <= 0d0)then
          
        else if ( (tt-T1_T2>0d0) .and. (tt-T1_T2<tpulse_2) ) then
          Ac_ext(iter,:)=Ac_ext(iter,:)+Epdir_2(:)*f0_2*(&
            &-(cos(omega_2*(tt-T1_T2)+phi_CEP_2*2d0*pi)-cos(phi_CEP_2*2d0*pi))/(2*omega_2)&
            &+(cos((omega_2+2*Pi/tpulse_2)*(tt-T1_T2)+phi_CEP_2*2d0*pi)-cos(phi_CEP_2*2d0*pi))/(4*(omega_2+2*Pi/tpulse_2))&
            &+(cos((omega_2-2*Pi/tpulse_2)*(tt-T1_T2)+phi_CEP_2*2d0*pi)-cos(phi_CEP_2*2d0*pi))/(4*(omega_2-2*Pi/tpulse_2))&
            &-0.5*(cos(omega_2*tpulse_2+phi_CEP_2*2d0*pi)-cos(phi_CEP_2*2d0*pi))/(omega_2*((omega_2*tpulse_2/(2*pi))**2-1.d0))&
            &*(3.d0*((tt-T1_T2)/tpulse_2)**2-2.d0*((tt-T1_T2)/tpulse_2)**3)&
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
        if (tt<tpulse_1) then
          Ac_ext(iter,:)=-Epdir_1(:)*f0_1/omega_1*(sin(pi*tt/tpulse_1))**2*cos(omega_1*tt+phi_CEP_1*2d0*pi)
        end if
      enddo
! probe laser
      do iter=0,Nt+1
        tt=iter*dt
        if ( (tt-T1_T2>0d0) .and. (tt-T1_T2<tpulse_2) ) then
          Ac_ext(iter,:)=Ac_ext(iter,:) &
            &-Epdir_2(:)*f0_2/omega_2*(sin(pi*(tt-T1_T2)/tpulse_2))**2*cos(omega_2*(tt-T1_T2)+phi_CEP_2*2d0*pi)
        endif
      enddo
    case('input')
      Ac_ext=0d0
      if(comm_is_root())then
        open(899,file='input_Ac.dat')
        do iter=0,Nt
          read(899,*)Ac_ext(iter,1),Ac_ext(iter,2),Ac_ext(iter,3)
        end do
        close(899)
      end if
      call comm_bcast(Ac_ext,proc_group(1))

    case('Asin2_cw')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
      do iter=0,Nt+1
        tt=iter*dt
        if (tt<tpulse_1*0.5d0) then
          Ac_ext(iter,:)=-Epdir_1(:)*f0_1/omega_1*(sin(pi*tt/tpulse_1))**2*cos(omega_1*tt+phi_CEP_1*2d0*pi)
        else
          Ac_ext(iter,:)=-Epdir_1(:)*f0_1/omega_1*cos(omega_1*tt+phi_CEP_1*2d0*pi)
        end if
      enddo

    end select

    Ac_ind=0.d0
  end select
  Ac_tot=Ac_ind+Ac_ext

  return
End Subroutine init_Ac
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
