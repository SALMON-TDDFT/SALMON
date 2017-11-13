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
Subroutine Fourier_tr
  use Global_Variables
  use salmon_file, only: open_filehandle
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  implicit none
  integer :: ihw,ixyz,iter
  real(8) :: hw,tt
  complex(8) :: jav_w(3),E_ext_w(3),E_tot_w(3)
  complex(8) :: jav_d(3),jav_s(3),smt_s
  complex(8) :: zsigma_w(3),zeps(3)
  integer :: fh_lr

  if (comm_is_root(nproc_id_global)) then
    fh_lr = open_filehandle(file_lr_data)
    if (ae_shape1 == 'impulse' .and. trans_Longi == 'lo') then
      write(fh_lr,"('#',99(1X,A))") &
        & "Frequency[au]", &
        & "Re(eps_x)", "Re(eps_y)", "Re(eps_z)", & 
        & "Im(eps_x)", "Im(eps_y)", "Im(eps_z)"
    else if (ae_shape1 == 'impulse' .and. Trans_Longi == 'tr') then
      write(fh_lr,"('#',99(1X,A))") &
        & "Frequency[au]", &
        & "Re(sigma_x)[au]", "Re(sigma_y)[au]", "Re(sigma_z)[au]", & 
        & "Im(sigma_x)[au]", "Im(sigma_y)[au]", "Im(sigma_z)[au]", &
        & "Re(eps_x)", "Re(eps_y)", "Re(eps_z)", & 
        & "Im(eps_x)", "Im(eps_y)", "Im(eps_z)"
    else
      write(fh_lr,"('#',99(1X,A))") &
        & "Frequency[au]", &
        & "Re(jmatter_x)[au]", "Re(jmatter_y)[au]", "Re(jmatter_z)[au]", &
        & "Im(jmatter_x)[au]", "Im(jmatter_y)[au]", "Im(jmatter_z)[au]", &
        & "Re(E_ext_x)[au]", "Re(E_ext_y)[au]", "Re(E_ext_z)[au]", &
        & "Im(E_ext_x)[au]", "Im(E_ext_y)[au]", "Im(E_ext_z)[au]", &
        & "Re(E_tot_x)[au]", "Re(E_tot_y)[au]", "Re(E_tot_z)[au]", &
        & "Im(E_tot_x)[au]", "Im(E_tot_y)[au]", "Im(E_tot_z)[au]"
    endif
  endif
  
  if (KbTev < 0d0) then
    ! sigma(omega=0) correcton
    jav_s=0d0; smt_s=0d0;
    do iter=0,Nt
      tt=iter*dt
      jav_s(:)=jav_s(:)+javt(iter,:)*smoothing_t(tt)
      smt_s=smt_s+smoothing_t(tt)
    end do
    jav_d(:)=jav_s(:)/smt_s
  else
    jav_d(:)=0d0
  end if

  do ihw=1,Nomega
    hw=ihw*domega
    jav_w=0.d0; E_ext_w=0.d0; E_tot_w=0.d0
    do iter=0,Nt
      tt=iter*dt
      jav_w(:)=jav_w(:)+(javt(iter,:)-jav_d(:))*exp(zI*hw*tt)*smoothing_t(tt)
      E_ext_w(:)=E_ext_w(:)+E_ext(iter,:)*exp(zI*hw*tt)*smoothing_t(tt)
      E_tot_w(:)=E_tot_w(:)+E_tot(iter,:)*exp(zI*hw*tt)*smoothing_t(tt)
    enddo
    jav_w(:)=jav_w(:)*dt; E_ext_w(:)=E_ext_w(:)*dt; E_tot_w(:)=E_tot_w(:)*dt
    if (ae_shape1 == 'impulse') then
      if (Trans_Longi == 'lo')  then 
        zeps(:)=1.d0/(1.d0-E_tot_w(:)/dAc)
      else if (Trans_Longi == 'tr') then
        zsigma_w(:)=jav_w(:)/dAc
        zeps=epdir_re1(:)+zI*4.d0*pi*zsigma_w(:)/hw
      end if
    end if
      write(fh, "(ES15.6E3,6(1X,ES22.14E3,1X))") &
    if (comm_is_root(nproc_id_global)) then
      if (ae_shape1 == 'impulse' .and. trans_Longi == 'lo') then
        write(fh_lr,'(ES15.6E3,99(1X,ES22.14E3))') hw &
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)
      else if (ae_shape1 == 'impulse' .and. Trans_Longi == 'tr') then
        write(fh_lr,'(ES15.6E3,99(1X,ES22.14E3))') hw &
             &,(real(zsigma_w(ixyz)),ixyz=1,3)&
             &,(aimag(zsigma_w(ixyz)),ixyz=1,3)&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)
      else
        write(fh_lr,'(ES15.6E3,99(1X,ES22.14E3))') hw &
             &,(real(jav_w(ixyz)),ixyz=1,3)&
             &,(aimag(jav_w(ixyz)),ixyz=1,3)&
             &,(real(E_ext_w(ixyz)),ixyz=1,3)&
             &,(aimag(E_ext_w(ixyz)),ixyz=1,3)&
             &,(real(E_tot_w(ixyz)),ixyz=1,3)&
             &,(aimag(E_tot_w(ixyz)),ixyz=1,3)
      endif
    endif
  enddo
 
  if (comm_is_root(nproc_id_global)) then
    close(fh_lr)
  endif

  return
Contains
!======
!======
Function smoothing_t(tt)
  implicit none
  real(8),intent(IN) :: tt
  real(8) :: smoothing_t
  
  smoothing_t=1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3

  return
End Function smoothing_t
!======
!======
End Subroutine Fourier_tr
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
