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
  use inputoutput, only: t_unit_current, t_unit_energy, t_unit_elec
  implicit none
  integer :: ihw,ixyz,iter
  real(8) :: hw,tt
  complex(8) :: jav_w(3),E_ext_w(3),E_tot_w(3)
  complex(8) :: jav_d(3),jav_s(3),smt_s
  complex(8) :: zsigma_w(3),zeps(3)
  integer :: fh_lr

  if (comm_is_root(nproc_id_global)) then
    fh_lr = open_filehandle(file_lr_data)

    write(fh_lr, '("#",1X,A)') "Fourier-transform spectra"

    if (ae_shape1 == 'impulse' .and. trans_Longi == 'lo') then
      write(fh_lr, '("#",1X,A,":",1X,A)') "eps", "Dielectric constant"
      write(fh_lr, '("#",1X,A,":",1X,A)') "sigma", "Conductivity"
      write(fh_lr, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Frequency", trim(t_unit_energy%name), &
        & 2, "Re(eps_x)", "none", &
        & 3, "Re(eps_y)", "none", &
        & 4, "Re(eps_z)", "none", &
        & 5, "Im(eps_x)", "none", &
        & 6, "Im(eps_y)", "none", &
        & 7, "Im(eps_z)", "none"
    else if (ae_shape1 == 'impulse' .and. Trans_Longi == 'tr') then
      write(fh_lr, '("#",1X,A,":",1X,A)') "sigma", "Conductivity"
      write(fh_lr, '("#",1X,A,":",1X,A)') "eps", "Dielectric constant"
      write(fh_lr, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Frequency", trim(t_unit_energy%name), &
        & 2, "Re(sigma_x)", "a.u.", &
        & 3, "Re(sigma_y)", "a.u.", &
        & 4, "Re(sigma_z)", "a.u.", &
        & 5, "Im(sigma_x)", "a.u.", &
        & 6, "Im(sigma_y)", "a.u.", &
        & 7, "Im(sigma_z)", "a.u.", &
        & 8, "Re(eps_x)", "none", &
        & 9, "Re(eps_y)", "none", &
        & 10, "Re(eps_z)", "none", &
        & 11, "Im(eps_x)", "none", &
        & 12, "Im(eps_y)", "none", &
        & 13, "Im(eps_z)", "none"
    else
      write(fh_lr, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
      write(fh_lr, '("#",1X,A,":",1X,A)') "E_ext", "External electric field"
      write(fh_lr, '("#",1X,A,":",1X,A)') "E_tot", "Total electric potential field"
      write(fh_lr, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Frequency", trim(t_unit_energy%name), &
        & 2, "Re(Jm_x)", trim(t_unit_current%name), &
        & 3, "Re(Jm_y)", trim(t_unit_current%name), &
        & 4, "Re(Jm_z)", trim(t_unit_current%name), &
        & 5, "Im(Jm_x)", trim(t_unit_current%name), &
        & 6, "Im(Jm_y)", trim(t_unit_current%name), &
        & 7, "Im(Jm_z)", trim(t_unit_current%name), &
        & 8, "Re(E_ext_x)", trim(t_unit_elec%name), &
        & 9, "Re(E_ext_y)", trim(t_unit_elec%name), &
        & 10, "Re(E_ext_z)", trim(t_unit_elec%name), &
        & 11, "Im(E_ext_x)", trim(t_unit_elec%name), &
        & 12, "Im(E_ext_y)", trim(t_unit_elec%name), &
        & 13, "Im(E_ext_z)", trim(t_unit_elec%name), &
        & 14, "Re(E_tot_x)", trim(t_unit_elec%name), &
        & 15, "Re(E_tot_y)", trim(t_unit_elec%name), &
        & 16, "Re(E_tot_z)", trim(t_unit_elec%name), &
        & 17, "Im(E_tot_x)", trim(t_unit_elec%name), &
        & 18, "Im(E_tot_y)", trim(t_unit_elec%name), &
        & 19, "Im(E_tot_z)", trim(t_unit_elec%name)
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
    if (comm_is_root(nproc_id_global)) then
      if (ae_shape1 == 'impulse' .and. trans_Longi == 'lo') then
        write(fh_lr,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)
      else if (ae_shape1 == 'impulse' .and. Trans_Longi == 'tr') then
        write(fh_lr,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zsigma_w(ixyz)),ixyz=1,3)&
             &,(aimag(zsigma_w(ixyz)),ixyz=1,3)&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)
      else
        write(fh_lr,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(jav_w(ixyz)) * t_unit_current%conv, ixyz=1,3)&
             &,(aimag(jav_w(ixyz)) * t_unit_current%conv, ixyz=1,3)&
             &,(real(E_ext_w(ixyz)) * t_unit_elec%conv ,ixyz=1,3)&
             &,(aimag(E_ext_w(ixyz)) * t_unit_elec%conv ,ixyz=1,3)&
             &,(real(E_tot_w(ixyz)) * t_unit_elec%conv ,ixyz=1,3)&
             &,(aimag(E_tot_w(ixyz)) * t_unit_elec%conv ,ixyz=1,3)
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

subroutine analysis_dns_trans(it)
  use Global_Variables
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use inputoutput, only: t_unit_energy
  implicit none
  integer :: it, fh_dns_trans1,fh_dns_trans2,i,j,ix,iy,iz
  real(8) :: omg,tt
  character(256) :: file_dns_trans1,file_dns_trans2

  if(.not. allocated(rho_trans)) then
     allocate(rho_trans(NL))
     rho_trans(:)=0d0
  endif

  omg = out_dns_trans_energy
  tt  = it*dt
  rho_trans(:) = rho_trans(:) + (rho(:)-rho_gs(:))*exp(zI*omg*tt)

  if(it==Nt)then  
      rho_trans(:) = rho_trans(:)/Nt

      ! Print
      fh_dns_trans1 = 503
      fh_dns_trans2 = 504
      if (comm_is_root(nproc_id_global)) then

      select case(format3d)
      case ('cube')
         file_dns_trans1 = "out_dns_trans_Re.cube"
         file_dns_trans2 = "out_dns_trans_Im.cube"
         open(fh_dns_trans1,file=file_dns_trans1,status="unknown")
         open(fh_dns_trans2,file=file_dns_trans2,status="unknown")

         write(fh_dns_trans1,8010) omg*t_unit_energy%conv,trim(t_unit_energy%name)
         write(fh_dns_trans1,8020) NI,  0d0, 0d0, 0d0
         write(fh_dns_trans1,8020) NLx, Hx,  0d0, 0d0
         write(fh_dns_trans1,8020) NLy, 0d0, Hy,  0d0
         write(fh_dns_trans1,8020) NLz, 0d0, 0d0, Hz

         write(fh_dns_trans2,8010) omg*t_unit_energy%conv,trim(t_unit_energy%name)
         write(fh_dns_trans2,8020) NI,  0d0, 0d0, 0d0
         write(fh_dns_trans2,8020) NLx, Hx,  0d0, 0d0
         write(fh_dns_trans2,8020) NLy, 0d0, Hy,  0d0
         write(fh_dns_trans2,8020) NLz, 0d0, 0d0, Hz

         do i=1, NI
            write(fh_dns_trans1,8030) Zatom(Kion(i)),0d0,Rion(1,i),Rion(2,i),Rion(3,i) 
            write(fh_dns_trans2,8030) Zatom(Kion(i)),0d0,Rion(1,i),Rion(2,i),Rion(3,i) 
         end do

         j=1
         do ix=0, NLx-1
         do iy=0, NLy-1
         do iz=0, NLz-1
            i=Lxyz(ix,iy,iz)
            if(mod(j,6)==0) then
               write(fh_dns_trans1,8040) real(rho_trans(i))
               write(fh_dns_trans2,8040) aimag(rho_trans(i))
            else
               write(fh_dns_trans1,8040,advance='no') real(rho_trans(i))
               write(fh_dns_trans2,8040,advance='no') aimag(rho_trans(i))
            endif
            j=j+1
         end do
         end do
         end do

         close(fh_dns_trans1)
         close(fh_dns_trans2)

         !format for cube
8010     format("# SALMON",/, &
         &      "# out_dns_trans option: energy in FT =",f18.10," [",a,"]" )
8020     format(I5,3(F12.6))
8030     format(I5,4(F12.6))
8040     format(ES12.4)

      case ('vtk')

         file_dns_trans1 = "out_dns_trans.cube"
         file_dns_trans2 = "out_dns_trans_ratio.cube"
         open(fh_dns_trans1,file=file_dns_trans1,status="unknown")
         open(fh_dns_trans2,file=file_dns_trans2,status="unknown")

         write(fh_dns_trans1,6010) 
         write(fh_dns_trans1,6020) NLx, NLy, NLz
         write(fh_dns_trans1,6030) 0d0, 0d0, 0d0
         write(fh_dns_trans1,6040) Hx, Hy, Hz
         write(fh_dns_trans1,6050) NLx * NLy * NLz
         write(fh_dns_trans1,6060) 

         write(fh_dns_trans2,6010) 
         write(fh_dns_trans2,6020) NLx, NLy, NLz
         write(fh_dns_trans2,6030) 0d0, 0d0, 0d0
         write(fh_dns_trans2,6040) Hx, Hy, Hz
         write(fh_dns_trans2,6050) NLx * NLy * NLz
         write(fh_dns_trans2,6060) 

         do ix=0, NLx-1
         do iy=0, NLy-1
         do iz=0, NLz-1
            i=Lxyz(ix,iy,iz)
            write(fh_dns_trans1,6070) real(rho_trans(i))
            write(fh_dns_trans2,6070) aimag(rho_trans(i))
         end do
         end do
         end do

         close(fh_dns_trans1)
         close(fh_dns_trans2)


         !format for vtk
6010     format("# vtk DataFile Version 3.0",/, &
         &      "vtk output",/,                 &
         &      "ASCII",/,                      &
         &      "DATASET STRUCTURED_POINTS" )

6020     format("DIMENSIONS",3(1X,I2)  )
6030     format("ORIGIN",    3(1X,F3.1))
6040     format("SPACING",   3(1X,F7.3))
6050     format("POINT_DATA",1X,I6     )
6060     format("SCALARS scalars float",/, &
         &      "LOOKUP_TABLE default" )
6070     format(ES12.5)

      end select
      endif

  endif

  return
end Subroutine analysis_dns_trans

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
