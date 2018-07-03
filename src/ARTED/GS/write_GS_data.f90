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
Subroutine write_GS_data
  use Global_Variables
  use salmon_global, only: out_dos, &
                         & out_dos_start, &
                         & out_dos_end, &
                         & iout_dos_nenergy, &
                         & out_dos_smearing, &
                         & out_dos_method, &
                         & out_dos_fshift
  use inputoutput, only: unit_length, t_unit_energy_inv, t_unit_energy
  use salmon_parallel, only: nproc_group_global, nproc_id_global, nproc_group_tdks
  use salmon_communication, only: comm_is_root,comm_summation, comm_bcast, comm_sync_all
  use salmon_file, only: open_filehandle
  implicit none
  integer ik,ib,ia,iter,j

  if (comm_is_root(nproc_id_global)) then
    open(403,file=file_GS)
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#grid information-----------------------------------------'
    write(403,*) '#aL =',aL
    write(403,*) '#al(1),al(2),al(3) =',al(1),al(2),al(3)
    write(403,*) '#aLx,aLy,aLz =',aLx,aLy,aLz
    write(403,*) '#bLx,bLy,bLz =',bLx,bLy,bLz
    write(403,*) '#Nd =',Nd
    write(403,*) '#NLx,NLy,NLz=',NLx,NLy,NLz
    write(403,*) '#NL =',NL
    write(403,*) '#Hx,Hy,Hz =',Hx,Hy,Hz
    write(403,*) '#(pi/max(Hx,Hy,Hz))**2 =',(pi/max(Hx,Hy,Hz))**2
    write(403,*) '#(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2 =',(pi/Hx)**2+(pi/Hy)**2+(pi/Hz)**2
    write(403,*) '#Hxyz =',Hxyz
    write(403,*) '#NKx,NKy,NKz=',NKx,NKy,NKz
    write(403,*) '#NKxyz =',NKxyz
    write(403,*) '#Sym=',Sym
    write(403,*) '#NK =',NK
    write(403,*) '#NEwald, aEwald =',NEwald, aEwald 
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS calc. option------------------------------------------'
    write(403,*) '#FSset_option =',FSset_option
    write(403,*) '#Ncg=',Ncg
    write(403,*) '#Nmemory_MB,alpha_MB =',Nmemory_MB,alpha_MB
    write(403,*) '#NFSset_start,NFSset_every =',NFSset_start,NFSset_every
    write(403,*) '#Nscf=',Nscf
    write(403,*) '#Nscf_conv=',Nscf_conv
    write(403,*) '#NI,NE=',NI,NE
    write(403,*) '#Zatom=',(Zatom(j),j=1,NE)
    write(403,*) '#Lref=',(Lref(j),j=1,NE)
    write(403,*) '#i,Kion(ia)','(Rion(j,a),j=1,3)'
    do ia=1,NI
      write(403,*) '#',ia,Kion(ia)
      write(403,*) '#',(Rion(j,ia),j=1,3)
    end do
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#GS information-------------------------------------------'
    write(403,*) '#NB,Nelec=',NB,Nelec
    write(403,*) '#Eall =',Eall
    write(403,*) '#ddns(iter = Nscf_conv)',ddns(Nscf_conv)
    write(403,*) '#ddns_abs_1e(iter = Nscf_conv)',ddns_abs_1e(Nscf_conv)
    write(403,*) '#esp_var_ave(iter = Nscf_conv)',esp_var_ave(Nscf_conv)
    write(403,*) '#esp_var_max(iter = Nscf_conv)',esp_var_max(Nscf_conv)
    write(403,*) '#NBoccmax is ',NBoccmax
    write(403,*) '#---------------------------------------------------------'
    write(403,*) '#band information-----------------------------------------'
    write(403,*) '#Bottom of VB',minval(esp_vb_min(:))
    write(403,*) '#Top of VB',maxval(esp_vb_max(:))
    write(403,*) '#Bottom of CB',minval(esp_cb_min(:))
    write(403,*) '#Top of CB',maxval(esp_cb_max(:))
    write(403,*) '#Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(403,*) '#Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*2.0*Ry
    write(403,*) '#BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(403,*) '#BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*2.0*Ry
    write(403,*) '#Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(403,*) '#Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(403,*) '#---------------------------------------------------------'
    write(403,'(a)') ' #iter     total-energy          ddns/nelec         esp_var_ave         esp_var_max'
    do iter=1,Nscf_conv
      !write(403,'(1x,i5,4e20.10)') iter, Eall_GS(iter),ddns(iter),esp_var_ave(iter),esp_var_max(iter)
      write(403,'(1x,i5,4e20.10)') iter, Eall_GS(iter),ddns_abs_1e(iter),esp_var_ave(iter),esp_var_max(iter)
    end do
    close(403)

    !! NOTE: We have changed the SYSNAME_band.out -> SYSNAME_eigen.data
    ! open(405,file=file_band)
    ! write(405,*) '#Bandmap at Ground State'
    ! write(405,*) '#(NK,NB)=','(',NK,NB,')'
    ! do ik=1,NK
    !   do ib=1,NB
    !     write(405,'(1x,2i5,5e26.16e3)') ik,ib,kAc(ik,1),kAc(ik,2),kAc(ik,3) &
    !       ,esp(ib,ik),occ(ib,ik)
    !   enddo
    ! enddo
    ! close(405)

  end if

  
  if(out_dos == 'y') call write_dos_data
  if(out_psi == 'y') call write_psi_data
  call write_k_data
  call write_eigen_data
  if(out_tm  == 'y') call write_tm_data

  return

  contains

    subroutine write_dos_data
      implicit none
      integer :: fh_dos
      real(8) :: vbmax, cbmin, emax, emin, efermi, eshift
      real(8) :: ww, fk, dw
      integer :: iw
      real(8) :: dos(iout_dos_nenergy), dos_l(iout_dos_nenergy)  
    
      if (comm_is_root(nproc_id_global)) then

        emin = minval(esp(:,:))
        emax = maxval(esp(:,:))
        cbmin = minval(esp_cb_min(:))
        vbmax = maxval(esp_vb_max(:))
      
        if (out_dos_fshift == 'y') then
          efermi = (vbmax + cbmin) * 0.5d0
          eshift = efermi
        else
          eshift = 0d0
        endif
      end if
      call comm_bcast(emin,nproc_group_global)
      call comm_bcast(emax,nproc_group_global)
      call comm_bcast(eshift,nproc_group_global)
      
      out_dos_start = max(out_dos_start, emin - 0.25d0 * (emax - emin))
      out_dos_end = min(out_dos_end, emax + 0.25d0 * (emax - emin))
      
      dos_l = 0d0
      dw = (out_dos_end - out_dos_start) / (iout_dos_nenergy - 1)
      
      select case (out_dos_method)
      case('lorentzian')
        !$omp parallel do private(ik,ib,fk,iw,ww) reduction(+:dos_l) collapse(2)
        do ik = NK_s,NK_e
          do ib = 1,NB
            fk = 2.d0/(NKxyz)*wk(ik)*(out_dos_smearing/pi)          
            do iw = 1, iout_dos_nenergy
              ww =  out_dos_start + (iw-1) * dw + eshift - esp(ib,ik) 
              dos_l(iw) = dos_l(iw) + fk/(ww**2 + out_dos_smearing**2)
            end do
          end do
        end do
        
      case('gaussian')
        !$omp parallel do private(ik,ib,fk,iw,ww) reduction(+:dos_l) collapse(2)
        do ik = NK_s,NK_e
          do ib = 1,NB
            fk = (2.d0 / (NKxyz * sqrt(2.d0*pi) * out_dos_smearing)) * wk(ik)
            do iw = 1, iout_dos_nenergy
              ww =  out_dos_start + (iw-1) * dw + eshift - esp(ib,ik) 
              dos_l(iw) = dos_l(iw) + fk * exp(-(0.5/out_dos_smearing**2)*ww**2)
            end do
          end do
        end do
      end select
      
      call comm_summation(dos_l,dos,iout_dos_nenergy,nproc_group_tdks)

      if (comm_is_root(nproc_id_global)) then
        fh_dos = open_filehandle(file_dos)
        write(fh_dos, '("#",1X,A)') "Density of States"
        
        write(fh_dos, '("#",1X,A,":",1X,A)') "Energy", "Electron energy"
        write(fh_dos, '("#",1X,A,":",1X,A)') "DoS", "Density of States"
        
        write(fh_dos, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "Energy", trim(t_unit_energy%name), &
          & 2, "DoS", trim(t_unit_energy_inv%name)
        do iw = 1, iout_dos_nenergy
          ww =  out_dos_start + (iw-1) * dw 
          write(fh_dos,'(F16.8,99(1X,E23.15E3))') &
            & ww * t_unit_energy%conv, &
            & dos(iw) * t_unit_energy_inv%conv
        end do
        close(fh_dos)
      end if
      call comm_sync_all
      
    end subroutine write_dos_data

  !--------------------------------------------------------------------------------
  !! export all orbital wave functions in cube or vtk format (multiplying phase factor)
    subroutine write_psi_data()
      use misc_routines
      implicit none
      integer :: fh_psi
      integer :: ik,ib,i,j,ix,iy,iz
      real(8) :: kr,psi
      character(256) :: gs_wfn_k_cube_vtk_dir

      fh_psi=502

      select case(format3d)
      case ('cube')
         write(gs_wfn_k_cube_vtk_dir,'(A,A)') trim(directory),'/gs_wfn_cube/'
         call create_directory(gs_wfn_k_cube_vtk_dir)

         do ik=NK_s,NK_e
         do ib=1,NB

            write(file_psi_gs,7000) trim(gs_wfn_k_cube_vtk_dir),ib,ik
            open(fh_psi,file=file_psi_gs,status="unknown")

            write(fh_psi,8010)
            write(fh_psi,8020) NI,  0d0, 0d0, 0d0
            write(fh_psi,8020) NLx, Hx,  0d0, 0d0
            write(fh_psi,8020) NLy, 0d0, Hy,  0d0
            write(fh_psi,8020) NLz, 0d0, 0d0, Hz

            do i=1, NI
               write(fh_psi,8030) Zatom(Kion(i)), 0d0, Rion(1,i),Rion(2,i),Rion(3,i) 
            end do

            j=1
            do ix=0, NLx-1
            do iy=0, NLy-1
            do iz=0, NLz-1
               i=Lxyz(ix,iy,iz)
               kr = kac0(ik,1)*Lx(i)*Hx + kac0(ik,2)*Ly(i)*Hy + kac0(ik,3)*Lz(i)*Hz
               psi= zu_GS0(i,ib,ik)*exp(zI*kr)
               if(mod(j,6)==0) then
                  write(fh_psi,8040) psi
               else
                  write(fh_psi,8040,advance='no') psi
               endif
               j=j+1
            end do
            end do
            end do

            close(fh_psi)
         enddo
         enddo

         !format for cube
7000     format(A,'/psi_gs_b',I5.5,'_k',I5.5,'.cube')
8010     format("# SALMON",/, &
         &      "# COMMENT" )
8020     format(I5,3(F12.6))
8030     format(I5,4(F12.6))
8040     format(ES12.4)


      case ('vtk')

         write(gs_wfn_k_cube_vtk_dir,'(A,A)') trim(directory),'/gs_wfn_vtk/'
         call create_directory(gs_wfn_k_cube_vtk_dir)

         do ik=NK_s,NK_e
         do ib=1,NB

            write(file_psi_gs,5000) trim(gs_wfn_k_cube_vtk_dir),ib,ik
            open(fh_psi,file=file_psi_gs,status="unknown")

            write(fh_psi,6010) 
            write(fh_psi,6020) NLx, NLy, NLz
            write(fh_psi,6030) 0d0, 0d0, 0d0
            write(fh_psi,6040) Hx, Hy, Hz
            write(fh_psi,6050) NLx * NLy * NLz
            write(fh_psi,6060) 

            do ix=0, NLx-1
            do iy=0, NLy-1
            do iz=0, NLz-1
               i=Lxyz(ix,iy,iz)
               kr = kac0(ik,1)*Lx(i)*Hx + kac0(ik,2)*Ly(i)*Hy + kac0(ik,3)*Lz(i)*Hz
               psi= zu_GS0(i,ib,ik)*exp(zI*kr)
               write(fh_psi,6070) psi
            end do
            end do
            end do

            close(fh_psi)
         enddo
         enddo

         !format for vtk
5000     format(A,'/psi_gs_b',I5.5,'_k',I5.5,'.vtk')
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


    end subroutine write_psi_data


  !--------------------------------------------------------------------------------
  !! export SYSNAME_k.data file
    subroutine write_k_data()
      implicit none
      integer :: fh_k
      integer :: ik

      if (comm_is_root(nproc_id_global)) then
        fh_k = open_filehandle(file_k_data, status="replace")
        write(fh_k, '("#",1X,A)') "k-point distribution"
        
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "Reduced coordinate of k-points"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "none", &
          & 3, "ky", "none", &
          & 4, "kz", "none", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, & 
            & kAc0(ik,1) / bLx, &
            & kAc0(ik,2) / bLy, &
            & kAc0(ik,3) / bLz, &
            & wk(ik)
        end do !ik
        close(fh_k)
      end if
      call comm_sync_all
    end subroutine write_k_data

  !--------------------------------------------------------------------------------
  !! export SYSNAME_eigen.data file
    subroutine write_eigen_data()
      implicit none
      integer :: fh_eigen
      integer :: ik, ib
   
      if (comm_is_root(nproc_id_global)) then
        fh_eigen = open_filehandle(file_eigen_data, status="replace")
        write(fh_eigen, '("#",1X,A)') "Ground state eigenenergies"
        
        write(fh_eigen, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "ib", "Band index"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "energy", "Eigenenergy"
        write(fh_eigen, '("#",1X,A,":",1X,A)') "occup", "Occupation"
        
        write(fh_eigen, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "ib", "none", &
          & 3, "energy", trim(t_unit_energy%name), &
          & 4, "occup", "none"
        do ik = 1, NK
          do ib = 1, NB
            write(fh_eigen, '(I6,1X,I6,99(1X,E23.15E3))') &
              & ik, ib, esp(ib,ik)*t_unit_energy%conv, occ(ib,ik)/wk(ik)*NKxyz
          end do !ib
        end do !ik
        close(fh_eigen)
      end if
      call comm_sync_all
    end subroutine write_eigen_data

  !! export SYSNAME_eigen.data file
    subroutine write_tm_data()
      use projector
      implicit none
      integer :: fh_tm
      integer :: i,j,ik,ib,ib1,ib2, ilma,ia,ix,iy,iz
      real(8) :: nabt(12), x,y,z
      complex(8) :: u_nab_u(3), upu(3,NB,NB,NK),upu_l(3,NB,NB,NK)
      complex(8) :: uVpsi(NB,NK),uVpsix(NB,NK),uVpsiy(NB,NK),uVpsiz(NB,NK)
      complex(8) :: uVpsixx(NB,NK),uVpsixy(NB,NK),uVpsixz(NB,NK)
      complex(8) :: uVpsiyy(NB,NK),uVpsiyz(NB,NK),uVpsizz(NB,NK)
      complex(8) :: u_rVnl_Vnlr_u(3,NB,NB,NK),u_rVnl_Vnlr_u_l(3,NB,NB,NK)
      complex(8) :: u_rVnl_u(3),u_Vnlr_u(3),veik
      complex(8) :: u_rVnlr_Vnlrr_u(3,3,NB,NK),u_rVnlr_Vnlrr_u_l(3,3,NB,NK)
      complex(8) :: ctmp1,ctmp2

      !calculate <u_nk|p_j|u_mk>  (j=x,y,z)
      nabt( 1: 4)=nabx(1:4)
      nabt( 5: 8)=naby(1:4)
      nabt( 9:12)=nabz(1:4)

      upu_l(:,:,:,:) = 0d0
      do ik=NK_s,NK_e
!$omp parallel
!$omp do private(ib1,ib2,u_nab_u) collapse(2)
      do ib1=1,NB
      do ib2=1,NB
         call u_nab_u_stencil(zu_GS0(:,ib1,ik),nabt,zu_GS0(:,ib2,ik),u_nab_u)
         upu_l(:,ib1,ib2,ik) = -zI*u_nab_u(:)*Hxyz
      enddo
      enddo
!$omp end do
!$omp end parallel
      enddo
      call comm_summation(upu_l,upu,3*NB*NB*NK,nproc_group_tdks)


      call update_projector(kac)

      !calculate <u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)
      u_rVnl_Vnlr_u_l = 0d0
      do ik=NK_s,NK_e
      do ilma=1,Nlma
         ia=a_tbl(ilma)
         uVpsi=0d0;  uVpsix=0d0;  uVpsiy=0d0;  uVpsiz=0d0

         do j=1,Mps(ia)
            i=Jxyz(j,ia)

            ix=Jxx(j,ia);  x=Lx(i)*Hx-ix*aLx
            iy=Jyy(j,ia);  y=Ly(i)*Hy-iy*aLy
            iz=Jzz(j,ia);  z=Lz(i)*Hz-iz*aLz

            veik = conjg(zproj(j,ilma,ik))
            do ib=1,NB
            uVpsi( ib,ik) =uVpsi( ib,ik)+ veik*    zu_GS0(i,ib,ik) !=<v|e^ik|u>
            uVpsix(ib,ik) =uVpsix(ib,ik)+ veik* x *zu_GS0(i,ib,ik) !=<v|e^ik*x|u>
            uVpsiy(ib,ik) =uVpsiy(ib,ik)+ veik* y *zu_GS0(i,ib,ik) !=<v|e^ik*y|u>
            uVpsiz(ib,ik) =uVpsiz(ib,ik)+ veik* z *zu_GS0(i,ib,ik) !=<v|e^ik*z|u>
            enddo
         end do

         uVpsi  = uVpsi *Hxyz *iuV(ilma)
         uVpsix = uVpsix*Hxyz
         uVpsiy = uVpsiy*Hxyz
         uVpsiz = uVpsiz*Hxyz

!$omp parallel
!$omp do private(ib1,ib2,u_rVnl_u,u_Vnlr_u) collapse(2)
         do ib1=1,NB
         do ib2=1,NB
            !<u|e^-ik*r|v><v|e^ik|u>
            u_rVnl_u(1)= conjg(uVpsix(ib1,ik))*uVpsi(ib2,ik) 
            u_rVnl_u(2)= conjg(uVpsiy(ib1,ik))*uVpsi(ib2,ik) 
            u_rVnl_u(3)= conjg(uVpsiz(ib1,ik))*uVpsi(ib2,ik) 
            !<u|e^-ik|v><v|e^ik*r|u>
            u_Vnlr_u(1)= conjg(uVpsi(ib1,ik))*uVpsix(ib2,ik) 
            u_Vnlr_u(2)= conjg(uVpsi(ib1,ik))*uVpsiy(ib2,ik) 
            u_Vnlr_u(3)= conjg(uVpsi(ib1,ik))*uVpsiz(ib2,ik) 

            u_rVnl_Vnlr_u_l(:,ib1,ib2,ik) = u_rVnl_Vnlr_u_l(:,ib1,ib2,ik)  &
            &                               + u_rVnl_u(:) - u_Vnlr_u(:)
         enddo
         enddo
!$omp end do
!$omp end parallel
      enddo  !ilma
      enddo  !ik
      call comm_summation(u_rVnl_Vnlr_u_l,u_rVnl_Vnlr_u,3*NB*NB*NK,nproc_group_tdks)


      !calculate <u_nk|[r_j,dVnl^(0)]r|u_nk>  (j=x,y,z)
      u_rVnlr_Vnlrr_u_l(:,:,:,:) = 0d0
      do ik=NK_s,NK_e
      do ilma=1,Nlma
         ia=a_tbl(ilma)
         uVpsi=0d0;  uVpsix=0d0;  uVpsiy=0d0;  uVpsiz=0d0
         uVpsixx=0d0;  uVpsixy=0d0;  uVpsixz=0d0
                       uVpsiyy=0d0;  uVpsiyz=0d0
                                     uVpsizz=0d0
         do j=1,Mps(ia)
            i=Jxyz(j,ia)

            ix=Jxx(j,ia);  x=Lx(i)*Hx-ix*aLx
            iy=Jyy(j,ia);  y=Ly(i)*Hy-iy*aLy
            iz=Jzz(j,ia);  z=Lz(i)*Hz-iz*aLz

            veik = conjg(zproj(j,ilma,ik))
            do ib=1,NB
            uVpsi(  ib,ik)=uVpsi(  ib,ik)+veik*    zu_GS0(i,ib,ik) !=<v|e^ik|u>
            uVpsix( ib,ik)=uVpsix( ib,ik)+veik* x *zu_GS0(i,ib,ik) !=<v|e^ik*x|u>
            uVpsiy( ib,ik)=uVpsiy( ib,ik)+veik* y *zu_GS0(i,ib,ik) !=<v|e^ik*y|u>
            uVpsiz( ib,ik)=uVpsiz( ib,ik)+veik* z *zu_GS0(i,ib,ik) !=<v|e^ik*z|u>
            uVpsixx(ib,ik)=uVpsixx(ib,ik)+veik*x*x*zu_GS0(i,ib,ik) !=<v|e^ik*xx|u>
            uVpsixy(ib,ik)=uVpsixy(ib,ik)+veik*x*y*zu_GS0(i,ib,ik) !=<v|e^ik*xy|u>
            uVpsixz(ib,ik)=uVpsixz(ib,ik)+veik*x*z*zu_GS0(i,ib,ik) !=<v|e^ik*xz|u>
            uVpsiyy(ib,ik)=uVpsiyy(ib,ik)+veik*y*y*zu_GS0(i,ib,ik) !=<v|e^ik*yy|u>
            uVpsiyz(ib,ik)=uVpsiyz(ib,ik)+veik*y*z*zu_GS0(i,ib,ik) !=<v|e^ik*yz|u>
            uVpsizz(ib,ik)=uVpsizz(ib,ik)+veik*z*z*zu_GS0(i,ib,ik) !=<v|e^ik*zz|u>
            enddo

         end do
         uVpsi  = uVpsi *Hxyz
         uVpsix = uVpsix*Hxyz
         uVpsiy = uVpsiy*Hxyz
         uVpsiz = uVpsiz*Hxyz

         uVpsixx = uVpsixx*Hxyz
         uVpsixy = uVpsixy*Hxyz
         uVpsixz = uVpsixz*Hxyz
         uVpsiyy = uVpsiyy*Hxyz
         uVpsiyz = uVpsiyz*Hxyz
         uVpsizz = uVpsizz*Hxyz

         do ib=1,NB
            !xx
            ctmp1 = conjg(uVpsix(ib,ik))*uVpsix( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsixx(ib,ik)
            u_rVnlr_Vnlrr_u_l(1,1,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(1,1,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !xy
            ctmp1 = conjg(uVpsix(ib,ik))*uVpsiy( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
            u_rVnlr_Vnlrr_u_l(1,2,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(1,2,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !xz
            ctmp1 = conjg(uVpsix(ib,ik))*uVpsiz( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
            u_rVnlr_Vnlrr_u_l(1,3,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(1,3,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !yx
            ctmp1 = conjg(uVpsiy(ib,ik))*uVpsix( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
            u_rVnlr_Vnlrr_u_l(2,1,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(2,1,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !yy
            ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiy( ib,ik)  
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyy(ib,ik)
            u_rVnlr_Vnlrr_u_l(2,2,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(2,2,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !yz
            ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiz( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
            u_rVnlr_Vnlrr_u_l(2,3,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(2,3,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !zx
            ctmp1 = conjg(uVpsiz(ib,ik))*uVpsix( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
            u_rVnlr_Vnlrr_u_l(3,1,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(3,1,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !zy
            ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiy( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
            u_rVnlr_Vnlrr_u_l(3,2,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(3,2,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)
            !zz
            ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiz( ib,ik)
            ctmp2 = conjg(uVpsi( ib,ik))*uVpsizz(ib,ik)
            u_rVnlr_Vnlrr_u_l(3,3,ib,ik) = &
            u_rVnlr_Vnlrr_u_l(3,3,ib,ik) + (ctmp1 - ctmp2)*iuV(ilma)

         enddo

      enddo  !ilma
      enddo  !ik
      call comm_summation(u_rVnlr_Vnlrr_u_l,u_rVnlr_Vnlrr_u,3*3*NB*NK,nproc_group_tdks)

    
      if (comm_is_root(nproc_id_global)) then
        fh_tm = open_filehandle(file_tm_data, status="replace")
        write(fh_tm, '("#",1X,A)') "#Transition Moment between occupied and unocupied orbitals in GS"
        write(fh_tm, '("#",1X,A)') "# (Separated analysis tool is available)"


        !Currently, TEST: print format is not decided

         !<u_nk|p_j|u_mk>  (j=x,y,z)
         write(fh_tm,*) "#<u_nk|p_j|u_mk>  (j=x,y,z)"
         do ik =1,NK
         do ib1=1,NB
         do ib2=1,NB
            write(fh_tm,9000) ik,ib1,ib2,(upu(j,ib1,ib2,ik),j=1,3)
         enddo
         enddo
         enddo
!9000     format(3i8,6e18.10)
9000     format(3i8,6e18.5)

         !<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)
         write(fh_tm,*) "#<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)"
         do ik=1,NK
         do ib=1,NB
            do i=1,3
               write(fh_tm,9000) ik,ib,i,(u_rVnlr_Vnlrr_u(i,j,ib,ik),j=1,3)
            enddo
         enddo
         enddo

         !<u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)
         write(fh_tm,*) "#<u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)"
         do ik =1,NK
         do ib1=1,NB
         do ib2=1,NB
            write(fh_tm,9000) ik,ib1,ib2,(u_rVnl_Vnlr_u(j,ib1,ib2,ik),j=1,3)
         enddo
         enddo
         enddo

      end if
        

      if (comm_is_root(nproc_id_global)) then
        close(fh_tm)
      end if
      call comm_sync_all
    end subroutine write_tm_data

    subroutine u_nab_u_stencil(A,B,C,D)
      use global_variables, only: NLx,NLy,NLz,zI
#ifndef ARTED_DOMAIN_POWER_OF_TWO
      use opt_variables, only: modx, mody, modz
#endif
      implicit none
      complex(8),intent(in)   :: A(0:NLz-1,0:NLy-1,0:NLx-1)
      real(8),   intent(in)   :: B(12)
      complex(8),intent(in)   :: C(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8),intent(out)  :: D(3)

      integer    :: ix,iy,iz
      complex(8) :: w1,w2,w3

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1)
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx)
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix
#endif

      D(:)=0d0

      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1

         w1=(B( 9)*(C(IDZ(1))-C(IDZ(-1))) &
         &  +B(10)*(C(IDZ(2))-C(IDZ(-2))) &
         &  +B(11)*(C(IDZ(3))-C(IDZ(-3))) &
         &  +B(12)*(C(IDZ(4))-C(IDZ(-4))))

         w2=(B( 5)*(C(IDY(1))-C(IDY(-1))) &
         &  +B( 6)*(C(IDY(2))-C(IDY(-2))) &
         &  +B( 7)*(C(IDY(3))-C(IDY(-3))) &
         &  +B( 8)*(C(IDY(4))-C(IDY(-4))))

         w3=(B( 1)*(C(IDX(1))-C(IDX(-1))) &
         &  +B( 2)*(C(IDX(2))-C(IDX(-2))) &
         &  +B( 3)*(C(IDX(3))-C(IDX(-3))) &
         &  +B( 4)*(C(IDX(4))-C(IDX(-4))))

         D(1) = D(1) + conjg(A(iz,iy,ix))*w3  !x
         D(2) = D(2) + conjg(A(iz,iy,ix))*w2  !y
         D(3) = D(3) + conjg(A(iz,iy,ix))*w1  !z
      end do
      end do
      end do

    end subroutine u_nab_u_stencil
    
End Subroutine write_GS_data
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
