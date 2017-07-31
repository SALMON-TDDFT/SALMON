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
!This file is "initialization.f90"
!This file contains initialization of solid state part
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module initialization
  implicit none
contains
  subroutine initialize
    use Global_Variables
    use timer
    use opt_variables
    use salmon_parallel
    use salmon_communication
    use environment
    use misc_routines
    use inputfile,only: transfer_input
    use restart,only: prep_restart_read
    implicit none
!$ integer :: omp_get_max_threads  

    call timer_initialize

    call load_environments

    if(comm_is_root(nproc_id_global)) then
       write(*,'(A)') "Welcome to SALMON-TDDFT"
       write(*,'(A)') "(Preliminary Developers Version)"
       write(*,'(2A)') "based on ARTED ver. = ",ARTED_ver
       call print_optimize_message
    end if

    NUMBER_THREADS=1
!$  NUMBER_THREADS=omp_get_max_threads()
!$  if(.true.) then
!$    if(comm_is_root(nproc_id_global))write(*,*)'parallel = Hybrid'
!$  else
    if(comm_is_root(nproc_id_global))write(*,*)'parallel = Flat MPI'
!$  end if

    if(comm_is_root(nproc_id_global))write(*,*)'NUMBER_THREADS = ',NUMBER_THREADS

    call timer_begin(LOG_ALL)

    call timer_begin(LOG_STATIC)
    Time_start=get_wtime() !reentrance
    call comm_bcast(Time_start,nproc_group_global)

    if(restart_option == 'restart') then
      if (comm_is_root(nproc_id_global)) call timer_show_current_hour('Restore...', LOG_ALL)
      call prep_restart_read
      return
    end if

    call transfer_input
    call Read_data

    call fd_coef
    call init

! initialize for optimization.
    call opt_vars_initialize_p1

    call input_pseudopotential_YS !shinohara

    call prep_ps_periodic('initial    ')

! initialize for optimization.
    call opt_vars_initialize_p2


  end subroutine initialize
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  Subroutine Read_data
    use salmon_global
    use Global_Variables
    use opt_variables
    use environment
    use salmon_parallel
    use salmon_communication
    use misc_routines
    use timer
    implicit none
    integer :: ia,i,j
    integer :: ix_m,iy_m
    
    if (comm_is_root(nproc_id_global)) then
       write(*,*) 'Nprocs=',nproc_size_global
       write(*,*) 'nproc_id_global=0:  ',nproc_id_global
       write(*,*) 'Time_shutdown=',Time_shutdown,'sec'
    end if
    
    if(comm_is_root(nproc_id_global))then
       
       need_backup = (backup_frequency > 0)
       write(*,*) 'need backup?',need_backup
       if (need_backup) write(*,*) '  frequency (# of iter) :',backup_frequency
       
       write(*,*) 'entrance_iter=',entrance_iter
       write(*,*) 'SYSname=',trim(SYSname)
       write(*,*) 'directory=',trim(directory)
       !yabana
       write(*,*) 'functional=',functional
       if(functional == 'TBmBJ') write(*,*) 'cvalue=',cval
       !yabana
       write(*,*) 'propagator=',propagator
       write(*,*) 'pseudo_file =',(trim(pseudo_file(i)),i=1,NE)
       write(*,*) 'PSmask_option =',PSmask_option !shinohara
       write(*,*) 'alpha_mask, gamma_mask, eta_mask =',alpha_mask, gamma_mask, eta_mask !shinohara
       file_GS=trim(directory)//trim(SYSname)//'_GS.out'
       file_RT=trim(directory)//trim(SYSname)//'_RT.out'
       file_epst=trim(directory)//trim(SYSname)//'_t.out'
       file_epse=trim(directory)//trim(SYSname)//'_e.out'
       file_force_dR=trim(directory)//trim(SYSname)//'_force_dR.out'
       file_j_ac=trim(directory)//trim(SYSname)//'_j_ac.out'
       file_DoS=trim(directory)//trim(SYSname)//'_DoS.out'
       file_band=trim(directory)//trim(SYSname)//'_band.out'
       file_dns=trim(directory)//trim(SYSname)//'_dns.out'
       file_ovlp=trim(directory)//trim(SYSname)//'_ovlp.out'
       file_nex=trim(directory)//trim(SYSname)//'_nex.out'
       write(*,*) 'al(1),al(2),al(3)=',al(1),al(2),al(3)
       write(*,*) 'Sym=',Sym,'crystal structure=',crystal_structure !sym
       write(*,*) 'Nd,NLx,NLy,NLz,NKx,NKy,NKz=',Nd,NLx,NLy,NLz,NKx,NKy,NKz
       write(*,*) 'FDTDdim=',FDTDdim
       write(*,*) 'TwoD_shape=',TwoD_shape 
       write(*,*) 'NX_m,NY_m=',NX_m,NY_m
       write(*,*) 'HX_m,HY_m=',HX_m,HY_m
       write(*,*) 'NKsplit,NXYsplit=',NKsplit,NXYsplit
       write(*,*) 'NXvacL_m,NXvacR_m=',NXvacL_m,NXvacR_m
       write(*,*) 'NEwald, aEwald =',NEwald, aEwald 
       write(*,*) 'KbTev=',KbTev ! sato
    end if
    
#ifdef ARTED_USE_FORTRAN2008
    write (process_directory,'(A,A,I5.5,A)') trim(directory),'/work_p',nproc_id_global,'/'
    call create_directory(process_directory)
#else
    process_directory = trim(directory)
#endif
    
    call comm_bcast(need_backup,nproc_group_global)
    call comm_bcast(file_GS,nproc_group_global)
    call comm_bcast(file_RT,nproc_group_global)
    call comm_bcast(file_epst,nproc_group_global)
    call comm_bcast(file_epse,nproc_group_global)
    call comm_bcast(file_force_dR,nproc_group_global)
    call comm_bcast(file_j_ac,nproc_group_global)
    call comm_bcast(file_DoS,nproc_group_global)
    call comm_bcast(file_band,nproc_group_global)
    call comm_bcast(file_dns,nproc_group_global)
    call comm_bcast(file_ovlp,nproc_group_global)
    call comm_bcast(file_nex,nproc_group_global)
    call comm_bcast(file_kw,nproc_group_global)

    if(use_ms_maxwell == 'y')then
       if(FDTDdim == '1D' .and. TwoD_shape /= 'periodic') then
          if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 1D calculation ! TwoD_shape is not good'
          TwoD_shape='periodic'
       end if
       if(FDTDdim == '1D' .and. NY_m /= 1) then
          if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 1D calculation ! NY_m is not good'
          NY_m=1
       end if
       if(FDTDdim == '2D' .and. TwoD_shape /= 'periodic') then
          if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 2D calculation ! TwoD_shape is not good'
          TwoD_shape='periodic'
       end if
    end if
    !sym ---
    select case(crystal_structure)
    case("diamond2")
       if(functional == "PZ" .or. functional == "PZM" &
            .or. functional == "TBmBJ" .or. functional == "BJ_PW")then
          if(Sym == 8)then
             if((mod(NLx,2)+mod(NLy,2)+mod(NLz,4)) /= 0)call err_finalize('Bad grid point')
             if(NLx /= NLy)call err_finalize('Bad grid point: NLx /= NLy')
             if(NKx /= NKy)call err_finalize('NKx /= NKy')
          else if(Sym /=1)then
             call err_finalize('Bad crystal structure')
          end if
       else
          if(Sym /= 1)call err_finalize('Bad crystal structure')
       end if
    case("diamond")
       if(functional == "PZ" .or. functional == "PZM" &
            .or. functional == "TBmBJ"  .or. functional == "BJ_PW")then
          if(Sym == 8)then
             if((mod(NLx,4)+mod(NLy,4)+mod(NLz,4)) /= 0)call err_finalize('Bad grid point')
             if(NLx /= NLy)call err_finalize('Bad grid point')
             if(NKx /= NKy) call err_finalize('NKx /= NKy')
          else if(Sym ==4 )then
             if(NLx /= NLy)call err_finalize('Bad grid point')
             if(NKx /= NKy) call err_finalize('NKx /= NKy')
          else if(Sym /= 1)then
             call err_finalize('Bad crystal structure')
          end if
       else
          if(Sym /= 1)call err_finalize('Bad crystal structure')
       end if
    case("tetragonal")
       if(functional == "PZ" .or. functional == "PZM")then
          if((mod(NLx,4)+mod(NLy,4)+mod(NLz,4)) /= 0)call err_finalize('Bad grid point')
          if(NLx /= NLy)call err_finalize('Bad grid point')
          if(NKx /= NKy) call err_finalize('NKx /= NKy')
       else if(Sym /= 1)then
          call err_finalize('Bad crystal structure')
       end if
    case default
       if(Sym /= 1)call err_finalize('Bad symmetry')
    end select
    !sym ---
    
    if(mod(NKx,2)+mod(NKy,2)+mod(NKz,2) /= 0) call err_finalize('NKx,NKy,NKz /= even')
    
    call comm_sync_all
    
    !  aLx=ax*aL;    aLy=ay*aL;    aLz=az*aL
    aLx=aL(1);    aLy=aL(2);    aLz=aL(3)
    aLxyz=aLx*aLy*aLz
    bLx=2*Pi/aLx; bLy=2*Pi/aLy; bLz=2*Pi/aLz
    Hx=aLx/NLx;   Hy=aLy/NLy;   Hz=aLz/NLz
    Hxyz=Hx*Hy*Hz
    NL=NLx*NLy*NLz
    NG=NL
    
    if (0<NKx .and. 0<NKy .and. 0<NKz) then
       NKxyz=NKx*NKy*NKz
       select case(Sym)
       case(1)
          NK=NKx*NKy*NKz
       case(4)
          NK=(NKx/2)*(NKy/2)*NKz
       case(8)
          NK=NKz*(NKx/2)*((NKx/2)+1)/2
       end select
    else
       if (comm_is_root(nproc_id_global)) then
          write(*,*) "Use non-uniform k-points distribution"
          write(*,*) "file_kw=", file_kw
          open(410, file=file_kw, status="old")
          read(410,*) NK, NKxyz
          close(410)
          write(*,*) "NK=", NK, "NKxyz=", NKxyz
       endif
       call comm_bcast(NK,nproc_group_global)
       call comm_bcast(NKxyz,nproc_group_global)
    endif
    

    if(use_ms_maxwell == 'y')then
       if(NXYsplit /= 1 .and. NKsplit /=1) call err_finalize('cannot respond your request')
       if(NX_m*NY_m*NKsplit/NXYsplit /= nproc_size_global) call err_finalize('NProcs is not good')
    
       NXY_s=NXYsplit*nproc_id_global/NKsplit
       NXY_e=(NXYsplit*(nproc_id_global+1)-1)/NKsplit
    
       allocate(NX_table(0:NX_m*NY_m-1),NY_table(0:NX_m*NY_m-1))
       i=-1
       do ix_m=1,NX_m
          do iy_m=1,NY_m
             i=i+1
             NX_table(i)=ix_m
             NY_table(i)=iy_m
          end do
       end do
    
       macRANK=NXY_s
       kRANK=mod(nproc_id_global,NKsplit)
    
       nproc_group_tdks = comm_create_group(nproc_group_global, macRANK, kRANK)
       call comm_get_groupinfo(nproc_group_tdks, nproc_id_tdks, nproc_size_tdks)
    end if
    
    !  NK_ave=NK/Nprocs; NK_remainder=NK-NK_ave*Nprocs
    !  NG_ave=NG/Nprocs; NG_remainder=NG-NG_ave*Nprocs
    
    NK_ave=NK/nproc_size_tdks; NK_remainder=NK-NK_ave*nproc_size_tdks
    NG_ave=NG/nproc_size_tdks; NG_remainder=NG-NG_ave*nproc_size_tdks
    
    if(is_symmetric_mode() == 1 .and. ENABLE_LOAD_BALANCER == 1) then
       call symmetric_load_balancing(NK,NK_ave,NK_s,NK_e,NK_remainder,nproc_id_tdks,nproc_size_tdks)
    else
       if (NK/nproc_size_tdks*nproc_size_tdks == NK) then
          NK_s=NK_ave*nproc_id_tdks+1
          NK_e=NK_ave*(nproc_id_tdks+1)
       else
          if (nproc_id_tdks < (nproc_size_tdks-1) - NK_remainder + 1) then
             NK_s=NK_ave*nproc_id_tdks+1
             NK_e=NK_ave*(nproc_id_tdks+1)
          else
             NK_s=NK-(NK_ave+1)*((nproc_size_tdks-1)-nproc_id_tdks)-NK_ave
             NK_e=NK-(NK_ave+1)*((nproc_size_tdks-1)-nproc_id_tdks)
          end if
       end if
       if(nproc_id_tdks == nproc_size_tdks-1 .and. NK_e /= NK) call err_finalize('prep. NK_e error')
    endif
    
    if (NG/nproc_size_tdks*nproc_size_tdks == NG) then
       NG_s=NG_ave*nproc_id_tdks+1
       NG_e=NG_ave*(nproc_id_tdks+1)
    else
       if (nproc_id_tdks < (nproc_size_tdks-1) - NG_remainder + 1) then
          NG_s=NG_ave*nproc_id_tdks+1
          NG_e=NG_ave*(nproc_id_tdks+1)
       else
          NG_s=NG-(NG_ave+1)*((nproc_size_tdks-1)-nproc_id_tdks)-NG_ave
          NG_e=NG-(NG_ave+1)*((nproc_size_tdks-1)-nproc_id_tdks)
       end if
    end if
    if(nproc_id_tdks == nproc_size_tdks-1 .and. NG_e /= NG) call err_finalize('prep. NG_e error')
    ! sato ---------------------------------------------------------------------------------------
    
    allocate(lap(-Nd:Nd),nab(-Nd:Nd))
    allocate(lapx(-Nd:Nd),lapy(-Nd:Nd),lapz(-Nd:Nd))
    allocate(nabx(-Nd:Nd),naby(-Nd:Nd),nabz(-Nd:Nd))
    allocate(Lx(NL),Ly(NL),Lz(NL),Gx(NG),Gy(NG),Gz(NG))
    allocate(Lxyz(0:NLx-1,0:NLy-1,0:NLz-1))
    allocate(ifdx(-Nd:Nd,1:NL),ifdy(-Nd:Nd,1:NL),ifdz(-Nd:Nd,1:NL))
    allocate(kAc(NK,3),kAc0(NK,3),kAc_new(NK,3))
    allocate(Vh(NL),Vexc(NL),Eexc(NL),rho(NL),Vpsl(NL),Vloc(NL),Vloc_GS(NL),Vloc_t(NL))
    allocate(Vloc_new(NL),Vloc_old(NL,2))
    !yabana
    allocate(tmass(NL),tjr(NL,3),tjr2(NL),tmass_t(NL),tjr_t(NL,3),tjr2_t(NL))
    !yabana
    allocate(rhoe_G(NG_s:NG_e),rhoion_G(NG_s:NG_e))
    allocate(rho_gs(NL))
    allocate(tpsi(NL),htpsi(NL),ttpsi(NL))
    allocate(tpsi_omp(NL,0:NUMBER_THREADS-1),htpsi_omp(NL,0:NUMBER_THREADS-1)) ! sato
    allocate(ttpsi_omp(NL,0:NUMBER_THREADS-1)) ! sato
    allocate(xk_omp(NL,0:NUMBER_THREADS-1),hxk_omp(NL,0:NUMBER_THREADS-1)) ! sato
    allocate(gk_omp(NL,0:NUMBER_THREADS-1),pk_omp(NL,0:NUMBER_THREADS-1)) ! sato
    allocate(pko_omp(NL,0:NUMBER_THREADS-1),txk_omp(NL,0:NUMBER_THREADS-1)) ! sato
    
    allocate(tau_s_l_omp(NL,0:NUMBER_THREADS-1),j_s_l_omp(NL,3,0:NUMBER_THREADS-1)) ! sato
    
    allocate(work(-Nd:NLx+Nd-1,-Nd:NLy+Nd-1,-Nd:NLz+Nd-1))
    allocate(zwork(-Nd:NLx+Nd-1,-Nd:NLy+Nd-1,-Nd:NLz+Nd-1))
    allocate(nxyz(0-NLx/2:NLx-1-NLx/2,0-NLy/2:NLy-1-NLy/2,0-NLz/2:NLz-1-NLz/2)) !Hartree
    allocate(rho_3D(0:NLx-1,0:NLy-1,0:NLz-1),Vh_3D(0:NLx-1,0:NLy-1,0:NLz-1))!Hartree
    allocate(rhoe_G_temp(1:NG),rhoe_G_3D(-NLx/2:NLx-1-NLx/2,-NLy/2:NLy-1-NLy/2,-NLz/2:NLz-1-NLz/2))!Hartree
    allocate(f1(0:NLx-1,0:NLy-1,-NLz/2:NLz-1-NLz/2),f2(0:NLx-1,-NLy/2:NLy-1-NLy/2,-NLz/2:NLz-1-NLz/2))!Hartree
    allocate(f3(-NLx/2:NLx-1-NLx/2,-NLy/2:NLy-1-NLy/2,0:NLz-1),f4(-NLx/2:Nlx-1-NLx/2,0:NLy-1,0:NLz-1))!Hartree
    allocate(eGx(-NLx/2:NLx-1-NLx/2,0:NLx-1),eGy(-NLy/2:NLy-1-NLy/2,0:NLy-1),eGz(-NLz/2:NLz-1-NLz/2,0:NLz-1))!Hartree
    allocate(eGxc(-NLx/2:NLx-1-NLx/2,0:NLx-1),eGyc(-NLy/2:NLy-1-NLy/2,0:NLy-1),eGzc(-NLz/2:NLz-1-NLz/2,0:NLz-1))!Hartree
    allocate(itable_sym(Sym,NL)) ! sym
    allocate(rho_l(NL),rho_tmp1(NL),rho_tmp2(NL)) !sym
    
    if (comm_is_root(nproc_id_global)) then
       write(*,*) 'NB,Nelec=',NB,Nelec
    endif
    if( kbTev < 0d0 )then ! sato
       NBoccmax=Nelec/2
    else 
       NBoccmax=NB
    end if
    

    call comm_bcast(NBoccmax,nproc_group_global)
    call comm_sync_all
    NKB=(NK_e-NK_s+1)*NBoccmax ! sato
    
    allocate(occ(NB,NK),wk(NK),esp(NB,NK))
    allocate(ovlp_occ_l(NB,NK),ovlp_occ(NB,NK))
    allocate(zu_GS(NL,NB,NK_s:NK_e),zu_GS0(NL,NB,NK_s:NK_e))
    allocate(zu_t(NL,NBoccmax,NK_s:NK_e))
    allocate(ik_table(NKB),ib_table(NKB)) ! sato
    allocate(esp_var(NB,NK))
    allocate(NBocc(NK)) !redistribution
    NBocc(:)=NBoccmax
    allocate(esp_vb_min(NK),esp_vb_max(NK)) !redistribution
    allocate(esp_cb_min(NK),esp_cb_max(NK)) !redistribution
    if (comm_is_root(nproc_id_global)) then
       write(*,*) 'FSset_option =',FSset_option
       write(*,*) 'Ncg=',Ncg
       write(*,*) 'Nmemory_MB,alpha_MB =',Nmemory_MB,alpha_MB
       write(*,*) 'NFSset_start,NFSset_every =',NFSset_start,NFSset_every
       write(*,*) 'Nscf=',Nscf
       !    write(*,*) 'ext_field =',ext_field
       !    write(*,*) 'Longi_Trans =',Longi_Trans
       write(*,*) 'use_ehrenfest_md =', use_ehrenfest_md
       write(*,*) 'projection_option =', projection_option
       write(*,*) 'Nt,dt=',Nt,dt
    endif
    call comm_sync_all
    !  if(ext_field /= 'LF' .and. ext_field /= 'LR' ) call err_finalize('incorrect option for ext_field')
    !  if(Longi_Trans /= 'Lo' .and. Longi_Trans /= 'Tr' ) call err_finalize('incorrect option for Longi_Trans')
    if(projection_option /= 'td' .and. projection_option /= 'gs' .and. &
         & projection_option /= 'no' ) &
         & call err_finalize('incorrect option for projection_option')
    
    call comm_sync_all
    
    allocate(javt(0:Nt,3))
    allocate(Ac_ext(-1:Nt+1,3),Ac_ind(-1:Nt+1,3),Ac_tot(-1:Nt+1,3))
    allocate(E_ext(0:Nt,3),E_ind(0:Nt,3),E_tot(0:Nt,3))
    
    NYvacB_m = 1
    NYvacT_m = NY_m  
    allocate(Ac_m(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
    allocate(Ac_old_m(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
    allocate(Ac_new_m(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
    allocate(g(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
    allocate(Elec(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(Bmag(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(j_m(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(jmatter_m(1:3,1:NX_m,1:NY_m))
    allocate(jmatter_m_l(1:3,1:NX_m,1:NY_m))
    jmatter_m_l=0d0;j_m=0d0
    
    allocate(zu_m(NL,NBoccmax,NK_s:NK_e,NXY_s:NXY_e))
    if(NXYsplit /= 1)then
       allocate(rho_m(NL,NXY_s:NXY_e))
       allocate(Vh_m(NL,NXY_s:NXY_e))
       allocate(Vexc_m(NL,NXY_s:NXY_e))
       allocate(Eexc_m(NL,NXY_s:NXY_e))
       allocate(Vloc_m(NL,NXY_s:NXY_e))
       allocate(Vloc_old_m(NL,2,NXY_s:NXY_e))
    end if
    
    allocate(energy_joule(NXvacL_m:NXvacR_m, NYvacB_m:NYvacT_m))
    allocate(energy_elec_Matter_l(1:NX_m,1:NY_m))
    allocate(energy_elec_Matter(1:NX_m,1:NY_m))
    allocate(energy_elec(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(energy_elemag(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(energy_total(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
    allocate(excited_electron_l(1:NX_m,1:NY_m))
    allocate(excited_electron(1:NX_m,1:NY_m))
    energy_elec_Matter_l(:,:)=0d0
    excited_electron_l=0d0
    Ndata_out = Nt / Nstep_write
    Ndata_out_per_proc = NData_out / nproc_size_global
    allocate(data_out(16,NXvacL_m:NXvacR_m,NY_m+1,0:Ndata_out_per_proc))
    allocate(data_local_Ac(3,NXY_s:NXY_e,0:Nt),data_local_jm(3,NXY_s:NXY_e,0:Nt))
    allocate(data_vac_Ac(3,2,0:Nt))
    ! sato ---------------------------------------------------------------------------------------
    
    call comm_sync_all
    
    
    allocate(Rps(NE),NRps(NE))
    allocate(Rion_eq(3,NI),dRion(3,NI,-1:Nt+1))
    allocate(Zps(NE),NRloc(NE),Rloc(NE),Mass(NE),force(3,NI))
    allocate(dVloc_G(NG_s:NG_e,NE),force_ion(3,NI))
    allocate(Mps(NI),Mlps(NE))
    allocate(anorm(0:Lmax,NE),inorm(0:Lmax,NE))
    allocate(rad(Nrmax,NE),vloctbl(Nrmax,NE),dvloctbl(Nrmax,NE))
    allocate(radnl(Nrmax,NE))
    allocate(udVtbl(Nrmax,0:Lmax,NE),dudVtbl(Nrmax,0:Lmax,NE))
    allocate(Floc(3,NI),Fnl(3,NI),Fion(3,NI))                         
    
    select case(iflag_atom_coor)
    case(ntype_atom_coor_cartesian)
       Rion_red(1,:)=Rion(1,:)/aLx
       Rion_red(2,:)=Rion(2,:)/aLy
       Rion_red(3,:)=Rion(3,:)/aLz
    case(ntype_atom_coor_reduced)
       Rion(1,:)=Rion_red(1,:)*aLx
       Rion(2,:)=Rion_red(2,:)*aLy
       Rion(3,:)=Rion_red(3,:)*aLz
    end select
    
    if (comm_is_root(nproc_id_global)) then
       write(*,*) 'Zatom=',(Zatom(j),j=1,NE)
       write(*,*) 'Lref=',(Lref(j),j=1,NE)
       write(*,*) 'i,Kion(ia)','(Rion_red(j,a),j=1,3)'
       do ia=1,NI
          write(*,*) ia,Kion(ia)
          write(*,'(3f12.8)') (Rion_red(j,ia),j=1,3)
       end do
    endif
    call comm_sync_all
    
    return
  End Subroutine Read_data
end module initialization
