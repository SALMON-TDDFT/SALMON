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
    use io_gs_wfn_k,only: modify_initial_guess_copy_1stk_to_all
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

    if(use_ehrenfest_md=='y') then
       call init_md
    endif

! modify initial guess if read and modify options = y
    if(read_initial_guess=='y') then
       if(modify_initial_guess=='copy_1stk_to_all') then
          call modify_initial_guess_copy_1stk_to_all
       endif
    endif

  end subroutine initialize
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  Subroutine Read_data
    use salmon_global
    use Global_Variables
    use opt_variables
    use environment
    use salmon_parallel
    use salmon_communication
    use salmon_file
    use misc_routines
    use timer
    implicit none
    integer :: ia,i,j
    integer :: ix_m,iy_m,iz_m
    integer :: fh
    
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
       write(*,*) 'PSmask_option =',PSmask_option
       write(*,*) 'alpha_mask, gamma_mask, eta_mask =',real(alpha_mask), real(gamma_mask), real(eta_mask)
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
       file_k_data=trim(directory) // trim(SYSname) // '_k.data'
       file_eigen_data=trim(directory) // trim(SYSname) // '_eigen.data'
       file_rt_data=trim(directory) // trim(SYSname) // '_rt.data'
       
       write(*,*) 'al(1),al(2),al(3)=',real(al(1)),real(al(2)),real(al(3))
       write(*,*) 'Sym=',Sym,'crystal structure=',crystal_structure !sym
       write(*,*) 'Nd,NLx,NLy,NLz,NKx,NKy,NKz=',Nd,NLx,NLy,NLz,NKx,NKy,NKz
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
    call comm_bcast(file_k_data,nproc_group_global)
    call comm_bcast(file_eigen_data,nproc_group_global)
    call comm_bcast(file_rt_data,nproc_group_global)
    call comm_bcast(file_kw,nproc_group_global)

    if(use_ms_maxwell == 'y')then
      !! TODO: Modify the conditions for present implementation
      !  if(FDTDdim == '1D' .and. TwoD_shape /= 'periodic') then
      !     if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 1D calculation ! TwoD_shape is not good'
      !     TwoD_shape='periodic'
      !  end if
      !  if(FDTDdim == '1D' .and. NY_m /= 1) then
      !     if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 1D calculation ! NY_m is not good'
      !     NY_m=1
      !  end if
      !  if(FDTDdim == '2D' .and. TwoD_shape /= 'periodic') then
      !     if(comm_is_root(nproc_id_global))write(*,*)'Warning !! 2D calculation ! TwoD_shape is not good'
      !     TwoD_shape='periodic'
      !  end if
    end if
      
    !sym ---
    select case(crystal_structure)
    case("diamond2")
       if(functional == "PZ" .or. functional == "PZM" &
            .or. functional == "TBmBJ" .or. functional == "BJ_PW")then
          if(Sym == 8)then
             if((mod(NLx,2)+mod(NLy,2)+mod(NLz,4)) /= 0)call err_finalize('Bad grid point')
             if(NLx /= NLy) call err_finalize('Bad grid point: NLx /= NLy')
             if(NKx /= NKy) call err_finalize('NKx /= NKy')
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
             if(NLx /= NLy) call err_finalize('Bad grid point')
             if(NKx /= NKy) call err_finalize('NKx /= NKy')
          else if(Sym ==4 )then
             if(NLx /= NLy) call err_finalize('Bad grid point')
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
          write(*,*) "file_kw=", trim(file_kw)
          open(410, file=file_kw, status="old")
          read(410,*) NK, NKxyz
          close(410)
          write(*,*) "NK=", NK, "NKxyz=", NKxyz
       endif
       call comm_bcast(NK,nproc_group_global)
       call comm_bcast(NKxyz,nproc_group_global)
    endif
  
   !! Assign the number of macropoints into "nmacro"
   if (use_ms_maxwell == 'y') then
     !! Number of the macropoint and bg_media in Multiscale grid
     if (len_trim(file_macropoint) > 0) then
       call set_macropoint_from_file()
     else
       call set_macropoint()
     end if
     !! Determine NXYsplit and NKsplit from the number of MPI processes
     call set_nksplit_nxysplit()
     !! Output calculation condition
     if (comm_is_root(nproc_id_global)) then
       write(*,*) 'FDTDdim=',FDTDdim
       write(*,*) 'TwoD_shape=',TwoD_shape 
       write(*,*) 'NX_m,NY_m=',NX_m,NY_m
       write(*,*) 'HX_m,HY_m=',HX_m,HY_m
       write(*,*) 'NXvacL_m,NXvacR_m=',NXvacL_m,NXvacR_m
       write(*,*) 'NKsplit,NXYsplit=',NKsplit,NXYsplit
     end if
   else
     nmacro = 1; NKsplit = 1; NXYsplit = 1
   end if

    ! Create communicator "nproc_group_tdks"
    kRANK = mod(nproc_id_global, NKsplit)
    macRANK = (nproc_id_global - kRANK) / NKsplit
    nproc_group_tdks = comm_create_group(nproc_group_global, macRANK, kRANK)
    call comm_get_groupinfo(nproc_group_tdks, nproc_id_tdks, nproc_size_tdks)

    NK_ave=NK/nproc_size_tdks; NK_remainder=NK-NK_ave*nproc_size_tdks
    NG_ave=NG/nproc_size_tdks; NG_remainder=NG-NG_ave*nproc_size_tdks
    
1    if(is_symmetric_mode() == 1 .and. ENABLE_LOAD_BALANCER == 1) then
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
    if(read_initial_guess=='y') iflag_gs_init_wf=2

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
    
    allocate(javt(0:Nt+1,3))
    allocate(Ac_ext(-1:Nt+1,3),Ac_ind(-1:Nt+1,3),Ac_tot(-1:Nt+1,3))
    allocate(E_ext(0:Nt,3),E_ind(0:Nt,3),E_tot(0:Nt,3))
    
    !! Maxwell+TDDFT Multiscale Calculation:
    if (use_ms_maxwell == 'y') then
       !! Allocate multiscale variables
       call allocate_multiscale_vars()
    end if
    
    ! sato ---------------------------------------------------------------------------------------
    
    call comm_sync_all
    
    allocate(Rps(NE),NRps(NE))
    allocate(Rion_eq(3,NI),dRion(3,NI,-1:Nt+1))
    dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0
    allocate(Zps(NE),NRloc(NE),Rloc(NE),Mass(NE),force(3,NI))
    allocate(dVloc_G(NG_s:NG_e,NE),force_ion(3,NI))
    allocate(Mps(NI),Mlps(NE))
    allocate(anorm(0:Lmax,NE),inorm(0:Lmax,NE))
    allocate(rad(Nrmax,NE),vloctbl(Nrmax,NE),dvloctbl(Nrmax,NE))
    allocate(radnl(Nrmax,NE))
    allocate(udVtbl(Nrmax,0:Lmax,NE),dudVtbl(Nrmax,0:Lmax,NE))
    allocate(Floc(3,NI),Fnl(3,NI),Fion(3,NI))                         
    if(use_ehrenfest_md=='y')then
       allocate(velocity(3,NI)) ; velocity(:,:)=0d0
       allocate(save_dVloc_G(NG_s:NG_e,NE))
    endif
    
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
  Subroutine init_md
    use Global_Variables
    implicit none    

    if(restart_option == 'new') then
       if(set_ini_velocity=='y' .or. step_velocity_scaling>=1) &
       call set_initial_velocity
    endif

  End Subroutine init_md

  Subroutine set_initial_velocity
    use salmon_global
    use Global_Variables
    use salmon_parallel
    use salmon_communication
    !use misc_routines
    use salmon_math
    implicit none
    integer :: ia,ixyz,iseed
    real(8) :: rnd1,rnd2,rnd, sqrt_kT_im, kB,mass_au
    real(8) :: v_com(3), sum_mass, Temperature_ion, scale_v

    write(*,*) "  Initial velocities with maxwell-bolthman distribution was set"
    write(*,*) "  Set temperature is ", real(temperature0_ion)

    kB = 8.6173303d-5 / 27.211396d0 ![au/K]

    iseed= 123
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
       sqrt_kT_im = sqrt( kB * temperature0_ion / mass_au )

       do ixyz=1,3
          call quickrnd(iseed,rnd1)
          call quickrnd(iseed,rnd2)
          rnd = sqrt(-2d0*log(rnd1))*cos(2d0*Pi*rnd2)
          velocity(ixyz,ia) = rnd * sqrt_kT_im
       enddo
       !write(*,*)"velocity1:",ia,real(velocity(:,ia))
    enddo

    !!(check temperature)
    !Tion=0d0
    !do ia=1,NI
    !   Tion = Tion + 0.5d0*umass*Mass(Kion(ia))*sum(velocity(:,ia)**2d0)
    !enddo
    !Temperature_ion = Tion * 2d0 / (3d0*NI) / kB
    !write(*,*)"  Temperature: random-vel",real(Temperature_ion)

    !velocity of center of mass is removed
    v_com(:)=0d0
    sum_mass=0d0
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
       v_com(:) = v_com(:) + mass_au * velocity(:,ia)
       sum_mass = sum_mass + mass_au
    enddo
    v_com(:) = v_com(:)/sum_mass
    do ia=1,NI
       velocity(:,ia) = velocity(:,ia) - v_com(:)
    enddo

    !(check velocity of center of mass)
    v_com(:)=0d0
    do ia=1,NI
       v_com(:) = v_com(:) + umass*Mass(Kion(ia)) * velocity(:,ia)
    enddo
    v_com(:) = v_com(:)/sum_mass
    write(*,*)"    v_com =",real(v_com(:))

    !rotation around center of mass is removed (do nothing now)


    !scaling: set temperature exactly to input value
    Tion=0d0
    do ia=1,NI
       Tion = Tion + 0.5d0 * umass*Mass(Kion(ia)) * sum(velocity(:,ia)**2d0)
    enddo
    Temperature_ion = Tion * 2d0 / (3d0*NI) / kB
    !write(*,*)"    Temperature: befor-scaling",real(Temperature_ion)

    scale_v = sqrt(temperature0_ion/Temperature_ion)
    velocity(:,:) = velocity(:,:) * scale_v

    !(check)
    Tion=0d0
    do ia=1,NI
       Tion = Tion + 0.5d0 * umass*Mass(Kion(ia)) * sum(velocity(:,ia)**2d0)
    enddo
    Temperature_ion = Tion * 2d0 / (3d0*NI) / kB
    write(*,*)"    Initial Temperature: after-scaling",real(Temperature_ion)

    call comm_bcast(velocity ,nproc_group_global)

  End Subroutine set_initial_velocity
  
  
  
  subroutine set_nksplit_nxysplit()
    use Global_Variables
    use salmon_parallel
    implicit none
    !! Assign the macropoint into the every MPI procs    
    if (nproc_size_global <= nmacro) then
      !! Parallization Case 1:
      !! Assign more than one TDDFT cells in every single MPI processe
      if (mod(nmacro, nproc_size_global) == 0) then
        nksplit = 1
        nxysplit = nmacro / nproc_size_global
        nmacro_s = nxysplit * nproc_id_global + 1
        nmacro_e = nxysplit * (nproc_id_global + 1)
      else
        call err_finalize('Error! Set nproc as mod(num_macripoint, nproc) == 0')
      end if
      
    else !! (nproc_size_global > nmacro)
      !! Parallization Case 2:
      !! Fork the single TDDFT cell into more than two MPI processes
      if (mod(nproc_size_global, nmacro) == 0) then
        nksplit = nproc_size_global / nmacro
        nxysplit = 1
        nmacro_s = nproc_id_global / nksplit + 1
        nmacro_e = nmacro_s
      else
        call err_finalize('Error! Set nproc as mod(nproc % num_macripoint) == 0')
      end if
      
    end if
    
  end subroutine
  
  
  ! TODO: Create deallocate variables for the finalization of the program
  subroutine allocate_multiscale_vars()
    use Global_Variables
    use salmon_parallel
    implicit none
    
    !! Set the size of macroscopic grid
    nx1_m = min(NXvacL_m, nx_origin_m) 
    nx2_m = max(NXvacR_m, nx_origin_m + (nx_m - 1)) 
    ny1_m = ny_origin_m 
    ny2_m = ny_origin_m + (ny_m - 1) 
    nz1_m = nz_origin_m 
    nz2_m = nz_origin_m + (nz_m - 1)
    !! Set the actual size of macroscpic grid (including overlap region)
    mx1_m = nx1_m - novlp_m; mx2_m = nx2_m + novlp_m
    my1_m = ny1_m - novlp_m; my2_m = ny2_m + novlp_m
    mz1_m = nz1_m - novlp_m; mz2_m = nz2_m + novlp_m
  
    !! Allocate macroscopic electromagnetic field variables
    !! NOTE: "Ac_(old|new)?_ms" are the vector potential Ac(r,t)
    !!       "Jm_(old|new)?_ms" are the matter current density Jm(r,t)
    !!       In the RT iteration, the variable with suffix "new" and "old" 
    !!       indicate the data at the time "iter+1" and "iter-1", respectively.
    allocate(Ac_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Ac_old_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    allocate(Ac_new_ms(1:3, mx1_m:mx2_m, my1_m:my2_m, mz1_m:mz2_m))
    Ac_ms = 0d0; Ac_old_ms = 0d0; Ac_new_ms = 0d0
    allocate(Jm_ms(1:3, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    allocate(Jm_old_ms(1:3, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    allocate(Jm_new_ms(1:3, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    Jm_ms = 0d0; Jm_old_ms = 0d0; Jm_new_ms = 0d0
    allocate(elec_ms(1:3, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    allocate(bmag_ms(1:3, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    elec_ms = 0d0; bmag_ms = 0d0
    allocate(energy_joule_ms(nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    allocate(energy_elemag_ms(nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    allocate(energy_elec_ms(nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))
    energy_joule_ms = 0d0; energy_elemag_ms = 0d0;  energy_elec_ms = 0d0;
    !! Total energy of Entire System
    total_energy_elemag_old = 0d0; total_energy_elemag = 0d0
    total_energy_absorb_old = 0d0; total_energy_absorb = 0d0
    total_energy_elec_old = 0d0; total_energy_elec = 0d0
    total_energy_em_old = 0d0; total_energy_em  = 0d0
    !! Allocate local macropoint field variables
    allocate(Ac_m(1:3, 1:nmacro))
    allocate(Ac_new_m(1:3, 1:nmacro))
    allocate(Jm_m(1:3, nmacro))
    allocate(Jm_new_m(1:3, nmacro))
    allocate(jm_new_m_tmp(1:3, nmacro))
    Ac_m = 0d0; Ac_new_m = 0d0;
    Jm_m = 0d0; Jm_new_m = 0d0; jm_new_m_tmp = 0d0

    allocate(energy_elec_Matter_new_m(1:nmacro))
    allocate(energy_elec_Matter_new_m_tmp(1:nmacro))
    allocate(excited_electron_new_m(1:nmacro))
    allocate(excited_electron_new_m_tmp(1:nmacro))
    energy_elec_Matter_new_m = 0d0; energy_elec_Matter_new_m_tmp = 0d0
    excited_electron_new_m = 0d0; excited_electron_new_m_tmp = 0d0
    
    ndata_out = (nt / out_ms_step) + 1
    ndata_out_per_proc = ndata_out / nproc_size_global + 1
    
    allocate(data_out(1:ndata_out_column, nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m, 0:ndata_out_per_proc))
    allocate(data_local_Ac(3, nmacro_s:nmacro_e, 0:Nt))
    allocate(data_local_jm(3, nmacro_s:nmacro_e, 0:Nt))
    allocate(data_vac_Ac(3, 2, 0:Nt))
    
    ! Temporal Storage of Microscopic System
    allocate(zu_m(NL,NBoccmax,NK_s:NK_e,nmacro_s:nmacro_e))
    if(NXYsplit /= 1)then
      allocate(rho_m(NL, nmacro_s:nmacro_e))
      allocate(Vh_m(NL, nmacro_s:nmacro_e))
      allocate(Vexc_m(NL, nmacro_s:nmacro_e))
      allocate(Eexc_m(NL, nmacro_s:nmacro_e))
      allocate(Vloc_m(NL, nmacro_s:nmacro_e))
      allocate(Vloc_old_m(NL, 2, nmacro_s:nmacro_e))
    end if
    

  end subroutine allocate_multiscale_vars
  
  
  
  subroutine set_macropoint()
    use salmon_global
    use Global_variables
    implicit none
    integer :: ix_m, iy_m, iz_m, icount
    
    select case (FDTDdim)
    case("1D", "1d")
      nmacro = nx_m
      ny_m = 1; ny_origin_m = 1;
      nz_m = 1; nz_origin_m = 1;
    
    case("2D", "2d")
      nmacro = nx_m * ny_m
      nz_m = 1; nz_origin_m = 1;
    
    case("3D", "3d")
      nmacro = nx_m * ny_m * nz_m
    
    case default
      call err_finalize("Unknown FDTDdim:" // trim(FDTDdim))
  
    end select
      
    allocate(macropoint(1:4, nmacro))
    
    icount = 1
    do ix_m = 0, nx_m - 1
      do iy_m = 0, ny_m - 1
        do iz_m = 0, nz_m - 1
          macropoint(1, icount) = ix_m + nx_origin_m
          macropoint(2, icount) = iy_m + ny_origin_m
          macropoint(3, icount) = iz_m + nz_origin_m
          macropoint(4, icount) = 0
          icount = icount + 1
        end do
      end do
    end do    
      
    end subroutine
    
    subroutine set_macropoint_from_file()
      use salmon_file
      use salmon_parallel
      use salmon_communication
      use global_variables
      implicit none
      integer :: fh, icount, itmp

      namelist/macroscopic_system/ &
        & nx_origin_m, nx_m, hx_m, &
        & ny_origin_m, ny_m, hy_m, &
        & nz_origin_m, nz_m, hz_m, &
        & fdtddim, TwoD_shape, & 
        & nmacro, nmacro_attr, &
        & nbg_media, nbg_media_attr, &
        & ninit_acfield
      
      nmacro = 1
      nmacro_attr = 0
      nbg_media = 0
      nbg_media_attr = 0
      ninit_acfield = 0

      if(comm_is_root(nproc_id_global)) then
        fh = open_filehandle(trim(directory) // trim(file_macropoint))
        read(fh, nml=macroscopic_system)
      end if
      
      call comm_bcast(nx_origin_m,nproc_group_global)
      call comm_bcast(nx_m,nproc_group_global)
      call comm_bcast(hx_m,nproc_group_global)
      call comm_bcast(ny_origin_m,nproc_group_global)
      call comm_bcast(ny_m,nproc_group_global)
      call comm_bcast(hy_m,nproc_group_global)
      call comm_bcast(nz_origin_m,nproc_group_global)
      call comm_bcast(nz_m,nproc_group_global)
      call comm_bcast(hz_m,nproc_group_global)
      call comm_bcast(fdtddim,nproc_group_global)
      call comm_bcast(TwoD_shape,nproc_group_global)
      call comm_bcast(nmacro,nproc_group_global)
      call comm_bcast(nmacro_attr,nproc_group_global)
      call comm_bcast(nbg_media,nproc_group_global)
      call comm_bcast(nbg_media_attr,nproc_group_global)
      call comm_bcast(ninit_acfield,nproc_group_global)

      allocate(macropoint(1:4, nmacro))
      allocate(macropoint_attr(1:nattr_column, nmacro_attr))
      allocate(bg_media_point(1:4, nbg_media))
      allocate(bg_media_attr(1:nattr_column, nbg_media_attr))
      allocate(init_acfield_point(1:3, ninit_acfield))
      allocate(init_acfield_val(1:6, ninit_acfield))

      if(comm_is_root(nproc_id_global)) then
        do icount = 1, nmacro
          read(fh, *) itmp, macropoint(1:4, icount)
        end do
        do icount = 1, nmacro_attr
          read(fh, *) itmp, macropoint_attr(1:nattr_column, icount)
        end do
        do icount = 1, nbg_media
          read(fh, *) itmp, bg_media_point(1:4, icount)
        end do
        do icount = 1, nbg_media_attr
          read(fh, *) itmp, bg_media_attr(1:nattr_column, icount)
        end do
        do icount = 1, ninit_acfield
          read(fh, *) itmp, init_acfield_point(1:3, icount), &
            & init_acfield_val(1:6, icount)
        end do
        close(fh)
      end if
      call comm_bcast(macropoint,nproc_group_global)
      call comm_bcast(macropoint_attr,nproc_group_global)
      call comm_bcast(bg_media_point,nproc_group_global)
      call comm_bcast(bg_media_attr,nproc_group_global)
      call comm_bcast(init_acfield_point,nproc_group_global)
      call comm_bcast(init_acfield_val,nproc_group_global)
      return
  end subroutine set_macropoint_from_file

  
  
    
    
end module initialization
