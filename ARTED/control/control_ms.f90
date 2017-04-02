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
!This file is "control_ms.f90"
!This file contains ms-mode program
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
module control_ms
  implicit none
contains
subroutine main
  use Global_Variables
  use timer
  use opt_variables
  use environment
  use performance_analyzer
  use communication
  use misc_routines, only: get_wtime
  implicit none
  integer :: iter,ik,ib,ia
  character(3) :: Rion_update
  character(10) :: functional_t
  integer :: ix_m,iy_m,ixy_m
  integer :: index, n

!$ integer :: omp_get_max_threads  

  call timer_initialize

  call load_environments

  if(comm_is_root(1)) then
    write(*,'(A)') "Welcome to SALMON-TDDFT multiscale mode"
    write(*,'(A)') "(Preliminary Developers Version)"
    write(*,'(2A)') "based on ARTED ver. = ",ARTED_ver
    call print_optimize_message
  end if

  NUMBER_THREADS=1
!$  NUMBER_THREADS=omp_get_max_threads()
!$  if(iter*0 == 0) then
!$    if(comm_is_root(1))write(*,*)'parallel = Hybrid'
!$  else
  if(comm_is_root(1))write(*,*)'parallel = Flat MPI'
!$  end if

  if(comm_is_root(1))write(*,*)'NUMBER_THREADS = ',NUMBER_THREADS

  call timer_begin(LOG_ALL)

  call timer_begin(LOG_STATIC)
  Time_start=get_wtime() !reentrance
  call comm_bcast(Time_start,proc_group(1))

  Rion_update='on'

  call Read_data
  if (entrance_option == 'reentrance' ) go to 2

  allocate(rho_in(1:NL,1:Nscf+1),rho_out(1:NL,1:Nscf+1))
  rho_in(1:NL,1:Nscf+1)=0.d0; rho_out(1:NL,1:Nscf+1)=0.d0
  allocate(Eall_GS(0:Nscf),esp_var_ave(1:Nscf),esp_var_max(1:Nscf),dns_diff(1:Nscf))
  call fd_coef
  call init
  call init_wf
  call Gram_Schmidt
  rho=0.d0; Vh=0.d0
!  call psi_rho_omp !sym

! initialize for optimization.
  call opt_vars_initialize_p1

  call psi_rho_GS !sym
  rho_in(1:NL,1)=rho(1:NL)
  call input_pseudopotential_YS !shinohara
!  call input_pseudopotential_KY
  call prep_ps_periodic('initial    ')

! initialize for optimization.
  call opt_vars_initialize_p2

  call Hartree
! yabana
  functional_t = functional
  if(functional_t == 'TBmBJ') functional = 'PZM'
  call Exc_Cor('GS')
  if(functional_t == 'TBmBJ') functional = 'TBmBJ'
! yabana
  Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
!  call Total_Energy(Rion_update,'GS')
  call Total_Energy_omp(Rion_update,'GS') ! debug
  call Ion_Force_omp(Rion_update,'GS')
  if (MD_option /= 'Y') Rion_update = 'off'
  Eall_GS(0)=Eall

  if(comm_is_root(1)) then
    write(*,*) 'This is the end of preparation for ground state calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

  call reset_gs_timer
  call timer_begin(LOG_GROUND_STATE)
  do iter=1,Nscf
    if (comm_is_root(1))  write(*,*) 'iter = ',iter
    if( kbTev < 0d0 )then ! sato
      if (FSset_option == 'Y') then
        if (iter/NFSset_every*NFSset_every == iter .and. iter >= NFSset_start) then
          do ik=1,NK 
            esp_vb_min(ik)=minval(esp(1:NBocc(ik),ik))
            esp_vb_max(ik)=maxval(esp(1:NBocc(ik),ik))
            esp_cb_min(ik)=minval(esp(NBocc(ik)+1:NB,ik))
            esp_cb_max(ik)=maxval(esp(NBocc(ik)+1:NB,ik))
          end do
          if (minval(esp_cb_min(:))-maxval(esp_vb_max(:))<0.d0) then
            call Occupation_Redistribution
          else
            if (comm_is_root(1)) then
              write(*,*) '======================================='
              write(*,*) 'occupation redistribution is not needed'
              write(*,*) '======================================='
            end if
          end if
        end if
      end if
    else if( iter /= 1 )then ! sato
      call Fermi_Dirac_distribution
      if((comm_is_root(1)).and.(iter == Nscf))then
        open(126,file='occ.out')
        do ik=1,NK
          do ib=1,NB
            write(126,'(2I7,e26.16E3)')ik,ib,occ(ib,ik)
          end do
        end do
        close(126)
      end if
    end if
    call Gram_Schmidt
    call diag_omp
    call Gram_Schmidt
    call CG_omp(Ncg)
    call Gram_Schmidt

!    call psi_rho_omp !sym
    call psi_rho_GS
    call Density_Update(iter) 
    call Hartree
! yabana
    functional_t = functional
    if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'PZM'
    call Exc_Cor('GS')
    if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'TBmBJ'
! yabana
    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    call Total_Energy_omp(Rion_update,'GS')
    call Ion_Force_omp(Rion_update,'GS')
    call sp_energy_omp
    call current_GS
    Eall_GS(iter)=Eall
    esp_var_ave(iter)=sum(esp_var(:,:))/(NK*Nelec/2)
    esp_var_max(iter)=maxval(esp_var(:,:))
    dns_diff(iter)=sqrt(sum((rho_out(:,iter)-rho_in(:,iter))**2))*Hxyz

    if (comm_is_root(1)) then
      write(*,*) 'Total Energy = ',Eall_GS(iter),Eall_GS(iter)-Eall_GS(iter-1)
      write(*,'(a28,3e15.6)') 'jav(1),jav(2),jav(3)= ',jav(1),jav(2),jav(3)
      write(*,'(4(i3,f12.6,2x))') (ib,esp(ib,1),ib=1,NB)
      do ia=1,NI
        write(*,'(1x,i7,3f15.6)') ia,force(1,ia),force(2,ia),force(3,ia)
      end do
      write(*,*) 'var_ave,var_max=',esp_var_ave(iter),esp_var_max(iter)
      write(*,*) 'dns. difference =',dns_diff(iter)
      if (iter/20*20 == iter) then
         write(*,*) '====='
         call timer_show_current_min('elapse time=',LOG_ALL)
      end if
      write(*,*) '-----------------------------------------------'
    end if
  end do
  call timer_end(LOG_GROUND_STATE)

  if(comm_is_root(1)) then
    call timer_show_hour('Ground State time  :', LOG_GROUND_STATE)
    call timer_show_min ('CG time            :', LOG_CG)
    call timer_show_min ('Gram Schmidt time  :', LOG_GRAM_SCHMIDT)
    call timer_show_min ('diag time          :', LOG_DIAG)
    call timer_show_min ('sp_energy time     :', LOG_SP_ENERGY)
    call timer_show_min ('hpsi time          :', LOG_HPSI)
    call timer_show_min (' - stencil time    :', LOG_HPSI_STENCIL)
    call timer_show_min (' - pseudo pt. time :', LOG_HPSI_PSEUDO)
    call timer_show_min ('psi_rho time       :', LOG_PSI_RHO)
    call timer_show_min ('Hartree time       :', LOG_HARTREE)
    call timer_show_min ('Exc_Cor time       :', LOG_EXC_COR)
    call timer_show_min ('current time       :', LOG_CURRENT)
    call timer_show_min ('Total_Energy time  :', LOG_TOTAL_ENERGY)
    call timer_show_min ('Ion_Force time     :', LOG_ION_FORCE)
    call timer_show_min ('Allreduce time     :', LOG_ALLREDUCE)
  end if
  if(comm_is_root(1)) then
    write(*,*) 'This is the end of GS calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

  zu_GS0(:,:,:)=zu_GS(:,:,:)

  zu(:,:,:)=zu_GS(:,1:NBoccmax,:)
  Rion_eq=Rion
  dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0

!  call psi_rho_omp !sym
  call psi_rho_GS
  call Hartree
! yabana
  call Exc_Cor('GS')
! yabana
  Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
  Vloc_GS(:)=Vloc(:)
  call Total_Energy_omp(Rion_update,'GS')
  Eall0=Eall
  if(comm_is_root(1)) write(*,*) 'Eall =',Eall

  call timer_end(LOG_STATIC)
  if (comm_is_root(1)) then
    write(*,*) '-----------------------------------------------'
    call timer_show_min('static time=',LOG_STATIC)
    write(*,*) '-----------------------------------------------'
  end if

  if (comm_is_root(1)) then
    write(*,*) '-----------------------------------------------'
    write(*,*) '----some information for Band map--------------'
    do ik=1,NK 
      esp_vb_min(ik)=minval(esp(1:NBocc(ik),ik))
      esp_vb_max(ik)=maxval(esp(1:NBocc(ik),ik))
      esp_cb_min(ik)=minval(esp(NBocc(ik)+1:NB,ik))
      esp_cb_max(ik)=maxval(esp(NBocc(ik)+1:NB,ik))
    end do
    write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
    write(*,*) 'Top of VB',maxval(esp_vb_max(:))
    write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
    write(*,*) 'Top of CB',maxval(esp_cb_max(:))
    write(*,*) 'The Bandgap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(*,*) 'Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(*,*) '-----------------------------------------------'
    write(*,*) '-----------------------------------------------'
  end if

  call write_GS_data

  deallocate(rho_in,rho_out)
  deallocate(Eall_GS,esp_var_ave,esp_var_max,dns_diff)
!====GS calculation============================

#ifdef ARTED_LBLK
  call opt_vars_init_t4ppt()
#endif

  if(comm_is_root(1)) then
    write(*,*) 'This is the end of preparation for Real time calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation============================

!  call init_Ac
  if (trim(FDTDdim) == '2DC') then
    call init_Ac_ms_2dc()
  else
    call init_Ac_ms
  endif
  
  rho_gs(:)=rho(:)

  Vloc_old(:,1) = Vloc(:); Vloc_old(:,2) = Vloc(:)
! sato ---------------------------------------
  if(NXYsplit /= 1)then
    do ixy_m=NXY_s,NXY_e
      zu_m(:,:,:,ixy_m)=zu(:,1:NBoccmax,:)
      rho_m(:,ixy_m)=rho(:)
      Vh_m(:,ixy_m)=Vh(:)
      Vexc_m(:,ixy_m)=Vexc(:)
      Eexc_m(:,ixy_m)=Eexc(:)
      Vloc_m(:,ixy_m)=Vloc(:)
      Vloc_old_m(:,:,ixy_m)=Vloc_old(:,:)
    end do
  end if
! sato ---------------------------------------

!reentrance
2 if (entrance_option == 'reentrance') then
    position_option='append'
  else
    position_option='rewind'
    entrance_iter=-1
  end if
  
  ! Output filename
  write(file_energy_transfer, "(A,'energy-transfer.out')") trim(directory)
  write(file_ac_vac, "(A,'Ac_Vac.out')") trim(directory)
  write(file_ac_vac_back, "(A,'Ac_Vac_back.out')") trim(directory)
  write(file_ac_m, "(A,'Ac_M',I6.6,'.out')") trim(directory), NXY_s
  
  if (comm_is_root(1)) then
!    open(7,file=file_epst,position = position_option)
!    open(8,file=file_dns,position = position_option)
!    open(9,file=file_force_dR,position = position_option)
    open(940,file=file_energy_transfer, position = position_option)
    open(941,file=file_ac_vac, position = position_option)
!    if (ovlp_option == 'yes') then 
!      open(404,file=file_ovlp,position = position_option) 
!      open(408,file=file_nex,position = position_option) 
!    end if
  endif
  if (procid(1) == ROOT_PROCID+1) then
    open(942,file=file_ac_vac_back, position = position_option)
  end if

  if(comm_is_root(2))then
    open(943,file=file_ac_m ,position = position_option)
  end if

  call comm_sync_all

!$acc enter data copyin(ik_table,ib_table)
!$acc enter data copyin(lapx,lapy,lapz)
!$acc enter data copyin(nabx,naby,nabz)
!$acc enter data copyin(modx,mody,modz)
!$acc enter data copyin(zJxyz,zKxyz)
!$acc enter data copyin(uV,iuV)

!$acc enter data create(kAc)

  call reset_rt_timer
  call timer_begin(LOG_DYNAMICS)
!$acc enter data copyin(zu)
  RTiteratopm : do iter=entrance_iter+1,Nt ! sato

    call dt_evolve_Ac ! sato
    Macro_loop : do ixy_m=NXY_s,NXY_e ! sato
      call timer_begin(LOG_OTHER)
      ix_m=NX_table(ixy_m)
      iy_m=NY_table(ixy_m)
      if(NXYsplit /= 1)then
        call get_macro_data(ixy_m)
      end if
      call timer_end(LOG_OTHER)

      call dt_evolve_KB_MS(ix_m,iy_m)

      call timer_begin(LOG_OTHER)
! sato ---------------------------------------
      if(NXYsplit /= 1)then
        call put_macro_data(ixy_m)
      end if
      kAc(:,1)=kAc0(:,1)+Ac_new_m(1,ix_m,iy_m)
      kAc(:,2)=kAc0(:,2)+Ac_new_m(2,ix_m,iy_m)
      kAc(:,3)=kAc0(:,3)+Ac_new_m(3,ix_m,iy_m)
!$acc update device(kAc)
! sato ---------------------------------------
      call timer_end(LOG_OTHER)

      call current_RT

      call timer_begin(LOG_OTHER)
! sato ---------------------------------------
      if(Sym /= 1)then
        jav(1)=0d0
        jav(2)=0d0
      end if
      if(comm_is_root(2))then
        jmatter_m_l(1:3,ix_m,iy_m)=jav(1:3)
      end if
! sato ---------------------------------------
      call timer_end(LOG_OTHER)

      javt(iter,:)=jav(:)
      if (MD_option == 'Y') then
!$acc update self(zu)
        call Ion_Force_omp(Rion_update,'RT')
        if (mod(iter, Nstep_write) == 0) then
          call Total_Energy_omp(Rion_update,'RT')
        end if
      else
        if (mod(iter, Nstep_write) == 0) then
!$acc update self(zu)
          call Total_Energy_omp(Rion_update,'RT')
          call Ion_Force_omp(Rion_update,'RT')
        end if
      end if
    
      call timer_begin(LOG_OTHER)
      if(comm_is_root(2))then ! sato
        energy_elec_Matter_l(ix_m,iy_m)=Eall-Eall0 ! sato
      end if ! sato
      call timer_end(LOG_OTHER)

      call timer_begin(LOG_K_SHIFT_WF)
!Adiabatic evolution
      if (AD_RHO /= 'No' .and. mod(iter,100) == 0) then
        call k_shift_wf(Rion_update,2)
        if(comm_is_root(2))then ! sato
          excited_electron_l(ix_m,iy_m)=sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if ! sato
      else if (iter == Nt ) then
        call k_shift_wf(Rion_update,2)
        if(comm_is_root(2))then ! sato
          excited_electron_l(ix_m,iy_m)=sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if ! sato
      end if
      call timer_end(LOG_K_SHIFT_WF)
      
    end do Macro_loop

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(jmatter_m_l,jmatter_m,3*NX_m*NY_m,proc_group(1))
    j_m(:,1:NX_m,1:NY_m)=jmatter_m(:,1:NX_m,1:NY_m)
    if(mod(iter,10) == 1) then
      call comm_bcast(reentrance_switch,proc_group(1))
    end if
    call timer_end(LOG_ALLREDUCE)

    call timer_begin(LOG_OTHER)
!write section ================================================================================
    if(comm_is_root(1)) then
      write(941,'(4e26.16E3)') iter*dt, Ac_new_m(1:3,0,1)
    end if
    if(procid(1) == ROOT_PROCID+1) then
      ix_m=min(NXvacR_m,NX_m+1)
      write(942,'(4e26.16E3)') iter*dt, Ac_new_m(1:3,ix_m,1)
    end if
    if(comm_is_root(2)) then
      ix_m=NX_table(NXY_s)
      write(943,'(7e26.16E3)') iter*dt, Ac_new_m(1:3,ix_m,1), j_m(1:3,ix_m,1)
    end if
    
    call calc_elec_field()
    call calc_bmag_field()
    call calc_energy_joule()
    call calc_energy_elemag()
    
    if (mod(iter, Nstep_write) == 0) then

      call timer_end(LOG_OTHER)
      
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(energy_elec_Matter_l,energy_elec_Matter,NX_m*NY_m,proc_group(1))
      call timer_end(LOG_ALLREDUCE)

      call timer_begin(LOG_OTHER)

      energy_elec(1:NX_m,1:NY_m)=energy_elec_Matter(1:NX_m,1:NY_m) 
      energy_total=energy_elemag+energy_elec
      
      n = iter / Nstep_write
      if (mod(n, nprocs(1)) == procid(1)) then
        index = (n - procid(1)) / nprocs(1)
        data_out(1,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Ac_new_m(1,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(2,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Ac_new_m(2,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Ac_new_m(3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(4,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Elec(1,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(5,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Elec(2,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(6,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Elec(3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(7,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Bmag(1,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(8,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Bmag(2,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(9,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=Bmag(3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(10,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=j_m(1,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(11,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=j_m(2,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(12,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=j_m(3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(13,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=energy_elemag(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(14,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=energy_joule(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(15,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=energy_elec(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
        data_out(16,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m,index)=energy_total(NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m)
      end if
      
      if(comm_is_root(1))then
        write(940,'(4e26.16E3)')iter*dt,sum(energy_elec)*HX_m*HY_m/aLxyz &
          &,sum(energy_elemag)*HX_m*HY_m/aLxyz,sum(energy_total)*HX_m*HY_m/aLxyz
      end if
    end if
    call timer_end(LOG_OTHER)

    if (AD_RHO /= 'No' .and. mod(iter,100) == 0 ) then 
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(excited_electron_l,excited_electron,NX_m*NY_m,proc_group(1))
      call timer_end(LOG_ALLREDUCE)
      if(comm_is_root(1))call write_excited_electron(iter)
    else if (iter == Nt ) then
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(excited_electron_l,excited_electron,NX_m*NY_m,proc_group(1))
      call timer_end(LOG_ALLREDUCE)
      if(comm_is_root(1))call write_excited_electron(iter)
    end if

    call timer_begin(LOG_OTHER)
    if (reentrance_switch == 1) then 
      call comm_sync_all
      write(*,*) procid(1),'iter =',iter
      iter_now=iter
!$acc update self(zu)
      call prep_Reentrance_write
      go to 1
    end if

!Timer
    if (iter/1000*1000 == iter.and.comm_is_root(1)) then
      write(*,*) 'iter =',iter
      call timer_show_current_hour('dynamics time      :', LOG_DYNAMICS)
    end if

!Timer for shutdown
    if (mod(iter,10) == 0) then
      Time_now=get_wtime()
      if (comm_is_root(1) .and. iter/100*100 == iter) then
        write(*,*) 'Total time =',(Time_now-Time_start)
      end if
      if ((Time_now - Time_start)>Time_shutdown) then 
        reentrance_switch=1
      end if
    end if
! sato ---------------------------------------
    call timer_end(LOG_OTHER)

  enddo RTiteratopm !end of RT iteraction========================
!$acc exit data copyout(zu)
  call timer_end(LOG_DYNAMICS)

  if(comm_is_root(1)) then
    call timer_show_hour('dynamics time      :', LOG_DYNAMICS)
    call timer_show_min ('dt_evolve_Ac time  :', LOG_DT_EVOLVE_AC)
    call timer_show_min ('dt_evolve time     :', LOG_DT_EVOLVE)
    call timer_show_min ('hpsi time          :', LOG_HPSI)
    call timer_show_min (' - init time       :', LOG_HPSI_INIT)
    call timer_show_min (' - stencil time    :', LOG_HPSI_STENCIL)
    call timer_show_min (' - pseudo pt. time :', LOG_HPSI_PSEUDO)
    call timer_show_min (' - update time     :', LOG_HPSI_UPDATE)
    call timer_show_min ('psi_rho time       :', LOG_PSI_RHO)
    call timer_show_min ('Hartree time       :', LOG_HARTREE)
    call timer_show_min ('Exc_Cor time       :', LOG_EXC_COR)
    call timer_show_min ('current time       :', LOG_CURRENT)
    call timer_show_min ('Total_Energy time  :', LOG_TOTAL_ENERGY)
    call timer_show_min ('Ion_Force time     :', LOG_ION_FORCE)
    call timer_show_min ('k_shift_wf time    :', LOG_K_SHIFT_WF)
    call timer_show_min ('Other time         :', LOG_OTHER)
    call timer_show_min ('Allreduce time     :', LOG_ALLREDUCE)
  end if
  call write_performance(trim(directory)//'ms_performance')

  if(comm_is_root(1)) write(*,*) 'This is the start of write section'
  call timer_begin(LOG_IO)
  call write_result_all
  call timer_end(LOG_IO)
  if(comm_is_root(1)) then
    write(*,*) 'This is the end of write section'
    call timer_show_min('write time =',LOG_IO)
  end if

  if(comm_is_root(1)) then
    write(*,*) 'This is the end of RT calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation===========================
  call comm_sync_all

  if (comm_is_root(1)) write(*,*) 'This is the end of all calculation'
  Time_now=get_wtime()
  call timer_end(LOG_ALL)
  if (comm_is_root(1)) call timer_show_hour('Total time =',LOG_ALL)

1 if(comm_is_root(1)) write(*,*)  'This calculation is shutdown successfully!'
  if(comm_is_root(1)) then
    close(940)
    close(941)
!    close(7)
!    close(8)
!    close(9)
!    if (ovlp_option == 'yes') then
!      close(404)
!      close(408)                                                      
!    end if
  endif
  if(procid(1) == ROOT_PROCID+1) then
    close(942)
  end if
  if(comm_is_root(2)) then
    close(943)
  end if
  call comm_finalize

contains
  subroutine get_macro_data(ixy_m)
    implicit none
    integer, intent(in) :: ixy_m
    integer :: il,ib,ik
!$omp parallel default(none) &
!$    shared(NK_s,NK_e,NBoccmax,NL,zu,zu_m,Vh,Vh_m,Vexc, &
!$           Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m,Vloc_old,Vloc_old_m) &
!$    firstprivate(ixy_m)

!$omp do collapse(2) private(ik,ib)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
      zu(:,ib,ik) = zu_m(:,ib,ik,ixy_m)
    end do
    end do
!$omp end do nowait

!$omp do private(il)
    do il=1,NL
      Vh(il)         = Vh_m(il,ixy_m)
      Vexc(il)       = Vexc_m(il,ixy_m)
      Eexc(il)       = Eexc_m(il,ixy_m)
      Vloc(il)       = Vloc_m(il,ixy_m)
      Vloc_old(il,:) = Vloc_old_m(il,:,ixy_m)
    end do
!$omp end do

!$omp end parallel
  end subroutine

  subroutine put_macro_data(ixy_m)
    implicit none
    integer, intent(in) :: ixy_m
    integer :: il,ib,ik
!$omp parallel default(none) &
!$    shared(NK_s,NK_e,NBoccmax,NL,zu,zu_m,Vh,Vh_m,Vexc,Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m) &
!$    firstprivate(ixy_m)

!$omp do collapse(2) private(ik,ib)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
      zu_m(:,ib,ik,ixy_m) = zu(:,ib,ik)
    end do
    end do
!$omp end do nowait

!$omp do private(il)
    do il=1,NL
      Vh_m(il,ixy_m)   = Vh(il)
      Vexc_m(il,ixy_m) = Vexc(il)
      Eexc_m(il,ixy_m) = Eexc(il)
      Vloc_m(il,ixy_m) = Vloc(il)
    end do
!$omp end do

!$omp end parallel
  end subroutine

  subroutine reset_gs_timer
    implicit none
    integer :: i
    do i = LOG_CG,LOG_GRAM_SCHMIDT
      call timer_reset(i)
    end do
    call reset_rt_timer
  end subroutine

  subroutine reset_rt_timer
    implicit none
    integer :: i
    do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
      call timer_reset(i)
    end do
  end subroutine
end subroutine main
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_data
  use Global_Variables
  use opt_variables
  use environment
  use communication
  implicit none
  integer :: ia,i,j
  integer :: ix_m,iy_m

  if (comm_is_root()) then
    write(*,*) 'Nprocs=',nprocs(1)
    write(*,*) 'procid(1)=0:  ',procid(1)
    write(*,*) 'entrance_option=',entrance_option
    write(*,*) 'Time_shutdown=',Time_shutdown,'sec'
  end if

  if(entrance_option == 'reentrance') then
    call err_finalize('Re-entrance function does not work now!')
!    call prep_Reentrance_Read
    return
  else if(entrance_option == 'new') then
  else 
    call err_finalize('entrance_option /= new or reentrance')
  end if


  if(comm_is_root())then
    write(*,*) 'entrance_iter=',entrance_iter
    write(*,*) SYSname
    write(*,*) directory
!yabana
    write(*,*) 'functional=',functional
    if(functional == 'TBmBJ') write(*,*) 'cvalue=',cval
!yabana
    write(*,*) 'propagator=',propagator
    write(*,*) 'ps_format =',ps_format !shinohara
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
    write(*,*) 'aL,ax,ay,az=',aL,ax,ay,az
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

  call comm_bcast(file_GS,proc_group(1))
  call comm_bcast(file_RT,proc_group(1))
  call comm_bcast(file_epst,proc_group(1))
  call comm_bcast(file_epse,proc_group(1))
  call comm_bcast(file_force_dR,proc_group(1))
  call comm_bcast(file_j_ac,proc_group(1))  
  call comm_bcast(file_DoS,proc_group(1))
  call comm_bcast(file_band,proc_group(1))
  call comm_bcast(file_dns,proc_group(1))
  call comm_bcast(file_ovlp,proc_group(1))
  call comm_bcast(file_nex,proc_group(1))
  call comm_bcast(file_kw,proc_group(1))

  if(FDTDdim == '1D' .and. TwoD_shape /= 'periodic') then
     if(comm_is_root())write(*,*)'Warning !! 1D calculation ! TwoD_shape is not good'
     TwoD_shape='periodic'
  end if
  if(FDTDdim == '1D' .and. NY_m /= 1) then
     if(comm_is_root())write(*,*)'Warning !! 1D calculation ! NY_m is not good'
     NY_m=1
  end if
  if(FDTDdim == '2D' .and. TwoD_shape /= 'periodic') then
     if(comm_is_root())write(*,*)'Warning !! 2D calculation ! TwoD_shape is not good'
     TwoD_shape='periodic'
  end if

!sym ---
  select case(crystal_structure)
  case("diamond2")
     if(functional == "PZ" .or. functional == "PZM" .or. functional == "TBmBJ")then
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
     if(functional == "PZ" .or. functional == "PZM" .or. functional == "TBmBJ")then
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
  if(mod(NLx,2)+mod(NLy,2)+mod(NLz,2) /= 0) call err_finalize('NLx,NLy,NLz /= even')

  call comm_sync_all

  aLx=ax*aL;    aLy=ay*aL;    aLz=az*aL
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
    if (comm_is_root()) then
      write(*,*) "Use non-uniform k-points distribution"
      write(*,*) "file_kw=", file_kw
      open(410, file=file_kw, status="old")
      read(410,*) NK, NKxyz
      close(410)
      write(*,*) "NK=", NK, "NKxyz=", NKxyz
    endif
    call comm_bcast(NK,proc_group(1))
    call comm_bcast(NKxyz,proc_group(1))
  endif

! sato ---------------------------------------------------------------------------------------
  if(NXYsplit /= 1 .and. NKsplit /=1) call err_finalize('cannot respond your request')
  if(NX_m*NY_m*NKsplit/NXYsplit /= nprocs(1)) call err_finalize('NProcs is not good')

  NXY_s=NXYsplit*procid(1)/NKsplit
  NXY_e=(NXYsplit*(procid(1)+1)-1)/NKsplit

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
  kRANK=mod(procid(1),NKsplit)

  call comm_set_level2_group(macRANK, kRANK)

!  NK_ave=NK/Nprocs; NK_remainder=NK-NK_ave*Nprocs
!  NG_ave=NG/Nprocs; NG_remainder=NG-NG_ave*Nprocs

  NK_ave=NK/nprocs(2); NK_remainder=NK-NK_ave*nprocs(2)
  NG_ave=NG/nprocs(2); NG_remainder=NG-NG_ave*nprocs(2)

  if(is_symmetric_mode() == 1 .and. ENABLE_LOAD_BALANCER == 1) then
    call symmetric_load_balancing(NK,NK_ave,NK_s,NK_e,NK_remainder,procid(2),nprocs(2))
  else
    if (NK/nprocs(2)*nprocs(2) == NK) then
      NK_s=NK_ave*procid(2)+1
      NK_e=NK_ave*(procid(2)+1)
    else
      if (procid(2) < (nprocs(2)-1) - NK_remainder + 1) then
        NK_s=NK_ave*procid(2)+1
        NK_e=NK_ave*(procid(2)+1)
      else
        NK_s=NK-(NK_ave+1)*((nprocs(2)-1)-procid(2))-NK_ave
        NK_e=NK-(NK_ave+1)*((nprocs(2)-1)-procid(2))
      end if
    end if
    if(procid(2) == nprocs(2)-1 .and. NK_e /= NK) call err_finalize('prep. NK_e error')
  endif

  if (NG/nprocs(2)*nprocs(2) == NG) then
    NG_s=NG_ave*procid(2)+1
    NG_e=NG_ave*(procid(2)+1)
  else
    if (procid(2) < (nprocs(2)-1) - NG_remainder + 1) then
      NG_s=NG_ave*procid(2)+1
      NG_e=NG_ave*(procid(2)+1)
    else
      NG_s=NG-(NG_ave+1)*((nprocs(2)-1)-procid(2))-NG_ave
      NG_e=NG-(NG_ave+1)*((nprocs(2)-1)-procid(2))
    end if
  end if
  if(procid(2) == nprocs(2)-1 .and. NG_e /= NG) call err_finalize('prep. NG_e error')
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
  allocate(nxyz(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1)) !Hartree
  allocate(rho_3D(0:NLx-1,0:NLy-1,0:NLz-1),Vh_3D(0:NLx-1,0:NLy-1,0:NLz-1))!Hartree
  allocate(rhoe_G_temp(1:NG),rhoe_G_3D(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1))!Hartree
  allocate(f1(0:NLx-1,0:NLy-1,-NLz/2:NLz/2-1),f2(0:NLx-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1))!Hartree
  allocate(f3(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,0:NLz-1),f4(-NLx/2:NLx/2-1,0:NLy-1,0:NLz-1))!Hartree
  allocate(eGx(-NLx/2:NLx/2-1,0:NLx-1),eGy(-NLy/2:NLy/2-1,0:NLy-1),eGz(-NLz/2:NLz/2-1,0:NLz-1))!Hartree
  allocate(eGxc(-NLx/2:NLx/2-1,0:NLx-1),eGyc(-NLy/2:NLy/2-1,0:NLy-1),eGzc(-NLz/2:NLz/2-1,0:NLz-1))!Hartree
  allocate(itable_sym(Sym,NL)) ! sym
  allocate(rho_l(NL),rho_tmp1(NL),rho_tmp2(NL)) !sym

  if (comm_is_root()) then
    write(*,*) 'NB,Nelec=',NB,Nelec
  endif
  if( kbTev < 0d0 )then ! sato
    NBoccmax=Nelec/2
  else 
    NBoccmax=NB
  end if


  call comm_bcast(NBoccmax,proc_group(1))
  call comm_sync_all
  NKB=(NK_e-NK_s+1)*NBoccmax ! sato

  allocate(occ(NB,NK),wk(NK),esp(NB,NK))
  allocate(ovlp_occ_l(NB,NK),ovlp_occ(NB,NK))
  allocate(zu_GS(NL,NB,NK_s:NK_e),zu_GS0(NL,NB,NK_s:NK_e))
  allocate(zu(NL,NBoccmax,NK_s:NK_e))
  allocate(ik_table(NKB),ib_table(NKB)) ! sato
  allocate(esp_var(NB,NK))
  allocate(NBocc(NK)) !redistribution
  NBocc(:)=NBoccmax
  allocate(esp_vb_min(NK),esp_vb_max(NK)) !redistribution
  allocate(esp_cb_min(NK),esp_cb_max(NK)) !redistribution
  if (comm_is_root()) then
    write(*,*) 'FSset_option =',FSset_option
    write(*,*) 'Ncg=',Ncg
    write(*,*) 'Nmemory_MB,alpha_MB =',Nmemory_MB,alpha_MB
    write(*,*) 'NFSset_start,NFSset_every =',NFSset_start,NFSset_every
    write(*,*) 'Nscf=',Nscf
!    write(*,*) 'ext_field =',ext_field
!    write(*,*) 'Longi_Trans =',Longi_Trans
    write(*,*) 'MD_option =', MD_option
    write(*,*) 'AD_RHO =', AD_RHO
    write(*,*) 'Nt,dt=',Nt,dt
  endif
  call comm_sync_all
!  if(ext_field /= 'LF' .and. ext_field /= 'LR' ) call err_finalize('incorrect option for ext_field')
!  if(Longi_Trans /= 'Lo' .and. Longi_Trans /= 'Tr' ) call err_finalize('incorrect option for Longi_Trans')
  if(AD_RHO /= 'TD' .and. AD_RHO /= 'GS' .and. AD_RHO /= 'No' ) call err_finalize('incorrect option for Longi_Trans')

  call comm_sync_all

  allocate(javt(0:Nt,3))
  allocate(Ac_ext(-1:Nt+1,3),Ac_ind(-1:Nt+1,3),Ac_tot(-1:Nt+1,3))
  allocate(E_ext(0:Nt,3),E_ind(0:Nt,3),E_tot(0:Nt,3))

  NYvacB_m = 1
  NYvacT_m = NY_m  
  allocate(Ac_m(1:3,NXvacL_m-1:NXvacR_m+1,0:NY_m+1))
  allocate(Ac_old_m(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
  allocate(Ac_new_m(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
  allocate(g(1:3,NXvacL_m-1:NXvacR_m+1,NYvacB_m-1:NYvacT_m+1))
  allocate(Elec(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
  allocate(Bmag(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
  allocate(j_m(1:3,NXvacL_m:NXvacR_m,NYvacB_m:NYvacT_m))
  allocate(jmatter_m(1:3,1:NX_m,1:NY_m))
  allocate(jmatter_m_l(1:3,1:NX_m,1:NY_m))
  jmatter_m_l=0d0;j_m=0d0

  if(NXYsplit /= 1)then
    allocate(zu_m(NL,NBoccmax,NK_s:NK_e,NXY_s:NXY_e))         
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
    Ndata_out_per_proc = NData_out / nprocs(1)
    allocate(data_out(16,NXvacL_m:NXvacR_m,NY_m+1,0:Ndata_out_per_proc))
! sato ---------------------------------------------------------------------------------------

  if (comm_is_root()) then
    write(*,*) 'dAc=',dAc
    write(*,*) 'Nomega,etep=',Nomega,domega
    write(*,*) 'AE_shape=',AE_shape
    write(*,*) 'IWcm2_1, tpulsefs_1, omegaev_1, phi_CEP_1 =',IWcm2_1,tpulsefs_1,omegaev_1,phi_CEP_1
    write(*,*) 'Epdir_1(1), Epdir_1(2), Epdir_1(3) =', Epdir_1(1),Epdir_1(2),Epdir_1(3)
    write(*,*) 'IWcm2_2, tpulsefs_2, omegaev_2, phi_CEP_2 =',IWcm2_2,tpulsefs_2,omegaev_2,phi_CEP_2
    write(*,*) 'Epdir_2(1), Epdir_2(2), Epdir_2(3) =', Epdir_2(1),Epdir_2(2),Epdir_2(3)
    write(*,*) 'T1_T2fs =', T1_T2fs
    write(*,*) ''
    write(*,*) '===========ion configuration================'
    write(*,*) 'NI,NE=',NI,NE
  endif
  call comm_sync_all

  if(AE_shape /= 'Asin2cos' .and. AE_shape /= 'Esin2sin' &
    &.and. AE_shape /= 'input' .and. AE_shape /= 'Asin2_cw' ) call err_finalize('incorrect option for AE_shape')


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

  if (comm_is_root()) then
    write(*,*) 'Zatom=',(Zatom(j),j=1,NE)
    write(*,*) 'Lref=',(Lref(j),j=1,NE)
    write(*,*) 'i,Kion(ia)','(Rion(j,a),j=1,3)'
    do ia=1,NI
      write(*,*) ia,Kion(ia)
      write(*,'(3f12.8)') (Rion(j,ia),j=1,3)
    end do
  endif
  call comm_sync_all

  Rion(1,:)=Rion(1,:)*aLx
  Rion(2,:)=Rion(2,:)*aLy
  Rion(3,:)=Rion(3,:)*aLz

  return
End Subroutine Read_data
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine prep_Reentrance_Read
  use Global_Variables
  use timer,       only: timer_reentrance_read
  use opt_variables, only: opt_vars_initialize_p1, opt_vars_initialize_p2
  use communication
  use misc_routines, only: get_wtime
  implicit none
  real(8) :: time_in,time_out

  call comm_sync_all
  time_in=get_wtime()

  write(cMyrank,'(I5.5)')procid(1)
  file_reentrance=trim(directory)//'tmp_re.'//trim(cMyrank)
  open(500,file=file_reentrance,form='unformatted')

  read(500) iter_now,entrance_iter

!======== read section ===========================
!== read data ===!
! constants
!  real(8),parameter :: Pi=3.141592653589793d0
!  complex(8),parameter :: zI=(0.d0,1.d0)
!  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
!  real(8),parameter :: umass=1822.9d0

! DFT parameters
!  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
!  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
!  real(8),parameter :: CU=0.002d0,DU=-0.0116d0

! grid
  read(500) NLx,NLy,NLz,Nd,NL,NG,NKx,NKy,NKz,NK,Sym,nGzero
  read(500) NKxyz 
  read(500) aL,ax,ay,az,aLx,aLy,aLz,aLxyz
  read(500) bLx,bLy,bLz,Hx,Hy,Hz,Hxyz

! pseudopotential
!  integer,parameter :: Nrmax=3000,Lmax=4
  read(500) ps_type
  read(500) Nps,Nlma

! material
  read(500) NI,NE,NB,NBoccmax
  read(500) Ne_tot

! physical quantities
  read(500) Eall,Eall0,jav(3),Tion
  read(500) Ekin,Eloc,Enl,Eh,Exc,Eion,Eelemag                      

  read(500) Nelec !FS set

! Bloch momentum,laser pulse, electric field
!  real(8) :: f0,Wcm2,pulseT,wave_length,omega,pulse_time,pdir(3),phi_CEP=0.00*2*pi
  read(500) AE_shape
  read(500) f0_1,IWcm2_1,tpulsefs_1,omegaev_1,omega_1,tpulse_1,Epdir_1(3),phi_CEP_1 ! sato
  read(500) f0_2,IWcm2_2,tpulsefs_2,omegaev_2,omega_2,tpulse_2,Epdir_2(3),phi_CEP_2 ! sato
  read(500) T1_T2fs,T1_T2

! control parameters
  read(500) NEwald                      !Ewald summation
  read(500) aEwald
  read(500) Ncg                        !# of conjugate gradient (cg)
  read(500) dt,dAc,domega
  read(500) Nscf,Nt,Nomega
  read(500) Nmemory_MB                   !Modified-Broyden (MB) method
  read(500) alpha_MB
  read(500) NFSset_start,NFSset_every !Fermi Surface (FS) set 

! file names, flags, etc
  read(500) SYSname,directory
  read(500) file_GS,file_RT
  read(500) file_epst,file_epse
  read(500) file_force_dR,file_j_ac
  read(500) file_DoS,file_band
  read(500) file_dns,file_ovlp,file_nex
  read(500) ext_field
  read(500) Longi_Trans
  read(500) FSset_option,MD_option
  read(500) AD_RHO !ovlp_option

  read(500) propagator

!  integer :: procid(1),Nprocs
!  integer :: NEW_COMM_WORLD,nprocs(2),procid(2) ! sato
  read(500) NK_ave,NG_ave,NK_s,NK_e,NG_s,NG_e
  read(500) NK_remainder,NG_remainder
! Timer
!  read(500) Time_shutdown
!  read(500) Time_start,Time_now
!  read(500) iter_now,entrance_iter  !read by top
!  read(500) entrance_option    !initial or reentrance        
  read(500) position_option

  read(500) NKB
  read(500) crystal_structure !sym
  read(500) KbTev

  allocate(Lx(NL),Ly(NL),Lz(NL),Gx(NG),Gy(NG),Gz(NG))
  allocate(Lxyz(0:NLx-1,0:NLy-1,0:NLz-1))
  allocate(ifdx(-Nd:Nd,1:NL),ifdy(-Nd:Nd,1:NL),ifdz(-Nd:Nd,1:NL))
  allocate(lap(-Nd:Nd),nab(-Nd:Nd))
  allocate(lapx(-Nd:Nd),lapy(-Nd:Nd),lapz(-Nd:Nd))
  allocate(nabx(-Nd:Nd),naby(-Nd:Nd),nabz(-Nd:Nd))

  read(500) Lx(:),Ly(:),Lz(:),Lxyz(:,:,:)
  read(500) ifdx(:,:),ifdy(:,:),ifdz(:,:)
  read(500) Gx(:),Gy(:),Gz(:)
  read(500) lap(:),nab(:)
  read(500) lapx(:),lapy(:),lapz(:)
  read(500) nabx(:),naby(:),nabz(:)

  allocate(Mps(NI),Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI))
  allocate(Mlps(NE),Lref(NE),Zps(NE),NRloc(NE))
  allocate(NRps(NE),inorm(0:Lmax,NE),iuV(Nlma),a_tbl(Nlma))
  allocate(rad(Nrmax,NE),Rps(NE),vloctbl(Nrmax,NE),udVtbl(Nrmax,0:Lmax,NE))
  allocate(radnl(Nrmax,NE))
  allocate(Rloc(NE),uV(Nps,Nlma),duV(Nps,Nlma,3),anorm(0:Lmax,NE))
  allocate(dvloctbl(Nrmax,NE),dudVtbl(Nrmax,0:Lmax,NE))

  read(500) Mps(:),Jxyz(:,:),Jxx(:,:),Jyy(:,:),Jzz(:,:)
  read(500) Mlps(:),Lref(:),Zps(:),NRloc(:)
  read(500) NRps(:),inorm(:,:),iuV(:),a_tbl(:)
  read(500) rad(:,:),Rps(:),vloctbl(:,:),udVtbl(:,:,:)
  read(500) radnl(:,:)
  read(500) Rloc(:),uV(:,:),duV(:,:,:),anorm(:,:)
  read(500) dvloctbl(:,:),dudVtbl(:,:,:)


  allocate(Zatom(NE),Kion(NI))
  allocate(Rion(3,NI),Mass(NE),Rion_eq(3,NI),dRion(3,NI,-1:Nt+1))
  allocate(occ(NB,NK),wk(NK))
  
  read(500) Zatom(:),Kion(:)
  read(500) Rion(:,:),Mass(:),Rion_eq(:,:),dRion(:,:,:)
  read(500) occ(:,:),wk(:)



  allocate(javt(0:Nt,3))
  allocate(Vpsl(NL),Vh(NL),Vexc(NL),Eexc(NL),Vloc(NL),Vloc_GS(NL),Vloc_t(NL))
  allocate(dVloc_G(NG_s:NG_e,NE))
  allocate(rho(NL),rho_gs(NL))
  allocate(rhoe_G(NG_s:NG_e),rhoion_G(NG_s:NG_e))
  allocate(force(3,NI),esp(NB,NK),force_ion(3,NI))
  allocate(Floc(3,NI),Fnl(3,NI),Fion(3,NI))                         
  allocate(ovlp_occ_l(NB,NK),ovlp_occ(NB,NK))

  read(500) javt(:,:)
  read(500) Vpsl(:),Vh(:),Vexc(:),Eexc(:),Vloc(:),Vloc_GS(:),Vloc_t(:)
  read(500) dVloc_G(:,:)
  read(500) rho(:),rho_gs(:)
!  real(8),allocatable :: rho_in(:,:),rho_out(:,:) !MB method
  read(500) rhoe_G(:),rhoion_G(:)
  read(500) force(:,:),esp(:,:),force_ion(:,:)
  read(500) Floc(:,:),Fnl(:,:),Fion(:,:)               
  read(500) ovlp_occ_l(:,:),ovlp_occ(:,:)



  allocate(NBocc(NK)) !redistribution
  allocate(esp_vb_min(NK),esp_vb_max(NK)) !redistribution
  allocate(esp_cb_min(NK),esp_cb_max(NK)) !redistribution
!  allocate(Eall_GS(0:Nscf),esp_var_ave(1:Nscf),esp_var_max(1:Nscf),dns_diff(1:Nscf))

  read(500) NBocc(:) !FS set
  read(500) esp_vb_min(:),esp_vb_max(:) !FS set
  read(500) esp_cb_min(:),esp_cb_max(:) !FS set
!  read(500) Eall_GS(:),esp_var_ave(:),esp_var_max(:),dns_diff(:)


  allocate(zu(NL,NBoccmax,NK_s:NK_e),zu_GS(NL,NB,NK_s:NK_e),zu_GS0(NL,NB,NK_s:NK_e))
  allocate(tpsi(NL),htpsi(NL),zwork(-Nd:NLx+Nd-1,-Nd:NLy+Nd-1,-Nd:NLz+Nd-1),ttpsi(NL))
  allocate(work(-Nd:NLx+Nd-1,-Nd:NLy+Nd-1,-Nd:NLz+Nd-1))
  allocate(esp_var(NB,NK))

! wave functions, work array
  read(500) zu(:,:,:),zu_GS(:,:,:),zu_GS0(:,:,:)
!  read(500) tpsi(:),htpsi(:),zwork(:,:,:),ttpsi(:)
!  read(500) work(:,:,:)
  read(500) esp_var(:,:)


  allocate(nxyz(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1)) !Hartree
  allocate(rho_3D(0:NLx-1,0:NLy-1,0:NLz-1),Vh_3D(0:NLx-1,0:NLy-1,0:NLz-1))!Hartree
  allocate(rhoe_G_temp(1:NG),rhoe_G_3D(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1))!Hartree
  allocate(f1(0:NLx-1,0:NLy-1,-NLz/2:NLz/2-1),f2(0:NLx-1,-NLy/2:NLy/2-1,-NLz/2:NLz/2-1))!Hartree
  allocate(f3(-NLx/2:NLx/2-1,-NLy/2:NLy/2-1,0:NLz-1),f4(-NLx/2:NLx/2-1,0:NLy-1,0:NLz-1))!Hartree
  allocate(eGx(-NLx/2:NLx/2-1,0:NLx-1),eGy(-NLy/2:NLy/2-1,0:NLy-1),eGz(-NLz/2:NLz/2-1,0:NLz-1))!Hartree
  allocate(eGxc(-NLx/2:NLx/2-1,0:NLx-1),eGyc(-NLy/2:NLy/2-1,0:NLy-1),eGzc(-NLz/2:NLz/2-1,0:NLz-1))!Hartree

! variables for 4-times loop in Fourier transportation
  read(500) nxyz(:,:,:)
  read(500) rho_3D(:,:,:),Vh_3D(:,:,:)
  read(500) rhoe_G_temp(:),rhoe_G_3D(:,:,:)
  read(500) f1(:,:,:),f2(:,:,:),f3(:,:,:),f4(:,:,:)
  read(500) eGx(:,:),eGy(:,:),eGz(:,:),eGxc(:,:),eGyc(:,:),eGzc(:,:)



  allocate(E_ext(0:Nt,3),E_ind(0:Nt,3),E_tot(0:Nt,3))
  allocate(kAc(NK,3),kAc0(NK,3),kAc_new(NK,3))
  allocate(Ac_ext(-1:Nt+1,3),Ac_ind(-1:Nt+1,3),Ac_tot(-1:Nt+1,3))

  read(500) E_ext(:,:),E_ind(:,:),E_tot(:,:)
  read(500) kAc(:,:),kAc0(:,:),kAc_new(:,:)             !k+A(t)/c (kAc)
  read(500) Ac_ext(:,:),Ac_ind(:,:),Ac_tot(:,:) !A(t)/c (Ac)


  allocate(ekr(Nps,NI)) ! sato
  allocate(ekr_omp(Nps,NI,NK_s:NK_e))
  allocate(tpsi_omp(NL,0:NUMBER_THREADS-1),htpsi_omp(NL,0:NUMBER_THREADS-1)) ! sato
  allocate(ttpsi_omp(NL,0:NUMBER_THREADS-1)) ! sato
  allocate(xk_omp(NL,0:NUMBER_THREADS-1),hxk_omp(NL,0:NUMBER_THREADS-1)) ! sato
  allocate(gk_omp(NL,0:NUMBER_THREADS-1),pk_omp(NL,0:NUMBER_THREADS-1)) ! sato
  allocate(pko_omp(NL,0:NUMBER_THREADS-1),txk_omp(NL,0:NUMBER_THREADS-1)) ! sato
  allocate(ik_table(NKB),ib_table(NKB)) ! sato

! sato
  read(500) ekr(:,:)  
! omp
!  integer :: NUMBER_THREADS
  read(500) ekr_omp(:,:,:)
  read(500) tpsi_omp(:,:),ttpsi_omp(:,:),htpsi_omp(:,:)
  read(500) xk_omp(:,:),hxk_omp(:,:),gk_omp(:,:),pk_omp(:,:),pko_omp(:,:),txk_omp(:,:)
  read(500) ik_table(:),ib_table(:)


  allocate(itable_sym(Sym,NL)) ! sym
  allocate(rho_l(NL),rho_tmp1(NL),rho_tmp2(NL)) !sym

  read(500) itable_sym(:,:) ! sym
  read(500) rho_l(:),rho_tmp1(:),rho_tmp2(:) !sym

  call timer_reentrance_read(500)
  call opt_vars_initialize_p1
  call opt_vars_initialize_p2

!== read data ===!  

  close(500)

  call comm_sync_all
  time_out=get_wtime()

  if(comm_is_root())write(*,*)'Reentrance time read =',time_out-time_in,' sec'

  return
end subroutine prep_Reentrance_Read
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine prep_Reentrance_write
  use Global_Variables
  use timer, only: timer_reentrance_write
  use communication
  use misc_routines, only: get_wtime
  implicit none
  real(8) :: time_in,time_out

  call comm_sync_all
  time_in=get_wtime()

  if (comm_is_root()) then
    open(501,file=trim(directory)//trim(SYSname)//'_re.dat')
    write(501,*) "'reentrance'"! entrance_option
    write(501,*) Time_shutdown
    write(501,*) "'"//trim(directory)//"'"
    close(501)
  end if


  write(cMyrank,'(I5.5)')procid(1)

  file_reentrance=trim(directory)//'tmp_re.'//trim(cMyrank)
  open(500,file=file_reentrance,form='unformatted')

  write(500) iter_now,iter_now !iter_now=entrance_iter

  
!======== write section ===========================
!== write data ===!
! grid
  write(500) NLx,NLy,NLz,Nd,NL,NG,NKx,NKy,NKz,NK,Sym,nGzero
  write(500) NKxyz 
  write(500) aL,ax,ay,az,aLx,aLy,aLz,aLxyz
  write(500) bLx,bLy,bLz,Hx,Hy,Hz,Hxyz

! pseudopotential
!  integer,parameter :: Nrmax=3000,Lmax=4
  write(500) ps_type
  write(500) Nps,Nlma

! material
  write(500) NI,NE,NB,NBoccmax
  write(500) Ne_tot

! physical quantities
  write(500) Eall,Eall0,jav(3),Tion
  write(500) Ekin,Eloc,Enl,Eh,Exc,Eion,Eelemag                      

  write(500) Nelec !FS set

! Bloch momentum,laser pulse, electric field
!  real(8) :: f0,Wcm2,pulseT,wave_length,omega,pulse_time,pdir(3),phi_CEP=0.00*2*pi
  write(500) AE_shape
  write(500) f0_1,IWcm2_1,tpulsefs_1,omegaev_1,omega_1,tpulse_1,Epdir_1(3),phi_CEP_1 ! sato
  write(500) f0_2,IWcm2_2,tpulsefs_2,omegaev_2,omega_2,tpulse_2,Epdir_2(3),phi_CEP_2 ! sato
  write(500) T1_T2fs,T1_T2

! control parameters
  write(500) NEwald                      !Ewald summation
  write(500) aEwald
  write(500) Ncg                        !# of conjugate gradient (cg)
  write(500) dt,dAc,domega
  write(500) Nscf,Nt,Nomega
  write(500) Nmemory_MB                   !Modified-Broyden (MB) method
  write(500) alpha_MB
  write(500) NFSset_start,NFSset_every !Fermi Surface (FS) set 

! file names, flags, etc
  write(500) SYSname,directory
  write(500) file_GS,file_RT
  write(500) file_epst,file_epse
  write(500) file_force_dR,file_j_ac
  write(500) file_DoS,file_band
  write(500) file_dns,file_ovlp,file_nex
  write(500) ext_field
  write(500) Longi_Trans
  write(500) FSset_option,MD_option
  write(500) AD_RHO !ovlp_option

  write(500) propagator

!  integer :: procid(1),Nprocs
!  integer :: NEW_COMM_WORLD,nprocs(2),procid(2) ! sato
  write(500) NK_ave,NG_ave,NK_s,NK_e,NG_s,NG_e
  write(500) NK_remainder,NG_remainder
! Timer
!  write(500) Time_shutdown
!  write(500) Time_start,Time_now
!  write(500) iter_now,entrance_iter  !write by top
!  write(500) entrance_option    !initial or reentrance        
  write(500) position_option

  write(500) NKB
  write(500) crystal_structure !sym
  write(500) KbTev

  write(500) Lx(:),Ly(:),Lz(:),Lxyz(:,:,:)
  write(500) ifdx(:,:),ifdy(:,:),ifdz(:,:)
  write(500) Gx(:),Gy(:),Gz(:)
  write(500) lap(:),nab(:)
  write(500) lapx(:),lapy(:),lapz(:)
  write(500) nabx(:),naby(:),nabz(:)

  write(500) Mps(:),Jxyz(:,:),Jxx(:,:),Jyy(:,:),Jzz(:,:)
  write(500) Mlps(:),Lref(:),Zps(:),NRloc(:)
  write(500) NRps(:),inorm(:,:),iuV(:),a_tbl(:)
  write(500) rad(:,:),Rps(:),vloctbl(:,:),udVtbl(:,:,:)
  write(500) radnl(:,:)
  write(500) Rloc(:),uV(:,:),duV(:,:,:),anorm(:,:)
  write(500) dvloctbl(:,:),dudVtbl(:,:,:)


  write(500) Zatom(:),Kion(:)
  write(500) Rion(:,:),Mass(:),Rion_eq(:,:),dRion(:,:,:)
  write(500) occ(:,:),wk(:)



  write(500) javt(:,:)
  write(500) Vpsl(:),Vh(:),Vexc(:),Eexc(:),Vloc(:),Vloc_GS(:),Vloc_t(:)
  write(500) dVloc_G(:,:)
  write(500) rho(:),rho_gs(:)
!  real(8),allocatable :: rho_in(:,:),rho_out(:,:) !MB method
  write(500) rhoe_G(:),rhoion_G(:)
  write(500) force(:,:),esp(:,:),force_ion(:,:)
  write(500) Floc(:,:),Fnl(:,:),Fion(:,:)               
  write(500) ovlp_occ_l(:,:),ovlp_occ(:,:)



  write(500) NBocc(:) !FS set
  write(500) esp_vb_min(:),esp_vb_max(:) !FS set
  write(500) esp_cb_min(:),esp_cb_max(:) !FS set
!  write(500) Eall_GS(:),esp_var_ave(:),esp_var_max(:),dns_diff(:)


! wave functions, work array
  write(500) zu(:,:,:),zu_GS(:,:,:),zu_GS0(:,:,:)
!  write(500) tpsi(:),htpsi(:),zwork(:,:,:),ttpsi(:)
!  write(500) work(:,:,:)
  write(500) esp_var(:,:)


! variables for 4-times loop in Fourier transportation
  write(500) nxyz(:,:,:)
  write(500) rho_3D(:,:,:),Vh_3D(:,:,:)
  write(500) rhoe_G_temp(:),rhoe_G_3D(:,:,:)
  write(500) f1(:,:,:),f2(:,:,:),f3(:,:,:),f4(:,:,:)
  write(500) eGx(:,:),eGy(:,:),eGz(:,:),eGxc(:,:),eGyc(:,:),eGzc(:,:)


  write(500) E_ext(:,:),E_ind(:,:),E_tot(:,:)
  write(500) kAc(:,:),kAc0(:,:),kAc_new(:,:)                  !k+A(t)/c (kAc)
  write(500) Ac_ext(:,:),Ac_ind(:,:),Ac_tot(:,:) !A(t)/c (Ac)


! sato
  write(500) ekr(:,:)  
! omp
!  integer :: NUMBER_THREADS
  write(500) ekr_omp(:,:,:)
  write(500) tpsi_omp(:,:),ttpsi_omp(:,:),htpsi_omp(:,:)
  write(500) xk_omp(:,:),hxk_omp(:,:),gk_omp(:,:),pk_omp(:,:),pko_omp(:,:),txk_omp(:,:)
  write(500) ik_table(:),ib_table(:)


  write(500) itable_sym(:,:) ! sym
  write(500) rho_l(:),rho_tmp1(:),rho_tmp2(:) !sym

  call timer_reentrance_write(500)

!== write data ===!  

  close(500)
  call comm_sync_all
  time_out=get_wtime()

  if(comm_is_root())write(*,*)'Reentrance time write =',time_out-time_in,' sec'

  return
end subroutine prep_Reentrance_write
end module control_ms
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
