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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module control_ms
  implicit none
contains
subroutine main
  use Global_Variables
  use timer
  use opt_variables
  use environment
  use performance_analyzer
  use salmon_parallel
  use salmon_communication
  use misc_routines
  use initialization
  implicit none
  integer :: iter,ik,ib,ia
  logical :: Rion_update
  character(10) :: functional_t
  integer :: ix_m,iy_m,ixy_m
  integer :: index, n
  character(len=128) :: fmt
  
  real(8) calc_pulse_xcenter

  call initialize




  if (entrance_option == 'reentrance' ) go to 2

  Rion_update = rion_update_on

  allocate(rho_in(1:NL,1:Nscf+1),rho_out(1:NL,1:Nscf+1))
  rho_in(1:NL,1:Nscf+1)=0.d0; rho_out(1:NL,1:Nscf+1)=0.d0
  allocate(Eall_GS(0:Nscf),esp_var_ave(1:Nscf),esp_var_max(1:Nscf),dns_diff(1:Nscf))

  call init_wf
  call Gram_Schmidt
  rho=0.d0; Vh=0.d0

  call psi_rho_GS !sym
  rho_in(1:NL,1)=rho(1:NL)


  call Hartree
! yabana
  functional_t = functional
  if(functional_t == 'TBmBJ' .or. functional_t == 'BJ_PW') functional = 'PZM'
  call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
  if(functional_t == 'TBmBJ') functional = 'TBmBJ'
  if(functional_t == 'BJ_PW') functional = 'BJ_PW'

! yabana
  Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
!  call Total_Energy(Rion_update,calc_mode_gs)
  call Total_Energy_omp(Rion_update,calc_mode_gs) ! debug
  call Ion_Force_omp(Rion_update,calc_mode_gs)
  if (use_ehrenfest_md /= 'y') Rion_update = rion_update_off
  Eall_GS(0)=Eall

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of preparation for ground state calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

  call reset_gs_timer
  call timer_begin(LOG_GROUND_STATE)
  do iter=1,Nscf
    if (comm_is_root(nproc_id_global))  write(*,*) 'iter = ',iter
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
            if (comm_is_root(nproc_id_global)) then
              write(*,*) '======================================='
              write(*,*) 'occupation redistribution is not needed'
              write(*,*) '======================================='
            end if
          end if
        end if
      end if
    else if( iter /= 1 )then ! sato
      call Fermi_Dirac_distribution
      if((comm_is_root(nproc_id_global)).and.(iter == Nscf))then
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
    if(functional_t == 'BJ_PW' .and. iter < 20) functional = 'PZM'
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
    if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'TBmBJ'
    if(functional_t == 'BJ_PW' .and. iter < 20) functional = 'BJ_PW'
! yabana
    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    call Total_Energy_omp(Rion_update,calc_mode_gs)
    call Ion_Force_omp(Rion_update,calc_mode_gs)
    call sp_energy_omp
    call current_GS
    Eall_GS(iter)=Eall
    esp_var_ave(iter)=sum(esp_var(:,:))/(NK*Nelec/2)
    esp_var_max(iter)=maxval(esp_var(:,:))
    dns_diff(iter)=sqrt(sum((rho_out(:,iter)-rho_in(:,iter))**2))*Hxyz

    if (comm_is_root(nproc_id_global)) then
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

  if(comm_is_root(nproc_id_global)) then
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
  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of GS calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

  zu_GS0(:,:,:)=zu_GS(:,:,:)

  zu_t(:,:,:)=zu_GS(:,1:NBoccmax,:)
  Rion_eq=Rion
  dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0

!  call psi_rho_omp !sym
  call psi_rho_GS
  call Hartree
! yabana
  call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
! yabana
  Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
  Vloc_GS(:)=Vloc(:)
  call Total_Energy_omp(Rion_update,calc_mode_gs)
  Eall0=Eall
  if(comm_is_root(nproc_id_global)) write(*,*) 'Eall =',Eall

  call timer_end(LOG_STATIC)
  if (comm_is_root(nproc_id_global)) then
    write(*,*) '-----------------------------------------------'
    call timer_show_min('static time=',LOG_STATIC)
    write(*,*) '-----------------------------------------------'
  end if

  if (comm_is_root(nproc_id_global)) then
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

  if(comm_is_root(nproc_id_global)) then
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
  do ixy_m=NXY_s,NXY_e
    zu_m(:,:,:,ixy_m)=zu_t(:,1:NBoccmax,:)
  end do

  if (NXYsplit /= 1) then
    do ixy_m=NXY_s,NXY_e
      rho_m(:,ixy_m)=rho(:)
      Vh_m(:,ixy_m)=Vh(:)
      Vexc_m(:,ixy_m)=Vexc(:)
      Eexc_m(:,ixy_m)=Eexc(:)
      Vloc_m(:,ixy_m)=Vloc(:)
      Vloc_old_m(:,:,ixy_m)=Vloc_old(:,:)
    end do
  end if

  deallocate(zu_t)
! sato ---------------------------------------

!reentrance
2 if (entrance_option == 'reentrance') then
    position_option='asis'
    if (use_ehrenfest_md /= 'y') Rion_update = rion_update_off
  else
    position_option='rewind'
    entrance_iter=-1
    call reset_rt_timer
  end if
  
  ! Output filename
  write(file_energy_transfer, "(A,'energy-transfer.out')") trim(directory)
  write(file_ac_vac, "(A,'Ac_Vac.out')") trim(directory)
  write(file_ac_vac_back, "(A,'Ac_Vac_back.out')") trim(directory)
  write(file_ac_m, "(A,'Ac_M',I6.6,'.out')") trim(process_directory), NXY_s
  
!  if (comm_is_root(nproc_id_global)) then
!    open(940,file=file_energy_transfer, position = position_option)
!  endif

!$acc enter data copyin(ik_table,ib_table)
!$acc enter data copyin(lapx,lapy,lapz)
!$acc enter data copyin(nabx,naby,nabz)
!$acc enter data copyin(modx,mody,modz)
!$acc enter data copyin(zJxyz,zKxyz)
!$acc enter data copyin(uV,iuV)

!$acc enter data create(kAc)

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

      call dt_evolve_KB_MS(ixy_m)

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

      call current_RT_MS(ixy_m)

      call timer_begin(LOG_OTHER)
! sato ---------------------------------------
      if(Sym /= 1)then
        jav(1)=0d0
        jav(2)=0d0
      end if
      if(comm_is_root(nproc_id_tdks))then
        jmatter_m_l(1:3,ix_m,iy_m)=jav(1:3)
      end if
! sato ---------------------------------------
      call timer_end(LOG_OTHER)

      javt(iter,:)=jav(:)
      if (use_ehrenfest_md == 'y') then
!$acc update self(zu)
        call Ion_Force_omp(Rion_update,calc_mode_rt,ixy_m)
        if (mod(iter, Nstep_write) == 0) then
          call Total_Energy_omp(Rion_update,calc_mode_rt,ixy_m)
        end if
      else
        if (mod(iter, Nstep_write) == 0) then
!$acc update self(zu)
          call Total_Energy_omp(Rion_update,calc_mode_rt,ixy_m)
          call Ion_Force_omp(Rion_update,calc_mode_rt,ixy_m)
        end if
      end if
    
      call timer_begin(LOG_OTHER)
      if(comm_is_root(nproc_id_tdks))then ! sato
        energy_elec_Matter_l(ix_m,iy_m)=Eall-Eall0 ! sato
      end if ! sato
      call timer_end(LOG_OTHER)

      call timer_begin(LOG_K_SHIFT_WF)
!Adiabatic evolution
      if (projection_option /= 'no' .and. mod(iter,100) == 0) then
        call k_shift_wf(Rion_update,2,zu_m(:,:,:,ixy_m))
        if(comm_is_root(nproc_id_tdks))then ! sato
          excited_electron_l(ix_m,iy_m)=sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if ! sato
      else if (iter == Nt ) then
        call k_shift_wf(Rion_update,2,zu_m(:,:,:,ixy_m))
        if(comm_is_root(nproc_id_tdks))then ! sato
          excited_electron_l(ix_m,iy_m)=sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if ! sato
      end if
      call timer_end(LOG_K_SHIFT_WF)
      
    end do Macro_loop

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(jmatter_m_l,jmatter_m,3*NX_m*NY_m,nproc_group_global)
    j_m(:,1:NX_m,1:NY_m)=jmatter_m(:,1:NX_m,1:NY_m)
    if(mod(iter,10) == 1) then
      call comm_bcast(reentrance_switch,nproc_group_global)
    end if
    call timer_end(LOG_ALLREDUCE)

    call timer_begin(LOG_OTHER)
!write section ================================================================================
    if(comm_is_root(nproc_id_global)) then
      ix_m=min(NXvacR_m,NX_m+1)
      data_vac_Ac(1:3,1,iter) = Ac_new_m(1:3,0,1)
      data_vac_Ac(1:3,2,iter) = Ac_new_m(1:3,ix_m,1)
    end if
    if(comm_is_root(nproc_id_tdks)) then
      do ixy_m=NXY_s,NXY_e
        ix_m=NX_table(ixy_m)
        iy_m=NY_table(ixy_m)

        data_local_Ac(1:3,ixy_m,iter) = Ac_new_m(1:3,ix_m,iy_m)
        data_local_jm(1:3,ixy_m,iter) = j_m(1:3,ix_m,iy_m)
      end do
    end if
    
    call calc_elec_field()
    call calc_bmag_field()
    call calc_energy_joule()
    call calc_energy_elemag()
    
    if (mod(iter, Nstep_write) == 0) then

      call timer_end(LOG_OTHER)
      
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(energy_elec_Matter_l,energy_elec_Matter,NX_m*NY_m,nproc_group_global)
      call timer_end(LOG_ALLREDUCE)

      call timer_begin(LOG_OTHER)

      energy_elec(1:NX_m,1:NY_m)=energy_elec_Matter(1:NX_m,1:NY_m) 
      energy_total=energy_elemag+energy_elec
      
      n = iter / Nstep_write
      if (mod(n, nproc_size_global) == nproc_id_global) then
        index = (n - nproc_id_global) / nproc_size_global
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
      
!      if(comm_is_root(nproc_id_global))then
!        write(940,'(4e26.16E3)')iter*dt,sum(energy_elec)*HX_m*HY_m/aLxyz &
!          &,sum(energy_elemag)*HX_m*HY_m/aLxyz,sum(energy_total)*HX_m*HY_m/aLxyz
!      end if
    end if
    call timer_end(LOG_OTHER)

    if (projection_option /= 'no' .and. mod(iter,100) == 0 ) then 
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(excited_electron_l,excited_electron,NX_m*NY_m,nproc_group_global)
      call timer_end(LOG_ALLREDUCE)
      if(comm_is_root(nproc_id_global))call write_excited_electron(iter)
    else if (iter == Nt ) then
      call timer_begin(LOG_ALLREDUCE)
      call comm_summation(excited_electron_l,excited_electron,NX_m*NY_m,nproc_group_global)
      call timer_end(LOG_ALLREDUCE)
      if(comm_is_root(nproc_id_global))call write_excited_electron(iter)
    end if

    call timer_begin(LOG_OTHER)
    if (reentrance_switch == 1) then 
      call comm_sync_all
      write(*,*) nproc_id_global,'iter =',iter
      iter_now=iter
!$acc update self(zu)
      call timer_end(LOG_DYNAMICS)
      call prep_Reentrance_write
      go to 1
    end if

!Timer
    if ((mod(iter, 1000) == 0) .and. comm_is_root(nproc_id_global)) then
      write(*,*) 'iter =', iter
      write(*,*) "pulse_xcenter =", calc_pulse_xcenter() 
      call timer_show_current_hour('dynamics time      :', LOG_DYNAMICS)
    end if
    

!Timer for shutdown
    if (mod(iter,10) == 0) then
      Time_now=get_wtime()
      if (comm_is_root(nproc_id_global) .and. iter/100*100 == iter) then
        write(*,*) 'Total time =',(Time_now-Time_start)
      end if
      if ((Time_now - Time_start)>Time_shutdown .and. Time_shutdown >= 0d0) then 
        reentrance_switch=1
      end if
    end if
! sato ---------------------------------------
    call timer_end(LOG_OTHER)

    ! backup for system failure
    if (need_backup .and. iter > 0 .and. mod(iter, backup_frequency) == 0) then
      if (comm_is_root(nproc_id_global)) call timer_show_current_hour('Backup...', LOG_ALL)
      call timer_end(LOG_DYNAMICS)
      call timer_end(LOG_ALL)
      iter_now=iter
      call prep_Reentrance_write
      call timer_begin(LOG_ALL)
      call timer_begin(LOG_DYNAMICS)
    end if
  enddo RTiteratopm !end of RT iteraction========================
!$acc exit data copyout(zu)
  call timer_end(LOG_DYNAMICS)

  if(comm_is_root(nproc_id_global)) then
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

  if(comm_is_root(nproc_id_global)) write(*,*) 'This is the start of write section'
  call timer_begin(LOG_IO)
  call write_result_all

  if (comm_is_root(nproc_id_global)) then
    open(941,file=file_ac_vac, position = position_option)
    do iter=0,Nt
       write(941,"(7e26.16e3)")iter*dt,data_vac_Ac(1:3,1,iter) &
            ,data_vac_Ac(1:3,2,iter)
    end do
    close(941)
  end if



  if(comm_is_root(nproc_id_tdks))then
    write (fmt,"(A,I2,A)")"(",(NXY_e-NXY_s+1)*6+1,"e26.16e3)"
    open(943,file=file_ac_m ,position = position_option)
    write(943,"(A,2x,I6,2x,A,2x,I6)")"# Data of macro points",NXY_s,"-",NXY_e
    do iter=0,Nt
       write(943,fmt)iter*dt,(data_local_Ac(1:3,ixy_m,iter) &
            ,data_local_jm(1:3,ixy_m,iter),ixy_m = NXY_s,NXY_e)
    end do
    close(943)
  end if


  call timer_end(LOG_IO)
  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of write section'
    call timer_show_min('write time =',LOG_IO)
  end if

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of RT calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation===========================
  call comm_sync_all

  if (comm_is_root(nproc_id_global)) write(*,*) 'This is the end of all calculation'
  Time_now=get_wtime()
  call timer_end(LOG_ALL)
  if (comm_is_root(nproc_id_global)) call timer_show_hour('Total time =',LOG_ALL)

1 if(comm_is_root(nproc_id_global)) write(*,*)  'This calculation is shutdown successfully!'
!  if(comm_is_root(nproc_id_global)) then
!    close(940)
!  endif
  !call comm_finalize

contains
  subroutine get_macro_data(ixy_m)
    implicit none
    integer, intent(in) :: ixy_m
    integer :: il
!$omp parallel do &
!$omp&    default(none) private(il) firstprivate(ixy_m) &
!$omp&    shared(NL,Vh,Vh_m,Vexc,Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m,Vloc_old,Vloc_old_m)
    do il=1,NL
      Vh(il)         = Vh_m(il,ixy_m)
      Vexc(il)       = Vexc_m(il,ixy_m)
      Eexc(il)       = Eexc_m(il,ixy_m)
      Vloc(il)       = Vloc_m(il,ixy_m)
      Vloc_old(il,:) = Vloc_old_m(il,:,ixy_m)
    end do
!$omp end parallel do
  end subroutine

  subroutine put_macro_data(ixy_m)
    implicit none
    integer, intent(in) :: ixy_m
    integer :: il
!$omp parallel do &
!$omp&    default(none) private(il) firstprivate(ixy_m) &
!$omp&    shared(NL,Vh,Vh_m,Vexc,Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m)
    do il=1,NL
      Vh_m(il,ixy_m)   = Vh(il)
      Vexc_m(il,ixy_m) = Vexc(il)
      Eexc_m(il,ixy_m) = Eexc(il)
      Vloc_m(il,ixy_m) = Vloc(il)
    end do
!$omp end parallel do
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
end module control_ms
