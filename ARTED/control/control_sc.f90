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
!This file is "control_sc.f90"
!This file contains sc-mode program
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
module control_sc
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
  integer :: iter,ik,ib,ia,i,ixyz
  logical :: Rion_update
  character(10) :: functional_t
!$ integer :: omp_get_max_threads  

  call timer_initialize

  call load_environments

  if(comm_is_root()) then
    write(*,'(A)') "Welcome to SALMON-TDDFT singlecell mode"
    write(*,'(A)') "(Preliminary Developers Version)"
    write(*,'(2A)') "based on ARTED ver. = ",ARTED_ver
    call print_optimize_message
  end if

  NUMBER_THREADS=1
!$  NUMBER_THREADS=omp_get_max_threads()
!$  if(iter*0 == 0) then
!$    if(comm_is_root())write(*,*)'parallel = Hybrid'
!$  else
  if(comm_is_root())write(*,*)'parallel = Flat MPI'
!$  end if

  if(comm_is_root())write(*,*)'NUMBER_THREADS = ',NUMBER_THREADS

  call timer_begin(LOG_ALL)

  call timer_begin(LOG_STATIC)
  Time_start=get_wtime() !reentrance
  call comm_bcast(Time_start,proc_group(1))

  Rion_update = rion_update_on

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
!  call input_pseudopotential_KY !shinohara
!  call comm_finalize !!debug shinohara
!  stop !!debug shinohara
  call prep_ps_periodic('initial    ')

! initialize for optimization.
  call opt_vars_initialize_p2

  call Hartree
! yabana
  functional_t = functional
  if(functional_t == 'TBmBJ') functional = 'PZM'
  call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
  if(functional_t == 'TBmBJ') functional = 'TBmBJ'
! yabana
  Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
  call Total_Energy_omp(Rion_update,calc_mode_gs) ! debug
  call Ion_Force_omp(Rion_update,calc_mode_gs)
  if (MD_option /= 'Y') Rion_update = rion_update_off
  Eall_GS(0)=Eall

  if(comm_is_root()) then
    write(*,*) 'This is the end of preparation for ground state calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

  call reset_gs_timer
  call timer_begin(LOG_GROUND_STATE)
  do iter=1,Nscf
    if (comm_is_root())  write(*,*) 'iter = ',iter
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
            if (comm_is_root()) then
              write(*,*) '======================================='
              write(*,*) 'occupation redistribution is not needed'
              write(*,*) '======================================='
            end if
          end if
        end if
      end if
    else if( iter /= 1 )then ! sato
      call Fermi_Dirac_distribution
      if((comm_is_root()).and.(iter == Nscf))then
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
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
    if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'TBmBJ'
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

    if (comm_is_root()) then
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

  if(comm_is_root()) then
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
  end if
  if(comm_is_root()) then
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
  if(comm_is_root()) write(*,*) 'Eall =',Eall

  call timer_end(LOG_STATIC)
  if (comm_is_root()) then
    write(*,*) '-----------------------------------------------'
    call timer_show_min('static time=',LOG_STATIC)
    write(*,*) '-----------------------------------------------'
  end if

  if (comm_is_root()) then
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

  if(comm_is_root()) then
    write(*,*) 'This is the end of preparation for Real time calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation============================

  call init_Ac
  do ixyz=1,3
    kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter,ixyz)
  enddo
  call current0(zu_t)
  javt(0,:)=jav(:)

  Vloc_old(:,1) = Vloc(:); Vloc_old(:,2) = Vloc(:)
! yabana
!  kAc0=kAc
! yabana
  rho_gs(:)=rho(:)

!reentrance
2 if (entrance_option == 'reentrance') then
    position_option='append'
  else
    position_option='rewind'
    entrance_iter=-1
    call reset_rt_timer
  end if

  if (comm_is_root()) then
    open(7,file=file_epst,position = position_option)
    open(8,file=file_dns,position = position_option)
    open(9,file=file_force_dR,position = position_option)
    if (AD_RHO /= 'No') then 
      open(404,file=file_ovlp,position = position_option) 
      open(408,file=file_nex,position = position_option) 
    end if
  endif

  call comm_sync_all

!$acc enter data copyin(ik_table,ib_table)
!$acc enter data copyin(lapx,lapy,lapz)
!$acc enter data copyin(nabx,naby,nabz)
!$acc enter data copyin(modx,mody,modz)
!$acc enter data copyin(zJxyz,zKxyz)
!$acc enter data copyin(uV,iuV)

!$acc enter data create(kAc)

#ifdef ARTED_USE_PAPI
  call papi_begin
#endif

  call timer_begin(LOG_DYNAMICS)
!$acc enter data copyin(zu)
  do iter=entrance_iter+1,Nt

    if (Longi_Trans == 'Lo') then 
      Ac_ind(iter+1,:)=2*Ac_ind(iter,:)-Ac_ind(iter-1,:)-4*Pi*javt(iter,:)*dt**2
      if (Sym /= 1) then
        Ac_ind(iter+1,1)=0.d0
        Ac_ind(iter+1,2)=0.d0
      end if
      Ac_tot(iter+1,:)=Ac_ext(iter+1,:)+Ac_ind(iter+1,:)
    else if (Longi_Trans == 'Tr') then 
      Ac_tot(iter+1,:)=Ac_ext(iter+1,:)
    end if

    call dt_evolve_KB(iter,zu_t)

    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter+1,ixyz)
    enddo
!$acc update device(kAc,kAc_new)
    call current_RT(zu_t)

    javt(iter+1,:)=jav(:)
    if (MD_option == 'Y') then
!$acc update self(zu)
      call Ion_Force_omp(Rion_update,calc_mode_rt)
      if (iter/10*10 == iter) then
        call Total_Energy_omp(Rion_update,calc_mode_rt)
      end if
    else
      if (iter/10*10 == iter) then
!$acc update self(zu)
        call Total_Energy_omp(Rion_update,calc_mode_rt)
        call Ion_Force_omp(Rion_update,calc_mode_rt)
      end if
    end if

    call timer_begin(LOG_OTHER)


    E_ext(iter,:)=-(Ac_ext(iter+1,:)-Ac_ext(iter-1,:))/(2*dt)
    E_ind(iter,:)=-(Ac_ind(iter+1,:)-Ac_ind(iter-1,:))/(2*dt)
    E_tot(iter,:)=-(Ac_tot(iter+1,:)-Ac_tot(iter-1,:))/(2*dt)

    Eelemag=aLxyz*sum(E_tot(iter,:)**2)/(8.d0*Pi)
    Eall=Eall+Eelemag
    do ia=1,NI
      force_ion(:,ia)=Zps(Kion(ia))*E_tot(iter,:)
    enddo
    force=force+force_ion
!pseudo potential update
    if (MD_option == 'Y') then
      Tion=0.d0
      do ia=1,NI
        dRion(:,ia,iter+1)=2*dRion(:,ia,iter)-dRion(:,ia,iter-1)+force(:,ia)*dt**2/(umass*Mass(Kion(ia)))
        Rion(:,ia)=Rion_eq(:,ia)+dRion(:,ia,iter+1)
        Tion=Tion+0.5d0*umass*Mass(Kion(ia))*sum((dRion(:,ia,iter+1)-dRion(:,ia,iter-1))**2)/(2*dt)**2
      enddo
      call prep_ps_periodic('not_initial')
    else
      dRion(:,:,iter+1)=0.d0
      Tion=0.d0
    endif
    Eall=Eall+Tion
!write section
    if (iter/10*10 == iter.and.comm_is_root()) then
      write(*,'(1x,f10.4,8f12.6,f22.14)') iter*dt,&
           &E_ext(iter,1),E_tot(iter,1),&
           &E_ext(iter,2),E_tot(iter,2),&
           &E_ext(iter,3),E_tot(iter,3),&
           &Eall,Eall-Eall0,Tion
      write(7,'(1x,100e16.6E3)') iter*dt,&
           &E_ext(iter,1),E_tot(iter,1),&
           &E_ext(iter,2),E_tot(iter,2),&
           &E_ext(iter,3),E_tot(iter,3),&
           &Eall,Eall-Eall0,Tion
      write(9,'(1x,100e16.6E3)') iter*dt,((force(ixyz,ia),ixyz=1,3),ia=1,NI),((dRion(ixyz,ia,iter),ixyz=1,3),ia=1,NI)
    endif
!Dynamical Density
    if (iter/100*100 == iter.and.comm_is_root()) then
      write(8,'(1x,i10)') iter
      do i=1,NL
        write(8,'(1x,2e16.6E3)') rho(i),(rho(i)-rho_gs(i))
      enddo
    endif


!j_Ac writing
    if(comm_is_root())then
      if (iter/1000*1000 == iter .or. iter == Nt) then
        open(407,file=file_j_ac)
        do i=0,Nt
          write(407,'(100e26.16E3)') i*dt,javt(i,1),javt(i,2),javt(i,3),&
            &Ac_ext(i,1),Ac_ext(i,2),Ac_ext(i,3),Ac_tot(i,1),Ac_tot(i,2),Ac_tot(i,3)
        end do
        write(407,'(100e26.16E3)') (Nt+1)*dt,javt(Nt,1),javt(Nt,2),javt(Nt,3),&
          &Ac_ext(Nt+1,1),Ac_ext(Nt+1,2),Ac_ext(Nt+1,3),Ac_tot(Nt+1,1),Ac_tot(Nt+1,2),Ac_tot(Nt+1,3)
        close(407)
      end if
    end if
!Adiabatic evolution
    if (AD_RHO /= 'No' .and. iter/100*100 == iter) then
      call k_shift_wf(Rion_update,5,zu_t)
      if (comm_is_root()) then
        do ia=1,NI
          write(*,'(1x,i7,3f15.6)') ia,force(1,ia),force(2,ia),force(3,ia)
        end do
        write(404,'(1x,i10)') iter
        do ik=1,NK
          write(404,'(1x,i5,500e16.6)')ik,(ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
        enddo
        write(408,'(1x,3e16.6E3)') iter*dt,sum(ovlp_occ(NBoccmax+1:NB,:)),sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        write(*,*) 'number of excited electron',sum(ovlp_occ(NBoccmax+1:NB,:)),sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        write(*,*) 'var_tot,var_max=',sum(esp_var(:,:))/(NK*Nelec/2),maxval(esp_var(:,:)) 
      end if
    end if


!Timer
    if (iter/1000*1000 == iter.and.comm_is_root()) then
      call timer_show_current_hour('dynamics time      :', LOG_DYNAMICS)
      call timer_show_min         ('dt_evolve time     :', LOG_DT_EVOLVE)
      call timer_show_min         ('Hartree time       :', LOG_HARTREE)
      call timer_show_min         ('current time       :', LOG_CURRENT)
    end if
!Timer for shutdown
    if (iter/10*10 == iter) then
      Time_now=get_wtime()
      call comm_bcast(Time_now,proc_group(1))
      if (comm_is_root() .and. iter/100*100 == iter) then
        write(*,*) 'Total time =',(Time_now-Time_start)
      end if
      if ((Time_now - Time_start)>Time_shutdown) then 
        call comm_sync_all
        write(*,*) procid(1),'iter =',iter
        iter_now=iter
!$acc update self(zu)
        call timer_end(LOG_DYNAMICS)
        call prep_Reentrance_write
        go to 1
      end if
    end if

    call timer_end(LOG_OTHER)

    ! backup for system failure
    if (need_backup .and. iter > 0 .and. mod(iter, backup_frequency) == 0) then
      call timer_end(LOG_DYNAMICS)
      call timer_end(LOG_ALL)
      iter_now=iter
      call prep_Reentrance_write
      call timer_begin(LOG_ALL)
      call timer_begin(LOG_DYNAMICS)
    end if
  enddo !end of RT iteraction========================
!$acc exit data copyout(zu)
  call timer_end(LOG_DYNAMICS)

#ifdef ARTED_USE_PAPI
  call papi_end
#endif

  if(comm_is_root()) then
#ifdef ARTED_USE_PAPI
    call papi_result(timer_get(LOG_DYNAMICS))
#endif
    call timer_show_hour('dynamics time      :', LOG_DYNAMICS)
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
    call timer_show_min ('Other time         :', LOG_OTHER)
    call timer_show_min ('Allreduce time     :', LOG_ALLREDUCE)
  end if
  call write_performance(trim(directory)//'sc_performance')

  if(comm_is_root()) then
    close(7)
    close(8)
    close(9)
    if (AD_RHO /= 'No') then
      close(404)
      close(408)                                                      
    end if
  endif

  if(comm_is_root()) then
    write(*,*) 'This is the end of RT calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation===========================
!====Analyzing calculation====================

!Adiabatic evolution
  call k_shift_wf_last(Rion_update,10,zu_t)

  call Fourier_tr

  call comm_sync_all

  if (comm_is_root()) write(*,*) 'This is the end of all calculation'
  Time_now=get_wtime()
  call timer_end(LOG_ALL)
  if (comm_is_root()) call timer_show_hour('Total time =',LOG_ALL)

1 if(comm_is_root()) write(*,*)  'This calculation is shutdown successfully!'
  call comm_finalize

contains
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
  use opt_variables, only: symmetric_load_balancing, is_symmetric_mode
  use environment
  use communication
  use misc_routines
  implicit none
  integer :: ia,i,j

  if (comm_is_root()) then
    write(*,*) 'nprocs(1)=',nprocs(1)
    write(*,*) 'procid(1)=0:  ',procid(1)
    write(*,*) 'entrance_option=',entrance_option
    write(*,*) 'Time_shutdown=',Time_shutdown,'sec'
  end if

  if(entrance_option == 'reentrance') then
    call prep_Reentrance_Read
    return
  else if(entrance_option == 'new') then
  else 
    call err_finalize('entrance_option /= new or reentrance')
  end if

  if(comm_is_root())then

    need_backup = (backup_frequency > 0)
    write(*,*) 'need backup?',need_backup
    if (need_backup) write(*,*) '  frequency (# of iter) :',backup_frequency

    write(*,*) 'entrance_iter=',entrance_iter
    write(*,*) SYSname
    write(*,*) directory
    write(*,*) 'functional=',functional
    if(functional == 'TBmBJ') write(*,*) 'cvalue=',cval
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
    write(*,*) 'NEwald, aEwald =',NEwald, aEwald 
    write(*,*) 'KbTev=',KbTev ! sato
  end if


  write (process_directory,'(A,A,I5.5,A)') trim(directory),'/work_p',procid(1),'/'
  call create_directory(process_directory)

  call comm_bcast(need_backup,proc_group(1))
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

  ! Added by M.Uemoto on 2016-07-07
  if (0<NKx .and. 0<NKy .and. 0<NKy) then
    ! Use uniform rectangular k-grid
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
    ! Use non-uniform k-points
    if (comm_is_root()) then
      write(*,*) "Use non-uniform k-points distribution"
      write(*,*) "file_kw=", file_kw
      open(410, file=file_kw, status="old")
      read(410, *) NK, NKxyz
      close(410)
      write(*,*) "NK=", NK, "NKxyz=", NKxyz      
    endif
    call comm_bcast(NK,proc_group(1))
    call comm_bcast(NKxyz,proc_group(1))
  endif

  NK_ave=NK/nprocs(1); NK_remainder=NK-NK_ave*nprocs(1)
  NG_ave=NG/nprocs(1); NG_remainder=NG-NG_ave*nprocs(1)

  if(is_symmetric_mode() == 1 .and. ENABLE_LOAD_BALANCER == 1) then
    call symmetric_load_balancing(NK,NK_ave,NK_s,NK_e,NK_remainder,procid(1),nprocs(1))
  else
    if (NK/nprocs(1)*nprocs(1) == NK) then
      NK_s=NK_ave*procid(1)+1
      NK_e=NK_ave*(procid(1)+1)
    else
      if (procid(1) < (nprocs(1)-1) - NK_remainder + 1) then
        NK_s=NK_ave*procid(1)+1
        NK_e=NK_ave*(procid(1)+1)
      else
        NK_s=NK-(NK_ave+1)*((nprocs(1)-1)-procid(1))-NK_ave
        NK_e=NK-(NK_ave+1)*((nprocs(1)-1)-procid(1))
      end if
    end if
    if(procid(1) == nprocs(1)-1 .and. NK_e /= NK) call err_finalize('prep. NK_e error')
  end if

  if (NG/nprocs(1)*nprocs(1) == NG) then
    NG_s=NG_ave*procid(1)+1
    NG_e=NG_ave*(procid(1)+1)
  else
    if (procid(1) < (nprocs(1)-1) - NG_remainder + 1) then
      NG_s=NG_ave*procid(1)+1
      NG_e=NG_ave*(procid(1)+1)
    else
      NG_s=NG-(NG_ave+1)*((nprocs(1)-1)-procid(1))-NG_ave
      NG_e=NG-(NG_ave+1)*((nprocs(1)-1)-procid(1))
    end if
  end if
  if(procid(1) == nprocs(1)-1 .and. NG_e /= NG) call err_finalize('prep. NG_e error')

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
  allocate(zu_t(NL,NBoccmax,NK_s:NK_e))
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
    write(*,*) 'ext_field =',ext_field
    write(*,*) 'Longi_Trans =',Longi_Trans
    write(*,*) 'MD_option =', MD_option
    write(*,*) 'AD_RHO =', AD_RHO
    write(*,*) 'Nt,dt=',Nt,dt
  endif
  call comm_sync_all

  if(ext_field /= 'LF' .and. ext_field /= 'LR' ) call err_finalize('incorrect option for ext_field')
  if(Longi_Trans /= 'Lo' .and. Longi_Trans /= 'Tr' ) call err_finalize('incorrect option for Longi_Trans')
  if(AD_RHO /= 'TD' .and. AD_RHO /= 'GS' .and. AD_RHO /= 'No' ) call err_finalize('incorrect option for Longi_Trans')

  call comm_sync_all

  allocate(javt(0:Nt+1,3))
  allocate(Ac_ext(-1:Nt+1,3),Ac_ind(-1:Nt+1,3),Ac_tot(-1:Nt+1,3))
  allocate(E_ext(0:Nt,3),E_ind(0:Nt,3),E_tot(0:Nt,3))

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
end module control_sc
