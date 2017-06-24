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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module control_sc
  implicit none
contains
subroutine main
  use Global_Variables
  use timer
  use opt_variables
  use environment
  use performance_analyzer
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use misc_routines, only: get_wtime
  use salmon_global, only: format3d, out_dns, out_dns_rt, out_dns_rt_step
  use inputoutput, only: t_unit_time, t_unit_current, t_unit_ac
  use initialization
  implicit none
  integer :: iter,ik,ib,ia,i,ixyz
  logical :: Rion_update
  character(10) :: functional_t


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

  if (comm_is_root(nproc_id_global)) then
    open(7,file=file_epst,position = position_option)
    open(8,file=file_dns,position = position_option)
    open(9,file=file_force_dR,position = position_option)
    if (projection_option /= 'no') then 
      open(404,file=file_ovlp,position = position_option) 
      open(408,file=file_nex,position = position_option) 
    end if
    
    ! Export electronic density
    if (out_dns == 'y') then
      select case(format3d)
      case ('cube')
        write(file_dns_gs, '(2A,"_dns_gs.cube")') trim(directory), trim(SYSname)
        open(502,file=file_dns_gs,position = position_option)
        call write_density_cube(502, .false.)
        close(502)
      case ('vtk')
        write(file_dns_gs, '(2A,"_dns_gs.vtk")') trim(directory), trim(SYSname)
        open(502,file=file_dns_gs,position = position_option)
        call write_density_vtk(502, .false.)
        close(502)
      end select
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

    if (trans_longi == 'lo') then 
      Ac_ind(iter+1,:)=2*Ac_ind(iter,:)-Ac_ind(iter-1,:)-4*Pi*javt(iter,:)*dt**2
      if (Sym /= 1) then
        Ac_ind(iter+1,1)=0.d0
        Ac_ind(iter+1,2)=0.d0
      end if
      Ac_tot(iter+1,:)=Ac_ext(iter+1,:)+Ac_ind(iter+1,:)
    else if (trans_longi == 'tr') then 
      Ac_tot(iter+1,:)=Ac_ext(iter+1,:)
    end if

    call dt_evolve_KB(iter,zu_t)

    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter+1,ixyz)
    enddo
!$acc update device(kAc,kAc_new)
    call current_RT(zu_t)

    javt(iter+1,:)=jav(:)
    if (use_ehrenfest_md == 'y') then
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
    if (use_ehrenfest_md == 'y') then
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
    if (iter/10*10 == iter.and.comm_is_root(nproc_id_global)) then
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
    if (iter/100*100 == iter.and.comm_is_root(nproc_id_global)) then
      write(8,'(1x,i10)') iter
      do i=1,NL
        write(8,'(1x,2e16.6E3)') rho(i),(rho(i)-rho_gs(i))
      enddo
    endif
    
    ! Export electronic density
    if (out_dns_rt == 'y') then
      if (mod(iter, out_dns_rt_step) == 0) then
        select case(format3d)
        case ('cube')
          write(file_dns_rt, '(2A,"_dns_rt_",I6.6,".cube")') trim(directory), trim(SYSname), iter
          write(file_dns_dlt, '(2A,"_dns_dlt_",I6.6,".cube")') trim(directory), trim(SYSname), iter
          open(501,file=file_dns_rt,position = position_option)
          call write_density_cube(501, .false.)
          close(501)
          open(501,file=file_dns_dlt,position = position_option)
          call write_density_cube(501, .true.)
          close(501)
        case ('vtk')          
          write(file_dns_rt, '(2A,"_dns_rt_",I6.6,".vtk")') trim(directory), trim(SYSname), iter
          write(file_dns_dlt, '(2A,"_dns_dlt_",I6.6,".vtk")') trim(directory), trim(SYSname), iter
          open(501,file=file_dns_rt,position = position_option)
          call write_density_vtk(501, .false.)
          close(501)
          open(501,file=file_dns_dlt,position = position_option)
          call write_density_vtk(501, .true.)
          close(501)
        end select
      end if
    end if
    


!j_Ac writing
    if(comm_is_root(nproc_id_global))then
      if (iter/1000*1000 == iter .or. iter == Nt) then
        open(407,file=file_j_ac)
        write(407,'(A)')'# J     : Matter current density'
        write(407,'(A)')'# Ac_ext: External vector-potential devided by light-velocity'
        write(407,'(A)')'# Ac_tot: Total Vector-potential devided by light-velocity'
        write(407,'(A1,A25,9A26)')'#','Time ['//trim(t_unit_time%name)//']', &
             'Jx ['//trim(t_unit_current%name)//']', &
             'Jy ['//trim(t_unit_current%name)//']', &
             'Jz ['//trim(t_unit_current%name)//']', &
             'Ac_ext_x ['//trim(t_unit_ac%name)//']', &
             'Ac_ext_y ['//trim(t_unit_ac%name)//']', &
             'Ac_ext_z ['//trim(t_unit_ac%name)//']', &
             'Ac_tot_x ['//trim(t_unit_ac%name)//']', &
             'Ac_tot_y ['//trim(t_unit_ac%name)//']', &
             'Ac_tot_z ['//trim(t_unit_ac%name)//']'
        do i=0,Nt
          write(407,'(100e26.16E3)') i*dt*t_unit_time%conv, &
                                     javt(i,1)*t_unit_current%conv, &
                                     javt(i,2)*t_unit_current%conv, &
                                     javt(i,3)*t_unit_current%conv, &
                                     Ac_ext(i,1)*t_unit_ac%conv, &
                                     Ac_ext(i,2)*t_unit_ac%conv, &
                                     Ac_ext(i,3)*t_unit_ac%conv, &
                                     Ac_tot(i,1)*t_unit_ac%conv, &
                                     Ac_tot(i,2)*t_unit_ac%conv, &
                                     Ac_tot(i,3)*t_unit_ac%conv
        end do
!        write(407,'(100e26.16E3)') (Nt+1)*dt,javt(Nt,1),javt(Nt,2),javt(Nt,3),&
!          &Ac_ext(Nt+1,1),Ac_ext(Nt+1,2),Ac_ext(Nt+1,3),Ac_tot(Nt+1,1),Ac_tot(Nt+1,2),Ac_tot(Nt+1,3)
        close(407)
      end if
    end if
!Adiabatic evolution
    if (projection_option /= 'no' .and. iter/100*100 == iter) then
      call k_shift_wf(Rion_update,5,zu_t)
      if (comm_is_root(nproc_id_global)) then
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
    if (iter/1000*1000 == iter.and.comm_is_root(nproc_id_global)) then
      call timer_show_current_hour('dynamics time      :', LOG_DYNAMICS)
      call timer_show_min         ('dt_evolve time     :', LOG_DT_EVOLVE)
      call timer_show_min         ('Hartree time       :', LOG_HARTREE)
      call timer_show_min         ('current time       :', LOG_CURRENT)
    end if
!Timer for shutdown
    if (iter/10*10 == iter) then
      Time_now=get_wtime()
      call comm_bcast(Time_now,nproc_group_global)
      if (comm_is_root(nproc_id_global) .and. iter/100*100 == iter) then
        write(*,*) 'Total time =',(Time_now-Time_start)
      end if
      if ((Time_now - Time_start)>Time_shutdown .and. Time_shutdown >= 0d0) then 
        call comm_sync_all
        write(*,*) nproc_id_global,'iter =',iter
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

  if(comm_is_root(nproc_id_global)) then
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

  if(comm_is_root(nproc_id_global)) then
    close(7)
    close(8)
    close(9)
    if (projection_option /= 'no') then
      close(404)
      close(408)                                                      
    end if
  endif

  if(comm_is_root(nproc_id_global)) then
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

  if (comm_is_root(nproc_id_global)) write(*,*) 'This is the end of all calculation'
  Time_now=get_wtime()
  call timer_end(LOG_ALL)
  if (comm_is_root(nproc_id_global)) call timer_show_hour('Total time =',LOG_ALL)

1 if(comm_is_root(nproc_id_global)) write(*,*)  'This calculation is shutdown successfully!'

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
end module control_sc
