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
!This file is "control_sc.f90"
!This file contains sc-mode program
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module control_sc
  implicit none
contains
subroutine tddft_sc
  use Global_Variables
  use timer
  use opt_variables
  use performance_analyzer
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use misc_routines, only: get_wtime
  use salmon_global, only: format3d, out_dns, out_dns_rt, out_dns_rt_step
  use salmon_file, only: open_filehandle
  use inputoutput, only: t_unit_time, t_unit_current, t_unit_ac,  t_unit_energy, t_unit_elec
  use restart, only: prep_restart_write
  use Ac_alocal_laser
  use md_ground_state
  implicit none
  integer :: iter,ia,i,ixyz  !,ib,ik
  real(8) :: Temperature_ion,mass_au,Ework,dt_h, aforce(3,NI)
  real(8) :: Enh, Enh_gkTlns, gkT, Qnh
  character(100) :: comment_line

#ifdef ARTED_LBLK
  call opt_vars_init_t4ppt()
#endif


!reentrance
2 if (restart_option == 'restart') then
    position_option='append'
  else if(restart_option == 'new')then
    position_option='rewind'
    entrance_iter=-1
  end if

  select case(use_ehrenfest_md)
  case('y')
    Rion_update_rt = rion_update_on
  case('n')
    Rion_update_rt = rion_update_off
  end select

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of preparation for Real time calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation============================
  if ( .not. (restart_option == 'restart')) call init_Ac
  iter=entrance_iter+1
  do ixyz=1,3
    kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter,ixyz)
  enddo
  call current0(zu_t)
  javt(0,:)=jav(:)

  if(use_ehrenfest_md=='y') then
     call current_RT_ion
     javt_ion(0,:)=jav_ion(:)  
  endif

  Vloc_old(:,1) = Vloc(:); Vloc_old(:,2) = Vloc(:)

  rho_gs(:)=rho(:)

  call reset_rt_timer


  ! Export electronic density
  if (out_dns == 'y') call write_density(iter,'gs')

  ! Export transition electronic density at specified energy
  if (out_dns_trans == 'y') call analysis_dns_trans(iter)

  if (comm_is_root(nproc_id_global)) then
    ! open(7,file=file_epst,    position = position_option) !! TODO: remove output of "_t.out" file future
    if (projection_option /= 'no') then 
      open(404, file=file_ovlp,         position = position_option) 
      open(408, file=file_nex,          position = position_option) 
      open(409, file=file_last_band_map,position = position_option) 
      if(projection_option=='gs' .and. projection_decomp=='atom') &
    & open(420, file=file_nex_atom,  position = position_option) 
      call write_projection_header
    end if

    if (out_old_dns == 'y') open(8,file=file_dns, position = position_option)
    
  endif
  call comm_sync_all

  ! Export to file_trj (initial step)
  if (out_rvf_rt=='y')then
       call Ion_Force_omp(Rion_update_rt,calc_mode_rt)
       write(comment_line,110) -1, 0.0d0
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       &  write(comment_line,112) trim(comment_line), xi_nh
       call write_xyz(comment_line,"new","rvf")
       call write_xyz(comment_line,"add","rvf")
  endif

  Temperature_ion= 0d0
  Ework          = 0d0
  if(use_ehrenfest_md=='y') then
     dt_h       = dt*0.5d0
     Enh_gkTlns = 0d0
     Enh        = 0d0
     if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
        gkT = 3d0*NI * kB/hartree2J*temperature0_ion
        Qnh = gkT * thermostat_tau**2d0
     endif
  endif

!$acc enter data copyin(ik_table,ib_table)
!$acc enter data copyin(lapx,lapy,lapz)
!$acc enter data copyin(nabx,naby,nabz)
!$acc enter data copyin(modx,mody,modz)
!$acc enter data copyin(zJxyz,zKxyz)
!$acc enter data copyin(uV,iuV)
!$acc enter data copyin(kAc)
!$acc enter data copyin(zproj)
!$acc enter data copyin(ik_table,ib_table)
!$acc enter data copyin(ekr_omp)
!$acc enter data copyin(a_tbl,jxyz,mps)

!$acc enter data create(kAc_new)
!$acc enter data create(ghtpsi)

  call timer_begin(LOG_DYNAMICS)
!$acc enter data copyin(zu_t)
  do iter=entrance_iter+1,Nt

    if (use_ehrenfest_md == 'y') then
       !(Velocity Verlet integrator for ion dynamics)
       !NHC act on velocity with dt/2
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call apply_nose_hoover_velocity(dt_h)
       endif
       !update ion velocity with dt/2
       if(iter==0) dRion(:,:,iter-1)= dRion(:,:,iter) - velocity(:,:)*dt
       do ia=1,NI
          mass_au = umass*Mass(Kion(ia))
          velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       enddo
       !velocity scaling
       if(step_velocity_scaling>=1 .and. mod(iter,step_velocity_scaling)==0) then
          call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
          call apply_velocity_scaling_ion(Temperature_ion,velocity)
       endif
       !update ion coordinate with dt
       do ia=1,NI
          mass_au = umass*Mass(Kion(ia))
          dRion(:,ia,iter+1) = dRion(:,ia,iter) + velocity(:,ia)*dt
          Rion(:,ia) = Rion_eq(:,ia) + dRion(:,ia,iter+1)
       enddo
       if (mod(iter,step_update_ps)==0 ) then
          call prep_ps_periodic('update_all       ')
       else if (mod(iter,step_update_ps2)==0 ) then
          call prep_ps_periodic('update_wo_realloc')
       endif

       !NHC act on thermostat with dt
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
          call apply_nose_hoover_thermostat(Temperature_ion,dt)
          Enh_gkTlns = Enh_gkTlns + gkT * xi_nh*dt
          Enh        = Enh_gkTlns + 0.5d0 * Qnh * xi_nh*xi_nh
       endif
    endif

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
    if(alocal_laser=='y') call prep_RT_Ac_alocal_laser(iter+1)

!$acc update device(kAc,kAc_new)
    call current_RT(zu_t)

    javt(iter+1,:)=jav(:)
    if (use_ehrenfest_md == 'y') then
      call current_RT_ion
      javt_ion(iter+1,:)=jav_ion(:)
!$acc update self(zu_t)
      aforce(:,:) = force(:,:)
      call Ion_Force_omp(Rion_update_rt,calc_mode_rt)
      call Total_Energy_omp(Rion_update_rt,calc_mode_rt)
    else
      if(mod(iter,nstep_energy_calc)==0)then
        call Total_Energy_omp(Rion_update_rt,calc_mode_rt)
      end if
    end if

    call timer_begin(LOG_OTHER)


    E_ext(iter,:)=-(Ac_ext(iter+1,:)-Ac_ext(iter-1,:))/(2*dt)
    E_ind(iter,:)=-(Ac_ind(iter+1,:)-Ac_ind(iter-1,:))/(2*dt)
    E_tot(iter,:)=-(Ac_tot(iter+1,:)-Ac_tot(iter-1,:))/(2*dt)

    Eelemag= aLxyz*sum(E_tot(iter,:)**2)/(8.d0*Pi)
    Eall   = Eall + Eelemag
    do ia=1,NI
      FionAc(:,ia)=Zps(Kion(ia))*E_tot(iter,:)
    enddo
    if(alocal_laser=='y') call get_Eelemag_FionAc_alocal_laser(iter)
    force=force+FionAc

    if (use_ehrenfest_md == 'y') then
       aforce(:,:) = 0.5d0*( aforce(:,:) + force(:,:) )

       !update ion velocity with dt/2
       Ework = 0d0
       do ia=1,NI
          mass_au = umass*Mass(Kion(ia))
          velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
          Ework = Ework - sum(aforce(:,ia)*(dRion(:,ia,iter+1)-dRion(:,ia,iter)))
       enddo

       !NHC act on velocity with dt/2
       if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
          call apply_nose_hoover_velocity(dt_h)
       endif

       if (stop_system_momt=='y') call remove_system_momentum(0)
       call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)

    else
      dRion(:,:,iter+1)=0.d0
      Tion=0.d0
    endif
    Eall=Eall+Tion
    
    Eall_t(iter) = Eall
    Tion_t(iter) = Tion
    Temperature_ion_t(iter) = Temperature_ion
    Ework_integ_fdR(iter) = Ework_integ_fdR(iter-1) + Ework

    if(use_ehrenfest_md=='y'.and.ensemble=="NVT".and.thermostat=="nose-hoover")then
       Enh_t(iter)  = Enh
       Hnvt_t(iter) = Eall + Enh
    endif

!---Write section---
    ! Export to file_rt_data
    if (mod(iter, 1000) == 0 .or. iter == Nt) then
      call write_rt_data(iter)
    end if

    ! Export to standard log file
    if (comm_is_root(nproc_id_global)) then
       if (iter/10*10==iter) then
       if (use_ehrenfest_md=='y') then
           write(*,120) iter*dt,&
           & (E_ext(iter,ixyz),E_tot(iter,ixyz),ixyz=1,3),&
           &  Eall, Eall-Eall0, Tion, Temperature_ion
         !! TODO: exclude _t.out file future implementation
         ! write(7,'(1x,100e16.6E3)') iter*dt,&
         !      & (E_ext(iter,ixyz),E_tot(iter,ixyz),ixyz=1,3),&
         !      &  Eall, Eall-Eall0, Tion, Temperature_ion
       else
           write(*,120) iter*dt,&
           & (E_ext(iter,ixyz),E_tot(iter,ixyz),ixyz=1,3),&
           &  Eall, Eall-Eall0, Tion
          ! write(7,'(1x,100e16.6E3)') iter*dt,&
          !      & (E_ext(iter,ixyz),E_tot(iter,ixyz),ixyz=1,3),&
          !      &  Eall, Eall-Eall0, Tion
       endif
       endif
    endif

    ! Export to file_trj
    if (out_rvf_rt=='y' .and. mod(iter,out_rvf_rt_step)==0)then
       if(use_ehrenfest_md=='n') &
       &  call Ion_Force_omp(Rion_update_rt,calc_mode_rt)
       write(comment_line,110) iter, iter*dt
110    format("#rt   step=",i8,"   time",e16.6)
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       &  write(comment_line,112) trim(comment_line), xi_nh
112    format(a,"  xi_nh=",e18.10)
       call write_xyz(comment_line,"add","rvf")
    endif

120 format(1x,f10.4,6f12.6,2e20.10E3,e14.6,f14.5)

    !! Export Dynamical Density (file_dns)
    !! NOTE: The support for _dns.file will be disabled in future release...
    if (out_old_dns == 'y' .and. mod(iter, out_dns_rt_step) == 0) then
      if (comm_is_root(nproc_id_global)) then
        write(8,'(1x,i10)') iter
        do i=1,NL
          write(8,'(1x,2e16.6E3)') rho(i),(rho(i)-rho_gs(i))
        enddo
      endif
      call comm_sync_all()
    end if
    
    ! Export electronic density (cube or vtk)
    if(out_dns_rt=='y' .and. mod(iter,out_dns_rt_step)==0) then
       call write_density(iter,'rt')
    end if

    ! Export transition electronic density at specified energy
    if (out_dns_trans == 'y') call analysis_dns_trans(iter)

    ! Export analysis data(Adiabatic evolution) to file_ovlp,file_nex
    if(projection_option /='no' .and. mod(iter,out_projection_step)==0)then
      call analysis_RT_using_GS(Rion_update_rt,Nscf,zu_t,iter,"projection")
    end if

    !! TODO: Remove the outpit of "j_ac.out" file future...
    ! if(comm_is_root(nproc_id_global))then
    !   if (iter/1000*1000 == iter .or. iter == Nt) then
    !     open(407,file=file_j_ac)
    !     write(407,'(A)')'# J     : Matter current density'
    !     write(407,'(A)')'# Ac_ext: External vector-potential devided by light-velocity'
    !     write(407,'(A)')'# Ac_tot: Total Vector-potential devided by light-velocity'
    !     write(407,'(A1,A25,9A26)')'#','Time ['//trim(t_unit_time%name)//']', &
    !          'Jx ['//trim(t_unit_current%name)//']', &
    !          'Jy ['//trim(t_unit_current%name)//']', &
    !          'Jz ['//trim(t_unit_current%name)//']', &
    !          'Ac_ext_x ['//trim(t_unit_ac%name)//']', &
    !          'Ac_ext_y ['//trim(t_unit_ac%name)//']', &
    !          'Ac_ext_z ['//trim(t_unit_ac%name)//']', &
    !          'Ac_tot_x ['//trim(t_unit_ac%name)//']', &
    !          'Ac_tot_y ['//trim(t_unit_ac%name)//']', &
    !          'Ac_tot_z ['//trim(t_unit_ac%name)//']'
    !     do i=0,Nt
    !       write(407,'(100e26.16E3)') i*dt*t_unit_time%conv, &
    !                                  javt(i,1)*t_unit_current%conv, &
    !                                  javt(i,2)*t_unit_current%conv, &
    !                                  javt(i,3)*t_unit_current%conv, &
    !                                  Ac_ext(i,1)*t_unit_ac%conv, &
    !                                  Ac_ext(i,2)*t_unit_ac%conv, &
    !                                  Ac_ext(i,3)*t_unit_ac%conv, &
    !                                  Ac_tot(i,1)*t_unit_ac%conv, &
    !                                  Ac_tot(i,2)*t_unit_ac%conv, &
    !                                  Ac_tot(i,3)*t_unit_ac%conv
    !     end do
    !     close(407)
    !   end if
    ! end if
    
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
!$acc update self(zu_t)
        call timer_end(LOG_DYNAMICS)
        call prep_restart_write
        go to 1
      end if
    end if

    call timer_end(LOG_OTHER)

    ! backup for system failure
    if (need_backup .and. iter > 0 .and. mod(iter, backup_frequency) == 0) then
      call timer_end(LOG_DYNAMICS)
      call timer_end(LOG_ALL)
      iter_now=iter
      call prep_restart_write
      call timer_begin(LOG_ALL)
      call timer_begin(LOG_DYNAMICS)
    end if
  enddo !end of RT iteraction========================
!$acc exit data copyout(zu_t)
  call timer_end(LOG_DYNAMICS)

  if(comm_is_root(nproc_id_global)) then
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
    ! close(7) !! TODO: Remove output of "_t.out" file future
    if (out_old_dns == 'y') then
      close(8) !! TODO: Disable output of "_dns.out" file for future release
    end if
    close(9)
    if (projection_option /= 'no') then
      close(404)
      close(408)
      close(409)
    end if
  endif

  if (use_ehrenfest_md == 'y') call print_restart_data_md_gs

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'This is the end of RT calculation'
    call timer_show_current_hour('elapse time=',LOG_ALL)
    write(*,*) '-----------------------------------------------------------'
  end if

!====RT calculation===========================
!====Analyzing calculation====================

  ! Export analysis data:Adiabatic evolution to file_last_band_map
  if (projection_option /= 'no') then
    call analysis_RT_using_GS(Rion_update_rt,Nscf,zu_t,iter,"proj_last ")
  end if

  call Fourier_tr

  call comm_sync_all

  if (comm_is_root(nproc_id_global)) write(*,*) 'This is the end of all calculation'
  Time_now=get_wtime()
  call timer_end(LOG_ALL)
  if (comm_is_root(nproc_id_global)) call timer_show_hour('Total time =',LOG_ALL)

1 if(comm_is_root(nproc_id_global)) write(*,*)  'This calculation is shutdown successfully!'

contains

  subroutine reset_rt_timer
    implicit none
    integer :: i
    do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
      call timer_reset(i)
    end do
  end subroutine
  
  
  subroutine write_rt_data(niter)
    implicit none
    integer, intent(in) :: niter
    integer :: fh_rt, fh_rt_energy, iiter

    if (comm_is_root(nproc_id_global)) then
      fh_rt = open_filehandle(file_rt_data)
      
      write(fh_rt, '("#",1X,A)') "Real time calculation"
      
      write(fh_rt, '("#",1X,A,":",1X,A)') "Ac_ext", "External vector potential field"
      write(fh_rt, '("#",1X,A,":",1X,A)') "E_ext", "External electric field"
      write(fh_rt, '("#",1X,A,":",1X,A)') "Ac_tot", "Total vector potential field"
      write(fh_rt, '("#",1X,A,":",1X,A)') "E_tot", "Total electric field"
      if(use_ehrenfest_md=='y') then
        write(fh_rt, '("#",1X,A,":",1X,A)') "Jm", "Matter current density(electrons)"
        write(fh_rt, '("#",1X,A,":",1X,A)') "Jmi","Matter current density(ions)"
      else
        write(fh_rt, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
      endif
      write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Ac_ext_x", trim(t_unit_ac%name), &
        & 3, "Ac_ext_y", trim(t_unit_ac%name), &
        & 4, "Ac_ext_z", trim(t_unit_ac%name), &
        & 5, "E_ext_x", trim(t_unit_elec%name), &
        & 6, "E_ext_y", trim(t_unit_elec%name), &
        & 7, "E_ext_z", trim(t_unit_elec%name), &
        & 8, "Ac_tot_x", trim(t_unit_ac%name), &
        & 9, "Ac_tot_y", trim(t_unit_ac%name), &
        & 10, "Ac_tot_z", trim(t_unit_ac%name), &
        & 11, "E_tot_x", trim(t_unit_elec%name), &
        & 12, "E_tot_y", trim(t_unit_elec%name), &
        & 13, "E_tot_z", trim(t_unit_elec%name), &
        & 14, "Jm_x", trim(t_unit_current%name), &
        & 15, "Jm_y", trim(t_unit_current%name), &
        & 16, "Jm_z", trim(t_unit_current%name)
      if(use_ehrenfest_md=='y') then
      write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 17, "Jmi_x", trim(t_unit_current%name), &
        & 18, "Jmi_y", trim(t_unit_current%name), &
        & 19, "Jmi_z", trim(t_unit_current%name)
      endif

      write(fh_rt,*)

      do iiter = 0, niter
        write(fh_rt, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & iiter * dt * t_unit_time%conv, &
          & Ac_ext(iiter, 1:3) * t_unit_ac%conv, &
          & E_ext(iiter, 1:3) * t_unit_elec%conv, &
          & Ac_tot(iiter, 1:3) * t_unit_ac%conv, &
          & E_tot(iiter, 1:3) * t_unit_elec%conv, &
          & javt(iiter, 1:3) * t_unit_current%conv
        if(use_ehrenfest_md=='y') then
        write(fh_rt, "(99(1X,E23.15E3))",advance='no') &
          & javt_ion(iiter, 1:3) * t_unit_current%conv
        endif
        write(fh_rt,*)
      end do
      close(fh_rt)

      fh_rt_energy = open_filehandle(file_rt_energy_data)

      write(fh_rt_energy, '("#",1X,A)') "Real time calculation"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Eall", "Total energy"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Eall0", "Initial energy"
      if(use_ehrenfest_md=='y') then
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Tion", "Kinetic energy of ions"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Temperature_ion", "Temperature of ions"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "E_work", "Work energy of ions(sum f*dr)"
      if(ensemble=="NVT".and.thermostat=="nose-hoover")then
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Enh", "NH thermostat energy (MD)"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat(MD)"
      write(fh_rt_energy, '("#",1X,A,":",1X,A)') "Hnvt'","Hnvt using E_work"
      endif
      endif

      write(fh_rt_energy, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Eall", trim(t_unit_energy%name), &
        & 3, "Eall-Eall0", trim(t_unit_energy%name)

      if(use_ehrenfest_md=='y') then
      write(fh_rt_energy, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 4, "Tion", trim(t_unit_energy%name), &
        & 5, "Temperature_ion", "K", &
        & 6, "E_work", trim(t_unit_energy%name)
      if(ensemble=="NVT".and.thermostat=="nose-hoover")then
      write(fh_rt_energy, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 7, "Enh",  trim(t_unit_energy%name), &
        & 8, "Hnvt", trim(t_unit_energy%name), &
        & 9, "Hnvt'",trim(t_unit_energy%name)
      endif
      endif

      write(fh_rt_energy,*)

      do iiter = 0, niter
        if( use_ehrenfest_md/='y' .and. mod(iiter,nstep_energy_calc)/=0)cycle
        write(fh_rt_energy, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & iiter * dt * t_unit_time%conv, &
          & Eall_t(iiter) * t_unit_energy%conv, &
          & (Eall_t(iiter) - Eall0) * t_unit_energy%conv
        if(use_ehrenfest_md=='y') then
        write(fh_rt_energy, "(99(1X,E23.15E3))",advance='no') &
          & Tion_t(iiter) * t_unit_energy%conv, &
          & Temperature_ion_t(iiter), &
          & Ework_integ_fdR(iiter) * t_unit_energy%conv
        if(ensemble=="NVT".and.thermostat=="nose-hoover")then
        write(fh_rt_energy, "(99(1X,E23.15E3))",advance='no') &
          & Enh_t(iiter) * t_unit_energy%conv, &
          & Hnvt_t(iiter) * t_unit_energy%conv, &
          & (Tion_t(iiter)+Ework_integ_fdR(iiter)+Enh_t(iiter)) * t_unit_energy%conv
        endif
        endif
        write(fh_rt_energy,*)
      end do


      close(fh_rt_energy)


    end if
    call comm_sync_all
    return
  end subroutine write_rt_data
  
end subroutine tddft_sc

  subroutine write_xyz(comment,action,rvf)
  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use Global_Variables
    use inputoutput, only: au_length_aa
    use salmon_global, only: SYSname,iflag_atom_coor,ntype_atom_coor_cartesian,ntype_atom_coor_reduced
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none
    integer :: ia,unit_xyz=200
    character(3) :: action,rvf
    character(1024) :: file_trj
    character(*) :: comment

    if(.not. comm_is_root(nproc_id_global)) return

    if(action=='new') then

       file_trj=trim(SYSname)//'_trj.xyz'
       open(unit_xyz,file=trim(file_trj),status="unknown")

    else if(action=='add') then

       write(unit_xyz,*) NI
       write(unit_xyz,*) trim(comment)
       do ia=1,NI
          if(      rvf=="r  " ) then
             write(unit_xyz,100) trim(atom_name(ia)),Rion(1:3,ia)*au_length_aa
          else if( rvf=="rv " ) then
             write(unit_xyz,110) trim(atom_name(ia)),Rion(1:3,ia)*au_length_aa,velocity(1:3,ia)
          else if( rvf=="rvf" ) then
             write(unit_xyz,120) trim(atom_name(ia)),Rion(1:3,ia)*au_length_aa,velocity(1:3,ia),force(1:3,ia)
          endif
       enddo

    else if(action=='end') then
       close(unit_xyz)
    endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

  end subroutine write_xyz

end module control_sc
