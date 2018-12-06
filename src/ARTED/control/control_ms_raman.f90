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
!This file is "control_ms.f90"
!This file contains ms-mode program
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module control_ms_raman
  implicit none
contains
subroutine raman_maxwell_ms
  use Global_Variables
  use timer
  use opt_variables
  use performance_analyzer
  use salmon_parallel
  use salmon_communication
  use salmon_file
  use misc_routines
  use inputoutput, only: t_unit_time, t_unit_current, t_unit_ac
  use restart, only: prep_restart_write

  implicit none
  integer :: iter
  integer :: ix_m, iy_m, iz_m
  integer :: imacro,igroup,i
  integer :: index
  logical :: flag_shutdown = .false.
  logical :: flg_out_ms_step, flg_out_ms_next_step
  logical :: flg_out_projection_step, flg_out_projection_next_step
  real(8) :: aforce(3,NI),Temperature_ion
 !real(8) :: Enh, Enh_gkTlns, gkT, Qnh   !NVT
  character(100) :: comment_line

#ifdef ARTED_LBLK
  call opt_vars_init_t4ppt()
#endif

  if (restart_option == 'restart') then
    position_option='asis'
  else if(restart_option == 'new')then

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
    if (trim(FDTDdim) == '2DC') then
      ! TODO: FIx the initialization routine
      call init_Ac_ms_2dc()
    else
      call init_Ac_ms
    endif

!xx    !(for restarting with read_rt_wfn_k_ms=y option)
!xx    if(read_rt_wfn_k_ms=='y') then
!xx       call add_init_Ac_ms
!xx       call assign_mp_variables_omp()
!xx       do imacro = nmacro_s, nmacro_e
!xx          call psi_rho_RT(zu_t)
!xx          call Hartree
!xx          call Exc_Cor(calc_mode_rt,NBoccmax,zu_t)
!xx          Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
!xx          call Total_Energy_omp(rion_update_on,calc_mode_rt,imacro)
!xx          call Ion_Force_omp(rion_update_on,calc_mode_rt,imacro)
!xx       enddo
!xx       do imacro = 1, nmacro
!xx          ix_m = macropoint(1,imacro)
!xx          iy_m = macropoint(2,imacro)
!xx          iz_m = macropoint(3,imacro)
!xx          Jm_new_m(1:3,imacro)=matmul(trans_mat(1:3,1:3),Jm_new_ms(1:3,ix_m,iy_m,iz_m))
!xx       end do
!xx    endif

!xx    rho_gs(:)=rho(:)

!xx    Vloc_old(:,1) = Vloc(:); Vloc_old(:,2) = Vloc(:)
    
!xx    do imacro = nmacro_s, nmacro_e
!xx      zu_m(:, :, :, imacro) =zu_t(:, 1:NBoccmax, :)
!xx    end do
    
!xx    if (NXYsplit /= 1) then
!xx      do imacro = nmacro_s, nmacro_e
!xx        rho_m(:,imacro)=rho(:)
!xx        Vh_m(:,imacro)=Vh(:)
!xx        Vexc_m(:,imacro)=Vexc(:)
!xx        Eexc_m(:,imacro)=Eexc(:)
!xx        Vloc_m(:,imacro)=Vloc(:)
!xx        Vloc_old_m(:,:,imacro)=Vloc_old(:,:)
!xx      end do
!xx    end if
    
    deallocate(zu_t)

!reentrance

    position_option='rewind'
    entrance_iter=-1
    call reset_rt_timer
  end if

  ! create directory for exporting
  call create_dir_ms

  ! Export to file_trj (initial step)
  if (out_rvf_rt=='y')then
      write(comment_line,110) -1, 0.0d0
     !if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
     !&  write(comment_line,112) trim(comment_line), xi_nh

      do imacro = nmacro_s, nmacro_e
        !call Ion_Force_omp(Rion_update_rt,calc_mode_rt,imacro)
         call Ion_Force_ForceField_MS(imacro)
         call Ion_Force_RamanTensor_MS(imacro)
      enddo

      do imacro = nmacro_s, nmacro_e
         if(flag_ms_ff_LessPrint)then
            if (imacro==1)then
            call write_xyz_ms(comment_line,"new","rvf",imacro)
            call write_xyz_ms(comment_line,"add","rvf",imacro)
            endif
         else
            call write_xyz_ms(comment_line,"new","rvf",imacro)
            call write_xyz_ms(comment_line,"add","rvf",imacro)
         endif
      enddo
      call comm_sync_all
  endif

  !(get ion current for initial step)
  if(use_ehrenfest_md=='y') then 
     do imacro = nmacro_s, nmacro_e
        call current_RT_ion_MS(imacro)
        if (comm_is_root(nproc_id_tdks)) then
            jm_ion_new_m_tmp(1:3, imacro) = jav_ion(1:3)
        end if
     enddo
     call comm_summation(jm_ion_new_m_tmp,Jm_ion_m,3*nmacro,nproc_group_global)
     do imacro = 1, nmacro
        ix_m = macropoint(1,imacro)
        iy_m = macropoint(2,imacro)
        iz_m = macropoint(3,imacro)
        Jm_ion_ms(1:3,ix_m,iy_m,iz_m) = matmul(trans_inv(1:3,1:3),Jm_ion_m(1:3,imacro))
     enddo
  endif

  !if(use_ehrenfest_md=='y') then
    !Enh_gkTlns = 0d0
    !Enh        = 0d0
    !if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
    !   gkT = 3d0*NI * kB/hartree2J*temperature0_ion
    !   Qnh = gkT * thermostat_tau**2d0
    !endif
  !endif
 
!xx  ! Output filename
!xx  write(file_energy_transfer, "(A,'energy-transfer.data')") trim(directory)
  
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
!$acc enter data copyin(zu_m)
  RTiteratopm : do iter=entrance_iter+1, Nt ! sato

    !! NOTE: flg_out_ms_step (the macroscopic field will exported in this step)
    flg_out_ms_step = .false.
    if (mod(iter, out_ms_step) == 0) then
      flg_out_ms_step = .true.
    end if

    !! NOTE: flg_out_projection_step (the projection will exported in this step)
    flg_out_projection_step = .false.
    if (projection_option /= 'no') then
      if ((mod(iter, out_projection_step) == 0) .or. (iter == Nt)) then 
        flg_out_projection_step = .true.
      end if
    end if

    !! Update of the Macroscopic System
    
    !===========================================================================
    call timer_begin(LOG_OTHER)
    !! NOTE: Update the macroscopic variables:
    !!       Ac_old_ms = Ac_ms; Ac_ms = Ac_new_ms; Ac_new_ms = 0
    !!       Jm_old_ms = Jm_ms; Jm_ms = Jm_new_ms; Jm_new_ms = 0
    !!       Energy_Elec_ms = map(Energy_Elec_Matter_new_m); Jm_m = Jm_new_m;
    call proceed_ms_variables_omp()
    call timer_end(LOG_OTHER)
    !===========================================================================

    !===========================================================================
    !! Calculate "iter+1" macroscopic field (Ac_new_ms, Jm_new_ms)
    call timer_begin(LOG_DT_EVOLVE_AC)
   !call dt_evolve_Ac() ! Timer: LOG_DT_EVOLVE_AC
    iter_save=iter
    call dt_evolve_Ac_1d_raman() ! Timer: LOG_DT_EVOLVE_AC
    call timer_end(LOG_DT_EVOLVE_AC)
    !===========================================================================
    !===========================================================================    
    !! Compute EM field, energy and other important quantities...
    call timer_begin(LOG_OTHER) ! TODO: Create new timer "LOG_EMFIELD" later
    call calc_energy_joule()
    if (flg_out_ms_step) then !! mod(iter, out_ms_step) == 0
      call calc_elec_field()
      call calc_bmag_field()
      call calc_energy_elemag()
      call calc_total_energy
    end if
    call timer_end(LOG_OTHER)
    !===========================================================================    


    !=========================================================================== 
    call timer_begin(LOG_OTHER)
    !! NOTE: Mapping between the macropoint and the gridsystem
    !!       Ac_m <- Ac_ms; Ac_new_m <- Ac_new_ms
    !!       data_local_Ac_m = Ac_m; data_local_Jm_m = Jm_m
    call assign_mp_variables_omp()
    !! Store data_vac_ac variables
    call store_data_local_ac_jm()
    call store_data_vac_ac()
    !! Store data_out variables
    if (flg_out_ms_step) then !! mod(iter, out_ms_step) == 0
      index = iter / out_ms_step
      call store_data_out_omp(index)
    end if
    call timer_end(LOG_OTHER)
    !===========================================================================

    !===========================================================================
    call timer_begin(LOG_OTHER)
    if (flg_out_ms_step .and. comm_is_root(nproc_id_global)) then
      call trace_ms_calculation()
    end if
    call timer_end(LOG_OTHER)
    !===========================================================================

    !! Update of the Microscopic System
    
    ! if (iter == Nt) break

    !! NOTE: Since the microscopic calculation at final step doesn't written in
    !!       any output file by the program. The calculation of the microscopic
    !!       system at iter=Nt is seemingfly meaningless.
    !!       However, in present implementation, the iter=Nt data is stored to 
    !!       the BACKUP, this checkpoint helps to restart from Nt+1 step.
    !!       If the backup is not required in the calculation, 
    !!       the use of the above "break" is prefereble... 
    !! TODO: Make it possible to export all data from 0 to Nt+1.
    
    !! NOTE: flg_out_ms_step (the macroscopic field will exported in next step)
    flg_out_ms_next_step = .false.
    if (mod(iter+1, out_ms_step) == 0) then
      flg_out_ms_next_step = .true.
    end if

    !! NOTE: flg_out_projection_step (the projection will exported in next step)
    flg_out_projection_next_step = .false.
    if (projection_option /= 'no') then
      if ((mod(iter+1, out_projection_step) == 0) .or. (iter+1 == Nt)) then 
        flg_out_projection_next_step = .true.
      end if
    end if
    
    Macro_loop : do imacro = nmacro_s, nmacro_e
      
      !===========================================================================
      !! Extract microscopic state of "imacro"-th macropoint
      call timer_begin(LOG_OTHER)
      if (NXYsplit /= 1) call get_macro_data(imacro)
      call timer_end(LOG_OTHER)
      !===========================================================================

      !===========================================================================
      !! update ion coordinate and velocity in MD option (part-1)
      if (use_ehrenfest_md == 'y') then
         call dt_evolve_MD_1_MS(iter,imacro)
      endif
      !===========================================================================

      !===========================================================================
!xx      call dt_evolve_KB_MS(imacro) ! Timer: LOG_DT_EVOLVE
      !===========================================================================

      !===========================================================================
      !! Store microscopic state of "imacro"-th macropoint
      call timer_begin(LOG_OTHER)
      if(NXYsplit /= 1)then
        call put_macro_data(imacro)
      end if
      call timer_end(LOG_OTHER)
      !===========================================================================

!xx      kAc(:,1) = kAc0(:,1) + Ac_new_m(1,imacro)
!xx      kAc(:,2) = kAc0(:,2) + Ac_new_m(2,imacro)
!xx      kAc(:,3) = kAc0(:,3) + Ac_new_m(3,imacro)

      !===========================================================================
!xx!$acc update device(kAc,kAc_new)
!xx      call current_RT_MS(imacro) ! Timer: LOG_CURRENT
!xx      if (use_ehrenfest_md == 'y') then
!xx         call current_RT_ion_MS(imacro)
!xx      endif
      !===========================================================================

      !===========================================================================
      call timer_begin(LOG_OTHER)
      if(Sym /= 1)then
        jav(1:2) = 0d0
      end if
      
      !===========================================================================
      !! Special rule for debug mode (macroscopic system does not have a matter)
      !! NOTE: This condition will be removed and replaced by FDTD mode
      !!       after merging the common Maxwell calculation routine.
      if (debug_switch_no_radiation) then
        jav(:) = 0d0
      end if
      !===========================================================================
      
!xx      if (comm_is_root(nproc_id_tdks)) then
!xx        jm_new_m_tmp(1:3, imacro) = jav(1:3)
!xx      end if
!xx      javt(iter+1,:) = jav(:)
!xx      if (use_ehrenfest_md == 'y') then
!xx         if (comm_is_root(nproc_id_tdks)) then
!xx            jm_ion_new_m_tmp(1:3, imacro) = jav_ion(1:3)
!xx         end if
!xx      end if
      call timer_end(LOG_OTHER)
      !===========================================================================

      if (use_ehrenfest_md == 'y') then

       aforce(:,:) = force_m(:,:,imacro)
       call Ion_Force_ForceField_MS(imacro)
       call Ion_Force_RamanTensor_MS(imacro)

!xx!$acc update self(zu_m(:,:,:,imacro))
!xx        aforce(:,:) = force_m(:,:,imacro)
!xx        call Ion_Force_omp(Rion_update_rt,calc_mode_rt,imacro)
!xx        call Ion_Force_Ac_MS(imacro)
!xx        if (flg_out_ms_next_step) then
!xx          call Total_Energy_omp(Rion_update_rt,calc_mode_rt,imacro)
!xx        end if
!xx      else
!xx        if (flg_out_ms_next_step) then
!xx          call Total_Energy_omp(Rion_update_rt,calc_mode_rt,imacro)
!xx!$acc update self(zu_m(:,:,:,imacro))
!xx          call Ion_Force_omp(Rion_update_rt,calc_mode_rt,imacro)
!xx        end if
      end if
    
      !===========================================================================
      ! Calculate + store electron energy (if required in the next iteration..)
      call timer_begin(LOG_OTHER)
      if (flg_out_ms_next_step) then !! mod(iter+1, out_ms_step) == 0
        if(comm_is_root(nproc_id_tdks)) then ! sato
          energy_elec_Matter_new_m_tmp(imacro) = Eall - Eall0_m(imacro)
        end if
      end if ! sato
      call timer_end(LOG_OTHER)
      !===========================================================================

      !===========================================================================
      !! update ion coordinate and velocity in MD option (part-2)
      if (use_ehrenfest_md == 'y') then
         call dt_evolve_MD_2_MS(aforce,iter,imacro)
      endif
      !(for exporting to file_trj later)
!xx      if (out_rvf_rt=='y' .and. mod(iter,out_rvf_rt_step)==0)then
!xx         if(use_ehrenfest_md=='n') &
!xx         &  call Ion_Force_omp(Rion_update_rt,calc_mode_rt,imacro)
!xx      endif
      !===========================================================================

      !===========================================================================
      ! Calculate + store excitation number (if required in the next iteration..)
!xx      call timer_begin(LOG_ANA_RT_USEGS)
!xx      if (flg_out_projection_next_step) then
!xx!$acc update self(zu_m(:,:,:,imacro))
!xx        call analysis_RT_using_GS(Rion_update_rt,Nscf,zu_m(:,:,:,imacro),iter,"projection")
!xx        if(comm_is_root(nproc_id_tdks))then ! sato
!xx          excited_electron_new_m_tmp(imacro) = sum(occ) - sum(ovlp_occ(1:NBoccmax,:))
!xx        end if ! sato
!xx      end if
!xx      call timer_end(LOG_ANA_RT_USEGS)
    
    end do Macro_loop !end of Macro_loop iteraction========================

    !===========================================================================
    call timer_begin(LOG_ALLREDUCE)
!xx    call comm_summation(jm_new_m_tmp, Jm_new_m, 3 * nmacro, nproc_group_global)
!xx    if(use_ehrenfest_md=='y') &
!xx    & call comm_summation(jm_ion_new_m_tmp, Jm_ion_new_m, 3*nmacro,nproc_group_global)
    if (flg_out_ms_next_step) then
      call comm_summation(energy_elec_Matter_new_m_tmp, energy_elec_Matter_new_m, nmacro, nproc_group_global)
    end if
    if (flg_out_projection_next_step) then
      call comm_summation(excited_electron_new_m_tmp, excited_electron_new_m, nmacro, nproc_group_global)
      if(comm_is_root(nproc_id_global)) call write_excited_electron(iter+1)
    end if
    call timer_end(LOG_ALLREDUCE)
    !===========================================================================

    
    !===========================================================================
!xx    call timer_begin(LOG_OTHER)
!xx!$omp parallel do default(shared) private(imacro, ix_m, iy_m, iz_m)
!xx    do imacro = 1, nmacro
!xx      !! NOTE: If the array "macropoint" is appropriately setted,
!xx      !!       every set of (ix_m, iy_m, iz_m) would be independent. Therefore,
!xx      !!       the OpenMP directives of 'reduction' is not required as above.
!xx      ix_m = macropoint(1, imacro)
!xx      iy_m = macropoint(2, imacro)
!xx      iz_m = macropoint(3, imacro)
!xx      !! Map the local macropoint current into the jm field
!xx      !Jm_new_ms(1:3, ix_m, iy_m, iz_m) = Jm_new_ms(1:3, ix_m, iy_m, iz_m) & 
!xx      !                               & + Jm_new_m(1:3, imacro)
!xx      Jm_new_ms(1:3, ix_m, iy_m, iz_m) = matmul(trans_inv(1:3,1:3), Jm_new_m(1:3, imacro))
!xx      if(use_ehrenfest_md=='y') &
!xx      Jm_ion_new_ms(1:3,ix_m,iy_m,iz_m) = matmul(trans_inv(1:3,1:3), Jm_ion_new_m(1:3,imacro))
!xx    end do
!xx!$omp end parallel do
!xx    call timer_end(LOG_OTHER)
    !===========================================================================

    ! Export to file_trj
    if (out_rvf_rt=='y' .and. mod(iter,out_rvf_rt_step)==0)then
       write(comment_line,110) iter, iter*dt
110    format("#rt   step=",i8,"   time",e16.6)
       !if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       !&  write(comment_line,112) trim(comment_line), xi_nh
       !112      format(a,"  xi_nh=",e18.10)
       if(flag_ms_ff_LessPrint) then  !AY
          do imacro = nmacro_s, nmacro_e
             if(imacro==1) call write_xyz_ms(comment_line,"add","rvf",imacro)
          enddo
          call comm_sync_all
       else   !AY
          do igroup=1,ndivide_macro
             do i=1,nmacro_write_group
                imacro = (igroup-1)*nmacro_write_group + i
                if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e) &
                     & call write_xyz_ms(comment_line,"add","rvf",imacro)
             enddo
             call comm_sync_all
          enddo
       endif
    endif
    !===========================================================================
    
    ! Shutdown sequence
    call timer_begin(LOG_OTHER)
    if(mod(iter,10) == 1) then
      call comm_bcast(flag_shutdown,nproc_group_global)
    end if
    if (flag_shutdown) then 
      call comm_sync_all
      write(*,*) nproc_id_global,'iter =',iter
      iter_now=iter
!$acc update self(zu_m)
      call timer_end(LOG_DYNAMICS)
      call prep_restart_write
      go to 1
    end if
    call timer_end(LOG_OTHER)
    
    ! Timer for shutdown
    if (mod(iter,10) == 0) then
      Time_now=get_wtime()
      if (comm_is_root(nproc_id_global) .and. iter/100*100 == iter) then
        write(*,*) 'Total time =',(Time_now-Time_start)
      end if
      if ((Time_now - Time_start)>Time_shutdown .and. Time_shutdown >= 0d0) then 
        flag_shutdown =.true.
      end if
    end if

    ! backup for system failure
    if (need_backup .and. iter > 0 .and. mod(iter, backup_frequency) == 0) then
      if (comm_is_root(nproc_id_global)) call timer_show_current_hour('Backup...', LOG_ALL)
      call timer_end(LOG_DYNAMICS)
      call timer_end(LOG_ALL)
      iter_now=iter
      call prep_restart_write
      call timer_begin(LOG_ALL)
      call timer_begin(LOG_DYNAMICS)
    end if
    
  enddo RTiteratopm !end of RT iteraction========================
!$acc exit data copyout(zu_m)
  call timer_end(LOG_DYNAMICS)

  if(comm_is_root(nproc_id_global)) then
    write(*,'(1X,A)') '-----------------------------------------------'

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
    call timer_show_min ('ana_RT_useGS time  :', LOG_ANA_RT_USEGS)
    call timer_show_min ('Other time         :', LOG_OTHER)
    call timer_show_min ('Allreduce time     :', LOG_ALLREDUCE)
  end if
  call write_performance(trim(directory)//'ms_performance')

  if(comm_is_root(nproc_id_global)) write(*,*) 'This is the start of write section'
  call timer_begin(LOG_IO)

  ! Write data out by using MPI
  if(.not. flag_ms_ff_LessPrint) then !AY
!XXXXX Currently do not write this one to reduce calculation time
!XXXXX If you want to write down, remove the following the commented out.
!     do index = 0, Ndata_out-1
!        call write_data_out(index)
!     end do
  endif  
  call write_data_local_ac_jm()
  call write_data_vac_ac()

  ! Export last atomic coordinate and velocity & Close file_trj
  if (use_ehrenfest_md=='y')then
     if(flag_ms_ff_LessPrint) then !AY
        do imacro = nmacro_s, nmacro_e
           if(imacro==1) then
              call write_xyz_ms(comment_line,"end","rvf",imacro)
              call write_ini_coor_vel_ms_for_restart(imacro)
           endif
        enddo
        call comm_sync_all
     else  !AY
        do igroup=1,ndivide_macro
           do i=1,nmacro_write_group
              imacro = (igroup-1)*nmacro_write_group + i
              if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e) then
                 call write_xyz_ms(comment_line,"end","rvf",imacro)
                 call write_ini_coor_vel_ms_for_restart(imacro)
              endif
           enddo
           call comm_sync_all
        enddo
     endif
  endif


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
  subroutine get_macro_data(jmacro)
    implicit none
    integer, intent(in) :: jmacro
    integer :: il
!$omp parallel do &
!$omp     default(none) private(il) firstprivate(jmacro) &
!$omp     shared(NL,Vh,Vh_m,Vexc,Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m,Vloc_old,Vloc_old_m)
    do il=1,NL
      Vh(il)         = Vh_m(il,jmacro)
      Vexc(il)       = Vexc_m(il,jmacro)
      Eexc(il)       = Eexc_m(il,jmacro)
      Vloc(il)       = Vloc_m(il,jmacro)
      Vloc_old(il,:) = Vloc_old_m(il,:,jmacro)
    end do
!$omp end parallel do
  end subroutine

  subroutine put_macro_data(jmacro)
    implicit none
    integer, intent(in) :: jmacro
    integer :: il
!$omp parallel do &
!$omp     default(none) private(il) firstprivate(jmacro) &
!$omp     shared(NL,Vh,Vh_m,Vexc,Vexc_m,Eexc,Eexc_m,Vloc,Vloc_m)
    do il=1,NL
      Vh_m(il,jmacro)   = Vh(il)
      Vexc_m(il,jmacro) = Vexc(il)
      Eexc_m(il,jmacro) = Eexc(il)
      Vloc_m(il,jmacro) = Vloc(il)
    end do
!$omp end parallel do
  end subroutine


  subroutine reset_rt_timer
    implicit none
    integer :: i
    do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
      call timer_reset(i)
    end do
  end subroutine
  
  
  subroutine proceed_ms_variables_omp()
    implicit none
    integer :: iimacro, iix_m, iiy_m, iiz_m
!$omp parallel do collapse(3) default(shared) private(iix_m, iiy_m, iiz_m)
    do iiz_m = mz1_m, mz2_m
      do iiy_m = my1_m, my2_m
        do iix_m = mx1_m, mx2_m
          Ac_old_ms(1:3, iix_m, iiy_m, iiz_m) = Ac_ms    (1:3, iix_m, iiy_m, iiz_m)
          Ac_ms    (1:3, iix_m, iiy_m, iiz_m) = Ac_new_ms(1:3, iix_m, iiy_m, iiz_m)
          Ac_new_ms(1:3, iix_m, iiy_m, iiz_m) = 0d0
        end do
      end do
    end do
!$omp end parallel do

!$omp parallel do collapse(3) default(shared) private(iix_m, iiy_m, iiz_m)
    do iiz_m = nz1_m, nz2_m
      do iiy_m = ny1_m, ny2_m
        do iix_m = nx1_m, nx2_m
          Jm_old_ms(1:3, iix_m, iiy_m, iiz_m) = Jm_ms    (1:3, iix_m, iiy_m, iiz_m)
          Jm_ms    (1:3, iix_m, iiy_m, iiz_m) = Jm_new_ms(1:3, iix_m, iiy_m, iiz_m)
          Jm_new_ms(1:3, iix_m, iiy_m, iiz_m) = 0d0
        end do
      end do
    end do
!$omp end parallel do

!$omp parallel do default(shared) private(iimacro, iix_m, iiy_m, iiz_m)
    do iimacro = 1, nmacro
      Jm_m(1:3, iimacro) = Jm_new_m(1:3, iimacro)
      if (flg_out_ms_step) then
        iix_m = macropoint(1, iimacro)
        iiy_m = macropoint(2, iimacro)
        iiz_m = macropoint(3, iimacro)
        energy_elec_ms(iix_m, iiy_m, iiz_m) = energy_elec_Matter_new_m(iimacro)
      end if
    end do
!$omp end parallel do

  if(use_ehrenfest_md=='y')then
!$omp parallel do collapse(3) default(shared) private(iix_m, iiy_m, iiz_m)
    do iiz_m = nz1_m, nz2_m
      do iiy_m = ny1_m, ny2_m
        do iix_m = nx1_m, nx2_m
          Jm_ion_old_ms(1:3,iix_m,iiy_m,iiz_m) = Jm_ion_ms    (1:3,iix_m,iiy_m,iiz_m)
          Jm_ion_ms    (1:3,iix_m,iiy_m,iiz_m) = Jm_ion_new_ms(1:3,iix_m,iiy_m,iiz_m)
          Jm_ion_new_ms(1:3,iix_m,iiy_m,iiz_m) = 0d0
        end do
      end do
    end do
!$omp end parallel do

!$omp parallel do default(shared) private(iimacro)
    do iimacro = 1, nmacro
      Jm_ion_m(1:3, iimacro) = Jm_ion_new_m(1:3, iimacro)
    end do
!$omp end parallel do
  endif

  end subroutine proceed_ms_variables_omp


  subroutine assign_mp_variables_omp()
    implicit none
    integer :: iix_m, iiy_m, iiz_m, iimacro
!$omp parallel do default(shared) private(iimacro, iix_m, iiy_m, iiz_m)
    do iimacro = 1, nmacro
      iix_m = macropoint(1, iimacro)
      iiy_m = macropoint(2, iimacro)
      iiz_m = macropoint(3, iimacro)
      !! Assign the vector potential into the local macropoint variables
      !Ac_m(1:3, iimacro) = Ac_ms(1:3, iix_m, iiy_m, iiz_m)
      !Ac_new_m(1:3, iimacro) = Ac_new_ms(1:3, iix_m, iiy_m, iiz_m)
      !Ac_old_m(1:3, iimacro) = Ac_old_ms(1:3, iix_m, iiy_m, iiz_m)
      Ac_m(1:3, iimacro)     = matmul(trans_mat(1:3,1:3),Ac_ms(    1:3, iix_m, iiy_m, iiz_m))
      Ac_new_m(1:3, iimacro) = matmul(trans_mat(1:3,1:3),Ac_new_ms(1:3, iix_m, iiy_m, iiz_m))
      Ac_old_m(1:3, iimacro) = matmul(trans_mat(1:3,1:3),Ac_old_ms(1:3, iix_m, iiy_m, iiz_m))
    end do
!$omp end parallel do
  end subroutine assign_mp_variables_omp
  
  
  subroutine store_data_local_ac_jm()
    implicit none
    integer :: iimacro
!$omp parallel do default(shared) private(iimacro)
    do iimacro = nmacro_s, nmacro_e
      !! Store data_local_Ac, data_local_Jm
      data_local_Ac(1:3, iimacro, iter) = Ac_m(1:3, iimacro)
      data_local_jm(1:3, iimacro, iter) = Jm_m(1:3, iimacro)
      if(use_ehrenfest_md=='y') &
      &  data_local_jm_ion(1:3,iimacro,iter) = Jm_ion_m(1:3,iimacro)
    end do
!$omp end parallel do
  end subroutine store_data_local_ac_jm
  

  subroutine store_data_out_omp(index)
    use Global_Variables
    implicit none
    integer, intent(in) :: index
    integer :: iproc, ipos
    integer :: iix_m, iiy_m, iiz_m
    
    iproc = mod(index, nproc_size_global)
    ipos = (index - iproc) / nproc_size_global
    if (iproc == nproc_id_global) then
!$omp parallel do collapse(3) default(shared) private(iiy_m, iix_m, iiz_m)
      do iiz_m = nz1_m, nz2_m
        do iiy_m = ny1_m, ny2_m
          do iix_m = nx1_m, nx2_m
            data_out(1:3, iix_m, iiy_m, iiz_m, ipos) = Ac_ms(1:3, iix_m, iiy_m, iiz_m)
            data_out(4:6, iix_m, iiy_m, iiz_m, ipos) = Elec_ms(1:3, iix_m, iiy_m, iiz_m)
            data_out(7:9, iix_m, iiy_m, iiz_m, ipos) = Bmag_ms(1:3, iix_m, iiy_m, iiz_m)
            data_out(10:12, iix_m, iiy_m, iiz_m, ipos) = Jm_ms(1:3, iix_m, iiy_m, iiz_m)
            data_out(13, iix_m, iiy_m, iiz_m, ipos) = Energy_elec_ms(iix_m, iiy_m, iiz_m)
            data_out(14, iix_m, iiy_m, iiz_m, ipos) = Energy_joule_ms(iix_m, iiy_m, iiz_m)
            data_out(15, iix_m, iiy_m, iiz_m, ipos) = Energy_elemag_ms(iix_m, iiy_m, iiz_m)
          end do
        end do
      end do
!$omp end parallel do
    end if
  end subroutine
  
  
  subroutine create_dir_ms
    implicit none
    integer :: imacro

    call create_directory(dir_ms)
    call create_directory(dir_ms_RT_Ac)

    do imacro = nmacro_s, nmacro_e
       call create_directory(dir_ms_M(imacro))
    enddo

  end subroutine

  subroutine write_data_out(index)
    implicit none
    integer, intent(in) :: index
    integer :: iix_m, iiy_m, iiz_m, fh_ac
    integer :: iproc, ipos
    
    iproc = mod(index, nproc_size_global)
    ipos = (index - iproc) / nproc_size_global
    if (iproc == nproc_id_global) then
      write(file_ac, "(A,A,'_Ac_',I6.6,'.data')") & 
        & trim(dir_ms_RT_Ac), trim(SYSname),  out_ms_step * index
      fh_ac = open_filehandle(file_ac)
      write(fh_ac, '("#",1X,A)') "Macroscopic field distribution"
      
      write(fh_ac, '("#",1X,A,":",1X,A)') "IX,IY,IZ", "Coordinate"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Ac", "Vector potential field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "E", "Electric field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "B", "Magnetic field"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eex", "Electron excitation energy"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eabs", "Absorbed energy"
      write(fh_ac, '("#",1X,A,":",1X,A)') "Eemf", "Total EM field energy"
      !if(use_ehrenfest_md=='y') then
      !   write(fh_ac, '("#",1X,A,":",1X,A)') "Jmi", "Matter current density of Ion"
      !   write(fh_ac, '("#",1X,A,":",1X,A)') "Tmp_ion", "Temperature of Ion"
      !endif
      write(fh_ac, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 1, "Ix", "none", &
        & 2, "Iy", "none", &
        & 3, "Iz", "none", &
        & 4, "Ac_x", "a.u.", & !!, trim(t_unit_ac%name), &
        & 5, "Ac_y", "a.u.", & !!, trim(t_unit_ac%name), &
        & 6, "Ac_z", "a.u.", & !!, trim(t_unit_ac%name), &
        & 7, "E_x", "a.u.", & !!, trim(t_unit_current%name), &
        & 8, "E_y", "a.u.", & !!, trim(t_unit_current%name), &
        & 9, "E_z", "a.u.", & !!, trim(t_unit_current%name), &
        & 10, "B_x", "a.u.", & !!, trim(t_unit_current%name), &
        & 11, "B_y", "a.u.", & !!, trim(t_unit_current%name), &
        & 12, "B_z", "a.u.", & !!, trim(t_unit_current%name), &
        & 13, "Jm_x", "a.u.", & !!, trim(t_unit_current%name), &
        & 14, "Jm_y", "a.u.", & !!, trim(t_unit_current%name), &
        & 15, "Jm_z", "a.u.", & !!, trim(t_unit_current%name), &
        & 16, "Eex", "a.u./unitcell", & !!, trim(t_unit_current%name), &
        & 17, "Eabs", "a.u./unitcell", & !!, trim(t_unit_current%name), &
        & 18, "Eemf", "a.u./unitcell" !!, trim(t_unit_current%name)
      !if(use_ehrenfest_md=='y') then
      !  write(fh_ac, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
      !  & 19, "Jmi_x", "a.u.", & !!, trim(t_unit_current%name), &
      !  & 20, "Jmi_y", "a.u.", & !!, trim(t_unit_current%name), &
      !  & 21, "Jmi_z", "a.u.", & !!, trim(t_unit_current%name), &
      !  & 22, "Tmp_ion", "K"
      !endif
      write(fh_ac,*)

      !! TODO: Support the automatic unit-system conversion of _ac.data files

      do iiz_m = nz1_m, nz2_m
        do iiy_m = ny1_m, ny2_m
          do iix_m = nx1_m, nx2_m
            write(fh_ac,'(I6,1X,I6,1X,I6,99(1X,E23.15E3))',advance='no')  &
              & iix_m, iiy_m, iiz_m, &
              & data_out(1:ndata_out_column, iix_m, iiy_m, iiz_m, ipos)
           !if(use_ehrenfest_md=='y') then
           !    !Add Jm_ion_xyz and Temperature_ion for MD (Later)
           !endif
            write(fh_ac,*)
          end do
        end do
      end do
      close(fh_ac)
    end if
   
    return
  end subroutine write_data_out


  subroutine store_data_vac_ac()
    implicit none
    ! Export the Ac field of detecting point
    if(comm_is_root(nproc_id_global)) then
      data_vac_Ac(1:3, 1, iter) = Ac_ms(1:3,ix_detect_l,iy_detect,iz_detect)
      data_vac_Ac(1:3, 2, iter) = Ac_ms(1:3,ix_detect_r,iy_detect,iz_detect)
    end if
  end subroutine store_data_vac_ac
  
  
  subroutine trace_ms_calculation()
    implicit none
    integer :: iix_m, iiy_m, iiz_m
    real(8) :: rrx, rry, rrz, eem, sem
    real(8) :: sx1, sy1, sz1
    real(8) :: sx2, sy2, sz2
    real(8) :: cx1, cy1, cz1
    real(8) :: cx2, cy2, cz2
    real(8) :: wx, wy, wz
    
    write(*,'(1X,A)') '-----------------------------------------------'
    write(*,'(1X,A,I6,A,I6)') 'Multiscale iter =', iter, '/', Nt
    
    write(*,'(1X,A)') 'Microscopic system:'
    
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Excitation energy:', &
    & "Eex =", total_energy_elec, &
    & "Diff =", total_energy_elec - total_energy_elec_old, &
    & "[au]"

    write(*,'(1X,A)') 'Macroscopic system:'
    
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Absorbed energy:', &
      & "Eabs =", total_energy_absorb, &
      & "Diff =", total_energy_absorb - total_energy_absorb_old, &
      & "[au]"
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Field energy:', &
      & "Eemf =", total_energy_elemag, &
      & "Diff =", total_energy_elemag - total_energy_elemag_old, &
      & "[au]"
    write(*,'(2X,A20,2(A10,ES15.6E3),A5)') 'Total EM energy:', &
      & "Etot =", total_energy_em, &
      & "Diff =", total_energy_em - total_energy_em_old, &
      & "[au]"
      
    sem = 0d0;
    sx1 = 0d0; sy1 = 0d0; sz1 = 0d0
    sx2 = 0d0; sy2 = 0d0; sz2 = 0d0
!$omp parallel do default(shared) collapse(3) &
!$omp private(iix_m, iiy_m, iiz_m, rrx, rry, rrz, eem) &
!$omp reduction(+: sx1, sy1, sz1, sx2, sy2, sz2, sem)
    do iiz_m = nz1_m, nz2_m
      do iiy_m = ny1_m, ny2_m
        do iix_m = nx1_m, nx2_m
          rrx = iix_m * HX_m
          rry = iiy_m * HY_m
          rrz = iiz_m * HZ_m
          
          eem = Energy_elemag_ms(iix_m, iiy_m, iiz_m)
          sem = sem + eem
          
          sx1 = sx1 + eem * rrx
          sy1 = sy1 + eem * rry
          sz1 = sz1 + eem * rrz
          
          sx2 = sx2 + eem * rrx ** 2
          sy2 = sy2 + eem * rry ** 2
          sz2 = sz2 + eem * rrz ** 2
        end do
      end do
    end do
!$omp end parallel do

    if (0d0 < sem) then
      cx1 = sx1 / sem; cy1 = sy1 / sem; cz1 = sz1 / sem 
      cx2 = sx2 / sem; cy2 = sy2 / sem; cz2 = sz2 / sem 
    else
      cx1 = 0d0; cy1 = 0d0; cz1 = 0d0
      cx2 = 0d0; cy2 = 0d0; cz2 = 0d0
    endif
    wx = sqrt(abs(cx2-cx1**2))
    wy = sqrt(abs(cy2-cy1**2))
    wz = sqrt(abs(cz2-cz1**2))
    write(*, '(1X,A)') "Position of wavepacket:"
    write(*,'(2X,A20,A10,3(ES12.3E3,","),A5)') &
      & 'Central position:', "<r> =", cx1, cy2, cz2, "[au]"
    write(*,'(2X,A20,A10,3(ES12.3E3,","),A5)') &
      & 'Spatial spreading:', "rms =", wx, wy, wz, "[au]"
    return 
  end subroutine trace_ms_calculation


  !===============================================================

  subroutine write_data_local_ac_jm()
    use salmon_file
    implicit none
    integer :: fh_ac_m
    integer :: iimacro, iiter
    character(100) :: file_ac_m

    !do iimacro = nmacro_s, nmacro_e
    if(flag_ms_ff_LessPrint) then !AY

       do iimacro = nmacro_s, nmacro_e

          if(iimacro.ne.1) cycle

          if(comm_is_root(nproc_id_tdks)) then

          if(iimacro.ge.nmacro_s .and. iimacro.le.nmacro_e) then

          write(file_ac_m,"(A,A,'_Ac_M.data')")trim(dir_ms_M(iimacro)),trim(SYSname)
          fh_ac_m = open_filehandle(file_ac_m)

          write(fh_ac_m, '("#",1X,A)') "Local variable at macro point"
        
          write(fh_ac_m, "('#',1X,A,':',3(1X,I6))") "Macropoint", macropoint(1:3, iimacro)
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Ac", "External vector potential field"
          if(use_ehrenfest_md=='y') then
            write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jmi", "Matter current density of Ion"
            write(fh_ac_m, '("#",1X,A,":",1X,A)') "Tmp_ion", "Temperature of Ion"
          endif

          write(fh_ac_m, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Ac_x", trim(t_unit_ac%name), &
          & 3, "Ac_y", trim(t_unit_ac%name), &
          & 4, "Ac_z", trim(t_unit_ac%name), &
          & 5, "Jm_x", trim(t_unit_current%name), &
          & 6, "Jm_y", trim(t_unit_current%name), &
          & 7, "Jm_z", trim(t_unit_current%name)
          if(use_ehrenfest_md=='y') then
            write(fh_ac_m, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
            & 8, "Jmi_x", trim(t_unit_current%name), &
            & 9, "Jmi_y", trim(t_unit_current%name), &
            &10, "Jmi_z", trim(t_unit_current%name), &
            &11, "Tmp_ion", "K"
          endif
          write(fh_ac_m,*)
          do iiter = 0, Nt
            write(fh_ac_m, "(F16.8,6(1X,ES22.14E3,1X))",advance='no') &
            & iiter * dt * t_unit_time%conv, &
            & data_local_Ac(1:3, iimacro, iiter) * t_unit_ac%conv, &
            & data_local_jm(1:3, iimacro, iiter) * t_unit_current%conv
            if(use_ehrenfest_md=='y') then
              write(fh_ac_m, "(F16.8,6(1X,ES22.14E3,1X))",advance='no') &
              data_local_jm_ion(1:3, iimacro, iiter) * t_unit_current%conv, &
              data_local_Tmp_ion(iimacro, iiter)
            endif
            write(fh_ac_m,*)
          end do !iiter

          close(fh_ac_m)

          endif !<--if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e)

          endif ! nproc_id_tdks

       enddo !iimacro

       call comm_sync_all

    else

    do igroup=1,ndivide_macro
       do i=1,nmacro_write_group

          if(comm_is_root(nproc_id_tdks)) then

          iimacro = (igroup-1)*nmacro_write_group + i
          if(iimacro.ge.nmacro_s .and. iimacro.le.nmacro_e) then

          write(file_ac_m,"(A,A,'_Ac_M.data')")trim(dir_ms_M(iimacro)),trim(SYSname)
          fh_ac_m = open_filehandle(file_ac_m)

          write(fh_ac_m, '("#",1X,A)') "Local variable at macro point"
        
          write(fh_ac_m, "('#',1X,A,':',3(1X,I6))") "Macropoint", macropoint(1:3, iimacro)
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jm", "Matter current density"
          write(fh_ac_m, '("#",1X,A,":",1X,A)') "Ac", "External vector potential field"
          if(use_ehrenfest_md=='y') then
            write(fh_ac_m, '("#",1X,A,":",1X,A)') "Jmi", "Matter current density of Ion"
            write(fh_ac_m, '("#",1X,A,":",1X,A)') "Tmp_ion", "Temperature of Ion"
          endif

          write(fh_ac_m, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Ac_x", trim(t_unit_ac%name), &
          & 3, "Ac_y", trim(t_unit_ac%name), &
          & 4, "Ac_z", trim(t_unit_ac%name), &
          & 5, "Jm_x", trim(t_unit_current%name), &
          & 6, "Jm_y", trim(t_unit_current%name), &
          & 7, "Jm_z", trim(t_unit_current%name)
          if(use_ehrenfest_md=='y') then
            write(fh_ac_m, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
            & 8, "Jmi_x", trim(t_unit_current%name), &
            & 9, "Jmi_y", trim(t_unit_current%name), &
            &10, "Jmi_z", trim(t_unit_current%name), &
            &11, "Tmp_ion", "K"
          endif
          write(fh_ac_m,*)
          do iiter = 0, Nt
            write(fh_ac_m, "(F16.8,6(1X,E23.15E3,1X))",advance='no') &
            & iiter * dt * t_unit_time%conv, &
            & data_local_Ac(1:3, iimacro, iiter) * t_unit_ac%conv, &
            & data_local_jm(1:3, iimacro, iiter) * t_unit_current%conv
            if(use_ehrenfest_md=='y') then
              write(fh_ac_m, "(F16.8,6(1X,E23.15E3,1X))",advance='no') &
              data_local_jm_ion(1:3, iimacro, iiter) * t_unit_current%conv, &
              data_local_Tmp_ion(iimacro, iiter)
            endif
            write(fh_ac_m,*)
          end do !iiter

          close(fh_ac_m)

          endif !<--if(imacro.ge.nmacro_s .and. imacro.le.nmacro_e)

          endif ! nproc_id_tdks

       enddo !i
       call comm_sync_all
    enddo !igroup
    endif

  end subroutine
  
  
  subroutine write_data_vac_ac()
    use salmon_file
    implicit none
    integer :: fh_ac_vac
    integer :: iiter
    character(100) :: file_ac_vac
    
    if (comm_is_root(nproc_id_global)) then
      write(file_ac_vac,"(A, A, '_Ac_vac.data')") trim(dir_ms_RT_Ac),trim(SYSname)
      fh_ac_vac = open_filehandle(file_ac_vac)
      write(fh_ac_vac, '("#",1X,A)') "Ac vacuum region"
      write(fh_ac_vac, '("#",1X,A)') "Data of Ac field at the end of media"
      
      write(fh_ac_vac, "('#',1X,A,':',3(1X,I6))") "L", ix_detect_l, iy_detect, iz_detect
      write(fh_ac_vac, "('#',1X,A,':',3(1X,I6))") "R", ix_detect_r, iy_detect, iz_detect
      
      write(fh_ac_vac, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Ac_x(L)", trim(t_unit_ac%name), &
        & 3, "Ac_y(L)", trim(t_unit_ac%name), &
        & 4, "Ac_z(L)", trim(t_unit_ac%name), &
        & 5, "Ac_x(R)", trim(t_unit_ac%name), &
        & 6, "Ac_y(R)", trim(t_unit_ac%name), &
        & 7, "Ac_z(R)", trim(t_unit_ac%name)
      do iiter = 0, Nt
        write(fh_ac_vac, "(F16.8,6(1X,E23.15E3,1X))") &
          & iiter * dt * t_unit_time%conv, &
          & data_vac_Ac(1:3, 1, iiter) * t_unit_ac%conv, &
          & data_vac_Ac(1:3, 2, iiter) * t_unit_ac%conv
      end do
      close(fh_ac_vac)
    end if
    call comm_sync_all  
  end subroutine write_data_vac_ac

  subroutine write_xyz_ms(comment,action,rvf,imacro)

  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use Global_Variables
    use inputoutput, only: au_length_aa
    use salmon_global, only: SYSname,iflag_atom_coor,ntype_atom_coor_cartesian,ntype_atom_coor_reduced
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file
    implicit none
    integer :: ia,unit_xyz,imacro 
    character(3) :: action,rvf
    character(1024) :: file_trj
    character(*) :: comment

    if(comm_is_root(nproc_id_tdks)) then

      unit_xyz = 1000 + imacro

      if(action=='new') then

        write(file_trj, "(A,A,'_trj.xyz')") trim(dir_ms_M(imacro)),trim(SYSname)
        open(unit_xyz,file=trim(file_trj),status="unknown")

      else if(action=='add') then

         write(unit_xyz,*) NI
         write(unit_xyz,*) trim(comment)
         do ia=1,NI
            if(      rvf=="r  " ) then
               write(unit_xyz,100) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa
            else if( rvf=="rv " ) then
               write(unit_xyz,110) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa,&
                           & velocity_m(1:3,ia,imacro)
            else if( rvf=="rvf" ) then
               write(unit_xyz,120) trim(atom_name(ia)), &
                           & Rion_m(1:3,ia,imacro)*au_length_aa,&
                           & velocity_m(1:3,ia,imacro), force_m(1:3,ia,imacro)
            endif
         enddo

      else if(action=='end') then

         close(unit_xyz)

      endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

    endif

!    call comm_sync_all
  end subroutine

  subroutine write_ini_coor_vel_ms_for_restart(imacro)
  !print atomic coordinate and velocity at the last step just to set_ini_coor_vel option
    use Global_Variables
    use inputoutput, only: au_length_aa
!    use salmon_global, only: SYSname,iflag_atom_coor,ntype_atom_coor_cartesian,ntype_atom_coor_reduced
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file
    implicit none
    integer :: ia,unit_ini_rv,imacro 
    character(1024) :: file_ini_rv
    real(8) :: ulength_to_au

    if(comm_is_root(nproc_id_tdks)) then

      unit_ini_rv = 5000 + imacro

      write(file_ini_rv, "(A,'last_coor_vel.data')") trim(dir_ms_M(imacro))
      open(unit_ini_rv,file=trim(file_ini_rv),status="unknown")

      select case(unit_length)
      case('au','a.u.')
         ulength_to_au   = 1d0
      case('AA','angstrom','Angstrom')
         ulength_to_au   = 1d0/au_length_aa
      end select
      
      do ia=1,NI
         write(unit_ini_rv,110)  &
              & Rion_m(1:3,ia,imacro)/ulength_to_au,&
              & velocity_m(1:3,ia,imacro)
      enddo
110   format(6e24.12)

      close(unit_ini_rv)

    endif

!    call comm_sync_all
  end subroutine

  subroutine dt_evolve_MD_1_MS(iter,imacro)
  ! Velocity Verlet integrator for ion dynamics in multi-scale
  ! Currentry, NVT ensemble is not supported
    use md_ground_state
    implicit none
    integer :: iter,imacro,ia
    real(8) :: dt_h,mass_au
  character(100) :: comment_line

    dt_h = dt/2d0
    velocity(:,:) = velocity_m(:,:,imacro)
    force(:,:)    = force_m(:,:,imacro)
    Rion(:,:)     = Rion_m(:,:,imacro)
    Rion_eq(:,:)  = Rion_eq_m(:,:,imacro)
    dRion(:,:,iter:iter+1)  = dRion_m(:,:,iter:iter+1,imacro)


    !!NHC act on velocity with dt/2
    !if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
    !   call apply_nose_hoover_velocity(dt_h)
    !endif
    !update ion velocity with dt/2

    if(iter==0) then
       dRion(:,:,iter-1)= dRion(:,:,iter) - velocity(:,:)*dt
       dRion_m(:,:,iter-1,imacro) = dRion(:,:,iter-1)
    endif

    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
      !velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h &
       &                               - friction/mass_au * velocity(:,ia) * dt_h
    enddo

    !velocity scaling
    if(step_velocity_scaling>=1 .and. mod(iter,step_velocity_scaling)==0)then
       call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity(:,:))
       call apply_velocity_scaling_ion(Temperature_ion,velocity(:,:))
    endif

    !update ion coordinate with dt
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
       dRion(:,ia,iter+1) = dRion(:,ia,iter) + velocity(:,ia)*dt
       Rion(:,ia) = Rion_eq(:,ia) + dRion(:,ia,iter+1)
    enddo

!xx    !update pseudopotential (check .... different in each macrogrid)
!xx    if (mod(iter,step_update_ps)==0 ) then
!xx       call prep_ps_periodic('update_all       ')
!xx    else if (mod(iter,step_update_ps2)==0 ) then
!xx      call prep_ps_periodic('update_wo_realloc')
!xx    endif

    !!NHC act on thermostat with dt
    !if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
    !   call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
    !   call apply_nose_hoover_thermostat(Temperature_ion,dt)
    !   Enh_gkTlns = Enh_gkTlns + gkT * xi_nh*dt
    !   Enh        = Enh_gkTlns + 0.5d0 * Qnh * xi_nh*xi_nh
    !endif

    velocity_m(:,:,imacro) = velocity(:,:)
    Rion_m(:,:,imacro)     = Rion(:,:)
    dRion_m(:,:,iter:iter+1,imacro) = dRion(:,:,iter:iter+1)
    
  end subroutine

  subroutine dt_evolve_MD_2_MS(aforce,iter,imacro)
    use md_ground_state
    implicit none
    integer :: iter,imacro, ia
    real(8) :: dt_h,Ework,aforce(3,NI),Temperature_ion, mass_au

    dt_h = dt/2d0
    velocity(:,:) = velocity_m(:,:,imacro)
    force(:,:)    = force_m(:,:,imacro)
    dRion(:,:,iter:iter+1)  = dRion_m(:,:,iter:iter+1,imacro)

    aforce(:,:) = 0.5d0*( aforce(:,:) + force(:,:) )

    !update ion velocity with dt/2
    Ework = 0d0
    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
      !velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
       velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h &
       &                               - friction/mass_au * velocity(:,ia) * dt_h
       Ework = Ework - sum(aforce(:,ia)*(dRion(:,ia,iter+1)-dRion(:,ia,iter)))
    enddo

    !!NHC act on velocity with dt/2
    !if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
    !   call apply_nose_hoover_velocity(dt_h)
    !endif

    if (stop_system_momt=='y') call remove_system_momentum(0)
    call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)


!    Eall=Eall+Tion
    
!    Eall_t(iter) = Eall
!    Tion_t(iter) = Tion
     data_local_Tmp_ion(imacro,iter) = Temperature_ion
!    Ework_integ_fdR(iter) = Ework_integ_fdR(iter-1) + Ework

    !if(ensemble=="NVT".and.thermostat=="nose-hoover")then
    !   Enh_t(iter)  = Enh
    !   Hnvt_t(iter) = Eall + Enh
    !endif

    velocity_m(:,:,imacro) = velocity(:,:)

  end subroutine

  subroutine Ion_Force_Ac_MS(imacro)
    implicit none
    integer :: ia,imacro
    real(8) :: E_tot(3),FionAc(3,NI)

    E_tot(:)= -(Ac_new_m(:,imacro)-Ac_old_m(:,imacro))/(2*dt)

!    Eelemag= aLxyz*sum(E_tot(iter,:)**2)/(8.d0*Pi)
!    Eall   = Eall + Eelemag

    do ia=1,NI
      FionAc(:,ia)=Zps(Kion(ia))*E_tot(:)
      force_m(:,ia,imacro) = force_m(:,ia,imacro) + FionAc(:,ia)
    enddo

  end subroutine

  subroutine Ion_Force_ForceField_MS(imacro)
    implicit none
    integer :: ia,imacro
    real(8) :: E_tot(3),Fion_FF(3,NI), mass_au, omg_au, k_au

    force_m(:,:,imacro) = 0d0

    omg_au  = Omg_dt/dt

    do ia=1,NI
       mass_au = umass*Mass(Kion(ia))
       k_au    = mass_au * omg_au**2d0
       Fion_FF(:,ia)= - k_au *(Rion_m(:,ia,imacro) - Rion_eq0(:,ia))
       force_m(:,ia,imacro) = force_m(:,ia,imacro) + Fion_FF(:,ia)
    enddo

  end subroutine

  subroutine Ion_Force_RamanTensor_MS(imacro)
    implicit none
    integer :: ia,imacro
    real(8) :: E_tot(3),Fion_E2(3,NI)

    E_tot(:)= -(Ac_new_m(:,imacro)-Ac_old_m(:,imacro))/(2*dt)

    !(in the case of diamond,Si)
    do ia=1,NI
      Fion_E2(1,ia) = 0.5d0 * aLxyz * c_pmode(ia) / dble(NI) &
                       * ( dchidq(2,3) * E_tot(2)*E_tot(3) &
                         + dchidq(3,2) * E_tot(3)*E_tot(2) )
      force_m(1,ia,imacro) = force_m(1,ia,imacro) + Fion_E2(1,ia)
    enddo

  end subroutine

end subroutine raman_maxwell_ms
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
end module control_ms_raman
