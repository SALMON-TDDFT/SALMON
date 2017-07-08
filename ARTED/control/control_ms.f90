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
subroutine tddft_maxwell_ms
  use Global_Variables
  use timer
  use opt_variables
  use environment
  use performance_analyzer
  use salmon_parallel
  use salmon_communication
  use misc_routines
  use restart, only: prep_restart_write

  implicit none
  integer :: iter
  integer :: ix_m,iy_m,ixy_m
  integer :: index, n
  character(len=128) :: fmt
  
  real(8) calc_pulse_xcenter
  logical :: flag_shutdown = .false.


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
      call init_Ac_ms_2dc()
    else
      call init_Ac_ms
    endif
  
    rho_gs(:)=rho(:)
    
    Vloc_old(:,1) = Vloc(:); Vloc_old(:,2) = Vloc(:)
    
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

!reentrance

    position_option='rewind'
    entrance_iter=-1
    call reset_rt_timer
  end if
  
  ! Output filename
  write(file_energy_transfer, "(A,'energy-transfer.out')") trim(directory)
  write(file_ac_vac, "(A,'Ac_Vac.out')") trim(directory)
  write(file_ac_vac_back, "(A,'Ac_Vac_back.out')") trim(directory)
  write(file_ac_m, "(A,'Ac_M',I6.6,'.out')") trim(process_directory), NXY_s
  
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
        call Ion_Force_omp(Rion_update_rt,calc_mode_rt,ixy_m)
        if (mod(iter, Nstep_write) == 0) then
          call Total_Energy_omp(Rion_update_rt,calc_mode_rt,ixy_m)
        end if
      else
        if (mod(iter, Nstep_write) == 0) then
!$acc update self(zu)
          call Total_Energy_omp(Rion_update_rt,calc_mode_rt,ixy_m)
          call Ion_Force_omp(Rion_update_rt,calc_mode_rt,ixy_m)
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
        call k_shift_wf(Rion_update_rt,2,zu_m(:,:,:,ixy_m))
        if(comm_is_root(nproc_id_tdks))then ! sato
          excited_electron_l(ix_m,iy_m)=sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if ! sato
      else if (iter == Nt ) then
        call k_shift_wf(Rion_update_rt,2,zu_m(:,:,:,ixy_m))
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
      call comm_bcast(flag_shutdown,nproc_group_global)
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
    if (flag_shutdown) then 
      call comm_sync_all
      write(*,*) nproc_id_global,'iter =',iter
      iter_now=iter
!$acc update self(zu)
      call timer_end(LOG_DYNAMICS)
      call prep_restart_write
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
        flag_shutdown =.true.
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
      call prep_restart_write
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


  subroutine reset_rt_timer
    implicit none
    integer :: i
    do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
      call timer_reset(i)
    end do
  end subroutine
end subroutine tddft_maxwell_ms
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
end module control_ms
