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
!This file is "md_gs.f90"
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

!-----------------------------------------------------------------
! Optimization in the ground state by conjugate gradient method
!-----------------------------------------------------------------
module md_ground_state
  use Global_Variables
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use salmon_parallel, only: nproc_id_global
  implicit none
contains

subroutine calc_md_ground_state
  use timer
  use opt_variables
  use performance_analyzer
  use misc_routines, only: get_wtime
  use salmon_global, only: format3d, out_dns, out_dns_rt, out_dns_rt_step
  use ground_state
  use io_gs_wfn_k
  implicit none
  integer :: entrance_it,it,ia
  real(8) :: Vion_gs, Ework, Ework_integ_fdR_gs, dt_h, aforce(3,NI), Eall00
  real(8) :: Temperature_ion, mass_au, Rion_prev(3,NI)
  real(8) :: Htot, Enh, Enh_gkTlns, gkT, Qnh
  character(100) :: comment_line


  iflag_gs_init_wf = 1  !flag to skip giving randam number initial guess
  PrLv_scf = 0

  position_option='rewind'
  entrance_it=-1
  Rion_update_rt = rion_update_on

  ! Export electronic density
  if(out_dns == 'y') call write_density(it,'gs')

  ! Export to file_trj (initial step)
  if (out_rvf_rt=='y')then
       call Ion_Force_omp(Rion_update_rt,calc_mode_gs)
       write(comment_line,110) -1, 0.0d0
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
       &  write(comment_line,112) trim(comment_line), xi_nh
       call write_xyz(comment_line,"new","rvf")
       call write_xyz(comment_line,"add","rvf")
  endif

  it=entrance_it+1
  dRion(:,:,it-1) = dRion(:,:,it) - velocity(:,:)*dt  !only for euler
  dt_h       = dt*0.5d0
  Eall00     = Eall0
  Enh_gkTlns = 0d0
  Enh        = 0d0

  if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
     gkT = 3d0*NI * kB/hartree2J*temperature0_ion
     Qnh = gkT * thermostat_tau**2d0
  endif

  if(comm_is_root(nproc_id_global)) then
     write(*,'(a)') "# time, Tion, Vion_gs, Ework, Enh, Eall, Eall-Eall0, Htot, Temperature_ion [all in a.u.]"
  endif

  call comm_sync_all

  ! === Main loop of time step ===
  do it = entrance_it+1, Nt

     ! Velocity Verlet integrator

     !NHC act on velocity with dt/2
     if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
        call apply_nose_hoover_velocity(dt_h)
     endif

     do ia=1,NI
        mass_au = umass*Mass(Kion(ia))
        velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
     enddo

     if(step_velocity_scaling>=1 .and. mod(it,step_velocity_scaling)==0) then
        call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
        call apply_velocity_scaling_ion(Temperature_ion,velocity)
     endif

     Rion_prev(:,:) = Rion(:,:)
     do ia=1,NI
        Rion(:,ia) = Rion(:,ia) + velocity(:,ia)*dt
     enddo

     !put SHAKE here in future

     !NHC act on thermostat with dt
     if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
        call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
        call apply_nose_hoover_thermostat(Temperature_ion,dt)
        Enh_gkTlns = Enh_gkTlns + gkT * xi_nh*dt
        Enh        = Enh_gkTlns + 0.5d0 * Qnh * xi_nh*xi_nh
     endif


     !update force (electric state) with updated coordinate
     aforce(:,:) = force(:,:)
     call prep_ps_periodic('update_all       ')
     call read_write_gs_wfn_k(iflag_read)
     call calc_ground_state
     call read_write_gs_wfn_k(iflag_write)

     call Ion_Force_omp(Rion_update_rt,calc_mode_gs)
     call Total_Energy_omp(Rion_update_rt,calc_mode_gs)

     aforce(:,:) = 0.5d0*( aforce(:,:) + force(:,:) )

     Ework = 0d0
     do ia=1,NI
        mass_au = umass*Mass(Kion(ia))
        velocity(:,ia) = velocity(:,ia) + force(:,ia)/mass_au * dt_h
        Ework = Ework - sum(aforce(:,ia)*(Rion(:,ia)-Rion_prev(:,ia)))
     enddo

     !NHC act on velocity with dt/2
     if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
        call apply_nose_hoover_velocity(dt_h)
     endif

     !remove system momentum
     if(stop_system_momt=='y') call remove_system_momentum(0)


     call cal_Tion_Temperature_ion(Tion,Temperature_ion,velocity)
     Vion_gs  = Eall
     Eall     = Eall + Tion
     Htot     = Eall + Enh
     Ework_integ_fdR_gs = Ework_integ_fdR_gs + Ework


     !---Write section---
     ! Export to file_rt_data
     call write_md_gs_data(it,Vion_gs,Tion,Temperature_ion,Eall,Eall00,Ework_integ_fdR_gs,Enh,Htot)

     ! Export to standard log file
     if(comm_is_root(nproc_id_global)) then
        write(*,120) it*dt,Tion,Vion_gs,Ework_integ_fdR_gs,Enh,Eall,Eall-Eall00,Htot,Temperature_ion
120     format(1x,f10.4, 7e20.10E3,f12.3)
     endif

     ! Export to file_trj
     if (out_rvf_rt=='y' .and. mod(it,out_rvf_rt_step)==0)then
        write(comment_line,110) it, it*dt
110     format("#md-gs  step=",i8,"   time",e16.6)
        if(ensemble=="NVT" .and. thermostat=="nose-hoover") &
        &  write(comment_line,112) trim(comment_line), xi_nh
112     format(a,"  xi_nh=",e18.10)
        call write_xyz(comment_line,"add","rvf")
     endif

     ! Export electronic density (cube or vtk)
     if(out_dns_rt=='y' .and. mod(it,out_dns_rt_step)==0) then
        call write_density(it,'rt')
     end if

  enddo  !main loop of it

  call print_restart_data_md_gs
  call comm_sync_all

  if(comm_is_root(nproc_id_global)) write(*,*) 'This calculation is shutdown successfully!'

end subroutine calc_md_ground_state

  subroutine apply_nose_hoover_velocity(dt_h) !assuming half dt
    implicit none
    integer :: ia
    real(8) :: tmp_exp,dt_h 

    tmp_exp = exp(-xi_nh*dt_h)
    do ia=1,NI
       velocity(:,ia) = velocity(:,ia) * tmp_exp
    enddo

    return
  end subroutine apply_nose_hoover_velocity

  subroutine apply_nose_hoover_thermostat(Temp_ion,dt_f) !assuming full dt
    implicit none
    real(8) :: Temp_ion,dt_f 

    xi_nh = xi_nh + (Temp_ion/temperature0_ion-1d0)/(thermostat_tau**2d0)*dt_f 

    return
  end subroutine apply_nose_hoover_thermostat

  subroutine write_md_gs_data(it,Vion_gs_t,Tion_gs_t,Temp_ion_t,Eall_gs_t,Eall_gs0,Ework_gs_t,Enh_t,Htot_gs_t)
    use inputoutput, only: t_unit_time,t_unit_energy
    implicit none
    integer :: fh_rt=202, it
    real(8) :: Vion_gs_t,Tion_gs_t,Temp_ion_t,Eall_gs_t,Ework_gs_t,Eall_gs0,Enh_t,Htot_gs_t

    if(comm_is_root(nproc_id_global)) then
      if(it==0) then
         open(fh_rt,file=trim(file_rt_data),status="unknown")
      
         write(fh_rt, '("#",1X,A)') "Real time calculation (Ground state MD)"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Tion", "Kinetic energy of ions"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Vion", "Potential energy of ions"
         write(fh_rt, '("#",1X,A,":",1X,A)') "E_work", "Work energy of ions(sum f*dr)"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Enh", "Energy of NH thermostat"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Eall", "Total energy"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Eall0", "Initial energy"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Hnvt'", "Hnvt using E_work"
         write(fh_rt, '("#",1X,A,":",1X,A)') "Temperature_ion", "Temperature of ions"
         write(fh_rt, '("#",99(1X,I0,":",A,"[",A,"]"))') &
              &  1, "Time",   trim(t_unit_time%name), &
              &  2, "Tion",   trim(t_unit_energy%name), &
              &  3, "Vion",   trim(t_unit_energy%name), &
              &  4, "E_work", trim(t_unit_energy%name), &
              &  5, "Enh",    trim(t_unit_energy%name), &
              &  6, "Eall",   trim(t_unit_energy%name), &
              &  7, "Eall-Eall0", trim(t_unit_energy%name), &
              &  8, "Hnvt",   trim(t_unit_energy%name), &
              &  9, "Hnvt'",  trim(t_unit_energy%name), &
              & 10, "Temperature_ion", "K"
      endif

      write(fh_rt, "(F16.8,99(1X,E23.15E3))") &
          & it * dt    * t_unit_time%conv, &
          & Tion_gs_t  * t_unit_energy%conv, &
          & Vion_gs_t  * t_unit_energy%conv, &
          & Ework_gs_t * t_unit_energy%conv, &
          & Enh_t      * t_unit_energy%conv, &
          & Eall_gs_t  * t_unit_energy%conv, &
          & (Eall_gs_t-Eall_gs0) * t_unit_energy%conv, &
          & Htot_gs_t  * t_unit_energy%conv, &
          & (Tion_gs_t+Ework_gs_t+Enh_t) * t_unit_energy%conv, &
          & Temp_ion_t
      flush(fh_rt)

      if(it==Nt) close(fh_rt)
    end if

    call comm_sync_all
    return
  end subroutine write_md_gs_data

  subroutine cal_Tion_Temperature_ion(Ene_ion,Temp_ion, vel)
     implicit none
     integer :: ia
     real(8) :: Ene_ion,Temp_ion, vel(3,NI)

     Ene_ion=0.d0
     do ia=1,NI
        Ene_ion = Ene_ion + 0.5d0 * umass*Mass(Kion(ia)) * sum(vel(:,ia)**2d0)
     enddo
     Temp_ion = Ene_ion * 2d0 / (3d0*NI) / (kB/hartree2J)

     return
  end subroutine cal_Tion_Temperature_ion
  
  subroutine apply_velocity_scaling_ion(Temp_ion,vel)
     implicit none
     real(8) :: Temp_ion, vel(3,NI), fac_vscaling

     fac_vscaling = sqrt(temperature0_ion/Temp_ion)
     vel(:,:) = vel(:,:) * fac_vscaling

     return
  end subroutine apply_velocity_scaling_ion

  subroutine print_restart_data_md_gs
    use salmon_communication, only: comm_is_root
    use salmon_parallel, only: nproc_id_global
    use inputoutput, only: au_length_aa
    implicit none
    integer :: ia,unit_atomic_coor_tmp=201,ik,j
    real(8) :: tmpr(3), uconv
    character(100)  :: char_atom

    if(comm_is_root(nproc_id_global)) then

       select case(iflag_atom_coor)
       case(ntype_atom_coor_cartesian)
          open(unit_atomic_coor_tmp,file='.atomic_coor.tmp',status='old')
       case(ntype_atom_coor_reduced)
          open(unit_atomic_coor_tmp,file='.atomic_red_coor.tmp',status='old')
       end select

       if(unit_length=='AA')then
          uconv = au_length_aa
       else  !au
          uconv = 1d0
       endif

       write(*,*) 
       write(*,9000) "##------ Restarting Data for MD-GS -------"
       write(*,9000) "# Coordinate in specified unit in input"
       write(*,9000) "# Copy into new input file"
       write(*,9000) "&atomic_coor"
       do ia = 1,NI
          read(unit_atomic_coor_tmp,*) char_atom, (tmpr(j),j=1,3),ik
          write(*,7000) trim(char_atom), Rion(1:3,ia)*uconv, ik
7000      format("     '",a,"'  ",3f18.10,i4)
       enddo
       write(*,9000) "/"

       write(*,*) 
       write(*,9000) "# Velocity (atoms, thermostat if nose-hoover ooption) in [au]"
       write(*,9000) "# Copy to separated file used in file_ini_vel option with set_ini_velocity='r'"
       do ia = 1,NI      
          write(*,8000) velocity(1:3,ia)
       enddo
8000   format(3e18.10)
       if(ensemble=="NVT" .and. thermostat=="nose-hoover") then
          write(*,8000) xi_nh
       endif
       write(*,*) 


       close(unit_atomic_coor_tmp)
9000   format(a)
    endif

    return
  end subroutine print_restart_data_md_gs

  subroutine remove_system_momentum(flag_print_check)
    ! remove center of mass and momentum of whole system
    use salmon_communication, only: comm_is_root
    use salmon_parallel, only: nproc_id_global
    implicit none
    integer :: ia, flag_print_check
    real(8) :: v_com(3), sum_mass, mass_au

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

    !rotation of system is removed
    !(this is only for isolated system)--- do nothing

    !(check velocity of center of mass)
    if(flag_print_check==1) then

       v_com(:)=0d0
       do ia=1,NI
          v_com(:) = v_com(:) + umass*Mass(Kion(ia)) * velocity(:,ia)
       enddo
       v_com(:) = v_com(:)/sum_mass
       if (comm_is_root(nproc_id_global)) &
       write(*,*) "    v_com =",real(v_com(:))

    endif

  end subroutine remove_system_momentum

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

end module md_ground_state
