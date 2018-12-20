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
!This file is "io_gs_wfn_k.f90"
!This file contains I/O routines for ground state wavefunctions
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module io_gs_wfn_k
  use salmon_global
  use salmon_communication
  use salmon_parallel
  use global_variables
  use misc_routines
  implicit none

  character(256) :: gs_wfn_directory
  character(256) :: gs_wfn_file, occ_file
  integer,parameter :: nfile_gs_wfn = 41
  integer,parameter :: nfile_occ    = 42

  integer,parameter :: iflag_read = 0
  integer,parameter :: iflag_write= 1

contains
  subroutine read_write_gs_wfn_k(iflag_read_write)
    implicit none
    integer,intent(in) :: iflag_read_write
    integer :: ik
  ! temporal communicator for multi-scale (same k-points at different macro point)
    integer :: nproc_group_kpoint_ms
    integer :: nproc_id_kpoint_ms
    integer :: nproc_size_kpoint_ms

    write (gs_wfn_directory,'(A,A)') trim(directory),'/gs_wfn_k/'
    if(iflag_read_write == iflag_write)call create_directory(gs_wfn_directory)

    if(comm_is_root(nproc_id_global))then
      occ_file = trim(gs_wfn_directory)//'occupation'
      open(nfile_occ,file=trim(occ_file),form='unformatted')
      select case(iflag_read_write)
      case(iflag_write); write(nfile_occ)occ
      case(iflag_read ); read(nfile_occ)occ
      end select
      close(nfile_occ)
    end if

    select case(iflag_read_write)
    case(iflag_read )
      call comm_bcast(occ,nproc_group_global)
    end select

    if(use_ms_maxwell == 'n' .or. (use_ms_maxwell == 'y'.and. nmacro_s == 1))then
      do ik=NK_s,NK_e
        
        write (gs_wfn_file,'(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik,'.wfn'
        open(nfile_gs_wfn,file=trim(gs_wfn_file),form='unformatted')
        select case(iflag_read_write)
        case(iflag_write); write(nfile_gs_wfn)zu_GS(:,:,ik)
        case(iflag_read ); read(nfile_gs_wfn)zu_GS(:,:,ik)
        end select
        close(nfile_gs_wfn)

      end do

    end if

    if(iflag_read_write==iflag_write)return

    if(use_ms_maxwell == 'y')then
      nproc_group_kpoint_ms = comm_create_group(nproc_group_global, NK_s, nmacro_s)
      call comm_get_groupinfo(nproc_group_kpoint_ms, nproc_id_kpoint_ms, nproc_size_kpoint_ms)
      call comm_bcast(zu_GS,nproc_group_kpoint_ms)
    end if


    zu_GS0(:,:,:)=zu_GS(:,:,:)
    
    zu_t(:,:,:)=zu_GS(:,1:NBoccmax,:)
    Rion_eq=Rion
    !dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0  !moved to initialization
    call psi_rho_GS
    call Hartree
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)

    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    Vloc_GS(:)=Vloc(:)
    call Total_Energy_omp(rion_update_on,calc_mode_gs)
    call Ion_Force_omp(rion_update_on,calc_mode_gs)
    Eall0=Eall
    if(use_ms_maxwell == 'y') Eall0_m(:)=Eall
    Eall_GS0=Eall
    if(PrLv_scf==3 .and. comm_is_root(nproc_id_global)) write(*,*)'Eall =',Eall

  end subroutine read_write_gs_wfn_k

  subroutine read_gs_wfn_k_ms_each_macro_grid
    use inputoutput, only: au_length_aa
    implicit none
    integer :: ik,imacro,nfile_occ_ms,nfile_gs_wfn_ms
    integer :: unit_ini_rv,ia,j
    character(1024) :: file_ini_rv
    real(8) :: ulength_to_au

    do imacro = nmacro_s, nmacro_e

       nfile_occ_ms    = 7000 + imacro
       nfile_gs_wfn_ms = nfile_occ_ms
       write (gs_wfn_directory,'(A,A)') trim(dir_ms_M(imacro)),'/gs_wfn_k/'
       
       if(comm_is_root(nproc_id_tdks)) then
          occ_file = trim(gs_wfn_directory)//'occupation'
          open(nfile_occ_ms,file=trim(occ_file),form='unformatted')
          read(nfile_occ_ms) occ
          close(nfile_occ_ms)

          !! maybe must be out of "if(comm_is_root(nproc_id_tdks)) then"
          do ik=NK_s,NK_e
             write(gs_wfn_file,'(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik,'.wfn'
             open(nfile_gs_wfn_ms,file=trim(gs_wfn_file),form='unformatted')
             read(nfile_gs_wfn_ms) zu_GS(:,:,ik)
             close(nfile_gs_wfn_ms)
          end do
       end if

       call comm_bcast(occ,nproc_group_tdks)

       zu_GS0(:,:,:) = zu_GS(:,:,:)
       zu_t(:,:,:) = zu_GS(:,1:NBoccmax,:)

       if(set_ini_coor_vel=='y') then
 
           select case(unit_length)
           case('au','a.u.')
              ulength_to_au   = 1d0
           case('AA','angstrom','Angstrom')
              ulength_to_au   = 1d0/au_length_aa
           end select

           ! file for initial atomic coordinate and velocity in each macro grid
           ! must be "(parent dir)/multiscale/M00XXXX/ini_coor_vel.dat"
           ! the file format is
           !  r_x  r_y  r_z  v_x  v_y  v_z  (for 1th atom)
           !  r_x  r_y  r_z  v_x  v_y  v_z  (for 2th atom)
           !    ..................
           !  r_x  r_y  r_z  v_x  v_y  v_z  (for NI-th atom)
           ! (unit of coordinate is A or a.u., but velocity must be a.u.)
           if(comm_is_root(nproc_id_tdks)) then
              unit_ini_rv = 5000 + imacro
              file_ini_rv = trim(dir_ms_M(imacro))//'ini_coor_vel.dat'
              open(unit_ini_rv,file=trim(file_ini_rv),status="old")
              do ia=1,NI
                 read(unit_ini_rv,*) (Rion(j,ia),j=1,3),(velocity(j,ia),j=1,3)
              enddo
              Rion(:,:)    = Rion(:,:)*ulength_to_au
              Rion_eq(:,:) = Rion(:,:)
              close(unit_ini_rv)
           endif
           
           call comm_bcast(Rion,    nproc_group_tdks)
           call comm_bcast(Rion_eq, nproc_group_tdks)
           call comm_bcast(velocity,nproc_group_tdks)

           Rion_m(:,:,imacro)     = Rion(:,:)
           Rion_eq_m(:,:,imacro)  = Rion_m(:,:,imacro)
           velocity_m(:,:,imacro) = velocity(:,:)
           
       endif

       call psi_rho_GS
       call Hartree
       call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)

       Vloc(1:NL) = Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
       Vloc_GS(:) = Vloc(:)
       call Total_Energy_omp(rion_update_on,calc_mode_gs)
       call Ion_Force_omp(rion_update_on,calc_mode_gs)
       Eall0_m(imacro) = Eall
       Eall_GS0 = Eall

    enddo !imacro

  end subroutine read_gs_wfn_k_ms_each_macro_grid

  subroutine modify_initial_guess_copy_1stk_to_all
    implicit none
    integer :: ik,ik1,ib
    !integer :: nproc_group_kpoint_ms
    !integer :: nproc_id_kpoint_ms
    !integer :: nproc_size_kpoint_ms
    real(8) :: tot_occ
    real(8),allocatable :: occ_tmp(:,:)
    complex(8),allocatable :: zu_GS1(:,:,:)

    write (gs_wfn_directory,'(A,A)') trim(directory),'/gs_wfn_k/'

    if(comm_is_root(nproc_id_global))then
      occ_file = trim(gs_wfn_directory)//'occupation'

      !read only gammma point
      allocate(occ_tmp(NB,1))
      open(nfile_occ,file=trim(occ_file),form='unformatted')
      read(nfile_occ)occ_tmp
      close(nfile_occ)

      !over-write with all k-points
      tot_occ=0d0
      do ib=1,NB
         tot_occ = tot_occ + occ_tmp(ib,1)
      enddo
      occ_tmp(:,1) = occ_tmp(:,1)/tot_occ*dble(nelec)/dble(NK)

      do ik=1,NK
         occ(:,ik) = occ_tmp(:,1)
      enddo
      open(nfile_occ,file=trim(occ_file),form='unformatted')
      write(nfile_occ)occ
      close(nfile_occ)

      deallocate(occ_tmp)

   end if
      
   if(use_ms_maxwell == 'n' .or. (use_ms_maxwell == 'y'.and. nmacro_s == 0))then

     allocate(zu_GS1(NL,NB,1))
     if(comm_is_root(nproc_id_global))then
        !read the first k-point data 
        ik1=1
        write (gs_wfn_file,'(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik1,'.wfn'
        open(nfile_gs_wfn,file=trim(gs_wfn_file),form='unformatted')
        read(nfile_gs_wfn) zu_GS1(:,:,1)
     endif
     call comm_bcast(zu_GS1,nproc_group_global)

     do ik=NK_s,NK_e
        zu_GS(:,:,ik) = zu_GS1(:,:,1)
        write (gs_wfn_file,'(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik,'.wfn'
        open(nfile_gs_wfn,file=trim(gs_wfn_file),form='unformatted')
        write(nfile_gs_wfn) zu_GS(:,:,ik)
        close(nfile_gs_wfn)
     end do

     deallocate(zu_GS1)

     if(comm_is_root(nproc_id_global))then
       write(*,*) "  Initial guess was modified:"
       write(*,*) "  Wave function at the first k-point was coppied to the other all."
      endif

    endif

  end subroutine modify_initial_guess_copy_1stk_to_all
end module io_gs_wfn_k

!-----------------------------------------------------------------

module io_rt_wfn_k
  use salmon_global
  use salmon_communication
  use salmon_parallel
  use global_variables
  use misc_routines
  implicit none

  character(256) :: rt_wfn_directory
  character(256) :: rt_wfn_file, occ_file, md_file, ae_file, al_file
  integer,parameter :: nfile_rt_wfn = 41
  integer,parameter :: nfile_occ    = 42
  integer,parameter :: nfile_md     = 43
  integer,parameter :: nfile_ae     = 44
  integer,parameter :: nfile_al     = 45

  integer,parameter :: iflag_read_rt = 0
  integer,parameter :: iflag_write_rt= 1

  character(1) :: alocal_laser_tmp

contains
  subroutine read_write_rt_wfn_k(iflag_read_write)
    implicit none
    integer,intent(in) :: iflag_read_write
    integer :: ik
  ! temporal communicator for multi-scale (same k-points at different macro point)
    integer :: nproc_group_kpoint_ms
    integer :: nproc_id_kpoint_ms
    integer :: nproc_size_kpoint_ms

    write (rt_wfn_directory,'(A,A)') trim(directory),'/rt_wfn_k/'
    if(iflag_read_write == iflag_write_rt)call create_directory(rt_wfn_directory)

    if(comm_is_root(nproc_id_global))then
      occ_file = trim(rt_wfn_directory)//'occupation'
      open(nfile_occ,file=trim(occ_file),form='unformatted')
      select case(iflag_read_write)
      case(iflag_write_rt); write(nfile_occ)occ
      case(iflag_read_rt ); read(nfile_occ )occ
      end select
      close(nfile_occ)

      md_file = trim(rt_wfn_directory)//'Rion_velocity'
      open(nfile_md,file=trim(md_file),form='unformatted')
      select case(iflag_read_write)
      case(iflag_write_rt); write(nfile_md)Rion,velocity
      case(iflag_read_rt ); read(nfile_md )Rion,velocity
      end select
      close(nfile_md)

      ae_file = trim(rt_wfn_directory)//'ae_field'
      open(nfile_ae,file=trim(ae_file),form='unformatted')
      select case(iflag_read_write)
      case(iflag_write_rt)
         t1_delay = t1_delay - Nt*dt
         write(nfile_ae) t1_delay, &
                         trans_longi, &
                         e_impulse, &
                         ae_shape1, ae_shape2, &
                         amplitude1,amplitude2, &
                         pulse_tw1, pulse_tw2, &
                         omega1,    omega2, &
                         epdir_re1, epdir_im1, & 
                         phi_cep1,  phi_cep2, &
                         epdir_re2, epdir_im2, &
                         rlaser_int_wcm2_1, rlaser_int_wcm2_2, &
                         t1_t2, &
                         quadrupole, quadrupole_pot, &
                         rlaserbound_sta, rlaserbound_end, &
                         alocal_laser
      case(iflag_read_rt )
         read(nfile_ae)  t1_delay, &
                         trans_longi, &
                         e_impulse, &
                         ae_shape1, ae_shape2, &
                         amplitude1,amplitude2, &
                         pulse_tw1, pulse_tw2, &
                         omega1,    omega2, &
                         epdir_re1, epdir_im1, &
                         phi_cep1,  phi_cep2, &
                         epdir_re2, epdir_im2, &
                         rlaser_int_wcm2_1, rlaser_int_wcm2_2, &
                         t1_t2, &
                         quadrupole, quadrupole_pot, &
                         rlaserbound_sta, rlaserbound_end, &
                         alocal_laser_tmp
      end select
      close(nfile_ae)

      if(alocal_laser=='y') then
         al_file = trim(rt_wfn_directory)//'alocal_laser_w'
         open(nfile_al,file=trim(al_file),form='unformatted')
         if(.not. allocated(weight_Ac_alocal)) &
         & allocate(weight_Ac_alocal(NL),weight_Ac_alocal_ion(NI), divA_al(NL))
         select case(iflag_read_write)
         case(iflag_write_rt);write(nfile_al)weight_Ac_alocal,weight_Ac_alocal_ion,divA_al
         case(iflag_read_rt );read(nfile_al )weight_Ac_alocal,weight_Ac_alocal_ion,divA_al
         end select
         close(nfile_al)
      endif

    end if

    select case(iflag_read_write)
    case(iflag_read_rt )
      call comm_bcast(occ,     nproc_group_global)
      call comm_bcast(Rion,    nproc_group_global)
      call comm_bcast(velocity,nproc_group_global)

      call comm_bcast(t1_delay,          nproc_group_global)
      call comm_bcast(trans_longi,       nproc_group_global)
      call comm_bcast(e_impulse,         nproc_group_global)
      call comm_bcast(ae_shape1,         nproc_group_global)
      call comm_bcast(ae_shape2,         nproc_group_global)
      call comm_bcast(amplitude1,        nproc_group_global)
      call comm_bcast(amplitude2,        nproc_group_global)
      call comm_bcast(pulse_tw1,         nproc_group_global)
      call comm_bcast(pulse_tw2,         nproc_group_global)
      call comm_bcast(omega1,            nproc_group_global)
      call comm_bcast(omega2,            nproc_group_global)
      call comm_bcast(epdir_re1,         nproc_group_global)
      call comm_bcast(epdir_im1,         nproc_group_global)
      call comm_bcast(phi_cep1,          nproc_group_global)
      call comm_bcast(phi_cep2,          nproc_group_global)
      call comm_bcast(epdir_re2,         nproc_group_global)
      call comm_bcast(epdir_im2,         nproc_group_global)
      call comm_bcast(rlaser_int_wcm2_1, nproc_group_global)
      call comm_bcast(rlaser_int_wcm2_2, nproc_group_global)
      call comm_bcast(t1_t2,             nproc_group_global)
      call comm_bcast(quadrupole,        nproc_group_global)
      call comm_bcast(quadrupole_pot,    nproc_group_global)
      call comm_bcast(rlaserbound_sta,   nproc_group_global)
      call comm_bcast(rlaserbound_end,   nproc_group_global)
      call comm_bcast(alocal_laser,      nproc_group_global)

      if(alocal_laser=='y') then
         if(.not. allocated(weight_Ac_alocal)) then
            allocate(weight_Ac_alocal(NL), weight_Ac_alocal_ion(NI), divA_al(NL))
         endif
         call comm_bcast(weight_Ac_alocal,     nproc_group_global)
         call comm_bcast(weight_Ac_alocal_ion, nproc_group_global)
         call comm_bcast(divA_al,              nproc_group_global)
      endif

    end select

    if(use_ms_maxwell == 'n' .or. (use_ms_maxwell == 'y'.and. nmacro_s == 1))then
      do ik=NK_s,NK_e
        
        write (rt_wfn_file,'(A,A,I7.7,A)') trim(rt_wfn_directory),'/wfn_rt_k',ik,'.wfn'
        open(nfile_rt_wfn,file=trim(rt_wfn_file),form='unformatted')
        select case(iflag_read_write)
        case(iflag_write_rt); write(nfile_rt_wfn)zu_t(:,:,ik)
        case(iflag_read_rt ); read( nfile_rt_wfn)zu_t(:,:,ik)
        end select
        close(nfile_rt_wfn)

      end do

    end if

    if(iflag_read_write==iflag_write_rt)return

    if(use_ms_maxwell == 'y')then
      nproc_group_kpoint_ms = comm_create_group(nproc_group_global, NK_s, nmacro_s)
      call comm_get_groupinfo(nproc_group_kpoint_ms, nproc_id_kpoint_ms, nproc_size_kpoint_ms)
      call comm_bcast(zu_t,nproc_group_kpoint_ms)
    end if

    Rion_eq=Rion
    call psi_rho_RT(zu_t)
    call Hartree
    call Exc_Cor(calc_mode_rt,NBoccmax,zu_t)

    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    call Total_Energy_omp(rion_update_on,calc_mode_rt)
    call Ion_Force_omp(rion_update_on,calc_mode_rt)
    if(PrLv_scf==3 .and. comm_is_root(nproc_id_global)) write(*,*)'Eall =',Eall

  end subroutine read_write_rt_wfn_k

  subroutine read_write_rt_wfn_k_ms_each_macro_grid(iflag_read_write_ms)
    implicit none
    integer,intent(in) :: iflag_read_write_ms
    integer :: ik,imacro
    integer :: nfile_rt_wfn_ms,nfile_occ_ms,nfile_md_ms,nfile_ae_ms,nfile_other_ms
    integer :: ix_m,iy_m,iz_m
    character(1024) :: md_file_ms, ae_file_ms, dir_ae_file_ms, other_file_ms
    real(8),allocatable :: energy_joule_ms_tmp(:,:,:)
    allocate(energy_joule_ms_tmp(nx1_m:nx2_m, ny1_m:ny2_m, nz1_m:nz2_m))

    if (comm_is_root(nproc_id_global)) then
       nfile_ae_ms = 7000
       dir_ae_file_ms = trim(dir_ms)//'rt_ae_field/'
       if(iflag_read_write_ms==iflag_write_rt) call create_directory(dir_ae_file_ms)
       ae_file_ms = trim(dir_ae_file_ms)//'ae_field'
       open(nfile_ae_ms,file=trim(ae_file_ms),form='unformatted')
       select case(iflag_read_write_ms)
       case(iflag_write_rt)
          write(nfile_ae_ms)     Ac_new_ms,    Ac_ms,Jm_new_ms,Jm_ms
       case(iflag_read_rt )
          read(nfile_ae_ms ) add_Ac_new_ms,add_Ac_ms,Jm_new_ms,Jm_ms
       end select
       close(nfile_ae_ms)
    end if

    do imacro = nmacro_s, nmacro_e

       if(comm_is_root(nproc_id_tdks)) then

          write (rt_wfn_directory,'(A,A)') trim(dir_ms_M(imacro)),'/rt_wfn_k/'
          if(iflag_read_write_ms==iflag_write_rt) call create_directory(rt_wfn_directory)

          nfile_occ_ms    = 7000 + imacro
          nfile_rt_wfn_ms = nfile_occ_ms
          nfile_md_ms     = nfile_occ_ms
          nfile_other_ms  = nfile_occ_ms

          occ_file = trim(rt_wfn_directory)//'occupation'
          open(nfile_occ_ms,file=trim(occ_file),form='unformatted')
          select case(iflag_read_write_ms)
          case(iflag_write_rt); write(nfile_occ_ms) occ
          case(iflag_read_rt ); read(nfile_occ_ms ) occ
          end select
          close(nfile_occ_ms)

          md_file_ms = trim(rt_wfn_directory)//'Rion_velocity'
          open(nfile_md_ms,file=trim(md_file_ms),form='unformatted')
          select case(iflag_read_write_ms)
          case(iflag_write_rt); write(nfile_md_ms) Rion,velocity
          case(iflag_read_rt ); read(nfile_md_ms ) Rion,velocity
          end select
          close(nfile_md)

          other_file_ms = trim(rt_wfn_directory)//'others'
          open(nfile_other_ms,file=trim(other_file_ms),form='unformatted')
          ix_m = macropoint(1,imacro)
          iy_m = macropoint(2,imacro)
          iz_m = macropoint(3,imacro)
          select case(iflag_read_write_ms)
          case(iflag_write_rt); write(nfile_other_ms) Eall0_m(imacro), energy_joule_ms(ix_m,iy_m,iz_m)
          case(iflag_read_rt ); read(nfile_other_ms ) Eall0_m(imacro), energy_joule_ms_tmp(ix_m,iy_m,iz_m)
          end select
          close(nfile_other_ms)

          !! maybe must be out of "if(comm_is_root(nproc_id_tdks)) then"
          do ik=NK_s,NK_e
             write(rt_wfn_file,'(A,A,I7.7,A)') trim(rt_wfn_directory),'/wfn_rt_k',ik,'.wfn'
             open(nfile_rt_wfn_ms,file=trim(rt_wfn_file),form='unformatted')
             select case(iflag_read_write_ms)
             case(iflag_write_rt); write(nfile_rt_wfn_ms) zu_m(:,:,ik,imacro)
             case(iflag_read_rt ); read( nfile_rt_wfn_ms) zu_t(:,:,ik)
             end select
             close(nfile_rt_wfn_ms)
          end do

       end if

       if(iflag_read_write_ms == iflag_read_rt) then

          call comm_bcast(occ,     nproc_group_tdks)
          call comm_bcast(Rion,    nproc_group_tdks)
          call comm_bcast(velocity,nproc_group_tdks)
          Rion_eq(:,:)           = Rion(:,:)
          Rion_m(:,:,imacro)     = Rion(:,:)
          Rion_eq_m(:,:,imacro)  = Rion_m(:,:,imacro)
          velocity_m(:,:,imacro) = velocity(:,:)

          call comm_summation(energy_joule_ms_tmp,energy_joule_ms,&
                             &(nx2_m-nx1_m+1)*(ny2_m-ny1_m+1)*(nz2_m-nz1_m+1),nproc_group_global)
          call comm_bcast(add_Ac_ms,     nproc_group_global)
          call comm_bcast(add_Ac_new_ms, nproc_group_global)
          call comm_bcast(Jm_new_ms,     nproc_group_global)
          call comm_bcast(Jm_ms,         nproc_group_global)

          call psi_rho_RT(zu_t)
          call Hartree
          call Exc_Cor(calc_mode_rt,NBoccmax,zu_t)

          Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
          call Total_Energy_omp(rion_update_on,calc_mode_rt,imacro)
          call Ion_Force_omp(rion_update_on,calc_mode_rt,imacro)

       endif

    enddo  !imacro

  end subroutine read_write_rt_wfn_k_ms_each_macro_grid

end module io_rt_wfn_k
