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
  integer,parameter :: nfile_occ = 42

  integer,parameter :: iflag_read = 0
  integer,parameter :: iflag_write = 1

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
    Eall_GS0=Eall
    if(PrLv_scf==3 .and. comm_is_root(nproc_id_global)) write(*,*)'Eall =',Eall


  end subroutine read_write_gs_wfn_k

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
  character(256) :: rt_wfn_file, occ_file, md_file, ae_file
  integer,parameter :: nfile_rt_wfn = 41
  integer,parameter :: nfile_occ    = 42
  integer,parameter :: nfile_md     = 43
  integer,parameter :: nfile_ae     = 44

  integer,parameter :: iflag_read_rt = 0
  integer,parameter :: iflag_write_rt= 1

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
                         alocal_laser
      end select
      close(nfile_ae)
    end if

    select case(iflag_read_write)
    case(iflag_read_rt )
      call comm_bcast(occ,     nproc_group_global)
      call comm_bcast(Rion,    nproc_group_global)
      call comm_bcast(velocity,nproc_group_global)
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

end module io_rt_wfn_k
