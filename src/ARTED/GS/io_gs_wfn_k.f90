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


    if(use_ms_maxwell == 'n' .or. (use_ms_maxwell == 'y'.and. NXY_s == 0))then
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
      nproc_group_kpoint_ms = comm_create_group(nproc_group_global, NK_s, NXY_s)
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
    if(PrLv_scf==3 .and. comm_is_root(nproc_id_global)) write(*,*)'Eall =',Eall


  end subroutine read_write_gs_wfn_k
end module io_gs_wfn_k

