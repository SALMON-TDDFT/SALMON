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
!This file is "main.f90"
!This file contains mainroutine program
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

Subroutine err_finalize(err_message)
  use Global_Variables
  use salmon_parallel, only: nproc_id_global, end_parallel
  use salmon_communication, only: comm_is_root
  implicit none
  character(*),intent(in) :: err_message
  if (comm_is_root(nproc_id_global)) then
    write(*,*) err_message
  endif
  call end_parallel

  stop
End Subroutine Err_finalize


subroutine arted
  use salmon_global, only:    use_ms_maxwell, restart_option
  use control_sc,       only: tddft_sc
  use control_ms,       only: tddft_maxwell_ms


  use salmon_parallel
  use initialization
  use ground_state
  use io_gs_wfn_k
  
  implicit none

  nproc_group_tdks = nproc_group_global
  nproc_id_tdks    = nproc_id_global
  nproc_size_tdks  = nproc_size_global

  call initialize

  if(restart_option == 'new')then
    select case(iflag_calc_mode)
    case(iflag_calc_mode_gs_rt)
      call calc_ground_state
    case(iflag_calc_mode_gs)
      call calc_ground_state
      call read_write_gs_wfn_k(iflag_write)
      return
    case(iflag_calc_mode_rt)
      call read_write_gs_wfn_k(iflag_read)
    end select
  else if(restart_option == 'restart')then
  else
    call Err_finalize("Invalid restart_option!")
  end if

  select case(use_ms_maxwell)
  case ('y')
    call tddft_maxwell_ms
  case ('n')
    call tddft_sc
  case default
    call Err_finalize("Invalid use_ms_maxwell parameter!")
  end select
  
  return
end subroutine arted
