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
  use salmon_global, only:    use_ms_maxwell
  use Global_Variables, only: calc_mode, &
                            & calc_mode_sc, &
                            & calc_mode_ms
  use control_sc,       only: main_sc => main
  use control_ms,       only: main_ms => main
  use inputfile,        only: read_arted => read_input, &
                            & dump_inputdata
  use salmon_parallel
  
  implicit none

!  nproc_group_maxwell = nproc_group_global
  nproc_group_tdks    = nproc_group_global
!  nproc_id_maxwell    = nproc_id_global
  nproc_id_tdks       = nproc_id_global
!  nproc_size_maxwell  = nproc_size_global
  nproc_size_tdks     = nproc_size_global

  call read_arted()
  !call dump_inputdata

  select case(use_ms_maxwell)
  case ('y')
    call main_ms
  case ('n')
    call main_sc
  case default
    call Err_finalize("Invalid use_ms_maxwell parameter!")
  end select
  
  return
end subroutine arted
