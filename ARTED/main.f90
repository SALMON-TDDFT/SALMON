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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130

Subroutine err_finalize(err_message)
  use Global_Variables
  use communication
  implicit none
  character(*),intent(in) :: err_message
  if (comm_is_root()) then
    write(*,*) err_message
  endif
  call comm_finalize

  stop
End Subroutine Err_finalize


Program main
  use Global_Variables, only: calc_mode, &
                            & calc_mode_sc, &
                            & calc_mode_ms
  use communication,    only: comm_init, &
                            & comm_is_root, &
                            & comm_bcast_character, &
                            & proc_group              
  use control_sc,       only: main_sc => main
  use control_ms,       only: main_ms => main
  use inputfile,        only: read_input, &
                            & dump_inputdata
  implicit none
  
  call comm_init()
  if (comm_is_root()) read(*,*) calc_mode
  call comm_bcast_character(calc_mode, proc_group(1))
  
  call read_input
  !call dump_inputdata

  select case(calc_mode)
  case (calc_mode_sc)
    call main_sc
  case (calc_mode_ms)
    call main_ms
  case default
    call Err_finalize("Invalid calc_mode parameter!")
  end select
  
  stop
End Program main
