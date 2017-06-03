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


subroutine arted(x_nprocs, x_myrank)
  use salmon_global, only:    use_ms_maxwell
  use Global_Variables, only: calc_mode, &
                            & calc_mode_sc, &
                            & calc_mode_ms
  use communication,    only: proc_group, &
                            & nprocs, &
                            & procid
  use control_sc,       only: main_sc => main
  use control_ms,       only: main_ms => main
  use inputfile,        only: read_arted => read_input, &
                            & dump_inputdata
  use mpi,              only: MPI_COMM_WORLD
  
  implicit none
  integer, intent(in) :: x_myrank
  integer, intent(in) :: x_nprocs
  
  proc_group(:) = MPI_COMM_WORLD
  nprocs(:) = x_nprocs
  procid(:) = x_myrank

  call read_arted()
  !call dump_inputdata

  select case(use_ms_maxwell)
  case ('y')
    call main_sc
  case ('n')
    call main_ms
  case default
    call Err_finalize("Invalid use_ms_maxwell parameter!")
  end select
  
  return
end subroutine arted
