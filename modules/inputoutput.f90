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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
module inputoutput
  use salmon_global
  implicit none
!Physical constant
  real(8),parameter :: au_time_fs = 0.02418884326505d0
  real(8),parameter :: au_energy_ev = 27.21138505d0
  real(8),parameter :: au_length_aa = 0.52917721067d0


  integer, parameter :: fh_namelist = 901
  integer, parameter :: fh_atomic_spiecies = 902
  integer, parameter :: fh_atomic_positions = 903
  integer, parameter :: fh_reentrance = 904

  integer :: inml_group_function
  integer :: inml_units
  integer :: inml_control
  integer :: inml_system
  integer :: inml_incident
  integer :: inml_propagation
  integer :: inml_rgrid
  integer :: inml_kgrid
  integer :: inml_tstep
  integer :: inml_electrons
  integer :: inml_pseudo
  integer :: inml_response
  integer :: inml_multiscale
  integer :: inml_group_atom

!Input/Output units
  real(8) :: utime_to_au, utime_from_au
  real(8) :: ulength_to_au, ulength_from_au
  real(8) :: uenergy_to_au, uenergy_from_au
  real(8) :: ucharge_to_au, ucharge_from_au

contains


  subroutine read_stdin(myrank)
    implicit none
    include 'mpif.h'
    integer,intent(in) :: myrank
    integer :: ierr

    integer :: cur = fh_namelist
    integer :: ret = 0
    character(100) :: buff, text


    
    if (myrank == 0) then
      open(fh_namelist, file='.namelist.tmp', status='replace')
!      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')
      open(fh_atomic_positions, file='.atomic_positions.tmp', status='replace')
      open(fh_reentrance, file='.reenetrance.tmp', status='replace')
      
      do while (.true.)
        read(*, '(a)', iostat=ret) buff
        if (ret < 0) then
          exit
        else
          text = trim(adjustl(buff))
          ! Comment lines
          if (text(1:1) == '!') cycle
!          ! Beginning of 'atomic_species' part
!          if (text == '&atomic_spiecies') then
!            cur = fh_atomic_spiecies
!            cycle
!          end if
          ! Beginning of 'atomic_positions' part
          if (text == '&atomic_positions') then
            cur = fh_atomic_positions
            cycle
          end if
          ! Beginning of 'atomic_species' part
          if (text == '&reentrance') then
            cur = fh_reentrance
            cycle
          end if
          ! End of 'atomic_(spiecies|positions)' part
          if ((text == '/') .and. (cur /= fh_namelist)) then
            cur = fh_namelist
            cycle
          end if
          
          write(cur, '(a)') text
        end if
      end do
      close(fh_namelist)
      close(fh_atomic_positions)
!      close(fh_atomic_spiecies)
      close(fh_reentrance)
    end if

!    call comm_sync_all()


    return
  end subroutine read_stdin

  subroutine read_input_common(myrank,cfunction)
    implicit none
    include 'mpif.h'
    integer,intent(in) :: myrank
    character(30),intent(out) :: cfunction
    integer :: ierr


    namelist/group_function/ cfunction
    namelist/units/ &
            & unit_time, &
            & unit_length, &
            & unit_energy


    if (myrank == 0)then
      open(fh_namelist, file='.namelist.tmp', status='old')

      read(fh_namelist, nml=group_function, iostat=inml_group_function)
      rewind(fh_namelist)


      unit_time='au'
      unit_length='au'
      unit_energy='au'
      unit_charge='au'
      read(fh_namelist, nml=units, iostat=inml_units)
      rewind(fh_namelist)

      close(fh_namelist)
    end if


    call mpi_bcast(cfunction,30,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_time,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_length,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_energy,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_charge,16,mpi_character,0,mpi_comm_world,ierr)


    call initialize_inputoutput_units

  end subroutine read_input_common

  subroutine initialize_inputoutput_units
    implicit none


! Unit for time
    select case(unit_time)
    case('au','a.u.')
      utime_to_au   = 1d0
      utime_from_au = 1d0
    case('fs','femtosecond')
      utime_to_au   = 1d0/au_time_fs
      utime_from_au = au_time_fs
    case default
      stop "Invalid unit for time."
    end select

! Unit for length
    select case(unit_length)
    case('au','a.u.')
      ulength_to_au   = 1d0
      ulength_from_au = 1d0
    case('AA','angstrom','Angstrom')
      ulength_to_au   = 1d0/au_length_aa
      ulength_from_au = au_length_aa
    case default
      stop "Invalid unit for length."
    end select

! Unit for energy
    select case(unit_energy)
    case('au','a.u.')
      uenergy_to_au   = 1d0
      uenergy_from_au = 1d0
    case('ev','eV')
      uenergy_to_au   = 1d0/au_energy_ev
      uenergy_from_au = au_energy_ev
    case default
      stop "Invalid unit for energy."
    end select

! Unit for charge
    select case(unit_charge)
    case('au','a.u.')
      ucharge_to_au   = 1d0
      ucharge_from_au = 1d0
    case default
      stop "Invalid unit for charge."
    end select

  end subroutine initialize_inputoutput_units
    
end module inputoutput



