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
module input

  implicit none
  integer, parameter :: fh_namelist = 901
  integer, parameter :: fh_atomic_spiecies = 902
  integer, parameter :: fh_atomic_positions = 903
  integer, parameter :: fh_reentrance = 904

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



contains


  subroutine read_stdin(myrank,cfunction)
    implicit none
    include 'mpif.h'
    character(30),intent(out) :: cfunction
    integer :: ierr
    integer :: myrank
    namelist/group_function/ cfunction

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

    if (myrank == 0) then
       open(fh_namelist, file='.namelist.tmp', status='old')
       read(fh_namelist, nml=group_function)
       close(fh_namelist)
    end if
    call mpi_bcast(cfunction,30,mpi_character,0,mpi_comm_world,ierr)

    return
  end subroutine read_stdin

end module input



