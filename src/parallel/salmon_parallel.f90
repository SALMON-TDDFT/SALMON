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
module salmon_parallel
  implicit none

  !!! multi process
  ! global
  integer, public :: nproc_group_global
  integer, public :: nproc_id_global
  integer, public :: nproc_size_global

  ! TDKS eq.
  integer, public :: nproc_group_tdks
  integer, public :: nproc_id_tdks
  integer, public :: nproc_size_tdks

  integer, public :: nproc_group_spin
  integer, public :: nproc_id_spin
  integer, public :: nproc_size_spin

  integer, public :: nproc_group_kgrid
  integer, public :: nproc_id_kgrid
  integer, public :: nproc_size_kgrid

  integer, public :: nproc_group_korbital
  integer, public :: nproc_id_korbital
  integer, public :: nproc_size_korbital

  integer, public :: nproc_group_rho
  integer, public :: nproc_id_rho
  integer, public :: nproc_size_rho

  integer, public :: nproc_group_k
  integer, public :: nproc_id_k
  integer, public :: nproc_size_k

  integer, public :: nproc_group_grid
  integer, public :: nproc_id_grid
  integer, public :: nproc_size_grid

  integer, public :: nproc_group_orbitalgrid
  integer, public :: nproc_id_orbitalgrid
  integer, public :: nproc_size_orbitalgrid

  integer, public :: nproc_group_h
  integer, public :: nproc_id_h
  integer, public :: nproc_size_h

  integer, public :: nproc_group_korbital_vhxc
  integer, public :: nproc_id_korbital_vhxc
  integer, public :: nproc_size_korbital_vhxc

  integer, public :: nproc_group_bound(3)
  integer, public :: nproc_id_bound(3)
  integer, public :: nproc_size_bound(3)

  ! FFTE
  integer, public :: nproc_group_icommy
  integer, public :: nproc_id_icommy
  integer, public :: nproc_size_icommy

  integer, public :: nproc_group_icommz
  integer, public :: nproc_id_icommz
  integer, public :: nproc_size_icommz

  integer, public :: nproc_group_icommw
  integer, public :: nproc_id_icommw
  integer, public :: nproc_size_icommw

  integer, public :: nproc_group_icommy_copy
  integer, public :: nproc_id_icommy_copy
  integer, public :: nproc_size_icommy_copy

  ! call once
  public :: setup_parallel
  public :: end_parallel

  ! util
  public :: get_thread_id
  public :: is_distributed_parallel

contains
  subroutine setup_parallel
    use salmon_communication
    implicit none
    call comm_init
    call comm_get_globalinfo(nproc_group_global, nproc_id_global, nproc_size_global)
  end subroutine

  subroutine end_parallel
    use salmon_communication
    implicit none
    call comm_finalize
  end subroutine

  function get_thread_id() result(nid)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: nid
#ifdef _OPENMP
    nid = omp_get_thread_num()
#else
    nid = 0
#endif
  end function

  function is_distributed_parallel() result(ret)
    implicit none
    logical :: ret
    ret = (nproc_size_global > 1)
  end function
end module
