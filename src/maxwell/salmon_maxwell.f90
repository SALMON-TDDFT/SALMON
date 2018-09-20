!
!  Copyright 2018 SALMON developers
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
!-----------------------------------------------------------------------------------------
module salmon_maxwell
  implicit none
  
  type fdtd_grid
    integer :: ng_sta(3), ng_end(3)          ! Size of Local Grid System
    integer :: lg_sta(3), lg_end(3)          ! Size of Global Grid System
    integer :: no_sta(3), no_end(3)          ! Size of Entire (Allocated) Variables
    real(8) :: dt                            ! Delta t
    integer :: iter_now                      ! Present iteration Number
    real(8) :: hgs(3)                        ! Grid Spacing
    real(8) :: origin(3)                     ! Coordinate of Origin Point (TBA)
    integer :: i_bc(3, 2)                    ! Boundary Condition for 1:x, 2:y, 3:z and 1:bottom and 2:top
    character(16) :: gauge                   ! Gauge Condition (TBD)
    integer, allocatable :: imaterial(:,:,:) ! Material information
  end type fdtd_grid

  type fdtd_field
    real(8), allocatable :: vec_e(:,:,:,:)             ! E
    real(8), allocatable :: vec_h(:,:,:,:)             ! H
    real(8), allocatable :: vec_a(:,:,:,:), phi(:,:,:) ! Vector and Scalar Potential
    real(8), allocatable :: vec_j_em(:,:,:,:)          ! Electromagnetic Current
    real(8), allocatable :: rho_em(:,:,:,:)            ! Electromagnetic Charge
  end type fdtd_field
  
  type fdtd_tmp
    !weyl
    
    !coulomb
    
    !eh
    integer :: pml_l !pml parameter
    real(8) :: pml_m !pml parameter
    real(8) :: pml_r !pml parameter
    real(8),allocatable :: ex_y(:,:,:),c1_ex_y(:,:,:),c2_ex_y(:,:,:),ex_z(:,:,:),c1_ex_z(:,:,:),c2_ex_z(:,:,:) !e
    real(8),allocatable :: ey_z(:,:,:),c1_ey_z(:,:,:),c2_ey_z(:,:,:),ey_x(:,:,:),c1_ey_x(:,:,:),c2_ey_x(:,:,:) !e
    real(8),allocatable :: ez_x(:,:,:),c1_ez_x(:,:,:),c2_ez_x(:,:,:),ez_y(:,:,:),c1_ez_y(:,:,:),c2_ez_y(:,:,:) !e
    real(8),allocatable :: hx_y(:,:,:),c1_hx_y(:,:,:),c2_hx_y(:,:,:),hx_z(:,:,:),c1_hx_z(:,:,:),c2_hx_z(:,:,:) !h
    real(8),allocatable :: hy_z(:,:,:),c1_hy_z(:,:,:),c2_hy_z(:,:,:),hy_x(:,:,:),c1_hy_x(:,:,:),c2_hy_x(:,:,:) !h
    real(8),allocatable :: hz_x(:,:,:),c1_hz_x(:,:,:),c2_hz_x(:,:,:),hz_y(:,:,:),c1_hz_y(:,:,:),c2_hz_y(:,:,:) !h
    real(8),allocatable :: hx_m(:,:,:),hy_m(:,:,:),hz_m(:,:,:)                                                 !h
  end type fdtd_tmp
  
  contains
  
  subroutine init_maxwell(grid,field,tmp)
    use inputoutput, only: theory, use_ms_maxwell
    implicit none
    type(fdtd_grid)  :: grid
    type(fdtd_field) :: field
    type(fdtd_tmp)   :: tmp
    
    select case(theory)
    case('Maxwell+TDDFT')
      !this selection is temporary. 
      !After removing use_ms_maxwell, this selection is revised.
      select case(use_ms_maxwell) 
      case('y')
        grid%gauge = 'weyl' 
      case('s')
        grid%gauge = 'coulomb' 
      end select
    case('Maxwell')
      grid%gauge = 'eh' 
    case default
      stop 'invalid theory'
    end select
    
    select case(grid%gauge)
    case('weyl')
      call init_weyl(grid,field,tmp)
    case('coulomb')
      call init_coulomb(grid,field,tmp)
    case('eh')
      call init_eh(grid,field,tmp)
    end select
  end subroutine init_maxwell
  
  subroutine finalize_maxwell(grid,field,tmp)
    implicit none
    type(fdtd_grid)  :: grid
    type(fdtd_field) :: field
    type(fdtd_tmp)   :: tmp
    
    select case(grid%gauge)
    case('weyl')
      call finalize_weyl(grid,field,tmp)
    case('coulomb')
      call finalize_coulomb(grid,field,tmp)
    case('eh')
      call finalize_eh(grid,field,tmp)
    end select
  end subroutine finalize_maxwell
  
  subroutine calc_maxwell(grid,field,tmp)
    implicit none
    type(fdtd_grid)  :: grid
    type(fdtd_field) :: field
    type(fdtd_tmp)   :: tmp
    
    select case(grid%gauge)
    case('weyl')
      call calc_weyl(grid,field,tmp)
    case('coulomb')
      call calc_coulomb(grid,field,tmp)
    case('eh')
      call calc_eh(grid,field,tmp)
    end select
  end subroutine calc_maxwell
  
end module salmon_maxwell
