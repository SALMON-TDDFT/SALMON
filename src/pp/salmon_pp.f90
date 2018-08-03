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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module salmon_pp

  implicit none

  type pp_info
    real(8) :: zion
    integer :: lmax,lmax0
    integer :: nrmax,nrmax0
    logical :: flag_nlcc
    real(8),allocatable :: rmass(:)
    integer,allocatable :: mr(:)
    integer,allocatable :: lref(:)
    integer,allocatable :: nrps(:)
    integer,allocatable :: mlps(:)
    integer,allocatable :: zps(:)
    integer,allocatable :: nrloc(:)
    real(8),allocatable :: rloc(:)
    real(8),allocatable :: rps(:)
    real(8),allocatable :: anorm(:,:)
    integer,allocatable :: inorm(:,:)
    real(8),allocatable :: rad(:,:)
    real(8),allocatable :: radnl(:,:)
    real(8),allocatable :: vloctbl(:,:)
    real(8),allocatable :: dvloctbl(:,:)
    real(8),allocatable :: udvtbl(:,:,:)
    real(8),allocatable :: dudvtbl(:,:,:)
    real(8),allocatable :: rho_nlcc_tbl(:,:)
    real(8),allocatable :: tau_nlcc_tbl(:,:)
    real(8),allocatable :: upp_f(:,:,:)
    real(8),allocatable :: vpp_f(:,:,:)
    real(8),allocatable :: upp(:,:)
    real(8),allocatable :: dupp(:,:)
    real(8),allocatable :: vpp(:,:)
    real(8),allocatable :: dvpp(:,:)
  end type

  contains

  subroutine init_pp(pp,nrmax,lmax,flag_nlcc)
    use salmon_global, only: nelem,lloc_ps
    implicit none
    type(pp_info) :: pp
    integer, parameter :: nrmax0=50000, lmax0=4
    integer,intent(in) :: nrmax,lmax
    logical,intent(in) :: flag_nlcc

    pp%lmax0=lmax0

    pp%nrmax0=nrmax0
  
    pp%nrmax=nrmax
    pp%lmax=lmax
  
    allocate(pp%rmass(1:nelem))
    
    allocate(pp%mr(1:nelem))

    allocate(pp%lref(1:nelem))
    pp%lref(1:nelem)=lloc_ps(1:nelem)

    allocate(pp%nrps(1:nelem))
    allocate(pp%rps(1:nelem))
    allocate(pp%mlps(1:nelem))
    allocate(pp%zps(1:nelem))
    allocate(pp%nrloc(1:nelem))
    allocate(pp%rloc(1:nelem))
    
    allocate(pp%anorm(0:lmax,nelem))
    allocate(pp%inorm(0:lmax,nelem))
  
    allocate(pp%rad(nrmax,nelem))
    allocate(pp%radnl(nrmax,nelem))
    
    allocate(pp%vloctbl(nrmax,nelem))
    allocate(pp%dvloctbl(nrmax,nelem))
    allocate(pp%udvtbl(nrmax,0:lmax,nelem))
    allocate(pp%dudvtbl(nrmax,0:lmax,nelem))
    
    allocate(pp%rho_nlcc_tbl(nrmax,nelem))
    allocate(pp%tau_nlcc_tbl(nrmax,nelem))
  
    allocate(pp%vpp(0:nrmax0,0:lmax0),pp%upp(0:nrmax0,0:lmax0))
    allocate(pp%dvpp(0:nrmax0,0:lmax0),pp%dupp(0:nrmax0,0:lmax0))
    allocate(pp%vpp_f(0:nrmax0,0:lmax0,nelem),pp%upp_f(0:nrmax0,0:lmax0,nelem))
  
    pp%flag_nlcc=flag_nlcc
  
  end subroutine init_pp

end module salmon_pp
