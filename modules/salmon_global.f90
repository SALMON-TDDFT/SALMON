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
module salmon_global
  implicit none

!Parameters for pseudo-potential
  integer, parameter :: maxMKI=10
  integer :: MI,MKI
   !shinohara
  integer :: ipsfileform(maxMKI)   ! file format for pseudo potential
  character(16)  :: ps_format(maxMKI)
! List of pseudopotential file formats
  integer,parameter :: n_Yabana_Bertsch_psformat = 1 !.rps
  integer,parameter :: n_ABINIT_psformat = 2 ! .pspnc
  integer,parameter :: n_FHI_psformat = 3 ! .cpi
  integer,parameter :: n_ABINITFHI_psformat = 4 ! .fhi


!Input variables
!! &calculation
  character(16)  :: calc_mode
  character(1)   :: use_ehrenfest_md
  character(1)   :: use_ms_maxwell
  character(1)   :: use_force
  character(1)   :: use_geometry_opt

!! &control      
  character(8)   :: restart_option
  integer        :: backup_frequency
  real(8)        :: time_shutdown
  character(256) :: sysname
  character(256) :: directory
                 
!! &units
  character(16)  :: unit_time
  character(16)  :: unit_length
  character(16)  :: unit_energy
  character(16)  :: unit_charge
                 
!! &parallel
  character(1)   :: domain_parallel
  integer        :: nproc_ob
  integer        :: nproc_domain(3)
  integer        :: nproc_domain_s(3)
  integer        :: num_datafiles_in
  integer        :: num_datafiles_out

!! &system
  integer        :: iperiodic
  integer        :: ispin
  real(8)        :: al(3)
  integer        :: isym
  character(32)  :: crystal_structure
  integer        :: nstate
  integer        :: nelec
  real(8)        :: temperature
  integer        :: nelem
  integer        :: natom
  character(256) :: file_atom

!! &pseudo
  character(256) :: pseudo_file(maxMKI)
  integer        :: Lmax_ps(maxMKI)
  integer        :: Lloc_ps(maxMKI)
  integer        :: iZatom(maxMKI)
  character(1)   :: psmask_option
  real(8)        :: alpha_mask
  real(8)        :: gamma_mask
  real(8)        :: eta_mask

!! &functional
  character(32)  :: xc
  real(8)        :: cval

!! &rgrid
  real(8)        :: dl(3)
  integer        :: num_rgrid(3)

!! &kgrid
  integer        :: num_kgrid(3)
  character(256) :: file_kw

!! &tgrid
  integer        :: nt
  real(8)        :: dt

!! &propagation
  integer        :: n_hamil
  character(16)  :: propagator

!! &scf
  integer        :: ncg
  integer        :: nmemory_mb
  real(8)        :: alpha_mb
  character(1)   :: fsset_option
  integer        :: nfsset_start
  integer        :: nfsset_every
  integer        :: nscf
  integer        :: ngeometry_opt
  character(1)   :: subspace_diagonalization
  character(16)  :: cmixing
  real(8)        :: rmixrate
  character(3)   :: convergence
  real(8)        :: threshold

!! &emfield
  character(2)   :: trans_longi
  character(16)  :: ae_shape1
  real(8)        :: amplitude1
  real(8)        :: rlaser_int1
  real(8)        :: pulse_tw1
  real(8)        :: omega1
  real(8)        :: epdir_re1(3)
  real(8)        :: epdir_im1(3)
  real(8)        :: phi_cep1
  character(16)  :: ae_shape2
  real(8)        :: amplitude2
  real(8)        :: rlaser_int2
  real(8)        :: pulse_tw2
  real(8)        :: omega2
  real(8)        :: epdir_re2(3)
  real(8)        :: epdir_im2(3)
  real(8)        :: phi_cep2
  real(8)        :: t1_t2
  character(1)   :: quadrupole
  character(8)   :: quadrupole_pot

!! &linear_response
  real(8)        :: e_impulse

!! &multiscale
  character(16)  :: fdtddim
  character(16)  :: twod_shape
  integer        :: nx_m
  integer        :: ny_m
  integer        :: nz_m
  real(8)        :: hx_m
  real(8)        :: hy_m
  real(8)        :: hz_m
  integer        :: nksplit
  integer        :: nxysplit
  integer        :: nxvacl_m
  integer        :: nxvacr_m

!! &analysis
  character(2)   :: projection_option
  integer        :: nenergy
  real(8)        :: de
  character(1)   :: out_psi
  character(1)   :: out_dos
  character(1)   :: out_pdos
  character(1)   :: out_dns
  character(1)   :: out_elf
  character(1)   :: out_dns_rt
  integer        :: out_dns_rt_step
  character(1)   :: out_elf_rt
  integer        :: out_elf_rt_step
  character(1)   :: out_estatic_rt
  integer        :: out_estatic_rt_step
  character(16)  :: format3d

!! &hartree
  integer        :: meo
  integer        :: num_pole_xyz(3)

!! &ewald
  integer        :: newald
  real(8)        :: aewald



end module salmon_global
