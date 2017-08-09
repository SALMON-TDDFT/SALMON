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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
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

! Flag for atomic coordinate type
  integer :: iflag_atom_coor
  integer,parameter :: ntype_atom_coor_none      = 0
  integer,parameter :: ntype_atom_coor_cartesian = 1
  integer,parameter :: ntype_atom_coor_reduced   = 2

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
  character(256) :: dump_filename
                 
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
  integer        :: nstate_spin(2)
  integer        :: nelec
  integer        :: nelec_spin(2)
  real(8)        :: temperature
  integer        :: nelem
  integer        :: natom
  character(256) :: file_atom_coor
  character(256) :: file_atom_red_coor

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
  character(8)   :: amin_routine
  integer        :: ncg
  character(8)   :: amixing
  real(8)        :: rmixrate
  integer        :: nmemory_mb
  real(8)        :: alpha_mb
  character(1)   :: fsset_option
  integer        :: nfsset_start
  integer        :: nfsset_every
  integer        :: nscf
  integer        :: ngeometry_opt
  character(1)   :: subspace_diagonalization
  character(8)   :: convergence
  real(8)        :: threshold
  real(8)        :: threshold_pot

!! &emfield
  character(2)   :: trans_longi
  character(16)  :: ae_shape1
  real(8)        :: e_impulse
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
  character(1)   :: alocal_laser
  real(8)        :: rlaserbound_sta(3)
  real(8)        :: rlaserbound_end(3)

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
  real(8)        :: out_dos_start
  real(8)        :: out_dos_end
  integer        :: iout_dos_nenergy
  real(8)        :: out_dos_smearing
  character(16)  :: out_dos_method
  character(1)   :: out_dos_fshift
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
  integer        :: numfiles_out_3d
  character(1)   :: timer_process

!! &hartree
  integer        :: meo
  integer        :: num_pole_xyz(3)

!! &ewald
  integer        :: newald
  real(8)        :: aewald

!! &group_fundamental
  integer        :: iditerybcg
  integer        :: iditer_nosubspace_diag
  integer        :: ntmg
  integer        :: idisnum(2)
  integer        :: iwrite_projection
  integer        :: itwproj
  integer        :: iwrite_projnum
  integer        :: itcalc_ene

!! &group_parallel
  integer        :: isequential
  integer        :: imesh_s_all
  integer        :: iflag_comm_rho

!! &group_hartree
  real(8)        :: hconv
  integer        :: lmax_meo

!! &group_file
  integer        :: ic
  integer        :: oc
  integer        :: ic_rt
  integer        :: oc_rt

!! &group_others
  integer        :: iparaway_ob
  integer        :: iscf_order
  integer        :: iswitch_orbital_mesh
  integer        :: iflag_psicube
  real(8)        :: lambda1_diis
  real(8)        :: lambda2_diis
  character(100) :: file_ini
  integer        :: num_projection
  integer        :: iwrite_projection_ob(200)
  integer        :: iwrite_projection_k(200)
  character(100) :: filename_pot
  integer        :: iwrite_external
  integer        :: iflag_dip2
  integer        :: iflag_intelectron
  integer        :: num_dip2
  real(8)        :: dip2boundary(100)
  real(8)        :: dip2center(100)
  integer        :: iflag_fourier_omega
  integer        :: num_fourier_omega
  real(8)        :: fourier_omega(200)
  integer        :: itotntime2
  integer        :: iwdenoption
  integer        :: iwdenstep
  integer        :: iflag_estatic

!! &atomic_coor
!! &atomic_red_coor
integer,allocatable :: Kion(:)    
real(8),allocatable :: Rion(:,:)  
real(8),allocatable :: Rion_red(:,:)  
character(1),allocatable :: flag_geo_opt_atom(:)


end module salmon_global
