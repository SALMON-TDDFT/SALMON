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
  integer, parameter :: maxmki=10
  integer :: mi,mki
   !shinohara
  integer :: ipsfileform(maxmki)   ! file format for pseudo potential
  character(16)  :: ps_format(maxmki)
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
  character(16)  :: theory
  character(16)  :: calc_mode
  character(1)   :: use_ehrenfest_md
  character(1)   :: use_adiabatic_md
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
  character(20)  :: modify_gs_wfn_k
  character(1)   :: read_gs_wfn_k
  character(1)   :: read_rt_wfn_k
  character(1)   :: write_gs_wfn_k
  character(1)   :: write_rt_wfn_k
  character(1)   :: read_gs_wfn_k_ms
  character(1)   :: read_rt_wfn_k_ms
  character(1)   :: write_gs_wfn_k_ms
  character(1)   :: write_rt_wfn_k_ms

!! &units
  character(16)  :: unit_system
  character(16)  :: unit_time
  character(16)  :: unit_length
  character(16)  :: unit_energy
  character(16)  :: unit_charge
                 
!! &parallel
  character(1)   :: domain_parallel
  integer        :: nproc_k
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
  real(8)        :: temperature_k
  integer        :: nelem
  integer        :: natom
  character(256) :: file_atom_coor
  character(256) :: file_atom_red_coor

!! &pseudo
  character(256) :: pseudo_file(maxmki)
  integer        :: lmax_ps(maxmki)
  integer        :: lloc_ps(maxmki)
  integer        :: izatom(maxmki)
  character(1)   :: psmask_option
  real(8)        :: alpha_mask
  real(8)        :: gamma_mask
  real(8)        :: eta_mask

!! &functional
  character(64)  :: xc !, xcname
  character(64)  :: xname
  character(64)  :: cname
  character(64)  :: alibx
  character(64)  :: alibc
  character(64)  :: alibxc
  real(8)        :: cval
  character(1)   :: no_update_func

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
  character(16)  :: convergence
  real(8)        :: threshold
  real(8)        :: threshold_norm_rho
  real(8)        :: threshold_norm_pot
  character(1)   :: omp_loop
  character(1)   :: skip_gsortho
  integer        :: iditer_notemperature
  character(1)   :: gscg

!! &emfield
  character(2)   :: trans_longi
  character(16)  :: ae_shape1
  real(8)        :: e_impulse
  real(8)        :: amplitude1
  real(8)        :: rlaser_int_wcm2_1
  real(8)        :: pulse_tw1
  real(8)        :: omega1
  real(8)        :: epdir_re1(3)
  real(8)        :: epdir_im1(3)
  real(8)        :: phi_cep1
  character(16)  :: ae_shape2
  real(8)        :: amplitude2
  real(8)        :: rlaser_int_wcm2_2
  real(8)        :: pulse_tw2
  real(8)        :: omega2
  real(8)        :: epdir_re2(3)
  real(8)        :: epdir_im2(3)
  real(8)        :: phi_cep2
  real(8)        :: t1_t2
  real(8)        :: t1_delay
  character(1)   :: quadrupole
  character(8)   :: quadrupole_pot
  character(1)   :: alocal_laser
  real(8)        :: rlaserbound_sta(3)
  real(8)        :: rlaserbound_end(3)
  integer        :: nump
  real(8)        :: vecp(3,2)
  real(8)        :: coop(3,2)
  real(8)        :: radp_diele

!! &multiscale
  character(16)  :: fdtddim
  character(16)  :: twod_shape
  integer        :: nx_m
  integer        :: ny_m
  integer        :: nz_m
  real(8)        :: hx_m
  real(8)        :: hy_m
  real(8)        :: hz_m
  integer        :: nksplit !! TODO: remove this variable
  integer        :: nxysplit !! TODO: remove this variable
  ! The input variables nxvac(l|r)_m do not recommend to use,
  ! However I tempolary remain them for the reason of the compatibility.
  ! Please use  n(x|y|z)_origin_m to provide the same functionality.
  integer        :: nxvacl_m 
  integer        :: nxvacr_m
  integer        :: nx_origin_m
  integer        :: ny_origin_m
  integer        :: nz_origin_m
  character(100) :: file_macropoint
  integer        :: num_macropoint
  character(1)   :: set_ini_coor_vel
  integer        :: nmacro_write_group
  !! TODO: remove num_macropoint later

!! &maxwell
  real(8)        :: al_em(3)
  real(8)        :: dl_em(3)
  real(8)        :: dt_em
  integer        :: nt_em
  integer        :: iboundary(3,2)
  character(16)  :: wave_input
  real(8)        :: ek_dir1(3)
  real(8)        :: source_loc1(3)
  real(8)        :: ek_dir2(3)
  real(8)        :: source_loc2(3)
  integer        :: iobs_num_em
  integer        :: iobs_samp_em
  real(8)        :: obs_loc_em(200,3)
  character(256) :: shape_file
  integer        :: imedia_num
  character(16)  :: type_media(0:200)
  real(8)        :: epsilon(0:200)
  real(8)        :: rmu(0:200)
  real(8)        :: sigma(0:200)
  real(8)        :: omega_p_d(0:200)
  real(8)        :: gamma_d(0:200)
  character(1)   :: smooth_d
  real(8)        :: weight_d
  character(1)   :: wf_em
  
!! &analysis
  character(2)   :: projection_option
  character(4)   :: projection_decomp
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
  character(1)   :: out_old_dns
  character(1)   :: out_dns_rt
  integer        :: out_dns_rt_step
  character(1)   :: out_dns_trans
  real(8)        :: out_dns_trans_energy
  character(1)   :: out_elf
  character(1)   :: out_elf_rt
  integer        :: out_elf_rt_step
  character(1)   :: out_estatic_rt
  integer        :: out_estatic_rt_step
  character(1)   :: out_rvf_rt
  integer        :: out_rvf_rt_step
  character(1)   :: out_tm
  integer        :: out_projection_step
  integer        :: out_ms_step
  character(16)  :: format3d
  integer        :: numfiles_out_3d
  character(1)   :: timer_process

!! &hartree
  integer        :: meo
  integer        :: num_pole_xyz(3)

!! &ewald
  integer        :: newald
  real(8)        :: aewald

!! &opt
  real(8)        :: cg_alpha_ini
  real(8)        :: cg_alpha_up
  real(8)        :: cg_alpha_down
  real(8)        :: convrg_scf_force
  real(8)        :: convrg_scf_ene
  real(8)        :: convrg_opt_fmax
  real(8)        :: convrg_opt_ene

!! &md
  character(10)  :: ensemble
  character(20)  :: thermostat
  integer        :: step_velocity_scaling
  integer        :: step_update_ps
  integer        :: step_update_ps2
  real(8)        :: temperature0_ion
  character(1)   :: set_ini_velocity
  character(256) :: file_ini_velocity
  character(256) :: file_set_shake
  real(8)        :: thermostat_tau
  real(8)        :: friction
  character(1)   :: stop_system_momt

!! &misc
  character(4)   :: fourier

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
integer,allocatable :: kion(:)    
real(8),allocatable :: rion(:,:)  
real(8),allocatable :: rion_red(:,:)  
character(1),allocatable :: flag_geo_opt_atom(:)
character(256),allocatable :: atom_name(:)


end module salmon_global
