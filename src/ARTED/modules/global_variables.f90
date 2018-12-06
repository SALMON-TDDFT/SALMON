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
Module Global_Variables
  use salmon_global
  use salmon_xc, only: xc_functional
  use salmon_pp, only: pp_info
!ARTED version
  character(50),parameter :: ARTED_ver='ARTED.1.6.0 (based on 2014.08.10.2)'

! constants
  real(8),parameter :: Pi=3.141592653589793d0
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
  real(8),parameter :: umass=1822.9d0
  real(8),parameter :: c_light=137.03953250d0 ! sato

!yabana
!! DFT parameters
!  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
!  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
!  real(8),parameter :: CU=0.002d0,DU=-0.0116d0
!yabana

  integer :: iter_now,entrance_iter
  real(8) :: Time_start,Time_now

! grid
  integer,parameter :: Nd = 4
  integer :: NLx,NLy,NLz,NL,NG,NKx,NKy,NKz,NK,Sym,nGzero
  integer :: NKxyz 
  real(8) :: aLx,aLy,aLz,aLxyz
  real(8) :: bLx,bLy,bLz,Hx,Hy,Hz,Hxyz
  integer,allocatable :: Lx(:),Ly(:),Lz(:),Lxyz(:,:,:)
  integer,allocatable :: ifdx(:,:),ifdy(:,:),ifdz(:,:)
  real(8),allocatable :: Gx(:),Gy(:),Gz(:)
  real(8),allocatable :: lap(:),nab(:)
  real(8),allocatable :: lapx(:),lapy(:),lapz(:)
  real(8),allocatable :: nabx(:),naby(:),nabz(:)

! pseudopotential
  type(pp_info) :: pp
  integer,parameter :: Nrmax=3000,Lmax=4
  character(2) :: ps_type
  integer :: Nps,Nlma
  integer,allocatable :: lma_tbl(:,:)
  integer,allocatable :: Mps(:),Jxyz(:,:),Jxx(:,:),Jyy(:,:),Jzz(:,:)
  integer,allocatable :: Mlps(:),Lref(:),Zps(:),NRloc(:)
  integer,allocatable :: NRps(:),inorm(:,:),iuV(:),a_tbl(:)
  real(8),allocatable :: rad(:,:),Rps(:),vloctbl(:,:),udVtbl(:,:,:)
  real(8),allocatable :: radnl(:,:)
  real(8),allocatable :: Rloc(:),uV(:,:),duV(:,:,:),anorm(:,:)
  real(8),allocatable :: dvloctbl(:,:),dudVtbl(:,:,:)
  real(8),allocatable :: save_udVtbl_a(:,:,:)
  real(8),allocatable :: save_udVtbl_b(:,:,:)
  real(8),allocatable :: save_udVtbl_c(:,:,:)
  real(8),allocatable :: save_udVtbl_d(:,:,:)

! material
  integer :: NI,NE,NB,NBoccmax
  real(8) :: Ne_tot
  integer,allocatable :: Zatom(:)
  real(8),allocatable :: Mass(:),Rion_eq(:,:),dRion(:,:,:)
  real(8),allocatable :: occ(:,:),wk(:)

! physical quantities
  real(8) :: Eall,Eall0,Eall_GS0,jav(3),jav_ion(3),Tion
  real(8) :: Ekin,Eloc,Enl,Eh,Exc,Eion,Eelemag                      
  real(8),allocatable :: javt(:,:), javt_ion(:,:)
  real(8),allocatable :: Vpsl(:),Vpsl_ia(:,:),Vh(:),Vexc(:),Eexc(:),Vloc(:),Vloc_GS(:),Vloc_t(:)!yabana
  real(8),allocatable :: Eall_t(:),Tion_t(:),Temperature_ion_t(:),Ework_integ_fdR(:)
  real(8),allocatable :: Enh_t(:),Hnvt_t(:)
  real(8),allocatable :: Vloc_new(:),Vloc_old(:,:)
  real(8),allocatable :: tmass(:),tjr(:,:),tjr2(:),tmass_t(:),tjr_t(:,:),tjr2_t(:)
!yabana
  complex(8),allocatable :: dVloc_G(:,:), save_dVloc_G(:,:) 
  real(8),allocatable :: rho(:),rho_gs(:),rho_gs_t(:)
  real(8),allocatable :: rho_in(:,:),rho_out(:,:) !MB method
  complex(8),allocatable :: rhoe_G(:),rhoion_G(:)
  real(8),allocatable :: velocity(:,:),force(:,:),esp(:,:)
  real(8),allocatable :: Floc(:,:),Fnl(:,:),Fion(:,:),FionAc(:,:)          
  real(8),allocatable :: ovlp_occ_l(:,:),ovlp_occ(:,:)
  integer,allocatable :: assign_grid_atom(:) ! projection + decomp analysis
  integer,allocatable :: NBocc(:) !FS set
  real(8),allocatable :: esp_vb_min(:),esp_vb_max(:) !FS set
  real(8),allocatable :: esp_cb_min(:),esp_cb_max(:) !FS set
  real(8),allocatable :: Eall_GS(:),esp_var_ave(:),esp_var_max(:),ddns(:),ddns_abs_1e(:)
  real(8) :: xi_nh, kB, hartree2J
  parameter( kB = 1.38064852d-23 ) ![J/K]
  parameter( hartree2J = 4.359744650d-18 )

  logical :: flag_use_grad_wf_on_force
  parameter(flag_use_grad_wf_on_force=.true.)
!Nonlinear core correction
  logical :: flag_nlcc = .false.
  real(8),allocatable :: rho_nlcc_tbl(:,:),tau_nlcc_tbl(:,:)
  real(8),allocatable :: rho_nlcc(:),tau_nlcc(:)

! wave functions, work array
#ifdef _OPENACC
  complex(8),allocatable,pinned :: zu_t(:,:,:)
#else
  complex(8),allocatable :: zu_t(:,:,:)
#endif
  complex(8),allocatable :: zu_GS(:,:,:),zu_GS0(:,:,:)
  complex(8),allocatable :: tpsi(:),htpsi(:),zwork(:,:,:),ttpsi(:)
  real(8),allocatable :: work(:,:,:)
  real(8),allocatable :: esp_var(:,:)
  integer :: iflag_gs_init_wf=0   !use random number(=0)

! variables for 4-times loop in Fourier transportation
  integer,allocatable :: nxyz(:,:,:)
  real(8),allocatable :: rho_3D(:,:,:),Vh_3D(:,:,:)
  complex(8),allocatable :: rhoe_G_temp(:),rhoe_G_3D(:,:,:)
  complex(8),allocatable :: f1(:,:,:),f2(:,:,:),f3(:,:,:),f4(:,:,:)
  complex(8),allocatable :: eGx(:,:),eGy(:,:),eGz(:,:),eGxc(:,:),eGyc(:,:),eGzc(:,:)

! Bloch momentum,laser pulse, electric field
!  character(8) :: AE_shape
  real(8) :: f0_1 !,IWcm2_1,tpulsefs_1,omegaev_1,omega_1,tpulse_1,Epdir_1(3),phi_CEP_1 ! sato
  real(8) :: f0_2 !,IWcm2_2,tpulsefs_2,omegaev_2,omega_2,tpulse_2,Epdir_2(3),phi_CEP_2 ! sato
  real(8) :: T1_T2fs
  real(8),allocatable :: E_ext(:,:),E_ind(:,:),E_tot(:,:)
  real(8),allocatable :: kAc(:,:),kAc0(:,:),kAc_new(:,:)     !k+A(t)/c (kAc)
  real(8),allocatable :: Ac_ext(:,:),Ac_ind(:,:),Ac_tot(:,:) !A(t)/c (Ac)
  !alocal_laser option
  logical :: flag_set_ini_Ac_alocal=.false.
  real(8),allocatable :: weight_Ac_alocal(:), weight_Ac_alocal_ion(:)
  real(8),allocatable :: Ac_ext_al(:,:),Ac_tot_al(:,:)
  real(8),allocatable :: Ac1x_al(:),Ac1y_al(:),Ac1z_al(:),Ac2_al(:,:),divA_al(:)
  real(8) :: Ac_al_amp(3),Ac2_al_amp,nabt_al(12)
  !analysis of transition density at omega
  complex(8),allocatable :: rho_trans(:)

! control parameters
  real(8) :: dAc,domega
  integer :: Nomega

! file names, flags, etc
  character(256) :: file_GS,file_RT
  character(256) :: file_epst,file_epse
  character(256) :: file_force_dR,file_j_ac
  character(256) :: file_DoS,file_band
  character(256) :: file_dns,file_ovlp,file_nex,file_last_band_map,file_nex_atom
  character(256) :: file_dns_gs
  character(256) :: file_dns_rt
  character(256) :: file_dns_dlt
  character(256) :: file_psi_gs
  character(256) :: file_energy_transfer
  character(256) :: file_ac 
  character(256) :: file_k_data
  character(256) :: file_eigen_data
  character(256) :: file_tm_data
  character(256) :: file_rt_data
  character(256) :: file_lr_data
  character(256) :: process_directory
  character(1024):: dir_ms, dir_ms_RT_Ac
  character(1024),allocatable :: dir_ms_M(:)
  integer        :: ndivide_macro

! energy computation
  character(256) :: file_rt_energy_data
  integer,parameter :: nstep_energy_calc = 10

!  character(2) :: ext_field ! this variable is removed
!  character(2) :: Longi_Trans ! this variable is replaced by trans_longi
!  character(1) :: MD_option ! this variable is replaced by use_ehrenfest_md
!  character(2) :: AD_RHO !ovlp_option ! This variable is replaced by projection_option
!yabana
  character(10) :: functional
!yabana

  integer :: NK_ave,NG_ave,NK_s,NK_e,NG_s,NG_e
  integer :: NK_remainder,NG_remainder
! Timer
  character(10) :: position_option

! sato
  complex(8),allocatable :: ekr(:,:)  

! omp
  integer :: NUMBER_THREADS
  complex(8),allocatable :: ekr_omp(:,:,:), zproj(:,:,:)
  complex(8),allocatable :: tpsi_omp(:,:),ttpsi_omp(:,:),htpsi_omp(:,:)
  complex(8),allocatable :: xk_omp(:,:),hxk_omp(:,:),gk_omp(:,:),pk_omp(:,:),pko_omp(:,:),txk_omp(:,:)
  integer :: NKB
  integer,allocatable :: ik_table(:),ib_table(:)

  real(8),allocatable :: tau_s_l_omp(:,:),j_s_l_omp(:,:,:)

! sym
  integer,allocatable :: itable_sym(:,:) ! sym
  real(8),allocatable :: rho_l(:),rho_tmp1(:),rho_tmp2(:) !sym

! Finite temperature
  real(8) :: KbTev

! multi scale
  integer :: macRANK,kRANK ! sato
  
  integer, parameter :: nattr_column = 9
  !! Macropoint coordinates and attributes
  integer, allocatable :: macropoint(:,:)
  real(8), allocatable :: macropoint_attr(:,:)
  integer :: nmacro, nmacro_s, nmacro_e, nmacro_attr
  !! bg_media: Background media coordinates and attributes  
  integer, allocatable :: bg_media_point(:,:)
  real(8), allocatable :: bg_media_attr(:,:)
  integer :: nbg_media, nbg_media_attr
  !! acinit: Initial vector field data
  integer :: ninit_acfield
  integer, allocatable :: init_acfield_point(:,:)
  real(8), allocatable :: init_acfield_val(:,:)

  
  ! Size of Macroscopic FDTD Grid 
  integer :: nx1_m, ny1_m, nz1_m
  integer :: nx2_m, ny2_m, nz2_m
  ! Array size of Macroscopic FDTD Grid (including overlaped region)
  integer :: mx1_m, my1_m, mz1_m ! 
  integer :: mx2_m, my2_m, mz2_m
  integer,parameter :: novlp_m = 1

  !! Macroscopic electromagnetic field variables
  !! NOTE: "Ac_(old|new)?_ms" are the vector potential Ac(r,t)
  !!       "Jm_(old|new)?_ms" are the matter current density Jm(r,t)
  !!       In the RT iteration, the variable with suffix "new" and "old" 
  !!       indicate the data at the time "iter+1" and "iter-1", respectively.
  real(8),allocatable :: Ac_old_ms(:,:,:,:),  Ac_ms(:,:,:,:), Ac_new_ms(:,:,:,:)
  real(8),allocatable :: add_Ac_ms(:,:,:,:),  add_Ac_new_ms(:,:,:,:)
  real(8),allocatable :: Jm_old_ms(:,:,:,:),  Jm_ms(:,:,:,:), Jm_new_ms(:,:,:,:)
  real(8),allocatable :: Jm_ion_old_ms(:,:,:,:),Jm_ion_ms(:,:,:,:),Jm_ion_new_ms(:,:,:,:)
  real(8),allocatable :: elec_ms(:,:,:,:)
  real(8),allocatable :: bmag_ms(:,:,:,:)
  real(8),allocatable :: energy_joule_ms(:,:,:)
  real(8),allocatable :: energy_elec_ms(:,:,:), Eall0_m(:)
  real(8),allocatable :: energy_elemag_ms(:,:,:)
  real(8) :: total_energy_elemag_old, total_energy_elemag
  real(8) :: total_energy_absorb_old, total_energy_absorb
  real(8) :: total_energy_elec_old, total_energy_elec
  real(8) :: total_energy_em_old, total_energy_em  

  real(8) :: ms_angle_x, ms_angle_y, ms_angle_z
  real(8) :: trans_mat(3, 3), trans_inv(3, 3)

  !! Calculate Pure FDTD Calculation without TDDFT:
  !! NOTE: This switch will be removed after marging the common FDTD routine..
  logical :: debug_switch_no_radiation
  
#ifdef _OPENACC
  complex(8),allocatable,pinned :: zu_m(:,:,:,:)
#else
  complex(8),allocatable :: zu_m(:,:,:,:)
#endif
  real(8),allocatable :: Rion_m(:,:,:),velocity_m(:,:,:),force_m(:,:,:)
  real(8),allocatable :: Rion_eq_m(:,:,:),dRion_m(:,:,:,:)
  real(8),allocatable :: Vh_m(:,:)
  real(8),allocatable :: Vexc_m(:,:)
  real(8),allocatable :: Eexc_m(:,:)
  real(8),allocatable :: Vloc_m(:,:),Vloc_old_m(:,:,:)
  real(8),allocatable :: rho_m(:,:)
  
  real(8),allocatable :: Ac_m(:,:), Ac_new_m(:,:), Ac_old_m(:,:)
  real(8),allocatable :: Jm_m(:,:), Jm_new_m(:,:)
  real(8),allocatable :: Jm_ion_m(:,:), Jm_ion_new_m(:,:)
  real(8),allocatable :: jm_new_m_tmp(:,:),jm_ion_new_m_tmp(:,:)
  real(8),allocatable :: energy_elec_Matter_new_m(:)
  real(8),allocatable :: energy_elec_Matter_new_m_tmp(:)
  real(8),allocatable :: excited_electron_new_m(:)
  real(8),allocatable :: excited_electron_new_m_tmp(:)

  real(8), allocatable :: data_out(:,:,:,:,:)
  integer, parameter :: ndata_out_column = 15
  integer :: ndata_out, ndata_out_per_proc
  real(8), allocatable :: data_local_ac(:,:,:)
  real(8), allocatable :: data_local_jm(:,:,:)
  real(8), allocatable :: data_local_jm_ion(:,:,:)
  real(8), allocatable :: data_vac_ac(:,:,:)
  real(8), allocatable :: data_local_Tmp_ion(:,:)
  
  integer :: ix_detect_l, ix_detect_r, iy_detect, iz_detect


  character(30), parameter :: calc_mode_sc = 'singlecell'
  character(30), parameter :: calc_mode_ms = 'multiscale'

  integer :: restart_switch = 0

  logical :: need_backup      = .FALSE.

  ! calculation mode flag
  integer :: iflag_calc_mode
  integer,parameter :: iflag_calc_mode_gs_rt = 0
  integer,parameter :: iflag_calc_mode_gs    = 1
  integer,parameter :: iflag_calc_mode_rt    = 2

  ! calculation mode 
  integer, parameter :: calc_mode_gs = 1000
  integer, parameter :: calc_mode_rt = 1100

  ! scf
  integer :: PrLv_scf = 3          !(no print=0, all print=3)
  integer :: Nscf_conv, it_last_update_zu_GS_proj=-99
  logical :: flag_update_only_zu_GS=.false.
  logical :: flag_scf_conv_ene_force=.false.

  ! Rion update flag
  logical :: Rion_update_rt
  logical, parameter :: rion_update_on  = .true.
  logical, parameter :: rion_update_off = .false.
  
  ! Exchange Correlation
  type(xc_functional) :: xc_func

  !AY trial for Raman: Maxwell + Force-Field type MD (theory=Raman)
  logical :: flag_ms_ff_LessPrint
  integer :: Nm_FDTD, iter_save, imode_FDTD_raman
  integer :: interval_step_trj_raman
  real(8) :: Omg_dt,v_mxmt
  real(8) :: eps_diag, dchidq(3,3)
  real(8),allocatable :: c_pmode(:), Rion_eq0(:,:) 
  real(8),allocatable :: Rion_m_next(:,:,:),velocity_m_next(:,:,:)
  character(1024) :: dir_ion_trj

  interface 
    subroutine total_Energy_omp(Rion_update,GS_RT,ixy_m)
      integer,intent(in) :: GS_RT
      logical,intent(in) :: Rion_update
      integer,intent(in),optional :: ixy_m
    end subroutine total_Energy_omp

    subroutine Ion_Force_omp(Rion_update,GS_RT,ixy_m)
      integer,intent(in) :: GS_RT
      logical,intent(in) :: Rion_update
      integer,intent(in),optional :: ixy_m
    end subroutine Ion_Force_omp
  
    subroutine hpsi_omp_KB_base(ik,tpsi,htpsi,ttpsi)
      integer,intent(in)              :: ik
      complex(8),intent(in)           ::  tpsi
      complex(8),intent(out)          :: htpsi
      complex(8),intent(out),optional :: ttpsi
    end subroutine hpsi_omp_KB_base
  end interface

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
#else
# define MEM_ALIGNED 32
#endif

!dir$ attributes align:MEM_ALIGNED :: zu_t,zu_GS
!dir$ attributes align:MEM_ALIGNED :: xk_omp,txk_omp,hxk_omp,pko_omp
!dir$ attributes align:MEM_ALIGNED :: tpsi_omp,ttpsi_omp,htpsi_omp

End Module Global_Variables
