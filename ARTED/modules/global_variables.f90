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
Module Global_Variables
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
  character(10) :: entrance_option    !initial or reentrance
  real(8) :: Time_shutdown
  real(8) :: Time_start,Time_now

! grid
  integer :: NLx,NLy,NLz,Nd,NL,NG,NKx,NKy,NKz,NK,Sym,nGzero
  integer :: NKxyz 
  real(8) :: aL,ax,ay,az,aLx,aLy,aLz,aLxyz
  real(8) :: bLx,bLy,bLz,Hx,Hy,Hz,Hxyz
  integer,allocatable :: Lx(:),Ly(:),Lz(:),Lxyz(:,:,:)
  integer,allocatable :: ifdx(:,:),ifdy(:,:),ifdz(:,:)
  real(8),allocatable :: Gx(:),Gy(:),Gz(:)
  real(8),allocatable :: lap(:),nab(:)
  real(8),allocatable :: lapx(:),lapy(:),lapz(:)
  real(8),allocatable :: nabx(:),naby(:),nabz(:)

! pseudopotential
  integer,parameter :: Nrmax=3000,Lmax=4
  character(2) :: ps_type
  character(10) :: ps_format !shinohara
  character(1) :: PSmask_option = 'n' !shinohara
  real(8) :: alpha_mask = 0.8d0 , gamma_mask = 1.8d0, eta_mask = 15.d0 !shinohara
  integer :: Nps,Nlma
  integer,allocatable :: Mps(:),Jxyz(:,:),Jxx(:,:),Jyy(:,:),Jzz(:,:)
  integer,allocatable :: Mlps(:),Lref(:),Zps(:),NRloc(:)
  integer,allocatable :: NRps(:),inorm(:,:),iuV(:),a_tbl(:)
  real(8),allocatable :: rad(:,:),Rps(:),vloctbl(:,:),udVtbl(:,:,:)
  real(8),allocatable :: radnl(:,:)
  real(8),allocatable :: Rloc(:),uV(:,:),duV(:,:,:),anorm(:,:)
  real(8),allocatable :: dvloctbl(:,:),dudVtbl(:,:,:)

! material
  integer :: NI,NE,NB,NBoccmax
  real(8) :: Ne_tot
  integer,allocatable :: Zatom(:),Kion(:)
  real(8),allocatable :: Rion(:,:),Mass(:),Rion_eq(:,:),dRion(:,:,:)
  real(8),allocatable :: occ(:,:),wk(:)

! physical quantities
  real(8) :: Eall,Eall0,jav(3),Tion
  real(8) :: Ekin,Eloc,Enl,Eh,Exc,Eion,Eelemag                      
  real(8),allocatable :: javt(:,:)
  real(8),allocatable :: Vpsl(:),Vh(:),Vexc(:),Eexc(:),Vloc(:),Vloc_GS(:),Vloc_t(:)!yabana
  real(8),allocatable :: Vloc_new(:),Vloc_old(:,:)
  real(8),allocatable :: tmass(:),tjr(:,:),tjr2(:),tmass_t(:),tjr_t(:,:),tjr2_t(:)
!yabana
  complex(8),allocatable :: dVloc_G(:,:)
  real(8),allocatable :: rho(:),rho_gs(:)
  real(8),allocatable :: rho_in(:,:),rho_out(:,:) !MB method
  complex(8),allocatable :: rhoe_G(:),rhoion_G(:)
  real(8),allocatable :: force(:,:),esp(:,:),force_ion(:,:)
  real(8),allocatable :: Floc(:,:),Fnl(:,:),Fion(:,:)               
  real(8),allocatable :: ovlp_occ_l(:,:),ovlp_occ(:,:)
  integer :: Nelec !FS set
  integer,allocatable :: NBocc(:) !FS set
  real(8),allocatable :: esp_vb_min(:),esp_vb_max(:) !FS set
  real(8),allocatable :: esp_cb_min(:),esp_cb_max(:) !FS set
  real(8),allocatable :: Eall_GS(:),esp_var_ave(:),esp_var_max(:),dns_diff(:)
!Nonlinear core correction
  logical :: flag_nlcc = .false.
  real(8),allocatable :: rho_nlcc_tbl(:,:),tau_nlcc_tbl(:,:)
  real(8),allocatable :: rho_nlcc(:),tau_nlcc(:)

! wave functions, work array
  complex(8),allocatable :: zu_t(:,:,:),zu_GS(:,:,:),zu_GS0(:,:,:)
  complex(8),allocatable :: tpsi(:),htpsi(:),zwork(:,:,:),ttpsi(:)
  real(8),allocatable :: work(:,:,:)
  real(8),allocatable :: esp_var(:,:)

! variables for 4-times loop in Fourier transportation
  integer,allocatable :: nxyz(:,:,:)
  real(8),allocatable :: rho_3D(:,:,:),Vh_3D(:,:,:)
  complex(8),allocatable :: rhoe_G_temp(:),rhoe_G_3D(:,:,:)
  complex(8),allocatable :: f1(:,:,:),f2(:,:,:),f3(:,:,:),f4(:,:,:)
  complex(8),allocatable :: eGx(:,:),eGy(:,:),eGz(:,:),eGxc(:,:),eGyc(:,:),eGzc(:,:)

! Bloch momentum,laser pulse, electric field
!  real(8) :: f0,Wcm2,pulseT,wave_length,omega,pulse_time,pdir(3),phi_CEP=0.00*2*pi
  character(8) :: AE_shape
  real(8) :: f0_1,IWcm2_1,tpulsefs_1,omegaev_1,omega_1,tpulse_1,Epdir_1(3),phi_CEP_1 ! sato
  real(8) :: f0_2,IWcm2_2,tpulsefs_2,omegaev_2,omega_2,tpulse_2,Epdir_2(3),phi_CEP_2 ! sato
  real(8) :: T1_T2fs,T1_T2
  real(8),allocatable :: E_ext(:,:),E_ind(:,:),E_tot(:,:)
  real(8),allocatable :: kAc(:,:),kAc0(:,:),kAc_new(:,:)                  !k+A(t)/c (kAc)
  real(8),allocatable :: Ac_ext(:,:),Ac_ind(:,:),Ac_tot(:,:) !A(t)/c (Ac)

! control parameters
  integer :: NEwald=4                      !Ewald summation
  real(8) :: aEwald=0.5d0
  integer :: Ncg=10                        !# of conjugate gradient (cg)
  real(8) :: dt,dAc,domega
  integer :: Nscf,Nt,Nomega
  integer :: Nmemory_MB=8                   !Modified-Broyden (MB) method
  real(8) :: alpha_MB=0.75
  integer :: NFSset_start=75,NFSset_every=25 !Fermi Surface (FS) set 

! file names, flags, etc
  character(128) :: SYSname,directory
  character(256) :: file_GS,file_RT
  character(256) :: file_epst,file_epse
  character(256) :: file_force_dR,file_j_ac
  character(256) :: file_DoS,file_band
  character(256) :: file_dns,file_ovlp,file_nex
  character(256) :: file_kw              ! non-uniform k-grid
  character(256) :: file_energy_transfer ! 940
  character(256) :: file_ac_vac          ! 941
  character(256) :: file_ac_vac_back     ! 942
  character(256) :: file_ac_m            ! 943
  character(256) :: file_ac              ! 902
  character(256) :: file_ac_init         ! 902
  character(256) :: process_directory

  character(2) :: ext_field
  character(2) :: Longi_Trans
  character(1) :: FSset_option,MD_option
  character(2) :: AD_RHO !ovlp_option
!yabana
  character(10) :: functional
  real(8) :: cval ! cvalue for TBmBJ. If cval<=0, calculated in the program
!yabana
  character(10) :: propagator = 'default' ! propagation scheme: default, or etrs

  integer :: NK_ave,NG_ave,NK_s,NK_e,NG_s,NG_e
  integer :: NK_remainder,NG_remainder
! Timer
  character(10) :: position_option

! sato
  complex(8),allocatable :: ekr(:,:)  

! omp
  integer :: NUMBER_THREADS
  complex(8),allocatable :: ekr_omp(:,:,:)
  complex(8),allocatable :: tpsi_omp(:,:),ttpsi_omp(:,:),htpsi_omp(:,:)
  complex(8),allocatable :: xk_omp(:,:),hxk_omp(:,:),gk_omp(:,:),pk_omp(:,:),pko_omp(:,:),txk_omp(:,:)
  integer :: NKB
  integer,allocatable :: ik_table(:),ib_table(:)

  real(8),allocatable :: tau_s_l_omp(:,:),j_s_l_omp(:,:,:)

! sym
  character(50) :: crystal_structure !sym
  integer,allocatable :: itable_sym(:,:) ! sym
  real(8),allocatable :: rho_l(:),rho_tmp1(:),rho_tmp2(:) !sym

! Finite temperature
  real(8) :: KbTev

! multi scale
  character(20) :: FDTDdim,TwoD_shape
  integer :: NXY_s,NXY_e ! sato
  integer :: NKsplit,NXYsplit ! sato
  integer :: macRANK,kRANK ! sato

  real(8) :: HX_m,HY_m
  integer :: NX_m,NXvacL_m,NXvacR_m
  integer :: NY_m,NYvacT_m,NYvacB_m
  real(8),allocatable :: Ac_m(:,:,:),Ac_new_m(:,:,:),Ac_old_m(:,:,:)
  real(8),allocatable :: Elec(:,:,:),Bmag(:,:,:)
  real(8),allocatable :: j_m(:,:,:)
  real(8),allocatable :: jmatter_m(:,:,:),jmatter_m_l(:,:,:)
  real(8),allocatable :: g(:,:,:)
  real(8) :: bcon
  integer,allocatable :: NX_table(:),NY_table(:)
  character(50) BC_my

  complex(8),allocatable :: zu_m(:,:,:,:)
  real(8),allocatable :: Vh_m(:,:)
  real(8),allocatable :: Vexc_m(:,:)
  real(8),allocatable :: Eexc_m(:,:)
  real(8),allocatable :: Vloc_m(:,:),Vloc_old_m(:,:,:)
  real(8),allocatable :: rho_m(:,:)
  real(8),allocatable :: energy_joule(:,:)
  real(8),allocatable :: energy_elec_Matter_l(:,:)
  real(8),allocatable :: energy_elec_Matter(:,:)
  real(8),allocatable :: energy_elec(:,:)
  real(8),allocatable :: energy_elemag(:,:)
  real(8),allocatable :: energy_total(:,:)
  real(8),allocatable :: excited_electron_l(:,:)
  real(8),allocatable :: excited_electron(:,:)

  real(8),allocatable :: data_out(:,:,:,:)
  real(8),allocatable :: data_local_Ac(:,:,:)
  real(8),allocatable :: data_local_jm(:,:,:)
  real(8),allocatable :: data_vac_Ac(:,:,:)
  integer :: Nstep_write=100
  integer :: Ndata_out, Ndata_out_per_proc
  

  character(30) :: calc_mode
  character(30), parameter :: calc_mode_sc = 'singlecell'
  character(30), parameter :: calc_mode_ms = 'multiscale'

  integer :: reentrance_switch = 0

  integer :: backup_frequency = 10000 ! # of iteration
  logical :: need_backup      = .FALSE.


  ! calculation mode
  integer, parameter :: calc_mode_gs = 1000
  integer, parameter :: calc_mode_rt = 1100

  ! Rion update flag
  logical, parameter :: rion_update_on  = .true.
  logical, parameter :: rion_update_off = .false.

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
