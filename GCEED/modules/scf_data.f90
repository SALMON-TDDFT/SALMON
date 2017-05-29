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
MODULE scf_data
use salmon_global
implicit none
include 'mpif.h'
!-------------------- Parameters
integer, parameter :: maxntmg=10

! Physical constants
real(8),parameter :: E2=14.4d0, H2M=3.81d0, a_B=0.529177d0
real(8),parameter :: Ry=13.6058d0, Pi=3.141592653589793d0
real(8),parameter :: fs2eVinv = 1.51925d0

! Coefficients of finite-difference
integer,parameter :: Nd=4       ! (2*Nd+1)-points finite-difference used in Taylor expansion
integer,parameter :: Ndh=4       ! (2*Nd+1)-points finite-difference used in Hartree routine
real(8),parameter :: cN0=-205.d0/72.d0                 ! for Laplacian
real(8),parameter :: cN1=8.d0/5.d0  , cN2=-1.d0/5.d0
real(8),parameter :: cN3=8.d0/315.d0, cN4=-1.d0/560.d0
real(8),parameter :: bN1=4.d0/5.d0  , bN2=-1.d0/5.d0   ! for Nabla
real(8),parameter :: bN3=4.d0/105.d0, bN4=-1.d0/280.d0

real(8),parameter :: cN0_Nd1=-2.d0
real(8),parameter :: cN1_Nd1=1.d0
real(8),parameter :: bN1_Nd1=1.d0/2.d0

real(8),parameter :: cN0_Nd2=-5.d0/2.d0
real(8),parameter :: cN1_Nd2=4.d0/3.d0, cN2_Nd2=-1.d0/12.d0
real(8),parameter :: bN1_Nd2=2.d0/3.d0, bN2_Nd2=-1.d0/12.d0

real(8),parameter :: cN0_Nd3=-49.d0/18.d0 
real(8),parameter :: cN1_Nd3=3.d0/2.d0, cN2_Nd3=-3.d0/20.d0
real(8),parameter :: cN3_Nd3=1.d0/90.d0
real(8),parameter :: bN1_Nd3=3.d0/4.d0, bN2_Nd3=-3.d0/20.d0 
real(8),parameter :: bN3_Nd3=1.d0/60.d0

real(8),parameter :: cN0_Nd5=-5269.d0/1800.d0 
real(8),parameter :: cN1_Nd5=5.d0/3.d0, cN2_Nd5=-5.d0/21.d0
real(8),parameter :: cN3_Nd5=5.d0/126.d0, cN4_Nd5=-5.d0/1008.d0
real(8),parameter :: cN5_Nd5=1.d0/3150.d0

real(8),parameter :: cN0_Nd6=-5369.d0/1800.d0 
real(8),parameter :: cN1_Nd6=12.d0/7.d0, cN2_Nd6=-15.d0/56.d0
real(8),parameter :: cN3_Nd6=10.d0/189.d0, cN4_Nd6=-1.d0/112.d0
real(8),parameter :: cN5_Nd6=2.d0/1925.d0, cN6_Nd6=-1.d0/16632.d0

real(8),parameter :: cN0_Nd7=-266681.d0/88200.d0
real(8),parameter :: cN1_Nd7=7.d0/4.d0, cN2_Nd7=-7.d0/24.d0
real(8),parameter :: cN3_Nd7=7.d0/108.d0, cN4_Nd7=-7.d0/528.d0
real(8),parameter :: cN5_Nd7=7.d0/3300.d0, cN6_Nd7=-7.d0/30888.d0
real(8),parameter :: cN7_Nd7=1.d0/84084.d0

real(8),parameter :: cN0_Nd8=-1077749.d0/352800.d0 
real(8),parameter :: cN1_Nd8=16.d0/9.d0, cN2_Nd8=-14.d0/45.d0
real(8),parameter :: cN3_Nd8=112.d0/1485.d0, cN4_Nd8=-7.d0/396.d0
real(8),parameter :: cN5_Nd8=112.d0/32175.d0, cN6_Nd8=-2.d0/3861.d0
real(8),parameter :: cN7_Nd8=16.d0/315315.d0, cN8_Nd8=-1.d0/411840.d0

real(8),parameter :: cN0_Nd9=-9778141.d0/3175200.d0
real(8),parameter :: cN1_Nd9=9.d0/5.d0, cN2_Nd9=-18.d0/55.d0
real(8),parameter :: cN3_Nd9=14.d0/165.d0, cN4_Nd9=-63.d0/2860.d0
real(8),parameter :: cN5_Nd9=18.d0/3575.d0, cN6_Nd9=-2.d0/2145.d0
real(8),parameter :: cN7_Nd9=9.d0/70070.d0, cN8_Nd9=-9.d0/777920.d0
real(8),parameter :: cN9_Nd9=1.d0/1969110.d0

real(8),parameter :: cN0_Nd10=-1968329.d0/635040.d0
real(8),parameter :: cN1_Nd10=20.d0/11.d0, cN2_Nd10=-15.d0/44.d0
real(8),parameter :: cN3_Nd10=40.d0/429.d0, cN4_Nd10=-15.d0/572.d0
real(8),parameter :: cN5_Nd10=24.d0/3575.d0, cN6_Nd10=-5.d0/3432.d0
real(8),parameter :: cN7_Nd10=30.d0/119119.d0, cN8_Nd10=-5.d0/155584.d0
real(8),parameter :: cN9_Nd10=10.d0/3741309.d0, cN10_Nd10=-1.d0/9237800.d0

real(8),parameter :: cN0_Nd11=-239437889.d0/76839840.d0
real(8),parameter :: cN1_Nd11=11.d0/6.d0, cN2_Nd11=-55.d0/156.d0
real(8),parameter :: cN3_Nd11=55.d0/546.d0, cN4_Nd11=-11.d0/364.d0
real(8),parameter :: cN5_Nd11=11.d0/1300.d0, cN6_Nd11=-11.d0/5304.d0
real(8),parameter :: cN7_Nd11=55.d0/129948.d0, cN8_Nd11=-55.d0/806208.d0
real(8),parameter :: cN9_Nd11=11.d0/1360476.d0, cN10_Nd11=-11.d0/17635800.d0
real(8),parameter :: cN11_Nd11=1.d0/42678636.d0

real(8),parameter :: cN0_Nd12=-240505109.d0/76839840.d0
real(8),parameter :: cN1_Nd12=24.d0/13.d0, cN2_Nd12=-33.d0/91.d0
real(8),parameter :: cN3_Nd12=88.d0/819.d0, cN4_Nd12=-99.d0/2912.d0
real(8),parameter :: cN5_Nd12=396.d0/38675.d0, cN6_Nd12=-11.d0/3978.d0
real(8),parameter :: cN7_Nd12=132.d0/205751.d0, cN8_Nd12=-33.d0/268736.d0
real(8),parameter :: cN9_Nd12=44.d0/2380833.d0, cN10_Nd12=-3.d0/1469650.d0
real(8),parameter :: cN11_Nd12=12.d0/81800719.d0, cN12_Nd12=-1.d0/194699232.d0

real(8) :: cNmat(0:12,0:12),bNmat(0:12,0:12)

!-------------------- Global variables

integer :: iflag_ps

integer :: inumcpu_check

integer :: ierr,nproc,myrank

integer,allocatable :: ob_sta_all_kgrid(:),ob_end_all_kgrid(:),iobnum_all_kgrid(:)

integer :: nproc_ob, nproc_Mxin(3),nproc_Mxin_mul
integer :: nproc_ob_spin(2)
integer :: nproc_Mxin_s(3), nproc_Mxin_mul_s
integer :: nproc_Mxin_s_dm(3), nproc_Mxin_mul_s_dm

integer :: ista_Mx_ori(3),iend_Mx_ori(3),inum_Mx_ori(3)
integer,allocatable :: ista_Mxin(:,:),iend_Mxin(:,:),inum_Mxin(:,:)
integer,allocatable :: ista_Mxin_old(:,:),iend_Mxin_old(:,:),inum_Mxin_old(:,:)
integer,allocatable :: ista_Mxin_s(:,:),iend_Mxin_s(:,:),inum_Mxin_s(:,:)
integer :: max_lg_num

integer :: Miter       ! Total number of Iteration for SCF calculation
integer :: Miter_rt    ! Total number of Iteration for RT calculation

real(8) :: mixrate     ! Mixing rate for updating charge density
integer :: Nmemory_MB  ! number of iteration to stock for Broyden's method

integer :: minroutine  ! Type of routine for minimization

integer :: iflag_diisjump
real(8) :: lambda1_diis, lambda2_diis

integer :: iwrite_external
integer :: iflag_writepsi
integer :: iflag_dip2
integer :: iflag_quadrupole
integer :: iflag_intelectron
integer :: num_dip2
real(8),allocatable :: rto(:)
real(8) :: dip2boundary(100)
real(8) :: dip2center(100)
integer ,allocatable:: idip2int(:)
real(8),allocatable :: rto_ix(:,:)
real(8),allocatable :: rbox_array_dip2(:,:)
real(8),allocatable :: rbox_array2_dip2(:,:)
real(8),allocatable :: rbox_array_dip2q(:,:,:)
real(8),allocatable :: rbox_array2_dip2q(:,:,:)
real(8),allocatable :: rbox_array_dip2e(:)
real(8),allocatable :: rbox_array2_dip2e(:)

integer :: ilsda

integer :: iflag_stopt
integer :: iter_stopt
integer :: istopt

integer :: ntmg

integer :: MST(2),ifMST(2),itotMST
integer :: itotfMST
integer :: MST0(2),itotMST0
integer :: Mx(3),Mxin(3),Mxin_old(3)

integer,allocatable :: Kion(:)    ! Label for ions
real(8),allocatable :: Rion(:,:)  ! Ion coordinates
character(8),allocatable :: AtomName(:)   
integer,allocatable :: iAtomicNumber(:)   
integer,allocatable :: istopt_a(:)    

real(8) :: Hgs(3)        ! Grid spacing
real(8) :: Hold(3)     ! Grid spacing
real(8) :: Hvol
real(8) :: Harray(3,maxntmg)  ! Grid spacing
real(8) :: rLsize(3,maxntmg)    ! size of the box
integer :: ithresholdVh(maxntmg)  ! threshold value for Vh iteration

integer :: MEO
integer :: maxMps
integer :: lmax_MEO

! Pseudopotential
integer,allocatable :: Jxyz(:,:,:),Mps(:),Jxxyyzz(:,:,:)
integer,allocatable :: Jxyz_tmp1(:,:,:)
integer,allocatable :: Jxyz_tmp2(:,:,:)
integer,allocatable :: Jxxyyzz_tmp1(:,:,:)
integer,allocatable :: Jxxyyzz_tmp2(:,:,:)
integer,allocatable :: Jxxyyzz2nd(:,:,:)
integer :: Mlmps
integer :: Mlps(maxMKI),Lref(maxMKI)
real(8),allocatable :: Vpsl(:,:,:)                 ! Local pseudopotential
real(8),allocatable :: Vpsl_atom(:,:,:,:)
real(8),allocatable :: uV(:,:,:),uVu(:,:)          ! Non-local

real(8),allocatable :: rocc(:,:)                    ! Occupation number
real(8),allocatable :: psi(:,:,:,:,:)              ! Single particle orbitals
real(8),allocatable :: psi_mesh(:,:,:,:,:)         ! Single particle orbitals
real(8),allocatable :: zpsi_mesh(:,:,:,:,:)        ! Single particle orbitals
real(8),allocatable :: psi_old(:,:,:,:,:)
complex(8),allocatable :: zpsi_old(:,:,:,:,:)

real(8),allocatable :: rho_in(:,:,:,:)
real(8),allocatable :: rho_out(:,:,:,:)
real(8),allocatable :: rho_s_in(:,:,:,:,:)
real(8),allocatable :: rho_s_out(:,:,:,:,:)

real(8),allocatable :: esp(:,:)         ! Single particle energy
real(8),allocatable :: esp2(:,:)        ! Single particle energy
real(8),allocatable :: rho(:,:,:)       ! Single particle density
real(8),allocatable :: rho0(:,:,:)      ! Single particle density
real(8),allocatable :: rho_diff(:,:,:)  ! Single particle density
real(8),allocatable :: rho_s(:,:,:,:)   ! Single particle density for each spin
real(8),allocatable :: Vh(:,:,:)        ! Hartree potential
real(8),allocatable :: Vxc(:,:,:)       ! Exchange-Correlation potential
real(8),allocatable :: Vxc_s(:,:,:,:)   ! Exchange-Correlation potential for each spin 
integer :: ihpsieff
real(8),allocatable :: elf(:,:,:)
complex(8),allocatable :: zpsi(:,:,:,:,:)
complex(8),allocatable :: zpsi_in(:,:,:,:,:),zpsi_out(:,:,:,:,:)
complex(8),allocatable :: zpsi_t0(:,:,:,:,:)
integer :: iSCFRT

real(8),allocatable :: Veff(:,:,:),Veff2(:,:,:,:),Vbox(:,:,:)
real(8),allocatable :: Veff3(:,:,:),Veff4(:,:,:),Veff5(:,:,:,:)

real(8),allocatable :: Ex_fast(:,:,:),Ec_fast(:,:,:)

real(8) :: elp3(3000)
real(8) :: elp5(3000)

integer, allocatable :: idiis_sd(:)

integer :: iwksta(3),iwkend(3),iwknum(3)
integer :: iwk2sta(3),iwk2end(3),iwk2num(3)
integer :: iwk3sta(3),iwk3end(3),iwk3num(3)
integer :: iwk_size
integer :: imesh_s_all

integer :: numcoo
integer,allocatable :: jMps_l(:,:)
integer,allocatable :: max_jMps_l(:)

integer :: numcoo_s
integer,allocatable :: jMps_l_s(:,:)
integer,allocatable :: max_jMps_l_s(:)

integer :: maxlm

integer :: imesh_oddeven
integer :: iswitch_orbital_mesh

integer :: version_num(2)

complex(8), allocatable :: zc(:)
real(8), allocatable :: Dp(:,:), Dp2(:,:,:)
real(8), allocatable :: Qp(:,:,:), Qp2(:,:,:,:)
real(8), allocatable :: rIe(:), rIe2(:,:)
real(8) :: vecDs(3)
real(8),allocatable :: vecDs2(:,:)
real(8) :: vecQs(3,3)
real(8),allocatable :: vecQs2(:,:,:)

integer :: iflag_psicube

integer :: num_pole
integer :: num_pole_xyz(3)

integer :: Ncg

integer :: itotNtime

integer :: num_datafiles_IN
integer :: num_datafiles_OUT
integer :: num_datafiles_OUT2

integer :: iflag_hartree

real(8),allocatable :: Gs(:,:,:),Gl(:,:,:)
complex(8),allocatable :: tx_exp(:,:),ty_exp(:,:),tz_exp(:,:)

integer :: lg_sta(3),lg_end(3),lg_num(3)
integer :: mg_sta(3),mg_end(3),mg_num(3)
integer :: ng_sta(3),ng_end(3),ng_num(3)

integer :: lg_old_sta(3),lg_old_end(3),lg_old_num(3)
integer :: mg_old_sta(3),mg_old_end(3),mg_old_num(3)
integer :: ng_old_sta(3),ng_old_end(3),ng_old_num(3)

real(8),allocatable :: gridcoo(:,:)

integer :: iobnum

integer :: kx_hock_sta(3),kx_hock_end(3),kx_hock_num(3)

real(8),allocatable :: wtk(:)

real(8),allocatable :: norm_diff_psi_stock(:,:)

real(8) :: Etot
real(8) :: Exc

integer :: imr(3),imrs(3),igroup

integer :: iflag_convergence
real(8) :: threshold_norm_diff_rho(maxntmg)
real(8) :: threshold_square_norm_diff_Vlocal(maxntmg)

real(8),allocatable :: rho_stock(:,:,:,:)
real(8),allocatable :: Vlocal_stock(:,:,:,:)
integer,parameter :: num_rho_stock=21

integer :: iflag_subspace_diag
integer :: iDiter_nosubspace_diag

real(8),allocatable :: Vh_stock1(:,:,:)
real(8),allocatable :: Vh_stock2(:,:,:)

real(8),allocatable :: Vlocal(:,:,:,:)
real(8),allocatable :: Vlocal2(:,:,:,:)

! use for hartree routine
integer :: iterVh
real(8) :: Hconv

real(8),allocatable :: rhobox(:,:,:)
real(8),allocatable :: rhobox_s(:,:,:,:)

integer :: lg_num_fmax(3)

integer :: iflag_comm_rho
real(8), allocatable :: rhobox1_all(:,:,:), rhobox2_all(:,:,:)

integer :: iDiter(maxntmg)

integer :: IC,OC

character(LEN=100) :: file_OUT_ini
integer :: num_mol
real(8) :: rlatcon
real(8), allocatable :: coo_mol_ini(:,:)

integer :: lg_sta_ini(3),lg_end_ini(3),lg_num_ini(3)
integer :: mg_sta_ini(3),mg_end_ini(3),mg_num_ini(3)
real(8) :: H_ini(3)
real(8) :: rLsize_ini(3)

integer :: img

integer :: itt
integer :: ikind_eext   !0:No external field, 1: dipoleApprox
integer :: N_hamil

character(3)  :: dir
character(2)  :: dir2 

real(8) :: dt,Fst,Fst2(2)
real(8) :: romega, romega2(2)
real(8) :: pulse_T, pulse_T2(2) 
real(8) :: rlaser_I, rlaser_I2(2) 
real(8) :: tau, tau2(2), delay, rcycle

integer       :: iwdenoption
integer       :: iwdenstep ! step for writing density
character(2)  :: denplane  ! plane for writing density (xy, yz, xz)
integer       :: idensum   ! whether density is summed up along direction
                           ! perpendicular to the plane
                           ! (0: not summed, 1: summed)
real(8)       :: posplane  ! position of the plane
                           ! (only for idensum = 0)
integer       :: numfile_movie

real(8) :: fourier_omega(200)
integer :: num_fourier_omega
integer :: iflag_fourier_omega

integer :: itotNtime2

character(100):: rtOutFile
character(100):: rtDiffOutFile
character(100):: rtELFOutFile
character(100):: file_Projection
character(40):: fileTmp
character(40):: fileTmp2
character(20):: fileTmp3
character(20):: fileNumber
character(20):: fileELF

integer,allocatable :: Jxyz2nd(:,:,:)
real(8),allocatable :: uV2nd(:,:,:)

integer,allocatable :: numatom_ps(:,:,:)
integer,allocatable :: iatomnum_ps(:,:,:,:)
integer,allocatable :: Jxyz_all(:,:)
integer :: maxMps_all
integer,allocatable :: Mps3rd(:,:,:,:)

integer :: jamax
integer,allocatable :: Jxyz_all_2nd(:,:)
integer,allocatable :: iatomnum_ps_2nd(:)
integer,allocatable :: Mps3rd_2nd(:)
integer,allocatable :: jja(:,:)
integer :: MImax
integer,allocatable :: numatom_ps_2nd(:)

integer :: iflag_Estatic
real(8), allocatable :: rho_n(:,:,:)
real(8), allocatable :: Vh_n(:,:,:)
real(8), allocatable :: Vh0(:,:,:)
complex(8), allocatable :: zpsi_n(:,:,:,:,:)

complex(8), allocatable :: Ex_static(:,:,:),Ey_static(:,:,:),Ez_static(:,:,:)

real(8) :: lasbound_sta(3),lasbound_end(3)
integer :: ilasbound_sta(3),ilasbound_end(3)
real(8) :: lascenter(3)

integer :: isequential

complex(8) :: cumnum

real(8) :: Eion

integer :: iblacsinit
integer :: CONTEXT, IAM, MYCOL, MYROW, NPCOL, NPROCS2, NPROW
integer :: DESCA( 50 ), DESCZ( 50 )

real(8) :: rho_region1(100)
real(8) :: rho_region2(100)
integer :: num_rho_region(100)
integer,allocatable :: rho_region_nx(:,:)
integer,allocatable :: rho_region_ny(:,:)
integer,allocatable :: rho_region_nz(:,:)

integer :: numspin

integer :: icalcforce
real(8),allocatable :: rforce(:,:)
integer :: iflag_md
real(8),allocatable :: dRion(:,:,:)
real(8),allocatable :: Rion_eq(:,:)
real(8),parameter :: umass=1822.9d0

real(8),allocatable :: Eeff_dip(:,:,:,:)  
real(8),allocatable :: Eeff_dip0(:,:,:,:)  

integer :: idisnum(2)
integer :: wmaxMI

real(8) :: rad_diele

real(8) :: fcN(0:12)
real(8) :: fbN(0:12)

integer,parameter :: NEwald=4

integer,allocatable :: oblist(:)

integer :: iwrite_projection
integer :: iwrite_projnum
integer :: itwproj

integer :: iparaway_ob
integer :: iscf_order

integer :: num_projection
integer :: iwrite_projection_ob(200)
integer :: iwrite_projection_k(200)

integer,allocatable::icoo1d(:,:)

character(100) :: filename_pot

integer :: MI_read

integer :: ik_oddeven

! variables for FFTE routine
integer,dimension(3) :: LNPU
integer :: NPUZ,NPUY

integer :: itcalc_ene

real(8) :: absorption(0:100000)
real(8) :: absorption_d(0:100000)
real(8) :: absorption_id(0:100000)

integer :: iflag_dos
integer :: iflag_pdos

integer :: iflag_ELF

CONTAINS

!=========================================================================
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!=========================================================================
subroutine snum_procs
implicit none
integer :: i,ibox,ipow

nproc_Mxin_mul=nproc/nproc_ob

ibox=1
ipow=0
do i=1,29
  if(ibox<nproc_Mxin_mul)then
    ipow=ipow+1
    ibox=ibox*2
  end if
end do

nproc_Mxin(1:3)=1
do i=1,ipow
  if(mod(i,3)==1)then
    nproc_Mxin(3)=nproc_Mxin(3)*2
  else if(mod(i,3)==2)then
    nproc_Mxin(2)=nproc_Mxin(2)*2
  else
    nproc_Mxin(1)=nproc_Mxin(1)*2
  end if
end do

ibox=1
ipow=0
do i=1,29
  if(ibox<nproc)then
    ipow=ipow+1
    ibox=ibox*2
  end if
end do

nproc_Mxin_s(1:3)=1
do i=1,ipow
  if(mod(i,3)==1)then
    nproc_Mxin_s(3)=nproc_Mxin_s(3)*2
  else if(mod(i,3)==2)then
    nproc_Mxin_s(2)=nproc_Mxin_s(2)*2
  else
    nproc_Mxin_s(1)=nproc_Mxin_s(1)*2
  end if
end do

end subroutine snum_procs

!======================================================================
subroutine init_mesh_s
implicit none

nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

allocate(ista_Mxin_s(3,0:nproc-1),iend_Mxin_s(3,0:nproc-1))
allocate(inum_Mxin_s(3,0:nproc-1))

call setng(ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s, &
           nproc,myrank,nproc_Mxin,nproc_Mxin_s_dm,ista_Mxin,iend_Mxin,inum_Mxin,isequential)

end subroutine init_mesh_s

!=====================================================================
subroutine make_iwksta_iwkend
implicit none

if(iwk_size==1)then
  iwksta(1:3)=ista_Mxin(1:3,myrank)
  iwkend(1:3)=iend_Mxin(1:3,myrank)
  iwk2sta(1:3)=ista_Mxin(1:3,myrank)-Nd
  iwk2end(1:3)=iend_Mxin(1:3,myrank)+Nd
  iwk3sta(1:3)=ista_Mxin(1:3,myrank)
  iwk3end(1:3)=iend_Mxin(1:3,myrank)
else if(iwk_size==2)then
  iwksta(1:3)=ista_Mxin(1:3,myrank)-Nd
  iwkend(1:3)=iend_Mxin(1:3,myrank)+Nd
  iwk2sta(1:3)=ista_Mxin(1:3,myrank)-Nd
  iwk2end(1:3)=iend_Mxin(1:3,myrank)+Nd
  iwk3sta(1:3)=ista_Mxin(1:3,myrank)
  iwk3end(1:3)=iend_Mxin(1:3,myrank)
else if(iwk_size==11.or.iwk_size==31)then
  iwksta(1:3)=ista_Mxin_s(1:3,myrank)
  iwkend(1:3)=iend_Mxin_s(1:3,myrank)
  iwk2sta(1:3)=ista_Mxin_s(1:3,myrank)-Ndh
  iwk2end(1:3)=iend_Mxin_s(1:3,myrank)+Ndh
  iwk3sta(1:3)=ista_Mxin_s(1:3,myrank)
  iwk3end(1:3)=iend_Mxin_s(1:3,myrank)
else if(iwk_size==12.or.iwk_size==32)then
  iwksta(1:3)=ista_Mxin_s(1:3,myrank)-Ndh
  iwkend(1:3)=iend_Mxin_s(1:3,myrank)+Ndh
  iwk2sta(1:3)=ista_Mxin_s(1:3,myrank)-Ndh
  iwk2end(1:3)=iend_Mxin_s(1:3,myrank)+Ndh
  iwk3sta(1:3)=ista_Mxin_s(1:3,myrank)
  iwk3end(1:3)=iend_Mxin_s(1:3,myrank)
end if
iwknum(1:3)=iwkend(1:3)-iwksta(1:3)+1
iwk2num(1:3)=iwk2end(1:3)-iwk2sta(1:3)+1
iwk3num(1:3)=iwk3end(1:3)-iwk3sta(1:3)+1

end subroutine make_iwksta_iwkend

END MODULE scf_data

