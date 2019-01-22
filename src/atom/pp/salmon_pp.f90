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
    character(2),allocatable :: atom_symbol(:)
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

  type pp_grid
    integer :: nps
    integer,allocatable :: mps(:)
    integer,allocatable :: jxyz(:,:,:)
    integer,allocatable :: jxx(:,:)
    integer,allocatable :: jyy(:,:)
    integer,allocatable :: jzz(:,:)
    real(8),allocatable :: uv(:,:)
    real(8),allocatable :: duv(:,:,:)
    integer :: nlma
    integer,allocatable :: lma_tbl(:,:)
    integer,allocatable :: ia_tbl(:)
    real(8),allocatable :: rinv_uvu(:)
  end type

  contains

  subroutine init_pp(pp,nrmax,lmax,flag_nlcc)
    use salmon_global,only : nelem,lloc_ps
    use salmon_global,only : pseudo_file
    use salmon_global,only : n_Yabana_Bertsch_psformat,n_ABINIT_psformat, &
                             n_ABINITFHI_psformat,n_FHI_psformat, &
                             ps_format,izatom,nelem
    use salmon_parallel, only: nproc_group_global, nproc_id_global
    use salmon_communication, only: comm_bcast, comm_is_root
    implicit none
    type(pp_info) :: pp
    integer, parameter :: nrmax0=50000, lmax0=4
    integer,intent(in) :: nrmax,lmax
    logical,intent(in) :: flag_nlcc
    integer :: ik
    character(256) :: ps_file
    integer :: ips_type,nlen_psf

    allocate(pp%atom_symbol(nelem))
    allocate(pp%rmass(nelem))
    allocate(pp%mr(nelem))

    if (comm_is_root(nproc_id_global)) then
  
      do ik=1,nelem
         
        ps_file = trim(pseudo_file(ik))
        nlen_psf = len_trim(ps_file)
        if(ps_file(nlen_psf+1-8:nlen_psf) == '_rps.dat')then
          ips_type = n_Yabana_Bertsch_psformat
          ps_format(ik) = 'KY'
          call read_mr_yb(pp,ik,ps_file)
        else if(ps_file(nlen_psf+1-6:nlen_psf) == '.pspnc')then
          ips_type = n_ABINIT_psformat
          ps_format(ik) = 'ABINIT'
          call read_mr_abinit(pp,ik,ps_file)
        else if(ps_file(nlen_psf+1-4:nlen_psf) == '.fhi')then
          ips_type = n_ABINITFHI_psformat
          ps_format(ik) = 'ABINITFHI'
          call read_mr_abinitfhi(pp,ik,ps_file)
        else if(ps_file(nlen_psf+1-4:nlen_psf) == '.cpi')then
          ips_type = n_FHI_psformat
          ps_format(ik) = 'FHI'
          call read_mr_fhi(pp,ik,ps_file)
        else
          stop 'Unprepared ps_format is required input_pseudopotential_YS'
        end if

!! --- Making prefix ---
!      select case (ipsfileform(ik))
!      case(n_Yabana_Bertsch_psformat)   ; ps_postfix = '_rps.dat'
!      case(n_ABINIT_psformat)           ; ps_postfix = '.pspnc'
!      case(n_ABINITFHI_psformat)        ; ps_postfix = '.fhi'
!      case(n_FHI_psformat)              ; ps_postfix = '.cpi'
!!    case('ATOM')      ; ps_postfix = '.psf' !Not implemented yet
!      case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
!      end select

! --- input pseudopotential and wave function ---
        select case (izatom(ik))
        case (1) ; pp%atom_symbol(ik) = 'H ' ; pp%rmass(ik)=1.d0
        case (2) ; pp%atom_symbol(ik) = 'He' ; pp%rmass(ik)=4.d0
        case (3) ; pp%atom_symbol(ik) = 'Li' ; pp%rmass(ik)=7.d0
        case (4) ; pp%atom_symbol(ik) = 'Be' ; pp%rmass(ik)=9.d0
        case (5) ; pp%atom_symbol(ik) = 'B ' ; pp%rmass(ik)=11.d0
        case (6) ; pp%atom_symbol(ik) = 'C ' ; pp%rmass(ik)=12.d0
        case (7) ; pp%atom_symbol(ik) = 'N ' ; pp%rmass(ik)=14.d0
        case (8) ; pp%atom_symbol(ik) = 'O ' ; pp%rmass(ik)=16.d0
        case (9) ; pp%atom_symbol(ik) = 'F ' ; pp%rmass(ik)=19.d0
        case(10) ; pp%atom_symbol(ik) = 'Ne' ; pp%rmass(ik)=20.d0
        case(11) ; pp%atom_symbol(ik) = 'Na' ; pp%rmass(ik)=23.d0
        case(12) ; pp%atom_symbol(ik) = 'Mg' ; pp%rmass(ik)=24.d0
        case(13) ; pp%atom_symbol(ik) = 'Al' ; pp%rmass(ik)=27.d0
        case(14) ; pp%atom_symbol(ik) = 'Si' ; pp%rmass(ik)=28.d0
        case(15) ; pp%atom_symbol(ik) = 'P ' ; pp%rmass(ik)=31.d0
        case(16) ; pp%atom_symbol(ik) = 'S ' ; pp%rmass(ik)=32.d0
        case(17) ; pp%atom_symbol(ik) = 'Cl' ; pp%rmass(ik)=35.d0
        case(18) ; pp%atom_symbol(ik) = 'Ar' ; pp%rmass(ik)=40.d0
        case(19) ; pp%atom_symbol(ik) = 'K ' ; pp%rmass(ik)=39.d0
        case(20) ; pp%atom_symbol(ik) = 'Ca' ; pp%rmass(ik)=40.d0
        case(21) ; pp%atom_symbol(ik) = 'Sc' ; pp%rmass(ik)=45.d0
        case(22) ; pp%atom_symbol(ik) = 'Ti' ; pp%rmass(ik)=48.d0
        case(23) ; pp%atom_symbol(ik) = 'V ' ; pp%rmass(ik)=51.d0
        case(24) ; pp%atom_symbol(ik) = 'Cr' ; pp%rmass(ik)=52.d0
        case(25) ; pp%atom_symbol(ik) = 'Mn' ; pp%rmass(ik)=55.d0
        case(26) ; pp%atom_symbol(ik) = 'Fe' ; pp%rmass(ik)=56.d0
        case(27) ; pp%atom_symbol(ik) = 'Co' ; pp%rmass(ik)=59.d0
        case(28) ; pp%atom_symbol(ik) = 'Ni' ; pp%rmass(ik)=59.d0
        case(29) ; pp%atom_symbol(ik) = 'Cu' ; pp%rmass(ik)=63.d0
        case(30) ; pp%atom_symbol(ik) = 'Zn' ; pp%rmass(ik)=65.d0
        case(31) ; pp%atom_symbol(ik) = 'Ga' ; pp%rmass(ik)=69.d0
        case(32) ; pp%atom_symbol(ik) = 'Ge' ; pp%rmass(ik)=73.d0
        case(33) ; pp%atom_symbol(ik) = 'As' ; pp%rmass(ik)=75.d0
        case(34) ; pp%atom_symbol(ik) = 'Se' ; pp%rmass(ik)=79.d0
        case(35) ; pp%atom_symbol(ik) = 'Br' ; pp%rmass(ik)=80.d0
        case(36) ; pp%atom_symbol(ik) = 'Kr' ; pp%rmass(ik)=84.d0
        case(37) ; pp%atom_symbol(ik) = 'Rb' ; pp%rmass(ik)=85.d0
        case(38) ; pp%atom_symbol(ik) = 'Sr' ; pp%rmass(ik)=88.d0
        case(39) ; pp%atom_symbol(ik) = 'Y ' ; pp%rmass(ik)=89.d0
        case(40) ; pp%atom_symbol(ik) = 'Zr' ; pp%rmass(ik)=91.d0
        case(41) ; pp%atom_symbol(ik) = 'Nb' ; pp%rmass(ik)=93.d0
        case(42) ; pp%atom_symbol(ik) = 'Mo' ; pp%rmass(ik)=96.d0
        case(43) ; pp%atom_symbol(ik) = 'Tc' ; pp%rmass(ik)=98.d0
        case(44) ; pp%atom_symbol(ik) = 'Ru' ; pp%rmass(ik)=101.d0
        case(45) ; pp%atom_symbol(ik) = 'Rh' ; pp%rmass(ik)=103.d0
        case(46) ; pp%atom_symbol(ik) = 'Pd' ; pp%rmass(ik)=106.d0
        case(47) ; pp%atom_symbol(ik) = 'Ag' ; pp%rmass(ik)=108.d0
        case(48) ; pp%atom_symbol(ik) = 'Cd' ; pp%rmass(ik)=112.d0
        case(49) ; pp%atom_symbol(ik) = 'In' ; pp%rmass(ik)=115.d0
        case(50) ; pp%atom_symbol(ik) = 'Sn' ; pp%rmass(ik)=119.d0
        case(51) ; pp%atom_symbol(ik) = 'Sb' ; pp%rmass(ik)=122.d0
        case(52) ; pp%atom_symbol(ik) = 'Te' ; pp%rmass(ik)=128.d0
        case(53) ; pp%atom_symbol(ik) = 'I ' ; pp%rmass(ik)=127.d0
        case(54) ; pp%atom_symbol(ik) = 'Xe' ; pp%rmass(ik)=131.d0
        case(55) ; pp%atom_symbol(ik) = 'Cs' ; pp%rmass(ik)=133.d0
        case(56) ; pp%atom_symbol(ik) = 'Ba' ; pp%rmass(ik)=137.d0
        case(57) ; pp%atom_symbol(ik) = 'La' ; pp%rmass(ik)=139.d0
        case(58) ; pp%atom_symbol(ik) = 'Ce' ; pp%rmass(ik)=140.d0
        case(59) ; pp%atom_symbol(ik) = 'Pr' ; pp%rmass(ik)=141.d0
        case(60) ; pp%atom_symbol(ik) = 'Nd' ; pp%rmass(ik)=144.d0
        case(61) ; pp%atom_symbol(ik) = 'Pm' ; pp%rmass(ik)=145.d0
        case(62) ; pp%atom_symbol(ik) = 'Sm' ; pp%rmass(ik)=150.d0
        case(63) ; pp%atom_symbol(ik) = 'Eu' ; pp%rmass(ik)=152.d0
        case(64) ; pp%atom_symbol(ik) = 'Gd' ; pp%rmass(ik)=157.d0
        case(65) ; pp%atom_symbol(ik) = 'Tb' ; pp%rmass(ik)=159.d0
        case(66) ; pp%atom_symbol(ik) = 'Dy' ; pp%rmass(ik)=164.d0
        case(67) ; pp%atom_symbol(ik) = 'Ho' ; pp%rmass(ik)=165.d0
        case(68) ; pp%atom_symbol(ik) = 'Er' ; pp%rmass(ik)=167.d0
        case(69) ; pp%atom_symbol(ik) = 'Tm' ; pp%rmass(ik)=169.d0
        case(70) ; pp%atom_symbol(ik) = 'Yb' ; pp%rmass(ik)=173.d0
        case(71) ; pp%atom_symbol(ik) = 'Lu' ; pp%rmass(ik)=175.d0
        case(72) ; pp%atom_symbol(ik) = 'Hf' ; pp%rmass(ik)=178.d0
        case(73) ; pp%atom_symbol(ik) = 'Ta' ; pp%rmass(ik)=181.d0
        case(74) ; pp%atom_symbol(ik) = 'W ' ; pp%rmass(ik)=184.d0
        case(75) ; pp%atom_symbol(ik) = 'Re' ; pp%rmass(ik)=186.d0
        case(76) ; pp%atom_symbol(ik) = 'Os' ; pp%rmass(ik)=190.d0
        case(77) ; pp%atom_symbol(ik) = 'Ir' ; pp%rmass(ik)=192.d0
        case(78) ; pp%atom_symbol(ik) = 'Pt' ; pp%rmass(ik)=195.d0
        case(79) ; pp%atom_symbol(ik) = 'Au' ; pp%rmass(ik)=197.d0
        case(80) ; pp%atom_symbol(ik) = 'Hg' ; pp%rmass(ik)=201.d0
        case(81) ; pp%atom_symbol(ik) = 'Tl' ; pp%rmass(ik)=204.d0
        case(82) ; pp%atom_symbol(ik) = 'Pb' ; pp%rmass(ik)=207.d0
        case(83) ; pp%atom_symbol(ik) = 'Bi' ; pp%rmass(ik)=209.d0
        case default ; stop 'Unprepared atomic data is called input_pseudopotential_YS'
        end select

      end do

    end if

    call comm_bcast(pp%atom_symbol,nproc_group_global)
    call comm_bcast(pp%rmass,nproc_group_global)
    call comm_bcast(ps_format,nproc_group_global)

    pp%lmax0=lmax0

    pp%nrmax0=nrmax0
  
    pp%nrmax=nrmax
    pp%lmax=lmax
  
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
!======================================================================
  subroutine read_mr_yb(pp,ik,ps_file)
    implicit none
    type(pp_info) :: pp
    integer :: ik
    character(256) :: ps_file
    
    open(4,file=ps_file,status='old')
    read(4,*) pp%mr(ik)
    close(4)
    return
  
  end subroutine read_mr_YB
!======================================================================
  subroutine read_mr_abinit(pp,ik,ps_file)
    implicit none
    type(pp_info) :: pp
    integer :: ik
    character(256) :: ps_file
    real(8) :: zatom, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    read(4,*) dummy_text
    read(4,*) zatom, zion, pspdat
    read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
    close(4)
    
    pp%mr(ik)=mmax
  
  end subroutine read_mr_ABINIT
!======================================================================
  subroutine read_mr_abinitfhi(pp,ik,ps_file)
    implicit none
    type(pp_info) :: pp
    integer :: ik
    integer :: i
    character(256) :: ps_file
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    do i=1,18
      read(4,*) dummy_text
    end do
    read(4,*) pp%mr(ik)
    close(4)
    
  end subroutine read_mr_abinitfhi

!======================================================================
  subroutine read_mr_fhi(pp,ik,ps_file)
    implicit none
    type(pp_info) :: pp
    integer :: ik
    integer :: i
    character(256) :: ps_file
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    do i=1,11
      read(4,*) dummy_text
    end do
    read(4,*) pp%mr(ik)
    close(4)
    
  end subroutine read_mr_fhi

end module salmon_pp
