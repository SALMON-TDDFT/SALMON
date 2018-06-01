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
module salmon_xc
  use builtin_pz, only: exc_cor_pz
  use builtin_pzm, only: exc_cor_pzm
  use builtin_pbe, only: exc_cor_pbe
  use builtin_tbmbj, only: exc_cor_tbmbj
#ifdef SALMON_USE_LIBXC
  use xc_f90_types_m
  use xc_f90_lib_m
#endif

  implicit none

! List of Exchange Correlation Functionals
  integer, parameter :: salmon_xctype_pz    = 1
  integer, parameter :: salmon_xctype_pzm   = 2
  integer, parameter :: salmon_xctype_pbe   = 3
  integer, parameter :: salmon_xctype_tbmbj = 4
  integer, parameter :: salmon_xctype_tpss  = 5
  integer, parameter :: salmon_xctype_vs98  = 6
#ifdef SALMON_USE_LIBXC
  integer, parameter :: salmon_xctype_libxc_nonmag = 101
  integer, parameter :: salmon_xctype_libxc_spin = 102
  integer, parameter :: salmon_xctype_libxc_spinor = 103
#endif

  type xc_functional
    character(32) :: xcname
    integer :: xctype
    integer :: ispin
    real(8) :: cval
    logical :: use_gradient
    logical :: use_laplacian
    logical :: use_kinetic_energy
    logical :: use_current
#ifdef SALMON_USE_LIBXC
    integer :: nfunc
    type(xc_f90_pointer_t) :: func(2)
    type(xc_f90_pointer_t) :: info(2)
#endif
  end type

contains

  subroutine print_xc_info()
    implicit none
    integer :: vmajor, vminor, vmicro
    
#ifdef SALMON_USE_LIBXC
    call xc_f90_version(vmajor, vminor, vmicro)
    print '(2X,A,1X,I1,".",I1,".",I1)', "Libxc: [enabled]", vmajor, vminor, vmicro
#else
    print '(2X,A)', "Libxc: [disabled]"
#endif
    return
  end subroutine print_xc_info
    

  subroutine init_xc(xcname, ispin, cval, xc)
    implicit none
    character(*), intent(in)         :: xcname
    integer, intent(in)              :: ispin
    real(8), intent(in)              :: cval
    type(xc_functional), intent(out) :: xc
    
    character(32) :: xcprefix, xcname2
    integer :: ipos
        
    xc%xcname = xcname
    xc%ispin = ispin
    xc%cval = cval

    ! Separate a prefix part of the xcname delimitted by ':' letter.
    ! e.g. 'LIBXC: LDA_X' -> xcprefix='LIBXC' xcname='LDA_X'
    ipos = index(xcname, ':')
    if (1 <= ipos) then
      xcprefix = trim(adjustl(xcname(1:ipos-1)))
      xcname2 = trim(adjustl(xcname(ipos+1:)))
    else
      xcprefix = ''
      xcname2 = trim(adjustl(xcname))
    end if


    if (0 < len_trim(xcprefix)) then
      select case(xcprefix)

#ifdef SALMON_USE_LIBXC
      ! NOTE: xcname is starts with 'libxc:' e.g. 'libxc:lda_x'
      case ('libxc', 'LIBXC')
        call setup_libxc_nonmag(xcname2)
        return
        
#endif
  
      case default
        print *, "Undefined functional group prefix:" // trim(xcprefix)
        stop
  
      end select
    end if

    ! Given xcname has non-prefix part:
    select case(xcname2)
    case ('pz', 'PZ')
      xc%xctype = salmon_xctype_pz
      xc%use_gradient = .false.
      xc%use_laplacian = .false.
      xc%use_kinetic_energy = .false.
      xc%use_current = .false.
      return
      
    case ('pzm', 'PZM')
      xc%xctype = salmon_xctype_pzm
      xc%use_gradient = .false.
      xc%use_laplacian = .false.
      xc%use_kinetic_energy = .false.
      xc%use_current = .false.
      return
      
    case ('pbe', 'PBE')
      xc%xctype = salmon_xctype_pbe
      xc%use_gradient = .true.
      xc%use_laplacian = .false.
      xc%use_kinetic_energy = .false.
      xc%use_current = .false.
      return

    case ('tbmbj', 'TBmBJ', 'TBMBJ')
      xc%xctype = salmon_xctype_tbmbj
      xc%use_gradient = .true.
      xc%use_laplacian = .true.
      xc%use_kinetic_energy = .true.
      xc%use_current = .true.
      return
      
    case ('bj_pw', 'BJ_PW')
      xc%xctype = salmon_xctype_tbmbj; xc%cval = 1d0
      xc%use_gradient = .true.
      xc%use_laplacian = .true.
      xc%use_kinetic_energy = .true.
      xc%use_current = .true.
      return
      
    case ('tpss', 'TPSS')
      xc%xctype = salmon_xctype_tpss
      xc%use_gradient = .true.
      xc%use_laplacian = .true.
      xc%use_kinetic_energy = .true.
      xc%use_current = .true.
      return
      
    case ('vs98', 'VS98')
      xc%xctype = salmon_xctype_vs98
      xc%use_gradient = .true.
      xc%use_laplacian = .true.
      xc%use_kinetic_energy = .true.
      xc%use_current = .true.
      return

#ifdef SALMON_USE_LIBXC
    ! Presently only for spin-unpolarized xc only...
    case('libxc_pz', 'LIBXC_PZ')
      ! LIBXC_PZ is equivalent to 'LIBXC: LDA_X + LDA_C_PZ'
      call setup_libxc_nonmag('LDA_X + LDA_C_PZ')
      return
      
    case('libxc_pzm', 'LIBXC_PZM')
      ! LIBXC_PZM is equivalent to 'LIBXC: LDA_X + LDA_C_PZ_MOD'
      call setup_libxc_nonmag('LDA_X + LDA_C_PZ_MOD')
      return

    case('libxc_pbe', 'LIBXC_PBE')
      ! LIBXC_PBE is equivalent to 'LIBXC: GGA_X_PBE + GGA_C_PBE'
      call setup_libxc_nonmag('GGA_X_PBE + GGA_C_PBE')
      return
      
#endif

    case default  
      print *, "Undefined functional type: " // trim(xcname2)
      stop
      
    end select

    return
  
  contains
    
    
    
#ifdef SALMON_USE_LIBXC
    subroutine setup_libxc_nonmag(code)
      implicit none
      character(*), intent(in) :: code
      character(32) :: item(2)
      integer :: jpos, ii, iixc
      
      xc%xctype = salmon_xctype_libxc_nonmag
      xc%use_gradient = .false.       ! initial setting
      xc%use_laplacian = .false.      ! initial setting
      xc%use_kinetic_energy = .false. ! initial setting
      xc%use_current = .false.        ! libxc does not uses current...
      
      jpos = index(code, '+')
      if (1 <= jpos) then
        xc%nfunc = 2
        item(1) = code(:jpos-1)
        item(2) = code(jpos+1:)
      else
        xc%nfunc = 1
        item(1) = code(:)
      end if
      
      do ii = 1, xc%nfunc
        iixc = xc_f90_functional_get_number(trim(adjustl(item(ii))))
        
        if (iixc <= 0) then
          print '(2X,"Error! Unknown Libxc functional:",A)', trim(adjustl(item(ii)))
          stop
        end if
        
        call xc_f90_func_init(xc%func(ii), xc%info(ii), iixc, XC_UNPOLARIZED)
        
        select case(xc_f90_info_family(xc%info(ii)))
        case (XC_FAMILY_LDA)
        case (XC_FAMILY_GGA)
          xc%use_gradient = .true.
        case ( XC_FAMILY_MGGA)
          xc%use_gradient = .true.
          xc%use_laplacian = .true.
          xc%use_kinetic_energy = .true.
        case (XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
          print '(2X,"Error! Hybrid is not available:",A)', trim(adjustl(item(ii)))
          stop
        end select
      end do
      
      return
    end subroutine setup_libxc_nonmag
#endif
        
        
        
  end subroutine init_xc



  subroutine finalize_xc(xc)
    implicit none
    type(xc_functional), intent(inout) :: xc
    integer :: ii
    
#ifdef SALMON_USE_LIBXC
    do ii = 1, xc%nfunc
      call xc_f90_func_end(xc%func(ii))
    end do
#endif
    return
  end subroutine



  subroutine calc_xc(xc, rho, rho_s, exc, eexc, vxc, vxc_s, &
      & grho, grho_s, rlrho, rlrho_s, tau, tau_s, rj, rj_s, &
      & rho_nlcc, &
      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz, Hxyz, aLxyz)
    implicit none
    type(xc_functional), intent(in) :: xc
    real(8), intent(in), optional :: rho(:, :, :) ! ispin = 1
    real(8), intent(in), optional :: rho_s(:, :, :, :) ! ispin = 2
    real(8), intent(out), optional :: exc(:, :, :) ! epsilon_xc[rho]
    real(8), intent(out), optional :: eexc(:, :, :) ! rho * epsilon_xc[rho]
    real(8), intent(out), optional :: vxc(:, :, :) ! v_xc[rho] for ispin=1
    real(8), intent(out), optional :: vxc_s(:, :, :, :) ! v_xc[rho] ispin=2
    !real(8), intent(out), optional :: gvxc(:, :, :) ! v_xc[rho] for ispin=1
    !real(8), intent(out), optional :: gvxc_s(:, :, :, :) ! v_xc[rho] ispin=2
    real(8), intent(in), optional :: grho(:, :, :, :)
    real(8), intent(in), optional :: grho_s(:, :, :, :, :) ! ispin = 2
    real(8), intent(in), optional :: rlrho(:, :, :)
    real(8), intent(in), optional :: rlrho_s(:, :, :, :) ! ispin = 2
    real(8), intent(in), optional :: rj(:, :, :, :)
    real(8), intent(in), optional :: rj_s(:, :, :, :) ! ispin = 2
    real(8), intent(in), optional :: tau(:, :, :)
    real(8), intent(in), optional :: tau_s(:, :, :, :) ! ispin = 2
  
    real(8), intent(in), optional :: rho_nlcc(:, :, :)
    
    !===============================================================
    ! NOTE:
    !   The following section (finite difference table) is required 
    !   in the GGA/mGGA functionals which is originally used by ARTED subroutine. 
    ! TODO:
    !   Prepare more simplified solution to call the built-in GGA/mGGA functionals
    integer, intent(in), optional :: nd
    integer, intent(in), optional :: ifdx(:, :) 
    integer, intent(in), optional :: ifdy(:, :)
    integer, intent(in), optional :: ifdz(:, :)
    real(8), intent(in), optional :: nabx(:)
    real(8), intent(in), optional :: naby(:)
    real(8), intent(in), optional :: nabz(:)
    real(8), intent(in), optional :: Hxyz
    real(8), intent(in), optional :: aLxyz
    !===============================================================
    
    integer :: nx, ny, nz, nl

    ! Detect size of 3-dimensional grid
    if (xc%ispin == 1) then
      nx = ubound(rho, 1) - lbound(rho, 1) + 1;
      ny = ubound(rho, 2) - lbound(rho, 2) + 1;
      nz = ubound(rho, 3) - lbound(rho, 3) + 1;
    else
      nx = ubound(rho_s, 1) - lbound(rho_s, 1) + 1;
      ny = ubound(rho_s, 2) - lbound(rho_s, 2) + 1;
      nz = ubound(rho_s, 3) - lbound(rho_s, 3) + 1;
    end if
    nl = nx * ny * nz

    select case (xc%xctype)
    case(salmon_xctype_pz) 
      call exec_builtin_pz()
      
    case(salmon_xctype_pzm) 
      call exec_builtin_pzm()
      
    case(salmon_xctype_pbe) 
      call exec_builtin_pbe()
    
    case(salmon_xctype_tbmbj) 
      call exec_builtin_tbmbj()
    
#ifdef SALMON_USE_LIBXC
    case(salmon_xctype_libxc_nonmag)
      call exec_libxc_nonmag()
#endif

#ifdef SALMON_DEBUG_XC
    case default
      print *, "Undefined functional type: ", xc%xctype
      stop
#endif

    end select
    return

  contains



    subroutine exec_builtin_pz()
      implicit none
      real(8) :: rho_s_1d(nl)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      
      rho_s_1d = reshape(rho, (/nl/)) * 0.5
      
      rho_s_1d = reshape(rho(:, :, :), (/nl/)) * 0.5
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif
      
      call exc_cor_pz(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)
      
      if (present(vxc)) then
         vxc = reshape(vexc_1d, (/nx, ny, nz/))
      endif
      
      if (present(exc)) then
         exc = reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = reshape(eexc_1d, (/nx, ny, nz/))
      endif
      
      return
    end subroutine exec_builtin_pz



    subroutine exec_builtin_pzm()
      implicit none
      real(8) :: rho_s_1d(nl)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      
      rho_s_1d = reshape(rho, (/nl/)) * 0.5
      
      rho_s_1d = reshape(rho(:, :, :), (/nl/)) * 0.5
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif
      
      call exc_cor_pzm(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)
      
      if (present(vxc)) then
         vxc = reshape(vexc_1d, (/nx, ny, nz/))
      endif
      
      if (present(exc)) then
         exc = reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = reshape(eexc_1d, (/nx, ny, nz/))
      endif
      
      return
    end subroutine exec_builtin_pzm
    


    subroutine exec_builtin_pbe()
      implicit none
      real(8) :: rho_1d(nl)
      real(8) :: grho_s_1d(nl, 3)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      
      rho_1d = reshape(rho, (/nl/))
      grho_s_1d = reshape(grho(:, :, :, :), (/nl, 3/)) * 0.5
      
      call exc_cor_pbe(nl, rho_1d, grho_s_1d, exc_1d, eexc_1d, vexc_1d, &
      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz)
      
      if (present(vxc)) then
         vxc = reshape(vexc_1d, (/nx, ny, nz/))
      endif
      
      if (present(exc)) then
         exc = reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = reshape(eexc_1d, (/nx, ny, nz/))
      endif
      
      return
    end subroutine exec_builtin_pbe



    subroutine exec_builtin_tbmbj()
      implicit none
      real(8) :: rho_1d(nl)
      real(8) :: rho_s_1d(nl)
      real(8) :: grho_s_1d(nl, 3)
      real(8) :: rlrho_s_1d(nl)
      real(8) :: tau_s_1d(nl)
      real(8) :: j_s_1d(nl, 3)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      integer :: ii

      rho_1d = reshape(rho, (/nl/))
      rho_s_1d = rho_1d * 0.5
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif

      grho_s_1d = reshape(grho(:, :, :, :), (/nl, 3/)) * 0.5
      rlrho_s_1d = reshape(rlrho(:, :, :), (/nl/)) * 0.5
      tau_s_1d = reshape(tau(:, :, :), (/nl/)) * 0.5
      j_s_1d = reshape(rj(:, :, :, :), (/nl, 3/)) * 0.5

      call exc_cor_tbmbj(nl, rho_1d, rho_s_1d,  grho_s_1d, rlrho_s_1d, tau_s_1d, j_s_1d, xc%cval, eexc_1d, vexc_1d, Hxyz, aLxyz)

      if (present(vxc)) then
        vxc = reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
        ! NOTE: Take care for "zero-division error"
        exc = reshape(eexc_1d, (/nx, ny, nz/)) / rho
      endif

      if (present(eexc)) then
        eexc = reshape(eexc_1d, (/nx, ny, nz/))
      endif

      return
    end subroutine exec_builtin_tbmbj



#ifdef SALMON_USE_LIBXC
    subroutine exec_libxc_nonmag()
      implicit none
      real(8) :: rho_1d(nl), grho_1d(nl, 3)
      real(8) :: sigma_1d(nl), rlrho_1d(nl), tau_1d(nl)
      real(8) :: exc_1d(nl), exc_tmp_1d(nl)
      real(8) :: vxc_1d(nl), vxc_tmp_1d(nl)
      real(8) :: gvxc_tmp_1d(nl), lvxc_tmp_1d(nl), tvxc_tmp_1d(nl)
      integer :: ii
      
      exc_1d = 0d0
      vxc_1d = 0d0
      
      rho_1d = reshape(rho, (/nl/))
      if (present(rho_nlcc)) then
        rho_1d = rho_1d + reshape(rho_nlcc, (/nl/))
      end if
      
      if (xc%use_gradient) then
        sigma_1d = reshape( &
          & grho(:,:,:,1)**2+grho(:,:,:,2)**2+grho(:,:,:,3)**2, (/nl/) &
          & )
      end if
      
      if (xc%use_laplacian) then
        rlrho_1d = reshape(rlrho, (/nl/))
      end if
      
      if (xc%use_kinetic_energy) then
        tau_1d = reshape(tau, (/nl/))
      end if
  
      do ii = 1, xc%nfunc
        select case (xc_f90_info_family(xc%info(ii)))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc_vxc( &
            & xc%func(ii), nl, rho_1d(1), &
            & exc_tmp_1d(1), vxc_tmp_1d(1) &
            & )

        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc_vxc( &
            & xc%func(ii), nl, rho_1d(1), sigma_1d(1), &
            & exc_tmp_1d(1), vxc_tmp_1d(1), gvxc_tmp_1d(1) &
            & )
        
        case(XC_FAMILY_MGGA)
          call xc_f90_mgga_exc_vxc( &
            & xc%func(ii), nl, rho_1d(1), sigma_1d(1), rlrho_1d(1), tau_1d(1), &
            & exc_tmp_1d(1), vxc_tmp_1d(1), &
            & gvxc_tmp_1d(1), lvxc_tmp_1d(1), tvxc_tmp_1d(1))
          
        case(XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
          print *, "Hybrid Functional is not implemented"

        case default
          print *, "Error! SALMON does not support this functional family." 
          stop
          
        end select
        
        exc_1d = exc_1d + exc_tmp_1d
        vxc_1d = vxc_1d + vxc_tmp_1d
        
      end do
      
      if (present(exc)) then
        exc = reshape(exc_1d, (/nx, ny, nz/))
      endif
      
      if (present(eexc)) then
        eexc = reshape(exc_1d * rho_1d, (/nx, ny, nz/))
      endif
      
      if (present(vxc)) then
         vxc = reshape(vxc_1d, (/nx, ny, nz/))
      endif
      
      return
    end subroutine exec_libxc_nonmag
    
#endif

  end subroutine calc_xc
end module
