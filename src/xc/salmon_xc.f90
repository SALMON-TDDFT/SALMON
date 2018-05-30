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
  integer, parameter :: salmon_xctype_pz   = 1
  integer, parameter :: salmon_xctype_pzm  = 2
  integer, parameter :: salmon_xctype_pbe  = 3
  integer, parameter :: salmon_xctype_bj   = 4
  integer, parameter :: salmon_xctype_tpss = 5
  integer, parameter :: salmon_xctype_vs98 = 6
#ifdef SALMON_USE_LIBXC
  integer, parameter :: salmon_xctype_libxc_lda = 101
  integer, parameter :: salmon_xctype_libxc_gga = 102 ! Future Plan
  integer, parameter :: salmon_xctype_libxc_mgga = 103 ! Future Plan
#endif

  type xc_functional
    character(16) :: xcname
    integer :: xctype
    integer :: ispin
    real(8) :: cval
#ifdef SALMON_USE_LIBXC
    type(xc_f90_pointer_t) :: x_func
    type(xc_f90_pointer_t) :: x_info
    type(xc_f90_pointer_t) :: c_func
    type(xc_f90_pointer_t) :: c_info
#endif
  end type

contains

  subroutine print_xc_info()
    implicit none
    integer :: major, minor, micro
    
#ifdef SALMON_USE_LIBXC
    call xc_version(major, minor, micro)
    print '(A,3(I3,"."))', 'Libxc version ', major, minor, micro
#endif
    print '(A)', 'SALMON version 1.0.0 builtin functionals'
    return
  end subroutine print_xc_info
    

  subroutine init_xc(xcname, ispin, cval, xc)
  implicit none
  character(*), intent(in)         :: xcname
  integer, intent(in)              :: ispin
  real(8), intent(in)              :: cval
  type(xc_functional), intent(out) :: xc

  xc%xcname = xcname
  xc%ispin = ispin
  xc%cval = cval

  select case(trim(xcname))
  case ('pz', 'PZ')
    xc%xctype = salmon_xctype_pz
    
  case ('pzm', 'PZM')
    xc%xctype = salmon_xctype_pzm
    
  case ('pbe', 'PBE')
    xc%xctype = salmon_xctype_pbe
    
  case ('tbmbj', 'TBmBJ', 'TBMBJ')
    xc%xctype = salmon_xctype_bj
    
  case ('bj_pw', 'BJ_PW')
    xc%xctype = salmon_xctype_bj; xc%cval = 1d0
    
  case ('tpss', 'TPSS')
    xc%xctype = salmon_xctype_tpss
    
  case ('vs98', 'VS98')
    xc%xctype = salmon_xctype_vs98

#ifdef SALMON_USE_LIBXC
  ! Presently only for spin-unpolarized xc only...
  case("libxc_pz")
    xc%xctype = salmon_xctype_libxc_lda
    call xc_f90_func_init(xc%x_func, xc%x_info, XC_LDA_X, XC_UNPOLARIZED)
    call xc_f90_func_init(xc%c_func, xc%c_info, XC_LDA_C_PZ, XC_UNPOLARIZED)
    
  case("libxc_pzm")
    xc%xctype = salmon_xctype_libxc_lda
    call xc_f90_func_init(xc%x_func, xc%x_info, XC_LDA_X, XC_UNPOLARIZED)
    call xc_f90_func_init(xc%c_func, xc%c_info, XC_LDA_C_PZ_MOD, XC_UNPOLARIZED)
#endif

#ifdef SALMON_DEBUG_XC
  case default
    print *, "Undefined functional type: " // trim(xcname)
    stop
#endif 
    
  end select

  return
  end subroutine init_xc


  subroutine finalize_xc(xc)
    implicit none
    type(xc_functional), intent(inout) :: xc
#ifdef SALMON_USE_LIBXC
    call xc_f90_func_end(xc%x_func)
    call xc_f90_func_end(xc%c_func)
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
    
#ifdef SALMON_USE_LIBXC
    case(salmon_xctype_libxc_lda)
          call exec_libxc_lda()
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
      real(8) :: grho_s_1d(nl)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      
      rho_1d = reshape(rho, (/nl/))
      grho_s_1d = reshape(grho(:, :, :, :), (/nl/)) * 0.5
      
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
      real(8) :: grho_s_1d(nl)
      real(8) :: rlrho_s_1d(nl)
      real(8) :: tau_s_1d(nl)
      real(8) :: j_s_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      integer :: ii
      
      rho_1d = reshape(rho, (/nl/))
#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif
      grho_s_1d = reshape(grho(:, :, :, :), (/nl/)) * 0.5
      rlrho_s_1d = reshape(rlrho(:, :, :), (/nl/)) * 0.5
      tau_s_1d = reshape(tau(:, :, :), (/nl/)) * 0.5
      j_s_1d = reshape(rj(:, :, :, :), (/nl/)) * 0.5

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
    subroutine exec_libxc_lda()
      implicit none
      real(8) :: rho_1d(nl)
      real(8) :: vx_1d(nl), vc_1d(nl)
      real(8) :: ex_1d(nl), ec_1d(nl)

      rho_1d = reshape(rho, (/nl/))
      if (present(rho_nlcc)) then
        rho_1d = rho_1d + reshape(rho_nlcc, (/nl/))
      endif
      
      call xc_f90_lda_exc_vxc(xc%x_func, nl, rho_1d(1), ex_1d(1), vx_1d(1))
      call xc_f90_lda_exc_vxc(xc%c_func, nl, rho_1d(1), ec_1d(1), vc_1d(1))

      if (present(vxc)) then
         vxc = reshape((vx_1d + vc_1d), (/nx, ny, nz/))
      endif
      
      if (present(exc)) then
         exc = reshape((ex_1d + ec_1d), (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = reshape((ex_1d + ec_1d) * rho_1d, (/nx, ny, nz/))
      endif

      return
    end subroutine exec_libxc_lda
#endif

  end subroutine calc_xc
end module
