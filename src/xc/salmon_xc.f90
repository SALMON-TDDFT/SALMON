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
  use builtin_pz_sp, only: exc_cor_pz_sp
  use builtin_pzm, only: exc_cor_pzm
  use builtin_pbe, only: exc_cor_pbe
  use builtin_tbmbj, only: exc_cor_tbmbj
#ifdef SALMON_USE_LIBXC
  use xc_f90_types_m
  use xc_f90_lib_m
#endif

  implicit none

! List of Exchange Correlation Functionals
  integer, parameter :: salmon_xctype_none  = 0
  integer, parameter :: salmon_xctype_pz    = 1
  integer, parameter :: salmon_xctype_pzm   = 2
  integer, parameter :: salmon_xctype_pbe   = 3
  integer, parameter :: salmon_xctype_tbmbj = 4
  integer, parameter :: salmon_xctype_tpss  = 5
  integer, parameter :: salmon_xctype_vs98  = 6
#ifdef SALMON_USE_LIBXC
  integer, parameter :: salmon_xctype_libxc = 101
#endif

  type xc_functional
    integer :: xctype(3)
    integer :: ispin
    real(8) :: cval
    logical :: use_gradient
    logical :: use_laplacian
    logical :: use_kinetic_energy
    logical :: use_current
#ifdef SALMON_USE_LIBXC
    type(xc_f90_pointer_t) :: func(3)
    type(xc_f90_pointer_t) :: info(3)
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



  subroutine init_xc(xc, ispin, cval, xcname, xname, cname)
    implicit none
    type(xc_functional), intent(inout) :: xc
    integer, intent(in)                :: ispin
    real(8), intent(in)                :: cval
    character(*), intent(in), optional :: xcname
    character(*), intent(in), optional :: xname
    character(*), intent(in), optional :: cname

    ! Initialization of xc variable
    xc%xctype(1:3) = salmon_xctype_none
    xc%ispin = ispin
    xc%cval = cval
    xc%use_gradient = .false.
    xc%use_laplacian = .false.
    xc%use_kinetic_energy = .false.
    xc%use_current = .false.

    ! Exchange correlation
    if (present(xcname) .and. (len_trim(xcname) > 0)) then
      call setup_xcfunc(xcname)
    end if

    ! Exchange only
    if (present(xname) .and. (len_trim(xname) > 0)) then
      call setup_xfunc(xname)
    end if

    ! Correlation only
    if (present(cname) .and. (len_trim(cname) > 0)) then
      call setup_cfunc(cname)
    end if

    return

  contains



    subroutine setup_xcfunc(name)
      use salmon_global, only: iperiodic, domain_parallel
      use inputoutput, only: stop_by_bad_input2
      implicit none
      character(*), intent(in) :: name

#ifdef SALMON_USE_LIBXC
      ! Libxc prefix is used...
      if (lower(name(1:6)) == 'libxc:') then
        xc%xctype(1) = salmon_xctype_libxc
        call init_libxc(name(7:), 1)
        return
      end if
#endif

      select case(lower(name))
      case('none')

        print '(A, A)', "Error! Exchange functional is not specified!"
        stop

        return
      
      case ('pz')
        xc%xctype(1) = salmon_xctype_pz
        return

      case ('pzm')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=pzm')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=pzm')
        end if

        xc%xctype(1) = salmon_xctype_pzm
        return

      case ('pbe')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=pbe')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=pbe')
        end if

        xc%xctype(1) = salmon_xctype_pbe
        xc%use_gradient = .true.
        return

      case ('tbmbj')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=tbmbj')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=tbmbj')
        end if

        xc%xctype(1) = salmon_xctype_tbmbj
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('bj_pw')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=bj_pw')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=bj_pw')
        end if

        xc%xctype(1) = salmon_xctype_tbmbj; xc%cval = 1d0
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('tpss')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=tpss')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=tpss')
        end if

        xc%xctype(1) = salmon_xctype_tpss
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('vs98')
        if(iperiodic==0) then
          call stop_by_bad_input2('iperiodic=0','xc=vs98')
        else if(domain_parallel=='y')then
          call stop_by_bad_input2('iperiodic=3','domain_parallel=y','xc=vs98')
        end if

        xc%xctype(1) = salmon_xctype_vs98
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

#ifdef SALMON_USE_LIBXC
      case('libxc_pz')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ', 3)
        return

      case('libxc_pzm')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ_MOD', 3)
        return

      case('libxc_pbe')
        xc%xctype(2) = salmon_xctype_libxc
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc('GGA_X_PBE', 2)
        call init_libxc('GGA_C_PBE', 3)
        return
#endif

      case default

        print '(A, A)', "Error! Undefined exchange functional:", trim(name)
        stop

      end select
      return
    end subroutine



    subroutine setup_xfunc(name)
      implicit none
      character(*), intent(in) :: name

#ifdef SALMON_USE_LIBXC
      if (lower(name(1:6)) == 'libxc:') then
        xc%xctype(2) = salmon_xctype_libxc
        call init_libxc(name(7:), 2)
        return
      end if
#endif

      select case(name)
      case('none')
        ! xc%xctype(2) = salmon_xctype_none ! default
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

      case default

        print '(A, A)', "Error! Undefined exchange functional:", trim(name)
        stop

      end select
      return
    end subroutine



    subroutine setup_cfunc(name)
      implicit none
      character(*), intent(in) :: name

#ifdef SALMON_USE_LIBXC
      if (lower(name(1:6)) == 'libxc:') then
        xc%xctype(3) = salmon_xctype_libxc
        call init_libxc(name(7:), 3)
        return
      end if
#endif

      select case(name)
      case('none')
        ! xc%xctype(3) = salmon_xctype_none ! default
        return

      ! Please insert additional functional here:
      ! e.g.
      ! case ('additional_functional')
      !   initialization_process_of_functional
      !   return

      case default

        print '(A, A)', "Undefined correlation functional:", trim(name)
        stop

      end select
      return
    end subroutine



#ifdef SALMON_USE_LIBXC
    subroutine init_libxc(libxc_name, ii)
      implicit none
      character(*), intent(in) :: libxc_name
      integer, intent(in) :: ii
      integer ::  ixc

      ixc = xc_f90_functional_get_number(trim(libxc_name))

      if (ixc <= 0) then
        print '(A, A)', "Undefined libxc functional:", trim(libxc_name)
        stop
      end if

      if (ispin > 0) then
        print '(A)', "Spin polarized is not available"
        stop
      end if

      call xc_f90_func_init( &
        & xc%func(ii), xc%info(ii), ixc, XC_UNPOLARIZED &
        & )

      select case(xc_f90_info_family(xc%info(ii)))
      case (XC_FAMILY_LDA)
      case (XC_FAMILY_GGA)
        xc%use_gradient = .true.
      case ( XC_FAMILY_MGGA)
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
      case (XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
        print '(A, A)', "Hybrid is not available:", trim(libxc_name)
        stop
      case default
        print '(A, A)', "Unknown Family:", trim(libxc_name)
        stop
      end select

      return
    end subroutine init_libxc
#endif

  end subroutine init_xc



  subroutine finalize_xc(xc)
    implicit none
    type(xc_functional), intent(inout) :: xc

#ifdef SALMON_USE_LIBXC
    if (xc%xctype(1) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(1))
    if (xc%xctype(2) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(2))
    if (xc%xctype(3) == salmon_xctype_libxc) call xc_f90_func_end(xc%func(3))
#endif

    return
  end subroutine



  subroutine calc_xc(xc, rho, rho_s, exc, eexc, vxc, vxc_s, &
      & grho, grho_s, rlrho, rlrho_s, tau, tau_s, rj, rj_s, &
      & rho_nlcc, &
      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz)
!      & nd, ifdx, ifdy, ifdz, nabx, naby, nabz, Hxyz, aLxyz)
    implicit none
    type(xc_functional), intent(in) :: xc
    real(8), intent(in), optional :: rho(:, :, :) ! ispin = 0
    real(8), intent(in), optional :: rho_s(:, :, :, :) ! ispin = 1
    real(8), intent(out), optional :: exc(:, :, :) ! epsilon_xc[rho]
    real(8), intent(out), optional :: eexc(:, :, :) ! rho * epsilon_xc[rho]
    real(8), intent(out), optional :: vxc(:, :, :) ! v_xc[rho] for ispin=0
    real(8), intent(out), optional :: vxc_s(:, :, :, :) ! v_xc[rho] ispin=1
    !real(8), intent(out), optional :: gvxc(:, :, :) ! v_xc[rho] for ispin=0
    !real(8), intent(out), optional :: gvxc_s(:, :, :, :) ! v_xc[rho] ispin=1
    real(8), intent(in), optional :: grho(:, :, :, :)
    real(8), intent(in), optional :: grho_s(:, :, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: rlrho(:, :, :)
    real(8), intent(in), optional :: rlrho_s(:, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: rj(:, :, :, :)
    real(8), intent(in), optional :: rj_s(:, :, :, :) ! ispin = 1
    real(8), intent(in), optional :: tau(:, :, :)
    real(8), intent(in), optional :: tau_s(:, :, :, :) ! ispin = 1

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
    ! real(8), intent(in), optional :: Hxyz
    ! real(8), intent(in), optional :: aLxyz
    !===============================================================

    integer :: nx, ny, nz, nl

    ! Detect size of 3-dimensional grid
    if (xc%ispin == 0) then
      nx = ubound(rho, 1) - lbound(rho, 1) + 1;
      ny = ubound(rho, 2) - lbound(rho, 2) + 1;
      nz = ubound(rho, 3) - lbound(rho, 3) + 1;
    else
      nx = ubound(rho_s, 1) - lbound(rho_s, 1) + 1;
      ny = ubound(rho_s, 2) - lbound(rho_s, 2) + 1;
      nz = ubound(rho_s, 3) - lbound(rho_s, 3) + 1;
    end if
    nl = nx * ny * nz

    ! Initialize output variables
    if (present(exc)) exc = 0d0
    if (present(eexc)) eexc = 0d0
    if (present(vxc)) vxc = 0d0
    if (present(vxc_s)) vxc_s = 0d0

    ! Exchange-Correlation
    select case (xc%xctype(1))
    case(salmon_xctype_pz)
      call exec_builtin_pz()
    case(salmon_xctype_pzm)
      call exec_builtin_pzm()
    case(salmon_xctype_pbe)
      call exec_builtin_pbe()
    case(salmon_xctype_tbmbj)
      call exec_builtin_tbmbj()
#ifdef SALMON_USE_LIBXC
    case(salmon_xctype_libxc)
      call exec_libxc(1)
#endif
    end select

    ! Exchange Only
    select case (xc%xctype(2))
#ifdef SALMON_USE_LIBXC
    case(salmon_xctype_libxc)
      call exec_libxc(2)
#endif
    end select

    ! Correlation Only
    select case (xc%xctype(3))
#ifdef SALMON_USE_LIBXC
    case(salmon_xctype_libxc)
      call exec_libxc(3)
#endif
    end select

    return

  contains



    subroutine exec_builtin_pz()
      implicit none
      real(8) :: rho_s_1d(nl)
      real(8) :: rho_s_sp_1d(nl,2)
      real(8) :: exc_1d(nl)
      real(8) :: eexc_1d(nl)
      real(8) :: vexc_1d(nl)
      real(8) :: vexc_sp_1d(nl,2)

      if (xc%ispin == 0) then
        rho_s_1d = reshape(rho, (/nl/)) * 0.5
      else if (xc%ispin == 1) then
        rho_s_sp_1d = reshape(rho_s, (/nl,2/))
      end if

#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif

      if (xc%ispin == 0) then
        call exc_cor_pz(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)
      else if (xc%ispin == 1) then
        call exc_cor_pz_sp(nl, rho_s_sp_1d, exc_1d, eexc_1d, vexc_sp_1d)
      end if

      if (xc%ispin == 0) then
        if (present(vxc)) then
           vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
        endif
      else if(xc%ispin == 1) then
        if (present(vxc_s)) then
           vxc_s = vxc_s + reshape(vexc_sp_1d, (/nx, ny, nz,2/))
        endif
      end if

      if (present(exc)) then
         exc = exc + reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
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

#ifndef SALMON_DEBUG_NEGLECT_NLCC
      if (present(rho_nlcc)) then
        rho_s_1d = rho_s_1d + reshape(rho_nlcc, (/nl/)) * 0.5
      endif
#endif

      call exc_cor_pzm(nl, rho_s_1d, exc_1d, eexc_1d, vexc_1d)

      if (present(vxc)) then
         vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
         exc = exc + reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
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
         vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
         exc = exc + reshape(exc_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
         eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
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

      !call exc_cor_tbmbj(nl, rho_1d, rho_s_1d,  grho_s_1d, rlrho_s_1d, tau_s_1d, j_s_1d, xc%cval, eexc_1d, vexc_1d, Hxyz, aLxyz)
      call exc_cor_tbmbj(nl, rho_1d, rho_s_1d,  grho_s_1d, rlrho_s_1d, tau_s_1d, j_s_1d, xc%cval, eexc_1d, vexc_1d)

      if (present(vxc)) then
        vxc = vxc + reshape(vexc_1d, (/nx, ny, nz/))
      endif

      if (present(exc)) then
        ! NOTE: Take care for "zero-division error"
        exc = exc + reshape(eexc_1d, (/nx, ny, nz/)) / rho
      endif

      if (present(eexc)) then
        eexc = eexc + reshape(eexc_1d, (/nx, ny, nz/))
      endif

      return
    end subroutine exec_builtin_tbmbj



#ifdef SALMON_USE_LIBXC
    subroutine exec_libxc(ii)
      implicit none
      integer, intent(in) :: ii
      ! character(256) :: name

      real(8) :: rho_1d(nl), grho_1d(nl, 3)
      real(8) :: sigma_1d(nl), rlrho_1d(nl), tau_1d(nl)
      real(8) :: exc_tmp_1d(nl), vxc_tmp_1d(nl)
      real(8) :: gvxc_tmp_1d(nl), lvxc_tmp_1d(nl), tvxc_tmp_1d(nl)
      
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
      end select

      if (present(exc)) then
        exc = exc + reshape(exc_tmp_1d, (/nx, ny, nz/))
      endif

      if (present(eexc)) then
        eexc = eexc + reshape(exc_tmp_1d, (/nx, ny, nz/)) * rho
      endif

      if (present(vxc)) then
         vxc = vxc + reshape(vxc_tmp_1d, (/nx, ny, nz/))
      endif

      return
    end subroutine exec_libxc
#endif

  end subroutine calc_xc


  
  function lower(s) result(lower_s)
     implicit none
     character(*), intent(in) :: s
     character(len(s)) :: lower_s
     integer :: i
     
     do i = 1, len(s)
       if (('A' <= s(i:i)) .and. (s(i:i) <= 'Z')) then
         lower_s(i:i) = char(ichar(s(i:i)) + 32)
       else
         lower_s(i:i) = s(i:i)
       end if 
     end do
     
     return
   end function lower

end module
